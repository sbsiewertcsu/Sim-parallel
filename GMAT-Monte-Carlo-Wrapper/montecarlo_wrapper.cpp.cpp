/**
 * @file montecarlo_wrapper.cpp
 * @brief Monte Carlo simulation wrapper for LEO satellite decay analysis using GMAT and MPI.
 * @author Angel Martinez-Sanchez
 * @date 10/27/2025
 *
 * Varies: initial altitude (h0), drag coefficent (Cd), area-to-mass (A2M)
 * Output per run: lifetime (days) until Altitude < 122 km, or until MaxDays cap.
 *
 * Build:   mpicxx -O3 -std=c++17 -o montecarlo_wrapper montecarlo_wrapper.cpp
 * Run MC:  mpirun -np 4 ./montecarlo_wrapper
 *
 * Options:
 *   --n = Total number of Monte Carlo simulations (default: 100) 
 *   --mass = Satellite mass in kg (default: 200.0)
 *   --capDays = Max days to propagate if no decay (default: 90.0)
*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <filesystem>

#define GMAT_EXECUTABLE "../GmatConsole"
#define GLOBAL_SEED 1234

using namespace std;

struct Trial {
    double h0_km;  // initial altitude
    double Cd;     // drag coefficient
    double A2M;    // area-to-mass (m^2/kg)
};

struct Result {
    int    id;
    double lifetime_days; // (end - start) A1ModJulian days
    double end_alt_km;    // final altitude (should be ~122 if decay hit; higher if capped)
    Trial  trial;
    bool   ok;
};

/******************************************
           Monte Carlo Generators
*******************************************/
static inline double urand() {
    return (rand() / (double)RAND_MAX); // [0,1.0)
}

void generateRandomLEO(Trial& t) {
    t.h0_km = 500.0 + (500.0 * urand());   // 500..1000 km 
    t.Cd    = 2.0   + (0.6   * urand());   // 2.0..2.6
    t.A2M   = 0.005 + (0.045 * urand());   // 0.005..0.05 m^2/kg
}

/******************************************
                 Parsers
*******************************************/

string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == string::npos) return "";
    return s.substr(a, b - a + 1);
}

bool parseTwoLineReport(const string &path, double &startA1, double &startAlt, double &endA1, double &endAlt) {
    ifstream in(path);
    if (!in) return false;
    // Lines contain: A1ModJulian  Altitude
    string line;
    if (!getline(in, line)) return false;
    {
        istringstream iss(line);
        if (!(iss >> startA1 >> startAlt)) return false;
    }
    if (!getline(in, line)) return false;
    {
        istringstream iss(line);
        if (!(iss >> endA1 >> endAlt)) return false;
    }
    return true;
}

string absPath(const string& name) {
    auto cwd = filesystem::current_path();  // returns the current working directory
    filesystem::path p = cwd / name;        // appends filename to cwd
    return p.string();                      // returns full path of file
}

/******************************************
                 Main Logic
*******************************************/

void generateDecayScript(int id, const Trial& t, double massKg, double maxDaysCap) {
    const double Re_km = 6378;          // earth radius in km
    double sma = Re_km + t.h0_km;       // semi-major axis km
    double area_m2 = t.A2M * massKg;

    // absolute paths for report + log
    string csvName = "traj_" + to_string(id) + ".csv";
    string csvAbs  = absPath(csvName);

    string fname = "trajectory_" + to_string(id) + ".script";
    ofstream f(fname);
    auto w = [&](const string& s){ f << s; };   // shorthand labda f() to write text line

    // Spacecraft 
    w("Create Spacecraft S;\n");
    w("S.DateFormat = UTCGregorian;\n");            // human-readable date format
    w("S.Epoch = '01 Jan 2025 12:00:00.000';\n");   // start date/time
    w("S.CoordinateSystem = EarthMJ2000Eq;\n");     // earth-centered inertial frame
    w("S.DisplayStateType = Keplerian;\n");         // method to describe orbit
    
    // 6 Keplerian elements
    f << "S.SMA = " << sma << ";\n";                // semi-major axis (Earth radius + h0)
    w("S.ECC = 0.0;\n");                            // circular orbital
    w("S.INC = 28.5;\n");                           // inclination (degrees)
    w("S.RAAN = 0;\nS.AOP = 0;\nS.TA = 0;\n\n");    // remaining elements (degrees)

    // Drag-related spacecraft properties
    f << "S.DryMass = " << massKg << ";\n";         // spacecraft dry mass (kg)
    f << "S.Cd = " << t.Cd << ";\n";                // drag coefficient
    f << "S.DragArea = " << area_m2 << ";\n\n";     // effective drag area (m^2)

    // --- Force model (drag ON via model name) ---
    w("Create ForceModel FM;\n");
    w("FM.CentralBody = Earth;\n");                 // earth gravity central body
    w("FM.PrimaryBodies = {Earth};\n");
    w("FM.Drag = Exponential;\n");                  // use simple exponential atmosphere model      
    w("FM.SRP  = Off;\n\n");                        // turn off solar radiation pressure

    // --- Propagator ---
    w("Create Propagator Prop;\n");
    w("Prop.FM = FM;\n");
    w("Prop.Type = RungeKutta89;\n");       // Runge-Kutta 8(9) integrator
    w("Prop.InitialStepSize = 60;\n");      // [sec] starting time step (1 minute)
    w("Prop.Accuracy = 1e-12;\n");          // relative integration tolerance
    w("Prop.MinStep = 0.001;\n");           // [sec] minimum time step
    w("Prop.MaxStep = 600;\n\n");           // [sec] maximum time step (10 minutes)

    // --- Report file ---
    w("Create ReportFile R;\n");
    f << "R.Filename = '" << csvAbs << "';\n";  // absolute path for output CSV
    w("R.Precision = 15;\n");                   // decimal precision for output
    w("R.WriteHeaders = false;\n\n");           // don’t write column headers

    // --- Mission sequence ---
    w("BeginMissionSequence;\n");
    // Start row
    w("Report R S.A1ModJulian S.Altitude;\n");  
    // Stop at decay altitude OR at cap days (whichever first)
    f << "Propagate Prop(S) { S.Altitude = 122, S.ElapsedDays = " << maxDaysCap << " };\n";
    // End row
    w("Report R S.A1ModJulian S.Altitude;\n");
}


Result runSingleTrajectory(int id, const Trial& t, double massKg, double maxDaysCap) {
    generateDecayScript(id, t, massKg, maxDaysCap);

    // Make file names and paths
    string script = "trajectory_" + to_string(id) + ".script";  // script filename
    string scriptAbs = absPath(script);                         // absolute path to script
    string log = "traj_" + to_string(id) + ".log";              // log filename

    // Build command executable: "../GmatConsole -r "<script>" > "<log>" 2>&1"
    ostringstream cmd;
    cmd << GMAT_EXECUTABLE << " -r " << "\"" << scriptAbs << "\""   // '-r' tells GMAT to run the script and exit
        << " > " << "\"" << log << "\" 2>&1";                       // '2>&1' redirects stderr to stdout so both go into the log file

    int rc = system(cmd.str().c_str()); // Run GMAT with the generated script

    // Initilizes result
    Result result{};
    result.id = id;
    result.trial = t;
    result.ok = false;

    // Parse 2-line report
    double startA1 = 0, startAlt = 0;   // Line 1 → start time, start altitude
    double endA1 = 0, endAlt = 0;      // Line 2 → end time, end altitude
    string csvAbs = absPath("traj_" + to_string(id) + ".csv");

    // Parses the report file (may fail if integarator did not reach accuracy range)
    if (!parseTwoLineReport(csvAbs, startA1, startAlt, endA1, endAlt)) {
        cerr << "[WARN] Could not parse " << csvAbs << "\n";
        return result;
    }

    result.lifetime_days = endA1 - startA1; // orbital lifetime in days
    result.end_alt_km = endAlt;             // final altitude in km
    result.ok = true;                       // run was successful
    return result;
}

void printResult(const Result& r, int rankTag) {
    cout << "[rank " << rankTag << "] run #" << r.id
         << " lifetime_days=" << r.lifetime_days
         << " end_alt_km=" << r.end_alt_km
         << " | h0=" << r.trial.h0_km
         << " Cd=" << r.trial.Cd
         << " A2M=" << r.trial.A2M
         << "\n";
}

inline bool isBetter(const Result& a, const Result& b, double eps = 1e-9) {
    if (!a.ok) return false;          // a can't beat anything if it's invalid
    if (!b.ok) return true;           // any valid a beats an invalid b

    if (a.lifetime_days > b.lifetime_days + eps) return true;
    if (b.lifetime_days > a.lifetime_days + eps) return false;

    // lifetimes effectively equal → break tie by higher initial altitude
    return a.trial.h0_km > b.trial.h0_km;
}

void runMonteCarloDecay(int numSimulations, int rank, int size, double massKg, double maxDaysCap) {
    Result localBest{};
    localBest.ok = false;

    // Loops through each trial assigned to current rank
    for (int i = rank; i < numSimulations; i += size) {
        srand(GLOBAL_SEED + i);
        Trial t; generateRandomLEO(t);              
        Result r = runSingleTrajectory(i, t, massKg, maxDaysCap);
        // If successful, print and check for best
        if (r.ok) {
            printResult(r, rank);
            if (isBetter(r, localBest)) localBest = r;
        }
    }

    // Gather all results to rank 0
    if (rank == 0) {
        Result globalBest = localBest;

        // Receive results from other processors
        for (int received = 1; received < size; received++) {
            Result rcv{};
            MPI_Status st;
            MPI_Recv(&rcv, sizeof(Result), MPI_BYTE, MPI_ANY_SOURCE, 42, MPI_COMM_WORLD, &st);            
            // Compare to find the global "best" result
            if (isBetter(rcv, globalBest)) globalBest = rcv;
        }

        // Print the overall best trajectory found
        if (globalBest.ok) {
            cout << "[RESULT] Best run was #"
                 << globalBest.id
                 << " lifetime_days=" << globalBest.lifetime_days
                 << " end_alt_km=" << globalBest.end_alt_km
                 << " (h0=" << globalBest.trial.h0_km
                 << ", Cd=" << globalBest.trial.Cd
                 << ", A2M=" << globalBest.trial.A2M << ")\n";
        } else {
            cout << "[RESULT] No successful runs parsed.\n";
        }
    } else {
        // All other ranks send their local best result to rank 0
        MPI_Send(&localBest, sizeof(Result), MPI_BYTE, 0, 42, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int numSim = 50;            // total trials across all ranks
    double massKg = 200.0;      // fixed mass
    double maxDaysCap = 90.0;   // cap if no decay

    // Parse Arguments
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a.rfind("--n=",0)==0) numSim = atoi(a.substr(4).c_str());
        else if (a.rfind("--mass=",0)==0) massKg = atof(a.substr(7).c_str());
        else if (a.rfind("--capDays=",0)==0) maxDaysCap = atof(a.substr(10).c_str());
    }

    // Initial Setup Info
    if (rank == 0) {
        cout << "[INFO] Monte Carlo LEO decay: n=" << numSim
             << " mass=" << massKg << " kg cap=" << maxDaysCap << " days, ranks=" << size << "\n";
    }

    // Run Monte Carlo and Records Time When all Process Finish
    auto t0 = chrono::steady_clock::now();
    runMonteCarloDecay(numSim, rank, size, massKg, maxDaysCap);
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = chrono::steady_clock::now();

    // Outputs Timing Info
    if (rank == 0) {
        double sec = chrono::duration<double>(t1 - t0).count();
        cout << "[TIMING] ranks=" << size << " runtime=" << sec << " s\n";
    }

    MPI_Finalize();
    return 0;
}