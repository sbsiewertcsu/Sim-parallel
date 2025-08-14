#include <SDL2/SDL.h>
#include <cmath>
#include <iostream>

const int WIDTH = 800;
const int HEIGHT = 600;
const double PI = 3.141592653589793;
const double g = 9.81;
const double l1 = 150.0;
const double l2 = 150.0;
const double m1 = 1.0;
const double m2 = 1.0;
const double dt = 0.01;

struct Pendulum {
    double theta1 = PI / 2, omega1 = 0;
    double theta2 = PI, omega2 = 0;

    void update() {
        double delta = theta2 - theta1;
        double den1 = (m1 + m2) * l1 - m2 * l1 * std::cos(delta) * std::cos(delta);
        double den2 = (l2 / l1) * den1;

        double a1 = (m2 * l1 * omega1 * omega1 * std::sin(delta) * std::cos(delta)
                   + m2 * g * std::sin(theta2) * std::cos(delta)
                   + m2 * l2 * omega2 * omega2 * std::sin(delta)
                   - (m1 + m2) * g * std::sin(theta1)) / den1;

        double a2 = (-m2 * l2 * omega2 * omega2 * std::sin(delta) * std::cos(delta)
                   + (m1 + m2) * g * std::sin(theta1) * std::cos(delta)
                   - (m1 + m2) * l1 * omega1 * omega1 * std::sin(delta)
                   - (m1 + m2) * g * std::sin(theta2)) / den2;

        omega1 += a1 * dt;
        omega2 += a2 * dt;
        theta1 += omega1 * dt;
        theta2 += omega2 * dt;
    }

    void get_positions(int &x1, int &y1, int &x2, int &y2) {
        x1 = WIDTH / 2 + static_cast<int>(l1 * std::sin(theta1));
        y1 = HEIGHT / 3 + static_cast<int>(l1 * std::cos(theta1));
        x2 = x1 + static_cast<int>(l2 * std::sin(theta2));
        y2 = y1 + static_cast<int>(l2 * std::cos(theta2));
    }
};

int main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Window* win = SDL_CreateWindow("Double Pendulum Simulation", 100, 100, WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (!win) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren) {
        SDL_DestroyWindow(win);
        SDL_Quit();
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    Pendulum p;
    bool quit = false;
    SDL_Event e;
    while (!quit) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = true;
        }

        p.update();

        int x1, y1, x2, y2;
        p.get_positions(x1, y1, x2, y2);

        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderClear(ren);

        SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
        SDL_RenderDrawLine(ren, WIDTH / 2, HEIGHT / 3, x1, y1);
        SDL_RenderDrawLine(ren, x1, y1, x2, y2);
        SDL_RenderDrawPoint(ren, x1, y1);
        SDL_RenderDrawPoint(ren, x2, y2);

        SDL_RenderPresent(ren);
        SDL_Delay(10);
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
