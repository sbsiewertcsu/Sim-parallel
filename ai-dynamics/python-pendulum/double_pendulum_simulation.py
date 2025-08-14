import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constants
g = 9.81
l1 = 1.0
l2 = 1.0
m1 = 1.0
m2 = 1.0

def double_pendulum_derivs(t, y):
    theta1, omega1, theta2, omega2 = y
    delta = theta2 - theta1
    den1 = (m1 + m2) * l1 - m2 * l1 * np.cos(delta)**2
    den2 = (l2 / l1) * den1

    d_omega1 = (m2 * l1 * omega1**2 * np.sin(delta) * np.cos(delta) +
                m2 * g * np.sin(theta2) * np.cos(delta) +
                m2 * l2 * omega2**2 * np.sin(delta) -
                (m1 + m2) * g * np.sin(theta1)) / den1

    d_omega2 = (-m2 * l2 * omega2**2 * np.sin(delta) * np.cos(delta) +
                (m1 + m2) * g * np.sin(theta1) * np.cos(delta) -
                (m1 + m2) * l1 * omega1**2 * np.sin(delta) -
                (m1 + m2) * g * np.sin(theta2)) / den2

    return [omega1, d_omega1, omega2, d_omega2]

# Initial conditions and integration
y0 = [np.pi / 2, 0, np.pi, 0]
t_eval = np.linspace(0, 20, 1000)
sol = solve_ivp(double_pendulum_derivs, (0, 20), y0, t_eval=t_eval)

theta1, theta2 = sol.y[0], sol.y[2]
x1 = l1 * np.sin(theta1)
y1 = -l1 * np.cos(theta1)
x2 = x1 + l2 * np.sin(theta2)
y2 = y1 - l2 * np.cos(theta2)

# XY plot
plt.figure(figsize=(10, 6))
plt.plot(t_eval, x1, label='x1 (mass 1)', linestyle='-')
plt.plot(t_eval, y1, label='y1 (mass 1)', linestyle='--')
plt.plot(t_eval, x2, label='x2 (mass 2)', linestyle='-')
plt.plot(t_eval, y2, label='y2 (mass 2)', linestyle='--')
plt.title("Double Pendulum XY Position Over Time")
plt.xlabel("Time (seconds)")
plt.ylabel("Position (meters)")
plt.legend()
plt.grid(True)
plt.savefig("xy_position_plot.png")

# Animation
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
line, = ax.plot([], [], 'o-', lw=2)

def update(i):
    line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    return line,

ani = animation.FuncAnimation(fig, update, frames=len(t_eval), interval=20, blit=True)
ani.save("double_pendulum_simulation.mp4", writer='ffmpeg', fps=30)
