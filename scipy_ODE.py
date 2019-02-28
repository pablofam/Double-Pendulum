#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

l1 = 1
m1 = 1
l2 = 0.5
m2 = 0.5
g = 9.81

# state is current state of system: [theta1, omega1, theta2, omega2]

def derivatives(state, stupid):
    theta1, omega1, theta2, omega2 = state
    dtheta1 = omega1
    dtheta2 = omega2

    domega1 = (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2) - 2 * np.sin(theta1 - theta2) * m2 * (omega2 ** 2 * l2 + omega1 ** 2 * l1 * np.cos(theta1 - theta2))) / (l2 * (2 * m1 + m2 - m2 * np.cos( 2 * theta1 - 2 * theta2)))
    domega2 = (2 * np.sin(theta1 - theta2) * (omega1 ** 2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1) + omega2 ** 2 * l2 * m2 * np.cos(theta1 - theta2))) / (l2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

    return[dtheta1, domega1, dtheta2, domega2]


# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 250

t = np.linspace(0, stoptime, numpoints)

# initial conditions

initial = (1,0,0,0)

# Call the ODE solver.
wsol = scipy.integrate.odeint(derivatives, initial, t, atol=abserr, rtol=relerr)

print(wsol)

plt.plot(t, [i[0] for i in wsol], "g")
plt.plot(t, [i[2] for i in wsol], "r")
plt.show()
