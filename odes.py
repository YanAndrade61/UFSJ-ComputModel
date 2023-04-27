import numpy as np


def ode_abc_molecules(t, y, k1, k2):
	A, B, C = y

	dABdt = C * k2 - A * B * k1
	dCdt = A * B * k1 - C * k2

	return np.array([ dABdt, dABdt, dCdt ])

def ode_competition(t, y, r1, r2, w11, w12, w21, w22):
	Sa, Sb = y

	dN1dt = Sa * (1 - Sa * w11 - Sb * w21) * r1
	dN2dt = Sb * (1 - Sb * w22 - Sa * w12) * r2

	return np.array([ dN1dt, dN2dt ])

def ode_predator_prey(t, y, a, b, r, m):
	Py, Pr = y

	dHdt = Py * r - Py * Pr * a
	dPdt = Py * Pr * b - Pr * m

	return np.array([ dHdt, dPdt ])

def ode_inhibitive_competition(t, y, k1, k2, k3, k4, k5) -> np.ndarray:
	S, E, I, ES, EI, P = y

	dSdt = ES * k2 \
		 - S * E * k1

	dEdt = ES * k2 \
		 + ES * k3 \
		 + EI * k5 \
		 - S * E * k1 \
		 - E * I * k4
	
	dESdt = S * E * k1 \
		  - ES * k2 \
		  - ES * k3

	dPdt = ES * k3

	dIdt = EI * k5 \
		 - E * I * k4

	dEIdt = E * I * k4 \
		  - EI * k5

	return np.array([ dSdt, dEdt, dIdt, dESdt, dEIdt, dPdt ])

def ode_sirs(t, y, alpha, beta, gamma) -> np.ndarray:
    S, I, R = y

    dSdt = -beta * S * I + alpha * R
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I - alpha * R

    return np.array([ dSdt, dIdt, dRdt ])
