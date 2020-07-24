import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np


def ode_system_factory(a, b, alpha, beta, gamma, sigma):
    return lambda t, x : np.array([
        a * np.sum(x) -b*x[0] - gamma * np.sum(x)*x[0] - beta * x[2] * x[0],
        beta * x[2] * x[0] - b * x[1] - gamma * np.sum(x) * x[1] - sigma * x[1],
        sigma * x[1] - b * x[2] - gamma * np.sum(x) * x[2] - alpha * x[2]
        ])
def r0(a, b, alpha, beta, gamma, sigma):
    return (a-b)*beta * sigma / (gamma *(a + sigma)*(a+alpha))

b = 0.5
a = 1
sigma = 2/13
alpha = 1/73
beta = 11.5/79.69
gamma = 1/100

r = r0(a, b, alpha, beta, gamma, sigma)

f = ode_system_factory(a = a, b = b, alpha = alpha, beta = beta, gamma = gamma, sigma = sigma)
tmax = 100
sol = solve_ivp(f, (0, tmax), [
		(a-b)/gamma,
		100,
		0
	], dense_output=True)
t = np.linspace(0, tmax, 3000)
z = sol.sol(t).T
plt.subplot(121)
plt.plot(t, z)
plt.subplot(122)
plt.plot(t, (z[:,2]+z[:,1]))
plt.show()

# tmax = 4000
# eps = 0.0008
# delta = 0.0005
# xrange = np.arange(0, 1, delta)
# yrange = np.arange(0, 1, delta)


# X, Y = np.meshgrid(xrange, yrange)


# N  = (X - b - alpha * Y) / gamma

# Z1 = 1+ Y * ((alpha - beta *N)/a)  - X
# Z2 = sigma  * (1 - X) / (sigma + alpha + X  - alpha * Y) -Y

# Xbar = np.mean(X[np.logical_and(np.abs(Z2) <eps, np.abs(Z1) <eps)])
# Ybar = np.mean(Y[np.logical_and(np.abs(Z2) <eps, np.abs(Z1) <eps)])
# Nbar = np.mean(N[np.logical_and(np.abs(Z2) <eps, np.abs(Z1) <eps)])

# print(Nbar)
# print(Xbar)
# print(Ybar)


# Xbar = 0.569069
# Ybar = 0.050311
# Nbar =  (Xbar - b - alpha * Ybar) / gamma

# sol = solve_ivp(f, (0, tmax), [
# 		Nbar * Xbar,
# 		Nbar * (1 - Xbar - Ybar),
# 		Nbar * Ybar
# 	], dense_output=True)
# t = np.linspace(0, tmax, 3000)
# z = sol.sol(t).T
# zN = np.sum(z, axis=1)
# zbar = z / np.transpose([zN,zN,zN])






# plt.subplot(121)
# plt.contour(X,Y, Z1, [0])
# plt.contour(X,Y, Z2 , [0])
# plt.plot(zbar[:,0], zbar[:,2])

# plt.subplot(122)
# plt.plot([0,tmax],[zN[-1], zN[-1]])
# plt.plot([0,tmax],[Nbar, Nbar])
# plt.plot(t, z)
# plt.show()



# X = zbar[-1,0]
# Y = zbar[-1,2]
# N  = (X - b - alpha * Y) / gamma
# Z1 = Y * ((beta * N - alpha)/a) - 1 
# Z2 = sigma  * (1 - X) / (sigma + alpha + X  - alpha * Y) -Y
