from matplotlib import pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlim([-0.15, 0.15])
ax.set_ylim([-0.15, 0.15])
ax.grid(True)

dee_opening_angle = np.deg2rad(42.5)

d1 = 0.05
d2 = 0.05

dth = np.deg2rad(10.0)

a1 = d1*np.array([np.cos(dee_opening_angle / 2.0), np.sin(dee_opening_angle / 2.0)])
a2 = np.array([a1[0]*np.cos(dth) - a1[1]*np.sin(dth), a1[0]*np.sin(dth) + a1[1]*np.cos(dth)])
# a2 = d2*np.array([np.cos(dee_opening_angle / 2.0), np.sin(dee_opening_angle / 2.0)])
b1 = d1*np.array([np.cos(-dee_opening_angle / 2.0), np.sin(-dee_opening_angle / 2.0)])
b2 = d2*np.array([np.cos(-dee_opening_angle / 2.0), np.sin(-dee_opening_angle / 2.0)])

zv = np.array([0.0, 0.0])

t = np.linspace(-10, 10, 100)
xa1t = np.array([-a1[1] * t + a1[0], a1[0] * t + a1[1]]).T
xa1 = np.array([-a1[1] * t, a1[0] * t]).T

xa2t = np.array([-a2[1] * t + a2[0] + a1[0], a2[0] * t + a2[1] + a1[1]]).T
xa2 = np.array([-a2[1] * t + a1[0], a2[0] * t + a1[1]]).T


xb1t = np.array([-b1[1] * t + b1[0], b1[0] * t + b1[1]]).T
xb1 = np.array([-b1[1] * t , b1[0] * t]).T

xb2t = np.array([-b2[1] * t + b2[0] + b1[0], b2[0] * t + b2[1] + b1[1]]).T
xb2 = np.array([-b2[1] * t + b1[0], b2[0] * t + b1[1]]).T


ax.plot(np.array([zv, a1])[:, 0], np.array([zv, a1])[:, 1], color='b', linewidth=2)
ax.plot(np.array([a1, a2 + a1])[:, 0], np.array([a1, a2 + a1])[:, 1], color='m', linewidth=2)
ax.plot(xa1[:, 0], xa1[:, 1], '--', color='b')
ax.plot(xa1t[:, 0], xa1t[:, 1], '--', color='b')
ax.plot(xa2[:, 0], xa2[:, 1], '--', color='m')
ax.plot(xa2t[:, 0], xa2t[:, 1], '--', color='m')

ax.plot(np.array([zv, b1])[:, 0], np.array([zv, b1])[:, 1], color='r', linewidth=2)
ax.plot(np.array([zv, b2+b1])[:, 0], np.array([zv, b2+b1])[:, 1], color='r', linewidth=2)
ax.plot(xb1[:, 0], xb1[:, 1], '--', color='r')
ax.plot(xb1t[:, 0], xb1t[:, 1], '--', color='r')
# ax.plot(xb2[:, 0], xb2[:, 1], '--', color='r')
ax.plot(xb2t[:, 0], xb2t[:, 1], '--', color='r')

ang = np.linspace(0.0, 359.0, 360)
ang = np.deg2rad(ang)
ax.plot(d1*np.cos(ang), d1*np.sin(ang), ':', color='g')
ax.plot((d1+d2)*np.cos(ang), (d1+d2)*np.sin(ang), ':', color='g')

ax.set_aspect(1)


ra = a1
rb = a2 + a1
dr = rb - ra
s = np.linalg.norm(dr)
dth = -np.arccos(dr[0] / s)
dx, dy = -ra
T = np.array([[1.0, 0.0, dx], [0.0, 1.0, dy], [0.0, 0.0, 1.0]])
S = np.array([[1.0 / s, 0.0, 0.0], [0.0, 1.0 / s, 0.0], [0.0, 0.0, 1.0]])
R = np.array([[np.cos(dth), -np.sin(dth), 0.0], [np.sin(dth), np.cos(dth), 0.0], [0.0, 0.0, 1.0]])

# r = np.array([0.05, 0.08, 1.0])
r = np.array([ra[0] + dr[0] * 0.5, ra[1] + dr[1] * 0.5, 1.0])
tr = np.matmul(S, np.matmul(R, np.matmul(T, r)))[:2]
# tr = np.matmul(R, np.matmul(T, r))[:2]

print(tr)
ax.scatter(tr[0], tr[1])
ax.scatter(r[0], r[1], marker='X')
# ax.scatter(a1[0], a1[1], marker='X')
plt.show()
