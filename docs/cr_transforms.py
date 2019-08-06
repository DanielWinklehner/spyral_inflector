from matplotlib import pyplot as plt
import numpy as np
from dans_pymodules import MyColors


def homogenize(v):
    return np.array([v[0], v[1], 1.0])


def operate(v, T):
    v = homogenize(v)
    vp = np.matmul(T, v)
    return vp[:2]


def make_bounds(v1, v2, ax, col):
    t = np.linspace(-10, 10, 100)
    va = np.array([-v2[1] * t + v2[0] + v1[0], v2[0] * t + v2[1] + v1[1]]).T
    vb = np.array([-v2[1] * t + v1[0], v2[0] * t + v1[1]]).T

    ax.plot(va[:, 0], va[:, 1], '--', color=col)
    ax.plot(vb[:, 0], vb[:, 1], '--', color=col)


colors = MyColors()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-0.15, 0.15])
ax.set_ylim([-0.15, 0.15])
ax.grid(True)
ax.set_aspect(1)

dee_opening_angle = np.deg2rad(42.5)
dee_angle = dee_opening_angle / 2.0

d1 = 0.05
d2 = 0.05
point = np.array([0.07, 0.01])
# dth = np.deg2rad(25.0)
dth = 0.0

a1 = d1*np.array([np.cos(dee_angle), np.sin(dee_angle)])
a2 = np.array([a1[0]*np.cos(dth) - a1[1]*np.sin(dth), a1[0]*np.sin(dth) + a1[1]*np.cos(dth)])

ra = a1
rb = a2 + a1
dr = rb - ra
s = np.linalg.norm(dr)
dth = -np.arccos(dr[0] / s) * np.sign(dr[1])
dx, dy = -ra
T = np.array([[1.0, 0.0, dx], [0.0, 1.0, dy], [0.0, 0.0, 1.0]])
S = np.array([[1.0 / s, 0.0, 0.0], [0.0, 1.0 / s, 0.0], [0.0, 0.0, 1.0]])
R = np.array([[np.cos(dth), -np.sin(dth), 0.0], [np.sin(dth), np.cos(dth), 0.0], [0.0, 0.0, 1.0]])

ax.plot(np.array([a1, a2 + a1])[:, 0], np.array([a1, a2 + a1])[:, 1], color=colors[0], linewidth=2)
make_bounds(a1, a2, ax, col=colors[0])
ax.scatter(point[0], point[1], marker='X', color=colors[1])

plt.title(r'$\vec{r}$')
# plt.savefig('no_transform.png')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-0.15, 0.15])
ax.set_ylim([-0.15, 0.15])
ax.grid(True)
ax.set_aspect(1)

a1 = operate(a1, T)
point = operate(point, T)

ax.plot(np.array([a1, a2 + a1])[:, 0], np.array([a1, a2 + a1])[:, 1], color=colors[0], linewidth=2)
make_bounds(a1, a2, ax, col=colors[0])
ax.scatter(point[0], point[1], marker='X', color=colors[1])

plt.title(r'$\mathbf{T}\vec{r}$')
# plt.savefig('translated.png')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-0.15, 0.15])
ax.set_ylim([-0.15, 0.15])
ax.grid(True)
ax.set_aspect(1)

a1 = operate(a1, R)
a2 = operate(a2, R)
point = operate(point, R)

ax.plot(np.array([a1, a2 + a1])[:, 0], np.array([a1, a2 + a1])[:, 1], color=colors[0], linewidth=2)
make_bounds(a1, a2, ax, col=colors[0])
ax.scatter(point[0], point[1], marker='X', color=colors[1])

plt.title(r'$\mathbf{RT}\vec{r}$')
# plt.savefig('rotated.png')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([-2.0, 2.0])
ax.set_ylim([-2.0, 2.0])
ax.grid(True)
ax.set_aspect(1)

a1 = operate(a1, S)
a2 = operate(a2, S)
point = operate(point, S)

ax.plot(np.array([a1, a2 + a1])[:, 0], np.array([a1, a2 + a1])[:, 1], color=colors[0], linewidth=2)
make_bounds(a1, a2, ax, col=colors[0])
ax.scatter(point[0], point[1], marker='X', color=colors[1])

plt.title(r'$\mathbf{SRT}\vec{r}$')
# plt.savefig('scaled.png')
plt.show()
