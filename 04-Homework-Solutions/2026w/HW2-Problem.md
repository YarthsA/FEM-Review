# Variation Method & FEA Course Homework 2

> **对应知识**：[1-3 变分法基础](../../01-Lecture-Notes/1-3-Variational-Methods.md)（§3.4 Euler 方程、§3.5 含高阶导数泛函、§3.4.9 Lagrange 乘子法）

Submitted before May 25, 2026

1. Identify the shortest curve between the point $A(0,0)$ and the point $B(2,1)$ by the variation method.

Hints: Let the curve be expressed as $y = y(x)$ . The length functional can be written as

$$
S [ y ] = \int_ {0} ^ {2} \sqrt {1 + (y ^ {\prime}) ^ {2}} d x,
$$

with boundary conditions

$$
y (0) = 0, \qquad y (2) = 1.
$$

2. Deduce the Euler equation for the extremum condition of the functional containing third-order derivative of its independent function.

3. Solve the maximum and minimum values of $z$ on the intersection curve of the ellipsoid

$$
9 x ^ {2} + 4 y ^ {2} + z ^ {2} = 3 6
$$

and the plane

$$
2 x - y + z = 3.
$$

Hints: Consider the unconditional extremum of the function

$$
F = z + \lambda_ {1} (2 x - y + z - 3) + \lambda_ {2} (9 x ^ {2} + 4 y ^ {2} + z ^ {2} - 3 6),
$$

where $\lambda_{1}$ and $\lambda_{2}$ are Lagrange multipliers.

4. Find the function $y(x)$ leading to the extremum of the functional

$$
Q [ y ] = \int_ {0} ^ {1} [ (y ^ {\prime}) ^ {2} + 4 y ^ {2} - 8 x y ] d x,
$$

with boundary conditions

$$
y (0) = 1, \qquad y (1) = 2.
$$

Hints: First write the Euler equation corresponding to the functional, then solve the resulting ordinary differential equation using the boundary conditions.