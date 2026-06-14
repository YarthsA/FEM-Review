# Variation Method & FEA Course — Homework 2 参考答案

> 交大 FEM 课程 · 变分法应用专题
> **对应知识**：[1-3 变分法基础](../../01-Lecture-Notes/1-3-Variational-Methods.md)（§3.4 Euler 方程、§3.5 高阶导数泛函、§3.4.9 Lagrange 乘子法）

---

## 1. 变分法求最短路径

> 📖 [§3.4.2 简单 Euler 方程](../../01-Lecture-Notes/1-3-Variational-Methods.md#342-简单-euler-方程) · [§3.4.5 最速降线（类似题型）](../../01-Lecture-Notes/1-3-Variational-Methods.md#345-最速降线问题的-euler-方程)
>
> 经典变分法问题：对长度泛函列 Euler 方程，解出来一定是直线——和几何直觉一致。

**题目**：求 $A(0,0)$ 到 $B(2,1)$ 的最短曲线。

### 解答

**Step 1**：写出泛函与边界条件

弧长泛函：
$$S[y] = \int_0^2 \sqrt{1 + (y')^2}\,dx$$

边界条件：
$$y(0) = 0,\quad y(2) = 1$$

**Step 2**：写出 Euler 方程

被积函数 $F(x,y,y') = \sqrt{1 + y'^2}$，$F$ 不显含 $x$ 和 $y$。

由 $\frac{\partial F}{\partial y} = 0$，Euler 方程变为：
$$-\frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) = 0$$

所以：
$$\frac{\partial F}{\partial y'} = \text{常数}$$

**Step 3**：计算 $\frac{\partial F}{\partial y'}$

$$\frac{\partial F}{\partial y'} = \frac{y'}{\sqrt{1 + y'^2}} = C$$

**Step 4**：求解 $y'$

$$y' = C\sqrt{1 + y'^2}$$
$$y'^2 = C^2(1 + y'^2)$$
$$y'^2(1 - C^2) = C^2$$
$$y'^2 = \frac{C^2}{1 - C^2}$$

令 $m = \frac{C}{\sqrt{1 - C^2}}$，则 $y' = m$（常数）。

**Step 5**：积分

$$y = mx + b$$

**Step 6**：代入边界条件

$$y(0) = 0 \Rightarrow b = 0$$
$$y(2) = 1 \Rightarrow m \cdot 2 = 1 \Rightarrow m = \frac{1}{2}$$

**结论**：
$$\boxed{y = \frac{1}{2}x}$$

即从 $(0,0)$ 到 $(2,1)$ 的最短路径是过这两点的直线。■

---

## 2. 含三阶导数的 Euler 方程

> 规律：含 n 阶导数的泛函，Euler 方程是 n 阶 ODE，符号交替 $(-1)^k$ 来自反复分部积分。

**题目**：推导含三阶导数泛函 $Q[y] = \int_{x_1}^{x_2} F(x, y, y', y'', y''')dx$ 的极值条件。

### 解答

**Step 1**：泛函的一阶变分

$$\delta Q = \int_{x_1}^{x_2} \left[\frac{\partial F}{\partial y}\delta y + \frac{\partial F}{\partial y'}\delta y' + \frac{\partial F}{\partial y''}\delta y'' + \frac{\partial F}{\partial y'''}\delta y'''\right]dx$$

极值条件 $\delta Q = 0$。

**Step 2**：分部积分

对 $\frac{\partial F}{\partial y'}\delta y'$ 项：
$$\int_{x_1}^{x_2}\frac{\partial F}{\partial y'}\delta y'\,dx = \left[\frac{\partial F}{\partial y'}\delta y\right]_{x_1}^{x_2} - \int_{x_1}^{x_2}\frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right)\delta y\,dx$$

对 $\frac{\partial F}{\partial y''}\delta y''$ 项（分部积分两次）：
$$\int_{x_1}^{x_2}\frac{\partial F}{\partial y''}\delta y''\,dx = \left[\frac{\partial F}{\partial y''}\delta y'\right]_{x_1}^{x_2} - \left[\frac{d}{dx}\left(\frac{\partial F}{\partial y''}\right)\delta y\right]_{x_1}^{x_2} + \int_{x_1}^{x_2}\frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y''}\right)\delta y\,dx$$

对 $\frac{\partial F}{\partial y'''}\delta y'''$ 项（分部积分三次）：
$$\int_{x_1}^{x_2}\frac{\partial F}{\partial y'''}\delta y'''\,dx = \left[\frac{\partial F}{\partial y'''}\delta y''\right]_{x_1}^{x_2} - \left[\frac{d}{dx}\left(\frac{\partial F}{\partial y'''}\right)\delta y'\right]_{x_1}^{x_2} + \left[\frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y'''}\right)\delta y\right]_{x_1}^{x_2} - \int_{x_1}^{x_2}\frac{d^3}{dx^3}\left(\frac{\partial F}{\partial y'''}\right)\delta y\,dx$$

**Step 3**：合并 $\delta y$ 项

假设边界条件满足 $\delta y, \delta y', \delta y''$ 在端点 $x_1, x_2$ 处为零，则边界项全部消失，得：

$$\delta Q = \int_{x_1}^{x_2} \left[\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) + \frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y''}\right) - \frac{d^3}{dx^3}\left(\frac{\partial F}{\partial y'''}\right)\right]\delta y\,dx$$

**Step 4**：由 $\delta y$ 的任意性

$$\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) + \frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y''}\right) - \frac{d^3}{dx^3}\left(\frac{\partial F}{\partial y'''}\right) = 0}$$

**规律**：含 $n$ 阶导数的泛函，Euler 方程为：
$$\sum_{k=0}^n (-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right) = 0$$

其中 $y^{(0)} = y$。■

---

## 3. Lagrange 乘子法求条件极值

> 把约束条件用乘子 λ 引入目标函数，转化为无条件极值问题。这是变分法中处理约束的经典手段，后续 Lagrange 乘子法在有限元中也常用于处理不可压缩约束等。

**题目**：在椭球面 $9x^2 + 4y^2 + z^2 = 36$ 与平面 $2x - y + z = 3$ 的交线上，求 $z$ 的最大值和最小值。

### 解答

**Step 1**：转化为无条件极值

构造 Lagrange 函数：
$$F = z + \lambda_1(2x - y + z - 3) + \lambda_2(9x^2 + 4y^2 + z^2 - 36)$$

**Step 2**：求偏导并令为零

$$\frac{\partial F}{\partial x} = 2\lambda_1 + 18\lambda_2 x = 0 \quad (1)$$
$$\frac{\partial F}{\partial y} = -\lambda_1 + 8\lambda_2 y = 0 \quad (2)$$
$$\frac{\partial F}{\partial z} = 1 + \lambda_1 + 2\lambda_2 z = 0 \quad (3)$$

**Step 3**：由 (1)(2) 解 $x, y$ 用 $\lambda$ 表示

由 (1)：$x = -\frac{\lambda_1}{9\lambda_2}$
由 (2)：$y = \frac{\lambda_1}{8\lambda_2}$

**Step 4**：从 (3) 解 $z$

$$z = -\frac{1 + \lambda_1}{2\lambda_2}$$

**Step 5**：代入约束条件

代入平面方程 $2x - y + z = 3$：
$$2\left(-\frac{\lambda_1}{9\lambda_2}\right) - \left(\frac{\lambda_1}{8\lambda_2}\right) + \left(-\frac{1 + \lambda_1}{2\lambda_2}\right) = 3$$

乘 $72\lambda_2$：
$$-16\lambda_1 - 9\lambda_1 - 36(1 + \lambda_1) = 216\lambda_2$$
$$-16\lambda_1 - 9\lambda_1 - 36 - 36\lambda_1 = 216\lambda_2$$
$$-61\lambda_1 - 36 = 216\lambda_2$$
$$\lambda_2 = -\frac{61\lambda_1 + 36}{216} \quad (4)$$

**Step 6**：代入椭球方程

$$9x^2 + 4y^2 + z^2 = 36$$

$$9\left(-\frac{\lambda_1}{9\lambda_2}\right)^2 + 4\left(\frac{\lambda_1}{8\lambda_2}\right)^2 + \left(-\frac{1 + \lambda_1}{2\lambda_2}\right)^2 = 36$$

$$\frac{\lambda_1^2}{9\lambda_2^2} + \frac{\lambda_1^2}{16\lambda_2^2} + \frac{(1 + \lambda_1)^2}{4\lambda_2^2} = 36$$

$$\frac{1}{\lambda_2^2}\left(\frac{\lambda_1^2}{9} + \frac{\lambda_1^2}{16} + \frac{(1 + \lambda_1)^2}{4}\right) = 36$$

将括号内通分（公分母 144）：
$$\frac{16\lambda_1^2}{144} + \frac{9\lambda_1^2}{144} + \frac{36(1 + \lambda_1)^2}{144} = \frac{25\lambda_1^2 + 36(1 + 2\lambda_1 + \lambda_1^2)}{144}$$
$$= \frac{25\lambda_1^2 + 36 + 72\lambda_1 + 36\lambda_1^2}{144} = \frac{61\lambda_1^2 + 72\lambda_1 + 36}{144}$$

所以：
$$\frac{61\lambda_1^2 + 72\lambda_1 + 36}{144\lambda_2^2} = 36$$

代入 (4)：
$$144 \cdot \frac{(61\lambda_1 + 36)^2}{216^2} = \frac{61\lambda_1^2 + 72\lambda_1 + 36}{36}$$

$$144 \cdot \frac{(61\lambda_1 + 36)^2}{216^2} = \frac{61\lambda_1^2 + 72\lambda_1 + 36}{36}$$

注意到 $\frac{144}{216^2} = \frac{144}{46656} = \frac{1}{324}$：
$$\frac{(61\lambda_1 + 36)^2}{324} = \frac{61\lambda_1^2 + 72\lambda_1 + 36}{36}$$

两边乘 324：
$$(61\lambda_1 + 36)^2 = 9(61\lambda_1^2 + 72\lambda_1 + 36)$$

展开：
$$3721\lambda_1^2 + 4392\lambda_1 + 1296 = 549\lambda_1^2 + 648\lambda_1 + 324$$

$$3172\lambda_1^2 + 3744\lambda_1 + 972 = 0$$

除以 4：
$$793\lambda_1^2 + 936\lambda_1 + 243 = 0$$

**Step 7**：解二次方程

$$\Delta = 936^2 - 4\cdot793\cdot243 = 876096 - 770796 = 105300$$

$$\lambda_1 = \frac{-936 \pm \sqrt{105300}}{2\cdot793} = \frac{-936 \pm 324.5}{1586}$$

$$\lambda_1^{(1)} \approx -0.3856, \quad \lambda_1^{(2)} \approx -0.7952$$

**Step 8**：代入求 $z$

$$z = -\frac{1 + \lambda_1}{2\lambda_2}$$

由 (4)：$\lambda_2 = -\frac{61\lambda_1 + 36}{216}$

$$\lambda_2^{(1)} = -\frac{61(-0.3856) + 36}{216} = -\frac{-23.52 + 36}{216} = -\frac{12.48}{216} \approx -0.05778$$

$$\lambda_2^{(2)} = -\frac{61(-0.7952) + 36}{216} = -\frac{-48.51 + 36}{216} = -\frac{-12.51}{216} \approx 0.05792$$

对第一组 $(\lambda_1^{(1)}, \lambda_2^{(1)})$：
$$z^{(1)} = -\frac{1 - 0.3856}{2(-0.05778)} = -\frac{0.6144}{-0.11556} \approx \boxed{5.317}$$

对第二组 $(\lambda_1^{(2)}, \lambda_2^{(2)})$：
$$z^{(2)} = -\frac{1 - 0.7952}{2(0.05792)} = -\frac{0.2048}{0.11584} \approx \boxed{-1.768}$$

**结论**：
- $z$ 的最大值为 $\approx 5.317$
- $z$ 的最小值为 $\approx -1.768$

（可以通过代入椭球方程验证精确性。）

---

## 4. 变分法求泛函极值函数

> 📖 [§3.4.2 简单 Euler 方程](../../01-Lecture-Notes/1-3-Variational-Methods.md#342-简单-euler-方程) · [§3.4.3 推导证明](../../01-Lecture-Notes/1-3-Variational-Methods.md#343-推导证明)
>
> 标准的 Euler 方程应用题：先写出 F 对 y, y' 的偏导，代入 Euler 方程得到 ODE，再用边界条件定常数。

**题目**：求泛函 $Q[y] = \int_0^1 [(y')^2 + 4y^2 - 8xy]dx$ 的极值函数，边界条件 $y(0)=1, y(1)=2$。

### 解答

**Step 1**：写出 Euler 方程

被积函数：$F = (y')^2 + 4y^2 - 8xy$

$$
\frac{\partial F}{\partial y} = 8y - 8x
$$
$$
\frac{\partial F}{\partial y'} = 2y'
$$

Euler 方程：
$$\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) = 0$$

代入：
$$8y - 8x - \frac{d}{dx}(2y') = 0$$
$$8y - 8x - 2y'' = 0$$

化简：
$$\boxed{y'' - 4y = -4x}$$

**Step 2**：求解 ODE

齐次方程 $y'' - 4y = 0$ 的特征方程：$r^2 - 4 = 0$，$r = \pm 2$

齐次解：
$$y_h = C_1e^{2x} + C_2e^{-2x}$$

**Step 3**：求特解

设特解形式 $y_p = Ax + B$：
$$y_p' = A, \quad y_p'' = 0$$

代入原方程：
$$0 - 4(Ax + B) = -4x$$
$$-4Ax - 4B = -4x$$

比较系数：
$$-4A = -4 \Rightarrow A = 1$$
$$-4B = 0 \Rightarrow B = 0$$

特解：
$$y_p = x$$

**Step 4**：通解

$$y = y_h + y_p = C_1e^{2x} + C_2e^{-2x} + x$$

**Step 5**：代入边界条件

$y(0) = 1$：
$$C_1 + C_2 + 0 = 1 \Rightarrow C_1 + C_2 = 1 \quad (1)$$

$y(1) = 2$：
$$C_1e^{2} + C_2e^{-2} + 1 = 2 \Rightarrow C_1e^2 + C_2e^{-2} = 1 \quad (2)$$

由 (1)：$C_2 = 1 - C_1$

代入 (2)：
$$C_1e^2 + (1 - C_1)e^{-2} = 1$$
$$C_1(e^2 - e^{-2}) + e^{-2} = 1$$
$$C_1 = \frac{1 - e^{-2}}{e^2 - e^{-2}} = \frac{e^2 - 1}{e^4 - 1}$$

化简：
$$C_1 = \frac{(e^2-1)(e^2+1)}{(e^2-1)(e^2+1)e^2}...$$

更简洁的表达式：
$$C_1 = \frac{e^2 - 1}{e^4 - 1} = \frac{e^2 - 1}{(e^2-1)(e^2+1)} = \frac{1}{e^2 + 1}$$

因此 $C_1 = \frac{1}{e^2 + 1}$，$C_2 = 1 - \frac{1}{e^2 + 1} = \frac{e^2}{e^2 + 1}$。

**结论**：
$$\boxed{y(x) = \frac{1}{e^2 + 1}e^{2x} + \frac{e^2}{e^2 + 1}e^{-2x} + x}$$

或写成双曲函数形式：
$$\boxed{y(x) = \frac{1}{\sinh 2}\left[\sinh(2-2x) + \sinh(2x)\right] + x}$$

> **验证**：检查边界条件 $y(0) = \frac{1}{e^2+1} + \frac{e^2}{e^2+1} + 0 = \frac{1+e^2}{e^2+1} = 1$ ✓；$y(1) = \frac{1}{e^2+1}e^2 + \frac{e^2}{e^2+1}e^{-2} + 1 = \frac{e^2+1}{e^2+1} + 1 = 1+1 = 2$ ✓

■
