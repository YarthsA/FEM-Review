# 习题二 解答

---

## 1. 证明 $Q[y(x)] = \int_{a}^{b} y^{2}(x) dx$ 不是线性泛函

**题目**：证明 $Q[y(x)] = \int_{a}^{b} y^{2}(x) dx$ 不是线性泛函。

**解答**：

**线性泛函的定义**：泛函 $Q[y]$ 是线性的，当且仅当对任意函数 $y_1, y_2$ 和常数 $\alpha, \beta$，满足：

$$Q[\alpha y_1 + \beta y_2] = \alpha Q[y_1] + \beta Q[y_2]$$

**验证**：

$$Q[\alpha y_1 + \beta y_2] = \int_a^b (\alpha y_1 + \beta y_2)^2 dx = \alpha^2 \int_a^b y_1^2 dx + 2\alpha\beta \int_a^b y_1 y_2 dx + \beta^2 \int_a^b y_2^2 dx$$

$$= \alpha^2 Q[y_1] + 2\alpha\beta \int_a^b y_1 y_2 dx + \beta^2 Q[y_2]$$

而线性泛函要求：$\alpha Q[y_1] + \beta Q[y_2] = \alpha \int_a^b y_1^2 dx + \beta \int_a^b y_2^2 dx$

比较可知 $Q[\alpha y_1 + \beta y_2] \neq \alpha Q[y_1] + \beta Q[y_2]$（存在交叉项和系数不匹配）。

因此 $Q[y]$ **不是线性泛函**，而是**二次泛函**。 $\blacksquare$

---

## 2. 用变分法求解连接两点之间最短曲线问题

**题目**：用变分的方法求解连接两点之间最短曲线的问题。

**解答**：

**弧长泛函**：$Q[y] = \int_{x_1}^{x_2} \sqrt{1 + y'^2} \, dx$

**Euler 方程**：$F = \sqrt{1 + y'^2}$

$$\frac{\partial F}{\partial y} = 0, \quad \frac{\partial F}{\partial y'} = \frac{y'}{\sqrt{1 + y'^2}}$$

Euler 方程：$0 - \frac{d}{dx}\left(\frac{y'}{\sqrt{1 + y'^2}}\right) = 0$

因此 $\frac{y'}{\sqrt{1 + y'^2}} = C$（常数），解得 $y' = C_1$（常数），即：

$$\boxed{y = C_1 x + C_2}$$

这是**直线方程**，连接两点的最短曲线是直线。

---

## 3. 求泛函 $J[y]=\int_{0}^{1}y^{3}y^{\prime2}dx$ 的一阶变分

**题目**：求泛函 $J[y]=\int_{0}^{1}y^{3}y^{\prime2}dx$ 的一阶变分，$y(0)=1$。

**解答**：

**已知**：$y(0) = 1$（本质边界条件），$y(1)$ 自由。

$$J[y + \varepsilon\delta y] = \int_0^1 (y + \varepsilon\delta y)^3 (y' + \varepsilon\delta y')^2 dx$$

对 $\varepsilon$ 求导并令 $\varepsilon = 0$：

$$\delta J = \int_0^1 \left[3y^2 y'^2 \delta y + 2y^3 y' \delta y'\right] dx$$

对第二项分部积分：

$$\int_0^1 2y^3 y' \delta y' dx = \left[2y^3 y' \delta y\right]_0^1 - \int_0^1 \frac{d}{dx}(2y^3 y') \delta y \, dx$$

**边界项处理**：$x=0$ 处 $\delta y(0) = 0$；$x=1$ 处 $\delta y(1)$ 自由。

**展开导数**：$\frac{d}{dx}(2y^3 y') = 6y^2(y')^2 + 2y^3 y''$

**合并**：

$$\boxed{\delta J = 2y(1)^3 y'(1) \delta y(1) + \int_0^1 \left[-3y^2(y')^2 - 2y^3 y''\right] \delta y \, dx}$$

---

## 4. 求泛函 $J[y] = \int_{0}^{\pi/2}(y'^{2} + 12xy)dx$ 的极值曲线

**题目**：求泛函 $J[y] = \int_{0}^{\pi/2}(y'^{2} + 12xy)dx$ 的极值曲线，边界条件为 $y(0) = 1$，$y(\pi/2) = 1$。

**解答**：

**边界条件**（题目给定）：$y(0) = 1$，$y(\pi/2) = 1$

**Euler 方程**：$F = y'^2 + 12xy$

$$\frac{\partial F}{\partial y} = 12x, \quad \frac{\partial F}{\partial y'} = 2y'$$

Euler 方程：$12x - 2y'' = 0 \Rightarrow y'' = 6x$

积分得：$y' = 3x^2 + C_1$，$y = x^3 + C_1 x + C_2$

**代入边界条件**：

$y(0) = 1$：$C_2 = 1$

$y(\pi/2) = 1$：$(\pi/2)^3 + C_1(\pi/2) + 1 = 1 \Rightarrow C_1 = -\frac{\pi^2}{4}$

**极值曲线**：$\boxed{y = x^3 - \frac{\pi^2}{4}x + 1}$

---

## 5. 求椭球体与平面相交的 $z$ 最大与最小值

**题目**：求椭球体 $16x^{2} + 4y^{2} + z^{2} = 16$ 与平面 $x + y + z = 1$ 相交的 $z$ 最大与最小值。

**解答**：

**Lagrange 乘子法**：构造 $F = z + \lambda_1(x + y + z - 1) + \lambda_2(16x^2 + 4y^2 + z^2 - 16)$

求偏导并令为零：$x = -\frac{\lambda_1}{32\lambda_2}$，$y = -\frac{\lambda_1}{8\lambda_2}$，$z = -\frac{1 + \lambda_1}{2\lambda_2}$

由 $x$ 和 $y$ 的关系：$y = 4x$

代入约束条件得：$21x^2 - 2x - 3 = 0$

解得 $x_1 = \frac{3}{7}$，$x_2 = -\frac{1}{3}$

对应 $z$ 值：$z_1 = -\frac{8}{7}$，$z_2 = \frac{8}{3}$

**结果**：$\boxed{z_{\max} = \frac{8}{3}, \quad z_{\min} = -\frac{8}{7}}$

---

## 6. 推导包含自变函数三阶导数的泛函极值条件的 Euler 方程

**题目**：推导包含自变函数三阶导数的泛函极值条件的 Euler 方程。

**解答**：

**泛函**：$J[y] = \int_{x_1}^{x_2} F(x, y, y', y'', y''') \, dx$

对泛函求变分并分部积分，令 $\delta y$ 的系数为零：

$$\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) + \frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y''}\right) - \frac{d^3}{dx^3}\left(\frac{\partial F}{\partial y'''}\right) = 0}$$

这是**Euler-Poisson 方程**，符号交替出现：$+ - + -$。

---

## 7. 弹性地基梁的微分方程和边界条件

**题目**：泛函 $\Pi_{(w)} = \int_{0}^{L} \left[ \frac{EI}{2} \left( \frac{d^{2}w}{dx^{2}} \right)^{2} + \frac{kw^{2}}{2} + qw \right] dx$，其中 E，I，k 是常数，q 是给定函数，w 是未知函数，试导出原问题的微分方程和边界条件。

**解答**：

**Euler 方程**：$F = \frac{EI}{2}(w'')^2 + \frac{k}{2}w^2 + qw$

$$\frac{\partial F}{\partial w} = kw + q, \quad \frac{\partial F}{\partial w'} = 0, \quad \frac{\partial F}{\partial w''} = EI w''$$

Euler 方程（含四阶导数）：

$$kw + q + \frac{d^2}{dx^2}(EI w'') = 0$$

若 $EI$ 为常数：$\boxed{EI\frac{d^4 w}{dx^4} + kw + q = 0}$

**自然边界条件**：在 $x = 0$ 和 $x = L$ 处，$w'' = 0$（弯矩为零）或 $w$ 给定；$w''' = 0$（剪力为零）或 $w'$ 给定。

---

## 8. 均匀弦的横振动方程

**题目**：试根据变分原理导出完全柔软的均匀弦的横振动方程。

**解答**：

**物理模型**：弦长 $L$，线密度 $\rho$，张力 $T$，横向位移 $u(x,t)$。

**动能**：$T = \frac{1}{2}\int_0^L \rho \left(\frac{\partial u}{\partial t}\right)^2 dx$

**势能**：$V = \frac{1}{2}\int_0^L T \left(\frac{\partial u}{\partial x}\right)^2 dx$

**Hamilton 原理**：$\delta \int_{t_1}^{t_2} (T - V) dt = 0$

对时间部分分部积分：$\int_{t_1}^{t_2} \rho \frac{\partial u}{\partial t} \delta\left(\frac{\partial u}{\partial t}\right) dt = \left[\rho \frac{\partial u}{\partial t} \delta u\right]_{t_1}^{t_2} - \int_{t_1}^{t_2} \rho \frac{\partial^2 u}{\partial t^2} \delta u \, dt$

对空间部分分部积分：$\int_0^L T \frac{\partial u}{\partial x} \delta\left(\frac{\partial u}{\partial x}\right) dx = \left[T \frac{\partial u}{\partial x} \delta u\right]_0^L - \int_0^L T \frac{\partial^2 u}{\partial x^2} \delta u \, dx$

合并并令 $\delta u$ 的系数为零：

$$\boxed{\rho \frac{\partial^2 u}{\partial t^2} = T \frac{\partial^2 u}{\partial x^2}}$$

即一维波动方程，波速 $c = \sqrt{T/\rho}$。
