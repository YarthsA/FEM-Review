# 第3章：变分法基础

> **对应 PDF**：[`Chapter 3 Variation theory and applications-1.pdf`](../06-References/pdfs-originals/Chapter%203%20Variation%20theory%20and%20applications-1.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW2 Q1-Q4](../04-Homework-Solutions/2026w/HW2-Problem.md) · [HW3 Q3（弹性地基梁变分推导）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：微积分（微分定义、分部积分公式、一阶/二阶常微分方程求解）、线性代数

---

## 3.1 从函数极值到泛函极值

### 3.1.1 问题的提出

变分法研究的是**如何找到一个函数，使得某个依赖于该函数的量取极值**。

如果你学过微积分，你熟悉的是"函数极值"问题：给定 $y = f(x)$，求 $x$ 使 $y$ 最大/最小。变分法将这个概念推广到函数空间：

| | **函数极值** | **泛函极值** |
|--|------------|------------|
| 输入 | 一个数 $x$ | 一个**函数** $y(x)$ |
| 输出 | 一个数 $y$ | 一个数 $Q$ |
| 问题 | $\frac{dy}{dx} = 0$ | $\delta Q = 0$ |

**例**：最速降线问题（Brachistochrone，Bernoulli 1696）

> 在垂直平面内有两点 $A$ 和 $B$（$A$ 高于 $B$）。一质点在重力作用下从 $A$ 滑到 $B$（无摩擦），问沿哪条曲线所需时间最短？

时间 $T$ 是路径 $y(x)$ 的函数：
$$T[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{2gy}}\,dx$$

我们要求的是**一个函数 $y(x)$**（一条路径），而不是一个数。这就是变分法要解决的问题。

### 3.1.2 发展简史

- **1696** — Bernoulli 提出最速降线问题，Newton、Leibniz、L'Hôpital 等都在数月内给出解答
- **1736** — Euler 给出一般性解法 → **Euler 方程**
- **1755** — Lagrange 改进简化，变分法正式成为数学分支
- **1908** — Ritz 提出直接法（近似解法），绕过困难的 PDE 求解
- **1943** — Courant 分片近似 → 为 FEM 奠定基础
- **1960** — FEM 诞生，变分法成为其理论基石

---

## 3.2 泛函（Functional）

### 3.2.1 定义

泛函是**从函数空间到实数域的映射**。通俗地说：泛函是"函数的函数"——它输入一个函数，输出一个实数。

**定义**：设 $D$ 是一个函数集合。若对任意 $f(x) \in D$，都有唯一的一个实数 $Q$ 与之对应，则称 $Q$ 是 $f$ 的**泛函**，记作 $Q = Q[f(x)]$。

> 与函数的区别：函数 $y = f(x)$ 的取值依赖于自变量 $x$ 在某个**点**的取值；泛函 $Q[f]$ 的取值依赖于函数 $f$ 的**整体形状**（整个函数关系，而非某一点的值）。

**例 1**：$Q[f] = \int_0^1 f(x)\,dx$
- $f(x) = x$ → $Q = 1/2$；$f(x) = \sqrt{x}$ → $Q = 2/3$
- 同样的自变量值 $x=0.5$ 在这两个函数下完全不同，但 $Q$ 关心的是曲线下的**总面积**

**例 2**：$Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx$ ← 这是变分法中最常见的泛函形式

### 3.2.2 线性泛函

若泛函 $Q$ 满足：
1. **齐次性**：$Q[cy] = c\,Q[y]$（$c$ 为任意常数）
2. **可加性**：$Q[y_1 + y_2] = Q[y_1] + Q[y_2]$

则 $Q$ 称为**线性泛函**。两条件可合并为：
$$Q[c_1y_1 + c_2y_2] = c_1Q[y_1] + c_2Q[y_2]$$

**判断练习**：

| 泛函 | 是否线性？ | 原因 |
|------|-----------|------|
| $Q[y] = \int_a^b y(x)dx$ | ✅ 是 | 积分是线性运算 |
| $Q[y] = \int_a^b y^2(x)dx$ | ❌ 不是 | 平方破坏了可加性 |
| $Q[y] = \int_a^b (y'^2 + 4y^2)dx$ | ❌ 不是 | $y$ 以非线性形式出现 |
| $Q[y] = \int_a^b (x^2y + 2y')dx$ | ✅ 是 | 每一项都是 $y$ 或 $y'$ 的线性形式 |

### 3.2.3 泛函的连续性与邻域

函数之间的"接近"有不同的定义方式：

- **零阶接近**：两条曲线的函数值接近（$|y(x)-y_1(x)|$ 小）
- **一阶接近**：函数值和斜率都接近（$|y(x)-y_1(x)|$ 和 $|y'(x)-y_1'(x)|$ 都小）
- **$k$ 阶接近**：直到 $k$ 阶导数都接近

泛函连续的定义：当函数在 $k$ 阶接近的意义下无限接近时，泛函值之差无限小。

---

## 3.3 变分（Variation）

### 3.3.1 自变函数的变分

两个函数之间的差称为**变分**：
$$\delta y = y(x) - y_1(x)$$
$$\delta y' = y'(x) - y_1'(x)$$

自变函数的变分 $\delta y$ 是一个微量函数（类似于自变量的微分 $dx$ 是一个微量）。

> 变分的几何意义：在函数空间中连接 $y(x)$ 和 $y_1(x)$ 的"微小位移"。而且 $\delta y$ 必须在边界上为零（如果端点固定），即 $\delta y(a) = \delta y(b) = 0$。

### 3.3.2 函数的微分回顾

为了类比理解变分，先回顾微分的两种定义：

**定义一（常规）**：增量分解为线性主部 + 高阶小量
$$\Delta y = A(x)\Delta x + \varphi(x,\Delta x)\Delta x^2$$
$$dy = A(x)dx = y'(x)dx$$

**定义二（Lagrange）**：引入参数 $\varepsilon$，在 $\varepsilon=0$ 处对 $\varepsilon$ 求导
$$\frac{\partial}{\partial\varepsilon}y(x + \varepsilon\Delta x)\bigg|_{\varepsilon=0} = y'(x)\Delta x = dy$$

### 3.3.3 泛函变分的两种定义

**定义一（常规）**：类比函数的微分

泛函的增量：
$$\Delta Q = Q[y+\delta y] - Q[y]$$

如果能将 $\Delta Q$ 分解为：
$$\Delta Q = T[y(x), \delta y] + o(\delta y)$$

其中 $T$ 对 $\delta y$ 是线性的（$y$ 固定时），则 $\delta Q = T$ 称为泛函的一阶变分。

**例**：对 $Q[y] = \int_a^b y^2(x)dx$
$$\Delta Q = \int_a^b [(y+\delta y)^2 - y^2]dx = \int_a^b (2y\delta y + \delta y^2)dx$$

第一项 $\int_a^b 2y\delta y\,dx$ 对 $\delta y$ 是线性的 → 这就是 $\delta Q$

**定义二（Lagrange 法，计算利器）**：
$$\delta Q = \left.\frac{\partial}{\partial\varepsilon} Q[y + \varepsilon\delta y]\right|_{\varepsilon=0}$$

> **为什么有用？** 因为将变分问题转化为普通微分问题——对 $\varepsilon$ 求导即可，避免了每次手动分解。

### 3.3.4 变分算子的运算法则

$$ \delta(y') = (\delta y)' $$

即**微分算子和变分算子可交换顺序**。这个性质在推导 Euler 方程时非常方便，因为它允许我们在分部积分时将 $\delta y'$ 视为 $(\delta y)'$。

---

## 3.4 泛函极值的必要条件

类比函数极值 $dy=0$，泛函极值的**必要条件**是：
$$\delta Q = 0$$
在 $y=y_0$ 处。

也就是说：若 $y_0(x)$ 使 $Q[y]$ 取极值，则在该函数处泛函的一阶变分为零。

> 这只是必要条件而非充分条件（类似于 $f'(x)=0$ 只是极值候选点）。对于实际问题，物理背景通常保证了充分性。

---

## 3.5 Euler 方程的完整推导

这是变分法的核心，也是考试的**必考点**。请仔细理解每一步。

### 3.5.1 问题设定

求泛函 $Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx$ 在边界条件 $y(x_1)=y_1, y(x_2)=y_2$ 下的极值函数。

### 3.5.2 推导过程

**Step 1：构造邻域函数**

设 $y(x)$ 是使 $Q$ 取极值的函数。构造"扰动"函数：
$$y_1(x) = y(x) + \alpha\eta(x)$$

其中 $\eta(x)$ 是任意满足 $\eta(x_1)=\eta(x_2)=0$ 的光滑函数，$\alpha$ 是一个小参数。

当 $\alpha=0$ 时，$y_1 = y$，$Q$ 取极值。

**Step 2：将泛函转化为关于 $\alpha$ 的普通函数**

$$\varphi(\alpha) = Q[y+\alpha\eta] = \int_{x_1}^{x_2} F(x, y+\alpha\eta, y'+\alpha\eta')\,dx$$

既然 $Q$ 在 $\alpha=0$ 取极值，则 $\varphi(\alpha)$ 在 $\alpha=0$ 取极值，所以：
$$\varphi'(0) = 0$$

**Step 3：求导**

在积分号内对 $\alpha$ 求导，应用链式法则：
$$\varphi'(0) = \int_{x_1}^{x_2} \left[ \frac{\partial F}{\partial y}\eta + \frac{\partial F}{\partial y'}\eta' \right]dx = 0$$

**Step 4：分部积分**

对第二项分部积分：
$$\int_{x_1}^{x_2} \frac{\partial F}{\partial y'}\eta'\,dx = \left[ \frac{\partial F}{\partial y'}\eta \right]_{x_1}^{x_2} - \int_{x_1}^{x_2} \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right)\eta\,dx$$

由于 $\eta(x_1)=\eta(x_2)=0$，边界项为零，得：
$$\int_{x_1}^{x_2} \left[ \frac{\partial F}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right) \right]\eta(x)\,dx = 0$$

**Step 5：应用变分法预备定理**

> **预备定理**：若 $f(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 的光滑函数 $\eta(x)$ 都有 $\int_a^b f(x)\eta(x)dx = 0$，则 $f(x) \equiv 0$ 在 $[a,b]$ 上。

令 $f(x) = \frac{\partial F}{\partial y} - \frac{d}{dx}(\frac{\partial F}{\partial y'})$，则由预备定理：

$$\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right) = 0}$$

这就是 **Euler-Lagrange 方程**（通常简称为 Euler 方程）。

### 3.5.3 推导要点重述

1. **构造扰动** $\to$ 变分问题转为 $\alpha$ 的函数
2. **求导 + 分部积分** $\to$ 将 $\delta y'$ 转移成 $\delta y$
3. **利用边界条件** $\to$ 取消边界项
4. **$\delta y$ 的任意性** $\to$ 被积函数必须为零

---

## 3.6 Euler 方程的应用（解题模板）

### 标准流程

```
① 写出被积函数 F(x, y, y')
② 计算 ∂F/∂y 和 ∂F/∂y'
③ 计算 d/dx(∂F/∂y')
④ 代入 Euler 方程 → 得到 ODE
⑤ 代入边界条件 → 求积分常数
```

### 例 1：最短路径

求 $A(0,0)$ 到 $B(1,1)$ 的最短曲线。

**解**：弧长泛函 $S[y] = \int_0^1 \sqrt{1+y'^2}\,dx$，$F = \sqrt{1+y'^2}$

$$\frac{\partial F}{\partial y} = 0,\quad \frac{\partial F}{\partial y'} = \frac{y'}{\sqrt{1+y'^2}}$$

Euler 方程：
$$0 - \frac{d}{dx}\left(\frac{y'}{\sqrt{1+y'^2}}\right) = 0 \Rightarrow \frac{y'}{\sqrt{1+y'^2}} = C$$

解得 $y' = \text{常数}$ → $y = C_1 x + C_2$，代入边界条件得 $y = x$。■

### 例 2：$F$ 不显含 $x$ 的首次积分

当 $F = F(y, y')$（不显含 $x$）时，有**首次积分**：
$$F - y'\frac{\partial F}{\partial y'} = C$$

**推导**：
$$\frac{d}{dx}(F - y'F_{y'}) = F_y y' + F_{y'}y'' - y''F_{y'} - y'(F_{y'y}y' + F_{y'y'}y'') = y'(F_y - \frac{d}{dx}F_{y'}) = 0$$

这个技巧在解最速降线等问题时特别有用。

### 例 3：完整 ODE 求解

求 $Q[y] = \int_0^1 (y'^2 + 4y^2 - 8xy)dx$ 的极值函数，$y(0)=1, y(1)=2$。

$F = y'^2 + 4y^2 - 8xy$，$F_y = 8y - 8x$，$F_{y'} = 2y'$

Euler 方程：$8y - 8x - \frac{d}{dx}(2y') = 0$ → $y'' - 4y = -4x$

这是一个非齐次二阶线性 ODE：
- 齐次解：$y_h = C_1e^{2x} + C_2e^{-2x}$（特征方程 $r^2-4=0$）
- 特解（待定系数法）：设 $y_p = Ax + B$，代入得 $A=1, B=0$
- 通解：$y = C_1e^{2x} + C_2e^{-2x} + x$

代入 $y(0)=1$：$C_1 + C_2 = 1$
代入 $y(1)=2$：$C_1e^2 + C_2e^{-2} + 1 = 2$ → $C_1e^2 + C_2e^{-2} = 1$

解得 $C_1 = \frac{1}{e^2+1}$，$C_2 = \frac{e^2}{e^2+1}$。

$$y(x) = \frac{e^{2x}}{e^2+1} + \frac{e^2 \cdot e^{-2x}}{e^2+1} + x = \frac{\sinh(2-2x) + \sinh 2x}{\sinh 2} + x$$

验证：$y(0) = \frac{e^0 + e^2}{e^2+1} + 0 = 1$ ✅，$y(1) = \frac{e^2 + e^0}{e^2+1} + 1 = 2$ ✅

---

## 3.7 Euler 方程的推广

### 3.7.1 含二阶导数的泛函

$Q[y] = \int F(x, y, y', y'')dx$

$$\delta Q = \int (F_y\delta y + F_{y'}\delta y' + F_{y''}\delta y'')dx = 0$$

对 $F_{y''}\delta y''$ 项需要分部积分**两次**：
$$\int F_{y''}\delta y'' dx = F_{y''}\delta y' - \frac{d}{dx}F_{y''}\delta y + \int\frac{d^2}{dx^2}F_{y''}\delta y\,dx$$

最终：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0}$$

### 3.7.2 含 $n$ 阶导数（通用形式）

$$Q[y] = \int F(x, y, y', \ldots, y^{(n)})dx$$

$$\boxed{\sum_{k=0}^n (-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right) = 0}$$

规律：$+F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} - \frac{d^3}{dx^3}F_{y'''} + \cdots$（交替正负）

### 3.7.3 多个独立函数的泛函

$$Q[y_1, \ldots, y_n] = \int F(x, y_1, \ldots, y_n, y_1', \ldots, y_n')dx$$

对每个 $y_i$ 有独立的 Euler 方程：
$$F_{y_i} - \frac{d}{dx}F_{y_i'} = 0 \quad (i=1,2,\ldots,n)$$

### 3.7.4 多元函数的泛函（多变量）

$$Q[z(x,y)] = \iint_D F(x, y, z, z_x, z_y)\,dxdy$$

$$\boxed{F_z - \frac{\partial}{\partial x}F_{z_x} - \frac{\partial}{\partial y}F_{z_y} = 0}$$

---

## 3.8 边界条件分类

### 本质边界条件（Essential BC）
边界值**预先给定**。在变分中，$\delta y$ 在边界上为零。
- 例：固定端 $y(a) = y_a$，简支端 $w(0) = 0$

### 自然边界条件（Natural BC）
边界值**未给定**，由分部积分产生的边界项必须单独为零：
$$\left.\frac{\partial F}{\partial y'}\right|_{x=a} = 0$$
这种边界条件是"自然"出现的，无需预先指定。

### 无条件极值 vs 条件极值
- **无条件极值**：函数满足边界条件外无其他约束
- **条件极值**：附加等周条件 $\int \varphi\,dx = \alpha$ → 用 **Lagrange 乘子法**

---

## 3.9 条件极值与 Lagrange 乘子法

在约束 $\int_a^b \varphi(x,y,y')dx = \alpha$ 下求 $Q[y] = \int_a^b F dx$ 的极值。

构造新泛函：
$$Q^*[y] = \int_a^b (F + \lambda\varphi)dx - \lambda\alpha = \int_a^b F^* dx - \lambda\alpha$$

其中 $F^* = F + \lambda\varphi$，$\lambda$ 为 Lagrange 乘子（额外的未知常数）。

对 $Q^*$ 使用 Euler 方程：$F^*_y - \frac{d}{dx}F^*_{y'} = 0$

通解含 2 个积分常数 + $\lambda$，由 2 个边界条件 + 1 个等周条件确定。

---

## 3.10 最速降线问题（完整求解）

### 物理模型
质点在重力作用下沿光滑曲线从 $A$ 到 $B$，求时间最短的路径。

### 建立泛函
能量守恒：$\frac12 mv^2 = mgy \Rightarrow v = \sqrt{2gy}$
弧长元：$dS = \sqrt{1+y'^2}dx$
时间元：$dt = dS/v = \sqrt{(1+y'^2)/(2gy)}\,dx$

总时间（略去常数 $1/\sqrt{2g}$）：
$$Q[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx,\quad F(y,y') = \sqrt{\frac{1+y'^2}{y}}$$

### 利用首次积分
由于 $F$ 不显含 $x$：
$$F - y'F_{y'} = C$$

计算 $F_{y'} = \frac{y'}{\sqrt{y(1+y'^2)}}$，代入：
$$\sqrt{\frac{1+y'^2}{y}} - \frac{y'^2}{\sqrt{y(1+y'^2)}} = \frac{1}{\sqrt{y(1+y'^2)}} = C$$

即：
$$y(1+y'^2) = D \quad (D = 1/C^2)$$

### 参数化求解
令 $y' = \tan\theta$，代入得：
$$y = \frac{D}{1+\tan^2\theta} = D\cos^2\theta = \frac{D}{2}(1+\cos 2\theta)$$

$$dx = \frac{dy}{y'} = \frac{-D\sin 2\theta d\theta}{\tan\theta} = -D(1+\cos 2\theta)d\theta$$

积分：
$$\begin{cases}
x = -\frac{D}{2}(2\theta + \sin 2\theta) + E \\
y = \frac{D}{2}(1 + \cos 2\theta)
\end{cases}$$

令 $2\theta = \pi - \phi$，得到旋轮线（cycloid）的标准形式：
$$\boxed{\begin{cases}
x = r(\phi - \sin\phi) \\
y = r(1 - \cos\phi)
\end{cases}},\quad r = \frac{D}{2}$$

这就是最速降线的解——一条旋轮线！

---

## 3.11 解题常见陷阱

| 陷阱 | 说明 |
|------|------|
| 混淆偏导与全导 | Euler 方程中 $\partial F/\partial y'$ 是偏导，但 $d/dx$ 是全导：$\frac{d}{dx}F_{y'} = F_{y'x} + F_{y'y}y' + F_{y'y'}y''$ |
| 忽略边界项 | 分部积分后边界项 $[F_{y'}\eta]_{x_1}^{x_2}$ 必须处理。$\eta(x_1)=\eta(x_2)=0$ 时消失；否则自然 BC |
| 线性泛函判断不全 | 必须同时检查齐次性和可加性 |
| 条件极值中 $\lambda$ 的处理 | $\lambda$ 是额外未知数，$F^* = F + \lambda\varphi$ 代入 Euler 方程，**同时**满足等周条件 |

---

## 检查你的理解

1. 用自己的话解释：什么是泛函？变分 $\delta Q$ 和微分 $dy$ 有什么类比关系？
2. 写出 Euler 方程，并解释推导过程中每一步的数学依据。
3. 对于 $Q[y] = \int (y''^2 + 2y^2)dx$，它的 Euler 方程是什么？
4. 什么情况下可以用首次积分 $F - y'F_{y'} = C$？
5. 本质边界条件和自然边界条件的区别是什么？分别举例。

---

> **对应作业**：[HW2 Q1（最短路径）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q2（三阶导数 Euler 方程）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q3（Lagrange 乘子）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q4（泛函极值函数）](../04-Homework-Solutions/2026w/HW2-Problem.md)
> **往年参考**：[past/HW2/homework 2](../04-Homework-Solutions/past/HW2/homework%202.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW2/Ans%20to%20HM2_LIU%20Sai_handed%20in.md)
