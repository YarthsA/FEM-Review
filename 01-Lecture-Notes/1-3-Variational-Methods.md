# 第3章：变分法基础（Variational Calculus）

> **对应 PDF**：[`Chapter 3 Variation theory and applications-1.pdf`](../06-References/pdfs-originals/Chapter%203%20Variation%20theory%20and%20applications-1.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW2 Q1-Q4](../04-Homework-Solutions/2026w/HW2-Problem.md) · [HW3 Q3（弹性地基梁）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：微积分（微分定义、分部积分、常微分方程求解）、线性代数

---

## 3.1 引言：从最速降线问题说起

### 3.1.1 一个改变数学史的问题

1696 年 6 月，瑞士数学家 Johann Bernoulli 在《教师学报》上向全欧洲数学家发起了挑战。他提出了这样一个问题：

> 设 $A$ 和 $B$ 是垂直平面内不在同一垂线上的两点（$A$ 高于 $B$）。一质点在重力作用下从 $A$ 无摩擦地滑到 $B$。问沿哪条曲线滑行所需时间最短？

这就是著名的**最速降线问题（Brachistochrone Problem）**，源自希腊语 "brachistos"（最短）和 "chronos"（时间）。

这个问题在当时引起了轰动。Newton、Leibniz、L'Hôpital 和 Jacob Bernoulli（Johann 的哥哥）都在数月内给出了解答。令人惊讶的是，答案不是人们直觉以为的直线或圆弧，而是一条**旋轮线（cycloid）**——即一个圆沿直线滚动时，圆周上某点所描绘的轨迹。

### 3.1.2 最速降线的数学表述

建立坐标系：$A$ 在原点，$y$ 轴向下为正。质点质量为 $m$，重力加速度为 $g$。

由能量守恒：
$$\frac12 mv^2 = mgy \quad\Rightarrow\quad v = \sqrt{2gy}$$

弧长微分：$dS = \sqrt{1 + y'^2}\,dx$

时间微元：
$$dt = \frac{dS}{v} = \frac{\sqrt{1+y'^2}}{\sqrt{2gy}}\,dx$$

总时间 $T$ 是路径 $y(x)$ 的函数（泛函）：
$$T[y] = \int_{x_A}^{x_B} \frac{\sqrt{1+y'^2}}{\sqrt{2gy}}\,dx = \frac{1}{\sqrt{2g}}\int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx$$

这就是一个典型的变分问题：在满足端点条件 $y(x_A)=0, y(x_B)=y_B$ 的所有函数 $y(x)$ 中，寻找使 $T[y]$ 取最小值的那个。

### 3.1.3 核心问题

变分法的核心问题可以概括为：

> 从合适的函数集中选择某个函数 $y(x)$，使某个依赖于该函数整体的量（即泛函）达到极值。

数学上，最常见的泛函形式是：
$$Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx \to \text{extremum}$$

这里的"自变量"不再是微积分中的数 $x$，而是**函数 $y(x)$ 本身**。

---

## 3.2 从函数到泛函

### 3.2.1 泛函的定义

**泛函（Functional）**是从函数空间到实数域的映射。换句话说，泛函是"函数的函数"——它输入一个函数，输出一个实数。

**定义**：设 $D$ 是一个函数集合。如果对于任意函数 $f(x) \in D$，都有唯一一个实数 $Q$ 与之对应，则称 $Q$ 是 $f$ 的**泛函**，记作：
$$Q = Q[f(x)]$$

### 3.2.2 函数 vs 泛函

理解泛函的关键在于区分它与普通函数：

| 概念 | 函数 $y = f(x)$ | 泛函 $Q = Q[f]$ |
|------|----------------|-----------------|
| 输入 | 一个数 $x$ | 一个**函数** $f(x)$ |
| 输出 | 一个数 $y$ | 一个数 $Q$ |
| 依赖关系 | 依赖于 $x$ 的**取值** | 依赖于 $f$ 的**形状** |

**例 1**：$Q[f] = \int_0^1 f(x)\,dx$
- 当 $f(x) = x$ 时，$Q = \frac12$
- 当 $f(x) = \sqrt{x}$ 时，$Q = \frac23$
- 同一个输入 $x=0.5$ 在两个函数下的值完全不同，但 $Q$ 关心的是整个曲线的"形状"——即曲线下的面积

**例 2**：$N[f] = $ 函数 $f$ 在 $[a,b]$ 上的零点个数
- $N[\cos x]$ 在 $[0,2\pi]$ 上等于 2
- $N[\sin x]$ 在 $[0,2\pi]$ 上等于 3
- 这也是泛函——输入函数、输出数

### 3.2.3 线性泛函

若泛函 $Q$ 同时满足以下两个条件：

1. **齐次性**：$Q[cy] = c\,Q[y]$（$c$ 为任意常数）
2. **可加性**：$Q[y_1 + y_2] = Q[y_1] + Q[y_2]$

则称 $Q$ 是**线性泛函**。两个条件可以合并为：
$$Q[c_1y_1 + c_2y_2] = c_1Q[y_1] + c_2Q[y_2]$$

**例 3**：判断 $Q[y] = \int_a^b y(x)\,dx$ 是否为线性泛函

$$\begin{aligned}
Q[c_1y_1 + c_2y_2] &= \int_a^b [c_1y_1(x) + c_2y_2(x)]\,dx \\
&= c_1\int_a^b y_1(x)\,dx + c_2\int_a^b y_2(x)\,dx \\
&= c_1Q[y_1] + c_2Q[y_2]
\end{aligned}$$

因此 $Q[y] = \int_a^b y(x)dx$ **是**线性泛函。✅

**例 4**：判断 $Q[y] = \int_a^b y^2(x)\,dx$ 是否为线性泛函

$$\begin{aligned}
Q[c_1y_1 + c_2y_2] &= \int_a^b (c_1y_1 + c_2y_2)^2\,dx \\
&= \int_a^b [c_1^2y_1^2 + c_2^2y_2^2 + 2c_1c_2y_1y_2]\,dx
\end{aligned}$$

而：
$$c_1Q[y_1] + c_2Q[y_2] = c_1\int_a^b y_1^2\,dx + c_2\int_a^b y_2^2\,dx$$

两式不相等（多出了交叉项 $2c_1c_2y_1y_2$），因此 $Q[y] = \int_a^b y^2dx$ **不是**线性泛函。❌

> **判断线性泛函的陷阱**：必须同时检查齐次性和可加性，缺一不可。工程中常见的泛函大多是**非线性**的（因为涉及 $y^2, y'^2$ 等项），但这不影响我们用变分法求解。

---

## 3.3 变分的概念

### 3.3.1 自变函数的变分

在微积分中，我们研究自变量的微小变化 $\Delta x$ 对函数值的影响。在变分法中，我们研究自变函数（即曲线）的微小变化对泛函值的影响。

**定义**：两个函数之间的差称为**变分**（variation），记作 $\delta y$：
$$\delta y = y(x) - y_1(x)$$

类似地，导数的变分为：
$$\delta y' = y'(x) - y_1'(x)$$

> 变分 $\delta y$ 类似于微分 $dx$：$dx$ 是一个**微小量**，$\delta y$ 是一个**微小函数**。在端点固定的问题中，$\delta y$ 必须满足 $\delta y(a) = \delta y(b) = 0$，因为边界值已经给定，不能变化。

### 3.3.2 函数的微分（回顾）

为了理解泛函的变分，先回顾函数微分的两种定义方式：

**定义一（常规）**：将函数增量分解为线性主部和高阶小量。
$$\Delta y = A(x)\Delta x + \varphi(x,\Delta x)\Delta x^2$$
其中 $A(x)$ 与 $\Delta x$ 无关。当 $\Delta x \to 0$ 时 $\varphi(x,\Delta x) \to 0$。称 $dy = A(x)dx = y'(x)dx$ 为函数 $y$ 的**微分**。

**定义二（Lagrange 法）**：引入参数 $\varepsilon$。
$$\left.\frac{\partial}{\partial\varepsilon}y(x + \varepsilon\Delta x)\right|_{\varepsilon=0} = y'(x)\Delta x = dy$$

### 3.3.3 泛函变分的常规定义

仿照函数微分的定义一。令 $\delta y$ 为自变函数的变分，则泛函的增量为：
$$\Delta Q = Q[y + \delta y] - Q[y]$$

如果能将 $\Delta Q$ 分解为：
$$\Delta Q = T[y, \delta y] + \beta[y, \delta y]$$

其中：
- $T[y, \delta y]$ 对 $\delta y$ 是**线性泛函**（当 $y$ 固定时）
- $\beta[y, \delta y] / \max|\delta y| \to 0$ 当 $\max|\delta y| \to 0$（高阶小量）

则称 $T[y, \delta y]$ 为泛函 $Q$ 在 $y$ 处的**一阶变分**，记作 $\delta Q$。

**例 5**：求 $Q[y] = \int_a^b y^2(x)dx$ 的变分 $\delta Q$

$$\begin{aligned}
\Delta Q &= \int_a^b (y + \delta y)^2 dx - \int_a^b y^2 dx \\
&= \int_a^b [y^2 + 2y\delta y + (\delta y)^2] dx - \int_a^b y^2 dx \\
&= \int_a^b 2y\,\delta y\,dx + \int_a^b (\delta y)^2 dx
\end{aligned}$$

第一项 $\int_a^b 2y\,\delta y\,dx$ 对 $\delta y$ 是线性的 → 这就是 $\delta Q$。
第二项 $\int_a^b (\delta y)^2 dx$ 是 $\delta y$ 的高阶小量（因为 $(\delta y)^2$ 比 $\delta y$ 更快趋于零）。

因此：
$$\delta Q = \int_a^b 2y(x)\,\delta y(x)\,dx$$

### 3.3.4 泛函变分的 Lagrange 定义（计算利器）

**定义**：固定 $y(x)$ 和变分 $\delta y$，构造辅助函数：
$$\varphi(\alpha) = Q[y + \alpha\delta y]$$

则泛函变分为：
$$\boxed{\delta Q = \left.\frac{\partial\varphi}{\partial\alpha}\right|_{\alpha=0} = \left.\frac{\partial}{\partial\alpha}Q[y + \alpha\delta y]\right|_{\alpha=0}}$$

> **为什么这个定义好用？** 它将变分问题转化为普通微分问题——只需要对参数 $\alpha$ 求导，然后令 $\alpha=0$，不需要记忆复杂的分解公式。

**例 6**：用 Lagrange 定义重做例 5

构造 $\varphi(\alpha) = \int_a^b (y + \alpha\delta y)^2 dx$，则：
$$\varphi'(\alpha) = \int_a^b 2(y + \alpha\delta y)\delta y\,dx$$
$$\varphi'(0) = \int_a^b 2y\,\delta y\,dx = \delta Q$$

与常规定义的结果一致。✅

### 3.3.5 变分算子的运算法则

一个重要性质将在推导 Euler 方程时发挥关键作用：

**微分算子与变分算子可交换顺序**：
$$\boxed{\frac{d}{dx}(\delta y) = \delta(y')}$$

即：对变分求导数，等于对导数求变分。

**证明**：
$$\frac{d}{dx}(\delta y) = \frac{d}{dx}[y_1(x) - y(x)] = y_1'(x) - y'(x) = \delta(y')$$

这个性质使我们可以在分部积分时自由地将 $\delta y'$ 替换为 $(\delta y)'$。

---

## 3.4 泛函的极值

### 3.4.1 极值的必要条件

类比函数的极值条件（$dy = 0$），泛函取极值的**必要条件**是：
$$\boxed{\delta Q = 0}$$

即在极值函数处，泛函的一阶变分为零。

### 3.4.2 强极值与弱极值

根据所考虑的"邻域"的不同，泛函极值分为两类：

- **强极值**：在零阶邻域（函数值本身接近）内成立。即对任意 $|y(x)-y_0(x)| < \varepsilon$ 的曲线，$Q[y] \geq Q[y_0]$。
- **弱极值**：在一阶邻域（函数值和一阶导数都接近）内成立。即要求 $|y(x)-y_0(x)| < \varepsilon$ 且 $|y'(x)-y_0'(x)| < \varepsilon$。

> 对于大多数工程问题，我们只关心极值的**必要条件** $\delta Q = 0$。充分性通常由物理背景来保证——我们知道系统应该取最小势能，而不需要验证二阶变分的正定性。

---

## 3.5 Euler 方程的完整推导

这是全书最重要的一节。**考试必考**。

### 3.5.1 问题设定

求泛函：
$$Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx$$

在边界条件 $y(x_1) = y_1, y(x_2) = y_2$ 下的极值曲线。

### 3.5.2 推导过程

**Step 1：构造邻域函数**

设 $y(x)$ 是使 $Q$ 取极值的待求函数。考虑在 $y(x)$ 附近的一条"扰动"曲线：
$$y_1(x) = y(x) + \alpha\eta(x)$$

其中：
- $\eta(x)$ 是任意满足 $\eta(x_1) = \eta(x_2) = 0$ 的光滑函数（端点固定，所以扰动在端点为0）
- $\alpha$ 是一个小参数，用来控制扰动的大小

当 $\alpha = 0$ 时，$y_1 = y$，此时 $Q$ 取极值。

**Step 2：转化为 $\alpha$ 的函数**

将 $y_1$ 代入泛函：
$$\varphi(\alpha) = Q[y + \alpha\eta] = \int_{x_1}^{x_2} F(x, y + \alpha\eta, y' + \alpha\eta')\,dx$$

由于 $Q$ 在 $y$ 处取极值，$\varphi(\alpha)$ 在 $\alpha = 0$ 处取极值，因此根据 Fermat 定理：
$$\varphi'(0) = 0$$

**Step 3：在积分号内求导**

应用链式法则，在积分号内对 $\alpha$ 求导：
$$\varphi'(\alpha) = \int_{x_1}^{x_2} \left[\frac{\partial F}{\partial y}\eta + \frac{\partial F}{\partial y'}\eta'\right]dx$$

令 $\alpha = 0$：
$$\varphi'(0) = \int_{x_1}^{x_2} \left[F_y\,\eta + F_{y'}\,\eta'\right]dx = 0 \tag{3.1}$$

其中 $F_y = \partial F/\partial y$，$F_{y'} = \partial F/\partial y'$，均在 $\alpha=0$ 处取值。

**Step 4：分部积分**

对式 (3.1) 的第二项进行分部积分：
$$\int_{x_1}^{x_2} F_{y'}\,\eta'\,dx = \left[F_{y'}\,\eta\right]_{x_1}^{x_2} - \int_{x_1}^{x_2} \frac{d}{dx}(F_{y'})\,\eta\,dx$$

由于 $\eta(x_1) = \eta(x_2) = 0$（端点固定），边界项为零：
$$\int_{x_1}^{x_2} F_{y'}\,\eta'\,dx = -\int_{x_1}^{x_2} \frac{d}{dx}(F_{y'})\,\eta\,dx$$

代入式 (3.1)：
$$\int_{x_1}^{x_2} \left[F_y - \frac{d}{dx}(F_{y'})\right]\eta(x)\,dx = 0 \tag{3.2}$$

**Step 5：利用变分法预备定理**

> **变分法预备定理**（Fundamental Lemma of Calculus of Variations）：
> 若 $f(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 的光滑函数 $\eta(x)$，都有 $\int_a^b f(x)\eta(x)dx = 0$，则 $f(x) \equiv 0$ 在 $[a,b]$ 上。

**证明概要（反证法）**：假设存在 $x_0 \in (a,b)$ 使得 $f(x_0) > 0$。由连续性，存在 $\delta > 0$ 使在 $(x_0-\delta, x_0+\delta)$ 上 $f(x) > 0$。构造一个在该区间内为正、其余为零的光滑函数 $\eta(x)$，则 $\int f\eta dx > 0$，与假设矛盾。

**应用**：令 $f(x) = F_y - \frac{d}{dx}F_{y'}$，它在 $[x_1,x_2]$ 上连续（假设 $y$ 有二阶连续导数）。式 (3.2) 对任意满足端点为零的 $\eta(x)$ 成立。由预备定理：

$$\boxed{F_y - \frac{d}{dx}F_{y'} = 0}$$

这就是 **Euler-Lagrange 方程**，通常简称为 **Euler 方程**。

### 3.5.3 推导回顾

Euler 方程的推导是变分法的核心，需要透彻理解每一步：

1. **构造扰动** → 将变分问题 $\delta Q=0$ 转化为普通函数求导 $\varphi'(0)=0$
2. **求导 + 分部积分** → 将 $\eta'$ 项转移到 $\eta$ 上，以便提取公因子
3. **利用边界条件** → 分部积分产生的边界项因 $\eta$ 在端点为 0 而消失
4. **应用预备定理** → 从积分恒等式导出微分方程

---

## 3.6 Euler 方程的应用

### 3.6.1 解题标准流程

```
① 写出被积函数 F(x, y, y')
② 计算偏导 ∂F/∂y 和 ∂F/∂y'
③ 计算全导 d/dx(∂F/∂y')
④ 代入 Euler 方程 → 得到一个 ODE
⑤ 解 ODE → 通解
⑥ 代入边界条件 → 确定积分常数
```

### 3.6.2 例 7：最短路径问题

**问题**：求 $A(0,0)$ 到 $B(2,1)$ 的最短曲线。

**解**：两点间曲线的弧长为 $S[y] = \int_0^2 \sqrt{1+y'^2}\,dx$。

$F = \sqrt{1+y'^2}$，$F_y = 0$，$F_{y'} = \frac{y'}{\sqrt{1+y'^2}}$

Euler 方程：
$$0 - \frac{d}{dx}\left(\frac{y'}{\sqrt{1+y'^2}}\right) = 0 \quad\Rightarrow\quad \frac{y'}{\sqrt{1+y'^2}} = C$$

解出 $y'$：
$$y'^2 = C^2(1+y'^2) \quad\Rightarrow\quad y'^2 = \frac{C^2}{1-C^2}$$

即 $y' = m$（常数），$y = mx + b$

代入 $y(0)=0 \Rightarrow b=0$，$y(2)=1 \Rightarrow m = \frac12$。

**结论**：$y = \frac12 x$ —— 连接两点的直线。■

### 3.6.3 $F$ 不显含 $x$ 时的首次积分

当 $F = F(y, y')$（即 $F$ 不显含 $x$）时，有一个重要的**首次积分**：
$$\boxed{F - y'F_{y'} = C}$$

**推导**：
$$\begin{aligned}
\frac{d}{dx}(F - y'F_{y'}) &= F_y y' + F_{y'}y'' - y''F_{y'} - y'(F_{y'y}y' + F_{y'y'}y'') \\
&= y'[F_y - F_{y'y}y' - F_{y'y'}y''] \\
&= y'\left[F_y - \frac{d}{dx}F_{y'}\right] = 0
\end{aligned}$$

因为 $F_y - \frac{d}{dx}F_{y'} = 0$（Euler 方程），所以括号内为零，故 $F - y'F_{y'} =$ 常数。

> 这个首次积分在求解最速降线等问题时非常有用——它将二阶 ODE 降为一阶，大大简化了计算。

### 3.6.4 例 8：完整的 ODE 求解

**问题**：求 $Q[y] = \int_0^1 [(y')^2 + 4y^2 - 8xy]\,dx$ 在 $y(0)=1, y(1)=2$ 下的极值函数。

**解**：$F = y'^2 + 4y^2 - 8xy$

$$\frac{\partial F}{\partial y} = 8y - 8x,\quad \frac{\partial F}{\partial y'} = 2y'$$

Euler 方程：
$$(8y - 8x) - \frac{d}{dx}(2y') = 0 \quad\Rightarrow\quad y'' - 4y = -4x$$

这是一个非齐次二阶线性常系数 ODE。

**齐次解**：特征方程 $r^2 - 4 = 0$，$r = \pm 2$
$$y_h = C_1e^{2x} + C_2e^{-2x}$$

**特解**（待定系数法）：设 $y_p = Ax + B$，代入：
$$0 - 4(Ax + B) = -4x \quad\Rightarrow\quad A = 1,\; B = 0$$
$$y_p = x$$

**通解**：
$$y = y_h + y_p = C_1e^{2x} + C_2e^{-2x} + x$$

**代入边界条件**：
$$\begin{cases}
y(0) = C_1 + C_2 = 1 \\
y(1) = C_1e^2 + C_2e^{-2} + 1 = 2
\end{cases}$$

解得：
$$C_1 = \frac{1}{e^2+1},\quad C_2 = \frac{e^2}{e^2+1}$$

**最终解**：
$$\boxed{y(x) = \frac{e^{2x}}{e^2+1} + \frac{e^2 e^{-2x}}{e^2+1} + x}$$

**验证**：
- $y(0) = \frac{1}{e^2+1} + \frac{e^2}{e^2+1} + 0 = 1$ ✅
- $y(1) = \frac{e^2}{e^2+1} + \frac{e^2 e^{-2}}{e^2+1} + 1 = \frac{e^2+1}{e^2+1} + 1 = 2$ ✅

---

## 3.7 Euler 方程的推广

### 3.7.1 含二阶导数的泛函

对于 $Q[y] = \int F(x, y, y', y'')\,dx$：

$$\delta Q = \int (F_y\delta y + F_{y'}\delta y' + F_{y''}\delta y'')\,dx = 0$$

对 $F_{y''}\delta y''$ 分部积分**两次**：
$$\int F_{y''}\delta y'' dx = F_{y''}\delta y'\Big|_{x_1}^{x_2} - \left[\frac{d}{dx}F_{y''}\delta y\right]_{x_1}^{x_2} + \int \frac{d^2}{dx^2}F_{y''}\,\delta y\,dx$$

假设 $\delta y$ 和 $\delta y'$ 在端点为零，得：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0}$$

### 3.7.2 含 $n$ 阶导数的泛函（通用形式）

$$Q[y] = \int F(x, y, y', y'', \ldots, y^{(n)})\,dx$$

$$\boxed{\sum_{k=0}^n (-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right) = 0}$$

符号规律：$+F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} - \frac{d^3}{dx^3}F_{y'''} + \cdots + (-1)^n\frac{d^n}{dx^n}F_{y^{(n)}} = 0$

### 3.7.3 多个独立函数的泛函

$$Q[y_1, y_2, \ldots, y_n] = \int F(x, y_1, \ldots, y_n, y_1', \ldots, y_n')\,dx$$

对每个 $y_i$ 独立地成立 Euler 方程：
$$F_{y_i} - \frac{d}{dx}F_{y_i'} = 0,\quad i = 1, 2, \ldots, n$$

### 3.7.4 多元函数的泛函

$$Q[z(x,y)] = \iint_D F(x, y, z, p, q)\,dxdy,\quad p = \frac{\partial z}{\partial x},\; q = \frac{\partial z}{\partial y}$$

$$\boxed{F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0}$$

**例 9**：$Q[z] = \iint_D (z_x^2 + z_y^2)\,dxdy$ 给出 Laplace 方程！

$F = p^2 + q^2$，$F_z = 0$，$F_p = 2p$，$F_q = 2q$，代入得：
$$0 - \frac{\partial}{\partial x}(2z_x) - \frac{\partial}{\partial y}(2z_y) = 0 \quad\Rightarrow\quad \frac{\partial^2 z}{\partial x^2} + \frac{\partial^2 z}{\partial y^2} = 0$$

这就是著名的 **Laplace 方程**。变分法提供了推导 PDE 的另一个视角！

---

## 3.8 边界条件分类

在变分推导中，分部积分产生的边界项必须妥善处理。

### 3.8.1 本质边界条件（Essential BC）

**定义**：边界值预先给定，与泛函的变分无关。在变分中，$\delta y$ 在边界上为零。

- 例：固定端 $y(a) = y_a$，或简支端 $w(0) = 0$
- 来源：问题陈述中直接给定的边界值

### 3.8.2 自然边界条件（Natural BC）

**定义**：当边界值未给定时，分部积分产生的边界项必须单独为零，由此自动导出的条件。

以 $\delta Q = \int (F_y - \frac{d}{dx}F_{y'})\delta y\,dx + [F_{y'}\delta y]_{x_1}^{x_2} = 0$ 为例：

若 $y(x_1)$ 未指定（即 $\delta y(x_1) \neq 0$），则必须有：
$$\left.\frac{\partial F}{\partial y'}\right|_{x=x_1} = 0$$

这就是**自然边界条件**。

### 3.8.3 两类边界条件的对比

| | 本质边界条件 | 自然边界条件 |
|--|------------|------------|
| 又称 | Dirichlet BC, Geometric BC | Neumann BC, Static BC |
| 来源 | 问题给定 | 变分极值自动导出 |
| 变分中 | $\delta y = 0$ | $\partial F/\partial y' = 0$ |
| 梁的例 | $w(0)=0$（固定端位移） | $EI w''(l)=0$（自由端弯矩） |

---

## 3.9 条件极值与 Lagrange 乘子法

实际中常遇到附加约束的泛函极值问题。例如：在曲线长度固定的条件下求面积最大，或者在满足某等周条件 $\int\varphi\,dx = \alpha$ 下求泛函极值。

**方法**（与普通函数条件极值的 Lagrange 乘子法完全类似）：

构造新泛函：
$$Q^*[y] = \int_a^b (F + \lambda\varphi)\,dx - \lambda\alpha = \int_a^b F^*\,dx - \lambda\alpha$$

其中 $\lambda$ 为 Lagrange 乘子（常数），$F^* = F + \lambda\varphi$。

对 $Q^*$ 应用 Euler 方程：
$$F^*_y - \frac{d}{dx}F^*_{y'} = 0$$

通解包含 2 个积分常数 + $\lambda$，由 2 个边界条件 + 1 个等周条件确定。

---

## 3.10 最速降线问题的完整求解

作为变分法的经典范例，我们完整求解最速降线问题。

### 3.10.1 泛函的简化

已得总时间泛函（略去常数 $1/\sqrt{2g}$）：
$$Q[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx,\quad F(y,y') = \sqrt{\frac{1+y'^2}{y}}$$

注意 $F$ 不显含 $x$！

### 3.10.2 应用首次积分

由 $F - y'F_{y'} = C$，计算 $F_{y'} = \frac{y'}{\sqrt{y(1+y'^2)}}$：
$$\sqrt{\frac{1+y'^2}{y}} - \frac{y'^2}{\sqrt{y(1+y'^2)}} = C$$

化简：
$$\frac{1}{\sqrt{y(1+y'^2)}} = C \quad\Rightarrow\quad y(1+y'^2) = D \quad(D = 1/C^2)$$

### 3.10.3 参数化求解

令 $y' = \tan\theta$，代入：
$$y = \frac{D}{1+\tan^2\theta} = D\cos^2\theta = \frac{D}{2}(1+\cos 2\theta)$$

$$dx = \frac{dy}{y'} = \frac{-D\sin 2\theta\,d\theta}{\tan\theta} = -D(1+\cos 2\theta)\,d\theta$$

积分得参数方程：
$$\begin{cases}
x = -\frac{D}{2}(2\theta + \sin 2\theta) + E \\[4pt]
y = \frac{D}{2}(1 + \cos 2\theta)
\end{cases}$$

### 3.10.4 化为旋轮线标准形式

令 $2\theta = \pi - \phi$，则：
$$\boxed{\begin{cases}
x = r(\phi - \sin\phi) \\[4pt]
y = r(1 - \cos\phi)
\end{cases}},\quad r = \frac{D}{2}$$

这就是**旋轮线（cycloid）**的参数方程——一个半径为 $r$ 的圆沿直线滚动时，圆周上某点的轨迹。

---

## 3.11 变分法记号与运算法则

变分法中有一些方便的计算规则，熟练掌握可大幅简化推导。

**规则 1**：变分算子与微分算子可交换
$$\frac{d}{dx}(\delta y) = \delta(y')$$

**规则 2**：变分算子与积分算子可交换
$$\delta\int_a^b F\,dx = \int_a^b \delta F\,dx$$

**规则 3**：复合函数的变分（与全微分形式相同）
$$\delta F(x,y,y') = \frac{\partial F}{\partial y}\delta y + \frac{\partial F}{\partial y'}\delta y'$$

---

## 3.12 经典变分法与 FEM 的关系

### Ritz 法
Ritz 法将变分问题近似求解：设 $u \approx \sum c_i\varphi_i(x)$，代入泛函，由 $\partial\Pi/\partial c_i = 0$ 得到线性方程组。

**局限**：$\varphi_i$ 必须满足全部边界条件，复杂几何下难以构造。

### FEM 的突破
FEM 将求解域划分为单元，在每个单元内独立建立简单的形函数：
- 形函数不必满足单元间的协调条件
- 程序标准化、自动化
- 是"在单元层面上应用 Ritz/Galerkin 法"

```
经典变分法 (Euler-Lagrange) → 精确但求解困难
    ↓
Ritz 法 (整体 trial function) → 边界条件难满足
    ↓
Courant 分片近似 (1943) → 克服边界困难
    ↓
FEM (单元级 trial function, 1960) → 程序化、通用化
```

---

## 检查你的理解

1. 什么是泛函？它与普通函数有什么本质区别？请举一个不是积分形式的泛函例子。
2. 线性泛函必须满足哪两个条件？$Q[y] = \int_a^b (y^2 + y'^2)dx$ 是否是线性泛函？为什么？
3. 写出 Euler 方程并解释其推导过程中分部积分和预备定理各自的作用。
4. 对于 $Q[y] = \int_a^b F(x,y,y',y'')dx$，它的 Euler 方程是什么？
5. 什么条件下可以使用首次积分 $F - y'F_{y'} = C$？
6. 本质边界条件和自然边界条件有什么区别？在梁问题中各对应于什么物理边界？

---

> **对应作业**：[HW2 Q1（最短路径）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q2（三阶导数 Euler 方程）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q3（Lagrange 乘子法）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q4（泛函极值函数）](../04-Homework-Solutions/2026w/HW2-Problem.md)
> **往年参考**：[past/HW2/homework 2](../04-Homework-Solutions/past/HW2/homework%202.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW2/Ans%20to%20HM2_LIU%20Sai_handed%20in.md)
