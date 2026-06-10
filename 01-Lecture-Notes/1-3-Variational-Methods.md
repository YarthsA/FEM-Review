# 变分法基础（Variational Calculus）

> **对应课件**：[`6 FEM_Element construction.pdf`](../06-References/pdfs-originals/6%20FEM_Element%20construction.pdf) 第1章 §IV · [原文MD](../06-References/../06-References/../md_output/6%20FEM_Element%20construction.md) §4.1-4.7
> **PDF 章节定位**：Chapter 1 → Theory of elasticity → IV. Basics of Variational Principles (Section 4.1–4.7)
> **相关作业**：[HW2 Q1-Q4](../04-Homework-Solutions/2026w/HW2-Problem.md) · [HW3 Q3（弹性地基梁变分推导）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：微积分（微分定义、分部积分、ODE 求解）、线性代数

---

## §4.1 引言

### 4.1.1 为什么需要近似方法？

弹性力学微分形式的方程很难获得精确解，原因在于：

■ 复杂的控制方程（偏微分方程组，15个方程耦合）
■ 复杂的边界条件（不规则几何形状、混合边界条件）

因此近似方法是必要的。**有限元分析（FEA）** 就是其中非常重要的一种。近似解可以非常接近精确解——实际上，在前面推导控制方程时，我们通常已经做了一些简化假设，这本身也是一种近似。

> 换言之，这个世界上不存在绝对精确的解。

### 4.1.2 最速降线问题（Brachistochrone Problem）

1696年，Johann Bernoulli（1667–1748）向全欧洲的数学家提出了一个挑战：

> 设在垂直平面内有任意两点 A 和 B（A 高于 B）。一质点在重力作用下从 A 沿某曲线无摩擦地滑到 B。问沿哪条曲线所需时间最短？

泛函表述为：
$$Q[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx$$

这个问题之所以著名，不仅在于求最大值或最小值，而在于要求出一个**未知函数（曲线）** 来满足所有条件。这个创新而有趣的问题迅速吸引了当时最著名的数学家。该期刊在1697年5月发表了 L'Hôpital（1661–1704）、Jacob Bernoulli（1654–1705）、Leibniz（1646–1716）和 Newton（1642–1727）的解答。他们都用不同的方法得到了相同的答案——最速降线是**旋轮线（cycloid）**。但仅凭这些解答还不能建立变分法的框架。

### 4.1.3 Euler 方程与变分法的建立

**Euler（1707–1783）** 在1736年解决了最速降线问题，并得到了解决这类问题的一般方法。他指出：若以下积分达到极值，
$$Q[y] = \int_{x_1}^{x_2} F[x, y(x), y'(x)]\,dx$$

则函数 $y(x)$ 必须满足：
$$\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) = 0$$

这就是著名的 **Euler 方程**。

**Lagrange（1736–1813）** 在1755年用数学分析的方法改进和简化了 Euler 的研究，变分法由此作为数学的一个新分支正式建立。

### 4.1.4 经典变分法的问题

经典变分法的主要问题可以概括为：在适当的函数集合内选取函数，使得某个积分达到极值，这可以通过求解相应的 Euler 方程来解决。

**但事实并非如此简单！**

几乎所有自然定律都可以用变分原理的形式来表达。变分法可以实现数学上的**统一**——数学物理定解问题通常可以转化为变分问题。变分法是解数学物理问题的近似方法，其基本思想就是把这样的问题转化为变分问题。

经典变分法促进了 Ritz 法等近似方法的发展，但 Ritz 法的试探函数选取仍然十分困难，而且系数矩阵的计算也非常复杂。直到计算机的发展，FEM 才迅速诞生。

---

## §4.2 泛函（Functional）

### 4.2.1 从例子理解泛函

变分法的研究对象是**泛函**，它是函数概念的推广。

**例**：设 $y=f(x)$ 是在闭区间上连续且非负的函数。则以下积分的几何意义是曲线下的面积：
$$Q[f] = \int_0^1 f(x)\,dx$$

当 $f(x)=x$ 时，$Q[f]=\int_0^1 x\,dx=\frac12$
当 $f(x)=\sqrt{x}$ 时，$Q[f]=\int_0^1\sqrt{x}\,dx=\frac23$

可见对于任意一个连续非负函数 $f(x)$，都有一个确定的 $Q$ 值与之对应。

### 4.2.2 函数 vs 泛函

**函数**：对于两个变量 $x$ 和 $y$，$x$ 属于数集 $D$，$y$ 属于数集 $R$。若对于 $D$ 中的每个 $x$，都有唯一的 $y$ 与之对应，则 $y$ 是 $x$ 的函数。

在之前的例子中，给定函数 $f(x)$ 后，$Q[f(x)]$ 被唯一确定。这意味着 $Q[f(x)]$ 可以看作是 $f(x)$ 的函数（记作 $Q[f]$）。这里 $f$ 可以看作一种自变量，$Q$ 就是它的函数。**泛函实际上是函数概念的推广**。

### 4.2.3 泛函的定义

考虑一个函数集 $D$。若对于 $D$ 中的任意函数 $f(x)$，都有唯一的一个 $Q$ 值与之对应，那么变量 $Q$ 可以称为函数 $f$ 的泛函，记作：
$$Q = Q[f(x)] \quad \text{或} \quad Q = Q[f]$$

> 注意：这只是一个粗略的定义，但足够用于后续讨论。

**用通俗的话说**：泛函就是"函数的函数"（不是复合函数那种）。函数通常依赖于自变量的**取值**，而泛函依赖于自变函数的**形状**。

### 4.2.4 更多例子

记 $N[f(x)]$ 为函数 $f(x)$ 在区间 $[a,b]$ 上零点的个数。则对任意实函数 $f$ 在 $[a,b]$ 上，都有一个确定的 $N$ 值。

例如取 $[a,b]=[0,2\pi]$：
$$N[\cos x] = 2,\quad N[\sin x] = 3,\quad N[x^2-1] = 1$$

根据前述定义，$N$ 是 $f$ 的泛函。

---

## §4.3 变分（Variation）

### 4.3.1 相关概念的推广

在深入讨论之前，需要将函数的相关概念推广到泛函：
- 函数极值 → **泛函极值**
- 函数连续性 → **泛函连续性**
- 函数微分 → **函数变分**

### 4.3.2 自变函数的变分

函数 $y(x)$ 自变量的增量为两个值的差：
$$\Delta x = x - x_1$$

类似地，泛函 $Q[y(x)]$ 的自变函数的增量为两个函数值的差：
$$\Delta y = y(x) - y_1(x)$$

当这个增量足够小时，称为**变分**，记作 $\delta y(x)$ 或 $\delta y$：
$$\delta y = y(x) - y_1(x),\quad \delta y' = y'(x) - y_1'(x)$$

自变函数通常需要满足特定的边界条件（如 $y(a)=y_a, y(b)=y_b$），所以自变函数的变分必须满足齐次边界条件：
$$\delta y\big|_{x=a} = 0,\quad \delta y\big|_{x=b} = 0$$

### 4.3.3 函数的接近度

所谓两条曲线"接近"有不同的含义：

**零阶接近度**：仅函数值相近，即 $\max|y(x)-y_1(x)|$ 很小，但导数可能相差很大。

**一阶接近度**：函数值和导数值都相近，即 $|y(x)-y_1(x)|$ 和 $|y'(x)-y_1'(x)|$ 都很小。

**K 阶接近度**：直到 $k$ 阶导数都很接近，即 $\delta y, \delta y', \ldots, \delta y^{(k)}$ 都很小。

显然 K 越高，曲线之间越接近。

### 4.3.4 连续泛函

对于函数：若对任意小的 $\varepsilon>0$，存在 $\delta>0$，使 $|x-x_1|<\delta$ 时总有 $|y(x)-y(x_1)|<\varepsilon$，则函数在 $x_1$ 处连续。

泛函的定义类似：若对任意小的 $\varepsilon>0$，存在 $\delta>0$，使当 $|y(x)-y_1(x)|<\delta, |y'(x)-y_1'(x)|<\delta, \ldots, |y^{(k)}(x)-y_1^{(k)}(x)|<\delta$ 时，总有 $|Q[y(x)]-Q[y_1(x)]|<\varepsilon$，则称泛函 $Q[y(x)]$ 是 K 阶连续的。

**例**：对于泛函 $Q[y(x)]=\int_a^b y(x)dx$，若 $y(x)\in C[a,b]$，则 $Q[y(x)]$ 是连续泛函。

**证明**：对任意 $y(x)$，若 $\max_{a\leq x\leq b}|y(x)-y_0(x)|<\frac{\varepsilon}{b-a}$，则：
$$|Q[y(x)]-Q[y_0(x)]| = \left|\int_a^b[y(x)-y_0(x)]dx\right| \leq \int_a^b\left|[y(x)-y_0(x)]\right|dx < \frac{\varepsilon}{b-a}\int_a^b dx = \varepsilon$$
得证。

### 4.3.5 线性泛函

若泛函 $Q[y(x)]$ 满足：
1. $Q[cy(x)] = cQ[y(x)]$（c 为任意常数）
2. $Q[y_1(x)+y_2(x)] = Q[y_1(x)]+Q[y_2(x)]$

则称 $Q$ 为线性泛函。两条件可合并为：
$$Q[c_1y_1(x)+c_2y_2(x)] = c_1Q[y_1(x)]+c_2Q[y_2(x)]$$

**例**：$Q[y(x)]=\int_a^b y(x)dx$ 是线性泛函 ✅
**例**：$Q[y(x)]=\int_a^b y^2(x)dx$ **不是**线性泛函 ❌

**验证**：
$$Q[c_1y_1+c_2y_2] = \int_a^b (c_1y_1+c_2y_2)^2dx = \int_a^b(c_1^2y_1^2 + c_2^2y_2^2 + 2c_1c_2y_1y_2)dx$$
$$c_1Q[y_1] + c_2Q[y_2] = c_1\int_a^b y_1^2dx + c_2\int_a^b y_2^2dx$$
两式因交叉项 $2c_1c_2y_1y_2$ 的存在而不等。

### 4.3.6 函数的微分（两种定义）

**定义一（常规）**：若函数增量可展开为 $\Delta y = A(x)\Delta x + \varphi(x,\Delta x)\Delta x^2$，其中 $A(x)$ 与 $\Delta x$ 无关，$\Delta x\to0$ 时 $\varphi\to0$。则微分为：
$$dy = A(x)dx = y'(x)dx$$

即：函数的微分是增量中的**线性主部**。

**定义二（Lagrange）**：
$$\left.\frac{\partial}{\partial\varepsilon}y(x+\varepsilon\Delta x)\right|_{\varepsilon=0} = y'(x)\Delta x = dy(x)$$

即：在 $\varepsilon=0$ 处对 $\varepsilon$ 求导。这个定义与 Lagange 给出的变分定义类似。

### 4.3.7 泛函的变分

**定义一（常规）**：仿照函数微分的定义。

泛函 $Q[y(x)]$ 的增量：
$$\Delta Q = Q[y(x)+\delta y] - Q[y(x)]$$

若 $\Delta Q$ 可表为：
$$\Delta Q = T[y(x),\delta y] + \beta[y(x),\delta y]$$

其中 $T$ 是 $\delta y$ 的线性泛函（$y$ 给定时），且 $\beta$ 满足 $\frac{\beta}{\max|\delta y|}\to0$（高阶小量），则 $T$ 称为泛函的**一阶变分**，记作 $\delta Q$。

**例**：$Q[y]=\int_a^b y^2(x)dx$ 的变分
$$\Delta Q = \int_a^b[(y+\delta y)^2 - y^2]dx = \int_a^b 2y\delta y\,dx + \int_a^b(\delta y)^2dx$$

第一项 $\int_a^b 2y\delta y\,dx$ 对 $\delta y$ 是线性的 → 这就是 $\delta Q$。第二项 $(\delta y)^2$ 比 $\delta y$ 更快趋于零，是高阶小量。

**定义二（Lagrange 法）**：
$$\varphi(\alpha) = Q[y(x)+\alpha\delta y]$$
则：
$$\delta Q = \left.\frac{\partial\varphi}{\partial\alpha}\right|_{\alpha=0} = \left.\frac{\partial}{\partial\alpha}Q[y+\alpha\delta y]\right|_{\alpha=0}$$

这个定义在计算中非常方便——对参数 $\alpha$ 求导再令 $\alpha=0$ 即可，不需要手动分解。

### 4.3.8 泛函极值

**函数极值**：若在 $x=x_0$ 附近，$dy=y(x)-y(x_0)\leq0(\geq0)$，则函数在 $x_0$ 处取极大（小）值。必要条件：$dy=0$。

**泛函极值**的定义类似。引入曲线 $y_0(x)$ 的 $\varepsilon$ 邻域概念——满足 $|y(x)-y_0(x)|\leq\varepsilon$ 的所有曲线的集合。

**强极值**：在零阶 $\varepsilon$ 邻域内成立，对任意 $|y(x)-y_0(x)|\leq\varepsilon$ 的曲线都有 $\Delta Q \leq 0$。

**弱极值**：在一阶 $\varepsilon$ 邻域内成立，要求 $|y(x)-y_0(x)|\leq\varepsilon$ **且** $|y'(x)-y_0'(x)|\leq\varepsilon$。

泛函极值的**必要条件**（类比函数极值的 $dy=0$）：
$$\delta Q = 0$$

即泛函的一阶变分为零。

---

## §4.4 最简单的 Euler 方程

### 4.4.1 变分法预备定理

若函数 $f(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 且 $|\eta(x)|\leq\varepsilon$ 的具有连续导数的非零函数 $\eta(x)$，都有：
$$\int_a^b f(x)\eta(x)dx = 0$$
则 $f(x)$ 在 $[a,b]$ 上恒等于零。

**描述性证明**（反证法）：若存在 $x_0\in(a,b)$ 使 $f(x_0)>0$，由连续性存在 $\delta>0$ 使当 $|x-x_0|<\delta$ 时 $f(x)>0$。构造一个在该区间内为正、其余为零的光滑函数 $\varphi(x)$，令 $\eta(x)=A\varphi(x)$（$A$ 为足够小的正数使 $|\eta(x)|\leq\varepsilon$）。则 $\int_a^b f(x)\eta(x)dx = \int_{x_0-\delta}^{x_0+\delta}f(x)\eta(x)dx > 0$，与假设矛盾。

### 4.4.2 简单 Euler 方程

求泛函 $Q[y] = \int_a^b F[x, y(x), y'(x)]dx$ 的极值曲线。

**定理**：若 $F(x,y,y')$ 是三个变量的连续函数，且 $F$ 及其一阶、二阶偏导在定义域内连续。若泛函 $Q[y]$ 在具有二阶连续导数的曲线 $y(x)$ 上取极值（满足 $y(a)=y_0$，$y(b)=y_1$，且位于平面有界域 B 内），则 $y(x)$ 必须满足：
$$F_y - \frac{d}{dx}F_{y'} = 0$$

这就是 **Euler 方程**。

### 4.4.3 推导证明

**Step 1**：设 $y(x)$ 是使 $Q$ 取极值的函数。任选一个函数 $\eta(x)$，满足 $\eta(a)=\eta(b)=0$，$\eta(x)\in C^1[a,b]$。则当 $|\alpha|$ 足够小时，曲线 $y_1(x)=y(x)+\alpha\eta(x)$ 在 $y(x)$ 的一阶 $\varepsilon$ 邻域内。

**Step 2**：由于 $y(x)$ 和 $\eta(x)$ 都给定，$Q[y(x)+\alpha\eta(x)]$ 是 $\alpha$ 的函数，记作 $\varphi(\alpha)$。

**Step 3**：$Q$ 在 $y(x)$ 取极值 → $\varphi(\alpha)$ 在 $\alpha=0$ 取极值 → $\varphi'(0)=0$。

$$\begin{aligned}
\varphi'(0) &= \left.\frac{d}{d\alpha}\int_a^b F[x, y+\alpha\eta, y'+\alpha\eta']dx\right|_{\alpha=0} \\
&= \int_a^b [F_y\eta + F_{y'}\eta']dx = 0
\end{aligned}$$

**Step 4**：分部积分：
$$\int_a^b F_{y'}\eta' dx = F_{y'}\eta\Big|_a^b - \int_a^b \eta\frac{d}{dx}F_{y'}dx = -\int_a^b \eta\frac{d}{dx}F_{y'}dx$$

代入得：
$$\int_a^b \left(F_y - \frac{d}{dx}F_{y'}\right)\eta(x)dx = 0$$

**Step 5**：由变分法预备定理，被积函数必须恒为零：
$$F_y - \frac{d}{dx}F_{y'} = 0$$

> 注意：Euler 方程只是泛函取极值的**必要条件**而非充分条件。但对于大多数工程问题，物理背景提供了直观判断，不需要验证充分性。

### 4.4.4 Euler 方程示例

**例**：求泛函 $Q[y(x)]=\int_0^1(x^2y + y'^2)dx$ 满足 $y(0)=0, y(1)=1/3$ 的驻值曲线。

**解**：$F = x^2y + y'^2$，$F_y = x^2$，$F_{y'} = 2y'$

Euler 方程：$x^2 - \frac{d}{dx}(2y') = 0$ → $2y'' = x^2$ → $y'' = \frac12x^2$

积分得：$y = \frac{1}{24}x^4 + C_1x + C_2$

代入边界条件：$y(0)=0 \Rightarrow C_2=0$，$y(1)=\frac{1}{24}+C_1=\frac13 \Rightarrow C_1=\frac{7}{24}$

**解**：$y(x) = \frac{1}{24}x^4 + \frac{7}{24}x$

### 4.4.5 最速降线问题的 Euler 方程

对于最速降线问题，$F(y,y') = \sqrt{\frac{1+y'^2}{y}}$

直接计算 $F_y$ 和 $F_{y'}$ 并代入 Euler 方程非常复杂。但因为 $F$ 不显含 $x$，可以利用首次积分。

**首次积分**：当 $F=F(y,y')$ 时：
$$F - y'F_{y'} = C$$

**证明**：
$$\frac{d}{dx}(F - y'F_{y'}) = F_yy' + F_{y'}y'' - y''F_{y'} - y'(F_{y'y}y' + F_{y'y'}y'') = y'(F_y - \frac{d}{dx}F_{y'}) = 0$$

代入 $F$：
$$\sqrt{\frac{1+y'^2}{y}} - \frac{y'^2}{\sqrt{y(1+y'^2)}} = C \quad\Rightarrow\quad \frac{1}{\sqrt{y(1+y'^2)}} = C$$

即 $y(1+y'^2) = D$（$D=1/C^2$）。

令 $y'=\tan\theta$：
$$y = \frac{D}{1+\tan^2\theta} = D\cos^2\theta = \frac{D}{2}(1+\cos2\theta)$$

$$dx = \frac{dy}{y'} = \frac{-D\sin2\theta d\theta}{\tan\theta} = -D(1+\cos2\theta)d\theta$$

积分得参数方程：
$$\begin{cases}
x = -\frac{D}{2}(2\theta + \sin2\theta) + E \\
y = \frac{D}{2}(1+\cos2\theta)
\end{cases}$$

令 $2\theta = \pi - \phi$，得旋轮线标准形式：
$$\boxed{\begin{cases}
x = r(\phi - \sin\phi) \\
y = r(1 - \cos\phi)
\end{cases}},\quad r = \frac{D}{2}$$

### 4.4.6 变分记号与运算法则

**基本法则**：$\delta(y') = (\delta y)'$，即变分算子与微分算子可交换。

利用这一法则，Euler 方程的推导可以写为：
$$\delta Q = \int_a^b [F_y\delta y + F_{y'}\delta y']dx = \int_a^b[F_y\delta y + F_{y'}d(\delta y)]dx$$
$$= \left.\delta y F_{y'}\right|_a^b + \int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx = 0$$

由于 $\delta y(a)=\delta y(b)=0$（边界条件），得：
$$\int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx = 0 \quad\Rightarrow\quad F_y - \frac{d}{dx}F_{y'} = 0$$

### 4.4.7 本质边界条件与自然边界条件

前面讨论的边界条件（$y(a)=y_a$，$y(b)=y_b$）是预先定义和固定的，因此边界值与泛函的一阶变分无关。这类边界条件称为**本质边界条件**。

现在考虑可变边界问题。变分结果为：
$$\delta Q = \int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx + \left.F_{y'}\delta y\right|_a^b$$

第一项→Euler 方程必须成立。对于第二项，当 $x=a$ 和 $x=b$ 时，必须有：
$$\left.\frac{\partial F}{\partial y'}\right|_{x=a} = 0,\quad \left.\frac{\partial F}{\partial y'}\right|_{x=b} = 0$$

否则可以找到 $\delta y$ 使 $\delta Q \neq 0$。这类边界条件称为**自然边界条件**。

### 4.4.8 泛函的二阶变分

类似函数的极值条件，一阶变分 $\delta Q = 0$ 只是极值的必要条件（驻值条件）。

二阶变分：
$$\delta^2 Q = \frac12\int_a^b\left(\frac{\partial^2F}{\partial y^2}\delta y^2 + 2\frac{\partial^2F}{\partial y\partial y'}\delta y\delta y' + \frac{\partial^2F}{\partial y'^2}\delta y'^2\right)dx$$

充分条件：
- $\delta Q = 0$ 且 $\delta^2 Q > 0$ → 极小值
- $\delta Q = 0$ 且 $\delta^2 Q < 0$ → 极大值

实际应用中通常只考虑一阶变分，充分性由物理背景保证。

### 4.4.9 条件极值与 Lagrange 乘子法

前面的泛函极值问题要求自变函数光滑且满足给定边界条件，无其他附加条件——称为**无条件极值问题**。

具有附加约束条件的泛函极值问题称为**条件极值问题**。类似于函数的条件极值，可以用 Lagrange 乘子法转化为等价的无条件极值问题。

**例**：求泛函 $Q[y(x)] = \int_a^b F(x,y,y')dx$ 在等周条件 $\int_a^b \varphi(x,y,y')dx = \alpha$ 下的极值。

构造新泛函：
$$Q^*[y] = \int_a^b F(x,y,y')dx + \lambda\left[\int_a^b\varphi(x,y,y')dx - \alpha\right] = \int_a^b F^*(x,y,y',\lambda)dx - \lambda\alpha$$

其中 $F^* = F + \lambda\varphi$，$\lambda$ 是 Lagrange 乘子。

对应的 Euler 方程为：
$$F_y^* - \frac{d}{dx}F_{y'}^* = 0$$

通解含 2 个积分常数和 $\lambda$，由 2 个边界条件 + 等周条件确定。

---

## §4.5 更一般的 Euler 方程

### 4.5.1 含高阶导数的泛函

对于 $Q = \int_a^b F(x, y, y', y'')dx$，变分为：
$$\delta Q = \int_a^b (F_y\delta y + F_{y'}\delta y' + F_{y''}\delta y'')dx = 0$$

对 $F_{y''}\delta y''$ 项分部积分两次：
$$\int_a^b F_{y''}\delta y'' dx = \left.\delta y'F_{y''}\right|_a^b - \left.\delta y\frac{d}{dx}F_{y''}\right|_a^b + \int_a^b \delta y\frac{d^2}{dx^2}F_{y''}dx$$

假设 $\delta y'(a) = \delta y'(b) = 0$ 和 $\delta y(a) = \delta y(b) = 0$，得：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0}$$

**推广到 $n$ 阶导数**（通用形式）：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} - \cdots + (-1)^n\frac{d^n}{dx^n}F_{y^{(n)}} = 0}$$

### 4.5.2 含多个独立函数的泛函

对于 $Q[y_1(x), \ldots, y_n(x)] = \int_a^b F(x, y_1, \ldots, y_n, y_1', \ldots, y_n')dx$：

变分后得：
$$\delta Q = \int_a^b \left[\delta y_1\left(F_{y_1} - \frac{d}{dx}F_{y_1'}\right) + \cdots + \delta y_n\left(F_{y_n} - \frac{d}{dx}F_{y_n'}\right)\right]dx = 0$$

由于 $\delta y_i$ 可独立任意选取，令 $\delta y_2 = \cdots = \delta y_n = 0$，得：
$$F_{y_i} - \frac{d}{dx}F_{y_i'} = 0,\quad i = 1, 2, \ldots, n$$

即每个 $y_i$ 独立满足 Euler 方程。

### 4.5.3 含多元函数的泛函

对于 $Q[z(x,y)] = \iint_D F[x, y, z(x,y), p, q]\,dxdy$，其中 $p=\partial z/\partial x$，$q=\partial z/\partial y$。

**Euler 方程**：
$$\boxed{F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0}$$

**推导要点**：利用 Green 公式将面积分转化为边界积分，利用 $\delta z$ 在边界上为零消去边界项。

**例**：$Q[z] = \iint_D (z_x^2 + z_y^2)dxdy$，边界条件 $z|_{\partial D} = \varphi(x,y)$

$F = p^2+q^2$，$F_p=2p$，$F_q=2q$，代入 Euler 方程：
$$\frac{\partial^2 z}{\partial x^2} + \frac{\partial^2 z}{\partial y^2} = 0$$

这就是著名的 **Laplace 方程**！

---

## §4.6 变分法在力学中的应用

### 4.6.1 Fermat 原理与广义相对论

- **Fermat 原理**：在均匀介质中，光线沿所需时间最少的光路传播。
- **Einstein 广义相对论**：光线在四维 Riemann 空间中沿所需时间最短的路线传播。

变分法的应用不仅限于经典物理和工程学科，而是延伸到量子场论、控制论、信息论等现代高科技领域。

### 4.6.2 Hamilton 原理

当 $t=t_1$ 和 $t=t_2$ 时，质点分别位于 A 和 B。则其真实运动路径应使以下**作用量积分**取极值：
$$S = \int_{t_1}^{t_2} (T - U)dt$$

其中 $T$ 是动能，$U$ 是势能，$T-U$ 称为 **Lagrange 函数**。

Hamilton 原理（1834）是力学的基本原理，可推出 Newton 三定律、能量守恒、动量守恒和角动量守恒。

### 4.6.3 Lagrange 方程

对于具有 $K$ 个自由度的系统，广义坐标为 $q_1, \ldots, q_K$，广义速度为 $\dot{q}_1, \ldots, \dot{q}_K$。

Lagrange 函数 $L = T - U = L(t, q, \dot{q})$

由 Hamilton 原理和 Euler 方程：
$$\frac{\partial L}{\partial q_i} - \frac{d}{dt}\frac{\partial L}{\partial \dot{q}_i} = 0,\quad i=1,2,\ldots,K$$

这就是著名的 **Lagrange 方程**——在广义坐标下描述运动基本规律的方程。

**Lagrange 方程的优势**：只需写出广义坐标下的动能和势能表达式，代入即可，不需要画受力图。

### 4.6.4 单摆

**例**：质量为 $m$ 的质点悬挂在长度为 $l$ 的轻绳末端，在重力场中运动。

取广义坐标为 $\theta$：
- 动能：$T = \frac12mv^2 = \frac12m(l\dot{\theta})^2$
- 势能（$\theta=0$ 设为零）：$U = mgl(1-\cos\theta)$

Lagrange 函数：$L = \frac12ml^2\dot{\theta}^2 - mgl(1-\cos\theta)$

代入 Lagrange 方程：
$$\frac{\partial L}{\partial\theta} - \frac{d}{dt}\frac{\partial L}{\partial\dot{\theta}} = -mgl\sin\theta - ml^2\ddot{\theta} = 0$$

$$\ddot{\theta} + \frac{g}{l}\sin\theta = 0$$

小角度时 $\sin\theta \approx \theta$，得简谐振动 $\ddot{\theta} + \frac{g}{l}\theta = 0$。

### 4.6.5 两端固定弦的自由横向振动

建立坐标系，设 $u(x,t)$ 为点 $x$ 在 $t$ 时刻的横向位移。

考虑弦上一微段变形后长度 $dl = \sqrt{1+u_x^2}dx$，伸长为 $(\sqrt{1+u_x^2}-1)dx \approx \frac12u_x^2dx$。

弹性势能：$U = \int_0^l \frac{K}{2}u_x^2 dx$
动能：$T = \int_0^l \frac12\rho(x)u_t^2 dx$

作用量泛函：
$$S = \int_{t_1}^{t_2}\int_0^l\left[\frac12\rho u_t^2 - \frac12K u_x^2\right]dx\,dt$$

Euler 方程给出：
$$\frac{\partial}{\partial x}(Ku_x) - \frac{\partial}{\partial t}(\rho u_t) = 0 \quad\Rightarrow\quad u_{tt} = \frac{K}{\rho(x)}u_{xx}$$

这就是经典的**弦振动方程**。

---

## §4.7 结论

从以上例子可以总结变分法的主要步骤：

1. **建立泛函**及其边界条件（基于物理问题）
2. **利用变分法预备定理**，通过泛函变分得到 Euler 方程
3. **求解 Euler 方程**——这是一个微分方程问题

注意：变分法和 Euler 方程描述的是同一个物理问题。从变分法和 Euler 方程出发获得近似解具有相同的效果。

Euler 方程往往很难求解。但有些问题已知微分方程但很难求解，如果能转化为相应的泛函极值问题，就可以用近似方法（如有限元分析）方便地求解。

> 然而，**并非每个微分方程都能找到合适的泛函**。这种情况下，我们必须使用其他方法获得近似解，如最小二乘法和 Galerkin 法——这些将在后续章节讨论。

---

## 检查你的理解

1. 什么是泛函？它与普通函数有什么本质区别？请举一个不是积分形式的泛函例子。
2. 泛函变分的Lagrange定义是什么？为什么这个定义在计算中更方便？
3. Euler 方程推导的五个步骤中，分部积分和预备定理各自起什么作用？
4. 为什么 $F$ 不显含 $x$ 时有首次积分 $F - y'F_{y'} = C$？
5. 本质边界条件和自然边界条件有什么区别？各举一例。

---

> **对应作业**：[HW2 Q1（最短路径）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q2（三阶导数Euler方程）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q3（Lagrange乘子法）](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q4（泛函极值函数）](../04-Homework-Solutions/2026w/HW2-Problem.md)
> **往年参考**：[past/HW2/homework 2](../04-Homework-Solutions/past/HW2/homework%202.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW2/Ans%20to%20HM2_LIU%20Sai_handed%20in.md)
