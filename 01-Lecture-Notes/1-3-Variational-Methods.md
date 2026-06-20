# 变分法基础（Variational Calculus）

> **对应课件**：[`Chapter 3 Variation theory and applications-1.pdf`](../06-References/pdfs-originals/Chapter%203%20Variation%20theory%20and%20applications-1.pdf) · [原文MD](../../md_output/Chapter%203%20Variation%20theory%20and%20applications-1.md)
> **相关作业**：[HW2 Q1-Q4](../04-Homework-Solutions/2026w/HW2-Problem.md) · [HW3 Q3（弹性地基梁变分推导）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：微积分（微分定义、分部积分、ODE 求解）、线性代数

> **📋 考试范围覆盖**
>
> | 本讲义章节 | 考试大纲考点 |
> |-----------|-------------|
> | §3.2 泛函 | [Var. Princ.] Functional |
> | §3.3 变分 | [Var. Princ.] Variation of functional; [Var. Princ.] Functional extremum |
> | §3.4 最简单的 Euler 方程 | [Var. Princ.] Euler equation |
> | §3.4.7 本质与自然边界条件 | [Var. Princ.] Essential and natural boundary conditions |
> | §3.4.9 条件极值与 Lagrange 乘子 | [Var. Princ.] Conditional extremum of functional |
> | §3.5 更一般的 Euler 方程 | [Var. Princ.] Extended forms of Euler equation |
> | §3.6 变分法在力学中的应用 | [Var. Princ.] Applications of variation method in mechanics |
> | §3.8.3 虚功原理 | [Var. Princ.] Principle of virtual work; Permissible/virtual displacement and stress |
> | §3.8.4 虚位移/虚应力原理 | [Var. Princ.] Principle of virtual displacement / virtual stress |
> | §3.8.5 最小势能原理 | [Var. Princ.] Principle of minimum potential energy |
> | §3.8.6 弹性力学 Euler 方程 | [Var. Princ.] Euler equations in elastic mechanics |
> | §3.8.7 直接法与间接法 | [Var. Princ.] Direct and indirect methods of variation problems |

---

## 3.1 引言

### 3.1.1 为什么需要近似方法？

弹性力学微分形式的方程很难获得精确解，原因在于：

■ 复杂的控制方程（偏微分方程组，15个方程耦合）
■ 复杂的边界条件（不规则几何形状、混合边界条件）

因此近似方法是必要的。**有限元分析（FEA）** 就是其中非常重要的一种。近似解可以非常接近精确解——实际上，在前面推导控制方程时，我们通常已经做了一些简化假设，这本身也是一种近似。

> 换言之，这个世界上不存在绝对精确的解。

> 💡 理解关键：整个变分法这一章的本质目的只有一个——为后续的 FEM 近似解法（Ritz 法、Galerkin 法）提供理论基础。变分法把"解微分方程"转化成"找泛函极值"，后者在工程上更容易用数值方法逼近。

### 3.1.2 最速降线问题（Brachistochrone Problem）

1696年，Johann Bernoulli（1667–1748）向全欧洲的数学家提出了一个挑战：

> 设在垂直平面内有任意两点 A 和 B（A 高于 B）。一质点在重力作用下从 A 沿某曲线无摩擦地滑到 B。问沿哪条曲线所需时间最短？

泛函表述为：
$$Q[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx$$

这个问题之所以著名，不仅在于求最大值或最小值，而在于要求出一个**未知函数（曲线）** 来满足所有条件。这个创新而有趣的问题迅速吸引了当时最著名的数学家。该期刊在1697年5月发表了 L'Hôpital（1661–1704）、Jacob Bernoulli（1654–1705）、Leibniz（1646–1716）和 Newton（1642–1727）的解答。他们都用不同的方法得到了相同的答案——最速降线是**旋轮线（cycloid）**。但仅凭这些解答还不能建立变分法的框架。


### 3.1.3 Euler 方程与变分法的建立

**Euler（1707–1783）** 在1736年解决了最速降线问题，并得到了解决这类问题的一般方法。他指出：若以下积分达到极值，
$$Q[y] = \int_{x_1}^{x_2} F[x, y(x), y'(x)]\,dx$$

则函数 $y(x)$ 必须满足：
$$\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) = 0$$

这就是著名的 **Euler 方程**。

**Lagrange（1736–1813）** 在1755年用数学分析的方法改进和简化了 Euler 的研究，变分法由此作为数学的一个新分支正式建立。

> 💡 理解关键：Euler 方程可以类比为"泛函版本的 $f'(x)=0$"。普通函数极值的必要条件是导数为零 $\frac{df}{dx}=0$；泛函极值的必要条件就是 Euler 方程成立。因为泛函的"自变量"是一整条曲线，所以"导数"变成一个微分方程而不是一个数。

### 3.1.4 经典变分法的问题

经典变分法的主要问题可以概括为：在适当的函数集合内选取函数，使得某个积分达到极值，这可以通过求解相应的 Euler 方程来解决。

**但事实并非如此简单！**

> ⚠️ 重难点：看似变分法把问题"简化"成了 Euler 方程，但 Euler 方程本身往往也是难解的偏微分方程！所以如果到这里就停了，等于什么都没解决。真正的突破在于后面的"直接法"——绕过 Euler 方程，直接对泛函求极值。

几乎所有自然定律都可以用变分原理的形式来表达。变分法可以实现数学上的**统一**——数学物理定解问题通常可以转化为变分问题。变分法是解数学物理问题的近似方法，其基本思想就是把这样的问题转化为变分问题。

经典变分法促进了 Ritz 法等近似方法的发展，但 Ritz 法的试探函数选取仍然十分困难，而且系数矩阵的计算也非常复杂。直到计算机的发展，FEM 才迅速诞生。

> 🔗 跨章连接：这里埋下的"直接法"种子，将在第 4 章开花——Ritz 法和 Galerkin 法就是两种最核心的直接法，它们分别从泛函和微分方程出发用试探函数逼近真实解，而 FEM 本质上就是分片定义的 Ritz/Galerkin 法。

---

## 3.2 泛函（Functional）

### 3.2.1 从例子理解泛函

变分法的研究对象是**泛函**，它是函数概念的推广。

**例**：设 $y=f(x)$ 是在闭区间上连续且非负的函数。则以下积分的几何意义是曲线下的面积：
$$Q[f] = \int_0^1 f(x)\,dx$$

当 $f(x)=x$ 时，$Q[f]=\int_0^1 x\,dx=\frac12$
当 $f(x)=\sqrt{x}$ 时，$Q[f]=\int_0^1\sqrt{x}\,dx=\frac23$

可见对于任意一个连续非负函数 $f(x)$，都有一个确定的 $Q$ 值与之对应。

> ⚠️ 重难点：这是泛函概念的第一道坎。函数是"数 → 数"的映射，泛函是"函数 → 数"的映射。不要把"泛函"和"复合函数"搞混——复合函数 $g(f(x))$ 的输出仍然是 $x$ 的函数，而泛函 $Q[f]$ 的输出是一个确定的数值。

### 3.2.2 函数 vs 泛函

**函数**：对于两个变量 $x$ 和 $y$，$x$ 属于数集 $D$，$y$ 属于数集 $R$。若对于 $D$ 中的每个 $x$，都有唯一的 $y$ 与之对应，则 $y$ 是 $x$ 的函数。

在之前的例子中，给定函数 $f(x)$ 后，$Q[f(x)]$ 被唯一确定。这意味着 $Q[f(x)]$ 可以看作是 $f(x)$ 的函数（记作 $Q[f]$）。这里 $f$ 可以看作一种自变量，$Q$ 就是它的函数。**泛函实际上是函数概念的推广**。

### 3.2.3 泛函的定义

考虑一个函数集 $D$。若对于 $D$ 中的任意函数 $f(x)$，都有唯一的一个 $Q$ 值与之对应，那么变量 $Q$ 可以称为函数 $f$ 的泛函，记作：
$$Q = Q[f(x)] \quad \text{或} \quad Q = Q[f]$$

> 注意：这只是一个粗略的定义，但足够用于后续讨论。

**用通俗的话说**：泛函就是"函数的函数"（不是复合函数那种）。函数通常依赖于自变量的**取值**，而泛函依赖于自变函数的**形状**。

> 💡 理解关键：记法上有个细节——方括号 $Q[f]$ 暗示"自变量"是一整个函数，圆括号 $f(x)$ 暗示自变量是某一点的数值。变分法中用 $Q[y(x)]$ 或 $Q[y]$ 的方括号记法强调这一点。考试中写错括号无所谓，但自己心里要清楚区别。

### 3.2.4 更多例子

记 $N[f(x)]$ 为函数 $f(x)$ 在区间 $[a,b]$ 上零点的个数。则对任意实函数 $f$ 在 $[a,b]$ 上，都有一个确定的 $N$ 值。

例如取 $[a,b]=[0,2\pi]$：
$$N[\cos x] = 2,\quad N[\sin x] = 3,\quad N[x^2-1] = 1$$

根据前述定义，$N$ 是 $f$ 的泛函。


---

## 3.3 变分（Variation）

### 3.3.1 相关概念的推广

在深入讨论之前，需要将函数的相关概念推广到泛函：
- 函数极值 → **泛函极值**
- 函数连续性 → **泛函连续性**
- 函数微分 → **函数变分**

> 💡 理解关键：整个 §3.3 做的就是一件事情——把微积分里你熟悉的每个概念（导数、微分、连续、极值）都对应地推广到泛函上。思维框架是"函数 → 泛函"的类比，一路对过去就行。如果某个定义看不懂，回头想：它的函数版本是什么？

### 3.3.2 自变函数的变分

函数 $y(x)$ 自变量的增量为两个值的差：
$$\Delta x = x - x_1$$

类似地，泛函 $Q[y(x)]$ 的自变函数的增量为两个函数值的差：
$$\Delta y = y(x) - y_1(x)$$

当这个增量足够小时，称为**变分**，记作 $\delta y(x)$ 或 $\delta y$：
$$\delta y = y(x) - y_1(x),\quad \delta y' = y'(x) - y_1'(x)$$

自变函数通常需要满足特定的边界条件（如 $y(a)=y_a, y(b)=y_b$），所以自变函数的变分必须满足齐次边界条件：
$$\delta y\big|_{x=a} = 0,\quad \delta y\big|_{x=b} = 0$$

> ❌ 易错点：$\delta y(x)$ 和 $dy(x)$ 是两个完全不同的概念！$dy = y'dx$ 是函数在某一固定曲线上的微分（沿曲线走一小步），而 $\delta y$ 是两条不同曲线在同一 $x$ 处的差值（换一条曲线）。视觉上 $\delta y$ 是纵向比较（不同曲线间），$dy$ 是横向比较（同一曲线上）。


### 3.3.3 函数的接近度

所谓两条曲线"接近"有不同的含义：

**零阶接近度**：仅函数值相近，即 $\max|y(x)-y_1(x)|$ 很小，但导数可能相差很大。

**一阶接近度**：函数值和导数值都相近，即 $|y(x)-y_1(x)|$ 和 $|y'(x)-y_1'(x)|$ 都很小。

**K 阶接近度**：直到 $k$ 阶导数都很接近，即 $\delta y, \delta y', \ldots, \delta y^{(k)}$ 都很小。

显然 K 越高，曲线之间越接近。

> 💡 理解关键：区分"强极值"和"弱极值"的精髓就在这里。零阶接近度的邻域比一阶大得多（允许导数乱跳），所以在零阶邻域中的极值更难达到→叫"强"极值。简单记：邻域越宽 → 极值越强。

### 3.3.4 连续泛函

对于函数：若对任意小的 $\varepsilon>0$，存在 $\delta>0$，使 $|x-x_1|<\delta$ 时总有 $|y(x)-y(x_1)|<\varepsilon$，则函数在 $x_1$ 处连续。

泛函的定义类似：若对任意小的 $\varepsilon>0$，存在 $\delta>0$，使当 $|y(x)-y_1(x)|<\delta, |y'(x)-y_1'(x)|<\delta, \ldots, |y^{(k)}(x)-y_1^{(k)}(x)|<\delta$ 时，总有 $|Q[y(x)]-Q[y_1(x)]|<\varepsilon$，则称泛函 $Q[y(x)]$ 是 K 阶连续的。

**例**：对于泛函 $Q[y(x)]=\int_a^b y(x)dx$，若 $y(x)\in C[a,b]$，则 $Q[y(x)]$ 是连续泛函。记号 $C[a,b]$ 表示区间 $[a,b]$ 上全体连续函数构成的集合；$C^1[a,b]$ 表示区间 $[a,b]$ 上全体一阶导数连续的函数构成的集合。

**证明**：对任意 $y(x)$，若 $\max_{a\leq x\leq b}|y(x)-y_0(x)|<\frac{\varepsilon}{b-a}$，则：
$$|Q[y(x)]-Q[y_0(x)]| = \left|\int_a^b[y(x)-y_0(x)]dx\right| \leq \int_a^b\left|[y(x)-y_0(x)]\right|dx < \frac{\varepsilon}{b-a}\int_a^b dx = \varepsilon$$
得证。

### 3.3.5 线性泛函

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

> 💡 理解关键：线性泛函定义跟线性函数完全一致，只是把"数"换成"函数"。考试时判断一个泛函是否线性，最快速的方法：① 看被积函数里 $y$ 是几次方——平方以上一定不线性；② 变量代换法：令测试函数 $\eta_1, \eta_2$，看 $Q[c_1\eta_1+c_2\eta_2] = c_1Q[\eta_1]+c_2Q[\eta_2]$ 是否成立。

### 3.3.6 函数的微分（两种定义）

**定义一（常规）**：若函数增量可展开为 $\Delta y = A(x)\Delta x + \varphi(x,\Delta x)\Delta x^2$，其中 $A(x)$ 与 $\Delta x$ 无关，$\Delta x\to0$ 时 $\varphi\to0$。则微分为：
$$dy = A(x)dx = y'(x)dx$$

即：函数的微分是增量中的**线性主部**。

**定义二（Lagrange）**：
$$\left.\frac{\partial}{\partial\varepsilon}y(x+\varepsilon\Delta x)\right|_{\varepsilon=0} = y'(x)\Delta x = dy(x)$$

即：在 $\varepsilon=0$ 处对 $\varepsilon$ 求导。这个定义与 Lagange 给出的变分定义类似。

> ❌ 易错点：这两种定义对函数微分是等价的，但对泛函变分来说 Lagrange 定义更方便——因为不需要手动分解 $\Delta Q$ 的线性部分和高阶部分。考试计算变分时推荐用 Lagrange 法：令 $\varphi(\alpha)=Q[y+\alpha\delta y]$，然后 $\delta Q = \varphi'(0)$。

### 3.3.7 泛函的变分

**定义一（常规）**：仿照函数微分的定义。

泛函 $Q[y(x)]$ 的增量：
$$\Delta Q = Q[y(x)+\delta y] - Q[y(x)]$$

若 $\Delta Q$ 可表为：
$$\Delta Q = T[y(x),\delta y] + \beta[y(x),\delta y]$$

其中 $T$ 是 $\delta y$ 的线性泛函（$y$ 给定时），且 $\beta$ 满足 $\frac{\beta}{\max|\delta y|}\to0$（高阶小量），则 $T$ 称为泛函的**一阶变分**，记作 $\delta Q$。

**例**：$Q[y]=\int_a^b y^2(x)dx$ 的变分
$$\Delta Q = \int_a^b[(y+\delta y)^2 - y^2]dx = \int_a^b 2y\delta y\,dx + \int_a^b(\delta y)^2dx$$

第一项 $\int_a^b 2y\delta y\,dx$ 对 $\delta y$ 是线性的 → 这就是 $\delta Q$。第二项 $(\delta y)^2$ 比 $\delta y$ 更快趋于零，是高阶小量。

**例**：$Q[y]=\int_a^b y y'^2 dx$，用两种方法分别计算。

**常规法**：

$$\Delta Q = \int_a^b (y+\delta y)(y'+\delta y')^2 dx - \int_a^b y(y')^2 dx$$

展开 $(y'+\delta y')^2 = (y')^2 + 2y'\delta y' + (\delta y')^2$：

$$= \int_a^b \left[y(y')^2 + 2yy'\delta y' + y(\delta y')^2 + \delta y(y')^2 + 2\delta y\cdot y'\delta y' + \delta y(\delta y')^2\right]dx - \int_a^b y(y')^2 dx$$

$$= \int_a^b \left[2yy'\delta y' + y(\delta y')^2 + (\delta y)(y')^2 + 2(\delta y)(y')(\delta y') + \delta y(\delta y')^2\right]dx$$

其中对 $\delta y$ 或 $\delta y'$ 线性的项：$\int_a^b (2yy'\delta y' + (y')^2\delta y)\,dx$

其余项含 $(\delta y)^2$、$(\delta y')^2$、$\delta y\cdot\delta y'$ → 高阶小量，丢掉。

$$\boxed{\delta Q = \int_a^b (2yy'\delta y' + (y')^2\delta y)\,dx}$$

**Lagrange 法**：

$$\varphi(\alpha) = Q[y+\alpha\delta y] = \int_a^b (y+\alpha\delta y)(y'+\alpha\delta y')^2 dx$$

$$\frac{d\varphi}{d\alpha} = \int_a^b \left[\delta y\cdot(y'+\alpha\delta y')^2 + (y+\alpha\delta y)\cdot 2(y'+\alpha\delta y')\delta y'\right]dx$$

令 $\alpha=0$：

$$\boxed{\delta Q = \int_a^b \left[(y')^2\delta y + 2yy'\delta y'\right]dx}$$

> 两种方法结果相同，但 Lagrange 法不需要手动展开和判断哪些是高阶小量，直接对 $\alpha$ 求导即可。

**定义二（Lagrange 法）**：
$$\varphi(\alpha) = Q[y(x)+\alpha\delta y]$$
则：
$$\boxed{\delta Q = \left.\frac{\partial\varphi}{\partial\alpha}\right|_{\alpha=0} = \left.\frac{\partial}{\partial\alpha}Q[y+\alpha\delta y]\right|_{\alpha=0}}$$

这个定义在计算中非常方便——对参数 $\alpha$ 求导再令 $\alpha=0$ 即可，不需要手动分解。

> ⚠️ **重难点**：Lagrange 法是整个变分法的**核心计算工具**。后续 Euler 方程推导（§3.4）的五步全部基于此定义。考试中计算变分时**必须用 Lagrange 法**，不要用常规法——常规法需要手动分解线性部分和高阶部分，容易出错且效率低。

**Lagrange 法操作步骤**（以 $Q[y]=\int_a^b F(x,y,y')dx$ 为例）：

**Step 1**：构造扰动函数 $y+\alpha\delta y$，代入泛函：
$$\varphi(\alpha) = Q[y+\alpha\delta y] = \int_a^b F(x,\; y+\alpha\delta y,\; y'+\alpha\delta y')dx$$

**Step 2**：对 $\alpha$ 求导（在积分号下求导，利用链式法则）：
$$\frac{d\varphi}{d\alpha} = \int_a^b \left[\frac{\partial F}{\partial y}\delta y + \frac{\partial F}{\partial y'}\delta y'\right]dx$$

**Step 3**：令 $\alpha=0$：
$$\delta Q = \left.\frac{d\varphi}{d\alpha}\right|_{\alpha=0} = \int_a^b \left[F_y\delta y + F_{y'}\delta y'\right]dx$$

> 💡 **理解关键**：Lagrange 法的精髓是**引入参数 $\alpha$ 将泛函变分转化为普通函数求导**。$\alpha=0$ 对应原始函数 $y$，$\alpha\neq 0$ 对应扰动后的函数。对 $\alpha$ 求导再令 $\alpha=0$，等价于提取 $\delta y$ 的线性部分——这就是一阶变分的定义。

**计算示例**：求 $Q[y]=\int_0^1 (y'^2 + xy)dx$ 的变分。

Step 1：$\varphi(\alpha) = \int_0^1 [(y'+\alpha\delta y')^2 + x(y+\alpha\delta y)]dx$

Step 2：$\dfrac{d\varphi}{d\alpha} = \int_0^1 [2(y'+\alpha\delta y')\delta y' + x\,\delta y]dx$

Step 3：令 $\alpha=0$：$\delta Q = \int_0^1 [2y'\delta y' + x\,\delta y]dx$


### 3.3.8 泛函极值

**函数极值**：若在 $x=x_0$ 附近，$dy=y(x)-y(x_0)\leq0(\geq0)$，则函数在 $x_0$ 处取极大（小）值。必要条件：$dy=0$。

**泛函极值**的定义类似。引入曲线 $y_0(x)$ 的 $\varepsilon$ 邻域概念——满足 $|y(x)-y_0(x)|\leq\varepsilon$ 的所有曲线的集合。

**强极值**：在零阶 $\varepsilon$ 邻域内成立，对任意 $|y(x)-y_0(x)|\leq\varepsilon$ 的曲线都有 $\Delta Q \leq 0$。

**弱极值**：在一阶 $\varepsilon$ 邻域内成立，要求 $|y(x)-y_0(x)|\leq\varepsilon$ **且** $|y'(x)-y_0'(x)|\leq\varepsilon$。

泛函极值的**必要条件**（类比函数极值的 $dy=0$）：
$$\boxed{\delta Q = 0}$$

即泛函的一阶变分为零。

> ⚠️ 重难点：强/弱极值的区别经常出现在概念题中。记住：强极值对曲线的限制少（零阶邻域，导数可以乱跳），所以"对付"的曲线更多→更难满足极值条件→更强。弱极值要求一阶邻域（导数也接近），能被弱极值排除的曲线较少→条件较弱。工程上几乎只用弱极值。

---

## 3.4 最简单的 Euler 方程

> ⚠️ 重难点：§3.4 是本章最核心的技术内容。Euler 方程的推导只有五个步骤，但每一步都有明确的作用。考试推导题大概率考察完整推导或至少关键步骤。建议把这五步背下来：① 构造扰动曲线 $y+\alpha\eta$，② 写成 $\alpha$ 的函数 $\varphi(\alpha)$，③ 极值条件 $\varphi'(0)=0$，④ 分部积分消去 $\eta'$，⑤ 变分法预备定理令被积函数为零。

### 3.4.1 变分法预备定理

**变分法预备定理**：若函数 $f(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 且 $|\eta(x)|\leq\varepsilon$ 的具有连续导数的非零函数 $\eta(x)$，都有：
$$\boxed{\int_a^b f(x)\eta(x)dx = 0}$$
则 $f(x)$ 在 $[a,b]$ 上恒等于零。

**描述性证明**（反证法）：若存在 $x_0\in(a,b)$ 使 $f(x_0)>0$，由连续性存在 $\delta>0$ 使当 $|x-x_0|<\delta$ 时 $f(x)>0$。构造一个在该区间内为正、其余为零的光滑函数 $\varphi(x)$，令 $\eta(x)=A\varphi(x)$（$A$ 为足够小的正数使 $|\eta(x)|\leq\varepsilon$）。则 $\int_a^b f(x)\eta(x)dx = \int_{x_0-\delta}^{x_0+\delta}f(x)\eta(x)dx > 0$，与假设矛盾。

> 💡 理解关键：预备定理的核心思想——如果 $f(x)$ 和任何你选择的光滑小扰动 $\eta(x)$ 的乘积积分都为零，那 $f(x)$ 自身必定处处为零。这相当于说"一个人如果跟全世界所有人握手都没感觉，那他一定是死人"。考试中预备定理通常简化写成"由 $\eta$ 的任意性，被积函数恒为零"即可。

> 🔗 跨章连接：变分法预备定理是 Euler 方程推导的最后一步，也是后续 Galerkin 法的基础思想。Galerkin 法中 $\int_V R \cdot w_i\,dV = 0$ 的合理性同样来自"权函数的任意性"——只不过那里权函数不是任意选取的，而是有限个基函数，所以是近似而非精确。

### 3.4.2 简单 Euler 方程

求泛函 $Q[y] = \int_a^b F[x, y(x), y'(x)]dx$ 的极值曲线。

> 💡 **记号说明**：本节证明中用 $\eta(x)$ 表示扰动函数，而 §3.3 中用 $\delta y$ 表示变分——它们是**同一个东西**，都代表"任意的微小扰动"。证明时用独立符号 $\eta$ 是为了更清晰地叙述边界条件 $\eta(a)=\eta(b)=0$，以及最后引用变分法预备定理时强调"$\eta$ 的任意性"。计算具体泛函的变分时则直接用 $\delta y$，更简洁。


**定理**：若 $F(x,y,y')$ 是三个变量的连续函数，且 $F$ 及其一阶、二阶偏导在定义域内连续。若泛函 $Q[y]$ 在具有二阶连续导数的曲线 $y(x)$ 上取极值（满足 $y(a)=y_0$，$y(b)=y_1$，且位于平面有界域 B 内），则 $y(x)$ 必须满足：
$$F_y - \frac{d}{dx}F_{y'} = 0$$

这就是 **Euler 方程**。

### 3.4.3 推导证明

**Step 1**：设 $y(x)$ 是使 $Q$ 取极值的函数。任选一个函数 $\eta(x)$，满足 $\eta(a)=\eta(b)=0$，$\eta(x)\in C^1[a,b]$。则当 $|\alpha|$ 足够小时，曲线 $y_1(x)=y(x)+\alpha\eta(x)$ 在 $y(x)$ 的一阶 $\varepsilon$ 邻域内。

**Step 2**：由于 $y(x)$ 和 $\eta(x)$ 都给定，$Q[y(x)+\alpha\eta(x)]$ 是 $\alpha$ 的函数，记作 $\varphi(\alpha)$。

**Step 3**：$Q$ 在 $y(x)$ 取极值 → $\varphi(\alpha)$ 在 $\alpha=0$ 取极值 → $\varphi'(0)=0$。

$$\begin{aligned}
\varphi'(0) &= \left.\frac{d}{d\alpha}\int_a^b F[x, y+\alpha\eta, y'+\alpha\eta']dx\right|_{\alpha=0} \\
&= \int_a^b [F_y\eta + F_{y'}\eta']dx = 0
\end{aligned}$$

> 💡 理解关键：Step 3 这一步是 Lagrange 法的直接应用。把 $\alpha$ 当参数，把泛函极值转化成普通函数的极值问题 $\varphi'(0)=0$。这本质上就是"降维打击"——用普通微积分处理泛函问题。

**Step 4**：分部积分：
$$\int_a^b F_{y'}\eta' dx = F_{y'}\eta\Big|_a^b - \int_a^b \eta\frac{d}{dx}F_{y'}dx = -\int_a^b \eta\frac{d}{dx}F_{y'}dx$$

> ❌ 易错点：Step 4 分部积分的边界项 $F_{y'}\eta|_a^b$ 之所以等于零，不是因为有自然边界条件，而是因为 $\eta(a)=\eta(b)=0$（这是 Step 1 中人为选择的）。考试中很多同学会漏写这句话——直接写"边界项为零"却不说明为什么，会被扣分。

**Step 5**：代入得：
$$\int_a^b \left(F_y - \frac{d}{dx}F_{y'}\right)\eta(x)dx = 0$$

由变分法预备定理，被积函数必须恒为零：
$$\boxed{F_y - \frac{d}{dx}F_{y'} = 0}$$

> 注意：Euler 方程只是泛函取极值的**必要条件**而非充分条件。但对于大多数工程问题，物理背景提供了直观判断，不需要验证充分性。

### 3.4.4 Euler 方程示例

**例**：求泛函 $Q[y(x)]=\int_0^1(x^2y + y'^2)dx$ 满足 $y(0)=0, y(1)=1/3$ 的驻值曲线。

**解**：$F = x^2y + y'^2$，$F_y = x^2$，$F_{y'} = 2y'$

Euler 方程：$x^2 - \frac{d}{dx}(2y') = 0$ → $2y'' = x^2$ → $y'' = \frac12x^2$

积分得：$y = \frac{1}{24}x^4 + C_1x + C_2$

代入边界条件：$y(0)=0 \Rightarrow C_2=0$，$y(1)=\frac{1}{24}+C_1=\frac13 \Rightarrow C_1=\frac{7}{24}$

**解**：$y(x) = \frac{1}{24}x^4 + \frac{7}{24}x$

> ❌ 易错点：求出 Euler 方程（一个二阶 ODE）后，解这个 ODE 会引入两个积分常数 $C_1, C_2$。这需要两个边界条件来确定。如果题目只给了一个边界条件，通常意味着另一个边界是"自由"的，需要用自然边界条件 $F_{y'}|_{\text{边界}} = 0$ 来补充。

### 3.4.5 最速降线问题的 Euler 方程

对于最速降线问题，$F(y,y') = \sqrt{\frac{1+y'^2}{y}}$

直接计算 $F_y$ 和 $F_{y'}$ 并代入 Euler 方程非常复杂。但因为 $F$ 不显含 $x$，可以利用首次积分。


**首次积分**：当 $F=F(y,y')$ 时：
$$\boxed{F - y'F_{y'} = C}$$

**证明**：
$$\frac{d}{dx}(F - y'F_{y'}) = F_yy' + F_{y'}y'' - y''F_{y'} - y'(F_{y'y}y' + F_{y'y'}y'') = y'(F_y - \frac{d}{dx}F_{y'}) = 0$$

> 💡 理解关键：首次积分本质上是对 Euler 方程的一次积分，将二阶 ODE 降为一阶。这跟经典力学中能量守恒（对运动方程的一次积分）是同一个数学结构。

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

### 3.4.6 变分记号与运算法则

**基本法则**：$\delta(y') = (\delta y)'$，即变分算子与微分算子可交换。

利用这一法则，Euler 方程的推导可以写为：
$$\delta Q = \int_a^b [F_y\delta y + F_{y'}\delta y']dx = \int_a^b[F_y\delta y + F_{y'}d(\delta y)]dx$$
$$= \left.\delta y F_{y'}\right|_a^b + \int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx = 0$$

由于 $\delta y(a)=\delta y(b)=0$（边界条件），得：
$$\int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx = 0 \quad\Rightarrow\quad F_y - \frac{d}{dx}F_{y'} = 0$$

> 💡 理解关键：$\delta(y') = (\delta y)'$ 这个交换律非常重要——它保证了变分操作和微分操作可以互换顺序。物理直觉：先变分再微分和先微分再变分，结果是同一个东西。证明很简单：$\delta(y') = (y+\alpha\eta)' - y' = \alpha\eta' = (\alpha\eta)' = (\delta y)'$。

### 3.4.7 本质边界条件与自然边界条件

前面讨论的边界条件（$y(a)=y_a$，$y(b)=y_b$）是预先定义和固定的，因此边界值与泛函的一阶变分无关。这类边界条件称为**本质边界条件**。

现在考虑可变边界问题。变分结果为：
$$\boxed{\delta Q = \int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx + \left.F_{y'}\delta y\right|_a^b}$$

第一项→Euler 方程必须成立。对于第二项，当 $x=a$ 和 $x=b$ 时，必须有：
$$\boxed{\left.\frac{\partial F}{\partial y'}\right|_{x=a} = 0,\quad \left.\frac{\partial F}{\partial y'}\right|_{x=b} = 0}$$

否则可以找到 $\delta y$ 使 $\delta Q \neq 0$。这类边界条件称为**自然边界条件**。

> ⚠️ 重难点：这是本章最容易混淆的概念之一！记法：**本质边界条件 = 变分前就定好的（如位移边界），变分时 $\delta y|_{\text{边界}} = 0$；自然边界条件 = 变分过程中自然"掉出来"的（如力边界 $F_{y'}|_{\text{边界}} = 0$）**。在力学中，位移边界 → 本质边界，力边界 → 自然边界。考试概念题大概率考这个区分。


### 3.4.8 泛函的二阶变分

类似函数的极值条件，一阶变分 $\delta Q = 0$ 只是极值的必要条件（驻值条件）。

二阶变分：
$$\delta^2 Q = \frac12\int_a^b\left(\frac{\partial^2F}{\partial y^2}\delta y^2 + 2\frac{\partial^2F}{\partial y\partial y'}\delta y\delta y' + \frac{\partial^2F}{\partial y'^2}\delta y'^2\right)dx$$

充分条件：
- $\delta Q = 0$ 且 $\delta^2 Q > 0$ → 极小值
- $\delta Q = 0$ 且 $\delta^2 Q < 0$ → 极大值

实际应用中通常只考虑一阶变分，充分性由物理背景保证。

> 🔗 跨章连接：第 4 章 Ritz 法中用到的 $\delta\Pi = 0$（总势能驻值）就是一阶变分。而 $\delta^2\Pi \geq 0$ 保证了 FEM 中总势能泛函在真实解处取极小值——这称为"位移元的下限性"，是 FEM 解收敛性质的理论保证。

### 3.4.9 条件极值与 Lagrange 乘子法

前面的泛函极值问题要求自变函数光滑且满足给定边界条件，无其他附加条件——称为**无条件极值问题**。具有附加约束条件的泛函极值问题称为**条件极值问题**。本节先回顾熟悉的多元函数条件极值，再推广到泛函情形。

#### 一、函数的条件极值（复习）

**问题**：求函数 $f(x,y,z)$ 在约束 $g_1(x,y,z)=0$、$g_2(x,y,z)=0$ 下的极值。

**方法**：引入 Lagrange 乘子 $\lambda_1, \lambda_2$，构造辅助函数：

$$F(x,y,z,\lambda_1,\lambda_2) = f + \lambda_1 g_1 + \lambda_2 g_2$$

对所有变量求偏导=0，得到 $n+m$ 个方程（$n$ 个原变量 + $m$ 个乘子），联立求解。

**例**（对应 HW2 Q4）：在椭球面 $16x^2+4y^2+z^2=16$ 与平面 $x+y+z=1$ 的交线上，求 $z$ 的极值。

- $f = z$，$g_1 = x+y+z-1$，$g_2 = 16x^2+4y^2+z^2-16$
- 构造 $F = z + \lambda_1(x+y+z-1) + \lambda_2(16x^2+4y^2+z^2-16)$
- 令 $\partial F/\partial x = \partial F/\partial y = \partial F/\partial z = 0$，联立约束求解

> 核心思想：**乘子法将约束问题转化为无条件极值问题**——这一思想可以推广到泛函。

#### 二、泛函的条件极值（推广）

将上述思想从有限维推广到无穷维：函数 → 泛函，代数约束 → 积分约束。

**问题**：求泛函 $Q[y] = \int_a^b F(x,y,y')dx$ 在等周条件 $\int_a^b \varphi(x,y,y')dx = \alpha$ 下的极值。

**方法**：引入乘子 $\lambda$，构造增广泛函：

$$\boxed{Q^*[y] = \int_a^b F^*(x,y,y',\lambda)\,dx - \lambda\alpha, \quad F^* = F + \lambda\varphi}$$

对 $Q^*$ 求变分=0，得到 Euler 方程：

$$\boxed{F_y^* - \frac{d}{dx}F_{y'}^* = 0}$$

通解含 2 个积分常数和 $\lambda$，由 2 个边界条件 + 等周条件确定。

**对比**：

| | 函数极值 | 泛函极值 |
|---|---|---|
| 对象 | $f(x_1,\ldots,x_n)$ | $Q[y] = \int F\,dx$ |
| 约束 | $g_i = 0$（代数方程） | $\int \varphi\,dx = \alpha$（积分约束） |
| 乘子 | $\lambda_i$（常数） | $\lambda$（常数） |
| 求解 | 偏导=0 → 代数方程组 | 变分=0 → Euler 方程（ODE） |
| 维度 | 有限维 | 无穷维 |

#### 三、具体算例：固定周长的最宽围栏

**问题**：在 $x=0$ 和 $x=1$ 处固定高度 $y(0)=y(1)=0$，围栏总弧长固定为 $L$。求曲线 $y(x)$ 使围出的面积最大。

**泛函**（面积）：$Q[y] = \int_0^1 y\,dx$

**约束**（弧长）：$\int_0^1 \sqrt{1+y'^2}\,dx = L$

**Step 1**：构造增广泛函

$$Q^*[y] = \int_0^1 y\,dx + \lambda\left[\int_0^1 \sqrt{1+y'^2}\,dx - L\right] = \int_0^1 \left(y + \lambda\sqrt{1+y'^2}\right)dx - \lambda L$$

其中 $F^* = y + \lambda\sqrt{1+y'^2}$。

**Step 2**：列 Euler 方程

$F^*_y = 1$，$F^*_{y'} = \dfrac{\lambda y'}{\sqrt{1+y'^2}}$

$$1 - \frac{d}{dx}\left(\frac{\lambda y'}{\sqrt{1+y'^2}}\right) = 0$$

**Step 3**：求解

由于 $F^*$ 不显含 $x$，可用首次积分：$F^* - y'F^*_{y'} = C$

$$y + \lambda\sqrt{1+y'^2} - y'\cdot\frac{\lambda y'}{\sqrt{1+y'^2}} = C$$

化简得：$y + \dfrac{\lambda}{\sqrt{1+y'^2}} = C$

令 $y' = \sinh t$，最终得解为**圆弧**——这正是等周定理的结论：固定周长下面积最大的封闭曲线是圆。

> 💡 **理解关键**：$\lambda$ 的物理意义——它是弧长约束的"代价"。如果 $L$ 增大（允许更长的围栏），$\lambda$ 会变化，解出的圆弧半径也跟着变。在力学中，$\lambda$ 通常对应约束反力。

---

## 3.5 更一般的 Euler 方程

### 3.5.1 含高阶导数的泛函

对于 $Q = \int_a^b F(x, y, y', y'')dx$，变分为：
$$\delta Q = \int_a^b (F_y\delta y + F_{y'}\delta y' + F_{y''}\delta y'')dx = 0$$

对 $F_{y''}\delta y''$ 项分部积分两次：
$$\int_a^b F_{y''}\delta y'' dx = \left.\delta y'F_{y''}\right|_a^b - \left.\delta y\frac{d}{dx}F_{y''}\right|_a^b + \int_a^b \delta y\frac{d^2}{dx^2}F_{y''}dx$$

假设 $\delta y'(a) = \delta y'(b) = 0$ 和 $\delta y(a) = \delta y(b) = 0$，得：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0}$$

**推广到 $n$ 阶导数**（通用形式）：
$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} - \cdots + (-1)^n\frac{d^n}{dx^n}F_{y^{(n)}} = 0}$$

> 💡 理解关键：高阶 Euler 方程的符号规律——正负交替，符号为 $(-1)^k$，每一项是第 $k$ 阶全导数作用于 $F_{y^{(k)}}$。记忆口诀："零阶正，一阶负，二阶正，三阶负……"，类似 Taylor 展开的交替符号。


### 3.5.2 含多个独立函数的泛函

对于 $Q[y_1(x), \ldots, y_n(x)] = \int_a^b F(x, y_1, \ldots, y_n, y_1', \ldots, y_n')dx$：

变分后得：
$$\delta Q = \int_a^b \left[\delta y_1\left(F_{y_1} - \frac{d}{dx}F_{y_1'}\right) + \cdots + \delta y_n\left(F_{y_n} - \frac{d}{dx}F_{y_n'}\right)\right]dx = 0$$

由于 $\delta y_i$ 可独立任意选取，令 $\delta y_2 = \cdots = \delta y_n = 0$，得：
$$\boxed{F_{y_i} - \frac{d}{dx}F_{y_i'} = 0,\quad i = 1, 2, \ldots, n}$$

即每个 $y_i$ 独立满足 Euler 方程。

> 💡 理解关键：多函数泛函得到的是 Euler 方程组，而不是单个方程。关键是 $\delta y_1, \ldots, \delta y_n$ 可以独立变化——这保证了每个方程必须独立成立。在平面弹性力学中，$u$ 和 $v$ 两个位移分量各自对应一个 Euler 方程。

### 3.5.3 含多元函数的泛函

对于 $Q[z(x,y)] = \iint_D F[x, y, z(x,y), p, q]\,dxdy$，其中 $p=\partial z/\partial x$，$q=\partial z/\partial y$。这里 $D$ 是二维平面上的有界区域，$\partial D$ 是它的边界（一条封闭曲线）。

**Euler 方程**：
$$\boxed{F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0}$$

**推导要点**：利用 Green 公式（此处指高斯散度定理 $\int_V(ab)_{,j}dV=\int_S ab\,l_j dS$，非 §3.8.2.3 的本构关系 Green 公式）将面积分转化为边界积分，利用 $\delta z$ 在边界上为零消去边界项。

**例**：$Q[z] = \iint_D (z_x^2 + z_y^2)dxdy$，边界条件 $z|_{\partial D} = \varphi(x,y)$

$F = p^2+q^2$，$F_p=2p$，$F_q=2q$，代入 Euler 方程：
$$\frac{\partial^2 z}{\partial x^2} + \frac{\partial^2 z}{\partial y^2} = 0$$

这就是著名的 **Laplace 方程**！

> 🔗 跨章连接：这是一个里程碑式的例子——它说明变分法可以把一个物理 PDE（Laplace 方程）等价转化为一个泛函极值问题。反过来想：如果你要解 Laplace 方程但边界复杂、解析解不存在，你可以转而去找那个泛函的极小值——这正是 FEM 做的事。Laplace 方程在 FEM 中对应最简单的 Poisson 问题（稳态热传导、静电场等），是最基础的 FEM 入门算例。

---

## 3.6 变分法在力学中的应用

### 3.6.1 Fermat 原理与广义相对论

- **Fermat 原理**：在均匀介质中，光线沿所需时间最少的光路传播。
- **Einstein 广义相对论**：光线在四维 Riemann 空间中沿所需时间最短的路线传播。

变分法的应用不仅限于经典物理和工程学科，而是延伸到量子场论、控制论、信息论等现代高科技领域。

### 3.6.2 Hamilton 原理

当 $t=t_1$ 和 $t=t_2$ 时，质点分别位于 A 和 B。则其真实运动路径应使以下**作用量积分**取极值：
$$\boxed{S = \int_{t_1}^{t_2} (T - U)dt}$$

其中 $T$ 是动能，$U$ 是势能，$T-U$ 称为 **Lagrange 函数**。

Hamilton 原理（1834）是力学的基本原理，可推出 Newton 三定律、能量守恒、动量守恒和角动量守恒。

> 💡 理解关键：Hamilton 原理中的 $T-U$（动能减势能）为什么不是 $T+U$？因为 $S = \int (T-U)dt$ 取极小值意味着系统尽量让动能小、势能大——这听起来反直觉。但物理上，"动能小"意味着系统"不愿"快速运动（惯性），"势能大"意味着系统"想"往低势能处走。这个 $T-U$ 的符号恰好产生了"惯性力和恢复力竞争"的效果。

> 🔗 跨章连接：Hamilton 原理和最小势能原理（§3.8.4）形式非常相似，但 Hamilton 原理是动力学的（含时间、动能），最小势能原理是静力学的（不含时间、只有势能）。FEM 的静力学分析基于最小势能原理，动力学分析（模态分析、瞬态响应）则基于 Hamilton 原理。

### 3.6.3 Lagrange 方程

对于具有 $K$ 个自由度的系统，广义坐标为 $q_1, \ldots, q_K$，广义速度为 $\dot{q}_1, \ldots, \dot{q}_K$。

Lagrange 函数 $L = T - U = L(t, q, \dot{q})$

由 Hamilton 原理和 Euler 方程：
$$\boxed{\frac{\partial L}{\partial q_i} - \frac{d}{dt}\frac{\partial L}{\partial \dot{q}_i} = 0,\quad i=1,2,\ldots,K}$$

这就是著名的 **Lagrange 方程**——在广义坐标下描述运动基本规律的方程。

**Lagrange 方程的优势**：只需写出广义坐标下的动能和势能表达式，代入即可，不需要画受力图。

### 3.6.4 单摆

**例**：质量为 $m$ 的质点悬挂在长度为 $l$ 的轻绳末端，在重力场中运动。

取广义坐标为 $\theta$：
- 动能：$T = \frac12mv^2 = \frac12m(l\dot{\theta})^2$
- 势能（$\theta=0$ 设为零）：$U = mgl(1-\cos\theta)$

Lagrange 函数：$L = \frac12ml^2\dot{\theta}^2 - mgl(1-\cos\theta)$

代入 Lagrange 方程：
$$\frac{\partial L}{\partial\theta} - \frac{d}{dt}\frac{\partial L}{\partial\dot{\theta}} = -mgl\sin\theta - ml^2\ddot{\theta} = 0$$

$$\ddot{\theta} + \frac{g}{l}\sin\theta = 0$$

小角度时 $\sin\theta \approx \theta$，得简谐振动 $\ddot{\theta} + \frac{g}{l}\theta = 0$。

### 3.6.5 两端固定弦的自由横向振动

建立坐标系，设 $u(x,t)$ 为点 $x$ 在 $t$ 时刻的横向位移。

考虑弦上一微段变形后长度 $dl = \sqrt{1+u_x^2}dx$，伸长为 $(\sqrt{1+u_x^2}-1)dx \approx \frac12u_x^2dx$。

弹性势能：$U = \int_0^l \frac{K}{2}u_x^2 dx$
动能：$T = \int_0^l \frac12\rho(x)u_t^2 dx$

作用量泛函：
$$S = \int_{t_1}^{t_2}\int_0^l\left[\frac12\rho u_t^2 - \frac12K u_x^2\right]dx\,dt$$

Euler 方程给出：
$$\frac{\partial}{\partial x}(Ku_x) - \frac{\partial}{\partial t}(\rho u_t) = 0 \quad\Rightarrow\quad u_{tt} = \frac{K}{\rho(x)}u_{xx}$$

这就是经典的**弦振动方程**。

> 💡 理解关键：注意这里的 Euler 方程涉及两个自变量 $(x, t)$，被积函数含偏导数 $u_x$ 和 $u_t$。这属于 §3.5.3 多元函数泛函的应用。从 Hamilton 原理（$S = \int(T-U)dt$）导出弦振动方程的过程，完美展示了变分法如何从"能量"出发自然地给出 PDE——不需要对微元体做受力分析。

---

## 3.7 结论

从以上例子可以总结变分法的主要步骤：

1. **建立泛函**及其边界条件（基于物理问题）
2. **利用变分法预备定理**，通过泛函变分得到 Euler 方程
3. **求解 Euler 方程**——这是一个微分方程问题

注意：变分法和 Euler 方程描述的是同一个物理问题。从变分法和 Euler 方程出发获得近似解具有相同的效果。

Euler 方程往往很难求解。但有些问题已知微分方程但很难求解，如果能转化为相应的泛函极值问题，就可以用近似方法（如有限元分析）方便地求解。

> 然而，**并非每个微分方程都能找到合适的泛函**。这种情况下，我们必须使用其他方法获得近似解，如最小二乘法和 Galerkin 法——这些将在后续章节讨论。

> 🔗 跨章连接：这一段给出了 FEM 最重要的一张路线图。**能找到泛函 → Ritz 法（对泛函求极小）；找不到泛函 → Galerkin 法（直接处理微分方程的加权残量）**。第 4 章将从这两个方向分别展开。考试大概率会问"Ritz 法和 Galerkin 法的区别是什么"——答案的源头就在这里。


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

---

## 3.8 弹性力学变分框架

> 本节是**变分法通向 FEM 的桥梁**——源文件 §5 的全部内容重新整理。考试重点：虚功原理（证明与含义）、最小势能原理、应变能与余能。

> ⚠️ 重难点：§3.8 是整章最有考试价值的部分。前面的泛函/Euler 方程是数学工具，这里才是工程核心。三个递进关系必须理清：① 虚功原理（最一般，不涉及本构）→ ② 虚位移/虚应力原理（分别是虚功原理取特定变分方向的特例）→ ③ 最小势能原理（引入线弹性本构后，从虚位移原理+Green公式导出）。考试证明题大概率考这三个原理之间的推导关系。

### 3.8.1 微分法 vs 变分法

弹性力学问题有两种解法路线：

**微分法**（前三章已详述）：
- 对微元体列平衡方程、几何方程、本构方程
- 得到 15 个偏微分方程，在给定边界条件下求解
- 难点：一旦边界稍复杂，解析精确解几乎不可能

**变分法**（本章核心）：
- 从整个弹性系统的**能量关系**出发
- 将弹性力学问题转化为**泛函极值问题**
- 在给定的约束条件下求解泛函的极值

> 由于泛函与弹性系统都涉及能量，弹性力学的变分原理也称为**能量原理（Energy Principles）**，变分法也称为**能量法**。

变分法的两种求解策略：

| 方法 | 路线 | 代表人物 | 意义 |
|------|------|----------|------|
| **Euler 法**（间接法） | 变分问题 → Euler 方程 → 微分方程求解 | Euler, Lagrange | 证明变分法与微分法的**等价性** |
| **直接法**（Direct Method） | 直接求泛函极值 → 近似解 | Ritz, Galerkin | FEM 的理论基础 |

> 关键认知：并非每个微分方程都能找到对应的泛函。遇到这种情况需要用加权残量法（如 Galerkin 法）直接处理微分方程。

> 💡 理解关键：Euler 法（间接法）的价值在于"证明两种路线等价"，而不是"实际解题"。实际工程中没有人会走 Euler 法——把变分问题变回微分方程再求解，等于白费劲。真正的实用路线是直接法：不经过微分方程，直接在泛函层面做极小化。

---

### 3.8.2 应变能与余能

#### 3.8.2.1 应变能密度 $W$

对小变形理想弹性体，定义**应变能密度** $W$（单位体积的应变能）：

$$\boxed{W = \int_{0}^{\varepsilon_{ij}} \sigma_{kl}\,d\varepsilon_{kl}}$$

**几何意义**（以单向拉伸 $\sigma_x$-$\varepsilon_x$ 曲线为例）：$W$ 是应力-应变曲线与横轴围成的面积。

对线弹性体，应力-应变呈线性关系：
$$\boxed{W = \frac{1}{2}\sigma_{kl}\varepsilon_{kl} = \frac{1}{2}C_{ijkl}\,\varepsilon_{ij}\varepsilon_{kl}}$$

#### 3.8.2.2 余能密度 $W^*$

$$\boxed{W^* = \int_{0}^{\sigma_{ij}} \varepsilon_{kl}\,d\sigma_{kl}}$$

**几何意义**：应力-应变曲线与纵轴围成的面积。

显然有：
$$\boxed{W + W^* = \sigma_{kl}\varepsilon_{kl}}$$

对线弹性体：$W = W^*$（数值相等，但概念不同），且：
$$W^* = \frac{1}{2}d_{ijkl}\,\sigma_{ij}\sigma_{kl}$$

> 注意：$W$ 和 $W^*$ 都是状态函数，与加载路径无关。这要求变形体是理想弹性体且初始状态无应力、无变形。

> ❌ 易错点：线弹性时 $W = W^*$ 数值相等，但定义不同——$W$ 用应变表示、$W^*$ 用应力表示。考试不要写"余能就是应变能"，而是要区分：$W(\varepsilon)$ 和 $W^*(\sigma)$ 是不同自变量的函数，恰好在 $\sigma = C\varepsilon$ 线性关系下数值相等。

#### 3.8.2.3 Green 公式（本构关系）

> ⚠️ 注意区分：本节的 Green 公式（$\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$）与推导中用到的 Green 公式（高斯散度定理 $\int_V(ab)_{,j}dV=\int_S ab\,l_j dS$）是两个不同的公式，都归功于 George Green。前者用于建立本构关系，后者用于分部积分。

从应变能密度定义可直接导出：

$$\boxed{\sigma_{ij} = \frac{\partial W}{\partial \varepsilon_{ij}}}$$

$$\boxed{\varepsilon_{ij} = \frac{\partial W^*}{\partial \sigma_{ij}}}$$

这两个公式极其重要——它们将**应力-应变关系**用能量函数的导数形式表示。这也是连接虚功原理与最小势能原理的关键。

> 💡 理解关键：Green 公式的美妙之处在于把 6 个独立的应力-应变关系塞进了一个能量函数 $W$ 里。换句话说，只要你知道 $W(\varepsilon)$ 的表达式，对 $\varepsilon_{ij}$ 求一次偏导就能得到 $\sigma_{ij}$。这在高维本构模型中极其强大——你不用背 21 个独立的弹性常数，只需知道能量函数的"形状"。而能量函数的对称性和不变性天然地保证了本构关系的对称性和客观性。

**物理直觉**：以一维杆件为例，$U = \int_0^\Delta N\,d\Delta$ 是外力功（= 应变能）。由能量守恒，$N = dU/d\Delta$——力的"驱动效应"等于能量对位移的变化率。Green 公式正是这一直觉的多维推广。


#### 3.8.2.4 弹性系统的总势能 $\Pi$

弹性系统 = **弹性体** + **载荷系统** + **支撑系统**。

载荷系统包括：
- 体力 $F_i$（在域 $V$ 内），其中 $F_i$ 表示单位体积上的体力分量（如重力），$i=1,2,3$ 对应 $x,y,z$ 方向。
- 面力 $\bar{p}_i$（在力边界 $S_\sigma$ 上）

对于**保守力系统**（力做功只与位移终值有关，与路径无关），外力势能 $V$ 定义为：
$$\boxed{V = -\int_V F_i u_i\,dV - \int_{S_\sigma} \bar{p}_i u_i\,dS}$$

> 保守力的物理意义：力对系统做正功 → 系统势能减少 → $A = -V$

**总势能**：应变能 + 外力势能
$$\boxed{\Pi = U + V = \int_V W\,dV - \int_V F_i u_i\,dV - \int_{S_\sigma} \bar{p}_i u_i\,dS}$$


> ❌ 易错点：外力势能 $V$ 的符号是负的（$V = -\int F_i u_i dV - \cdots$）！因为"外力做正功 = 系统势能减少"。如果写 $V = \int F_i u_i dV + \cdots$（正号），那 $\Pi = U - V$ 求极小意味着外力尽量多做功——系统会发散到无穷远，完全错误。

---

### 3.8.3 虚功原理（Principle of Virtual Work）

> **考试核心考点**——必须掌握证明过程和结论含义。

> ⚠️ 重难点：虚功原理是"整个弹性力学变分框架的总根"。它不依赖任何本构假设（线弹性、非线弹性、塑性都可以），是纯粹从平衡和几何条件导出的。理解这一点才能理解为什么虚位移原理等价于平衡方程、虚应力原理等价于几何方程——因为虚功原理本身就把平衡（通过应力）和几何（通过应变）"耦合"在一起了。

#### 3.8.3.1 可能位移与可能应力

首先约定弹性体的边界划分：$S = S_u \cup S_\sigma$，其中 $S_u$ 是位移已知的边界，$S_\sigma$ 是外力已知的边界（两者互不重叠）。符号 $\bar{u}_i$ 和 $\bar{p}_i$ 分别表示边界上给定的已知位移和已知面力。

FEM 思想的核心逻辑：在"满足部分条件的解"中，找到"满足剩余条件的解"。

**可能位移（几何容许位移）**：满足几何方程 + 位移边界条件，但不一定满足平衡条件。
$$u_i \in \{u_i \;|\; \varepsilon_{ij} = \tfrac12(u_{i,j}+u_{j,i})\;\text{在}\;V\;\text{内},\; u_i = \bar{u}_i\;\text{在}\;S_u\;\text{上}\}$$

**可能应力（静力容许应力）**：满足平衡方程 + 力边界条件，但不一定满足变形协调。
$$\sigma_{ij} \in \{\sigma_{ij}\;|\; \sigma_{ij,j}+F_i = 0\;\text{在}\;V\;\text{内},\; \sigma_{ij}l_j = \bar{p}_i\;\text{在}\;S_\sigma\;\text{上}\}$$

其中 $l_j$ 是边界外法线方向与 $j$ 坐标轴夹角的方向余弦（$l_1,l_2,l_3$ 即 $n_x,n_y,n_z$）。

> 核心理解：可能位移有无穷多种，可能应力也有无穷多种。只有同时满足所有条件的才是**真实解**。变分法的思路就是：先从大量"可能"解中选取，再通过极值条件筛选出真实解。

> 💡 理解关键："可能位移"和"可能应力"这两个概念是理解后续一切的基础。它们分别只满足三组方程中的两组（可能位移 → 几何+位移边界；可能应力 → 平衡+力边界），第三组靠变分原理来"选"。这就是 FEM 分两步走的思想：① 构造满足几何条件的位移场（形函数）；② 用最小势能原理确定未知参数。第一步保证"可能位移"，第二步从可能位移中选出真实位移。

#### 3.8.3.2 虚位移与虚应力

**虚位移** $\delta u_i$：从一个可能位置到邻近可能位置的**无限小位移**。
- 满足变形连续条件：$\delta\varepsilon_{ij} = \frac12(\delta u_{i,j} + \delta u_{j,i})$
- 在 $S_u$ 上：$\delta u_i = 0$（位移边界已给定，不可变）

**虚应力** $\delta\sigma_{ij}$：从一种可能应力到另一种可能应力的**无限小变分**。
- 在 $V$ 内：$\delta\sigma_{ij,j} = 0$（平衡条件不变）
- 在 $S_\sigma$ 上：$\delta\sigma_{ij}l_j = 0$（力边界已给定，不可变）

> ❌ 易错点："虚"字的含义经常被误解。虚位移不是"假想的、不存在的位移"——它是真实可能的位移之间的差！叫"虚"是因为我们想象在某一瞬间"冻结"载荷，让系统从一个可能位移走到另一个可能位移。类似地，虚功也不是"不存在的功"，而是真实力在"虚位移"上做的功。这个"虚"字更接近"变分"（$\delta u_i$）的物理含义。

#### 3.8.3.3 虚功方程

对于在外力 $(F_i, \bar{p}_i)$ 作用下处于平衡的可变形体，可能应力 $\sigma_{ij}$ 和可能位移 $u_i$（及相应 $\varepsilon_{ij}$）满足：

$$\boxed{\int_V F_i u_i\,dV + \int_S p_i u_i\,dS = \int_V \sigma_{ij}\varepsilon_{ij}\,dV}$$

分解边界 $S = S_u \cup S_\sigma$：
$$\int_V F_i u_i\,dV + \int_{S_u} p_i \bar{u}_i\,dS + \int_{S_\sigma} \bar{p}_i u_i\,dS = \int_V \sigma_{ij}\varepsilon_{ij}\,dV$$

> 关键注意：虚功原理**不涉及本构关系**——适用于一切小变形体的任何材料（线弹性、非线性弹性、弹塑性均可）。


#### 3.8.3.4 充分性证明（正向）

**已知**：$\sigma_{ij}$ 满足平衡条件，$u_i$ 满足几何条件。**证**：虚功方程恒成立。

**Step 1** — 利用 Cauchy 公式和高斯公式：
$$\begin{aligned}
\int_S p_i u_i\,dS &= \int_S \sigma_{ij}l_j u_i\,dS = \int_V (\sigma_{ij}u_i)_{,j}\,dV \\
&= \int_V \sigma_{ij,j}u_i\,dV + \int_V \sigma_{ij}u_{i,j}\,dV
\end{aligned}$$

**Step 2** — 利用应力张量对称性和几何方程：
$$\int_V \sigma_{ij}u_{i,j}\,dV = \int_V \sigma_{ij} \cdot \frac12(u_{i,j} + u_{j,i})\,dV = \int_V \sigma_{ij}\varepsilon_{ij}\,dV$$

**Step 3** — 代入并利用平衡方程 $\sigma_{ij,j} + F_i = 0$：
$$\begin{aligned}
\int_S p_i u_i\,dS &= \int_V (-F_i)u_i\,dV + \int_V \sigma_{ij}\varepsilon_{ij}\,dV \\
\int_V F_i u_i\,dV + \int_S p_i u_i\,dS &= \int_V \sigma_{ij}\varepsilon_{ij}\,dV \quad\blacksquare
\end{aligned}$$


#### 3.8.3.5 虚位移原理（等价于平衡方程）

**表述**：对于给定外力的可变形体，若某应力场对**一切可能位移**都使虚功方程成立，则此应力场是平衡容许的。

**证明核心步骤**（仅需理解思路）：
- 任取两个相邻的可能位移 $u_i^1$、$u_i^2$，令 $\delta u_i = u_i^1 - u_i^2$（虚位移）
- 虚功方程相减得：
  $$\int_V F_i\delta u_i\,dV + \int_{S_\sigma}\bar{p}_i\delta u_i\,dS = \int_V \sigma_{ij}\delta\varepsilon_{ij}\,dV \tag{*}$$

- **对右边分部积分**（关键推导）：

  **Step 1**：代入几何关系 $\delta\varepsilon_{ij} = \frac12(\delta u_{i,j} + \delta u_{j,i})$

  $$\int_V \sigma_{ij}\delta\varepsilon_{ij}\,dV = \int_V \sigma_{ij}\cdot\frac12(\delta u_{i,j} + \delta u_{j,i})\,dV$$

  **Step 2**：利用应力对称性 $\sigma_{ij} = \sigma_{ji}$，两项合并

  $$= \int_V \sigma_{ij}\,\delta u_{i,j}\,dV$$

  **Step 3**：分部积分（高斯公式 $\int_V (a\,b)_{,j}\,dV = \int_S a\,b\,l_j\,dS$）

  $$= \int_S \sigma_{ij}l_j\,\delta u_i\,dS - \int_V \sigma_{ij,j}\,\delta u_i\,dV$$

  **Step 4**：分解边界 $S = S_u \cup S_\sigma$，利用 $\delta u_i|_{S_u} = 0$

  $$= \int_{S_\sigma} \sigma_{ij}l_j\,\delta u_i\,dS - \int_V \sigma_{ij,j}\,\delta u_i\,dV$$

- 代回 $(*)$ 式，移项整理：

  $$\int_V (F_i + \sigma_{ij,j})\delta u_i\,dV + \int_{S_\sigma}(\bar{p}_i - \sigma_{ij}l_j)\delta u_i\,dS = 0$$

- 由变分法预备定理，被积函数恒为零 → 平衡方程 $F_i + \sigma_{ij,j} = 0$ + 力边界 $\bar{p}_i = \sigma_{ij}l_j$

**结论**：
$$\boxed{\text{虚位移原理} \;\Longleftrightarrow\; \text{平衡方程} \;+\; \text{力边界条件}}$$

> 💡 理解关键：虚位移原理的逻辑是——如果一个应力场对"所有可能位移"都满足虚功方程，那它必然是满足平衡方程和力边界的。换句话说，用"无数个位移"去检验一个应力场，如果每次外力功都等于内力功，说明这个应力场在内部是自平衡且边界力匹配的。FEM 的求解逻辑正是反过来的：从位移出发求应力，靠虚位移原理保证求出的应力"最接近"平衡。

#### 3.8.3.6 虚应力原理（等价于几何方程）

**表述**：对于给定外力的可变形体，若某位移场对**一切可能应力**都使虚功方程成立，则此位移场是几何容许的。

**证明核心**（对偶于虚位移原理）：
- 任取两个相邻的可能应力 $\sigma_{ij}^1$、$\sigma_{ij}^2$，令 $\delta\sigma_{ij} = \sigma_{ij}^1 - \sigma_{ij}^2$
- 虚功方程相减，利用 $\delta F_i = 0$、$\delta\bar{p}_i = 0$（外力给定）：
  $$\int_{S_u} \delta p_i\,\bar{u}_i\,dS = \int_V \delta\sigma_{ij}\,\varepsilon_{ij}\,dV \tag{**}$$

- **对左边利用高斯定理**（关键推导，与 3.8.3.5 对偶）：

  **Step 1**：代入 Cauchy 公式 $\delta p_i = \delta\sigma_{ij}l_j$

  $$\int_{S_u} \delta p_i\,\bar{u}_i\,dS = \int_{S_u} \delta\sigma_{ij}l_j\,\bar{u}_i\,dS$$

  **Step 2**：高斯散度定理，将面积分化为体积分

  $$= \int_V (\delta\sigma_{ij}\bar{u}_i)_{,j}\,dV$$

  **Step 3**：乘积求导法则展开

  $$= \int_V \delta\sigma_{ij,j}\,\bar{u}_i\,dV + \int_V \delta\sigma_{ij}\,\bar{u}_{i,j}\,dV$$

  **Step 4**：利用虚应力的平衡约束 $\delta\sigma_{ij,j} = 0$，第一项消失

  $$= \int_V \delta\sigma_{ij}\,\bar{u}_{i,j}\,dV$$

  **Step 5**：利用应力对称性 $\delta\sigma_{ij} = \delta\sigma_{ji}$，改写为应变形式

  $$= \int_V \delta\sigma_{ij}\cdot\frac12(\bar{u}_{i,j} + \bar{u}_{j,i})\,dV$$

- 代回 $(**)$ 式，移项整理：

  $$\int_V \delta\sigma_{ij}\left[\varepsilon_{ij} - \frac12(u_{i,j}+u_{j,i})\right]dV + \int_{S_u}\delta p_i(u_i - \bar{u}_i)dS = 0$$

- 由变分法预备定理，被积函数恒为零 → 几何方程 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ + 位移边界 $u_i = \bar{u}_i$

**结论**：
$$\boxed{\text{虚应力原理} \;\Longleftrightarrow\; \text{几何方程} \;+\; \text{位移边界条件}}$$

> 💡 理解关键：虚应力原理和虚位移原理是完美的"对偶关系"——一个判应力是否平衡、一个判位移是否协调。把表格倒过来看：虚位移原理：独立变分 $\delta u_i$ → 检验 $\sigma_{ij}$ → 输出平衡方程+力边界；虚应力原理：独立变分 $\delta\sigma_{ij}$ → 检验 $u_i$ → 输出几何方程+位移边界。考试大概率会问你这两个原理的"对偶性"体现在哪里。


#### 3.8.3.7 功的互等定理（Betti 公式）

对**线弹性体**（需要用到本构关系的对称性），取两种不同真实状态：

$$\boxed{\int_V F_i^{(1)} u_i^{(2)}\,dV + \int_S p_i^{(1)} u_i^{(2)}\,dS = \int_V F_i^{(2)} u_i^{(1)}\,dV + \int_S p_i^{(2)} u_i^{(1)}\,dS}$$

**文字表述**：状态一的力在状态二的位移上做的功 = 状态二的力在状态一的位移上做的功。

**推导**：由虚功方程 + Hooke 定律的对称性 $C_{ijkl}=C_{klij}$：
$$\sigma_{ij}^{(1)}\varepsilon_{ij}^{(2)} = C_{ijkl}\varepsilon_{kl}^{(1)}\varepsilon_{ij}^{(2)} = C_{klij}\varepsilon_{ij}^{(2)}\varepsilon_{kl}^{(1)} = \sigma_{kl}^{(2)}\varepsilon_{kl}^{(1)} = \sigma_{ij}^{(2)}\varepsilon_{ij}^{(1)}$$

> 💡 理解关键：Betti 互等定理的推导只有一行，但这一步用到了一个极深的事实：弹性张量 $C_{ijkl}$ 满足 $C_{ijkl} = C_{klij}$（主对称性）。这说明 $\sigma_{ij}^{(1)}\varepsilon_{ij}^{(2)} = \sigma_{ij}^{(2)}\varepsilon_{ij}^{(1)}$——状态一的应力在状态二应变上做的功 = 状态二的应力在状态一应变上做的功。互等定理在 FEM 中保证刚度矩阵 $K$ 是对称矩阵——这正是有限元计算高效的数学根源。

**应用示例**：任意形状弹性体两端受压力 $P$（作用线长度为 $L$，横截面积 $A$），求体积变化 $\Delta V$。

**方法一：直接由应变求**（推荐）

- 实际状态：$\sigma_x = -P/A$，$\sigma_y = \sigma_z = 0$
- 由 Hooke 定律：$\varepsilon_x = -\frac{P}{EA}$，$\varepsilon_y = \varepsilon_z = \frac{\nu P}{EA}$
- 体积应变：$\varepsilon_V = \varepsilon_x + \varepsilon_y + \varepsilon_z = \frac{P}{EA}(-1+2\nu)$

$$\boxed{\Delta V = V\varepsilon_V = -\frac{(1-2\nu)}{E}PL}$$

**方法二：用互等定理**（演示定理用法）

核心思路：选一个虚拟状态，使得 Betti 定理的一端直接给出 $\Delta V$。

**Step 1**：定义两个状态

| | 实际状态（State 1） | 虚拟状态（State 2） |
|---|---|---|
| **载荷** | 轴向压力 $P$ | 全表面单位压力 $q=1$（静水压力） |
| **应力** | $\sigma_x^{(1)}=-P/A$，其余为 0 | $\sigma_x^{(2)}=\sigma_y^{(2)}=\sigma_z^{(2)}=-1$ |
| **应变** | $\varepsilon_x^{(1)}=-P/(EA)$，$\varepsilon_y^{(1)}=\varepsilon_z^{(1)}=\nu P/(EA)$ | $\varepsilon_x^{(2)}=\varepsilon_y^{(2)}=\varepsilon_z^{(2)}=-(1-2\nu)/E$ |

**Step 2**：关键观察

体积变化 $\Delta V = V \cdot (\varepsilon_x^{(1)}+\varepsilon_y^{(1)}+\varepsilon_z^{(1)})$，而虚拟状态的应力满足 $\sigma_x^{(2)}=\sigma_y^{(2)}=\sigma_z^{(2)}=-1$，所以：

$$\int_V \sigma_{ij}^{(2)}\varepsilon_{ij}^{(1)}\,dV = \int_V (-1)\cdot(\varepsilon_x^{(1)}+\varepsilon_y^{(1)}+\varepsilon_z^{(1)})\,dV = -\Delta V$$

右端恰好就是 $-\Delta V$！这就是选择静水压力作为虚拟状态的原因。

**Step 3**：计算左端

只有 $\sigma_x^{(1)} \neq 0$，所以：

$$\int_V \sigma_{ij}^{(1)}\varepsilon_{ij}^{(2)}\,dV = \int_V \sigma_x^{(1)}\cdot\varepsilon_x^{(2)}\,dV = \left(-\frac{P}{A}\right)\cdot\left(-\frac{1-2\nu}{E}\right)\cdot V = \frac{P(1-2\nu)}{AE}\cdot V$$

**Step 4**：由 Betti 定理求 $\Delta V$

$$\frac{P(1-2\nu)}{AE}\cdot V = -\Delta V$$

$$\boxed{\Delta V = -\frac{(1-2\nu)}{E}PL}$$

> 💡 **为什么选静水压力？** 因为我们要求的是 $\varepsilon_V = \varepsilon_x+\varepsilon_y+\varepsilon_z$（体积应变），而静水压力 $\sigma_x=\sigma_y=\sigma_z=-1$ 恰好能在 Betti 定理中"提取"出这个和。这是 Betti 定理的精髓——**虚拟状态的选择取决于你想求什么量**。


---

### 3.8.4 最小势能原理（Principle of Minimum Potential Energy）

> FEM 的理论基石——"FEM 就是用最小势能原理在所有可能位移中找到真实位移"

> ⚠️ 重难点：最小势能原理是本章的最高峰——它把前面所有概念（虚功原理、虚位移原理、Green 公式（本构关系 $\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$）、总势能）串在一起，最终给出 $\delta\Pi = 0$。理解了这个推导，就理解了为什么 FEM 能工作。考试如果出一道大题让你"从虚功原理推导最小势能原理"，本节就是标准答案。

#### 3.8.4.1 从虚位移原理到最小势能原理

将 Green 公式 $\sigma_{ij} = \partial W/\partial\varepsilon_{ij}$ 代入虚位移原理方程：

$$\int_V F_i\delta u_i\,dV + \int_{S_\sigma} \bar{p}_i\delta u_i\,dS = \int_V \frac{\partial W}{\partial\varepsilon_{ij}}\delta\varepsilon_{ij}\,dV = \delta\int_V W\,dV = \delta U$$

左端恰好是 $-\delta V$（外力势能变分的负值）：
$$\delta V = \delta\left(-\int_V F_i u_i\,dV - \int_{S_\sigma}\bar{p}_i u_i\,dS\right) = -\left(\int_V F_i\delta u_i\,dV + \int_{S_\sigma}\bar{p}_i\delta u_i\,dS\right)$$

因此：
$$\boxed{\delta\Pi = \delta(U + V) = 0}$$

这就是总势能的驻值条件。可进一步证明：

$$\boxed{\delta^2\Pi \geq 0}$$

因此这是**极小值**而非一般驻值。

> 💡 理解关键：$\delta^2\Pi \geq 0$（而不是 $\delta^2\Pi \leq 0$ 或正负不定）是由应变能的正定性保证的。物理上，弹性体储存的应变能总是正的（$W = \frac{1}{2}\sigma_{ij}\varepsilon_{ij} \geq 0$），所以总势能取极小值而不是极大值。这也对应于力学中的"最小势能"——稳定平衡状态对应势能最小。

> ❌ 易错点：推导中的一个常见错误是在做 $\delta\int_V W dV$ 时忘记变分和积分的交换。正确步骤是：$\int_V \frac{\partial W}{\partial\varepsilon_{ij}}\delta\varepsilon_{ij} dV = \int_V \delta W dV = \delta\int_V W dV$。$\delta$ 可以和 $\int$ 交换（线性算子），但前提是被积函数满足一定连续性条件——在线弹性范围内总是满足的。

#### 3.8.4.2 完整表述

> 对于小变形弹性系统，在一切满足几何方程和位移边界条件的可能位移场中，真实的位移场使系统总势能取最小值。反过来说，使总势能取最小值的可能位移场就是真实位移场。


#### 3.8.4.3 与 FEM 的直接联系

有限元法的基本思想链条：
```
可能位移场 → 构造试探函数 → 代入总势能泛函
→ δΠ = 0 → 刚度方程 Ku = f → 求解节点位移
```

这就是后续第 4、5 章的核心内容。

> 🔗 跨章连接：这条链条的每一步对应后续章节的一个主题：① "可能位移场 + 构造试探函数" → 第 6 章（形函数）；② "代入总势能泛函" → 第 5 章（单刚矩阵）；③ "$\delta\Pi = 0$ → 刚度方程" → 第 4 章（Ritz 法）和第 5 章（总体集成）；④ "求解节点位移" → 第 4 章（线性方程组求解）。考试综合题可能要求你把这个链条串起来解释。

---

### 3.8.5 弹性力学的 Euler 方程（从总势能变分导出平衡 PDE）

#### 3.8.5.1 基本思路

从最小势能原理 $\delta\Pi = 0$ 出发，经过变分运算自动导出：
- **域内**平衡微分方程（Euler 方程）
- **力边界**上的自然边界条件
- **位移边界**上的本质边界条件

这证明了变分法与微分法的**等价性**。

> 💡 理解关键：这一节的目的是"反着走"——不是从微分方程推变分，而是从变分推微分方程。逻辑链条是：最小势能原理 $\delta\Pi = 0$ → 变分运算 → Euler 方程（平衡方程）+ 自然边界（力边界）+ 本质边界（位移边界）。这完整证明了"变分法描述"和"微分方程描述"是等价的两个语言。

#### 3.8.5.2 平面应力问题的完整推导

以平面应力问题为例（为清晰起见不用张量记号）。

**Step 1 — 导出总势能**

应变能（用位移表示）：
$$U = \frac{E}{2(1-\mu^2)}\iint_\Omega \left[\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + 2\mu\frac{\partial u}{\partial x}\frac{\partial v}{\partial y} + \frac{1-\mu}{2}\left(\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\right)^2\right]dxdy$$

> 注意：此处 $\mu$ 是泊松比，与前文使用的 $\nu$ 含义相同。全文对同一物理量混用了 $\mu$ 和 $\nu$ 两个符号，阅读时请注意。
总势能：
$$\Pi = U - \iint_\Omega (F_x u + F_y v)dxdy - \int_{\Gamma_\sigma} (\bar{p}_x u + \bar{p}_y v)ds$$

此处 $\Omega$ 为平面应力问题的求解域，$\Gamma_\sigma$ 为给定面力的力边界，$\Gamma_u$ 为给定位移的位移边界（记号 $\Gamma_\sigma,\Gamma_u$ 分别等同于前文的 $S_\sigma,S_u$）。

**Step 2 — 取一阶变分，分部积分**

对 $\delta U$ 中的各项利用 Green 公式（高斯散度定理）将面积分转化为边界积分，整理后得：

**Step 3 — 结果**

由变分法预备定理，$\delta u$ 和 $\delta v$ 的系数分别恒为零。

**域内（在 $\Omega$ 内）——Euler 方程（即平衡方程）**：
$$\boxed{\frac{E}{1-\mu^2}\left(\frac{\partial^2 u}{\partial x^2} + \frac{1-\mu}{2}\frac{\partial^2 u}{\partial y^2} + \frac{1+\mu}{2}\frac{\partial^2 v}{\partial x\partial y}\right) + F_x = 0}$$
$$\boxed{\frac{E}{1-\mu^2}\left(\frac{\partial^2 v}{\partial y^2} + \frac{1-\mu}{2}\frac{\partial^2 v}{\partial x^2} + \frac{1+\mu}{2}\frac{\partial^2 u}{\partial x\partial y}\right) + F_y = 0}$$

**力边界（在 $\Gamma_\sigma$ 上）——自然边界条件**：
$$\frac{E}{1-\mu^2}\left[l\left(\frac{\partial u}{\partial x} + \mu\frac{\partial v}{\partial y}\right) + m\frac{1-\mu}{2}\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)\right] = \bar{p}_x$$
$$\frac{E}{1-\mu^2}\left[m\left(\frac{\partial v}{\partial y} + \mu\frac{\partial u}{\partial x}\right) + l\frac{1-\mu}{2}\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)\right] = \bar{p}_y$$

**位移边界（在 $\Gamma_u$ 上）——本质边界条件**：
$$u = \bar{u},\quad v = \bar{v}$$

> 本质边界/自然边界的关键区分：本质边界条件是变分前就预设好的（位移边界给定值），变分时边界变分 $\delta u|_{\Gamma_u}=0$；而自然边界条件是从变分计算中"掉出来"的产物。


---

### 3.8.6 直接法算例

#### 3.8.6.1 一边固定三边自由薄板（Ritz 法）

**题目**：薄板一边固定（$x=0$），三边自由，受均匀剪应力 $\tau$（$y=0$ 和 $y=b$ 方向相反，$x=a$ 处竖向）。忽略体力。

**分析**：固定边条件 $(u)_{x=0}=0$，$(v)_{x=0}=0$ 是**本质边界条件**，试探函数必须满足。

选择位移函数（与 §3.8.5 的坐标系一致，$x$ 为固定边法向）：
$$u = x(A_1 + A_2x + A_3y + \cdots),\quad v = x(B_1 + B_2x + B_3y + \cdots)$$

**Ritz 法步骤**：

1. 取**一阶近似**（只保留 $A_1$、$B_1$）：
   $$u = A_1x,\quad v = B_1x$$

2. 计算应变能：
   $$U = \frac{E}{2(1-\nu^2)}\int_0^a\int_0^b\left[\left(\frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2 + 2\nu\frac{\partial u}{\partial x}\frac{\partial v}{\partial y} + \frac{1-\nu}{2}\left(\frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}\right)^2\right]dxdy$$
   $$= \frac{Eab}{2(1-\nu^2)}\left(A_1^2 + \frac{1-\nu}{2}B_1^2\right)$$

3. 外力势能（注意剪应力方向）：
   $$V = - \int_{\Gamma} (\bar{X}u + \bar{Y}v)d\Gamma = -\tau B_1 ab$$

4. 极值条件：
   $$\frac{\partial\Pi}{\partial A_1} = 0,\quad \frac{\partial\Pi}{\partial B_1} = 0$$
   得：
   $$\boxed{A_1 = 0,\quad B_1 = \frac{2(1+\nu)}{E}\tau}$$

5. 最终位移：
   $$u = 0,\quad v = \frac{2(1+\nu)}{E}\tau x$$

> 注意：一阶近似只是近似解。取更多项可提高精度——这正是 FEM 的基本思想。

> 💡 理解关键：这个 Ritz 法算例完整展示了"直接法"的全流程。关键步骤 4 中 $\partial\Pi/\partial A_1 = 0$、$\partial\Pi/\partial B_1 = 0$——注意这里要求的是**对未知系数**求偏导，而非对位移函数变分！因为 $\Pi(A_1, B_1)$ 代入试探函数后已经变成了 $A_1, B_1$ 的普通函数，求极值就降级为普通多元函数求极值。这是 Ritz 法的精髓：把泛函极值转化成代数方程组。

> ❌ 易错点：本题中三边是自由的（力边界 $\bar{p}=0$）——这不要求试探函数满足！Ritz 法只需要试探函数满足位移边界条件（本质边界），力边界条件在极值条件中自然满足。如果错误地试图让试探函数满足所有边界条件，会增加不必要的复杂度。

#### 3.8.6.2 简支梁均布荷载（Galerkin 法）

**题目**：简支梁受均布荷载 $q$，用 Galerkin 法求解挠度。

**关键区别**：Galerkin 法要求试探函数**既要满足位移边界条件，又要满足力边界条件**。

位移边界条件（简支）：$w|_{x=0}=w|_{x=l}=0$

力边界条件（简支端弯矩为零）：$w''|_{x=0}=w''|_{x=l}=0$（因为 $M = -EIw''$，其中 $EI$ 为梁的抗弯刚度，$E$ 为杨氏模量，$I$ 为截面惯性矩）

**Galerkin 法步骤**：

1. 构造满足所有边界条件的多项式：
   $$w = a_1x^4 + a_2x^3 + a_3x^2 + a_4x + a_5$$
   代入四个边界条件得：
   $$a_3 = a_5 = 0,\quad a_2 = -2a_1l,\quad a_4 = a_1l^3$$
   $$\Rightarrow w = a_1(x^4 - 2lx^3 + l^3x)$$

2. 令 $w_1 = x^4 - 2lx^3 + l^3x$。Galerkin 方程（权函数 = 基函数）：
   $$\int_0^l (EIw'''' - q) \cdot w_1\,dx = 0$$

3. 代入 $w = a_1w_1$，$w'''' = 24a_1$，积分得：
   $$a_1 = \frac{q}{24EI}$$

4. 最终解：
   $$\boxed{w = \frac{q}{24EI}(x^4 - 2lx^3 + l^3x)}$$

> 注意：这正是经典梁理论的**精确解**！因为基函数恰好包含了真实解的形式。

> ❌ 易错点：Ritz 法和 Galerkin 法对试探函数的要求不同！Ritz 法只需要满足**本质边界条件**（位移边界），Galerkin 法需要满足**全部边界条件**（位移边界 + 力边界）。这道题如果不小心用 Ritz 法的试探函数（只满足 $w=0$）来做 Galerkin，就会出问题。


#### 3.8.6.3 矩形薄板位移边值问题（Galerkin 法）

**题目**：矩形薄板（$2a$ 宽 $\times b$ 高），左、右、下边固定，上边给定抛物线形压缩位移 $v = -\eta(1-x^2/a^2)$，忽略体力。

**选择位移函数**（基于变分原理的对称性考虑）：
$$u = A_1\left(1-\frac{x^2}{a^2}\right)\frac{x}{a}\frac{y}{b}\left(1-\frac{y}{b}\right)$$
$$v = -\eta\left(1-\frac{x^2}{a^2}\right)\frac{y}{b} + B_1\left(1-\frac{x^2}{a^2}\right)\frac{y}{b}\left(1-\frac{y}{b}\right)$$

验证边界条件：
- $x=\pm a$：$u=0$，$v$ 中两项均有 $(1-x^2/a^2)$ 因子 → $v=0$ ✓
- $y=0$：$u=0$（因子 $y/b$），$v=0$ ✓
- $y=b$：$v = -\eta(1-x^2/a^2)$（$B_1$ 项有 $(1-y/b)$ 因子为零）✓

由 Galerkin 方程解得：
$$\boxed{A_1 = \frac{35(1+\nu)\eta}{42\frac{b}{a} + 20(1-\nu)\frac{a}{b}},\quad B_1 = \frac{5(1-\nu)\eta}{16\frac{a^2}{b^2} + 2(1-\nu)}}$$

> 💡 理解关键：构造试探函数是 Ritz/Galerkin 法最"艺术"的一步——需要同时满足：① 所有边界条件；② 尽量简单（项数少）；③ 能反映解的定性特征（对称性、奇偶性）。这道题中 $(1-x^2/a^2)$ 因子巧妙地在 $x=\pm a$ 处自动归零；$y/b$ 因子在 $y=0$ 处归零；$(1-y/b)$ 因子在 $y=b$ 处归零。这种"因子化构造法"是考试和实际 FEM 中构造形函数的标准技巧。

#### 3.8.6.4 简支梁变分推导平衡方程与边界条件

**题目**：简支梁受均布荷载 $q$，从总势能出发，导出平衡微分方程和支座的力边界条件。

**总势能**（忽略剪力变形）：
$$\Pi = \frac{EI}{2}\int_0^l \left(\frac{d^2w}{dx^2}\right)^2 dx - \int_0^l q w\,dx$$

**$\delta\Pi = 0$**：
$$\delta\Pi = EI\int_0^l \frac{d^2w}{dx^2}\,\delta\!\left(\frac{d^2w}{dx^2}\right)dx - \int_0^l q\delta w\,dx = 0$$

**分部积分两次**（关键技巧）：
$$\begin{aligned}
EI\int_0^l \frac{d^2w}{dx^2}\,\delta\!\left(\frac{d^2w}{dx^2}\right)dx &= EI\int_0^l \frac{d^2w}{dx^2}\,d\!\left[\delta\!\left(\frac{dw}{dx}\right)\right] \\
&= \left.EI\frac{d^2w}{dx^2}\delta\!\left(\frac{dw}{dx}\right)\right|_0^l - EI\int_0^l \frac{d^3w}{dx^3}\,\delta\!\left(\frac{dw}{dx}\right)dx \\
&= \left.EI\left[\frac{d^2w}{dx^2}\delta\!\left(\frac{dw}{dx}\right) - \frac{d^3w}{dx^3}\delta w\right]\right|_0^l + EI\int_0^l \frac{d^4w}{dx^4}\delta w\,dx
\end{aligned}$$

**结果**：
$$\boxed{EI\frac{d^4w}{dx^4} - q = 0}\quad\text{——域内平衡方程（Euler 方程）}$$
$$\boxed{\left.\frac{d^2w}{dx^2}\right|_{x=0,l} = 0}\quad\text{——力边界条件（简支端弯矩为零，自然边界）}$$

> $\delta w|_{x=0,l}=0$ 不是推导结果，而是**本质边界条件**（简支约束，预设的）。力边界条件 $\frac{d^2w}{dx^2}=0$ 是变分过程中"掉出来"的**自然边界条件**。


> 💡 理解关键：分部积分两次后会产生两个边界项——$\delta w'$ 的系数 $w''$（弯矩）和 $\delta w$ 的系数 $w'''$（剪力）。简支端：$\delta w = 0$（本质边界）→ $w'''$ 的边界项自动消失；$\delta w' \neq 0$（转角可自由变化）→ $w''=0$（弯矩为零）是自然边界条件。固定端：$\delta w = 0$ 且 $\delta w' = 0$ → 两个边界项都自动消失，弯矩和剪力都由支反力提供，不作为自然边界出现。

---

### 3.8.7 总结：变分原理体系图

```
                   弹性力学基本方程
                        |
          +-------------+-------------+
          |                           |
    微分法路线                      变分法路线
  (平衡+几何+本构)              (能量泛函极值)
          |                           |
    15个PDE求解              +--------+--------+
    (几乎不可能)             |                  |
                      Euler法(间接法)    直接法(Ritz/Galerkin)
                     变分→微分方程        直接求泛函极值
                          |                  |
                    证明等价性            FEM的基础
```


| 原理 | 独立变分对象 | 等价于 | 涉及的场 |
|------|------------|--------|---------|
| 虚位移原理 | $\delta u_i$（位移变分） | 平衡方程 + 力边界条件 | 可能位移场 |
| 虚应力原理 | $\delta\sigma_{ij}$（应力变分） | 几何方程 + 位移边界条件 | 可能应力场 |
| 最小势能原理 | $u_i$（位移） | 真实位移使 $\Pi$ 最小 | 几何容许位移 |
| 最小余能原理 | $\sigma_{ij}$（应力） | 真实应力使 $\Pi^*$ 最小 | 静力容许应力 |

> 核心记忆口诀：**虚位移判应力，虚应力判位移；最小势能找位移，最小余能找应力。**


> 🔗 跨章连接：这个表格中的四个原理在整个 FEM 课程中以不同方式出现。最小势能原理 → 第 5 章位移元（$Kd = f$），最小余能原理 → 应力杂交元（更高级的内容），虚位移原理 → 第 4 章推导 Galerkin 弱形式时反复使用的技巧，虚应力原理 → 第 6 章分片试验的理论基础。

---

## 检查你的理解（续）

6. 应变能密度 $W$ 和余能密度 $W^*$ 的几何意义分别是什么？对线弹性体，$W$ 和 $W^*$ 的关系如何？
7. 写出总势能 $\Pi$ 的完整表达式。每一项的物理含义是什么？
8. 虚功原理的证明中，哪一步使用了"应力张量对称性"？哪一步使用了"Gauss 积分公式"？
9. 虚位移原理和虚应力原理各自等价于什么？为什么它们都不需要本构关系？
10. 试述功的互等定理并给出一维情形的一个应用实例。
11. 最小势能原理中，$\delta\Pi = 0$ 给了我们什么信息？$\delta^2\Pi \geq 0$ 又给了什么信息？
12. 在平面应力问题的总势能变分推导中，本质边界条件和自然边界条件分别是如何出现的？
