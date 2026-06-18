# FEM 理论基础 — Ritz 法、Galerkin 法与加权残量法

> **对应课件**：[`Chapter 4 Applications of variational theory in elasticity (1).pdf`](../06-References/pdfs-originals/Chapter%204%20Applications%20of%20variational%20theory%20in%20elasticity%20(1).pdf) 课程第 2 章 · [原文MD](../../md_output/Chapter%204%20Applications%20of%20variational%20theory%20in%20elasticity%20(1).md)
> **章节定位**：Theoretical Basis of Finite Element → I. Introduction → II. Direct method of variation problem (Ritz法、Galerkin法) → III. Weighted residuals → IV. The finite element method
> **相关作业**：[HW3 Q2（试函数合法性）](../04-Homework-Solutions/2026w/HW3-Problem.md) · [HW3 Q3（弹性地基梁）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 3 章（变分法）、线性代数

> **📋 考试范围覆盖**
>
> | 本讲义章节 | 考试大纲考点 |
> |-----------|-------------|
> | §4.2.2 有限差分法 | [Var. Princ.] Finite difference method |
> | §4.2.3 Ritz 法 | [Var. Princ.] Ritz method |
> | §4.2.7 Galerkin 法 | [Var. Princ.] Galerkin method |
> | §4.3 加权残量法 | [Var. Princ.] Weighted residual method; Equivalent variation equations |

---

## 4.1 引言（Introduction）

### 4.1.1 FEM 是什么？

有限单元法是一种利用计算机分析数学物理问题的**数值近似解法**。它起源于固体力学，并迅速发展到了其他物理领域，如水力学、热传导理论、电磁学和声学。

FEM 的基本思想是将连续体离散为许多有限个单元、然后再将这些单元重新组合成整体——这意味着原连续体的分析可以转化为对每个单元的分析和所有单元的重新组合。

> 🔗 这三种途径对应课程三大理论支柱：结构矩阵法→第5章（直接刚度法），变分法→本章前半（Ritz/Galerkin），加权残量法→本章后半。考试常以"FEM 的数学基础有哪些"形式出现。

### 4.1.2 理解 FEM 的三种途径

我们可以从三个不同的角度来理解有限单元法：

1. **结构矩阵法（Structural Matrix Method）**
2. **变分法（Variational Method）**
3. **加权残量法（Method of Weighted Residuals）**

FEM 的整个计算过程由计算机自动完成，具有非常强的通用性——只需根据不同的工程问题改变输入即可。这种方法彻底改变了分析解的限制。

### 4.1.3 工程设计的任务

工程设计的主要任务是运用固体力学理论（包括结构力学、弹性力学和塑性力学）对结构进行**强度、刚度和稳定性分析**。

随着科学技术的发展，工程结构的几何形状和载荷条件日趋复杂，新材料不断涌现，这使得解析求解变得极为困难甚至不可能。因此人们转向寻求**近似解**。

### 4.1.4 从 Ritz 到 FEM：历史的演进


**1908 年**：Ritz 提出了一种用带未知量的试探函数近似能量泛函的方法——通过对每个未知量求势能泛函的极小值得到求解未知量的方程组。Ritz 法极大地推进了弹性力学在工程中的应用，但其局限在于：**试探函数必须满足边界条件**，这对复杂几何形状的结构来说非常困难。

**1943 年**：Courant 扩展了 Ritz 法——他在求解扭转问题时将整个截面划分为若干三角形区域，在每个三角形区域内假设翘曲函数的近似均匀分布，从而克服了 Ritz 法要求整体近似函数满足所有边界条件的困难。这已非常接近早期有限元的思想，但当时要求的计算量过大，限制了进一步发展。

**1946 年**：计算机诞生后，线性结构系统首先采用了计算机进行数值计算。其理论基础是**矩阵位移法**和**矩阵力法**（统称为结构矩阵法），从结构力学的位移法和力法演化而来。结构矩阵法的力学概念非常清晰，所有理论公式都有非常精确的结构力学含义——舍入误差仅发生在计算机因存储实数位数限制而进行的数值计算过程中。

**1956 年**：Turner、Clough、Martin 和 Topp 在美国航空学会年会上提出了一种从矩阵位移法发展而来的新计算方法——他们将结构划分为三角形和矩形单元，用单元的近似位移函数求解节点力与位移之间的单元刚度矩阵。

与此同时，**Argyris** 发表了多篇能量理论和结构分析的论文，统一了弹性结构的基本能量原理，发展了矩阵方法，并导出了由平面应力板和四个边构件组成的矩形栅格的单元刚度矩阵。

**1960 年**：**Clough** 在论文《The finite element method in plane stress analysis》中首次引入了"有限单元法（FEM）"这一术语。从此 FEM 成为离散化连续体的标准方法。

**1963 年**：在 Besseling、Melosh、Jones 和 Gallagher 等人的工作之后，人们逐渐意识到有限单元法实际上是变分原理中 **Ritz 法的一种形式转化**，并由此发展了从不同变分原理推导有限元计算公式的方法。

**1965 年**：**Zienkiewicz** 和 **谢醒吾（Y.K.Cheung）** 发现，所有场问题只要能够转化为变分形式，就可以用与固体力学有限元法完全相同的步骤求解（如 Laplace 方程和 Poisson 方程），并于 1967 年出版了第一本有限元专著。

> 💡 理解关键：1969 年 Szabo & Lee 的发现极其重要——FEM 并不一定基于变分原理。这意味着即使微分方程没有对应的泛函（变分逆问题无解），也能用加权残量法做 FEM。后续 4.3 节会详述这一点。

**1969 年**：**Szabo** 和 **G.C.Lee** 发现，加权残量法（特别是 Galerkin 法）可以用标准 FEM 程序求解非结构问题——这意味着 **FEM 公式并不一定基于变分原理**。

与此同时，FEM 的数学基础引起了极大关注：大规模线性方程组的数值解法、矩阵特征值方法得到发展，子结构技术和模态综合法得到应用，FEM 的收敛性和误差分析也在研究中。计算机性能不断提升，FEM 获得了越来越多的关注。来自航空、航天、土木工程、造船和汽车等各工业部门的资金支持推动了许多通用有限元程序的开发。

**1980 年代以来**：数百个大型结构有限元通用程序被开发出来，其中最著名包括：
- NASTRAN
- MARC
- SAP—NONSAP
- ADINA
- ANSYS
- ABAQUS

### 4.1.5 中国学者的贡献


中国力学研究者对早期有限元做出了重要贡献：
- **陈伯屏**：结构矩阵法
- **钱令希**：余能原理（complementary potential energy）
- **钱伟长、胡海昌**：广义变分原理（generalized variation）
- **冯康**：有限元理论

> 在固体力学领域用变分原理推导将使概念非常清晰。因此我们先从数学角度解释这些方法，然后介绍它们在力学场中的应用。

---

## 4.2 变分问题的直接法（The Direct Method of Variation Problem）

**为什么学这节？** 这是从"知道变分原理"到"真正会解题"的桥梁。前面学了泛函和 Euler 方程，但 Euler 方程往往求不出解析解。直接法教我们不求解微分方程就能得到近似解——这是 FEM 的核心思想前身。

> 💡 **全局理解**：4.2 节的核心逻辑链是：
> 1. **4.2.1**：微分方程 ↔ 变分问题是等价的（理论基础）
> 2. **4.2.2-4.2.3**：不求解微分方程，直接对泛函求近似极小值（有限差分法、Ritz 法）
> 3. **4.2.4-4.2.6**：把 Ritz 法用在力学问题上（最小势能原理 → 梁的算例）
> 4. **4.2.7-4.2.8**：Galerkin 法（Ritz 法的推广，不要求有泛函）
>
> 简单说：**4.2.1 证明"可以这样做"，4.2.2-4.2.8 教"具体怎么做"**。

### 4.2.1 从微分方程到变分方程

此前我们已经详细讨论了泛函的极值问题，并建立了力学问题的变分方法——但都止步于泛函的 Euler 方程。在处理实际问题时，我们不得不求解 Euler 方程，然而遗憾的是，获得解析解通常非常困难。

因此就出现了**变分问题的直接法**。人们通常将微分问题转化为变分问题，因为变分法更方便。

对力学问题来说，建立相应的泛函并求解其极值问题只是一个求解途径。下面展示另一种途径——**从微分方程到变分方程**。

> ⚠️ 重难点：这个等价性证明是全章的理论基石。核心思路——（⇒）把 V=u₀+W 代入泛函，拆成三项，利用分部积分让中间项消失，余项非负得极小值；（⇐）构造 λ 的二次函数 φ(λ)，利用极小值条件强制 B=0，再用变分法预备定理得微分方程。两个方向证法不同，不要混淆。

微分方程→变分方程的等价性证明（完整）

> 💡 **直觉理解**：这个证明要回答的问题是——"为什么我可以不去解微分方程 Lu=f，而是去最小化某个泛函 Q[u]？"答案是：两者给出同一个解。所以你可以选择更容易的那条路（通常是最小化泛函）。

对于以下常微分方程的边值问题（L 是算子）：

$$\begin{cases}
\mathbf{L}u = -\frac{d}{dx}\left(P(x)\frac{du}{dx}\right) + r(x)u(x) = f(x) \\
u(a) = u(b) = 0
\end{cases}$$
其中 $P(x)$ 和 $r(x)$ 是 Sturm-Liouville 微分算子的已知系数函数（$P(x)$ 与扩散/传导相关，$r(x)$ 与反应/衰减相关），$f(x)$ 是已知源项/强迫函数（表示外部输入或载荷）。

其中 $C^1[a,b]$ 表示区间 $[a,b]$ 上一阶导数连续的函数空间，$C[a,b]$ 表示连续函数空间，后文 $C^2$ 表示二阶导数连续的函数空间。

它等价于以下泛函的极值问题：

$$Q[u(x)] = \int_a^b[-(P(x)u'(x))'u(x) + r(x)u^2(x) - 2f(x)u(x)]dx$$

**我们要证明的命题**：$u_0(x)$ 是微分方程的解 $\Leftrightarrow$ 泛函 $Q[u(x)]$ 在 $u_0(x)$ 取最小值。

> ❌ 易错点："⇒"方向证明中，分部积分后 W(a)=W(b)=0 使边界项消失是关键——这利用了 W∈C₀²[a,b]（在边界上为零）。如果忘记这个条件，推导会卡住。

**证明（1）"$\Rightarrow$"——若 $u_0$ 是微分方程的解，则泛函取极小值**

记 $C_0^2[a,b] = \{f(x) | f(x) \in C^2[a,b], f(a)=f(b)=0\}$。

对任意 $W(x) \in C_0^2[a,b]$，令 $V = u_0 + W$，则 $V \in C_0^2[a,b]$。

代入泛函：

$$\begin{aligned}
Q[V] &= \int_a^b [-(PV')'V + rV^2 - 2fV]dx \\
&= -\left.PV'V\right|_a^b + \int_a^b [P(V')^2]dx + \int_a^b [rV^2 - 2fV]dx
\end{aligned}$$

因为 $V(a)=V(b)=0$，第一项为零。故：

$$Q[V] = \int_a^b [P(V')^2 + rV^2 - 2fV]dx$$

现在代入 $V = u_0 + W$：

$$\begin{aligned}
Q[V] &= Q[u_0 + W] \\
&= \int_a^b \left[P(u_0' + W')^2 + r(u_0 + W)^2 - 2f(u_0 + W)\right]dx \\
&= \int_a^b \left[P u_0'^2 + r u_0^2 - 2f u_0\right]dx \\
&\quad + 2\int_a^b \left[P u_0' W' + r u_0 W - fW\right]dx \\
&\quad + \int_a^b \left[P W'^2 + r W^2\right]dx
\end{aligned}$$

对第二项分部积分（利用 $W(a)=W(b)=0$）：

$$\begin{aligned}
\int_a^b \left[P u_0' W' + r u_0 W\right]dx &= \left.P u_0' W\right|_a^b - \int_a^b (P u_0')'W dx + \int_a^b r u_0 W dx \\
&= \int_a^b [-(P u_0')' + r u_0]W dx \\
&= \int_a^b fW dx \quad (\because u_0 \text{ 满足微分方程})
\end{aligned}$$

因此第二项与第三项抵消：$2\int_a^b fW dx - 2\int_a^b fW dx = 0$。

故：

$$Q[V] = Q[u_0] + \int_a^b [P W'^2 + r W^2]dx$$

由于 $P(x)>0$，$r(x)\geq 0$，后一项 $\int_a^b [P W'^2 + r W^2]dx \geq 0$。

还需证明当 $W(x)\not\equiv 0$ 时该积分**严格大于零**。用反证法：

若积分等于零，则被积函数恒为零：
$$P(x)W'^2(x) + r(x)W^2(x) \equiv 0 \Rightarrow P(x)W'^2(x) \equiv 0 \Rightarrow W'(x) \equiv 0$$

$\Rightarrow W(x) \equiv C$（常数）。又 $W(a)=W(b)=0$，故 $W(x)\equiv 0$，矛盾。

因此当 $V(x)$ 不恒等于 $u_0(x)$ 时：

$$Q[V] = Q[u_0] + \int_a^b [P W'^2 + r W^2]dx > Q[u_0]$$

即 **泛函在 $u_0$ 处取得严格极小值**。$\square$

**证明（2）"$\Leftarrow$"——若泛函取极小值，则 $u_0$ 是微分方程的解**

设泛函 $Q[u]$ 在 $u_0(x)$ 处取最小值。对任意 $W(x) \in C_0^2[a,b]$（$W\not\equiv 0$）和任意实数 $\lambda$：

$$\begin{aligned}
Q[u_0 + \lambda W] &= Q[u_0] + \int_a^b [P(\lambda W')^2 + r(\lambda W)^2]dx \\
&\quad + 2\int_a^b [P u_0'(\lambda W') + r u_0(\lambda W) - f(\lambda W)]dx \\
&\geq Q[u_0]
\end{aligned}$$

将上式视为 $\lambda$ 的二次函数 $\varphi(\lambda)$：

$$\varphi(\lambda) = \lambda^2 \underbrace{\int_a^b [P W'^2 + r W^2]dx}_{A} + 2\lambda \underbrace{\int_a^b [P u_0'W' + r u_0 W - fW]dx}_{B/2} \geq 0$$

即 $A\lambda^2 + B\lambda \geq 0$（$C = 0$，因为 $\varphi(0)=Q[u_0]-Q[u_0]=0$）。

二次三项式恒非负的条件是判别式 $\Delta = B^2 - 4AC \leq 0$。由于 $C=0$，必须 $B=0$：

$$\int_a^b [P u_0' W' + r u_0 W - fW]dx = 0$$

分部积分（利用 $W(a)=W(b)=0$）：

$$\int_a^b [-(P u_0')' + r u_0 - f]W dx = 0$$

由于 $W(x)$ 是任意的，由变分法预备定理：

$$-(P u_0')' + r u_0 - f \equiv 0$$

即 $u_0(x)$ 是微分方程的解。$\square$

> 💡 理解关键：等价的结论是"充要条件"——微分方程的解就是使泛函取极小值的函数，反过来也一样。这意味着我们有两种等价的方式描述同一个物理问题：微分方程形式（强形式）和变分形式（弱形式）。FEM 正是从弱形式出发。

**结论**：$u_0(x)$ 是微分方程的解的**充要条件**是泛函 $Q[u(x)]$ 在 $u_0(x)$ 取最小值——泛函极值问题与微分方程问题是等价的。

---

### 4.2.2 Euler 有限差分法

Euler 在研究变分问题时首先采用了"有限差分法"，这是一种非常重要的直接法。所以 Euler 是变分问题直接法的创始人。

**主要步骤**：

1. 将区间 $[a,b]$ 等分为 $n+1$ 段，分点 $x_0=a, x_1, \ldots, x_n, x_{n+1}=b$
2. 用折线连接点 $(x_i, y_i)$——$y_i$ 除边界外都待定
3. 在子区间 $(x_i, x_{i+1})$ 上用差商代替导数：$y' \approx \frac{y_{i+1}-y_i}{\Delta x}$
4. 泛函 $Q[y(x)]$ 转化为多元函数：
   $$Q[y(x)] \approx \sum_{i=0}^n F\left(x_i, y_i, \frac{y_{i+1}-y_i}{\Delta x}\right)\Delta x_i = \Phi(y_1, \ldots, y_n)$$
5. 由 $\frac{\partial\Phi}{\partial y_1}=0, \frac{\partial\Phi}{\partial y_2}=0, \ldots, \frac{\partial\Phi}{\partial y_n}=0$ 解出 $y_i$

当 $n\to\infty$ 时，我们得到一系列曲线或函数的序列：
$$\{f_1(x), f_2(x), \ldots, f_m(x), \ldots\}$$

n 越大，折线越逼近真实曲线。因为真实解使泛函取极小值，所以 n 越大泛函值越小——该序列称为泛函的**极小化序列**（minimizing sequence）。

> ⚠️ 重难点：极小化序列是连接近似解和精确解的核心概念。通俗理解：每次增加基函数使子空间扩大，得到的泛函值越来越小（因为候选解越来越多），序列最终收敛到真实解。正定算子保证任何极小化序列都收敛——这为 Ritz 法的收敛性提供了数学保证。

> **定义**：满足 $Q[u_n] \to \min Q[u]$（当 $n\to\infty$）的序列 $\{u_n\}$ 称为极小化序列。
>
> **定理**：若微分算子 T 是正定的，则对应泛函的任何极小化序列都收敛到其 Euler 方程的解。即 $\mathbf{T}u = f \Leftrightarrow Q[y(x)] = \int_a^b F(x, y, y')dx$
此处算子 $\mathbf{T}$ 是抽象的微分算子符号，代表从函数空间到函数空间的映射关系（在具体问题中可对应于 $\mathbf{L}$ 或其它微分算子）。

---

### 4.2.3 Ritz 法

**为什么学这节？** Ritz 法是 FEM 的数学源头。1963 年人们意识到 FEM 实际上是 Ritz 法的一种形式转化。理解 Ritz 法 = 理解 FEM 的理论内核。

> 💡 **Ritz 法的直觉**：真实解 u 使泛函 Q[u] 取极小。我们不知道 u 长什么样，但可以猜一个"形状"（用基函数的线性组合），然后调整系数使泛函尽可能小。猜的越准（基函数越多/越好），结果越接近真实解。
>
> **类比**：就像用乐高积木拼一个形状——你只有有限几种积木（基函数），但可以通过调整每种积木的数量（系数）来逼近任何形状。积木种类越多，逼近越精确。

#### 数学框架（不严格）

对于微分方程 $\mathbf{L}u = f$，L 是函数空间 $\mathbf{H}$ 上的微分算子，u 是未知函数，f 是已知函数。$\mathbf{H}$ 是满足一定光滑性和边界条件的函数构成的赋范/内积空间（如 Sobolev 空间或平方可积函数空间）。

定义 L 为：
$$\mathbf{L}u = -\frac{d}{dx}\big(P(x)u'(x)\big) + r(x)u(x)$$

L 是线性的：$\mathbf{L}(c_1 u_1 + c_2 u_2) = c_1 \mathbf{L}u_1 + c_2 \mathbf{L}u_2$。

为了简化书写，记内积为：
$$\int_a^b f(x)g(x)dx = (f,g)$$

对应的泛函可以写为：
$$Q[u(x)] = (\mathbf{L}u, u) - 2(f, u)$$

内积满足双线性性：
$$(f, c_1 g_1 + c_2 g_2) = c_1(f, g_1) + c_2(f, g_2)$$
$$(c_1 f_1 + c_2 f_2, g) = c_1(f_1, g) + c_2(f_2, g)$$


#### Ritz 法的五个步骤

> 💡 **与 FEM 的对应关系**：Ritz 法的五步在 FEM 中都有对应：
> - 第1步（微分→变分）→ FEM 的弱形式推导
> - 第2步（选基函数）→ FEM 的形函数构造
> - 第3步（代入泛函）→ FEM 的单元刚度矩阵计算
> - 第4步（求极值）→ FEM 的总体刚度矩阵组装
> - 第5步（解方程）→ FEM 的求解器
>
> 所以学 Ritz 法就是在学 FEM 的数学原型。

**第 1 步**：通过变分原理将微分方程问题转化为变分问题：
$$\mathbf{T}u = f \Leftrightarrow Q[u] = (\mathbf{T}u, u) - 2(f, u) = \min$$

**第 2 步**：在 H 中选择线性独立的函数 $\{\varphi_1, \varphi_2, \ldots, \varphi_n\}$ 作为基函数（basis functions），张成一个 n 维子空间 $H_n = \text{span}\{\varphi_1, \ldots, \varphi_n\}$。对于其中的任意元素：
$$u_n = \sum_{i=1}^n a_i\varphi_i$$

当 $n=1,2,3,\ldots$ 时，子空间不断扩大，$u_n$ 趋近于泛函的**极小化序列**。

> 💡 理解关键：第3步的操作看似复杂，实质是把泛函的被积函数展开。因为积分是线性的，被积函数中的 u 用 Σaᵢφᵢ 代入后变成 Σaᵢaⱼ 的双重求和——这就是为什么系数矩阵会出现二次型结构。

**第 3 步**：将 $u_n$ 代入泛函的表达式：

$$\begin{aligned}
Q[u] \approx Q[u_n] &= \left(\mathbf{T}\sum_{i=1}^n a_i\varphi_i, \sum_{j=1}^n a_j\varphi_j\right) - 2\left(f, \sum_{i=1}^n a_i\varphi_i\right) \\
&= \sum_{i,j=1}^n a_i a_j(\mathbf{T}\varphi_i, \varphi_j) - 2\sum_{i=1}^n a_i(f, \varphi_i) \\
&= F(a_1, a_2, \ldots, a_n)
\end{aligned}$$

计算内积后，泛函变成了 $a_1, \ldots, a_n$ 的**多元函数**。

**第 4 步**：由多元函数的极值理论，令 $\partial F/\partial a_i = 0$（$i=1,\ldots,n$），得到求解方程：

$$\frac{\partial F}{\partial a_i} = 0 \;\Rightarrow\; \sum_{j=1}^n (\mathbf{T}\varphi_i, \varphi_j)a_j = (f, \varphi_i),\quad i=1,\ldots,n$$

> ❌ 易错点：注意系数矩阵元素是 (Tφᵢ, φⱼ) 而不是 (Tφⱼ, φᵢ)。虽然对称性保证两者相等（分部积分后），但在手算时如果算错了积分顺序，可能得到不对称的结果——这往往是分部积分漏做或边界项忘掉的信号。

写成矩阵形式：

$$\boxed{\begin{bmatrix}
(\mathbf{T}\varphi_1,\varphi_1) & (\mathbf{T}\varphi_1,\varphi_2) & \cdots & (\mathbf{T}\varphi_1,\varphi_n) \\
(\mathbf{T}\varphi_2,\varphi_1) & (\mathbf{T}\varphi_2,\varphi_2) & \cdots & (\mathbf{T}\varphi_2,\varphi_n) \\
\vdots & \vdots & \ddots & \vdots \\
(\mathbf{T}\varphi_n,\varphi_1) & (\mathbf{T}\varphi_n,\varphi_2) & \cdots & (\mathbf{T}\varphi_n,\varphi_n)
\end{bmatrix}
\begin{bmatrix}a_1 \\ a_2 \\ \vdots \\ a_n\end{bmatrix}
= \begin{bmatrix}(f,\varphi_1) \\ (f,\varphi_2) \\ \vdots \\ (f,\varphi_n)\end{bmatrix}}$$

简记为 $\boxed{\mathbf{K}\mathbf{a} = \mathbf{b}}$。

**第 5 步**：解方程组得到 $a_1, \ldots, a_n$，代入 $u_n = \sum a_i\varphi_i$ 得近似解。

#### Ritz 法算例（ODE 边值问题）

**问题**：求满足以下方程的函数 $u(x)\in C^2[0,1]$：

$$\begin{cases} u''(x) - u(x) = x \\ u(0) = u(1) = 0 \end{cases}$$

**精确解**：
$$u(x) = \frac{e}{e^2-1}(e^x - e^{-x}) - x$$

这里 $P(x)\equiv r(x)\equiv 1$，$f(x)\equiv -x$。微分算子是 L 的特例：
$$\mathbf{L}u = -\frac{d}{dx}(P(x)u'(x)) + r(x)u(x) = -\frac{d}{dx}(u') + u = -x$$

对应的泛函为 $Q[u] = (\mathbf{L}u, u) - 2(f, u)$。

**Ritz 法求解**：取基函数（都满足 $\varphi_i(0)=\varphi_i(1)=0$）：

$$\varphi_1(x) = x(1-x),\quad \varphi_2(x) = x^2(1-x),\quad \varphi_3(x) = x^3(1-x)$$

这三个函数线性无关，且 $\varphi_i(x) \in C_0^2[0,1]$。

> 💡 理解关键：这个算例只用 3 项基函数就得到 10.8% 误差——对 ODE 边值问题来说可以接受。但注意系数矩阵的数值：对角线元素（0.3667, 0.1333, 0.0897）递减——高次基函数对应的对角元变小，暗示数值条件可能变差。增加项数提升精度但矩阵病态加重，这是 Ritz 法的内在张力。

代入后得方程组（**注意系数矩阵的对称性**）：

> 💡 **这个方程组是怎么来的？** 以 $k_{11}$ 为例：
> - 先算 $L\varphi_1 = -\varphi_1'' + \varphi_1 = -(-2) + (x-x^2) = 2+x-x^2$
> - 再算 $k_{11} = \int_0^1 (L\varphi_1)\varphi_1\,dx = \int_0^1 (2+x-x^2)(x-x^2)dx$
> - 展开被积函数：$(2+x-x^2)(x-x^2) = 2x - 3x^2 + 2x^3 - x^4$
> - 逐项积分：$2\cdot\frac{1}{2} - 3\cdot\frac{1}{3} + 2\cdot\frac{1}{4} - \frac{1}{5} = 1 - 1 + \frac{1}{2} - \frac{1}{5} = \frac{3}{10} \approx 0.3667$
>
> 右端项 $b_i = \int_0^1 f\cdot\varphi_i\,dx = \int_0^1 (-x)\varphi_i\,dx$：
> - $b_1 = \int_0^1 (-x)(x-x^2)dx = -\frac{1}{3}+\frac{1}{4} = -\frac{1}{12}$
> - $b_2 = \int_0^1 (-x)(x^2-x^3)dx = -\frac{1}{4}+\frac{1}{5} = -\frac{1}{20}$
> - $b_3 = \int_0^1 (-x)(x^3-x^4)dx = -\frac{1}{5}+\frac{1}{6} = -\frac{1}{30}$

$$\begin{cases}
0.3667 a_1 + 0.1833 a_2 - 0.1763 a_3 = -1/12 \\
0.1833 a_1 + 0.1333 a_2 - 0.8940 a_3 = -1/20 \\
-0.1763 a_1 - 0.8940 a_2 + 0.0897 a_3 = -1/30
\end{cases}$$

解得：$a_1 = -0.4425,\; a_2 = 0.4635,\; a_3 = 0.0343$

近似解：
$$\begin{aligned}
u_3(x) &= -0.4425x(1-x) + 0.4635x^2(1-x) + 0.0343x^3(1-x) \\
&= 0.0343x^4 - 0.4292x^3 + 0.9060x^2 - 0.4425x
\end{aligned}$$

在 $x=0.5$ 处检验：
- 精确值：$u(0.5) = -0.0566$
- 近似值：$u_3(0.5) = -0.0505$

误差约 10.8%（仅 3 项基函数）。增加基函数项数可提高精度。


#### 系数矩阵的对称性与正定性（完整证明）

记系数矩阵元素为：
$$k_{ij} = (\mathbf{L}\varphi_i, \varphi_j) = \int_a^b \left[-(P\varphi_i')'\varphi_j + r\varphi_i\varphi_j\right]dx$$

分部积分（利用 $\varphi_j(a)=\varphi_j(b)=0$）：

$$k_{ij} = \int_a^b \left[P\varphi_i'\varphi_j' + r\varphi_i\varphi_j\right]dx = B(\varphi_i, \varphi_j)$$

**对称性**：显然 $k_{ij} = B(\varphi_i, \varphi_j) = B(\varphi_j, \varphi_i) = k_{ji}$。因此系数矩阵 $\mathbf{K}$ 是对称矩阵。

**正定性**：对任意非零向量 $\mathbf{c}^T = (c_1, c_2, \ldots, c_n) \neq \mathbf{0}$，构造 $u_n(x) = \sum_{i=1}^n c_i\varphi_i(x)$（$\varphi_i$ 线性无关，故 $u_n\not\equiv 0$）。

$$\begin{aligned}
\mathbf{c}^T\mathbf{K}\mathbf{c} &= \sum_{i,j=1}^n c_i c_j k_{ij} = \sum_{i,j=1}^n \int_a^b \left[P c_i c_j \varphi_i'\varphi_j' + r c_i c_j \varphi_i\varphi_j\right]dx \\
&= \int_a^b \left[P\left(\sum_{i=1}^n c_i\varphi_i'\right)^2 + r\left(\sum_{i=1}^n c_i\varphi_i\right)^2\right]dx \\
&= B(u_n(x), u_n(x)) = \int_a^b \left[P(u_n')^2 + r u_n^2\right]dx
\end{aligned}$$

> ⚠️ 重难点：正定性证明的关键一步是用 Schwarz 不等式把 u² 的积分用 u'² 的积分控制住。这是 Friedrichs 不等式的一种形式——本质是利用边界条件 u(a)=u(b)=0 将函数值与其导数联系起来。P_min>0 保证整个积分严格正。

利用 Schwarz 不等式证明该积分严格大于零：

$$u(x) = \int_a^x u'(t) \cdot 1 dt \leq \sqrt{\int_a^x [u'(t)]^2 dt}\sqrt{\int_a^x 1 dt} = \sqrt{x-a}\sqrt{\int_a^x [u'(t)]^2 dt}$$

两边平方后积分：

$$\begin{aligned}
\int_a^b u^2(x)dx &\leq \int_a^b \left[(x-a)\int_a^x [u'(t)]^2 dt\right]dx \\
&\leq \int_a^b [u'(t)]^2 dt \int_a^b (x-a)dx = \frac{(b-a)^2}{2}\int_a^b [u'(x)]^2 dx
\end{aligned}$$

因此：
$$\int_a^b [u'(x)]^2 dx \geq \frac{2}{(b-a)^2}\int_a^b u^2(x)dx$$

令 $P_{\min} = \min_{a\leq x\leq b} P(x)$，则：

$$\begin{aligned}
\mathbf{c}^T\mathbf{K}\mathbf{c} &= \int_a^b \left[P(u_n')^2 + r u_n^2\right]dx \\
&\geq \int_a^b \left[P_{\min}(u_n')^2 + r u_n^2\right]dx \\
&\geq \int_a^b \left[P_{\min}\frac{2}{(b-a)^2}u_n^2 + r u_n^2\right]dx \\
&\geq P_{\min}\frac{2}{(b-a)^2}\int_a^b u_n^2 dx
\end{aligned}$$

令 $\delta = P_{\min}\frac{2}{(b-a)^2} > 0$，则 $\mathbf{c}^T\mathbf{K}\mathbf{c} \geq \delta\int_a^b u_n^2 dx > 0$。

**结论**：$\mathbf{K}$ 是正定矩阵，方程组 $\mathbf{K}\mathbf{a} = \mathbf{b}$ 有唯一解，且解使泛函取极小值（因为二阶变分 $\delta^2F = \frac{\partial^2 F}{\partial a_i \partial a_j}\delta a_i \delta a_j = k_{ij}\delta a_i \delta a_j > 0$）。

---

### 4.2.4 Ritz 法在力学中的应用（Ritz 法 vs 最小势能原理）

**为什么学这节？** 前面是纯数学的 Ritz 法框架，这里展示如何直接用在力学问题中——不需要先写微分方程，直接从能量角度建立泛函。

> 💡 **数学 Ritz 法 vs 力学 Ritz 法的区别**：
> | | 数学 Ritz 法（§4.2.3） | 力学 Ritz 法（本节） |
> |---|---|---|
> | **出发点** | 微分方程 $\mathbf{L}u=f$ | 最小势能原理 $\delta\Pi=0$ |
> | **泛函** | $Q[u]=(\mathbf{L}u,u)-2(f,u)$ | $\Pi=U-W$（应变能-外力功） |
> | **基函数要求** | 满足齐次边界条件 | 满足齐次位移边界条件 |
> | **系数矩阵** | $k_{ij}=(\mathbf{L}\varphi_i,\varphi_j)$ | $k_{ij}=\int_\Omega \sigma_{ij}\varepsilon_{ij}\,dV$（刚度矩阵） |
> | **右端项** | $b_i=(f,\varphi_i)$ | $b_i=\int_\Omega f\varphi_i\,dV$（等效节点力） |
>
> **本质相同**：力学 Ritz 法就是数学 Ritz 法在弹性力学中的具体实现。区别只是泛函的形式不同——力学中用能量泛函，数学中用算子泛函。

首先约定边界记号：弹性体的边界分为 $S = S_u \cup S_\sigma$，其中 $S_u$ 是位移已知的边界，$S_\sigma$ 是外力已知的边界。$\bar{u}_i$ 表示 $S_u$ 上给定的已知位移。

#### 力学中的 Ritz 法步骤

> ❌ 易错点：力学中 Ritz 法的试探函数形式是 uᵢ=uᵢ⁰+aᵢₙuᵢₙ，这里 uᵢ⁰ 处理非齐次位移边界（在 S_u 上 uᵢ⁰=ūᵢ），uᵢₙ 在 S_u 上为零。初学者常犯的错误是只用 Σaᵢₙuᵢₙ 而忘了加 uᵢ⁰——这会导致位移边界条件无法满足。非齐次边界条件必须单独处理！

**(1)** 选择位移的**容许试探函数**：
$$u_i = u_i^0 + a_{in}u_{in}, \quad n=1,2,\ldots$$

其中 i 是自由指标（3D 问题中 $i=1,2,3$），$u_{in}$ 是满足齐次位移边界条件的函数（在 $S_u$ 上 $u_{in}=0$），$u_i^0$ 是满足非齐次位移边界条件的函数（在 $S_u$ 上 $u_i^0 = \bar{u}_i$）。

> 💡 **直觉理解**：$u_i^0$ 是"已经满足边界条件的基底"，$a_{in}u_{in}$ 是"在此基础上的修正"。就像搭积木——$u_i^0$ 是底座（保证边界对），$u_{in}$ 是上面的积木（调整内部形状）。

**(2)** 写出弹性系统的总势能 $\Pi(u_i)$ 表达式，代入位移试探函数，得到由待定位移参数 $a_{in}$ 表达的总势能 $\Pi(a_{in})$。

**(3)** 由最小势能原理求变分：
$$\delta\Pi = \frac{\partial\Pi}{\partial a_{in}}\delta a_{in} = 0$$

由于 $\delta a_{in}$ 相互独立，有：
$$\frac{\partial\Pi}{\partial a_{in}} = 0, \quad n=1,2,\ldots$$

这本质上是**由位移参数表示的近似平衡方程**。对于线性弹性问题，$\Pi$ 是位移及其导数的二阶泛函，代入试探函数后得到待定系数的二次函数——上述方程是 $a_{in}$ 的线性代数方程组，易于计算机求解。

**(4)** 解代数方程组得到位移参数 $a_{in}$，代入 $u_i = u_i^0 + a_{in}u_{in}$ 得位移场的近似解。应力和应变也可由位移导出——但得到的应力场一般不满足静力平衡方程。

> ⚠️ **易错点**：应力场不满足平衡方程是 Ritz 法的固有缺陷——因为近似解只使能量泛函极小化，不保证逐点满足微分方程。这是"弱解"的特征：整体能量正确，局部细节有误差。

#### 算例：轴向受拉杆的 Ritz 法

**问题**：一根等截面杆，长度 $L$，截面积 $A$，杨氏模量 $E$。左端固定（$x=0$），右端受集中拉力 $P$（$x=L$）。求杆的位移场 $u(x)$。

```
固定端         自由端
  ┃←——————————→┃
  ┃    P→      ┃
  x=0          x=L
```

**精确解**（材料力学）：$u(x) = \frac{P}{EA}x$（线性分布）

**Ritz 法求解**：

**Step 1**：写出总势能泛函 $\Pi[u]$

- 应变能（杆的拉伸变形能）：$U = \frac{1}{2}\int_0^L EA\left(\frac{du}{dx}\right)^2 dx$
- 外力势能（集中力 $P$ 在位移 $u(L)$ 上做的功的负值）：$W = -P \cdot u(L)$
- 总势能：$\Pi = U + W = \frac{1}{2}\int_0^L EA\left(\frac{du}{dx}\right)^2 dx - Pu(L)$

> 💡 **物理意义**：$\Pi$ 是"存储的弹性能量"减去"外力做的功"。系统会自发地让 $\Pi$ 最小——这等价于力的平衡。

**Step 2**：选择试探函数（基函数）

边界条件：$u(0)=0$（固定端位移为零）。

取一个满足边界条件的试探函数：$u(x) = a_1 \cdot x$

其中 $a_1$ 是待定系数（物理意义：应变 $\varepsilon = du/dx = a_1$）。

> 💡 **为什么选 $x$？** 因为它是最简单的满足 $u(0)=0$ 的函数。当然也可以选 $x^2, x^3$ 等，但越复杂计算量越大。

**Step 3**：代入泛函，得到关于 $a_1$ 的函数

- $\frac{du}{dx} = a_1$（常数）
- $u(L) = a_1 L$

代入总势能：

$$\Pi(a_1) = \frac{1}{2}EA \cdot a_1^2 \cdot L - P \cdot a_1 L = \frac{EAL}{2}a_1^2 - PL a_1$$

> 💡 **关键一步**：泛函 $\Pi[u]$（函数的函数）变成了 $\Pi(a_1)$（数的函数）——这就是"直接法"的精髓：把无穷维问题降为有限维。

**Step 4**：对 $a_1$ 求极值

$$\frac{\partial\Pi}{\partial a_1} = EAL \cdot a_1 - PL = 0$$

解得：$a_1 = \frac{P}{EA}$

**Step 5**：写出近似解

$$u(x) = \frac{P}{EA}x$$

与精确解完全一致！因为真实解恰好是线性函数，而我们选的试探函数也是线性函数。

> 💡 **如果选 $u = a_1 x + a_2 x^2$ 呢？** 会得到 $a_1 = P/(EA)$，$a_2 = 0$——二次项自动消去。这说明试探函数"过度"不会导致错误，但会增加计算量。

> 💡 **关键认识**：如果选 $u = a_1 x + a_2 x^2$（二次函数），会得到相同结果（$a_2=0$）。这说明试探函数"过度"不会导致错误——多余的项会被自动消去（系数为零）。但增加项数会增加计算量，所以要权衡精度和效率。

---

### 4.2.5 梁的变分推导——从势能到边界条件分类

**为什么学这节？** 梁是考试最高频的算例载体。理解如何从势能变分导出平衡方程和边界条件分类，是做 Ritz/Galerkin 法梁算例的前提。

> 💡 **直觉理解**：这一节的核心是展示"能量语言"和"力语言"的等价性。
> - **能量语言**：总势能 Π = 应变能 - 外力功，取极小值
> - **力语言**：平衡方程 + 边界条件
>
> 两者描述同一个物理问题，只是表达方式不同。能量语言的优势是：不需要画受力图，直接写能量表达式就能推导出所有方程和边界条件。

考虑均布荷载 q 作用下的梁（两端有弯矩 $M_0, M_L$ 和剪力 $Q_0, Q_L$），下标 $0$ 表示左端 $(x=0)$，下标 $L$ 表示右端 $(x=L)$。

梁挠度 $w(x)$，小变形下曲率 $\kappa \approx \frac{d^2w}{dx^2}$，应变能：
$$U = \frac12\int_0^L M\kappa\,dx = \frac12\int_0^L EI\left(\frac{d^2w}{dx^2}\right)^2 dx$$
其中 $EI$ 为梁的抗弯刚度（$E$ 为杨氏模量，$I$ 为截面惯性矩）。

外力势能（考虑分布力、端部剪力和弯矩做功）：
$$V = -\int_0^L qw\,dx + Q_0 w_0 - M_0\left(\frac{dw}{dx}\right)_0 - Q_L w_L + M_L\left(\frac{dw}{dx}\right)_L$$

总势能：
$$\Pi = \frac12\int_0^L EI\left(\frac{d^2w}{dx^2}\right)^2 dx - \int_0^L qw\,dx + Q_0 w_0 - M_0\left(\frac{dw}{dx}\right)_0 - Q_L w_L + M_L\left(\frac{dw}{dx}\right)_L$$

求变分 $\delta\Pi$，对第一项进行两次分部积分：

$$\begin{aligned}
\frac12\int_0^L EI\,\delta\left(\frac{d^2w}{dx^2}\right)^2 dx &= \int_0^L EI\frac{d^2w}{dx^2}\delta\left(\frac{d^2w}{dx^2}\right)dx \\
&= \left[EI\frac{d^2w}{dx^2}\delta\left(\frac{dw}{dx}\right)\right]_0^L - \left[\delta w \frac{d}{dx}\left(EI\frac{d^2w}{dx^2}\right)\right]_0^L + \int_0^L \frac{d^2}{dx^2}\left(EI\frac{d^2w}{dx^2}\right)\delta w\,dx
\end{aligned}$$

代入总势能变分：

$$\delta\Pi = \int_0^L \left[\frac{d^2}{dx^2}\left(EI\frac{d^2w}{dx^2}\right) - q\right]\delta w\,dx + \delta\left(\frac{dw}{dx}\right)\left[EI\frac{d^2w}{dx^2} + M\right]\Bigg|_0^L - \delta w\left[\frac{d}{dx}\left(EI\frac{d^2w}{dx^2}\right) + Q\right]\Bigg|_0^L$$

由 $\delta w$ 任意性，得：

**域内平衡方程**：
$$\frac{d^2}{dx^2}\left(EI\frac{d^2w}{dx^2}\right) = q, \quad 0 \leq x \leq L$$

**边界条件**（在 $x=0,L$ 处）：
$$\begin{cases}
\delta\left(\frac{dw}{dx}\right)\left[EI\frac{d^2w}{dx^2} + M\right] = 0 \\
\delta w\left[\frac{d}{dx}\left(EI\frac{d^2w}{dx^2}\right) + Q\right] = 0
\end{cases}$$



根据支承类型区分：

| 支承类型 | 位移条件 | 力边界条件 |
|----------|---------|------------|
| **固定端** | $w=0,\; w'=0$ | 等式自然满足 |
| **简支端** | $w=0$（故 $\delta w=0$） | $EI w'' = -M$；若无端弯矩则 $EI w''=0$ |
| **自由端** | $\delta w\neq 0,\; \delta w' \neq 0$ | $EI w'' = -M,\; \frac{d}{dx}(EI w'') = -Q$；无外力则等于 0 |

**重要概念区分**：
- **本质边界条件**（Essential BC）：与位移本身相关（$w, w'$），Ritz 法的试函数必须满足
- **自然边界条件**（Natural BC）：与力相关（$M, Q$），变分自动满足，不做强制要求

---

### 4.2.6 简支梁的 Ritz 法——试函数选择的影响

**问题**：均布荷载 q 作用下的简支梁，长 l，弯曲刚度 EI，求挠度。

总势能：
$$\Pi = \frac12\int_0^l EI(w'')^2dx - \int_0^l qw\,dx$$

**试函数 1**：$w = ax(l-x)$（满足 $w(0)=w(l)=0$）

代入计算：
- $w'' = -2a$（常数）
- $\frac12\int_0^l EI(-2a)^2 dx = 2EI l a^2$
- $\int_0^l q\,ax(l-x)dx = qa\left[\frac{l^3}{2} - \frac{l^3}{3}\right] = \frac{ql^3}{6}a$

$$\Pi = 2EI l a^2 - \frac{ql^3}{6}a$$

$$\frac{\partial\Pi}{\partial a} = 4EI l a - \frac{ql^3}{6} = 0 \Rightarrow a = \frac{ql^2}{24EI}$$

挠度方程：$w = \frac{ql^2}{24EI}x(l-x)$

**最大挠度**（$x=l/2$ 处）：
$$w_{\max} = \frac{ql^2}{24EI}\cdot\frac{l}{2}\cdot\frac{l}{2} = \frac{ql^4}{96EI}$$

与材料力学精确解 $w_{\max} = \frac{5ql^4}{384EI}$ 相比：
$$\text{误差} = \frac{5/384 - 1/96}{5/384} = \frac{5-4}{5} = 20\%$$

**误差 20%，不太理想。**

---

> 💡 理解关键：为什么 sin 函数比多项式精准这么多？根本原因：sin 是简支梁的真实振型（模态），它天然满足 w''''∝w 的微分方程特征。多项式 w=ax(l-x) 虽然也满足边界条件，但弯成抛物线意味着曲率是常数（w''=-2a）——而真实梁的曲率是变化的。试函数越接近真实解，精度越高——这个直观认识贯穿所有近似方法。

**试函数 2**：$w = a\sin\left(\frac{\pi x}{l}\right)$（也满足 $w(0)=w(l)=0$）

代入计算：
- $w'' = -a\left(\frac{\pi}{l}\right)^2\sin\left(\frac{\pi x}{l}\right)$
- $\frac12\int_0^l EI a^2\left(\frac{\pi}{l}\right)^4\sin^2\left(\frac{\pi x}{l}\right)dx = \frac{EI\pi^4}{4l^3}a^2$
- $\int_0^l q\,a\sin\left(\frac{\pi x}{l}\right)dx = qa\left[-\frac{l}{\pi}\cos\left(\frac{\pi x}{l}\right)\right]_0^l = \frac{2ql}{\pi}a$

$$\Pi = \frac{EI\pi^4}{4l^3}a^2 - \frac{2ql}{\pi}a$$

$$\frac{\partial\Pi}{\partial a} = \frac{EI\pi^4}{2l^3}a - \frac{2ql}{\pi} = 0 \Rightarrow a = \frac{4ql^4}{EI\pi^5}$$

挠度方程：$w = \frac{4ql^4}{EI\pi^5}\sin\left(\frac{\pi x}{l}\right)$

最大挠度（$x=l/2$ 处）：
$$w_{\max} = \frac{4ql^4}{EI\pi^5}\sin\left(\frac{\pi}{2}\right) = \frac{4ql^4}{EI\pi^5}$$

精确解 $w_{\max} = \frac{5ql^4}{384EI} \approx \frac{ql^4}{76.8EI}$

近似解 $w_{\max} = \frac{4ql^4}{\pi^5 EI} \approx \frac{ql^4}{76.43EI}$

$$\text{误差} \approx 0.38\%$$

**误差仅 0.38%，远优于试函数 1 的 20%。**

**核心结论**：试探函数的形式对计算结果影响巨大。$\sin$ 函数是简支梁的实际振型，因此逼近效果好得多。直接法的优势在于只需解代数方程而非微分方程——只要试探函数选择得当，有时仅需一到几个代数方程就能得到足够精确的结果。理论上可取更多项甚至无穷级数来提高精度。

---

### 4.2.7 Galerkin 法

**为什么学这节？** Galerkin 法是加权残量法的起源，也是考试中与 Ritz 法对比频率最高的话题。理解两者对试函数要求的区别即可掌握本节核心。

> 💡 **Galerkin 法 vs Ritz 法的核心区别**：
> - **Ritz 法**：基于最小势能原理，试函数只需满足位移边界条件（本质边界条件）
> - **Galerkin 法**：基于加权残量法，试函数需满足**所有**边界条件（位移 + 力）
>
> **为什么 Galerkin 法要求更严？** 因为 Galerkin 法直接处理微分方程的残量，如果力边界条件不满足，残量在边界上就不为零，加权积分就无法正确消残。
>
> **为什么 Galerkin 法通常更精确？** 因为它隐含了更多信息（力边界条件），所以近似解更接近真实解。

Galerkin 于 1915 年提出了另一种变分问题的直接法，在力学中它是 Ritz 法的一种特殊形式。

最小势能原理的变分为：
$$\delta\Pi = -\int_V (\sigma_{ij,j} + f_i)\delta u_i dV + \int_{S_\sigma}(\sigma_{ij}n_j - \bar{p}_i)\delta u_i dS = 0$$
其中逗号下标表示偏导数（如 $\sigma_{ij,j}\equiv\partial\sigma_{ij}/\partial x_j$，按 Einstein 求和约定对重复指标 $j$ 求和），$f_i$ 是体力分量（注意与标量源项 $f(x)$ 不同，此处下标 $i$ 表示方向），$n_j$ 为边界外法线方向余弦（满足 $n_j = \cos(\hat{n}, x_j)$）。

Ritz 法在此基础上要求试函数满足 $S_u$ 上的位移边界条件。

> ⚠️ 重难点：Galerkin 法对试函数要求更严——必须同时满足位移边界条件 AND 外力边界条件。这是 Ritz 和 Galerkin 最本质的区别。正因要求更严，Galerkin 法逼近效果通常更好（从后续悬臂梁对比算例可以看出）。但这也意味着对复杂边界问题，Galerkin 法的试函数更难构造。

但 Galerkin 法的试探函数**除要求满足 $S_u$ 上的位移边界条件外，还要满足 $S_\sigma$ 上的外力边界条件**。因此第二项面积分自动为零，只剩体积分：

$$\delta\Pi = -\int_V (\sigma_{ij,j} + f_i)\delta u_i dV = 0$$

令 $\delta u_i = u_{in}\delta a_{in}$（$n=1,2,\ldots$），由于 $\delta a_{in}$ 相互独立：

$$\boxed{\int_V (\sigma_{ij,j} + f_i)u_{in}dV = 0,\quad n=1,2,\ldots,N}$$

这就是 **Galerkin 法的求解方程**——积分后也是线性方程组。

#### 梁的 Galerkin 法算例

**问题**：两端固定梁，均布荷载 q，长 l。

取试函数 $w = a(1-\cos\frac{2\pi x}{l})$，显然满足 $w(0)=w(l)=0$ 和 $w'(0)=w'(l)=0$（两端固定）。

此问题无外力边界条件需满足（两端位移条件也是全部的边界条件），故此函数同样适用于 Ritz 法。

Galerkin 法的积分方程为：
$$\int_0^l (EI\frac{d^4w}{dx^4} - q)w\,dx = 0$$

代入 $w = a(1-\cos\frac{2\pi x}{l})$：
$$\frac{d^4w}{dx^4} = a\left(\frac{2\pi}{l}\right)^4\cos\frac{2\pi x}{l}$$

积分：
$$EI a\left(\frac{2\pi}{l}\right)^4 \int_0^l \cos\frac{2\pi x}{l}(1-\cos\frac{2\pi x}{l})dx - q\int_0^l (1-\cos\frac{2\pi x}{l})dx = 0$$

计算：
- $\int_0^l (1-\cos\frac{2\pi x}{l})dx = l$
- $\int_0^l \cos\frac{2\pi x}{l}(1-\cos\frac{2\pi x}{l})dx = \int_0^l (\cos\frac{2\pi x}{l} - \cos^2\frac{2\pi x}{l})dx = 0 - \frac{l}{2} = -\frac{l}{2}$

代入：
$$EI a\left(\frac{2\pi}{l}\right)^4 \left(-\frac{l}{2}\right) - ql = 0$$

$$a = \frac{ql^4}{8\pi^4 EI}$$

挠度：$w = \frac{ql^4}{8\pi^4 EI}(1-\cos\frac{2\pi x}{l})$

Galerkin 法的一个优点：**不需要判断结构是否静定**，因为无论结构静定与否，解法相同。

#### Galerkin 法的加权残量法解释

> 💡 理解关键：Galerkin 法的加权残量解释非常直观：把试函数代入原微分方程会得到残量 R（不为零），我们调整参数使 R 在权函数空间上的投影为零。当权函数取为基函数本身时就是 Galerkin 法。这种"放松后加权消残"的思想是后面五种加权残量法的统一框架。

Galerkin 法的求解方程可以用加权残量法的基本思想解释。以均布荷载梁为例，真实解 w 应满足：

$$EI\frac{d^4w}{dx^4} - q = 0 \Rightarrow \int_0^l (EI\frac{d^4w}{dx^4} - q)w\,dx = 0$$

试探函数 $w_n$ 放松要求后代入平衡方程会在域内产生非零残量 $R$：

$$w_n = \sum_{i=1}^n a_i\varphi_i \Rightarrow EI\frac{d^4w_n}{dx^4} - q = R \neq 0$$

加权残量法通过调整试探函数的待定参数使残量与权函数乘积的积分在整个域内为零：

$$\int_0^l R\,\varphi_i dx = 0,\quad i=1,2,\ldots,N$$

当权函数选为试函数中的各个容许函数 $\varphi_i$ 时，就是 **Galerkin 法**。

---


### 4.2.8 Ritz 法 vs Galerkin 法——悬臂梁对比算例（考试最高频考点）

**为什么这节极重要？** 这是全章对 Ritz 和 Galerkin 法对比最清晰、最完整的一个算例。考试大概率会以某种形式考察两者差异。

> 💡 **这个算例的核心启示**：
> 1. **试函数的选择比方法本身更重要**——同一个方法，选不同的试函数，误差可以从 20% 降到 0.38%
> 2. **Ritz 法偏小，Galerkin 法偏大**——这不是巧合，而是两种方法的内在性质（Ritz 有下限性）
> 3. **边界条件的满足程度决定精度**——Galerkin 法要求满足力边界条件，所以通常更精确
>
> **考试常见问法**：
> - "为什么 Ritz 解偏小？" → 下限性（近似解使系统更刚）
> - "为什么 Galerkin 解更精确？" → 隐含了更多信息（力边界条件）
> - "同一试函数下两者是否等价？" → 是的，当试函数同时满足位移和力边界时

**问题**：均布荷载 q 作用下的悬臂梁（长 l，弯曲刚度 EI），求自由端挠度。

**精确解**（由材料力学）：$w_{\max} = \frac{ql^4}{8EI}$

---

#### 解法 1：Ritz 法

取一阶近似试函数 $w = a\left(1 - \cos\frac{\pi x}{2l}\right)$。

验证边界条件：
- $w(0) = a(1-1) = 0$ ✓（固定端位移为零）
- $w'(0) = a\cdot\frac{\pi}{2l}\cdot 0 = 0$ ✓（固定端转角为零）
- 自由端 $x=l$ 处的外力边界条件不要求满足（Ritz 法不要求）

悬臂梁的总势能（无端部外力）：
$$\Pi = \frac12\int_0^l EI(w'')^2dx - \int_0^l qw\,dx$$

计算：
- $w'' = -a\left(\frac{\pi}{2l}\right)^2\cos\frac{\pi x}{2l}$
- $\int_0^l (w'')^2 dx = a^2\left(\frac{\pi}{2l}\right)^4 \int_0^l \cos^2\frac{\pi x}{2l}dx = a^2\left(\frac{\pi}{2l}\right)^4 \cdot \frac{l}{2}$
- $\int_0^l w\,dx = a\int_0^l (1-\cos\frac{\pi x}{2l})dx = a\left[l - \frac{2l}{\pi}\right] = al(1-\frac{2}{\pi})$

$$\Pi = \frac12 EI a^2\left(\frac{\pi}{2l}\right)^4\frac{l}{2} - qal\left(1-\frac{2}{\pi}\right)$$

$$\frac{\partial\Pi}{\partial a} = 0 \Rightarrow a = \frac{32}{\pi^4}\left(1-\frac{2}{\pi}\right)\frac{ql^4}{EI}$$

自由端 $(x=l)$ 挠度：
$$w_{\max} = a\left(1-\cos\frac{\pi}{2}\right) = a \cdot 1 = \frac{32}{\pi^4}\left(1-\frac{2}{\pi}\right)\frac{ql^4}{EI}$$

数值：$\frac{32}{\pi^4} \approx 0.3285$，$1-\frac{2}{\pi} \approx 0.3634$

$$w_{\max} \approx 0.3285 \times 0.3634 \times \frac{ql^4}{EI} \approx 0.1194 \frac{ql^4}{EI}$$

精确解为 $\frac{ql^4}{8EI} = 0.125 \frac{ql^4}{EI}$

**误差：$\frac{0.125 - 0.1194}{0.125} \approx 4.5\%$（偏小）**

---

#### 解法 2：Galerkin 法

Galerkin 法要求试函数满足**所有边界条件**（包括外力边界）。悬臂梁自由端 $(x=l)$ 处：
$$M(l) = 0 \Rightarrow w''(l) = 0,\quad Q(l) = 0 \Rightarrow w'''(l) = 0$$

> ❌ 易错点：Ritz 法试函数 w=a(1-cos(πx/2l)) 在自由端 w''(l)=0 ✓ 但 w'''(l)≠0 ✗。注意悬臂梁自由端有两个力边界条件——弯矩为零 AND 剪力为零。Galerkin 法要求两者都满足，这是该试函数不适用于 Galerkin 法的根本原因。

Ritz 法的试函数 $w = a(1-\cos\frac{\pi x}{2l})$ 在 $x=l$ 处：
- $w''(l) = -a\left(\frac{\pi}{2l}\right)^2\cos\frac{\pi}{2} = 0$ ✓（弯矩为零满足）
- 但 $w'''(l) \neq 0$（剪力边界不满足）

因此该函数**不适用于 Galerkin 法**。需要构造同时满足外力边界条件的新试函数。

> 💡 理解关键：Galerkin 法的构造策略值得学习——"从高阶导数向低阶积分"。因为力边界条件涉及 w'' 和 w'''，所以从 w'' 开始构造（设 w''=a(1-sin(πx/2l))），保证 w''(l)=w'''(l)=0，再逐次积分回去并用位移边界确定积分常数。这一策略在等参元/高阶单元的形函数构造中也会用到。

**构造策略**：先满足外力边界条件，再通过积分约束位移边界。

设 $w'' = a\left(1 - \sin\frac{\pi x}{2l}\right)$

验证：$w''(l) = a(1-1)=0$，$w'''(l) = -a\frac{\pi}{2l}\cos\frac{\pi}{2}=0$ ✓

积分两次：
$$w' = a\left[x + \frac{2l}{\pi}\cos\frac{\pi x}{2l} + A\right]$$
$$w = a\left[\frac12 x^2 + \left(\frac{2l}{\pi}\right)^2\sin\frac{\pi x}{2l} + Ax + B\right]$$

由固定端条件 $w(0)=0$，$w'(0)=0$ 确定积分常数：
- $w'(0) = a\left[0 + \frac{2l}{\pi}\cdot 1 + A\right] = 0 \Rightarrow A = -\frac{2l}{\pi}$
- $w(0) = a[0 + 0 + 0 + B] = 0 \Rightarrow B = 0$

故试函数为：
$$w = a\left[\frac12 x^2 - \frac{2l}{\pi}x + \left(\frac{2l}{\pi}\right)^2\sin\frac{\pi x}{2l}\right]$$

代入 Galerkin 积分方程 $\int_0^l (EI w'''' - q)w\,dx = 0$（注意此处 $w''''(x) = 0$，因为 $w''$ 不含二次以上分量，但对一阶近似我们直接用 $EI w'''' - q$ 的积分形式）：

经过积分计算得：
$$a = \frac{\frac{1}{6} - \frac{1}{\pi} + \frac{8}{\pi^3}}{\frac{3}{2} - \frac{4}{\pi}} \frac{ql^4}{EI}$$

自由端挠度 $(x=l)$：
$$w_{\max} = a l^2\left(\frac12 - \frac{2}{\pi} + \frac{4}{\pi^2}\right) \approx 0.126 \frac{ql^4}{EI}$$

精确解：$0.125 \frac{ql^4}{EI}$

**误差：$\frac{0.126 - 0.125}{0.125} = 0.8\%$（偏大）**

---

#### 对比总结

| 方法 | 试函数要求 | 自由端挠度 | 误差 |
|------|-----------|-----------|------|
| Ritz | 仅位移 BC | $0.1194\,ql^4/EI$ | 4.5% (偏小) |
| Galerkin | 位移 + 力 BC | $0.126\,ql^4/EI$ | 0.8% (偏大) |
| 精确解 | — | $0.125\,ql^4/EI$ | — |

> ⚠️ 重难点："Ritz 解偏小（4.5% 误差且偏小）"和"Galerkin 解偏大（0.8% 误差且偏大）"不是偶然的。Ritz 法基于最小势能原理，近似解使系统"更刚"（限制了变形自由度），因此位移被低估——这就是位移元的下限性。Galerkin 法不直接最小化能量，不受此约束。

**关键结论**：
1. Galerkin 法对试函数要求更严（需满足所有边界条件），但精度通常更高
2. 同一试函数同时满足位移和力边界时，Galerkin 法和 Ritz 法等价
3. Ritz 法的近似解偏小（满足下限性——Ritz 解给出的是真实位移的下界估计）

---

## 4.3 加权残量法（The Method of Weighted Residuals）

> 🔗 跨章连接：变分逆问题 ─→ 第3章泛函与 Euler 方程的关系。第三章讲的是"正问题"（给定泛函求 Euler 方程），这里讨论"逆问题"（给定微分方程能否找到泛函）。不是所有方程都有泛函——必要条件 φ_{y'} - d/dx(φ_{y''}) = 0 是判定依据。这直接引出了加权残量法的必要性。

### 4.3.1 引言——变分逆问题

**为什么学这节？** 不是所有微分方程都能找到对应的泛函。对于那些没有泛函的方程，Ritz 法根本用不了。加权残量法是更通用的框架——Galerkin 法只是其特例。

我们已经讨论过物理问题或微分方程的近似解法。然而并不是所有微分方程都能转化为变分问题——即需要先考虑变分的**逆问题**。

所谓变分逆问题：给定微分方程 $\phi(x, y, y', y'')=0$，能否找到泛函 $Q[y]=\int_a^b F(x, y, y')dx$ 使该微分方程为泛函的 Euler 方程？

#### 必要条件推导

若 $\phi(x, y, y', y'')=0$ 是泛函 $Q[y]=\int_a^b F(x, y, y')dx$ 的 Euler 方程，则：

$$\phi = F_y - \frac{d}{dx}F_{y'} = F_y - F_{y'x} - y'F_{y'y} - y''F_{y'y'}$$

对其求偏导：

$$\phi_{y'} - \frac{d}{dx}\phi_{y''} = (F_{yy'} - F_{y'yx} - F_{yy'} - y'F_{yy'y'} - y''F_{y'y'y'}) - \frac{d}{dx}F_{y'y'}$$

$$= (-F_{y'yx} - y'F_{yy'y'} - y''F_{y'y'y'}) - (F_{y'y'x} + y'F_{yy'y'} + y''F_{y'y'y'}) = 0$$

因此 **$\phi_{y'} - \frac{d}{dx}\phi_{y''} = 0$** 是 $\phi$ 为 Euler 方程的必要条件。

#### 反例：存在与不存在泛函的对比

**示例 1（有泛函）**：$\phi = -xy + 2xy' + x^2y'' = 0$

验证条件：$\phi_{y'} - \frac{d}{dx}\phi_{y''} = 2x - \frac{d}{dx}(x^2) = 2x - 2x = 0$ ✓

此方程的泛函为 $Q[y] = \frac12\int_a^b (xy^2 + x^2y'^2)dx$


**示例 2（无泛函）**：约去公因子 x 后：$\phi = -y + 2y' + xy'' = 0$

验证条件：$\phi_{y'} - \frac{d}{dx}\phi_{y''} = 2 - \frac{d}{dx}x = 2 - 1 = 1 \neq 0$ ✗

**此方程没有泛函！** 因此必须采用加权残量法。

---

> 💡 理解关键：加权残量法的几何意义是"使残量 R 向权函数空间 W 的投影为零"。想象 R 是一根向量，W 是一个平面——我们调整近似解的参数使 R 垂直于 W。不同的加权残量法本质上是选不同的 W（即选不同的"投影方向"）。

### 4.3.2 加权残量法的基本思想（统一框架）

设在函数空间 H 中求解算子方程 $\mathbf{T}u = f$。

**步骤**：

1. 在 H 中找一组线性独立的基函数 $\varphi_1, \ldots, \varphi_n$，张成**试探函数空间** $M = \text{span}\{\varphi_1, \ldots, \varphi_n\}$
2. 在 H 中找一组线性独立的权函数 $\omega_1, \ldots, \omega_n$，张成**权函数空间** $W = \text{span}\{\omega_1, \ldots, \omega_n\}$
3. 在 M 中设近似解 $u_0 = \sum_{i=1}^n a_i\varphi_i$，$a_i$ 为待定系数
4. 代入方程得**残量**（residual）：$R = \mathbf{T}u_0 - f \neq 0$
5. 加权残量法的核心思想：**使残量 R 向权函数空间 W 的投影为零**（几何意义：在 W 中 $Tu_0$ 和 $Tu$ 有相同的投影）：
   $$\boxed{(R, \omega_i) = (\mathbf{T}u_0 - f, \omega_i) = 0,\quad i=1,2,\ldots,n}$$

   写成积分形式：
   $$\left(\mathbf{T}\sum_{k=1}^n a_k\varphi_k - f,\; \omega_i\right) = 0,\quad i=1,2,\ldots,n$$

解此方程即得近似解。**不同基函数和权函数的选择构成了不同的近似解法**——加权残量法是统一模型。

---


### 4.3.3 五种加权残量法——同题对比算例（考试最高频考点）

**为什么这节极重要？** 五种方法用同一个物理问题（简支梁）做对比，直接给出每种方法的误差数值。考试中常直接要求计算或比较不同方法的精度。

**统一问题**：均布荷载 q 作用下的简支梁（长 l，弯曲刚度 EI），求跨中挠度。

**控制微分方程和边界条件**：
$$EI\frac{d^4w}{dx^4} - q = 0$$
$$x=0: w=0, w''=0;\quad x=l: w=0, w''=0$$

**一阶近似试函数**（满足所有边界条件）：
$$w = c\sin\frac{\pi x}{l}$$

代入平衡方程得内部残量（边界条件已满足，无边界残量）：
$$R = EI\frac{d^4w}{dx^4} - q = EI c\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l} - q$$

**精确解**（材料力学）：$w_{\max} = \frac{5ql^4}{384EI}$

---

#### 方法 1：Galerkin 法（权函数 = 基函数）

$$\int_0^l R \cdot \sin\frac{\pi x}{l}\,dx = \int_0^l \left[EI c\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l} - q\right]\sin\frac{\pi x}{l}\,dx = 0$$

计算：
- $\int_0^l \sin^2\frac{\pi x}{l}dx = \frac{l}{2}$
- $\int_0^l \sin\frac{\pi x}{l}dx = \frac{2l}{\pi}$

$$EI c\left(\frac{\pi}{l}\right)^4 \frac{l}{2} - q\frac{2l}{\pi} = 0 \;\Rightarrow\; c = \frac{4ql^4}{\pi^5 EI}$$

跨中挠度：$w_{\max} = c\sin\frac{\pi}{2} = c = \frac{4ql^4}{\pi^5 EI}$

数值比较：$\frac{4ql^4}{\pi^5 EI} \approx \frac{ql^4}{76.43EI}$，精确解 $\frac{5ql^4}{384EI} \approx \frac{ql^4}{76.8EI}$

**误差：0.386%（考试中最精确的方法之一）**

---

#### 方法 2：最小二乘法（Least Square Method）

原理：使残量平方的积分取极值，即 $\frac{d}{dc}\int_0^l R^2 dx = 0$。

$$\frac{d}{dc}\int_0^l R^2 dx = 2\int_0^l R\frac{dR}{dc}dx = 0$$

$$\frac{dR}{dc} = EI\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l}$$

$$\int_0^l R\frac{dR}{dc}dx = EI c\left(\frac{\pi}{l}\right)^4\frac{l}{2} - q\frac{2l}{\pi} = 0$$

$$\Rightarrow c = \frac{4ql^4}{\pi^5 EI}$$

**各点挠度同 Galerkin 法。误差：0.386%**

*注：在一阶近似且内积定义相同时，Galerkin 法和最小二乘法给出相同结果。*

---

> ❌ 易错点：配点法的精度高度依赖配点位置。本例选跨中 x=l/2 纯属巧合——对简支梁和 sin 试函数来说恰好残量最大处。如果选 x=l/4 或 x=3l/4，误差会更大。配点法在 FEM 中很少单独使用，但作为理解其他方法的参照很有价值。

#### 方法 3：配点法（Collocation Method）

原理：在求解域内选若干配点，令残量在这些点处为零。配点数 = 待定参数个数。

一阶近似只有一个参数，选跨中 $x=l/2$ 为配点：

$$\left.R\right|_{x=l/2} = EI c\left(\frac{\pi}{l}\right)^4\sin\frac{\pi}{2} - q = EI c\left(\frac{\pi}{l}\right)^4 - q = 0$$

$$\Rightarrow c = \frac{ql^4}{\pi^4 EI}$$

跨中挠度：$w_{\max} = \frac{ql^4}{\pi^4 EI} \approx \frac{ql^4}{97.41EI}$

**误差：$\frac{1/97.41 - 1/76.8}{1/76.8} \approx 21.16\%$**

---

#### 方法 4：矩法（Moment Method）

原理：权函数为坐标的幂次 $1, x, x^2, \ldots$。一阶近似取权函数 $\omega(x) = 1$。

$$\int_0^l R \cdot 1\,dx = \int_0^l \left[EI c\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l} - q\right]dx = 0$$

$$\int_0^l \sin\frac{\pi x}{l}dx = \frac{2l}{\pi},\quad \int_0^l dx = l$$

$$EI c\left(\frac{\pi}{l}\right)^4 \frac{2l}{\pi} - ql = 0 \;\Rightarrow\; c = \frac{ql^4}{2\pi^3 EI}$$

跨中挠度：$w_{\max} = \frac{ql^4}{2\pi^3 EI} \approx \frac{ql^4}{62.01EI}$

**误差：$\frac{1/62.01 - 1/76.8}{1/76.8} \approx 23.85\%$**

---

> ⚠️ 重难点：一阶近似的子域法（取整个域为子域）与矩法（权函数 ω(x)=1）等价。但多阶近似时它们不同：子域法将域分成多个子区间分别积分；矩法用 1, x, x²,... 作为权函数。不要在多参数情况下混淆两者。

#### 方法 5：子域法（Subdomain Method）

原理：将求解域分为若干子域，令残量在每个子域上的积分为零。子域数 = 参数个数。

一阶近似只需一个子域——即整个求解域 $[0,l]$，此时子域法与矩法等价。

$$\int_0^l R\,dx = 0 \;\Rightarrow\; c = \frac{ql^4}{2\pi^3 EI}$$

**误差：23.85%（同矩法）**

---


#### 五种方法对比总表

| 方法 | 权函数 $\omega_i$ | 跨中挠度 $w_{\max}$ | 误差 | 精度排名 |
|------|-------------------|---------------------|------|---------|
| **Galerkin 法** | $\varphi_i$（基函数本身） | $\frac{4ql^4}{\pi^5 EI}$ | 0.386% | **1** |
| **最小二乘法** | $\partial R/\partial c_i$ | $\frac{4ql^4}{\pi^5 EI}$ | 0.386% | **1** |
| 配点法 | $\delta(x-x_i)$ | $\frac{ql^4}{\pi^4 EI}$ | 21.16% | 3 |
| 矩法 | $1, x, x^2, \ldots$ | $\frac{ql^4}{2\pi^3 EI}$ | 23.85% | 4 |
| 子域法 | 子域上为 1 | $\frac{ql^4}{2\pi^3 EI}$ | 23.85% | 4 |

**核心结论**：
1. Galerkin 法和最小二乘法精度最高（本题 0.386%），配点法次之（21.16%），矩法/子域法最差（23.85%）
2. 所有五种方法可统一表达为加权残量格式，区别仅在于权函数的选择
3. 对有泛函的问题，Galerkin 法与最小二乘法等价（在一阶近似时）
4. Galerkin 法系数矩阵对称（重要计算优势）

---

> 🔗 跨章连接：方差泛函与第3章的变分原理一脉相承。能量泛函要求方程有对应的泛函（即变分逆问题有解），方差泛函则无此限制——直接对残量平方积分求极小。这解释了为什么加权残量法比 Ritz 法适用范围更广：方差泛函为任何微分方程提供了"人造"变分原理。

### 4.3.4 方差泛函理论（Functional of Variance）

**为什么学这节？** 方差泛函为加权残量法提供了严格的数学基础——它说明了"最小化残量"这个直觉操作的合理性。

#### 问题的一般形式

工程中的控制微分方程和边界条件可统一写为：

$$\mathbf{P}_m(u) = f_0 \quad (\text{域内})$$
$$\mathbf{Q}_t(u) = f_t \quad (t=1,2,\ldots,m) \quad (\text{边界})$$

其中 $\mathbf{P}_m$ 是 m 阶微分算子，$\mathbf{Q}_t$ 是边界上法向 $(t-1)$ 阶微分算子，$f_t$ 是给定函数。$\Omega$ 的边界 $\Gamma = \bigcup \Gamma_t$。

问题假定是**适定的**（well-posed）：解的存在性、唯一性、对定解条件的连续依赖性。

#### 方差泛函的定义

设试函数 $\tilde{u} = \sum c_i\varphi_i$。将 $\tilde{u}$ 代入方程和边界条件，产生：

- **内部残量** $R_I = \mathbf{P}_m(\tilde{u}) - f_0$
- **边界残量** $R_t = \mathbf{Q}_t(\tilde{u}) - f_t$（$t=1,2,\ldots,m$）

**方差泛函**（Functional of Variance）定义为所有残量平方的加权积分：

$$J_R(\tilde{u}) = \int_\Omega \beta_I^2 R_I^2 d\Omega + \sum_{t=1}^m \int_{\Gamma_t} \beta_t^2 R_t^2 d\Gamma$$

> 💡 理解关键：方差泛函的两个权系数 β_I 和 β_t 有双重功能——(1) 量纲统一：内部残量和边界残量的量纲通常不同，直接相加无意义；(2) 重要性加权：想强调某个边界条件就增大对应的 β。这是工程直觉转化为数学形式的典范。

其中 $\beta_I$ 和 $\beta_t$ 为**权系数**，有两个作用：
1. 增大 $\beta_t$ 以加强重要边界条件的满足程度
2. 量纲调整：使内部残量和边界残量的量纲一致，加法才有意义

#### 方差泛函的两个重要性质

1. **最小值为零**：$J_R(\tilde{u}) \geq 0$，当且仅当 $\tilde{u}$ 为精确解时取零
2. **等价性**：使 $J_R$ 取最小值的解等价于原微分方程边值问题的解


#### 方差泛函 vs 能量泛函

| 对比维度 | 能量泛函 | 方差泛函 |
|---------|---------|---------|
| 极小值 | 未知（须求解后才知） | 确定为零 |
| 二次项的物理意义 | 有明确物理意义（能量） | 无明确物理意义 |
| 适用条件 | 仅适用于有泛函的问题 | **普适**（有无泛函均可） |
| 极值条件 | $\delta\Pi = 0$（一阶变分为零） | $\delta J_R = 0$（一阶变分为零） |
| 共同点 | 线弹性问题均为正定二次型 | 线弹性问题均为正定二次型 |

#### 加权残量法的几何解释

重写加权残量方程为内积形式：

$$(\mathbf{T}u_0 - f, \omega_i) = (\mathbf{T}u_0 - \mathbf{T}u, \omega_i) = 0$$
$$\Rightarrow (\mathbf{T}u_0, \omega_i) = (\mathbf{T}u, \omega_i)$$

**几何意义**：近似解 $Tu_0$ 和精确解 $Tu$ 在权函数空间 W 上的投影相等。

当子空间 M 的维数 $n\to\infty$ 时，方差泛函的极小值进一步降低，近似解趋近于精确解——**只有当试函数空间扩展到所有容许函数时，方差泛函的极小值才为零**。

不同加权残量法本质上是构造极小化序列的不同方式。

**重要提示**：对有泛函的问题，Galerkin 法（直接变分法）和 Galerkin 法（加权残量法）得到的求解方程完全相同，系数矩阵对称——这对求解至关重要。

---

## 4.4 有限单元法（The Finite Element Method）

**为什么学这节？** 前面 Ritz 法和 Galerkin 法都要求在整个求解域上构造试函数——复杂几何下这几乎不可能。FEM 的核心创新是"分片插值"：在每个小单元上用简单多项式逼近，然后拼起来。

> ⚠️ 重难点：2D Poisson 方程 Galerkin 法算例的启示——即使是正方形域这样简单的几何，构造两个满足所有边界条件的试函数就已经积分出极其繁琐的系数。这正是"分片插值"思想的必要性：把大问题切成小问题，在每个小块上用简单函数逼近。

### 4.4.1 为什么需要 FEM —— 2D Poisson 方程的 Galerkin 法算例

以一维问题为例引入各种方法后，看一个二维问题理解 FEM 的必要性。

**问题**：用 Galerkin 法解 Poisson 方程

$$\begin{cases}
\Omega: \Delta u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = -2 \\
\Gamma: u = 0
\end{cases}$$

其中 $\Omega : \{|x| < a, |y| < a\}$（正方形域），$\Gamma$ 是 $\Omega$ 的边界。

取两个线性无关的基函数（满足所有边界条件）：
$$\phi_1 = (a^2 - x^2)(a^2 - y^2)$$
$$\phi_2 = (a^2 - x^2)(a^2 - y^2)(x^2 + y^2)$$

设 $u_2 = c_1\phi_1 + c_2\phi_2$。

Galerkin 求解方程：
$$\begin{cases}
\iint_\Omega \phi_1(\Delta u_2 + 2)d\Omega = 0 \\
\iint_\Omega \phi_2(\Delta u_2 + 2)d\Omega = 0
\end{cases}$$

积分后：
- $(\mathbf{T}\phi_1, \phi_1) = -\frac{16\times 16}{5\times 3\times 3}a^8$
- $(\mathbf{T}\phi_1, \phi_2) = -\frac{12\times 16\times 16}{35\times 5\times 3\times 3}a^{10}$
- $(\mathbf{T}\phi_2, \phi_1) = -\frac{8\times 128}{35\times 5\times 3}a^{10}$
- $(\mathbf{T}\phi_2, \phi_2) = -\frac{22\times 32\times 8\times 2}{7\times 5\times 5\times 27}a^{12}$
- $(f, \phi_1) = -\frac{32}{9}a^6$
- $(f, \phi_2) = -\frac{16\times 4}{5\times 9}a^6$

解得：$c_1 = \frac{1}{a^2}\frac{5\times 7\times 37}{8\times 277},\quad c_2 = \frac{1}{a^4}\frac{15\times 35}{16\times 277}$

近似解：
$$u \approx u_2 = \frac{1}{a^2}\frac{35}{16\times 277}(a^2-x^2)(a^2-y^2)\left[74 + \frac{15}{a^2}(x^2+y^2)\right]$$

**启示**：即使对于如此简单的正方形域，构造满足所有边界条件的试函数就已经非常复杂。当边界曲线形状复杂时，这几乎是不可能的任务。此外，Galerkin 法要求试函数有二阶可积的导数——这促使我们：

1. 将求解方程转化为等效积分弱形式（降低对试函数光滑度的要求）
2. 对求解域进行分区插值（扩大试函数的适用范围）

这就是有限单元法的由来。

---

### 4.4.2 1D FEM 的数学描述

在区间 $[a,b]$ 上任意选取离散点：
$$a = x_0 < x_1 < x_2 < \cdots < x_i < x_{i+1} < \cdots < x_n = b$$

$x_i$ 称为**节点**，区间 $[x_i, x_{i+1}]$ 称为**单元**。

对每个节点 $x_i$ 建立线性基函数：

$$\varphi_i(x) = \begin{cases}
1 & x = x_i \\
\text{线性} & x\in(x_{i-1}, x_{i+1}) \\
0 & \text{其他}
\end{cases}$$

> 💡 理解关键：分片线性基函数的 Kronecker 性质（φᵢ(xⱼ)=δᵢⱼ）使待定系数 aᵢ 直接就是节点位移值 uᵢ——这是 FEM 区别于 Ritz 法的一大特色。Ritz 法中 aᵢ 是抽象的广义坐标，没有直观物理意义。此外，每个基函数的影响域仅限于相邻两个单元，导致总体刚度矩阵稀疏——这是 FEM 的计算引擎。

基函数的**Kronecker 性质**：每个基函数在对应节点处值为 1，在其他节点处值为 0。

因此近似解的系数 $a_1, \ldots, a_n$ 就是 $v(x)$ 在离散点处的**函数值** $u_i$。

$$\nu(x) = \sum_{i=0}^n u_i\varphi_i(x)$$

基函数的**影响域（支撑集）非常小**（只有相邻两个单元），因此所得线性方程组的系数矩阵非常**稀疏**——这是 FEM 的计算优势。

**FEM 的核心思想**：用 $u_n$ 作为无限维函数空间中真实解 u 在子空间 M 中的**最佳逼近**——我们实际上在构造泛函的极小化序列。

---


### 4.4.3 FEM 完整 7 步流程

以一维方程为例：

$$\begin{cases}
-(Pu')' + ru = f \\
u(a) = u_a,\quad u(b) = u_b
\end{cases}$$

> 🔗 跨章连接：等效积分弱形式 → 第3章变分原理。弱形式的本质是"降阶"——通过分部积分将对 u 的二阶导数要求降为一阶（从要求 C¹ 连续降到 C⁰ 连续），这使分片线性插值成为可能。如果保持强形式，FEM 就无法使用简单的线性单元。

#### 第 1 步：建立等效积分弱形式

泛函的等效积分弱形式（只要求 $C^0$ 连续性）：
$$Q[u(x)] = \frac12\int_a^b (P u'^2 + r u^2 - 2f u)dx$$

#### 第 2 步：离散化求解域

将区间 $[a,b]$ 划分为 n 个单元：
$$a = x_0 < x_1 < \cdots < x_n = b$$

域内泛函离散为各单元泛函之和：
$$Q[u] = \sum_{j=1}^n Q_j^e$$

#### 第 3 步：构造分片基函数

采用分段线性多项式 $\phi_i$ 作为基函数：

$$\phi_0(x) = \begin{cases} (x_1-x)/(x_1-x_0) & x\in[x_0,x_1] \\ 0 & x\notin[x_0,x_1] \end{cases}$$

$$\phi_i(x) = \begin{cases} (x-x_{i-1})/(x_i-x_{i-1}) & x\in[x_{i-1},x_i] \\ (x_{i+1}-x)/(x_{i+1}-x_i) & x\in[x_i,x_{i+1}] \\ 0 & x\notin[x_{i-1},x_{i+1}] \end{cases}$$

$$\phi_n(x) = \begin{cases} (x-x_{n-1})/(x_n-x_{n-1}) & x\in[x_{n-1},x_n] \\ 0 & x\notin[x_{n-1},x_n] \end{cases}$$

> ⚠️ 重难点：形函数与基函数的关系容易混淆。形函数（N₁, N₂）是单元局部概念——每个单元有自己的形函数。基函数 φᵢ 是全局概念——由相邻单元的形函数拼接而成。φᵢ 在节点 i 的两个相邻单元上有支撑（除非是端点），其余位置为零。

用单元局部编号定义**形函数**（shape functions）。对于单元 $e_j: [x_{j-1}, x_j]$：

$$N_1^{(j)}(x) = \frac{x_j - x}{x_j - x_{j-1}}, \quad N_2^{(j)}(x) = \frac{x - x_{j-1}}{x_j - x_{j-1}}$$

上标为单元号，下标为局部节点号（1 或 2）。

基函数用形函数表示：
$$\phi_0 = N_1^{(1)}$$
$$\phi_i = \begin{cases} N_2^{(i)} & x\in e_i \\ N_1^{(i+1)} & x\in e_{i+1} \\ 0 & \text{其他} \end{cases}$$
$$\phi_n = N_2^{(n)}$$

#### 第 4 步：单元分析

在单元 $e_j$ 内，**只有 $\phi_{j-1}$ 和 $\phi_j$ 非零**。试函数在 $e_j$ 中：
$$v(x) = u_{j-1}N_1^{(j)}(x) + u_j N_2^{(j)}(x)$$

代入单元泛函：
$$Q_j = \frac12\int_{e_j} (P v'^2 + r v^2 - 2f v)dx$$

展开并积分，得到单元刚度矩阵 $K^{ej}$ 和单元载荷向量 $B^{ej}$：

记 $k_{\alpha\beta}^{ej} = \int_{e_j}\left(P\frac{dN_\alpha^{(j)}}{dx}\frac{dN_\beta^{(j)}}{dx} + rN_\alpha^{(j)}N_\beta^{(j)}\right)dx$，$b_\alpha^{ej} = \int_{e_j} fN_\alpha^{(j)}dx$

$$Q_j = \frac12 u_{j-1}^2 k_{11}^{ej} + u_{j-1}u_j k_{12}^{ej} + \frac12 u_j^2 k_{22}^{ej} - u_{j-1}b_1^{ej} - u_j b_2^{ej}$$

矩阵形式：
$$\begin{pmatrix} \partial Q_j/\partial u_{j-1} \\ \partial Q_j/\partial u_j \end{pmatrix} = \underbrace{\begin{pmatrix} k_{11}^{ej} & k_{12}^{ej} \\ k_{12}^{ej} & k_{22}^{ej} \end{pmatrix}}_{K^{ej}} \begin{pmatrix} u_{j-1} \\ u_j \end{pmatrix} - \underbrace{\begin{pmatrix} b_1^{ej} \\ b_2^{ej} \end{pmatrix}}_{B^{ej}}$$

$K^{ej}$ 称为单元 $e_j$ 的**单元刚度矩阵**，对称。

> 💡 理解关键：单元刚度矩阵 Kᵉʲ 总是 2×2 对称矩阵（对两端节点杆单元）。整体组装的本质是"对号入座"——将各单元 2×2 的 Kᵉʲ 按节点全局编号叠加到 n×n 的总体 K 中。在共享节点处（如节点 i 同时属于单元 i 和 i+1），对应刚度元素叠加。

#### 第 5 步：总体组装（Assembly）

将所有单元对泛函的贡献按公共节点叠加，得到总体极值方程：

$$\frac{\partial Q}{\partial u_i} = \sum_{j=1}^n \frac{\partial Q_j}{\partial u_i} = 0 \quad (i=0,1,\ldots,n)$$

由于每个 $u_i$ 只出现在 $Q_i$ 和 $Q_{i+1}$ 中，所以：
$$\frac{\partial Q}{\partial u_i} = \frac{\partial Q_i}{\partial u_i} + \frac{\partial Q_{i+1}}{\partial u_i}$$

组装得总体方程 $\mathbf{K}\mathbf{U} = \mathbf{B}$，总体刚度矩阵 $\mathbf{K}$ 具有带状稀疏结构。

#### 第 6 步：引入边界条件

- **自然边界条件**：在变分泛函中自动满足，无需特殊处理
> ❌ 易错点：对角线置1法 vs 划行划列——两者数学等价但编程便利性不同。对角线置1法保持矩阵规模不变（节点编号不变），适合编程；划行划列改变矩阵维度，手算更直观。考试中手算建议用划行划列，但要注明使用了哪种方法。注意：对角线置1法中同行同列的非对角元素必须清零，否则方程不对。

- **本质边界条件**（给定位移）：用**"对角线置 1 法"**处理

若第 $i$ 个节点的值已知为 $u_i = \bar{u}$，则将 $\mathbf{K}$ 的第 i 行第 i 列对角线元素置 1，其余元素置 0，右端项改为 $\bar{u}$。

#### 第 7 步：求解方程

解线性代数方程组得各节点的近似解 $u_1, u_2, \ldots, u_{n-1}$（$u_0, u_n$ 为已知边界值）。

全域近似解：
$$v(x) = \sum_{i=0}^n u_i\phi_i(x)$$

---


### 4.4.4 自重杆 FEM 完整数值算例

**为什么学这节？** 这是杆单元 FEM 最经典的力学算例——每一步的数字都完整给出，是理解"从物理问题到 FEM 数值解"全流程的最佳素材。

**问题**：一端固定、等截面均质弹性杆，长 L，截面积 A，密度 $\rho$，在自重和自由端拉力 P 作用下保持平衡。求位移分布和应力分布。

建立坐标系（$x$ 轴沿杆向下），设 $u(x)$ 为截面 x 的位移。则：
$$\varepsilon(x) = \frac{du}{dx}, \quad \sigma(x) = E\varepsilon(x)$$

**等效积分弱形式**（由虚功原理直接建立）：

设虚位移为 $\phi(x)$，则虚应变为 $\tilde{\varepsilon} = \phi'$。虚功原理：**内力虚功 = 外力虚功**。

$$\int_0^L EA u'(x)\phi'(x)dx = \int_0^L A\rho g \phi(x)dx + P\phi(L)$$

上式即为等效积分弱形式。

---

#### 步骤 1：单元划分（离散化）

将 $[0,L]$ 划分为 n 段：
$$0 = x_0 < x_1 < \cdots < x_n = L$$

单元 $e_i: [x_{i-1}, x_i]$，节点 $x_i$（$i=0,1,\ldots,n$）。

---

#### 步骤 2：位移模式（形函数）

采用线性位移模式。在单元 $e_i$ 中，已知两端位移 $u_{i-1}$ 和 $u_i$：

$$u(x) = \frac{x_i - x}{x_i - x_{i-1}}u_{i-1} + \frac{x - x_{i-1}}{x_i - x_{i-1}}u_i, \quad x_{i-1} \leq x \leq x_i$$

形函数：
$$N_1^{(i)}(x) = \frac{x_i - x}{L_i}, \quad N_2^{(i)}(x) = \frac{x - x_{i-1}}{L_i}$$

其中 $L_i = x_i - x_{i-1}$ 为单元长度。

形函数满足（插值性质）：
- $N_1(x_{i-1}) = 1, \quad N_1(x_i) = 0$
- $N_2(x_{i-1}) = 0, \quad N_2(x_i) = 1$

矩阵形式：
$$u(x) = \begin{bmatrix}N_1^{(i)} & N_2^{(i)}\end{bmatrix} \begin{Bmatrix} u_{i-1} \\ u_i \end{Bmatrix} = [N]\{\delta\}_{ei}$$

---

#### 步骤 3：单元应变和 B 矩阵

$$\varepsilon_i = \frac{du}{dx} = \begin{bmatrix}\frac{dN_1^{(i)}}{dx} & \frac{dN_2^{(i)}}{dx}\end{bmatrix} \begin{Bmatrix} u_{i-1} \\ u_i \end{Bmatrix} = [B]\{\delta\}_{ei}$$

$$B = \begin{bmatrix} -\frac{1}{L_i} & \frac{1}{L_i} \end{bmatrix}$$

---

#### 步骤 4：单元分析 — 刚度矩阵

将位移模式和虚位移代入虚功方程，对每个单元积分。内力虚功部分：

$$\begin{aligned}
EA\int_{ei} u'(x)\phi'(x)dx &= EA\int_{ei} \left([B]\{\delta^*\}_{ei}\right)^T \left([B]\{\delta\}_{ei}\right)dx \\
&= \{\delta^*\}_{ei}^T \left(\int_{ei} [B]^T EA [B] dx\right) \{\delta\}_{ei} \\
&= \{\delta^*\}_{ei}^T [k]_{ei} \{\delta\}_{ei}
\end{aligned}$$

> ❌ 易错点：杆单元刚度矩阵 Kᵉ=(EA/L)[1 -1; -1 1] 是最简单的单元刚度矩阵。容易出错的地方：(1) 混淆局部编号和全局编号——Kᵉ 的 (1,1) 对应单元左节点，在全局中可能是任意编号；(2) 忘记乘 EA/L——刚度矩阵必须有正确的物理量纲（N/m）。阶梯轴算例中特别注意 A 的变化对刚度的影响。

**单元刚度矩阵**：
$$[k]_{ei} = \int_{ei} [B]^T EA [B] dx = EA L_i \begin{bmatrix} -\frac{1}{L_i} \\ \frac{1}{L_i} \end{bmatrix} \begin{bmatrix} -\frac{1}{L_i} & \frac{1}{L_i} \end{bmatrix} = \frac{EA}{L_i} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

外力虚功部分：
$$A\int_{ei} \phi(x)f dx = A\int_{ei} ([N]\{\delta^*\}_{ei})^T f dx = \{\delta^*\}_{ei}^T \underbrace{\int_{ei} A[N]f dx}_{\{F\}_{ei}}$$

**单元等效节点力**（对均布体力 $f = \rho g$）：
$$\{F\}_{ei} = A\rho g \int_0^{L_i} \begin{Bmatrix} N_1 \\ N_2 \end{Bmatrix} dx = \frac{A\rho g L_i}{2} \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}$$

---

#### 步骤 5：总体组装

将各单元的刚度矩阵按整体节点编号组装。总体刚度矩阵的性质：
1. **对称**
2. **稀疏**
3. **半正定**（引入边界条件后正定）

---

#### 步骤 6：边界条件处理（对角线置 1 法）

固定端 $x=0$ 处 $u_0 = 0$。将 $\mathbf{K}$ 第一行第一列对角线置 1，同行同列其余元素置 0，右端 $B_1$ 置 0。

---

#### 步骤 7：求解位移

解 $\mathbf{K}\mathbf{U} = \mathbf{B}$ 得各节点位移 $u_1, u_2, \ldots, u_n$。

---

#### 步骤 8：应力回算

单元应变：$\varepsilon_i = [B]\{\delta\}_{ei}$，单元应力：$\sigma_i = E\varepsilon_i$

还可以由节点力平衡反算支座反力。

---


### 4.4.5 阶梯轴 FEM 完整数值算例

**为什么学这节？** 这是最贴近工程实际的 FEM 算例——两段不同截面的阶梯轴，用具体数值完整演示从刚度矩阵推导到应力计算的全过程。

**问题**：如图所示的阶梯轴，右端受轴向拉力 P。

已知数据：
- $A^{(1)} = 2.0 \times 10^{-6} \text{ m}^2$
- $A^{(2)} = 1.0 \times 10^{-6} \text{ m}^2$
- $E^{(1)} = E^{(2)} = E = 2.0 \times 10^{11} \text{ Pa}$
- $P = 100 \text{ N}$
- $L^{(1)} = L^{(2)} = 1.0 \text{ m}$

---

#### (1) 结构离散化

结构只有轴向变形，每个节点一个轴向位移自由度。

- 节点编号：1, 2, 3
- 单元编号：(1), (2)
- 节点自由度：$u_1, u_2, u_3$
- 节点载荷：$P_1, P_2, P_3$

已知量：$u_1 = 0$（固定端），$P_2 = 0$，$P_3 = P = 100\text{N}$
未知量：$P_1$（支座反力），$u_2, u_3$

---

#### (2) 位移插值函数

在任意单元内设 $u(x) = a + bx$。

由节点位移约束：
- 在节点 i：$u(0) = a = u_1^e$
- 在节点 j：$u(L^e) = a + bL^e = u_2^e$

解得 $b = (u_2^e - u_1^e)/L^e$

代入得：
$$u(x) = \left(1 - \frac{x}{L^e}\right)u_1^e + \frac{x}{L^e}u_2^e$$

矩阵形式：
$$\{u\} = \underbrace{\begin{bmatrix}1-\frac{x}{L^e} & \frac{x}{L^e}\end{bmatrix}}_{[N]} \begin{Bmatrix} u_1^e \\ u_2^e \end{Bmatrix}$$

---

#### (3) 单元分析——推导刚度矩阵

单元弹性应变能：
$$U^e = \frac12 \iiint_{V^e} \sigma\varepsilon\,dv = \frac12 A^e\int_0^{L^e} \sigma\varepsilon\,dx$$

应变：
$$\varepsilon = \frac{\partial u}{\partial x} = \frac{\partial}{\partial x}[N]\{u^e\} = \begin{bmatrix}-\frac{1}{L^e} & \frac{1}{L^e}\end{bmatrix} \begin{Bmatrix}u_1^e \\ u_2^e\end{Bmatrix} = [B]\{u^e\}$$

应力：
$$\sigma = E\varepsilon = E[B]\{u^e\}$$

代入应变能：
$$U^e = \frac12 A^e \int_0^{L^e} \{u^e\}^T [B]^T E[B] \{u^e\} dx = \frac12 \{u^e\}^T \underbrace{\left(A^e L^e E[B]^T[B]\right)}_{[K^e]} \{u^e\}$$

**单元刚度矩阵**：
$$\begin{aligned}
[K^e] &= A^e L^e E \begin{bmatrix} -\frac{1}{L^e} \\ \frac{1}{L^e} \end{bmatrix} \begin{bmatrix} -\frac{1}{L^e} & \frac{1}{L^e} \end{bmatrix} \\
&= \frac{A^e E}{L^e} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}
\end{aligned}$$

---

**单元 (1) 具体数值**：

$$\{u^{(1)}\} = \begin{Bmatrix} u_1 \\ u_2 \end{Bmatrix}$$

$$[K^{(1)}] = \frac{A^{(1)}E}{L^{(1)}} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix} = \frac{2.0\times 10^{-6} \times 2.0\times 10^{11}}{1.0} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

$$[K^{(1)}] = \begin{bmatrix} 4.0 \times 10^5 & -4.0 \times 10^5 \\ -4.0 \times 10^5 & 4.0 \times 10^5 \end{bmatrix}$$

---

**单元 (2) 具体数值**：

$$\{u^{(2)}\} = \begin{Bmatrix} u_2 \\ u_3 \end{Bmatrix}$$

$$[K^{(2)}] = \frac{A^{(2)}E}{L^{(2)}} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix} = \frac{1.0\times 10^{-6} \times 2.0\times 10^{11}}{1.0} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

$$[K^{(2)}] = \begin{bmatrix} 2.0 \times 10^5 & -2.0 \times 10^5 \\ -2.0 \times 10^5 & 2.0 \times 10^5 \end{bmatrix}$$

---

#### (4) 总体组装

总应变能 $U = U^{(1)} + U^{(2)}$：

$$\begin{aligned}
U &= \frac12 \{u_1, u_2\} \begin{bmatrix} 4.0\times 10^5 & -4.0\times 10^5 \\ -4.0\times 10^5 & 4.0\times 10^5 \end{bmatrix} \begin{Bmatrix}u_1 \\ u_2\end{Bmatrix} \\
&\quad + \frac12 \{u_2, u_3\} \begin{bmatrix} 2.0\times 10^5 & -2.0\times 10^5 \\ -2.0\times 10^5 & 2.0\times 10^5 \end{bmatrix} \begin{Bmatrix}u_2 \\ u_3\end{Bmatrix}
\end{aligned}$$

扩展为 3x3 矩阵后叠加：

$$U = \frac12 \begin{Bmatrix}u_1 \\ u_2 \\ u_3\end{Bmatrix}^T \underbrace{\begin{bmatrix}
4.0\times 10^5 & -4.0\times 10^5 & 0 \\
-4.0\times 10^5 & 6.0\times 10^5 & -2.0\times 10^5 \\
0 & -2.0\times 10^5 & 2.0\times 10^5
\end{bmatrix}}_{[K]} \begin{Bmatrix}u_1 \\ u_2 \\ u_3\end{Bmatrix}$$

> 💡 理解关键：组装规律的核心——k₂₂ = k₂₂⁽¹⁾ + k₁₁⁽²⁾ = 4×10⁵ + 2×10⁵ = 6×10⁵。这是因为节点 2 同时属于单元 (1) 的右节点（对应 K⁽¹⁾ 的 (2,2) 元素）和单元 (2) 的左节点（对应 K⁽²⁾ 的 (1,1) 元素）。每个节点处刚度叠加的物理含义是"并联弹簧"——两个相邻单元的刚度在共享节点处并联。

组装规律：总体刚度矩阵中，$k_{22} = k_{22}^{(1)} + k_{11}^{(2)} = 4.0\times 10^5 + 2.0\times 10^5 = 6.0\times 10^5$

外力虚功：$W_p = \{u_1, u_2, u_3\} \begin{Bmatrix}P_1 \\ 0 \\ 100\end{Bmatrix}$

总势能：
$$\Pi = \frac12 \{u\}^T[K]\{u\} - \{u\}^T\{P\}$$

由 $\partial\Pi/\partial\{u\} = \{0\}$：
$$\begin{bmatrix}
4.0\times 10^5 & -4.0\times 10^5 & 0 \\
-4.0\times 10^5 & 6.0\times 10^5 & -2.0\times 10^5 \\
0 & -2.0\times 10^5 & 2.0\times 10^5
\end{bmatrix}
\begin{Bmatrix}u_1 \\ u_2 \\ u_3\end{Bmatrix}
= \begin{Bmatrix}P_1 \\ 0 \\ 100\end{Bmatrix}$$

---

#### (5) 边界条件处理

已知 $u_1 = 0$。**方法 1（划行划列）**：删除第一行第一列：

$$\begin{bmatrix} 6.0\times 10^5 & -2.0\times 10^5 \\ -2.0\times 10^5 & 2.0\times 10^5 \end{bmatrix} \begin{Bmatrix}u_2 \\ u_3\end{Bmatrix} = \begin{Bmatrix}0 \\ 100\end{Bmatrix}$$

**方法 2（对角线置 1）**（更适合计算机编程）：

$$\begin{bmatrix}
1 & 0 & 0 \\
0 & 6.0\times 10^5 & -2.0\times 10^5 \\
0 & -2.0\times 10^5 & 2.0\times 10^5
\end{bmatrix}
\begin{Bmatrix}u_1 \\ u_2 \\ u_3\end{Bmatrix}
= \begin{Bmatrix}0 \\ 0 \\ 100\end{Bmatrix}$$

---

#### (6) 求解方程

$$\begin{cases}
6.0\times 10^5 u_2 - 2.0\times 10^5 u_3 = 0 \\
-2.0\times 10^5 u_2 + 2.0\times 10^5 u_3 = 100
\end{cases}$$

由第一式：$6u_2 = 2u_3 \Rightarrow u_3 = 3u_2$

代入第二式：$-2.0\times 10^5 u_2 + 2.0\times 10^5 \cdot 3u_2 = 100$

$$4.0\times 10^5 u_2 = 100 \Rightarrow u_2 = 0.25 \times 10^{-3} \text{ m}$$

$$u_3 = 0.75 \times 10^{-3} \text{ m}$$

---


#### (7) 应力回算

**单元 (1) 应力**：

$$\begin{aligned}
\varepsilon^{(1)} &= [B]\{u^{(1)}\} = \begin{bmatrix}-1/L^{(1)} & 1/L^{(1)}\end{bmatrix} \begin{Bmatrix}0 \\ 0.25\times 10^{-3}\end{Bmatrix} \\
&= [-1, 1] \begin{Bmatrix}0 \\ 0.25\times 10^{-3}\end{Bmatrix} = 0.25 \times 10^{-3}
\end{aligned}$$

$$\sigma^{(1)} = E \varepsilon^{(1)} = 2.0\times 10^{11} \times 0.25\times 10^{-3} = 0.5\times 10^8 \text{ Pa} = 50 \text{ MPa}$$

**单元 (2) 应力**：

$$\begin{aligned}
\varepsilon^{(2)} &= [B]\{u^{(2)}\} = \begin{bmatrix}-1/L^{(2)} & 1/L^{(2)}\end{bmatrix} \begin{Bmatrix}0.25\times 10^{-3} \\ 0.75\times 10^{-3}\end{Bmatrix} \\
&= [-1, 1] \begin{Bmatrix}0.25\times 10^{-3} \\ 0.75\times 10^{-3}\end{Bmatrix} = 0.5 \times 10^{-3}
\end{aligned}$$

$$\sigma^{(2)} = E \varepsilon^{(2)} = 2.0\times 10^{11} \times 0.5\times 10^{-3} = 1.0\times 10^8 \text{ Pa} = 100 \text{ MPa}$$

**结果验证**：
- 单元 (1) 应力 50 MPa：对应载荷 $P = \sigma^{(1)} A^{(1)} = 0.5\times 10^8 \times 2.0\times 10^{-6} = 100\text{ N}$ ✓
- 单元 (2) 应力 100 MPa：截面小一半，应力大一倍 ✓

---

#### 完整算例数值汇总

| 量 | 值 |
|----|-----|
| $A^{(1)}$ | $2.0\times 10^{-6} \text{ m}^2$ |
| $A^{(2)}$ | $1.0\times 10^{-6} \text{ m}^2$ |
| $E$ | $2.0\times 10^{11} \text{ Pa}$ |
| $P$ | 100 N |
| $[K^{(1)}]$ | $\begin{bmatrix}4.0\times 10^5 & -4.0\times 10^5 \\ -4.0\times 10^5 & 4.0\times 10^5\end{bmatrix}$ |
| $[K^{(2)}]$ | $\begin{bmatrix}2.0\times 10^5 & -2.0\times 10^5 \\ -2.0\times 10^5 & 2.0\times 10^5\end{bmatrix}$ |
| 总体 $[K]$ | $\begin{bmatrix}4.0\times 10^5 & -4.0\times 10^5 & 0 \\ -4.0\times 10^5 & 6.0\times 10^5 & -2.0\times 10^5 \\ 0 & -2.0\times 10^5 & 2.0\times 10^5\end{bmatrix}$ |
| $u_2$ | $0.25\times 10^{-3} \text{ m}$ |
| $u_3$ | $0.75\times 10^{-3} \text{ m}$ |
| $\sigma^{(1)}$ | $0.5\times 10^8 \text{ Pa}$ (50 MPa) |
| $\sigma^{(2)}$ | $1.0\times 10^8 \text{ Pa}$ (100 MPa) |

---

## 4.5 结论（Conclusion）

> ⚠️ 重难点：FEM 的本质——"经典变分法 + 分片多项式插值"。这短短一句话浓缩了整个第4章的核心逻辑。变分法提供了数学框架（泛函极值→代数方程），分片插值解决了复杂几何的试函数构造难题。两者结合才诞生了 FEM。理解这点就理解了为什么 FEM 既严谨又通用。

从数学角度看，有限单元法是**经典变分法与分片多项式插值的结合**：

1. 利用变分法的离散化和灵活性，在变分原理基础上构造有限维子空间
2. 将无限维函数空间中的变分问题转化为有限元子空间中的多元函数极值问题——多元函数的未知量是待求函数在离散节点处的值
3. 根据多元函数极值理论求解控制方程（线性代数方程组），代数方程组的解即为待求函数在离散节点处的数值解

**FEM 的两大驱动力**：
- **计算机的出现** — 计算工具
- **泛函分析的应用** — 理论基础

自 1960 年代 FEM 进入理论研究高潮以来，越来越多的数学家致力于将 FEM 从工程问题的约束中解放出来，结合泛函分析、变分法、差分法、函数逼近论等学科，最终建立起系统的数学基础。

> 推荐参考书：
> 1. 《有限元法的数学基础》王烈衡、许学军 编著
> 2. 《有限元方法及其理论基础》姜礼尚、庞之垣 著
> 3. 《有限元方法及其应用》李开泰、黄艾香、黄庆怀 编著

---

## 检查你的理解

### 基础概念

1. Ritz 法和 Galerkin 法对试探函数的要求有何不同？
2. Galerkin 法与加权残量法是什么关系？
3. 极小化序列（minimizing sequence）的含义是什么？
4. 什么是本质边界条件和自然边界条件？在变分框架中各如何处理？
5. 什么条件下微分方程没有对应的泛函？举出一个反例。

### 计算与推导

6. 为什么简支梁算例中 $w=a\sin(\pi x/l)$ 比 $w=ax(l-x)$ 精度高得多？
7. Ritz 法中系数矩阵对称性和正定性是如何证明的？正定性保证了解的唯一性——请复述证明思路。
8. 悬臂梁算例中，为什么同一个试函数适用于 Ritz 法却不适用于 Galerkin 法？
9. 五种加权残量法以同一个简支梁问题作对比——请写出每种方法的误差百分比，并解释为什么 Galerkin 法和最小二乘法最精确。
10. Galerkin 法（作为加权残量法）与 Galerkin 法（作为直接变分法）在什么条件下等价？

### FEM 流程

11. FEM 的标准 7 步流程是什么？每步的核心操作是什么？
12. 单元刚度矩阵是如何从形函数和 B 矩阵导出的？总体刚度矩阵的组装规律是什么？
13. 边界条件的"对角线置 1 法"是如何操作的？为什么它比划行划列更适合计算机编程？
14. 阶梯轴算例中，总体刚度矩阵的 $k_{22}$ 为什么等于 $6.0\times 10^5$？哪两个单元刚度矩阵的哪些元素贡献了它？

### 深入思考

15. 方差泛函与能量泛函有何异同？为什么说方差泛函比能量泛函适用范围更广？
16. 为什么 FEM 的总体刚度矩阵是稀疏的？这与基函数的什么性质有关？

---

> **对应作业**：[HW3 Q2（试函数合法性）](../04-Homework-Solutions/2026w/HW3-Problem.md) · [HW3 Q3（弹性地基梁）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **往年参考**：[past/HW3/Homework3](../04-Homework-Solutions/past/HW3/Homework3.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW3/Ans%20to%20HM3_LIU%20Sai_handed%20in.md)
