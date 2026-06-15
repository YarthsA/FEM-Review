# FEM 理论基础 — Ritz 法、Galerkin 法与加权残量法

> **对应课件**：[`Chapter 4 Applications of variational theory in elasticity (1).pdf`](../06-References/pdfs-originals/Chapter%204%20Applications%20of%20variational%20theory%20in%20elasticity%20(1).pdf) 课程第 2 章
> **章节定位**：Theoretical Basis of Finite Element → I. Introduction → II. Direct method of variation problem (Ritz法、Galerkin法) → III. Weighted residuals → IV. The finite element method
> **相关作业**：[HW3 Q2（试函数合法性）](../04-Homework-Solutions/2026w/HW3-Problem.md) · [HW3 Q3（弹性地基梁）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 3 章（变分法）、线性代数

---

## 4.1 引言（Introduction）

### 4.1.1 FEM 是什么？

有限单元法是一种利用计算机分析数学物理问题的**数值近似解法**。它起源于固体力学，并迅速发展到了其他物理领域，如水力学、热传导理论、电磁学和声学。

FEM 的基本思想是将连续体离散为许多有限个单元、然后再将这些单元重新组合成整体——这意味着原连续体的分析可以转化为对每个单元的分析和所有单元的重新组合。

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

**1946 年**：计算机诞生后，线性结构系统首先采用了计算机进行数值计算。其理论基础是**矩阵位移法**和**矩阵力法**（统称为结构矩阵法），从结构力学的位移法和力法演化而来。采用矩阵代数表达方程更为简洁，也便于计算机编程实现。

**1956 年**：Turner、Clough、Martin 和 Topp 在美国航空学会年会上提出了一种从矩阵位移法发展而来的新计算方法——他们将结构划分为三角形和矩形单元，用单元的近似位移函数求解节点力与位移之间的单元刚度矩阵。

与此同时，**Argyris** 发表了多篇能量理论和结构分析的论文，统一了弹性结构的基本能量原理，发展了矩阵方法，并导出了由平面应力板和四个边构件组成的矩形栅格的单元刚度矩阵。

**1960 年**：**Clough** 在论文《The finite element method in plane stress analysis》中首次引入了"有限单元法（FEM）"这一术语。从此 FEM 成为离散化连续体的标准方法。

**1963 年**：在 Besseling、Melosh、Jones 和 Gallagher 等人的工作之后，人们逐渐意识到有限单元法实际上是变分原理中 **Ritz 法的一种形式转化**，并由此发展了从不同变分原理推导有限元计算公式的方法。

**1965 年**：**Zienkiewicz** 和 **谢醒吾（Y.K.Cheung）** 发现，所有场问题只要能够转化为变分形式，就可以用与固体力学有限元法完全相同的步骤求解（如 Laplace 方程和 Poisson 方程），并于 1967 年出版了第一本有限元专著。

**1969 年**：**Szabo** 和 **G.C.Lee** 发现，加权残量法（特别是 Galerkin 法）可以用标准 FEM 程序求解非结构问题——这意味着 FEM 公式并不一定基于变分原理。

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

### 4.2.1 引言

此前我们已经详细讨论了泛函的极值问题，并建立了力学问题的变分方法——但都止步于泛函的 Euler 方程。在处理实际问题时，我们不得不求解 Euler 方程，然而遗憾的是，获得解析解通常非常困难。

因此就出现了**变分问题的直接法**。人们通常将微分问题转化为变分问题，因为变分法更方便。

对力学问题来说，建立相应的泛函并求解其极值问题只是一个求解途径。下面展示另一种途径——**从微分方程到变分方程**。

#### 从微分方程到变分方程

对于以下常微分方程的边值问题（L 是算子）：
$$\begin{cases}
\mathbf{L}u = -\frac{d}{dx}(P(x)\frac{du}{dx}) + r(x)u(x) = f(x) \\
u(a) = u(b) = 0
\end{cases}$$

其中 $P(x)\in C^1[a,b],\; r(x)\in C[a,b]$，$P(x)>0,\; r(x)\geq 0$ 且 $r(x)$ 不恒为零。

它等价于以下泛函的极值问题：
$$Q[u(x)] = \int_a^b[-(P(x)u'(x))'u(x) + r(x)u^2(x) - 2f(x)u(x)]dx$$

**证明其等价性的思路**：
1. 若 $u_0(x)$ 是微分方程的解，则它使泛函取极小值
2. 若泛函在 $u_0(x)$ 处取极小值，则 $u_0(x)$ 是微分方程的解

**证明（1）**：记 $C_0^2[a,b] = \{f(x)|f\in C^2[a,b], f(a)=f(b)=0\}$。对任意 $W\in C_0^2$，令 $V=u_0+W$。代入泛函并分部积分，由于 $V(a)=V(b)=0$，可得：
$$Q[V] = Q[u_0] + \int_a^b[PW'^2 + rW^2]dx \geq Q[u_0]$$

当 $W\not\equiv0$ 时，因 $P(x)>0, r(x)\geq0$，积分项为正，故 $Q[V] > Q[u_0]$。

**证明（2）**：若 $Q[u]$ 在 $u_0$ 处取最小值，对任意 $W$ 和实数 $\lambda$，二次型 $\varphi(\lambda)=Q[u_0+\lambda W] \geq Q[u_0]$。展开并令判别式 $\Delta = B^2 - 4AC \leq 0$，得 $B=0$。分部积分后由变分法预备定理导出微分方程。

**结论**：$u_0(x)$ 是微分方程的解的**充要条件**是泛函 $Q[u(x)]$ 在 $u_0(x)$ 取最小值——泛函极值问题与微分方程问题是等价的。

#### Euler 有限差分法

Euler 在研究变分问题时首先采用了"有限差分法"，这是一种非常重要的直接法。所以 Euler 是变分问题直接法的创始人。

主要步骤：
1. 将区间 $[a,b]$ 等分为 $n+1$ 段，分点 $x_0=a, x_1, \ldots, x_n, x_{n+1}=b$
2. 用折线连接点 $(x_i, y_i)$，$y_i$ 待定
3. 用差商代替导数：$y' \approx (y_{i+1}-y_i)/\Delta x$
4. 泛函 $Q[y]$ 转化为多元函数 $\Phi(y_1,\ldots,y_n)$
5. 由 $\partial\Phi/\partial y_i = 0$ 解出 $y_i$

当 $n\to\infty$ 时，我们得到一系列曲线或函数的序列：
$$\{f_1(x), f_2(x), \ldots, f_m(x), \ldots\}$$

n 越大，折线越逼近真实曲线。因为真实解使泛函取极小值，所以 $n$ 越大泛函值越小——该序列称为泛函的**极小化序列**。

> **定理**：若微分算子 T 是正定的，则对应泛函的任何极小化序列都收敛到其 Euler 方程的解。

---

### 4.2.2 Ritz 法

#### 回顾

从以上讨论可以确认微分问题可以转化为变分问题。现在用数学的方式（虽然不严格）来解释。

对于微分方程 $\mathbf{L}u = f$，L 是函数空间 H 上的微分算子，u 是未知函数，f 是已知函数。

为了简化书写，记 $\int_a^b f(x)g(x)dx = (f,g)$，称为函数 f 和 g 的**内积**。

对应的泛函可以写为：
$$Q[u(x)] = (\mathbf{L}u, u) - 2(f, u)$$

#### Ritz 法的步骤

**(1)** 通过变分原理将微分方程问题转化为变分问题：
$$\mathbf{T}u = f \Leftrightarrow Q[u] = (\mathbf{T}u, u) - 2(f, u) = \min$$

**(2)** 在 H 中选择线性独立的函数 $\{\varphi_1, \varphi_2, \ldots, \varphi_n\}$ 作为基函数，张成一个 n 维子空间 $H_n = \text{span}\{\varphi_1, \ldots, \varphi_n\}$。对于其中的任意元素：
$$u_n = \sum_{i=1}^n a_i\varphi_i$$

当 $n=1,2,3,\ldots$ 时，子空间不断扩大，$u_n$ 趋近于泛函的**极小化序列**。

**(3)** 将 $u_n$ 代入泛函的表达式：
$$Q[u] \approx Q[u_n] = \sum_{i,j=1}^n a_i a_j(\mathbf{T}\varphi_i, \varphi_j) - 2\sum_{i=1}^n a_i(f, \varphi_i) = F(a_1, \ldots, a_n)$$

计算内积后，泛函变成了 $a_1, \ldots, a_n$ 的**多元函数**。

**(4)** 由多元函数的极值理论，得到求解方程：
$$\frac{\partial F}{\partial a_i} = 0 \quad\Rightarrow\quad \sum_{j=1}^n (\mathbf{T}\varphi_i, \varphi_j)a_j = (f, \varphi_i),\quad i=1,\ldots,n$$

写成矩阵形式：
$$\begin{pmatrix}
(\mathbf{T}\varphi_1,\varphi_1) & (\mathbf{T}\varphi_1,\varphi_2) & \cdots \\
(\mathbf{T}\varphi_2,\varphi_1) & (\mathbf{T}\varphi_2,\varphi_2) & \cdots \\
\vdots & \vdots & \ddots
\end{pmatrix}
\begin{pmatrix}a_1 \\ a_2 \\ \vdots\end{pmatrix}
= \begin{pmatrix}(f,\varphi_1) \\ (f,\varphi_2) \\ \vdots\end{pmatrix}$$

**(5)** 解方程组得到 $a_1, \ldots, a_n$，代入 $u_n = \sum a_i\varphi_i$ 得近似解。

#### 算例

**问题**：求满足 $u''(x) - u(x) = x$，$u(0)=u(1)=0$ 的函数 $u(x)\in C^2[0,1]$。

**精确解**：
$$u(x) = \frac{e}{e^2-1}(e^x - e^{-x}) - x$$

这里 $P(x)\equiv r(x)\equiv 1$，$f(x)\equiv -x$。对应的泛函为 $Q[u] = (\mathbf{L}u, u) - 2(f, u)$。

**Ritz 法**：取基函数：
$$\varphi_1(x) = x(1-x),\; \varphi_2(x) = x^2(1-x),\; \varphi_3(x) = x^3(1-x)$$

它们显然满足 $\varphi_i(0) = \varphi_i(1) = 0$，且线性独立。

代入后得方程组（**注意系数矩阵的对称性**）：
$$\begin{cases}
0.3667a_1 + 0.1833a_2 - 0.1763a_3 = -1/12 \\
0.1833a_1 + 0.1333a_2 - 0.8940a_3 = -1/20 \\
-0.1763a_1 - 0.8940a_2 + 0.0897a_3 = -1/30
\end{cases}$$

解得 $a_1=-0.4425, a_2=0.4635, a_3=0.0343$。

近似解：
$$u_3(x) = -0.4425x(1-x) + 0.4635x^2(1-x) + 0.0343x^3(1-x)$$

在 $x=0.5$ 处：$u_{\text{exact}}(0.5)=-0.0566$，$u_3(0.5)=-0.0505$

误差明显（因仅取 3 项基函数），但增加基函数项数可提高精度。

#### 系数矩阵的进一步讨论

对称性：
$$k_{ij} = (\mathbf{L}\varphi_i,\varphi_j) = \int_a^b[P\varphi_i'\varphi_j' + r\varphi_i\varphi_j]dx = k_{ji}$$

正定性：对任意 $\mathbf{c}\neq\mathbf{0}$，
$$\mathbf{c}^T\mathbf{K}\mathbf{c} = B(u_n, u_n) = \int_a^b[P(u_n')^2 + r u_n^2]dx > 0$$

这确保了方程 $\mathbf{K}\mathbf{a} = \mathbf{b}$ 有唯一解，且解使泛函取极小值。

#### Ritz 法的力学应用

Ritz 法直接使用最小势能原理：
1. 选择位移的**容许试探函数**
2. 写出弹性系统的总势能 $\Pi(u_i)$ 表达式，代入位移试探函数
3. 由 $\partial\Pi/\partial a_{in} = 0$ 得到求解方程——本质上是由位移参数表示的**近似平衡方程**

#### 简支梁的 Ritz 法

**问题**：均布荷载 q 作用下的简支梁，求挠度。

总势能：$\Pi = \frac12\int_0^L EI(w'')^2dx - \int_0^L qw\,dx$

**试函数 1**：$w = ax(l-x)$（满足 $w(0)=w(l)=0$）

代入得 $\Pi = 2EIla^2 - \frac16ql^3a$，$\partial\Pi/\partial a = 0 \Rightarrow a = \frac{ql^2}{24EI}$

$$w_{\max} = \frac{ql^4}{96EI}$$

与材料力学精确解 $5ql^4/(384EI)$ 相比，**误差约 20%** → 不太理想。

**试函数 2**：$w = a\sin(\pi x/l)$（也满足边界条件）

代入得 $\Pi = \frac{EI\pi^4}{4l^3}a^2 - \frac{2ql}{\pi}a$

$a = \frac{4ql^4}{EI\pi^5}$，$w_{\max} = \frac{4ql^4}{EI\pi^5}$，**误差仅 0.38%**。

**结论**：试探函数的形式对计算结果影响巨大。采用更多项甚至无穷级数可以得到更精确的结果。直接法的优势在于只需解代数方程而非微分方程——只要试探函数选择得当，有时仅需少数几个（甚至一个）代数方程就能得到足够精确的结果。

---

### 4.2.3 Galerkin 法

Galerkin 于 1915 年提出了另一种变分问题的直接法，在力学中它是 Ritz 法的一种特殊形式。

变分方程 $\delta\Pi = 0$ 是 Ritz 法的基础——建立满足 $S_u$ 上位移边界条件的试探函数并代入。

但 Galerkin 法的试探函数**除要求满足 $S_u$ 上的位移边界条件外，还要满足 $S_\sigma$ 上的外力边界条件**。因此势能变分可以写为：
$$\delta\Pi = -\int_V (\sigma_{ij,j}+f_i)\delta u_i dV = 0$$

由 $\delta u_i = u_{in}\delta a_{in}$ 且 $\delta a_{in}$ 相互独立，得：
$$\int_V (\sigma_{ij,j}+f_i)u_{in}dV = 0,\quad n=1,2,\ldots,N$$

这就是 Galerkin 法的求解方程——积分后也是线性方程组。

#### 梁的 Galerkin 法示例

**问题**：两端固定梁，均布荷载 q。

取试函数 $w = a(1-\cos\frac{2\pi x}{l})$，显然满足 $w(0)=w(l)=0$，$w'(0)=w'(l)=0$。

Galerkin 法的积分方程为：
$$\int_0^l (EI\frac{d^4w}{dx^4} - q)w\,dx = 0$$

代入计算得 $a = \frac{ql^4}{8\pi^4 EI}$。

Galerkin 法的一个优点：**不需要判断结构是否静定**，因为无论结构静定与否，解算结果一样。

#### Galerkin 法的加权残量法解释

Galerkin 法的求解方程可以用加权残量法的基本思想解释。以均布荷载梁为例，真实解 w 应满足：
$$EI\frac{d^4w}{dx^4} - q = 0$$

试探函数 $w_n$ 放松要求后代入平衡方程会在域内产生非零残量 $R_i$。加权残量法通过调整试探函数的待定参数使残量与权函数乘积的积分在整个域内为零，从而获得近似解：
$$\int_0^l R\,\varphi_i dx = 0,\quad i=1,2,\ldots,N$$

当权函数选为试函数中的各个容许函数 $\varphi_i$ 时，就是 **Galerkin 法**。

---

## 4.3 加权残量法（The Method of Weighted Residuals）

### 4.3.1 引言

我们已经讨论过物理问题或微分方程的近似解法。然而并不是所有微分方程都能转化为变分问题——即需要先考虑变分的**逆问题**。

并不是所有微分方程都有使该微分方程自身为 Euler 方程的泛函。以二阶常微分方程为例，若它是泛函 $Q[y]=\int_a^b F(x,y,y')dx$ 的 Euler 方程，则必须满足 $\varphi_{y'} - \frac{d}{dx}\varphi_{y''}=0$。这个条件并非对所有微分方程都成立。

因此这类问题必须采用**加权残量法**。

### 4.3.2 加权残量法的基本思想

设在函数空间 H 中求解算子方程 $\mathbf{T}u = f$。

1. 在 H 中找一组线性独立的元素 $\varphi_1, \ldots, \varphi_n$，张成子空间 $M$
2. 找一组线性独立的权函数 $\omega_1, \ldots, \omega_n$，张成权函数空间 $W$
3. 在 M 中设 $u_0 = \sum a_i\varphi_i$
4. 代入方程得残量 $R = \mathbf{T}u_0 - f \neq 0$
5. 加权残量法的基本思想是使残量 R 向权函数空间 W 的投影为零：
$$(R, \omega_i) = 0,\quad i=1,2,\ldots,n$$

不同基函数和权函数的选择构成了不同的近似解法。

### 4.3.3 常见加权残量法

| 方法 | 权函数 | 特点 |
|------|--------|------|
| **Galerkin 法** | $\varphi_i$（基函数本身） | 最常用，系数对称 |
| **最小二乘法** | $\partial R/\partial c_i$ | 数学意义明确 |
| 配点法 | $\delta(x-x_i)$ | 实现简单 |
| 子域法 | 子域上为 1 | 物理解释清晰 |
| 矩法 | $1, x, x^2, \ldots$ | 最早的方法之一 |

---

## 4.4 有限单元法（The Finite Element Method）

### 4.4.1 概述

以一维问题为例。对于方程 $- (P u')' + r u = f$，其等效积分弱形式为：
$$Q[u(x)] = \frac12\int_a^b (P u'^2 + r u^2 - 2f u)dx$$

只要求试函数 $u(x)$ 具有 $C^0$ 连续性。

### 4.4.2 FEM 的数学描述

在区间 $[a,b]$ 上任意选取离散点 $a=x_0 < x_2 < \ldots < x_i < x_{i+1} < \ldots < x_n = b$。$x_i$ 称为节点，区间 $[x_i, x_{i+1}]$ 称为单元。

对每个节点 $x_i$ 建立线性基函数：
$$\varphi_i = \begin{cases}
1 & x = x_i \\
\text{线性} & x\in(x_{i-1}, x_{i+1}) \\
0 & \text{其他}
\end{cases}$$

由于每个基函数的值在对应的节点处为 1、在其他节点处为 0（**Kronecker 性质**），近似解的系数 $a_1, \ldots, a_n$ 就是 $v(x)$ 在离散点 $x_1, \ldots, x_n$ 处的函数值 $u_i$。基函数的**影响域（支撑集）非常小**，因此所得线性方程组的系数矩阵非常稀疏。

**FEM 的核心思想**：用 $u_n$ 作为无限维函数空间中真实解 u 在子空间 M 中的最佳逼近。

---

## 4.5 结论

经典 FEM 首先通过变分原理找到微分方程所对应的变分问题（即对应的泛函）。对于力学问题，可以从变分方法直接建立泛函。

但一些复杂微分方程的泛函不易找到，此时可以采用**加权残量法**——直接计算基函数与方程两端的内积，从而得到离散的求解方程组。

---

## 检查你的理解

1. Ritz 法和 Galerkin 法对试探函数的要求有何不同？
2. Galerkin 法与加权残量法是什么关系？
3. 极小化序列（minimizing sequence）的含义是什么？
4. 为什么简支梁算例中 $w=a\sin(\pi x/l)$ 比 $w=ax(l-x)$ 精度高得多？
5. 加权残量法有哪几种常见方法？各有什么特点？

---

> **对应作业**：[HW3 Q2（试函数合法性）](../04-Homework-Solutions/2026w/HW3-Problem.md) · [HW3 Q3（弹性地基梁）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **往年参考**：[past/HW3/Homework3](../04-Homework-Solutions/past/HW3/Homework3.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW3/Ans%20to%20HM3_LIU%20Sai_handed%20in.md)
