# 考前冲刺速查 — 考试范围章节分布 + 核心公式

> 按考试大纲四大模块排列，每个考点附核心公式。配合 [概念速查](../02-Concepts-Formulas/2-1-Concepts-Glossary.md) 和 [01 讲义笔记](../01-Lecture-Notes/) 使用。

---

## 一、Tensor（张量）

### 1.1 定义 Definition (scalar, vector, tensor)

- **标量** (scalar)：只有大小，无方向。不随坐标变换改变。零阶张量。
- **向量** (vector)：有大小和方向。变换律 $a_{i'}' = l_{i'j}a_j$。一阶张量。
- **张量** (tensor)：高阶推广。二阶张量变换律 $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$。

### 1.2 证明某量是张量 Proving tensors

> 详见下方证明步骤

1. 写出旧坐标系中的定义
2. 坐标变换 $x_{i'} = l_{i'j}x_j$，导数变换 $\dfrac{\partial}{\partial x_{j'}} = l_{mj'}\dfrac{\partial}{\partial x_m}$
3. 代入整理为 $l_{i'i}l_{j'j}\cdots \times \text{原分量}$ 形式
4. 与张量变换律对比 → 得证

**关键判据**：若某量满足 $\boxed{a_{i'j'\cdots}' = l_{i'i}l_{j'j}\cdots a_{ij\cdots}}$，则它是张量。

### 1.3 特殊符号 Special symbols

| 符号 | 定义 | 关键性质 |
|------|------|----------|
| Kronecker $\delta_{ij}$ | $\boxed{\delta_{ij}=\begin{cases}1&i=j\\0&i\neq j\end{cases}}$ | $\boxed{\delta_{ii}=3}$（3D），$\boxed{\delta_{ij}a_j=a_i}$ |
| Levi-Civita $e_{ijk}$ | $\boxed{\text{置换符号：偶排列+1，奇排列-1，有重复0}}$ | $\boxed{e_{123}=1}$，$\boxed{e_{ijk}e_{ist}=\delta_{js}\delta_{kt}-\delta_{jt}\delta_{ks}}$，$\boxed{e_{ijk}e_{ijt}=2\delta_{kt}}$ |

### 1.4 代数与张量的转换 Conversion between algebra and tensors

Einstein 求和约定：$a_i b_i = a_1b_1+a_2b_2+a_3b_3$（哑指标自动求和）

- **自由指标**：不重复出现 → 决定结果的阶数
- **哑指标**：重复出现 → 隐含求和

### 1.5 坐标旋转 Coordinate rotation

方向余弦矩阵 $[l]$：

$$l_{ij} = \cos(x_i', x_j),\quad \boxed{l_{ij}l_{ik}=\delta_{jk}},\quad [l]^T[l]=[I]$$

### 1.6 向量/矩阵在不同坐标系下的变换 Vector/Matrix transformation

| 量 | 变换公式 |
|----|----------|
| 向量 | $\boxed{a_{i'}' = l_{i'j}a_j}$ |
| 二阶张量 | $\boxed{a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}}$ |
| 三阶赝张量 | $\boxed{e_{i'j'k'}' = \det(\boldsymbol{L})\,e_{ijk}}$ |

### 1.7 点积与叉积 Scalar and cross products

$$\boxed{\boldsymbol{a}\cdot\boldsymbol{b}=a_ib_i},\quad \boxed{(\boldsymbol{a}\times\boldsymbol{b})_i=e_{ijk}a_jb_k}$$

**$\varepsilon$-$\delta$ 恒等式**（利用 Lagrange 公式）：

$$\boxed{e_{ijk}e_{ist}=\delta_{js}\delta_{kt}-\delta_{jt}\delta_{ks}}$$

---

## 二、Elasticity（弹性力学）

### 2.1 Cauchy 公式 Cauchy formula

斜截面上的应力分量：

$$\boxed{T_i^n = \sigma_{ji}l_j = \sigma_{ij}l_j}$$

矩阵形式：

$$\begin{pmatrix} T_1^n \\ T_2^n \\ T_3^n \end{pmatrix} = \begin{pmatrix} \sigma_{11}&\sigma_{21}&\sigma_{31} \\ \sigma_{12}&\sigma_{22}&\sigma_{32} \\ \sigma_{13}&\sigma_{23}&\sigma_{33} \end{pmatrix} \begin{pmatrix} l_1 \\ l_2 \\ l_3 \end{pmatrix}$$

由力矩平衡得**剪应力互等**：$\sigma_{ij}=\sigma_{ji}$（应力张量对称）

### 2.2 应力变换 Stress transformation

坐标旋转下，新旧坐标系应力关系：

$$\boxed{\sigma_{i'j'}' = l_{i'm}l_{j'n}\sigma_{mn}}$$

### 2.3 主应力与主应变 Principal stress/strain

**主应力特征方程**：

$$\boxed{|\sigma_{ij}-\sigma\delta_{ij}|=0}$$

展开：

$$\boxed{\sigma^3 - I_1\sigma^2 + I_2\sigma - I_3 = 0}$$

主应力 $\sigma_1,\sigma_2,\sigma_3$ 为特征方程的三个实根，对应主方向两两正交。

主应变 $\varepsilon_1,\varepsilon_2,\varepsilon_3$ 由同样方法求得。

### 2.4 应力/应变张量的不变量 Invariant of stress/strain tensors

**应力不变量**：

$$\boxed{I_1 = \sigma_{11}+\sigma_{22}+\sigma_{33} = \sigma_1+\sigma_2+\sigma_3}$$

$$\boxed{I_2 = \begin{vmatrix}\sigma_{11}&\sigma_{12}\\\sigma_{21}&\sigma_{22}\end{vmatrix}+\begin{vmatrix}\sigma_{22}&\sigma_{23}\\\sigma_{32}&\sigma_{33}\end{vmatrix}+\begin{vmatrix}\sigma_{33}&\sigma_{31}\\\sigma_{13}&\sigma_{11}\end{vmatrix} = \sigma_1\sigma_2+\sigma_2\sigma_3+\sigma_3\sigma_1}$$

$$\boxed{I_3 = \begin{vmatrix}\sigma_{11}&\sigma_{21}&\sigma_{31}\\\sigma_{12}&\sigma_{22}&\sigma_{32}\\\sigma_{13}&\sigma_{23}&\sigma_{33}\end{vmatrix} = \sigma_1\sigma_2\sigma_3}$$

应变不变量形式相同，将 $\sigma$ 换为 $\varepsilon$。

### 2.5 平衡方程 Equilibrium equation

$$\boxed{\sigma_{ij,j}+f_i=0 \quad (i=1,2,3)}$$

展开（无体力分量 $b_i$ 写法）：

$$\frac{\partial\sigma_{11}}{\partial x_1}+\frac{\partial\sigma_{21}}{\partial x_2}+\frac{\partial\sigma_{31}}{\partial x_3}+f_1=0$$

（其余 $i=2,3$ 类似，共 3 个方程）

### 2.6 几何方程 Geometric equation

**小变形应变张量**：

$$\boxed{\varepsilon_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i})}$$

展开：$\varepsilon_{11}=\dfrac{\partial u_1}{\partial x_1}$，$\varepsilon_{12}=\dfrac{1}{2}\left(\dfrac{\partial u_1}{\partial x_2}+\dfrac{\partial u_2}{\partial x_1}\right)=\dfrac{1}{2}\gamma_{12}$

### 2.7 Green 应变 Green strain (大变形)

$$\boxed{E_{ij}=\frac{1}{2}\left(\frac{\partial\xi_k}{\partial x_i}\frac{\partial\xi_k}{\partial x_j}-\delta_{ij}\right)}$$

利用位移 $u_i=\xi_i-x_i$，化简为：

$$\boxed{E_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i}+u_{k,i}u_{k,j})}$$

小变形时忽略二阶项 → 回退到 $\varepsilon_{ij}$

### 2.8 本构关系 Constitutive relation

**Green 公式**（本构关系，适用于任意弹性材料；注意与高斯散度定理区分）：

$$\boxed{\sigma_{ij}=\frac{\partial W}{\partial \varepsilon_{ij}}}$$

**广义 Hooke 定律**（各向同性线弹性）：

$$\boxed{\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2G\varepsilon_{ij}}$$

或写为：

$$\sigma_{ij}=2G\varepsilon_{ij}+\lambda\varepsilon_{kk}\delta_{ij}$$

**Lamé 常数**：

$$\boxed{\lambda=\frac{\nu E}{(1+\nu)(1-2\nu)}},\quad \boxed{G=\frac{E}{2(1+\nu)}=\mu}$$

### 2.9 各向同性线弹性材料的本构矩阵 Constitutive matrix of isotropic linear elastic materials

**三维一般形式**（用 Lamé 常数）：

$$\boxed{[D]=\begin{pmatrix}\lambda+2G&\lambda&\lambda&0&0&0\\\lambda&\lambda+2G&\lambda&0&0&0\\\lambda&\lambda&\lambda+2G&0&0&0\\0&0&0&G&0&0\\0&0&0&0&G&0\\0&0&0&0&0&G\end{pmatrix}}$$

其中 $\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}$，$G = \frac{E}{2(1+\nu)}$

**平面应力**（薄板 $\sigma_z=0$，从三维退化为 3×3）：

$$\boxed{[D]=\frac{E}{1-\nu^2}\begin{pmatrix}1&\nu&0\\\nu&1&0\\0&0&\frac{1-\nu}{2}\end{pmatrix}}$$

**平面应变**（长柱体 $\varepsilon_z=0$，从三维退化为 3×3）：

$$\boxed{[D]=\frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}1&\frac{\nu}{1-\nu}&0\\\frac{\nu}{1-\nu}&1&0\\0&0&\frac{1-2\nu}{2(1-\nu)}\end{pmatrix}}$$

> 💡 **关系**：三维 6×6 矩阵 → 施加约束（$\sigma_z=0$ 或 $\varepsilon_z=0$）→ 退化为 2D 3×3 矩阵。平面应力和平面应变不是并列的，都是三维的特例。

### 2.10 边界条件 Boundary conditions

- **位移边界条件**（本质/Dirichlet）：$\boxed{u_i=\bar{u}_i}$ 在 $\Gamma_u$ 上
- **力边界条件**（自然/Neumann）：$\boxed{\sigma_{ij}n_j=\bar{T}_i}$ 在 $\Gamma_\sigma$ 上

> 两类边界互补：$\Gamma_u \cup \Gamma_\sigma = \Gamma$（整个边界），$\Gamma_u \cap \Gamma_\sigma = \varnothing$

### 2.11 三组方程推导中的假设 Assumptions in the deduction of three sets of equations

1. **连续性**（Continuity）：物体内部由连续介质充满
2. **均匀性**（Homogeneity）：物体各部分材料性质相同
3. **各向同性**（Isotropy）：材料性质与方向无关
4. **完全弹性**（Perfect elasticity）：应力应变满足 Hooke 定律且可逆
5. **小变形**（Small deformation）：位移远小于物体尺寸
6. **无初始应力**（No initial stress）：外载荷施加前无预应力

---

## 三、Variational Principles（变分原理）

### 3.1 泛函 Functional

$$Q[y]:D\to\mathbb{R}$$

**线性泛函**：$Q[c_1y_1+c_2y_2]=c_1Q[y_1]+c_2Q[y_2]$

### 3.2 泛函的变分 Variation of functional

**定义一（常规法）**：泛函增量 $\Delta Q = Q[y+\delta y] - Q[y]$ 中，对 $\delta y$ 线性的部分即为一阶变分。

例：$Q[y]=\int y^2 dx$，$\Delta Q = \int[(y+\delta y)^2 - y^2]dx = \underbrace{\int 2y\delta y\,dx}_{\delta Q} + \underbrace{\int(\delta y)^2 dx}_{\text{高阶小量}}$

**定义二（Lagrange 法，推荐）**：

$$\boxed{\delta Q = \dfrac{\partial}{\partial\alpha}Q[y+\alpha\,\delta y]\bigg|_{\alpha=0}}$$

即：引入参数 $\alpha$，对 $\alpha$ 求导再令 $\alpha=0$。不需要手动判断线性项和高阶项。


**基本性质**：变分与微分可交换 $\boxed{\delta\left(\dfrac{dy}{dx}\right)=\dfrac{d}{dx}(\delta y)}$

### 3.3 泛函极值 Functional extremum

极值必要条件（驻值条件）：$\boxed{\delta Q=0}$

- $\delta Q=0$ 且 $\delta^2 Q>0$ → 极小值
- $\delta Q=0$ 且 $\delta^2 Q<0$ → 极大值

实际应用中通常只考虑一阶变分，充分性由物理背景保证。

### 3.4 Euler 方程 Euler equation

泛函极值问题 $\delta Q=0$ 转化为微分方程的工具。给定泛函 $Q[y]=\int F(x,y,y')dx$，Euler 方程就是其极值函数必须满足的微分方程。FEM 中反过来用：先建立泛函，变分得 Euler 方程（=控制方程），再用 Ritz/Galerkin 法近似求解。

$$\boxed{\frac{\partial F}{\partial y}-\frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right)=0}$$

**含高阶导数**：

$$\boxed{\sum_{k=0}^{n}(-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right)=0}$$

**多元函数**：$\boxed{F_{y_i}-\dfrac{d}{dx}F_{y_i'}=0\quad(i=1,\ldots,n)}$

**多变量**（$z=z(x,y)$）：$\boxed{F_z-\dfrac{\partial}{\partial x}F_p-\dfrac{\partial}{\partial y}F_q=0}$（$p=z_x,q=z_y$）

**$F$ 不显含 $x$**：首次积分 $\boxed{F-y'F_{y'}=C}$

### 3.5 本质边界条件与自然边界条件 Essential and natural boundary conditions

变分后边界项：

$$\boxed{\delta Q = \int_a^b\left(F_y - \frac{d}{dx}F_{y'}\right)\delta y\,dx + \left.F_{y'}\delta y\right|_a^b}$$

- **本质（Essential）**：预先给定 $\delta y|_{\Gamma}=0$，如固定端 $u=0$
- **自然（Natural）**：由变分自动导出 $\boxed{\left.\dfrac{\partial F}{\partial y'}\right|_{\text{边界}}=0}$，如自由端 $M=0$

### 3.6 泛函的条件极值 Conditional extremum of functional

**Lagrange 乘子法**：引入 $\lambda$，构造新泛函：

$$\boxed{\hat{Q}=Q[y]+\lambda\cdot C[y]}$$

其中 $C[y]=0$ 为约束条件。对 $\hat{Q}$ 取变分 $\delta\hat{Q}=0$，同时解出 $y$ 和 $\lambda$。

### 3.7 Euler 方程的推广形式 Extended forms of Euler equation

**含高阶导数**（$Q=\int_a^b F(x,y,y',y'')dx$）：

$$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0}$$

推广到 $n$ 阶：$\boxed{F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} - \cdots + (-1)^n\frac{d^n}{dx^n}F_{y^{(n)}} = 0}$

**多个独立函数**（$Q=\int_a^b F(x,y_1,\ldots,y_n,y_1',\ldots,y_n')dx$）：

$$\boxed{F_{y_i} - \frac{d}{dx}F_{y_i'} = 0,\quad i=1,2,\ldots,n}$$

**多元函数**（$Q=\iint_D F(x,y,z,p,q)\,dxdy$，$p=z_x$，$q=z_y$）：

$$\boxed{F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0}$$

### 3.8 变分法在力学中的应用 Applications of variation method in mechanics

核心联系：**变分问题 ↔ 微分方程边值问题**（当泛函存在时）

$$\boxed{\delta Q = 0 \;\Longleftrightarrow\; \text{Euler 方程} + \text{自然边界条件}}$$

力学中：$\delta\Pi=0$（最小势能原理）$\Longleftrightarrow$ 平衡方程 + 边界条件

### 3.9 虚功原理 Principle of virtual work

**可能位移**：满足几何方程和位移 BC 的位移场
**可能应力**：满足平衡方程和力 BC 的应力场
**虚位移** $\delta u_i$：满足位移 BC 的任意微小位移变分，$\delta\varepsilon_{ij}=\frac12(\delta u_{i,j}+\delta u_{j,i})$
**虚应力** $\delta\sigma_{ij}$：满足力 BC 的任意微小应力变分，$\delta\sigma_{ij,j}=0$

**虚功方程**（核心公式）：

$$\boxed{\int_V F_i\,\delta u_i\,dV + \int_{S_\sigma} \bar{p}_i\,\delta u_i\,dS = \int_V \sigma_{ij}\,\delta\varepsilon_{ij}\,dV}$$

左端 = 外力虚功，右端 = 内力虚功（应变能变分）。适用于任何小变形体，不依赖本构关系。

**虚位移原理**（$\delta u_i$ 任意 → 检验 $\sigma_{ij}$）：

$$\boxed{\int_V \sigma_{ij}\delta\varepsilon_{ij}\,dV = \int_V F_i\delta u_i\,dV + \int_{S_\sigma}\bar{p}_i\delta u_i\,dS}$$

$\Longleftrightarrow$ 平衡方程 $F_i+\sigma_{ij,j}=0$ + 力边界 $\bar{p}_i=\sigma_{ij}l_j$

> 💡 **虚功方程 vs 虚位移原理**：数学表达式相同，逻辑方向相反。
> - **虚功方程**：已知应力满足平衡 → 等式自动成立（从因到果）
> - **虚位移原理**：等式对所有 $\delta u$ 成立 → 应力必然平衡（从果到因）
> - FEM 中主要用虚位移原理：用它导出弱形式，再离散为 $[K]\{u\}=\{F\}$

**虚应力原理**（$\delta\sigma_{ij}$ 任意 → 检验 $u_i$）：

$$\boxed{\int_V \varepsilon_{ij}\delta\sigma_{ij}\,dV = \int_{S_u}u_i\delta p_i\,dS}$$

$\Longleftrightarrow$ 几何方程 $\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$ + 位移边界 $u_i=\bar{u}_i$

> 💡 **虚位移原理 vs 虚应力原理**：对偶。
> - **虚位移原理**：独立变分 $\delta u_i$ → 检验 $\sigma_{ij}$ → 输出平衡方程+力边界
> - **虚应力原理**：独立变分 $\delta\sigma_{ij}$ → 检验 $u_i$ → 输出几何方程+位移边界
> - FEM 中主要用虚位移原理（从位移出发求应力）；虚应力原理用于应力分析和混合元

**推导**（从虚功原理到两个原理）：

**虚功原理**：$\int_V F_i\delta u_i\,dV + \int_{S_\sigma} \bar{p}_i\delta u_i\,dS = \int_V \sigma_{ij}\delta\varepsilon_{ij}\,dV$

---

**虚位移原理**（取 $\delta u_i$ 为独立变分）：

右边利用 $\delta\varepsilon_{ij}=\frac12(\delta u_{i,j}+\delta u_{j,i})$ 和 $\sigma_{ij}$ 对称性：

$$\int_V \sigma_{ij}\delta\varepsilon_{ij}\,dV = \int_V \sigma_{ij}\delta u_{i,j}\,dV$$

分部积分（高斯公式）：

$$= \int_S \sigma_{ij}l_j\delta u_i\,dS - \int_V \sigma_{ij,j}\delta u_i\,dV$$

分解边界 $S = S_u \cup S_\sigma$，利用 $\delta u_i|_{S_u}=0$：

$$= \int_{S_\sigma} \sigma_{ij}l_j\delta u_i\,dS - \int_V \sigma_{ij,j}\delta u_i\,dV$$

代回虚功原理，移项：

$$\int_V (F_i+\sigma_{ij,j})\delta u_i\,dV + \int_{S_\sigma}(\bar{p}_i-\sigma_{ij}l_j)\delta u_i\,dS = 0$$

由 $\delta u_i$ 任意 → 平衡方程 + 力边界。

---

**虚应力原理**（取 $\delta\sigma_{ij}$ 为独立变分）：

取两个都满足平衡的应力场 $\sigma_{ij}^1$、$\sigma_{ij}^2$，令 $\delta\sigma_{ij}=\sigma_{ij}^1-\sigma_{ij}^2$。虚功方程对两个状态分别成立，相减得（$F_i$、$\bar{p}_i$ 给定，故左端只剩 $S_u$ 项）：

$$\int_{S_u} \delta p_i\,u_i\,dS = \int_V \delta\sigma_{ij}\,\varepsilon_{ij}\,dV \tag{*}$$

**左边**：

$\delta p_i = \delta\sigma_{ij}l_j$（Cauchy 公式），高斯公式化为体积分：

$$\int_{S_u} \delta\sigma_{ij}l_j\,u_i\,dS = \int_V (\delta\sigma_{ij}u_i)_{,j}\,dV = \int_V \delta\sigma_{ij,j}\,u_i\,dV + \int_V \delta\sigma_{ij}\,u_{i,j}\,dV$$

利用 $\delta\sigma_{ij,j}=0$（虚应力满足平衡），第一项消失。利用 $\delta\sigma_{ij}$ 对称性：

$$= \int_V \delta\sigma_{ij}\cdot\frac12(u_{i,j}+u_{j,i})\,dV$$

**代回 (*) 式**，移项得完整方程：

$$\int_V \delta\sigma_{ij}\left[\varepsilon_{ij} - \frac12(u_{i,j}+u_{j,i})\right]dV = 0$$

由变分法预备定理，$\delta\sigma_{ij}$ 在 $V$ 内任意：

$$\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$$

**补充位移边界条件**：在 $S_u$ 上，$\delta p_i$ 任意，由 (*) 式的边界项可得 $u_i = \bar{u}_i$。
。

### 3.10 功的互等定理 Reciprocal theorem (Betti's formula)

$$\boxed{\int_V F_i^{(1)} u_i^{(2)}\,dV + \int_S p_i^{(1)} u_i^{(2)}\,dS = \int_V F_i^{(2)} u_i^{(1)}\,dV + \int_S p_i^{(2)} u_i^{(1)}\,dS}$$

状态一的力在状态二的位移上做的功 = 状态二的力在状态一的位移上做的功。

### 3.11 最小势能原理 Principle of minimum potential energy

在一切**可能位移场**中，真实位移场使总势能 $\Pi$ 取最小值：

$$\boxed{\Pi=\underbrace{\frac{1}{2}\int_\Omega \sigma_{ij}\varepsilon_{ij}\,dV}_{\text{应变能 }U}+\underbrace{\left(-\int_\Omega f_iu_i\,dV-\int_{\Gamma_\sigma}T_iu_i\,dS\right)}_{\text{外力势 }W}}$$

极值条件 $\dfrac{\partial\Pi}{\partial u_i}=0$ → 平衡方程

### 3.12 弹性力学中的 Euler 方程 Euler equations in elastic mechanics

从总势能 $\delta\Pi=0$ 出发，对平面应力问题分部积分后得：

**域内（Euler 方程 = 平衡方程，以平面应力为例）**：

$$\boxed{\frac{E}{1-\mu^2}\left(\frac{\partial^2 u}{\partial x^2}+\frac{1-\mu}{2}\frac{\partial^2 u}{\partial y^2}+\frac{1+\mu}{2}\frac{\partial^2 v}{\partial x\partial y}\right)+F_x=0}$$

$$\boxed{\frac{E}{1-\mu^2}\left(\frac{\partial^2 v}{\partial y^2}+\frac{1-\mu}{2}\frac{\partial^2 v}{\partial x^2}+\frac{1+\mu}{2}\frac{\partial^2 u}{\partial x\partial y}\right)+F_y=0}$$

$\Longrightarrow$ 等价于平衡方程 $\sigma_{ij,j}+f_i=0$ + 力边界 $\sigma_{ij}n_j=\bar{T}_i$

### 3.13 变分问题的直接法与间接法 Direct and indirect methods of variation problems

- **间接法**：Euler 方程法（将变分问题转化为微分方程）
- **直接法**：Ritz 法、Galerkin 法、有限差分法（直接近似求解泛函极值）

### 3.14 微分方程的等效变分方程 Equivalent variation equations of differential equations

**等效积分形式**（Galerkin 弱形式）：

$$\boxed{\int_\Omega w_i\left(\mathbf{T}u - f\right)d\Omega = 0}$$

分部积分降低光滑度要求：$C^1 \to C^0$（导数阶数减一），使低阶单元可用。

### 3.15 有限差分法 Finite difference method

将区间 $[a,b]$ 等分为 $n+1$ 段，用差商代替导数：

$$\boxed{y' \approx \frac{y_{i+1}-y_i}{\Delta x}}$$

泛函转化为多元函数：

$$\boxed{Q[y] \approx \sum_{i=0}^n F\!\left(x_i,\,y_i,\,\frac{y_{i+1}-y_i}{\Delta x}\right)\Delta x = \Phi(y_1,\ldots,y_n)}$$

由 $\dfrac{\partial\Phi}{\partial y_i}=0\;(i=1,\ldots,n)$ 解出节点值。$n\to\infty$ 得**极小化序列**。

### 3.16 Ritz 法 Ritz method

1. 选基函数 $\varphi_i$（需满足**本质边界条件**）
2. 构造 $\boxed{u_n=\sum_{i=1}^n a_i\varphi_i}$
3. 代入泛函 $Q$，极值条件 $\boxed{\dfrac{\partial Q}{\partial a_i}=0\;(i=1,\ldots,n)}$
4. 解线性方程组 $\boxed{\boldsymbol{Ka}=\boldsymbol{b}}$，其中 $K_{ij}=\dfrac{\partial^2 Q}{\partial a_i\partial a_j}$

> **解题套路**：给定泛函 $Q[u]$，边界条件 $u|_{\Gamma}=g$
> 1. 若 $g=0$（齐次 BC）：直接选 $u_n = \sum a_i\varphi_i$，$\varphi_i$ 满足齐次 BC
> 2. 若 $g\neq 0$（非齐次 BC）：令 $u_n = u_0 + \sum a_i\varphi_i$，其中 $u_0$ 满足非齐次 BC（如 $u_0=x$ 满足 $u(1)=1$），$\varphi_i$ 满足齐次 BC（如 $\varphi_i$ 在边界为 0）
> 3. 代入泛函 $Q$，化为 $Q(a_1,\ldots,a_n)$
> 4. 列方程：$\partial Q/\partial a_i = 0$（$i=1,\ldots,n$）
> 5. 解出 $a_1,\ldots,a_n$

> **齐次 vs 非齐次边界条件**：
> - **齐次 BC**：$u|_{\Gamma}=0$（边界值为零），如简支梁 $w(0)=w(l)=0$
> - **非齐次 BC**：$u|_{\Gamma}=g\neq 0$（边界值不为零），如 $u(1)=1$
> - Ritz 法中，试函数只需满足齐次 BC 部分，非齐次部分用 $u_0$ 单独处理

### 3.17 Galerkin 法 Galerkin method

加权残量法中取**权函数 = 试函数本身**：

$$\boxed{\int_\Omega N_i\cdot R\,d\Omega=0\quad(i=1,\ldots,n)}$$

- 试函数需满足**全部边界条件**（位移 + 力）
- 当泛函存在时，与 Ritz 法等价

> **边界条件要求对比**：
>
> | | Ritz 法 | Galerkin 法 |
> |---|---|---|
> | **本质边界条件**（位移 BC） | 必须满足 | 必须满足 |
> | **自然边界条件**（力 BC） | 不要求 | **必须满足** |
> | **原因** | 变分时 delta u 在 S_u 上为 0，自动处理 | 残量 R 需在全部边界上正交 |
>
> **为什么 Galerkin 要求更严？** 因为 Galerkin 法直接处理微分方程的残量，如果力边界条件不满足，残量在边界上就不为零，加权积分就无法正确消残。
>
> **举例**：悬臂梁自由端 $w''(l)=0$（弯矩=0）和 $w'''(l)=0$（剪力=0）都是力边界条件。Galerkin 法的试函数必须同时满足这两个条件，而 Ritz 法只需满足固定端的位移条件。

> **解题套路**：给定微分方程 $Lu = f$，边界条件 $B(u)=0$
> 1. 选试函数 $u_n = \sum a_i\varphi_i$（满足全部 BC）
> 2. 算残量 $R = Lu_n - f$
> 3. 列方程：$\int_0^l R\cdot\varphi_i\,dx = 0$（$i=1,\ldots,n$）
> 4. 解出 $a_1,\ldots,a_n$

### 3.18 加权残量法 Weighted residual method

基本思想：使残差 $R=L(\tilde{u})-f$ 在权函数空间中投影为零：

$$\boxed{\int_\Omega w_i\cdot R\,d\Omega = 0\quad(i=1,\ldots,n)}$$

| 方法 | 权函数 $w_i$ | 精度（一阶近似，简支梁） |
|------|--------|--------------------------|
| Galerkin | $N_i$（基函数本身） | **0.38%** |
| 最小二乘 | $\partial R/\partial c_i$ | **0.38%** |
| 配点法 | $\delta(x-x_i)$ | 21.16% |
| 子域法 | $1$（子域内） | 23.85% |

> **解题套路**：三种方法的区别只在权函数 $w_i$ 的选择
> - Galerkin：$w_i = \varphi_i$（基函数本身），精度最高
> - 最小二乘：$w_i = \partial R/\partial c_i$，精度同 Galerkin
> - 配点法：$w_i = \delta(x-x_i)$，在配点处令 $R=0$
> - 子域法：$w_i = 1$（子域内），在子域上令 $\int R\,d\Omega=0$

---

## 四、Basics of FEA（有限元分析基础）

### 4.1 单元、节点、DOF 的概念 Concepts of elements, nodes, DOFs

- **单元 (element)**：连续体的离散子域
- **节点 (node)**：单元的顶点/边中点/形心
- **自由度 (DOF)**：每个节点的独立参数数（Lagrange 型：$n$ 维问题有 $n$ DOF）

> **举例**：
> - **一维杆单元**：2 个节点，每节点 1 个 DOF（轴向位移 $u$）→ 共 2 DOF
> - **二维平面 CST**：3 个节点，每节点 2 个 DOF（$u, v$）→ 共 6 DOF
> - **三维四面体**：4 个节点，每节点 3 个 DOF（$u, v, w$）→ 共 12 DOF
> - **梁单元（Hermite）**：2 个节点，每节点 2 个 DOF（挠度 $w$ + 转角 $\theta$）→ 共 4 DOF

### 4.2 形函数 Shape functions

**定义**：形函数 $N_i$ 是定义在单元上的插值函数，用于将节点值 $u_i$ 插值为单元内任意点的值：

$$u(x,y) = \sum_{i=1}^n N_i(x,y)\,u_i$$

**核心性质**：

| 性质 | 表达式 | 含义 |
|------|--------|------|
| **Kronecker 性** | $N_i(\text{节点}j) = \delta_{ij}$ | 在自身节点为 1，其他节点为 0 |
| **单位分解性** | $\sum_{i=1}^n N_i = 1$ | 所有形函数之和恒等于 1 |
| **完备性** | 含刚体位移和常应变模态 | 保证收敛 |

**几何含义**：形函数就是自然坐标的推广——一维的 $\lambda_i$、二维的面积坐标 $L_i$、三维的体积坐标，都是形函数本身。

**各单元形函数**：

- **一维线性**：$\boxed{N_1=\dfrac{x_2-x}{L}}$，$\boxed{N_2=\dfrac{x-x_1}{L}}$

- **三角形（CST）**：$\boxed{N_i=\dfrac{1}{2\Delta}(a_i+b_ix+c_iy)}$，其中 $b_i=y_j-y_m$，$c_i=x_m-x_j$，$a_i=x_jy_m-x_my_j$

- **矩形（双线性，自然坐标）**：$\boxed{N_i=\dfrac{1}{4}(1+\xi\xi_i)(1+\eta\eta_i)}$

### 4.3 1D 和 2D FEM 的完整流程 Formulation of 1D and 2D FEA

```
结构离散 → 形函数 → [B]矩阵 → [D]矩阵 → 单元刚度矩阵[k]e
→ 总刚组装[K] = Σ[k]e → 引入边界条件 → 求解[K]{δ}={F} → 应力回算
```

**单元刚度矩阵**：

$$\boxed{[k]^e = \int_{\Omega_e} [B]^T[D][B]\,d\Omega}$$

**总体刚度矩阵**：

$$\boxed{[K] = \sum_e [k]^e}$$

**求解**：

$$\boxed{[K]\{\delta\} = \{F\}}$$

> **各矩阵的来源和含义**：
>
> | 矩阵 | 来源 | 含义 | 举例 |
> |------|------|------|------|
> | $[N]$ | 形函数 | 节点值→单元内场 | $N_i = \frac{1}{2\Delta}(a_i+b_ix+c_iy)$ |
> | $[B]$ | $[N]$ 对坐标求导 | 位移→应变 | $B_{ij} = \frac{\partial N_i}{\partial x_j}$ |
> | $[D]$ | 材料本构 | 应变→应力 | 1D: $D=EA$；平面应力: $\frac{E}{1-\nu^2}\begin{pmatrix}1&\nu&0\\\nu&1&0\\0&0&\frac{1-\nu}{2}\end{pmatrix}$ |
> | $[k]^e$ | $[B]^T[D][B]$ 积分 | 节点力→节点位移 | $[k]^e = t\Delta_e[B]^T[D][B]$（CST） |
> | $[K]$ | 总装 | 整体刚度 | $[K]=\sum[k]^e$ |
> | $\{F\}$ | 外力 | 节点力向量 | 体力+面力 |
>
> **一维杆单元示例**（最简单）：
> - $[B] = \frac{1}{L}\begin{pmatrix}-1&1\end{pmatrix}$，$[D] = EA$，$t=1$
> - $[k]^e = \frac{EA}{L}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$
>
> **二维 CST 单元示例**：
> - $[B] = \frac{1}{2\Delta}\begin{pmatrix}b_1&0&b_2&0&b_3&0\\0&c_1&0&c_2&0&c_3\\c_1&b_1&c_2&b_2&c_3&b_3\end{pmatrix}$
> - $[D]$ 见 §2.9（平面应力/应变）
> - $[k]^e = t\Delta[B]^T[D][B]$（$B$ 是常数矩阵，无需积分）

### 4.4 单元刚度矩阵的特性 Characteristics of element stiffness matrix

- **对称性**：$[k]^e=[k]^e{}^T$
- **奇异性**（无约束时）→ 引入边界条件后消失
- **主元恒正**
- **带状稀疏性**（总刚）

### 4.5 单元刚度矩阵元素的物理意义 Physical representation of each element in the stiffness matrix

$$\boxed{K_{ij} = \text{第 }j\text{ 个 DOF 产生单位位移（其余为零）时，第 }i\text{ 个 DOF 上的节点力}}$$

即 $k_{ij}$ = 第 $j$ 个自由度方向施加单位位移时，在第 $i$ 个自由度方向产生的反力。

### 4.6 收敛准则 Convergence criteria of FEA

> **什么时候用？** 设计新单元或判断单元是否可用时。考试常以"为什么 CST 能收敛"或"判断某单元是否满足收敛条件"的形式出现。

| 准则 | 含义 | 实际体现 | 反例 |
|------|------|---------|------|
| **完备性** | 试函数必须包含刚体位移和常应变状态 | $\sum N_i=1$（单位分解）保证刚体位移；线性形函数保证常应变 | 若 $\sum N_i \neq 1$，单元受均匀应力时会"缩小"或"放大" |
| **协调性** | 相邻单元交界面上位移必须连续 | CST 用线性函数，边上线性插值 → 位移连续 | 用二次函数但只在节点处连续，边上会出现间隙 |
| **分片试验** | 给单元施加常应变状态，看能否精确再现 | 若能 → 单元合格；若不能 → 单元有问题 | 用于检验新开发的单元类型 |

> 💡 **直觉理解**：
> - **完备性**：单元必须能表示"什么都不动"（刚体位移）和"均匀拉伸"（常应变），否则连基本状态都表示不了
> - **协调性**：相邻单元不能"脱开"或"重叠"，否则能量会从缝隙泄漏
> - **分片试验**：用简单已知答案的算例验证单元是否正确

### 4.7 长度坐标 Length coordinate (一维)

$$\boxed{\lambda_1=\frac{x_{i+1}-x}{L_i},\quad\lambda_2=\frac{x-x_i}{L_i},\quad\lambda_1+\lambda_2=1}$$

$$\boxed{L_f(x)=f(x_i)\lambda_1+f(x_{i+1})\lambda_2}$$

**长度坐标积分公式**：

$$\boxed{\int_{Q_1}^{Q_2} \lambda_1^{\alpha_1} \lambda_2^{\alpha_2}\,dx = L_i \frac{\alpha_1!\,\alpha_2!}{(\alpha_1 + \alpha_2 + 1)!}}$$

### 4.8 Lagrange 插值 Linear and high-order Lagrange interpolation (一维)

**线性**：同 4.7

**二次**：$L_f(x)=\sum_{k=0}^{2}f(x_k)\ell_k(x)$，其中 $\ell_k(x)=\prod_{j\neq k}\dfrac{x-x_j}{x_k-x_j}$

**$n$ 次**：$\boxed{L_f(x)=\sum_{k=0}^{n}f(x_k)\ell_k(x)}$，其中 $\boxed{\ell_k(x)=\prod_{j\neq k}\dfrac{x-x_j}{x_k-x_j}}$

### 4.9 Hermite 三次插值 Hermite cubic interpolation (Euler-Bernoulli 梁单元)

$w(\xi)=N_1w_1+N_2\theta_1+N_3w_2+N_4\theta_2$，$\xi=\dfrac{x-x_1}{L}$

$$\boxed{N_1=1-3\xi^2+2\xi^3,\quad N_2=L(\xi-2\xi^2+\xi^3)}$$
$$\boxed{N_3=3\xi^2-2\xi^3,\quad N_4=L(-\xi^2+\xi^3)}$$

$$\boxed{[k]_e=\frac{EI}{L^3}\begin{pmatrix}12&6L&-12&6L\\6L&4L^2&-6L&2L^2\\-12&-6L&12&-6L\\6L&2L^2&-6L&4L^2\end{pmatrix}}$$

### 4.10 面积坐标 Area coordinate (三角形)

$$\boxed{L_1=\frac{\Delta_1}{\Delta},\quad L_2=\frac{\Delta_2}{\Delta},\quad L_3=\frac{\Delta_3}{\Delta}}$$

$$L_1+L_2+L_3=1$$

与直角坐标转换：$x=L_1x_1+L_2x_2+L_3x_3$，$y=L_1y_1+L_2y_2+L_3y_3$

线性形函数：$N_i=L_i$

**面积坐标积分公式**（极其实用）：

$$\boxed{\iint_{\Delta_e} L_1^{\alpha_1} L_2^{\alpha_2} L_3^{\alpha_3}\,dxdy = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!}{(\alpha_1 + \alpha_2 + \alpha_3 + 2)!} \cdot 2\Delta_e}$$

### 4.11 Lagrange 插值 Lagrange interpolation for triangular elements (线性与高阶)

**线性**：$N_i=L_i$

**二次**（6 节点）：

$$\boxed{\begin{aligned}
N_i &= L_i(2L_i-1) & \text{(顶点 } i=1,2,3\text{)} \\
N_4 &= 4L_1L_2,\quad N_5 = 4L_2L_3,\quad N_6 = 4L_3L_1 & \text{(边中点)}
\end{aligned}}$$

### 4.12 划线法 Method of "scraping line"

构造形函数的几何方法：令不含某个节点坐标的线族为常数，利用归一化条件确定形函数。适用于非规则单元（如 serendipity 元）。

### 4.13 矩形单元形函数 Shape functions for rectangular elements (自然坐标 $\xi,\eta\in[-1,1]$)

**4 节点双线性**（Lagrange 型）：

$$\boxed{N_i=\frac{1}{4}(1+\xi\xi_i)(1+\eta\eta_i)\quad(i=1,2,3,4)}$$

### 4.14 Serendipity 单元形函数 Shape functions for serendipity element

**8 节点 serendipity**（角点）：$\boxed{N_i=\frac{1}{4}(1+\xi\xi_i)(1+\eta\eta_i)(\xi\xi_i+\eta\eta_i-1)}$

**边中点**：$\boxed{N_i=\frac{1}{2}(1-\xi^2)(1+\eta\eta_i)}$ 或 $\boxed{N_i=\frac{1}{2}(1+\xi\xi_i)(1-\eta^2)}$

### 4.15 等参单元 Isoparametric elements

坐标变换与位移插值使用**相同形函数**和**相同节点**：

$$\boxed{x=\sum N_i(\xi,\eta)x_i},\quad \boxed{u=\sum N_i(\xi,\eta)u_i}$$

### 4.16 Jacobian 矩阵 Jacobian matrix

$$\boxed{[J]=\begin{pmatrix}\dfrac{\partial x}{\partial\xi}&\dfrac{\partial y}{\partial\xi}\\[6pt]\dfrac{\partial x}{\partial\eta}&\dfrac{\partial y}{\partial\eta}\end{pmatrix}=\begin{pmatrix}\sum\dfrac{\partial N_i}{\partial\xi}x_i&\sum\dfrac{\partial N_i}{\partial\xi}y_i\\\sum\dfrac{\partial N_i}{\partial\eta}x_i&\sum\dfrac{\partial N_i}{\partial\eta}y_i\end{pmatrix}}$$

单元刚度矩阵中的积分变换：$\boxed{\displaystyle\int\int f(x,y)\,dxdy=\int_{-1}^{1}\int_{-1}^{1}f(\xi,\eta)|\det[J]|\,d\xi d\eta}$

### 4.17 Gauss 数值积分 Numerical integration

**核心思想**：$n$ 个积分点（不等距），代数精度达 $(2n-1)$ 次。

**正交条件**确定积分点 $\xi_i$（$n$ 阶 Legendre 多项式 $P_n(\xi)$ 的根）：

$$\boxed{\int_{-1}^1 \xi^i P(\xi)\,d\xi = 0, \quad i=0,1,\ldots,n-1}$$

**积分公式**：

$$\boxed{\int_{-1}^1 f(\xi)\,d\xi \approx \sum_{i=1}^n A_i f(\xi_i)}$$

**常用一维 Gauss 积分点和权系数**：

| $n$ | 积分点 $\xi_i$ | 权系数 $A_i$ | 精度 |
|:---:|---------------|:-----------:|:---:|
| 1 | $0$ | $2$ | 1 |
| 2 | $\pm 1/\sqrt{3} \approx \pm 0.577350$ | $1, 1$ | 3 |
| 3 | $0,\;\pm\sqrt{3/5} \approx \pm 0.774597$ | $8/9,\;5/9,\;5/9$ | 5 |

**二维嵌套法**：

$$\boxed{\int_{-1}^1\int_{-1}^1 f(\xi,\eta)\,d\xi d\eta \approx \sum_{i=1}^n\sum_{j=1}^n A_i A_j f(\xi_i,\eta_j)}$$

> 权系数是乘积 $A_iA_j$（非 $A_i+A_j$）。4 节点等参元用 $2\times2$（4 点），8/9 节点用 $3\times3$（9 点）。

### 4.18 常用积分公式 Common integral formulas (Ritz/Galerkin 法高频)

**三角函数积分**（简支梁/悬臂梁算例必备）：

$$\boxed{\int_0^l \sin^2\frac{n\pi x}{l}\,dx = \frac{l}{2}} \quad (n=1,2,3,\ldots)$$

$$\boxed{\int_0^l \sin\frac{m\pi x}{l}\sin\frac{n\pi x}{l}\,dx = \begin{cases}l/2 & m=n \\ 0 & m\neq n\end{cases}}$$

$$\boxed{\int_0^l \cos^2\frac{n\pi x}{l}\,dx = \frac{l}{2}}$$

$$\boxed{\int_0^l \sin\frac{n\pi x}{l}\,dx = \frac{2l}{n\pi}(1-(-1)^n) = \begin{cases}4l/(n\pi) & n\text{ 奇数} \\ 0 & n\text{ 偶数}\end{cases}}$$

**推导依据**：半角公式 $\sin^2\theta = \frac{1-\cos 2\theta}{2}$，$\cos^2\theta = \frac{1+\cos 2\theta}{2}$

**多项式积分**（Ritz 法常用）：

$$\boxed{\int_0^l x^m(l-x)^n\,dx = \frac{m!\,n!}{(m+n+1)!}l^{m+n+1}}$$

**正交性**：三角函数族 $\{\sin(n\pi x/l)\}$ 在 $[0,l]$ 上正交——这是 Galerkin 法中系数矩阵为对角阵的根本原因。

> 💡 **记忆技巧**：$\int_0^l \sin^2(\cdot)dx = l/2$ 是最常用的，记住即可。其他公式可由正交性和半角公式推导。

> 祝考试顺利！
