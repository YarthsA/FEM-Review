# 第1章：FEA 概述与弹性力学基础

> **对应 PDF**：[`Chapter 1 Introduction to FEA.pdf`](../06-References/pdfs-originals/Chapter%201%20Introduction%20to%20FEA.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §5
> **相关作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：高等数学（多元微积分、偏微分、多重积分）、线性代数（矩阵运算、行列式）、大学物理（静力学平衡、Hooke 定律）

---

## 1.1 为什么需要有限元法？

### 1.1.1 工程问题的数学描述

几乎所有工程问题都可以用**偏微分方程（PDE）**来描述。以固体力学为例，我们需要求解一个弹性体在外力作用下的位移、应变和应力分布。这一问题的数学表述是一个**边值问题（Boundary Value Problem, BVP）**：

$$
\begin{cases}
\text{平衡方程：} & \sigma_{ij,j} + f_i = 0 \quad \text{在 }\Omega\text{ 内}\\
\text{几何方程：} & \varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})\\
\text{本构方程：} & \sigma_{ij} = D_{ijkl}\varepsilon_{kl}\\
\text{边界条件：} & u_i = \bar{u}_i \text{ 在 }S_u\text{ 上}\\
& \sigma_{ij}n_j = \bar{T}_i \text{ 在 }S_\sigma\text{ 上}
\end{cases}
$$

这 15 个方程（6 个应力 + 6 个应变 + 3 个平衡/几何方程）+ 边界条件，共同决定了场变量 $\{\sigma_{ij}, \varepsilon_{ij}, u_i\}$。

### 1.1.2 解析解的困难

对于简单的几何形状（如矩形梁、圆板），我们可以求得上述方程的精确解。例如，材料力学中悬臂梁自由端受集中力 $P$ 的挠度公式：

$$w_{\max} = \frac{PL^3}{3EI}$$

但一旦出现以下任一情况，解析解就不再可行：

| 困难来源 | 例子 | 为什么难 |
|----------|------|----------|
| 复杂几何 | 有孔的板、带缺口的轴 | 边界条件无法用简单函数描述 |
| 复杂材料 | 非均匀、各向异性、非线性 | 本构方程 $D_{ijkl}$ 不是常数 |
| 复杂载荷 | 接触、冲击、热应力 | 载荷分布无法用初等函数表示 |
| 多场耦合 | 热-力-电耦合 | 方程数量翻倍，耦合项复杂 |

### 1.1.3 数值方法概览

| 方法 | 核心思想 | 优点 | 缺点 |
|------|----------|------|------|
| **有限差分法 (FDM)** | 用差商代替导数 | 实现简单 | 不规则边界困难 |
| **有限元法 (FEM)** | 分片插值 + 变分原理 | 几何适应性强 | 理论较深 |
| **边界元法 (BEM)** | 只在边界上离散 | 降维处理 | 适用于线性问题 |
| **无网格法 (Meshfree)** | 基于点的插值 | 无网格生成成本 | 计算量大 |

FEM 是工程中最广泛使用的数值方法，因为它在几何适应性、理论完备性和计算效率之间取得了最佳平衡。

---

## 1.2 有限元法的核心思想

### 1.2.1 离散化（Discretization）

FEM 的第一要义是**离散化**：将一个具有无限自由度的连续体（continuum），分割成有限个几何简单的**单元（element）**，单元之间通过**节点（node）**连接。

$$
\text{连续体（无限自由度）} \xrightarrow{\text{网格划分}} \text{离散体（有限自由度）}
$$

### 1.2.2 分片近似（Piecewise Approximation）

在每个单元内部，我们用一个简单的**多项式函数**（称为**形函数**或**插值函数**）来近似该单元内的位移场：

$$u^{(e)}(x,y) = \sum_{i=1}^{n} N_i(x,y) u_i^{(e)}$$

其中：
- $u_i^{(e)}$ 是单元 $e$ 中第 $i$ 个节点的位移值（未知数）
- $N_i(x,y)$ 是第 $i$ 个节点的形函数，满足 $N_i(\text{节点}_j) = \delta_{ij}$

> **关键洞察**：如果我们在整个求解域上用一个全局多项式来近似，那就是 **Ritz 法**。但 Ritz 法的 trial function 必须满足**所有边界条件**，对于复杂几何几乎不可能构造。FEM 的突破在于：**将求解域分割成小单元，每个单元独立构造低阶多项式**，这样形函数极其简单，且不必满足单元间的全局协调条件。

### 1.2.3 为什么分片近似有效？

想象一条曲线：在整个区间上用单一高次多项式拟合可能需要 10 阶以上多项式，但在每个小区间上用线性函数拟合（也就是用许多折线段）却非常简单。当分段足够细时，折线无限逼近曲线。

这就是 FEM 的精髓：**化整为零、积零为整**。

### 1.2.4 FEM 的基本流程

```
┌─────────────────────────────────────────────────────┐
│ 1. 前处理 (Preprocessing)                          │
│    └── 几何建模 → 网格划分 → 材料属性 → 边界条件    │
├─────────────────────────────────────────────────────┤
│ 2. 求解 (Solution)                                 │
│    ├── 对每个单元 e:                               │
│    │   a. 确定形函数 N_i(x,y)                       │
│    │   b. 计算 [B] 矩阵（应变-位移矩阵）              │
│    │   c. 计算 [D] 矩阵（弹性矩阵，由材料确定）        │
│    │   d. 计算单元刚度矩阵 [k]_e = ∫[B]^T[D][B]dV   │
│    │   e. 计算等效节点力 {f}_e                       │
│    ├── 组装总体刚度矩阵 [K] = Σ[k]_e                │
│    ├── 引入边界条件                                 │
│    └── 求解 [K]{δ} = {F} → 得到节点位移              │
├─────────────────────────────────────────────────────┤
│ 3. 后处理 (Postprocessing)                         │
│    └── 计算单元应变 ε = [B]{δ}_e                    │
│        → 计算单元应力 σ = [D]ε                      │
│        → 可视化（云图、变形图）                       │
└─────────────────────────────────────────────────────┘
```

---

## 1.3 弹性力学基础

### 1.3.1 弹性力学解决什么问题？

弹性力学（Theory of Elasticity）是固体力学的一个分支，研究弹性体在外力作用下的**应力（stress）**和**变形（deformation）**分布。它与材料力学的关系：

| 对比维度 | 材料力学 (Mechanics of Materials) | 弹性力学 (Theory of Elasticity) |
|----------|-----------------------------------|--------------------------------|
| 研究对象 | 杆、梁、柱等**结构构件** | 板、壳、实体等**一般弹性体** |
| 基本假设 | 较多简化假设（平截面假设等） | 假设少，更接近真实 |
| 求解精度 | 工程近似 | 更精确 |
| 求解方法 | 初等公式（$\sigma=My/I$ 等） | PDE 求解或数值方法 |
| 适用范围 | 规则构件 | 任意几何 + 任意载荷 |

> **举例**：材料力学中矩形截面杆的扭转问题简单地给出 $\tau_{\max} = T/(\alpha b h^2)$，但无法处理非圆截面、变截面或带有键槽的扭转问题。弹性力学则可以处理所有这些情况——代价是需要求解更复杂的方程。

### 1.3.2 弹性力学的基本假设

为了让问题在数学上可处理，弹性力学引入了以下假设：

1. **连续性**（Continuity）：物体内部没有空隙——应力、应变是空间中的连续函数
2. **均匀性**（Homogeneity）：材料属性在各点相同（$E$, $\nu$ 为常数）
3. **各向同性**（Isotropy）：材料属性在各方向相同
4. **完全弹性**（Perfect Elasticity）：应力-应变关系满足 Hooke 定律，加载和卸载路径一致
5. **小变形**（Small Deformation）：变形远小于物体尺寸，平衡方程可在变形前的几何上建立
6. **无初始应力**（No Initial Stress）：无外力时物体内部应力为零

> **偏离假设的后果**：当这些假设不成立时，问题进入"非线性"范畴——材料非线性（塑性、粘弹性）、几何非线性（大变形）、接触非线性。FEM 也能处理这些情况，但需要更复杂的本构模型和求解算法。

### 1.3.3 三大基本变量

在直角坐标系 $(x,y,z)$ 中，弹性力学涉及三类基本未知量：

**① 位移场（3 个分量）**
$$ \mathbf{u} = \begin{pmatrix} u \\ v \\ w \end{pmatrix} $$
其中 $u, v, w$ 分别表示沿 $x, y, z$ 方向的位移，它们都是 $(x,y,z)$ 的函数。

**② 应变场（6 个分量）**

$$ \boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x & \varepsilon_y & \varepsilon_z & \gamma_{xy} & \gamma_{yz} & \gamma_{zx} \end{pmatrix}^T $$

- $\varepsilon_x, \varepsilon_y, \varepsilon_z$ —— **正应变**（线应变），表示单位长度的伸长/缩短
- $\gamma_{xy}, \gamma_{yz}, \gamma_{zx}$ —— **剪应变**（角应变），表示直角的角度变化

**③ 应力场（6 个分量）**

$$ \boldsymbol{\sigma} = \begin{pmatrix} \sigma_x & \sigma_y & \sigma_z & \tau_{xy} & \tau_{yz} & \tau_{zx} \end{pmatrix}^T $$

- $\sigma_x, \sigma_y, \sigma_z$ —— **正应力**，垂直于截面
- $\tau_{xy}, \tau_{yz}, \tau_{zx}$ —— **剪应力**，平行于截面

> **为什么是 6 个而不是 9 个？** 应力张量和应变张量都是**对称的**：$\sigma_{ij} = \sigma_{ji}$，$\varepsilon_{ij} = \varepsilon_{ji}$。这是因为从微元体的力矩平衡可以证明剪应力互等 $\tau_{xy} = \tau_{yx}$。所以独立分量数从 9 减少到 6。

---

## 1.4 弹性力学的三类基本方程

### 1.4.1 几何方程（Kinematic/Strain-Displacement Relations）

几何方程描述了**位移与应变**之间的关系。在小变形假设下：

$$ \varepsilon_x = \frac{\partial u}{\partial x},\quad \varepsilon_y = \frac{\partial v}{\partial y},\quad \varepsilon_z = \frac{\partial w}{\partial z} $$

$$ \gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x},\quad \gamma_{yz} = \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y},\quad \gamma_{zx} = \frac{\partial w}{\partial x} + \frac{\partial u}{\partial z} $$

写成矩阵形式：
$$ \boldsymbol{\varepsilon} = [\partial] \mathbf{u} $$

其中微分算子矩阵为 $6\times3$：
$$
[\partial] = \begin{pmatrix}
\frac{\partial}{\partial x} & 0 & 0 \\
0 & \frac{\partial}{\partial y} & 0 \\
0 & 0 & \frac{\partial}{\partial z} \\
\frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
\frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x}
\end{pmatrix}
$$

**物理意义**：正应变 $\varepsilon_x = \partial u/\partial x$ 表示 $x$ 方向上一段微小线段 $dx$ 的伸长率；剪应变 $\gamma_{xy} = \partial u/\partial y + \partial v/\partial x$ 表示原来在 $xy$ 平面内互相垂直的两条微小线段夹角的变化量。

### 1.4.2 物理方程（Constitutive Relations）

物理方程描述了**应力与应变**之间的关系，也就是材料的本构行为。对线弹性各向同性材料，这就是广义 Hooke 定律：

$$ \boldsymbol{\sigma} = \mathbf{D} \boldsymbol{\varepsilon} $$

弹性矩阵 $\mathbf{D}$ 是一个 $6\times6$ 的对称矩阵。对各向同性材料，它由两个独立常数（Lame 常数 $\lambda$ 和 $G$）决定：

$$
\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}
$$

Lame 常数与工程常数 $E$（弹性模量）和 $\nu$（泊松比）之间的关系：
$$ \lambda = \frac{\nu E}{(1+\nu)(1-2\nu)},\quad G = \frac{E}{2(1+\nu)} $$

> **直观理解**：$G$ 就是剪切模量，衡量材料抵抗剪切变形的能力。$\lambda$ 则是体积模量的"同伴"——它与 $\kappa = E/[3(1-2\nu)]$（体积模量）有直接关系 $\kappa = \lambda + 2G/3$。

### 1.4.3 平衡方程（Equilibrium Equations）

平衡方程描述了**应力与外力**之间的关系。考虑弹性体内一个微元体 $dxdydz$，对其 $x$ 方向列平衡方程 $\sum F_x = 0$：

$$
\left(\sigma_x + \frac{\partial\sigma_x}{\partial x}dx\right)dydz - \sigma_x dydz
+ \left(\tau_{yx} + \frac{\partial\tau_{yx}}{\partial y}dy\right)dxdz - \tau_{yx}dxdz
+ \left(\tau_{zx} + \frac{\partial\tau_{zx}}{\partial z}dz\right)dxdy - \tau_{zx}dxdy
+ f_x dxdydz = 0
$$

化简、消去 $dxdydz$，并利用 $\tau_{yx}=\tau_{xy}, \tau_{zx}=\tau_{xz}$：

$$ \frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{xz}}{\partial z} + f_x = 0 $$

类似地得到 $y$ 和 $z$ 方向。写成矩阵形式：
$$ [\partial]^T \boldsymbol{\sigma} + \mathbf{f} = \mathbf{0} $$

其中 $\mathbf{f} = (f_x, f_y, f_z)^T$ 是单位体积上的体力（如重力）。

### 1.4.4 边界条件

上述方程只在弹性体内部成立。在边界上，我们需要指定：

**位移边界**（Dirichlet/Essential BC）：
$$ \mathbf{u}|_{S_u} = \bar{\mathbf{u}} $$
在 $S_u$ 部分边界上，位移是已知的（如固定端 $u=0$）。

**力边界**（Neumann/Natural BC）：
$$ [\mathbf{n}]\boldsymbol{\sigma}|_{S_\sigma} = \bar{\mathbf{T}} $$
在 $S_\sigma$ 部分边界上，外力是已知的（如自由表面 $\bar{\mathbf{T}}=0$）。

### 1.4.5 三类方程的关系图

```
位移 u_i        
  │ 几何方程 ε_ij = (u_i,j + u_j,i)/2
  ▼
应变 ε_ij      
  │ 物理方程 σ_ij = D_ijkl ε_kl
  ▼
应力 σ_ij      
  │ 平衡方程 σ_ij,j + f_i = 0
  ▼
外力 f_i
```

15 个方程（6 几何 + 6 物理 + 3 平衡）对应 15 个未知量（3 位移 + 6 应变 + 6 应力），加上合适的边界条件，弹性力学问题在数学上是**适定**的——解存在且唯一。

---

## 1.5 平面问题

三维问题往往过于复杂。当结构在一个方向上的尺寸远小于或远大于另外两个方向时，可以简化：

### 1.5.1 平面应力（Plane Stress）

**适用条件**：薄板（$z$ 方向厚度很小），载荷在 $xy$ 平面内，板面自由。

由于板面自由，$\sigma_z = \tau_{zx} = \tau_{zy} = 0$。简化的弹性矩阵：

$$ \mathbf{D} = \frac{E}{1-\nu^2} \begin{pmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{pmatrix} $$

**工程实例**：含孔薄板的拉伸、薄壁压力容器、齿轮齿根应力分析。

### 1.5.2 平面应变（Plane Strain）

**适用条件**：长结构（$z$ 方向很长），截面沿 $z$ 方向不变。

由于 $z$ 方向变形被约束，$\varepsilon_z = \gamma_{zx} = \gamma_{zy} = 0$（但 $\sigma_z \neq 0$）。简化的弹性矩阵：

$$ \mathbf{D} = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)} \begin{pmatrix}
1 & \frac{\nu}{1-\nu} & 0 \\
\frac{\nu}{1-\nu} & 1 & 0 \\
0 & 0 & \frac{1-2\nu}{2(1-\nu)}
\end{pmatrix} $$

**工程实例**：水坝的横向截面分析、隧道围岩应力、长墙的横向受力。

> 🔑 **记忆方法**："薄板无应力" → 平面应力（$\sigma_z=0$）；"厚板无应变" → 平面应变（$\varepsilon_z=0$）。

---

## 1.6 最小势能原理

### 1.6.1 总势能的组成

弹性系统的总势能 $\Pi$ 由两部分组成：

$$ \Pi = U + V $$

**① 应变能 (Strain Energy)** $U$：
$$ U = \frac12 \int_\Omega \boldsymbol{\varepsilon}^T \boldsymbol{\sigma} \, dV = \frac12 \int_\Omega \boldsymbol{\varepsilon}^T \mathbf{D} \boldsymbol{\varepsilon} \, dV $$

**② 外力势能 (Potential of External Forces)** $V$：
$$ V = -\int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\bar{\mathbf{T}}\,dS $$

负号表示外力做正功时势能减少。

### 1.6.2 最小势能原理的表述

> **在所有满足位移边界条件的可能位移场中，真实位移场使总势能取最小值。**
> $$ \delta\Pi = 0 $$

这等价于平衡方程和力边界条件。也就是说，求解 PDE 平衡方程的问题，可以转化为求解泛函 $\Pi$ 的极值问题——这正是有限元法的变分基础。

---

## 1.7 FEA 术语参考

| 术语 | 英文 | 定义 |
|------|------|------|
| 节点 | Node | 网格中连接单元的离散点，存储自由度值 |
| 单元 | Element | 网格的基本构建块，有确定几何形状（三角形、四边形等） |
| 自由度 | DOF | 节点上的独立未知量（位移元中为位移分量个数） |
| 网格 | Mesh | 节点和单元的集合，是求解域的离散表示 |
| 形函数 | Shape Function | 单元内插值函数中对应每个节点的基函数 |
| 刚度矩阵 | Stiffness Matrix | 联系节点力和节点位移的矩阵 $[k]_e$ |
| 总装 | Assembly | 将各单元刚度矩阵叠加到总体刚度矩阵的过程 |

### 单元类型一览

```
1D 单元（线单元）：杆单元(Truss) → 仅轴向力
                  梁单元(Beam)  → 弯曲+轴向

2D 单元（平面单元）：三角形(Tri) → CST(3节点)、LST(6节点)
                     四边形(Quad) → Q4(4节点)、Q8(8节点)

3D 单元（体单元）：四面体(Tet) → 4节点、10节点
                   六面体(Hex)  → 8节点、20节点
```

### 解的收敛性

FEM 近似解随网格细化收敛到真实解的两种策略：

- **$h$ 方法**：加密网格（减小单元尺寸），不改变单元阶数
- **$p$ 方法**：提高单元多项式阶数，不改变网格
- **$hp$ 方法**：两者结合

---

## 1.8 弹性力学发展简史

| 年份 | 人物 | 贡献 |
|------|------|------|
| 1678 | Hooke | 发现弹性体的位移与外力成正比 |
| 1821 | Navier | 推导弹性体平衡方程 |
| 1823 | Cauchy | 给出线弹性边值问题完整表述 |
| 1855 | Saint-Venant | 扭转和弯曲理论（圣维南原理） |
| 1908 | Ritz | 提出用 trial function 近似能量泛函 |
| 1943 | Courant | 分片三角形区域近似，预示 FEM 思想 |
| 1956 | Turner, Clough 等 | 矩阵位移法解平面应力问题 |
| **1960** | **Clough** | **首次提出"Finite Element Method"** |

---

## 检查你的理解

1. FEM 和 Ritz 法的本质区别是什么？为什么 FEM 能解决 Ritz 法无法处理的复杂几何问题？
2. 弹性力学的 6 个基本假设分别是什么？如果其中一个不成立，会发生什么？
3. 三类基本方程（几何、物理、平衡）分别描述了什么物理量之间的关系？
4. 平面应力和平面应变的区别是什么？如何判断一个实际问题属于哪种情况？
5. 总势能 $\Pi$ 由哪两部分组成？最小势能原理的物理意义是什么？

> **建议**：这些问题的答案应能在 2-3 句话内清晰表达。如果不能，说明对应知识点还未完全理解。

---

> **对应作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) — 需要理解 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 的张量性质
> [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md) — 需要熟练使用 $T_{ii}, T_{ij}T_{ij}$ 等运算
