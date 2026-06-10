# 第1章：FEA 概述与弹性力学基础

> **对应 PDF**：[`Chapter 1 Introduction to FEA.pdf`](../06-References/pdfs-originals/Chapter%201%20Introduction%20to%20FEA.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §5
> **相关作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：高等数学（多元微积分、偏微分）、线性代数（矩阵、行列式）、大学物理（静力学平衡）

---

## 1.1 为什么要用有限元法？

### 1.1.1 工程问题的数学本质

所有的工程问题本质上都可以归结为在特定条件下求解**偏微分方程（PDE）**。以固体力学为例：给定一个弹性体，知其几何形状、材料属性、所受外力和约束条件，求体内各点的位移、应变和应力分布。

这构成了一个**边值问题（BVP）**：
$$\begin{cases}
\text{平衡方程：} & \sigma_{ij,j} + f_i = 0 \quad \text{(在域 } \Omega \text{ 内)} \\
\text{几何方程：} & \varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i}) \quad \text{(在域 } \Omega \text{ 内)} \\
\text{物理方程：} & \sigma_{ij} = D_{ijkl}\varepsilon_{kl} \quad \text{(在域 } \Omega \text{ 内)} \\
\text{边界条件：} & u_i = \bar{u}_i \text{ 在 } S_u \text{ 上},\quad \sigma_{ij}n_j = \bar{T}_i \text{ 在 } S_\sigma \text{ 上}
\end{cases}$$

这 15 个方程 + 边界条件，理论上可以求解 15 个未知量（3 位移 + 6 应变 + 6 应力）。但问题在于——**解析求解 PDE 极其困难**。

### 1.1.2 何时无法解析求解？

| 困难来源 | 具体例子 | 为什么难 |
|----------|----------|----------|
| 复杂几何 | 含孔板、带缺口轴、圆角过渡 | 边界形状无法用初等函数描述 |
| 复杂载荷 | 接触、冲击、非均匀压力 | 载荷分布无法简单表达 |
| 复杂材料 | 各向异性、非均匀、非线性 | 本构张量 $D_{ijkl}$ 不是常数 |
| 多场耦合 | 热-力耦合、流-固耦合、压电 | 方程数量翻倍，耦合项复杂 |

材料力学通过引入大量简化假设（如平截面假设）来绕过这些困难，但这牺牲了精度和适用范围。弹性力学可以处理更一般的问题，但其 PDE 在工程实际中几乎无法解析求解。

### 1.1.3 FEM 的定位

**有限单元法（FEM/FEA）**是一种**数值求解场问题的近似方法**。它不追求精确的闭合解，而是在可接受的计算成本内给出足够精确的近似解。

与 FEM 并列的其他数值方法：

| 方法 | 核心理念 | 优势 | 局限 |
|------|----------|------|------|
| 有限差分法（FDM） | 用差商代替微分 | 实现简单、效率高 | 不规则边界处理困难 |
| **有限元法（FEM）** | 分片插值 + 变分原理 | **几何适应性强** | 理论较深 |
| 边界元法（BEM） | 只在边界上离散 | 降维、精度高 | 适用于线性齐次问题 |
| 无网格法 | 基于点的函数逼近 | 无需网格 | 计算量大、边界条件处理复杂 |

FEM 在几何适应性、理论完备性和计算效率之间取得了最佳平衡，因此成为工程数值分析的事实标准。

---

## 1.2 有限元法的核心思想

### 1.2.1 离散化

FEM 的第一步是将一个具有**无限自由度**的连续体（continuum）分割成有限个几何简单的子区域——**单元（element）**，单元之间通过**节点（node）**连接。这个过程称为**网格划分（meshing）**。

$$\text{连续体（无限自由度）} \xrightarrow{\text{网格划分}} \text{离散体系（有限自由度）}$$

例如：一根连续的悬臂梁可以划分为 100 个杆单元，每个单元 2 个节点，总共 101 个节点。原来连续梁上的每一**点**都有位移（无限自由度），现在只有 101 个节点处的位移是未知的（有限自由度）。

### 1.2.2 分片近似

在每个单元内部，我们用一个简单的**多项式（形函数）**来近似该单元内的位移场：
$$u^{(e)}(x) = \sum_{i=1}^{n} N_i(x)\,u_i^{(e)}$$

其中 $u_i^{(e)}$ 是节点位移（真正的未知数），$N_i(x)$ 是形函数。

> **这就是 FEM 的精髓：化整为零，积零为整。**
> 想象一条复杂曲线，用单一高次多项式需要 10 阶以上才能拟合；但如果分成 100 小段，每段用一条直线就能近似得很好。分段越多，精度越高。

### 1.2.3 为什么分片近似是突破？

在 FEM 之前，Ritz 法已经知道可以用 trial function 的线性组合来近似解。但 Ritz 法的 trial function 必须满足**所有边界条件**——对于有孔的板、有缺口的轴这类复杂形状，构造这样的全局函数几乎不可能。

FEM 的突破在于：将求解域**分割成小单元**，每个单元独立构造低阶形函数。这些形函数不需要在全局满足边界条件，只需要满足单元间的大致连续即可。这不仅使形函数构造变得极其简单，而且使整个流程可以**标准化、程序化**。

### 1.2.4 FEM 基本流程

```
前处理（Preprocessing）
├── 几何建模：创建或导入 CAD 模型
├── 网格划分：选择合适的单元类型，生成节点和单元
├── 材料属性：定义弹性模量、泊松比、密度等
├── 边界条件：施加约束（固定端、简支等）
└── 载荷：施加外力（集中力、分布力、重力等）

求解（Solution）
├── 对每个单元 e：
│   ① 确定形函数 N_i
│   ② 计算 [B] 矩阵（应变-位移）
│   ③ 确定 [D] 矩阵（应力-应变，由材料决定）
│   ④ 计算单元刚度矩阵 [k]_e = ∫[B]^T[D][B] dV
│   ⑤ 计算等效节点力 {f}_e
├── 组装总体刚度矩阵 [K] = Σ[k]_e
├── 引入边界条件
└── 求解 [K]{δ} = {F} → 得到节点位移

后处理（Postprocessing）
├── 计算单元应变：ε = [B]{δ}_e
├── 计算单元应力：σ = [D]ε
└── 结果可视化：云图、变形图、动画
```

### 1.2.5 解的收敛性

随着网格细化，FEM 近似解收敛到真实解。两种策略：

- **h 方法**：加密网格（减小 $h$，即单元尺寸），不改变单元阶数
- **p 方法**：提高单元多项式阶数（增大 $p$），不改变网格
- **hp 方法**：两者结合

---

## 1.3 弹性力学基础

### 1.3.1 什么是弹性力学？

**弹性力学（Theory of Elasticity）**是固体力学的重要分支，研究弹性体在外力作用下内部**应力**和**变形**分布的规律。

- **弹性（Elasticity）**：外力去除后，材料完全恢复到原始状态
- **弹性体（Elastic body）**：以弹性为唯一材料特性的物理对象

### 1.3.2 弹性力学 vs 材料力学

| 对比维度 | 材料力学 | 弹性力学 |
|----------|----------|----------|
| 研究对象 | 杆、梁、柱等结构构件 | 板、壳、块体等一般弹性体 |
| 基本假设 | 多（平截面假设、圣维南原理等） | 少（仅 6 个基本假设） |
| 求解方法 | 初等公式（$\sigma=My/I$ 等） | PDE 求解或数值方法 |
| 精度 | 工程近似 | 更精确 |
| 适用范围 | 规则构件、简单载荷 | 任意几何、任意载荷 |

> **举例**：材料力学中矩形截面杆扭转给出 $\tau_{\max}=T/(\alpha b h^2)$，但无法处理非圆截面、变截面或带键槽的扭转问题。弹性力学则可以处理所有这些情况——代价是需要求解更复杂的偏微分方程。

### 1.3.3 弹性力学的基本假设

1. **连续性（Continuity）**：物体内部紧密无隙，应力/应变是空间中的连续函数
2. **均匀性（Homogeneity）**：材料属性在各点相同（$E, \nu$ 为常数）
3. **各向同性（Isotropy）**：材料属性在各方向相同
4. **完全弹性（Perfect elasticity）**：应力-应变满足 Hooke 定律，加载/卸载路径一致
5. **小变形（Small deformation）**：变形远小于物体尺寸，可在变形前几何上列平衡方程
6. **无初始应力（No initial stress）**：无外力时内部应力为零

### 1.3.4 三大基本变量

在直角坐标系 $(x,y,z)$ 下：

**① 位移（3 个分量）**
$$\mathbf{u} = \begin{pmatrix} u \\ v \\ w \end{pmatrix}$$
描述弹性体内任一点在 $x, y, z$ 方向的位移分量。

**② 应变（6 个分量）**
$$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \varepsilon_z \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{pmatrix}$$

- $\varepsilon_x, \varepsilon_y, \varepsilon_z$：**正应变**（线应变），表示单位长度的伸长量
- $\gamma_{xy}, \gamma_{yz}, \gamma_{zx}$：**剪应变**（角应变），表示直角的角度变化量

> **直观理解**：取一个无限小的立方体，正应变描述它各边的伸缩比，剪应变描述它各面角度的变化。

**③ 应力（6 个分量）**
$$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x \\ \sigma_y \\ \sigma_z \\ \tau_{xy} \\ \tau_{yz} \\ \tau_{zx} \end{pmatrix}$$

- $\sigma_x, \sigma_y, \sigma_z$：**正应力**，垂直于截面，拉为正、压为负
- $\tau_{xy}, \tau_{yz}, \tau_{zx}$：**剪应力**，平行于截面

> 为什么是 6 个独立分量不是 9 个？因为 $\sigma_{ij} = \sigma_{ji}$（剪应力互等定理，由微元体力矩平衡可证），所以 $3\times3$ 应力张量只有 6 个独立分量。

---

## 1.4 三类基本方程

### 1.4.1 几何方程（应变-位移关系）

描述**位移**如何产生**应变**。在小变形假设下：

$$\varepsilon_x = \frac{\partial u}{\partial x},\quad \varepsilon_y = \frac{\partial v}{\partial y},\quad \varepsilon_z = \frac{\partial w}{\partial z}$$

$$\gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x},\quad \gamma_{yz} = \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y},\quad \gamma_{zx} = \frac{\partial w}{\partial x} + \frac{\partial u}{\partial z}$$

矩阵形式：
$$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$$

其中微分算子 $[\partial]$ 是 $6\times3$ 矩阵：
$$[\partial] = \begin{pmatrix}
\frac{\partial}{\partial x} & 0 & 0 \\
0 & \frac{\partial}{\partial y} & 0 \\
0 & 0 & \frac{\partial}{\partial z} \\
\frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
\frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x}
\end{pmatrix}$$

### 1.4.2 物理方程（应力-应变关系）

描述**应力**与**应变**之间的关系，即材料的本构行为。对各向同性线弹性材料（广义 Hooke 定律）：
$$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$$

其中弹性矩阵 $\mathbf{D}$ 由两个 Lame 常数 $\lambda$ 和 $G$ 表达：
$$\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}$$

Lame 常数与工程常数的关系：
$$\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)},\quad G = \frac{E}{2(1+\nu)}$$

### 1.4.3 平衡方程（应力-外力关系）

考虑弹性体内一个微元体 $dx\,dy\,dz$。以 $x$ 方向为例列平衡 $\sum F_x = 0$：

$$\left(\sigma_x + \frac{\partial\sigma_x}{\partial x}dx\right)dydz - \sigma_x dydz$$
$$+ \left(\tau_{yx} + \frac{\partial\tau_{yx}}{\partial y}dy\right)dxdz - \tau_{yx}dxdz$$
$$+ \left(\tau_{zx} + \frac{\partial\tau_{zx}}{\partial z}dz\right)dxdy - \tau_{zx}dxdy$$
$$+ f_x\,dxdydz = 0$$

化简并利用 $\tau_{yx}=\tau_{xy}, \tau_{zx}=\tau_{xz}$：
$$\frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{xz}}{\partial z} + f_x = 0$$

类似地得到 $y, z$ 方向，矩阵形式：
$$[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$$

### 1.4.4 边界条件

- **位移边界**（$S_u$ 上）：$\mathbf{u} = \bar{\mathbf{u}}$（如固定端位移为零）
- **力边界**（$S_\sigma$ 上）：$[\mathbf{n}]\boldsymbol{\sigma} = \bar{\mathbf{T}}$（如自由表面力为零）

### 1.4.5 三类方程的关系

```
u_i → ε_ij = (u_i,j+u_j,i)/2 → ε_ij → σ_ij = D_ijkl ε_kl → σ_ij → σ_ij,j + f_i = 0
位移    几何方程               应变    物理方程（本构）       应力    平衡方程
```

15 个方程 + 边界条件 ⇒ 问题适定。

---

## 1.5 平面问题

### 1.5.1 平面应力（Plane Stress）

适用场景：**薄板**（$z$ 向厚度远小于 $x,y$ 向尺寸），载荷在 $xy$ 平面内。

由于板面自由，$\sigma_z = \tau_{zx} = \tau_{zy} = 0$。

$$[D] = \frac{E}{1-\nu^2}\begin{pmatrix}1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2}\end{pmatrix}$$

**例**：含孔薄板拉伸、齿轮齿根应力。

### 1.5.2 平面应变（Plane Strain）

适用场景：**长结构**（$z$ 向远大于 $x,y$ 向尺寸），截面沿 $z$ 不变。

由于 $z$ 向变形被约束，$\varepsilon_z = \gamma_{zx} = \gamma_{zy} = 0$。

$$[D] = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}1 & \frac{\nu}{1-\nu} & 0 \\ \frac{\nu}{1-\nu} & 1 & 0 \\ 0 & 0 & \frac{1-2\nu}{2(1-\nu)}\end{pmatrix}$$

**例**：水坝、隧道、长墙。

> **口诀**："薄板无应力，厚板无应变"

---

## 1.6 最小势能原理

弹性系统总势能 $\Pi$：
$$\Pi = \underbrace{\frac12\int_\Omega \boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV}_{\text{应变能 }U} - \underbrace{\int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\bar{\mathbf{T}}\,dS}_{\text{外力势能 }V}$$

**最小势能原理**：在所有满足位移边界条件的可能位移场中，真实位移场使总势能取最小值。即 $\delta\Pi = 0$。

这个原理等价于平衡方程和力边界条件，是 FEM 的变分基础。

---

## 1.7 FEA 术语与单元类型

| 术语 | 英文 | 含义 |
|------|------|------|
| 节点 | Node | 离散模型中连接单元的点 |
| 单元 | Element | 基本构建块，有确定几何形状 |
| 自由度 | DOF | 节点上独立未知量的个数 |
| 网格 | Mesh | 节点和单元的集合 |
| 形函数 | Shape Function | 单元内插值函数 |
| 刚度矩阵 | Stiffness Matrix | 联系节点力和节点位移的矩阵 $[k]_e$ |

**单元类型**：
- 1D：杆单元（truss，仅轴向）、梁单元（beam，弯曲+轴向）
- 2D：三角形 CST/LST、四边形 Q4/Q8
- 3D：四面体 Tet、六面体 Hex

---

## 1.8 弹性力学发展简史

| 年份 | 人物 | 贡献 |
|------|------|------|
| 1678 | Hooke | 发现位移与外力成正比 |
| 1821 | Navier | 推导弹性体平衡方程 |
| 1823 | Cauchy | 线弹性边值问题完整表述 |
| 1855 | Saint-Venant | 扭转和弯曲理论 |
| 1908 | Ritz | 能量法近似求解 |
| 1943 | Courant | 分片三角形近似 |
| 1956 | Turner-Clough | 矩阵位移法解平面应力 |
| **1960** | **Clough** | **提出"Finite Element Method"** |

---

> **对应作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
