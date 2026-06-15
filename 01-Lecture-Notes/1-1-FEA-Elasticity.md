# FEA 概述与弹性力学基础

> **对应课件**：[`Chapter 1 Introduction to FEA.pdf`](../06-References/pdfs-originals/Chapter%201%20Introduction%20to%20FEA.pdf) 绪论部分 · [原文MD](../../md_output/Chapter%201%20Introduction%20to%20FEA.md)
> **相关作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：高等数学、线性代数、大学物理

---

## 1.1 什么是有限元法？

> The finite element method (FEM), also called as finite element analysis (FEA), is a method for numerical solution of field problems.
> — 课程原始定义

> 💡 **理解关键**：FEM 的精髓不是"算"而是"拆"——把算不动的连续体拆成算得动的小块（单元），每块上假设简单的位移模式，再拼回去。这个"拆→算→拼"的策略贯穿整门课。

### 1.1.1 物理背景

FEM 解决的核心问题：**一个连续体结构在载荷和边界约束下，其内部的位移、应变和应力如何分布？**

对于简单几何形状和简单载荷，弹性力学可以给出解析解（exact solution）。但对于绝大多数实际工程问题——复杂的几何形状、多种材料属性、任意分布的载荷——解析解不存在。FEM 通过将连续体"切开"成有限个简单形状的单元（elements），用数值手段逼近真实解。

### 1.1.2 FEM 的基本原理

FEM 将连续体（continuum）离散为有限个单元（elements），通过在每个单元上假设位移函数（形函数），利用变分原理建立代数方程组来求解全场问题。将连续体离散化后，整个问题的分析就转化为对每个单元的分析和所有单元的重新组合。


### 1.1.3 FEM 的三种理解途径

1. **结构矩阵法（Structural Matrix Method）**：从结构力学的矩阵位移法发展而来
2. **变分法（Variational Method）**：基于最小势能原理或虚位移原理
3. **加权残量法（Method of Weighted Residuals）**：直接处理微分方程

FEM 的通用性极强——整个计算过程由计算机自动完成，只需根据不同的工程问题改变输入即可。这种方法彻底改变了分析解的限制。

### 1.1.4 其他数值方法对比

在 FEM 出现之前和同时期，存在几种替代数值方法：

| 方法 | 基本原理 | 优点 | 缺点 |
|------|---------|------|------|
| **有限差分法 (FDM)** | 用差分方程近似微分方程中的导数 | 概念简单；适合热传导/流体力学问题；平行边界的二维区域效果好 | 弯曲边界处理繁琐 |
| **配点法 (Collocation)** | 在选定的离散点上令残差为零 | 实现简单 | 精度依赖于配点选择 |
| **子域法 (Subdomain)** | 在每个子域上令残差积分平均为零 | — | — |
| **最小二乘法 (Least Squares)** | 最小化残差平方积分 | 与数据拟合兼容 | — |
| **Galerkin 法** | 加权函数取试探函数本身 | 对称矩阵；等价于变分法 | — |
| **变分法 (Variational)** | 基于能量泛函驻值原理 | 自然边界条件自动满足 | 要求存在泛函 |

> 前三种方法在 FEA 尚未普及的年代常用于解决连续体问题，但它们分别存在几何局限、或边界条件难以处理的问题。FEM 统一了这些方法的优点，同时克服了它们的缺点。

> 🔗 **连接**：这个表格里最需要记住的是 Galerkin 法和变分法的关系——当泛函存在时，Galerkin 法和变分法等价。第 4 章的 Ritz 法（变分路线）和 Galerkin 法（加权残量路线）实际上是从不同的门走进同一间屋子。

---

## 1.4 FEM 的八大优势

FEM 之所以成为主流工程分析方法的根本原因：

1. **通用性（Versatile）**：同一套计算框架适用于各种物理问题——固体力学、流体力学、热传导、电磁场等。

2. **无几何限制（No geometric restriction）**：不同于解析法需要简化几何，FEM 可以处理任意形状的物体。

3. **边界条件和加载无限制（Boundary conditions and loading are not restricted）**：可在任意位置施加任意类型的边界条件和载荷。

4. **材料属性无限制（Material properties are not restricted）**：各向同性、各向异性、复合材料、非均质材料均可处理。

5. **可处理不同行为的部件组合（Components with different behaviours）**：同一模型中可以包含不同本构关系的部件（如弹性和塑性区域）。

6. **模型可高保真逼近实际结构（FE models closely resemble the actual body）**：可以通过精确的几何建模和网格划分来忠实反映真实结构。

7. **精度可通过网格细化或高阶单元提高（Approximation is easily improved by refining the mesh or increasing the order of elements）**：FEM 的解可以通过**h-细化**（减小单元尺寸）或**p-细化**（提高形函数阶次）系统性地逼近精确解。

8. **可处理多物理场问题（Can deal with multi-physics problems）**：热-力耦合、流-固耦合、电-力耦合等多物理场问题可以在同一平台求解。

> 💡 **理解关键**：这八条中第 7 条"精度可系统性提高"是 FEM 的理论基石——它保证了 FEM 不仅仅"能用"而且"可靠"。离散→收敛→精确，这条逻辑链让 FEM 从"工程师的经验工具"变成了"可数学证明的数值方法"。

---

## 1.5 关键术语

FEM 中最核心的几个概念：

### Interpolation（插值）
通过节点值构造单元内部场的数学方法。例如，单元内位移 $u(x) = \sum N_i(x) u_i$，其中 $N_i$ 为形函数，$u_i$ 为节点位移。

### Elements（单元）
将连续体划分后得到的小块域。常见类型：1-D 线单元、2-D 平面单元、3-D 体单元。

### Nodes（节点）
单元的角点或边界上的特定位置。节点是未知变量（如位移）的承载点。

### DOF（自由度）
每个节点的独立未知变量个数。对三维弹性问题，每节点 3 个位移自由度（$u, v, w$）；对杆单元，每节点 1 个轴向位移自由度。

### Mesh（网格）
所有单元和节点的集合。网格质量直接影响计算精度和收敛性。

物理意义总结：**单元是"骨头"（定义拓扑结构），节点是"关节"（搭载未知量），D.O.F 是"运动维度"（每节点的未知量数量）。**

> ❌ **易错**：DOF 不是单元的属性，是节点的属性！很多初学者会把"单元的自由度"和"节点的自由度"混为一谈。正确说法：每节点有几个 DOF，系统总 DOF = 节点数 × 每节点 DOF。CST 单元有 3 个节点，每节点 2 个 DOF（u, v），所以一个 CST 单元的"单元自由度"是 6——但它本质上是 3 个节点的 DOF 拼起来的。

---

## 1.6 FEM 的发展简史

- **1908**：Ritz 提出用带未知量的试探函数近似能量泛函，得到求解未知量的方程组——但试探函数必须满足边界条件，复杂几何下极为困难
- **1943**：**Courant** 在解扭转问题时将截面划分为三角形区域，在每个三角形内假设翘曲函数均匀分布——**克服了 Ritz 法要求整体函数满足边界条件的困难**，预示了 FEM 思想
- **1954**：Turner 等人在波音公司开始用三角形单元分析飞机结构——内部报告于 1954 年完成
- **1955**：**Argyris** 发表能量理论和结构分析的多篇论文，统一了弹性结构的基本能量原理
- **1956**：**Turner, Clough, Martin, Topp** 在航空学会年会上正式发表论文"Stiffness and Deflection Analysis of Complex Structures"——介绍了用三角形和矩形单元求解平面应力的方法，即"直接刚度法"
- **1960**：**Clough** 在论文"The finite element method in plane stress analysis"中首次命名"Finite Element Method"——**学科诞生标志**
- **1960s**：Zienkiewicz, Argyris, Clough 三者并称 FEM 三大奠基人；冯康在中国独立提出了基于变分原理的有限元方法，但在国际上知名度较低
- **1970s**：在大中型计算机上得到广泛应用
- **1980s**：微型计算机普及，前/后处理软件成熟
- **1990s**：大规模结构体系分析成为现实

### 重要人物速查

| 人物 | 年份 | 贡献 |
|------|------|------|
| Courant | 1943 | 分片多项式近似——FEA 思想先驱 |
| Argyris | 1955 | 能量定理与结构分析统一理论 |
| Turner, Clough, Martin, Topp | 1956 | 直接刚度法——FEA 实用化开端 |
| Clough | 1960 | 命名"FEM"——学科诞生 |
| Zienkiewicz | 1960s-70s | FEM 理论体系化，写就经典教材 |
| 冯康 | 1960s | 中国独立发现变分形式的有限元方法 |


---

## 1.7 商业 FEM 软件

主流通用 FEM 软件：

| 软件 | 全称/开发方 | 主要应用领域 |
|------|------------|-------------|
| **ANSYS** | ANSYS Inc. | 通用多物理场（结构、热、流体、电磁） |
| **ABAQUS** | Dassault Systemes | 非线性结构分析（尤其擅长接触和大变形） |
| **NASTRAN** | NASA / MSC Software | 航空/航天结构分析 |
| **COMSOL** | COMSOL Inc. | 多物理场耦合（基于 PDE 直接输入） |
| **LS-DYNA** | LSTC / Ansys | 显式动力分析（碰撞、冲击、爆炸） |
| **HyperMesh** | Altair | 前处理/网格划分 |
| **SDRC/I-DEAS** | Siemens (原 SDRC) | 完整 CAD/CAM/CAE 一体化 |
| **ALGOR** | Autodesk | 通用结构/热分析 |
| **PATRAN** | MSC Software | 前后处理 |

此外还有 **MARC**、**SAP**、**NONSAP**、**ADINA** 等在特定历史时期有重要影响的软件。

---

## 1.8 弹性力学（Theory of Elasticity）

### 1.8.1 基本定义

- **弹性（Elasticity）**：外力去除后，物体恢复到原始状态的性质
- **弹性体（Elastic body）**：以弹性为唯一材料特性的物理对象
- **弹性力学**：研究弹性体在外载荷下应力和变形分布规律的学科

### 1.8.2 弹性力学 vs 材料力学

弹性力学是固体力学的重要分支，而材料力学可能是固体力学中历史最悠久的分支。材料力学关注结构构件的强度、刚度和稳定性问题，但材料力学只能处理杆、梁、柱等结构构件，板、壳和实体结构则很难处理。即使对杆梁构件，仍有一些问题留待解决。

弹性力学在原则上对弹性体的几何形状和载荷形式没有限制，允许比材料力学更真实的假设。例如梁弯曲理论中的"平截面假设"在某些条件下并不适用。

> 💡 **理解关键**：弹性力学 vs 材料力学的核心区别不在研究对象，而在假设的严格程度。材料力学需要额外假设（如平截面假设）来简化问题，弹性力学只依赖六条基本假设。FEM 建立在弹性力学基础上，所以天然比材料力学方法更通用。

### 1.8.3 弹性力学的基本假设

1. **连续性（Continuity）**：物体内部由连续介质充满，无空隙
2. **均匀性（Homogeneity）**：物体各部分的材料性质相同
3. **各向同性（Isotropy）**：材料性质与方向无关
4. **完全弹性（Perfect elasticity — Hooke's law）**：应力应变满足线性关系，且可逆
5. **小变形（Small deformation）**：位移远小于物体尺寸，忽略几何非线性
6. **无初始应力（No initial stress）**：外载荷施加前物体内部无预应力

> 这些假设确保控制方程为**线性偏微分方程**，从而叠加原理成立。当任一假设不满足时，问题变为非线性，FEM 求解需特殊处理。

> ⚠️ **难点**：这六条假设看似简单，但每一条被打破都对应一类非线性问题——连续性破坏→断裂力学、均匀性/各向同性破坏→复合材料、Hooke 定律破坏→塑性/超弹性、小变形假设破坏→几何非线性（大变形）。FEM 处理非线性问题时，本质上就是在逐个放宽这些假设。

---

## 1.9 弹性力学的三类基本变量与三类基本方程

> ⚠️ **难点**：这里是全章最核心的数学框架。"15 个未知量 = 15 个方程 = 封闭体系"是理解一切后续推导的前提。必须牢记三类变量的数量和物理含义——位移 3（u, v, w）、应变 6（为什么是 6？因为剪应变互等：γxy=γyx 等）、应力 6（同理，剪应力互等：τxy=τyx 等），3+6+6=15。

### 1.9.1 三类基本变量

三维弹性力学的基本变量：

**① 位移（3个分量）**：
$$\mathbf{u} = \begin{pmatrix} u \\ v \\ w \end{pmatrix}$$

**② 应变（6个分量）**：
$$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \varepsilon_z \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{pmatrix}$$

**③ 应力（6个分量）**：
$$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x \\ \sigma_y \\ \sigma_z \\ \tau_{xy} \\ \tau_{yz} \\ \tau_{zx} \end{pmatrix}$$

> **注意**：理论上应力/应变张量各 9 个分量，但剪应力互等（$\tau_{xy} = \tau_{yx}$，$\tau_{yz} = \tau_{zy}$，$\tau_{zx} = \tau_{xz}$）使得独立分量各为 6 个，共 15 个未知量。

### 1.9.2 三类基本方程

三类方程（几何 + 物理 + 平衡）形成**15 个方程 = 15 个未知量**的封闭体系。

#### 几何方程（位移-应变关系）
$$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$$

> ⚠️ **难点**：微分算子矩阵 $[\partial]$ 的写法是这个方程最容易出错的地方。注意这里 $[\partial]$ 是 6×3 矩阵（因为应变 6 行 × 位移 3 列），等式左边是 6×1 的应变向量，右边是 6×3 的矩阵乘 3×1 的位移向量。对应地，平衡方程中 $[\partial]^T$ 是 3×6 矩阵。两个方程是"共轭"关系——几何方程告诉你"位移的梯度=应变"，平衡方程告诉你"应力的散度=体力"。

其中微分算子矩阵：
$$[\partial] = \begin{pmatrix} 
\frac{\partial}{\partial x} & 0 & 0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial x} & \frac{\partial}{\partial z} & 0 \\
0 & 0 & \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
\end{pmatrix}^T$$

展开为：
$$\varepsilon_x = \frac{\partial u}{\partial x},\; \varepsilon_y = \frac{\partial v}{\partial y},\; \varepsilon_z = \frac{\partial w}{\partial z}$$

$$\gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x},\; \gamma_{yz} = \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y},\; \gamma_{zx} = \frac{\partial w}{\partial x} + \frac{\partial u}{\partial z}$$

物理意义：**应变是位移的空间梯度**。"正应变"描述线段长度的相对变化，"剪应变"描述夹角的变化。

> ❌ **易错**：工程剪应变 $\gamma_{xy} = \partial u/\partial y + \partial v/\partial x$，而不是 $\partial u/\partial y$ 单独一项——这是两个位移分量的交叉导数之和。在看应变矩阵时必须记住每个 $\gamma$ 是两个偏导数之和，单独一个偏导只描述了"剪切角的一部分"。

#### 物理方程（应力-应变关系 / 本构方程）

$$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$$

弹性矩阵 $\mathbf{D}$（Lame 常数形式）：
$$\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}$$

Lame 常数：$\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}$，$G = \frac{E}{2(1+\nu)}$

> 对简单拉伸：$\sigma_x = E\varepsilon_x$。但在三维中，一个方向的正应力会影响其他方向的正应变（泊松效应），因此需要完整的刚度矩阵。

> 💡 **理解关键**：$\mathbf{D}$ 矩阵的块状结构有明确的物理含义——左上 3×3 块是"正应力-正应变"耦合（含泊松效应），右下 3×3 块是"剪应力-剪应变"（对角，彼此不耦合）。非对角零块意味着"正应力不产生剪应变，剪应力不产生正应变"——这是各向同性材料特有的性质。

#### 平衡方程

$$[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$$

展开为：
$$\frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{zx}}{\partial z} + f_x = 0$$
$$\frac{\partial\tau_{xy}}{\partial x} + \frac{\partial\sigma_y}{\partial y} + \frac{\partial\tau_{yz}}{\partial z} + f_y = 0$$
$$\frac{\partial\tau_{zx}}{\partial x} + \frac{\partial\tau_{yz}}{\partial y} + \frac{\partial\sigma_z}{\partial z} + f_z = 0$$

物理意义：**微元体上的应力散度与体力之和为零**——即微元体在每一方向受力平衡（Newton 第二定律在静力学下的形式）。

> 🔗 **连接**：注意到几何方程 $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$ 和平衡方程 $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$ 使用了**同一个微分算子及其转置**。这不是巧合——它反映了"位移→应变"和"应力→平衡"之间的能量共轭关系，在第 4 章推导虚功原理时会反复利用这个对称性。

### 1.9.3 边界条件

- 位移边界（Dirichlet / essential）：$\mathbf{u}|_{S_u} = \bar{\mathbf{u}}$
- 外力边界（Neumann / natural）：$[\mathbf{n}]\boldsymbol{\sigma}|_{S_\sigma} = \mathbf{T}$

其中 $[\mathbf{n}]$ 是边界外法线方向的矩阵表示。

> ❌ **易错**：Dirichlet 边界又叫"本质边界"（essential），Neumann 边界又叫"自然边界"（natural）。这两个名字的来源——变分法中 Dirichlet 边界必须显式施加（否则不满足），Neumann 边界会被泛函的驻值条件"自动"满足。在 FEM 代码里，Dirichlet 边界通过"划行划列"施加，Neumann 边界通过等效节点力施加——千万不要搞反了！

### 1.9.4 最小势能原理

总势能泛函（包含应变能 U 和外力功 W）：

$$\Pi = \int_\Omega \frac12\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\mathbf{T}\,dS$$

- 第一项：应变能
- 第二项：体力势能
- 第三项：面力势能

> 在一切可能位移场中，真实位移场使总势能取最小值，即 $\delta\Pi = 0$。

这是 FEM 的核心理论基石——变分表述将求解偏微分方程问题转化为求解"最小化能量"问题。

> 💡 **理解关键**：为什么是"最小"而不是"最大"或"驻值"？物理直觉——弹性体在载荷下变形，它总是"选择"最省能量的变形方式。就像水往低处流，弹性体往能量最低的状态走。数学上 $\delta\Pi=0$ 只保证驻值，但第二变分 $\delta^2\Pi = \int \delta\boldsymbol{\varepsilon}^T\mathbf{D}\delta\boldsymbol{\varepsilon}\,dV > 0$（因为 $\mathbf{D}$ 正定）保证了它是极小值。

---

## 1.10 FEM 的求解步骤

以结构分析为例，完整流程为：

```
① 定义问题（Defining problems）
   明确物理模型：几何形状、材料属性、载荷、边界条件
               ↓
② 模型理想化（Model idealization）
   将实际结构简化为可分析的计算模型
               ↓
③ 前处理（Preprocessing）
   - 结构离散：按几何特点和精度要求划分单元和网格
   - 形成单元刚度矩阵：kₑ = ∫ BᵀEB dV
   - 形成等效节点荷载列阵
               ↓
④ 数值分析（Numerical analysis）
   - 集成总体刚度矩阵和荷载列阵
   - 引入强制边界条件
   - 求解方程 Ku = F，得到节点位移
               ↓
⑤ 后处理（Postprocessing）
   - 计算单元应变 ε = Bu
   - 计算单元应力 σ = Dε
   - 可视化变形、应力云图等
               ↓
⑥ 检查与诠释（Check and interpret results）
   变形的合理性、外力平衡检验、数量级检验
               ↓
⑦ 迭代（Reiteration）
   如需提高精度 → 加密网格 / 提高单元阶次 / 修正模型
```

关键判别：**网格加密后解的变化如果趋于稳定，说明解已经收敛**；如果变化剧烈，说明模型可能存在刚体位移、单元畸形或边界条件不正确。

---

## 1.11 单元类型

按空间维度分类：

| 维度 | 单元类型 | 典型形状 | 典型应用 |
|------|---------|---------|---------|
| 1-D (line) | 杆/梁单元 | 直线段 (2 节点) | 桁架、框架 |
| 2-D (plane) | 平面应力/应变 | 三角形 (3 节点)、四边形 (4 节点) | 薄板、深切剖面 |
| 3-D (solid) | 体单元 | 四面体 (4 节点)、六面体 (8 节点) | 通用三维结构 |

> 从 1-D 到 3-D，每个方向的单元有自己的形函数和刚度矩阵推导公式，但**推导逻辑完全一致**：位移插值 → 应变 → 应力 → 能量 → 变分 → 刚度矩阵。

> ❌ **易错**：平面应力与平面应变不是一回事！薄板（厚度方向自由伸缩）→平面应力（σz=0）；厚体/长体（厚度方向被约束）→平面应变（εz=0）。判据是厚度方向是否有约束，而不是几何外观。考试中看到"处于平面应力/应变状态"一定要先判断错了哪个假设。

---

## 1.12 变截面杆 FEM 完整算例


> 这是本章最核心的教学案例，展示了 FEM 从物理模型到数值求解的全过程。

### 1.12.1 问题描述

考虑一**变截面杆**，划分为 3 个单元（每个单元等截面），承受 3 个集中力 $F_1$、$F_2$、$F_3$ 和 $F_4$，如图所示。

**已知**：各段截面面积 $A_1$、$A_2$、$A_3$，长度 $l_1$、$l_2$、$l_3$，弹性模量 $E$（常数）
**求**：各加载点的位移 $u_2$、$u_3$、$u_4$（$u_1 = 0$ 为固定端）

**离散化**：共 3 个单元、4 个节点。

```
单元编号:  ①       ②       ③
节点编号:  1 ——— 2 ——— 3 ——— 4
截面面积:  A₁      A₂      A₃
杆段长度:  l₁      l₂      l₃
```

### 1.12.2 单元刚度矩阵

对轴向拉压杆单元，定义单元刚度系数 $k_i = \frac{EA_i}{l_i}$，单元平衡方程为：

**单元 1**（节点 1-2）：
$$\begin{bmatrix} k_1 & -k_1 \\ -k_1 & k_1 \end{bmatrix} \begin{Bmatrix} u_1 \\ u_2 \end{Bmatrix} = \begin{Bmatrix} P_{11} \\ P_{12} \end{Bmatrix}$$

**单元 2**（节点 2-3）：
$$\begin{bmatrix} k_2 & -k_2 \\ -k_2 & k_2 \end{bmatrix} \begin{Bmatrix} u_2 \\ u_3 \end{Bmatrix} = \begin{Bmatrix} P_{21} \\ P_{22} \end{Bmatrix}$$

**单元 3**（节点 3-4）：
$$\begin{bmatrix} k_3 & -k_3 \\ -k_3 & k_3 \end{bmatrix} \begin{Bmatrix} u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} P_{31} \\ P_{32} \end{Bmatrix}$$

其中 $P_{ij}$ 表示单元 $i$ 在节点 $j$ 处的内力（符号以节点对单元的作用方向为正向）。

### 1.12.3 节点平衡

对每个节点列平衡方程：

$$\begin{aligned}
\text{节点 1:} &\quad P_{11} = F_1 \\
\text{节点 2:} &\quad P_{12} + P_{21} = F_2 \\
\text{节点 3:} &\quad P_{22} + P_{31} = F_3 \\
\text{节点 4:} &\quad P_{32} = F_4
\end{aligned}$$

> 物理意义：每个节点上，相邻单元传给该节点的内力之和 = 作用在该节点上的外力。

### 1.12.4 总体刚度矩阵组装

将各单元的 2×2 刚度矩阵"扩展"到 4×4 全局矩阵中，在对应节点位置叠加：

**单元 1 贡献**（节点 1-2 → 矩阵行/列 1-2）：
$$\begin{bmatrix} 
k_1 & -k_1 & 0 & 0 \\ 
-k_1 & k_1 & 0 & 0 \\ 
0 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 0 
\end{bmatrix} \begin{Bmatrix} u_1 \\ u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} P_{11} \\ P_{12} \\ 0 \\ 0 \end{Bmatrix}$$

**单元 2 贡献**（节点 2-3 → 矩阵行/列 2-3）：
$$\begin{bmatrix} 
0 & 0 & 0 & 0 \\ 
0 & k_2 & -k_2 & 0 \\ 
0 & -k_2 & k_2 & 0 \\ 
0 & 0 & 0 & 0 
\end{bmatrix} \begin{Bmatrix} u_1 \\ u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} 0 \\ P_{21} \\ P_{22} \\ 0 \end{Bmatrix}$$

**单元 3 贡献**（节点 3-4 → 矩阵行/列 3-4）：
$$\begin{bmatrix} 
0 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 0 \\ 
0 & 0 & k_3 & -k_3 \\ 
0 & 0 & -k_3 & k_3 
\end{bmatrix} \begin{Bmatrix} u_1 \\ u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \\ P_{31} \\ P_{32} \end{Bmatrix}$$

将三个方程叠加，并代入节点平衡条件 $P_{11}=F_1$、$P_{12}+P_{21}=F_2$、$P_{22}+P_{31}=F_3$、$P_{32}=F_4$：

$$\boxed{\begin{bmatrix} 
k_1 & -k_1 & 0 & 0 \\ 
-k_1 & k_1+k_2 & -k_2 & 0 \\ 
0 & -k_2 & k_2+k_3 & -k_3 \\ 
0 & 0 & -k_3 & k_3 
\end{bmatrix} \begin{Bmatrix} u_1 \\ u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} F_1 \\ F_2 \\ F_3 \\ F_4 \end{Bmatrix}}$$

> 组装规则总结：**对角位置的刚度系数为连接到该节点的所有单元的 $k$ 之和**（例如 $K_{22} = k_1 + k_2$，因为节点 2 连接单元 1 和 2）；**非对角位置的非零项为 $-k$**（对应连接两节点的单元刚度）。

> 💡 **理解关键**：组装时各单元刚度矩阵"叠加"的物理含义是什么？因为节点上的合力 = 各单元传来的内力之和。每个扩展后的单元刚度矩阵实际上在说"这个单元对哪些节点有力贡献"——把三个单元的贡献加起来就是总合力。这就是"组装"的物理本质。

### 1.12.5 施加边界条件并求解

固定端 $u_1 = 0$，且 $F_1$ 变为未知反力 $R$。

$$\begin{bmatrix} 
k_1 & -k_1 & 0 & 0 \\ 
-k_1 & k_1+k_2 & -k_2 & 0 \\ 
0 & -k_2 & k_2+k_3 & -k_3 \\ 
0 & 0 & -k_3 & k_3 
\end{bmatrix} \begin{Bmatrix} 0 \\ u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} R \\ F_2 \\ F_3 \\ F_4 \end{Bmatrix}$$

删去第一行第一列（已知 $u_1=0$），得缩减方程组：

$$\begin{bmatrix} 
k_1+k_2 & -k_2 & 0 \\ 
-k_2 & k_2+k_3 & -k_3 \\ 
0 & -k_3 & k_3 
\end{bmatrix} \begin{Bmatrix} u_2 \\ u_3 \\ u_4 \end{Bmatrix} = \begin{Bmatrix} F_2 \\ F_3 \\ F_4 \end{Bmatrix}$$

求解得到单元应变：
$$\varepsilon_1 = \frac{u_2-u_1}{l_1} = \frac{u_2}{l_1},\quad \varepsilon_2 = \frac{u_3-u_2}{l_2},\quad \varepsilon_3 = \frac{u_4-u_3}{l_3}$$

单元应力：
$$\sigma_1 = E\varepsilon_1,\quad \sigma_2 = E\varepsilon_2,\quad \sigma_3 = E\varepsilon_3$$

反力 $R$ 从第一行回代：$R = -k_1 u_2$。

> ❌ **易错**：施加边界条件时常见的两个错误——①只划掉了一行而忘记划对应的列（必须行+列都划掉），②划掉后忘记把已知位移对应的项移到右边。这里我们采用了"划行划列法"（直接删除），等价但更简洁。另外注意 $R = -k_1 u_2$ 中的负号——杆被拉伸时 $u_2 > 0$，固定端反力方向朝左（与 x 正向相反），所以 $R < 0$ 是正确的。

---

## 1.13 总体刚度矩阵的性质

$$\begin{bmatrix} 
k_1 & -k_1 & 0 & 0 \\ 
-k_1 & k_1+k_2 & -k_2 & 0 \\ 
0 & -k_2 & k_2+k_3 & -k_3 \\ 
0 & 0 & -k_3 & k_3 
\end{bmatrix}$$

### 性质 1：非负对角元（Nonnegative $K_{ii}$）
$K_{ii} \ge 0$，对角元表示连接在节点 $i$ 上所有单元刚度系数之和。物理上，对角元代表要使节点 $i$ 产生单位位移（其他节点固定）所需的力——这是一个非负量。

### 性质 2：对称性（Symmetry）
$K_{ij} = K_{ji}$。来源于 **Betti 互等定理**——节点 $i$ 处单位位移在节点 $j$ 处产生的力，等于节点 $j$ 处单位位移在节点 $i$ 处产生的力。对称性意味着只需存储上三角部分，节省约一半内存。

### 性质 3：稀疏性（Sparsity）
大规模的总体刚度矩阵中大量元素为零。上例中 4×4 矩阵有 8 个零元素。对于大型三维结构（节点数 $10^4\text{--}10^7$），非零元素通常仅占总元素数的 **1% 或更少**。稀疏性来自"一个节点只通过单元与邻近节点连接"的局部特性。这正是 FEM 能在计算机上高效处理大型问题的关键——**稀疏矩阵存储和求解技术**。

### 性质 4：奇异性（Singularity — 无约束时）
在施加边界条件之前，$\det(\mathbf{K}) = 0$，即 $\mathbf{K}$ 不可逆。物理原因：未约束的系统可以发生**刚体位移**（平移和转动），存在无穷多组解。**消除奇异性的唯一方法是施加足够的本质边界条件**，消除所有刚体位移模式。

> 关键理解：$\mathbf{K}$ 的奇异性不是 bug 而是 feature——它反映了物理系统的真实特性。只有当约束足够时，系统才具有唯一解。


---

## 1.14 形函数法推导刚度矩阵（系统方法）

> 🔗 **连接**：这一节的推导链条是整个课程的"主干道"——形函数→B 矩阵→单元刚度→总体组装。第 5 章的 CST 平面单元、第 6 章的等参元都沿着同一条路走，只是形函数从线性的 $N_i=(1-\xi)/2$ 变成了三角形面积坐标或四边形双线性。这一节读透了，后面每章都是"换汤不换药"。

> 上一节的杆单元刚度矩阵 $k_i = EA_i/l_i$ 来自结构力学直觉。本节介绍**通用、系统化**的推导方法——形函数法，这是所有单元类型的统一推导框架。

### 1.14.1 坐标变换与形函数

将杆单元节点坐标从 $x_i$、$x_j$ 映射到自然坐标 $\xi \in [-1, 1]$：

$$\xi = \frac{2x}{l} \quad \text{（原点在杆中点）}$$

其中 $l = x_j - x_i$ 为杆长。

形函数 $N_i$、$N_j$ 定义为在自然坐标下满足插值条件的线性函数：

$$N_i(\xi) = \frac{1-\xi}{2},\quad N_j(\xi) = \frac{1+\xi}{2}$$

满足：$N_i(-1)=1,\; N_i(1)=0$；$N_j(-1)=0,\; N_j(1)=1$。

### 1.14.2 位移插值

单元内任意点的位移由节点位移插值得到：

$$\boxed{u(x) = N_i(\xi)u_i + N_j(\xi)u_j = \begin{bmatrix} N_i & N_j \end{bmatrix} \begin{Bmatrix} u_i \\ u_j \end{Bmatrix}}$$

几何坐标同样通过形函数插值：
$$x = N_i(\xi)x_i + N_j(\xi)x_j$$

> 这是**等参表示**——位移和几何使用相同的形函数（即等参元，详见第 6 章）。

### 1.14.3 应变-位移关系（应变矩阵 B）

$$\varepsilon = \frac{du}{dx} = \begin{bmatrix} \dfrac{dN_i}{dx} & \dfrac{dN_j}{dx} \end{bmatrix} \begin{Bmatrix} u_i \\ u_j \end{Bmatrix}$$

计算形函数对 $x$ 的导数（链式法则）：

$$\frac{dx}{d\xi} = \frac{dN_i}{d\xi}x_i + \frac{dN_j}{d\xi}x_j = -\frac{x_i}{2} + \frac{x_j}{2} = \frac{x_j - x_i}{2} = \frac{l}{2}$$

$$\frac{dN_i}{dx} = \frac{dN_i}{d\xi}\frac{d\xi}{dx} = \left(-\frac{1}{2}\right)\left(\frac{2}{l}\right) = -\frac{1}{l}$$

$$\frac{dN_j}{dx} = \frac{dN_j}{d\xi}\frac{d\xi}{dx} = \left(\frac{1}{2}\right)\left(\frac{2}{l}\right) = \frac{1}{l}$$

因此：

$$\varepsilon = \begin{bmatrix} -\dfrac{1}{l} & \dfrac{1}{l} \end{bmatrix} \begin{Bmatrix} u_i \\ u_j \end{Bmatrix} = \mathbf{Bu}$$

**应变矩阵** $\mathbf{B} = \begin{bmatrix} -1/l & 1/l \end{bmatrix}$。

> 注：对于线性形函数的杆单元，**应变和应力在单元内为常数**（因为 $\mathbf{B}$ 矩阵为常数矩阵）。

> ⚠️ **难点**：链式法则 $\frac{dN}{dx} = \frac{dN}{d\xi} \frac{d\xi}{dx}$ 在推导 B 矩阵时是核心操作。在 1-D 杆单元中 $\frac{d\xi}{dx} = \frac{2}{l}$ 很简单，但在第 6 章的 2-D/3-D 等参元中，链式法则会变成一个 Jacobian 矩阵求逆的问题——那是全课程计算最复杂的部分。这里的一维链式法则就是那个复杂版的简化预演。

### 1.14.4 应力

$$\sigma = E\varepsilon = E\mathbf{Bu} = \begin{bmatrix} -\dfrac{E}{l} & \dfrac{E}{l} \end{bmatrix} \begin{Bmatrix} u_i \\ u_j \end{Bmatrix}$$

### 1.14.5 单元刚度矩阵（能量法）

单元应变能：

$$\begin{aligned}
U &= \frac{1}{2}\int_V \sigma\varepsilon\,dV 
   = \frac{1}{2}\int_V (\mathbf{Bu})^T E (\mathbf{Bu})\,dV \\
  &= \frac{1}{2}\mathbf{u}^T \left(\int_V \mathbf{B}^T E\mathbf{B}\,dV\right) \mathbf{u}
\end{aligned}$$

其中体积分（对等截面杆元，$E$ 为常数，$\mathbf{B}$ 为常数矩阵）：

$$\begin{aligned}
\int_V \mathbf{B}^T E\mathbf{B}\,dV 
&= E \cdot (A l) \cdot \mathbf{B}^T \mathbf{B} \\
&= EA l \begin{bmatrix} -1/l \\ 1/l \end{bmatrix} \begin{bmatrix} -1/l & 1/l \end{bmatrix} \\
&= \frac{EA}{l} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}
\end{aligned}$$

令 $\mathbf{k} = \int_V \mathbf{B}^T E\mathbf{B}\,dV$，则：

$$\boxed{\mathbf{k} = \frac{EA}{l} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}}$$

**这就是杆单元的刚度矩阵的通用表达式。** 对比第 1.12.2 节中直接写出的 $k_i = EA_i/l_i$，此处 $\mathbf{k}$ 的推导完全由形函数出发，**适用于任意单元类型**。

> 💡 **理解关键**：$\mathbf{k} = \int_V \mathbf{B}^T \mathbf{D} \mathbf{B}\,dV$ 是整个 FEM 中最重要的一个公式，没有之一。它的逻辑链：$\mathbf{B}$ 把节点位移映射为应变，$\mathbf{D}$ 把应变映射为应力，$\mathbf{B}^T$ 把应力映射回节点力。所以 $\mathbf{B}^T\mathbf{D}\mathbf{B}$ 的物理含义是"节点位移→节点力"的刚度映射。为什么体积积分？因为单元内每一点都在产生应变能，需要在整个体积上累积。

### 1.14.6 最小势能原理导出平衡方程

外力功（含分布力 $q$ 和集中力 $P_i$、$P_j$）：

$$W = -\int_l q u\,dx - (P_i u_i + P_j u_j)$$

代入 $u = N_i u_i + N_j u_j$，积分得：

$$W = -\begin{bmatrix} u_i & u_j \end{bmatrix} \begin{Bmatrix} ql/2 + P_i \\ ql/2 + P_j \end{Bmatrix}$$

总势能 $\Pi = U + W$：

$$\Pi = \frac{1}{2}\mathbf{u}^T \mathbf{k} \mathbf{u} - \mathbf{u}^T \mathbf{f}$$

其中 $\mathbf{f} = \begin{Bmatrix} ql/2 + P_i \\ ql/2 + P_j \end{Bmatrix}$ 为等效节点力。

对 $\mathbf{u}^T$ 求导（$\partial\Pi/\partial\mathbf{u}^T = \mathbf{0}$）：

$$\boxed{\mathbf{ku} = \mathbf{f}}$$

> 这就是**单个单元的平衡方程**。所有单元的刚度矩阵组装后，即得整体结构方程组 $\mathbf{KU} = \mathbf{F}$。

> ❌ **易错**：分布力 $q$ 的等效节点力是 $ql/2$ 分配到每个节点——不是 $ql$（总载荷的一半分到每个节点）。物理原因是线性形函数 $N_i$ 从 $\xi=-1$ 处的 1 线性降到 $\xi=1$ 处的 0，$\int_{-1}^1 N_i q \cdot (l/2)d\xi = ql/2$。如果是二次形函数，分布就不均匀了。

### 1.14.7 推导链条总结

```
形函数 N(ξ) → 坐标变换 dx/dξ
    ↓
位移插值 u = Nu
    ↓
应变矩阵 B = dN/dx（链式法则）
    ↓
应力 σ = EBu
    ↓
单元刚度 k = ∫ BᵀEB dV
    ↓
应变能 U = (1/2) uᵀ k u
    ↓
外功 W = -uᵀ f
    ↓
总势能 Π = U + W → 变分 δΠ = 0
    ↓
单元平衡方程 ku = f
    ↓
总体组装 KU = F
```

**这条推导链条适用于所有单元类型（1-D / 2-D / 3-D），只需替换形函数和应变矩阵的具体形式。**

> 💡 **理解关键**：这条链条有 10 步，但核心只有 3 个"转换"——①位移→应变（B 矩阵，用链式法则）、②应变→应力（D 矩阵，用本构关系）、③应力→节点力（B^T 积分，用能量原理）。其他步骤是这 3 个转换的包装和重组。当你学新的单元类型时，只需要弄清楚这 3 个转换在那个单元中的具体形式。

---

## 1.15 网格质量与收敛性

### 1.15.1 网格质量

- **避免大长宽比（aspect ratio）**：理想长宽比接近 1:1
- **避免过于尖锐的内角**：三角形内角避免 < 15° 或 > 165°
- **单元间正确连接**：相邻单元的节点必须严格对齐（不允许"悬空节点"）

### 1.15.2 对称性利用

当结构和载荷满足对称条件时，可仅分析一半或四分之一模型：

- **镜像对称（Reflective symmetry）**：几何和载荷对称于一个面
- **旋转对称（Cyclic symmetry）**：结构在不同角度重复
- **轴对称（Axisymmetry）**：结构绕一个轴旋转对称

> 利用对称性可将问题规模减少 50%~75%，显著节省计算成本。

### 1.15.3 收敛性

FEM 解向精确解收敛的两种基本策略：

1. **h-收敛（h-convergence）**：保持单元阶次不变，减小单元尺寸（加密网格）
2. **p-收敛（p-convergence）**：保持网格不变，提高形函数多项式阶次（使用高阶单元）

> ⚠️ **难点**：h-收敛和 p-收敛的收敛速度不同。如果真解足够光滑，p-收敛的收敛速度可按指数级加速（谱收敛）；如果真解有奇异性（如裂纹尖端），h-收敛更有效因为加密网格可以更好地捕捉局部梯度。考试中常见的陷阱：加密网格后解仍不收敛→可能是因为单元阶次太低（需要 p-收敛）或者存在应力集中（需要局部 h-细化）。

### 1.15.4 结果检验

完成 FEM 分析后应检查：
1. **变形形状**：是否符合物理直觉（如约束点位移为零）
2. **外力平衡**：支反力之和是否等于外载荷之和
3. **数量级**：位移/应力是否在合理范围内

---

## 1.16 动力方程与非线性问题（概述）

### 1.16.1 结构动力方程

当惯性力和阻尼力不可忽略时，分时段静力方程推广为：

$$\mathbf{M}\ddot{\mathbf{D}} + \mathbf{C}\dot{\mathbf{D}} + \mathbf{K}\mathbf{D} = \mathbf{R}^{\text{ext}}$$

其中 $\mathbf{M}$ 为质量矩阵，$\mathbf{C}$ 为阻尼矩阵，$\mathbf{D}$ 为位移向量。

#### 逐步积分方法

**显式积分**（如中心差分法）：直接通过前两个时刻的解外推下一时刻；无需迭代，但时间步长受稳定条件限制。

**隐式积分**（如 Newmark-$\beta$ 法）：需在每个时间步求解联立方程组；无条件稳定（适当参数下），步长只取决于精度要求。

$$\mathbf{K}^{\text{eff}} = \frac{1}{\beta \Delta t^2}\mathbf{M} + \frac{\gamma}{\beta \Delta t}\mathbf{C} + \mathbf{K}$$

### 1.16.2 非线性问题分类

FEA 中三大类非线性：

| 类型 | 原因 | 示例 |
|------|------|------|
| **材料非线性** | 应力-应变关系不再满足 Hooke 定律 | 塑性变形、超弹性、蠕变 |
| **几何非线性** | 变形足够大，不能忽略几何变化对平衡的影响 | 大变形、屈曲、金属成形 |
| **接触非线性** | 不同部件之间的接触状态随加载而变化 | 碰撞、密封、装配分析 |

求解方法：**Newton-Raphson 迭代**（包括修正的半切线刚度法和**弧长法**处理 snap-through / snap-back 现象）。

> 💡 **理解关键**：三种非线性的区分——材料非线性是 D 矩阵不再恒定（σ=σ(ε) 非线性），几何非线性是 B 矩阵不再恒定（B=B(u) 依赖于位移），接触非线性是边界条件本身是未知的（接触状态待定）。每一种非线性对应 FEM 公式中不同的"变数"来源。

---

## 1.17 参考教材

1. Zienkiewicz O.C., Taylor, R.L. and Zhu, J.Z. *The Finite Element Method: Its Basis and Fundamentals* (7th ed.), Elsevier Butterworth-Heinemann, 2013.
2. Reddy, J.N. *An Introduction to the Finite Element Method*. McGraw-Hill, New York, 1993.
3. Cook. R.D., Malkus, D.S., Plesha, M.E. Witt, R.J. *Concepts and Applications of Finite Element Analysis* (4th ed.), John Wiley & Sons Inc. New York, 2002.
4. Bathe, K.J. *Finite Element Procedures*. Prentice Hall, Englewood Cliffs, NJ. 1996.
5. Rao, S.S. *The Finite Element Method in Engineering* (5th ed.) Elsevier Butterworth-Heinemann. 2011.
6. 王勖成，*有限单元法*，清华大学出版社，2003.
7. 朱伯芳，*有限单元法原理与应用*（第四版），水利水电出版社，2018.

---

## 检查你的理解

1. **变截面杆算例**：给定 $E = 200$ GPa，$A_1 = 100$ mm²，$A_2 = 50$ mm²，$A_3 = 25$ mm²，$l_1 = l_2 = l_3 = 100$ mm，$F_2 = 1000$ N，$F_3 = F_4 = 0$（自由端自由拉伸），$u_1 = 0$。请组装总体刚度矩阵，并求解 $u_2$、$u_3$、$u_4$。

2. **刚度矩阵性质**：为什么未施加边界条件的总体刚度矩阵是奇异的？给出物理原因和数学证据。

3. **形函数法推导链**：给定形函数 $N_i = (1-\xi)/2$、$N_j = (1+\xi)/2$，请从位移插值出发，经应变 → 应力 → 能量 → 变分，完整导出杆单元的刚度矩阵表达式。

4. **应变常数**：对于使用线性形函数的杆单元，为什么单元内的应变和应力是常数？这对网格划分有什么启示？

5. **显式 vs 隐式积分**：什么情况下选择显式时间积分，什么情况下选择隐式时间积分？

6. **收敛性**：h-收敛和 p-收敛有什么区别？如果加密网格后解仍在不断变化说明什么问题？

7. **对称性利用**：一段简支梁在中点承受集中力。请问可以利用何种对称性来减少计算量？边界条件应如何修改？

---

> **对应作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
