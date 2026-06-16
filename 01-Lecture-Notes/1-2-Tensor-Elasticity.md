# 张量分析与弹性力学

> **对应课件**：[`Chapter 2 Elastic theory.pdf`](../06-References/pdfs-originals/Chapter%202%20Elastic%20theory.pdf) 课程第 2 章 · [原文MD](../../md_output/Chapter%202%20Elastic%20theory.md)
> **章节定位**：Solid mechanics: theory of elasticity → Introduction to tensor → Theory of elasticity: differential context
> **相关作业**：[HW1 Q1-Q4](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：线性代数（向量、矩阵、行列式）、高等数学（偏导数、多元微积分）

---

> 🔗 **跨章连接**：本章是 FEM 课程的第 2 章，核心任务有两个——(1) 建立张量语言，这是后续所有公式推导的"通用语法"；(2) 导出弹性力学的 15 个方程，这是第 3 章变分法和第 4 章 FEM 求解的物理基础。学本章时重点关注"指标运算熟练度"，因为 HW1 全部是张量运算题，后面的 HW2-HW4 也会持续用到。

> **📋 考试范围覆盖**
>
> | 本讲义章节 | 考试大纲考点 |
> |-----------|-------------|
> | §2.1 指标与求和 | [Tensor] Conversion between algebra and tensors |
> | §2.2 常用符号 | [Tensor] Special symbols: Kronecker delta, Levi-Civita |
> | §2.3 线性代数复习 | [Tensor] Scalar and cross products; [Tensor] Vector/Matrix transformation |
> | §2.4 笛卡尔张量 | [Tensor] Definition (scalar, vector, tensor); [Tensor] Proving tensors |
> | §3.1 内力分析 | [Elasticity] Cauchy formula; [Elasticity] Stress transformation; [Elasticity] Principal stress/strain; [Elasticity] Invariant of stress/strain tensors |
> | §3.1.5 平衡方程 | [Elasticity] Equilibrium equation |
> | §3.2 变形分析 | [Elasticity] Geometric equation; [Elasticity] Green strain |
> | §3.3 本构关系 | [Elasticity] Constitutive relation; [Elasticity] Constitutive relation of isotropic linear elastic materials |
> | §3.4 弹性力学方程汇总 | [Elasticity] Boundary conditions; [Elasticity] Assumptions in deducing three sets of equations |

## 目录

- [1. 弹性力学发展简介](#1-弹性力学发展简介)
- [2. 张量引论（Introduction to Tensor）](#2-张量引论introduction-to-tensor)
  - [2.1 Index and summation（指标与求和）](#21-index-and-summation指标与求和)
  - [2.2 Common symbols（常用符号）](#22-common-symbols常用符号)
  - [2.3 Revision of linear algebra（线性代数复习）](#23-revision-of-linear-algebra线性代数复习)
  - [2.4 Cartesian tensor（笛卡尔张量）](#24-cartesian-tensor笛卡尔张量)
- [3. 弹性力学：微分表述](#3-弹性力学微分表述)
  - [3.1 内力分析（Internal Forces）](#31-内力分析internal-forces)
  - [3.2 变形分析（Deformation Analysis）](#32-变形分析deformation-analysis)
  - [3.3 本构关系（Constitutive Relations）](#33-本构关系constitutive-relations)
  - [3.4 弹性力学微分方程汇总](#34-弹性力学微分方程汇总)
  - [3.5 补充内容（Supplementary）](#35-补充内容supplementary)
  - [3.6 结论](#36-结论)

---

## 1. 弹性力学发展简介

弹性力学是一门具有 300 多年历史的学科。弹性（elasticity）定义为：撤去外力后，弹性材料制成的物体会恢复到其原始状态。

- **1678 年**：Hooke 发现弹性体的位移与外力成正比。
- **1821-1823 年**：Navier 和 Cauchy 分别推导出线弹性边值问题的控制平衡方程（governing equilibrium equations），首次揭示了弹性体力学行为的数学描述。
- **1855 年**：Saint-Venant 发表关于扭转（torsion）和弯曲（bending）的经典工作。
- **1933 年**：Muskhelishvili 提出复变函数方法（complex variable formulation），可以处理矩形截面扭转等问题。
- **20 世纪末**：弹性力学进一步发展到考虑与其他物理因素（热、粘性等）的耦合作用。


### 弹性力学 vs. 材料力学

| 对比维度 | 材料力学 | 弹性力学 |
|---------|--------|--------|
| 研究对象 | 主要针对杆、梁等杆件结构 | 无几何形状限制（杆、板、壳、实体） |
| 基本假设 | 较多简化假设（如"平截面假定"） | 更少的人为假设，更一般化 |
| 数学工具 | 初等微积分 | 偏微分方程、张量分析、复变函数等 |
| 解的精度 | 近似解 | 精确解（对可解问题）或系统化的近似 |

> 💡 **理解关键**：弹性力学 = 材料力学的"升级版"。材料力学用简化假设减少未知量（杆件→一维），弹性力学放弃很多假设后未知量暴增到 15 个，必须引入张量工具来维持表达简洁。

弹性力学需要强有力的数学工具，而**张量分析（tensor analysis）**正是描述弹性力学问题的最简洁、最系统的语言。张量形式的优势有三：
1. **简化书写**——一组方程可压缩为一个指标方程；
2. **揭示物理本质**——张量方程在坐标变换下形式不变（坐标不变性），直接反映了物理规律的客观性；
3. **便于文献阅读**——现代连续介质力学文献普遍采用张量表达。

---

## 2. 张量引论（Introduction to Tensor）

张量形式的核心思想：**物理定律不依赖于坐标系的选择**。在坐标系 $K$ 下某物理量分量为 $a,b,c,\ldots$ 的物理定律，在坐标系 $K'$ 下分量变为 $a',b',c',\ldots$，但定律的形式保持不变。张量为这种坐标变换下的物理量分量提供了清晰简洁的表达。

> 💡 **理解关键**：张量的核心不是什么高深数学，而是一个很朴素的思想——"物理定律不应该因为你换了个坐标系就变了"。就像你站在房间里测身高，不管你是面朝东还是面朝北，身高数值、物理规律都不应该变。张量就是确保这一点的数学工具。

我们按以下顺序引入张量：
1. 指标与求和记号
2. 两个常用符号：$\delta_{ij}$ 和 $e_{ijk}$
3. 线性代数中的坐标变换回顾
4. 基于坐标变换严格定义标量、向量和张量

### 2.1 Index and summation（指标与求和）

#### Einstein 求和约定

在张量运算中，**在同一项中重复出现的下标表示对该下标求和**（哑指标，dummy index）。求和范围通常为 1 到 3（三维空间）：

$$a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^3 a_i b_i$$

> ⚠️ **重难点**：Einstein 求和约定是整个张量运算的基础语法，初学最容易混淆的是"哪些指标求和、哪些指标不求和"。一条铁律——**同一项中出现两次的字母自动求和（哑指标），只出现一次的字母是自由指标代表一组方程**。出现三次或以上？表达式无意义。

#### 自由指标（Free index）

在单项中只出现一次的指标，代表一般的分量子编号（可取 1, 2, 3）。方程两边的自由指标必须一致。

**例 1**：
$$a_{ij} b_j = a_{i1}b_1 + a_{i2}b_2 + a_{i3}b_3 \quad (i = 1,2,3)$$

其中 $i$ 是自由指标（表示 3 个方程），$j$ 是哑指标（对 $j=1,2,3$ 求和）。

> ❌ **易错点**：看表达式 $a_{ij}b_j = c_i$。学生常犯的错误是以为 $j$ 是自由指标所以方程中还有 $j$。正确理解：$j$ 已经被求和掉了（消掉了），只剩下 $i$ ——这个方程是 3 个（$i=1,2,3$）。核心判断：**哑指标在等号右边看不见**。

**例 2**：线性方程组的简写
$$\left\{\begin{array}{l} a_{11}b_1 + a_{12}b_2 + a_{13}b_3 = c_1 \\ a_{21}b_1 + a_{22}b_2 + a_{23}b_3 = c_2 \\ a_{31}b_1 + a_{32}b_2 + a_{33}b_3 = c_3 \end{array}\right. \quad\Longrightarrow\quad a_{ij}b_j = c_i \quad (i=1,2,3)$$

#### 求和指标的可替换性

哑指标可以使用任意字母替换，不影响求和结果：
$$a_{ii} = a_{11} + a_{22} + a_{33} = a_{kk}$$

> **规则**：在单向中同一哑指标最多出现两次。出现两次代表求和；出现一次是自由指标；出现超过两次则该表达式无意义。

> 💡 **理解关键**：哑指标可替换性就像积分变量——$\int f(x)dx$ 和 $\int f(t)dt$ 是一样的。$a_{ii}$ 和 $a_{kk}$ 完全相同。这个性质在推导中极其有用：当你有两个求和项需要合并时，可能需要将某个哑指标换成另一字母以避免冲突。

### 2.2 Common symbols（常用符号）

#### Kronecker 符号 $\delta_{ij}$

$\delta_{ij}$ 有 9 个分量，定义如下：

$$\boxed{\delta_{ij} = \begin{cases} 1, & i=j \\ 0, & i \neq j \end{cases}}$$

矩阵形式为 $3 \times 3$ 的单位矩阵：
$$\delta_{ij} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$$


**重要性质**：
- **迹**：$\delta_{ii} = \delta_{11} + \delta_{22} + \delta_{33} = 3$
- **指标替换（substitution property）**：$\delta_{ij}A_i = A_j$
  - 推导：$\delta_{ij}A_i = \delta_{1j}A_1 + \delta_{2j}A_2 + \delta_{3j}A_3$，当 $j=1$ 时只剩 $A_1$，$j=2$ 时只剩 $A_2$，$j=3$ 时只剩 $A_3$。所以 $\delta_{ij}$ 将 $A_i$ 中的下标 $i$ 替换为 $j$。
- $\delta_{ij}\delta_{jk} = \delta_{ik}$
- **正交基条件**：标准正交基向量满足 $\mathbf{e}_i \cdot \mathbf{e}_j = \delta_{ij}$

**应用示例**：将方程组
$$\begin{cases} a_{11} = b + c_{11} \\ a_{22} = b + c_{22} \\ a_{33} = b + c_{33} \\ a_{12} = c_{12},\; a_{21} = c_{21} \\ a_{23} = c_{23},\; a_{32} = c_{32} \\ a_{31} = c_{31},\; a_{13} = c_{13} \end{cases}$$
简化为一项：
$$a_{ij} = b\delta_{ij} + c_{ij}$$

> 💡 **理解关键**：上面的例子很好地展示了 $\delta_{ij}$ 的本质功能——它是一个"二分类处理器"：当 $i=j$ 时加 $b$，当 $i \neq j$ 时不加。这是后来 Hooke 定律 $\sigma_{ij} = \lambda e \delta_{ij} + 2\mu \varepsilon_{ij}$ 中 $\lambda e \delta_{ij}$ 项的来源：$\delta_{ij}$ 确保只有对角线（$i=j$）上涉及体积变形的影响。

#### 置换符号（Permutation symbol）$e_{ijk}$

$e_{ijk}$ 有 27 个分量，也常称为 Levi-Civita 符号：

$$\boxed{e_{ijk} = \begin{cases}
0 & \text{任意两指标相等} \\
+1 & (i,j,k)\text{ 为偶排列：}(1,2,3),(2,3,1),(3,1,2)\text{（顺时针循环）} \\
-1 & (i,j,k)\text{ 为奇排列：}(1,3,2),(2,1,3),(3,2,1)\text{（逆时针循环）}
\end{cases}}$$

> ⚠️ **重难点**：$e_{ijk}$ 有三个指标共 27 个分量，不是矩阵。记忆方法：把 $(1,2,3)$ 的循环看作"顺时针"，偶排列（+1）就是顺时针方向数过去——$(1,2,3)$, $(2,3,1)$, $(3,1,2)$；奇排列（-1）是逆时针——$(1,3,2)$, $(2,1,3)$, $(3,2,1)$。任意两个指标相等时为 0（所以非零分量只有 6 个）。

**核心恒等式**（$\delta$-$e$ 关系）：

$$\boxed{e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}}$$


由此可推出：
- $e_{ijk}e_{ijt} = 2\delta_{kt}$
- $e_{ijk}e_{ijk} = 6$

**向量叉积的张量形式**：
$$(\mathbf{a}\times\mathbf{b})_i = e_{ijk}a_j b_k$$

即展开为行列式：
$$\mathbf{a} \times \mathbf{b} = \begin{vmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \\ a_1 & a_2 & a_3 \\ b_1 & b_2 & b_3 \end{vmatrix} = e_{ijk}a_i b_j \mathbf{e}_k$$

> 🔗 **跨章连接**：向量叉积的张量形式是后续推导 $\delta$-$e$ 恒等式时反复用到的工具。HW1 Q3 要求证明 $e_{ijk}$ 是三阶张量，解题时需要用行列式性质 $\det(\mathbf{L}) = 1$（右手正交变换）。$\delta$-$e$ 恒等式的推导也出现在 HW1 Q1。

**行列式的 $e_{ijk}$ 表示**：$3 \times 3$ 矩阵的行列式可表示为
$$a = |a_{ij}| = e_{ijk}a_{i1}a_{j2}a_{k3} = e_{ijk}a_{1i}a_{2j}a_{3k}$$

更一般地，列/行交换与行列式符号的关系：
$$a e_{rst} = e_{ijk}a_{ir}a_{js}a_{kt} \quad \text{（列交换）}$$
$$a e_{rst} = e_{ijk}a_{ri}a_{sj}a_{tk} \quad \text{（行交换）}$$

**本小节恒等式速查表**：

| 恒等式 | 用途 |
|-------|------|
| $\delta_{ij}A_i = A_j$ | 指标替换 |
| $\delta_{ii} = 3$ | 求迹 |
| $e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$ | $\delta$-$e$ 核心恒等式 |
| $e_{ijk}e_{ijt} = 2\delta_{kt}$ | 双重缩并 |
| $e_{ijk}e_{ijk} = 6$ | 三重缩并 |

---

> #### 检查你的理解：2.1-2.2
> 1. 在表达式 $a_{ij}b_{jk} = c_{ik}$ 中，哪些是自由指标，哪些是哑指标？它代表了多少个方程？
> 2. $\delta_{ij}\delta_{ij}$ 等于多少？为什么？
> 3. 利用 $e_{ijk}$ 写出 $( \mathbf{a} \times \mathbf{b} ) \cdot ( \mathbf{c} \times \mathbf{d} )$ 的张量形式并用 $\delta$-$e$ 恒等式证明它在坐标变换下不变。

---

### 2.3 Revision of linear algebra（线性代数复习）

#### 向量与坐标

在笛卡尔坐标系 $\{O; x_1,x_2,x_3\}$ 中，基向量 $\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3$ 是三个相互垂直的单位向量（标准正交基）。任意向量 $\mathbf{r}$ 可表示为：

$$\mathbf{r} = x_1\mathbf{e}_1 + x_2\mathbf{e}_2 + x_3\mathbf{e}_3 = x_i \mathbf{e}_i$$

向量 $\mathbf{r}$ 本身**不依赖**于坐标系的选择，但其分量 $x_i$ 在不同坐标系下取不同的值。

> 💡 **理解关键**：这里有一个深刻的物理直觉——向量本身是一个几何量（一根有长度有方向的箭头），不管你怎么旋转坐标系，那根箭头本身没变。变的只是它在各坐标轴上的"投影数值"（分量）。张量分析把这种"量本身不变，分量按规律变"的思想推广到更高阶的量。

#### 基向量变换

设新坐标系 $\{O; x_1',x_2',x_3'\}$，记方向余弦 $l_{i'j} = \cos(x_i', x_j)$（新坐标轴 $x_i'$ 与旧坐标轴 $x_j$ 之间的夹角余弦）。

新基向量用旧基向量表达：
$$\mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j \quad\text{矩阵形式}\quad \mathbf{e}' = \mathbf{L}\mathbf{e}$$

其中 $\mathbf{L} = [l_{i'j}]$ 是 $3 \times 3$ 方向余弦矩阵。

方向余弦表：

| | $x_1$ | $x_2$ | $x_3$ |
|--|------|------|------|
| $x_1'$ | $l_{1'1}$ | $l_{1'2}$ | $l_{1'3}$ |
| $x_2'$ | $l_{2'1}$ | $l_{2'2}$ | $l_{2'3}$ |
| $x_3'$ | $l_{3'1}$ | $l_{3'2}$ | $l_{3'3}$ |

$\mathbf{L}$ 是**正交矩阵**（orthogonal matrix），满足 $\mathbf{L}^{-1} = \mathbf{L}^T$。因此逆变换为：

$$\mathbf{e} = \mathbf{L}^T\mathbf{e}' \quad\text{即}\quad \mathbf{e}_j = l_{ij'}\mathbf{e}_i'$$

> ⚠️ **重难点**：为什么 $\mathbf{L}$ 是正交阵（$\mathbf{L}^T = \mathbf{L}^{-1}$）？因为新基向量 $\mathbf{e}_1', \mathbf{e}_2', \mathbf{e}_3'$ 也是标准正交基，所以 $\mathbf{e}_i' \cdot \mathbf{e}_j' = \delta_{ij}$。代入 $\mathbf{e}_i' = l_{ik}\mathbf{e}_k$ 可得 $l_{ik}l_{jk} = \delta_{ij}$，即 $\mathbf{L}\mathbf{L}^T = \mathbf{I}$。另一种理解：方向余弦矩阵的每一行是新坐标轴在旧坐标系中的方向分量，行向量是单位向量且两两正交→正交矩阵。

正交性的分量形式：
$$\boxed{l_{i'k}l_{j'k} = \delta_{i'j'}} \qquad \boxed{l_{k'i}l_{k'j} = \delta_{ij}}$$

> ❌ **易错点**：注意方向余弦指标中 $l_{i'j}$ 的指标顺序——$i'$（新坐标轴）是第一个指标，$j$（旧坐标轴）是第二个。初学者经常搞混旧到新和新到旧的公式，其实记住一个就好：$a_{i'} = l_{i'j}a_j$（新=方向余弦乘以旧并求和）。$i'$ 是新指标、$j$ 是旧指标、$j$ 是哑指标被消掉——用"新旧指标配对规则"核对：哑指标（j）在右边出现，自由指标（i'）在等号左边和右边都出现。

#### 向量分量变换

**旧到新**：
$$\boxed{a_{i'} = l_{i'j}a_j} \quad (i' = 1,2,3)$$

**新到旧**：
$$\boxed{a_i = l_{ij'}a_{j'}} \quad (i = 1,2,3)$$

**推导**：由向量不变性 $\mathbf{a} = a_i\mathbf{e}_i = a_{i'}\mathbf{e}_{i'}'$，代入 $\mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j$：
$$a_i\mathbf{e}_i = a_{i'}l_{i'j}\mathbf{e}_j = a_{i'}l_{i'i}\mathbf{e}_i \;\Rightarrow\; a_i = l_{i'i}a_{i'}$$
取逆即得 $a_{i'} = l_{i'j}a_j$。

#### 矩阵变换（线性变换的矩阵表示）

设线性变换算子 $\Psi$ 在某基下的矩阵表示为 $\mathbf{A}$：
$$\Psi\mathbf{a} = \mathbf{e}^T\mathbf{A}\mathbf{a}$$

在另一组基下：
$$\Psi\mathbf{a} = \mathbf{e}'^T\mathbf{A}'\mathbf{a}'$$

代入基变换 $\mathbf{e}^T = \mathbf{e}'^T\mathbf{L}$ 和分量变换 $\mathbf{a} = \mathbf{L}^T\mathbf{a}'$：
$$\mathbf{e}^T\mathbf{A}\mathbf{a} = \mathbf{e}'^T\mathbf{L}\mathbf{A}\mathbf{L}^T\mathbf{a}'$$

与 $\mathbf{e}'^T\mathbf{A}'\mathbf{a}'$ 对比，由于 $\mathbf{a}'$ 任意，得：

$$\boxed{\mathbf{A}' = \mathbf{L}\mathbf{A}\mathbf{L}^T}$$

用指标记号写出：
$$\boxed{a_{i'j'} = a_{mn}l_{i'm}l_{j'n}}$$


> 这就是二阶张量分量变换律的来源！可以说，**满足此变换律的矩阵就是二阶张量**。

#### 标量积与向量积

**标量积（点积）**：
$$\boxed{\mathbf{a}\cdot\mathbf{b} = a_i b_i}$$

**正交基条件**：$\mathbf{e}_i \cdot \mathbf{e}_j = \delta_{ij}$，由此可进一步证明
$$\boxed{l_{k'i}l_{k'j} = \delta_{ij}}$$

**向量积（叉积）**：
$$\boxed{(\mathbf{a}\times\mathbf{b})_i = e_{ijk}a_j b_k}$$

#### $\delta_{ij}$-$e_{ijk}$ 恒等式的证明

利用向量双重叉积公式 $\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = (\mathbf{a}\cdot\mathbf{c})\mathbf{b} - (\mathbf{a}\cdot\mathbf{b})\mathbf{c}$ 来证明 $e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}$。

**证明步骤**：

1. 将 $\mathbf{a} = a_i\mathbf{e}_i$，$\mathbf{b} = b_k\mathbf{e}_k$，$\mathbf{c} = c_s\mathbf{e}_s$ 代入双重叉积：
   $$\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = a_i\mathbf{e}_i \times (b_k c_s e_{ksp}\mathbf{e}_p) = a_i b_k c_s e_{ksp} e_{ipj} \mathbf{e}_j$$

2. 将相同向量代入右侧 $(\mathbf{a}\cdot\mathbf{c})\mathbf{b} - (\mathbf{a}\cdot\mathbf{b})\mathbf{c}$：
   $$= a_i c_s \delta_{is} b_k \mathbf{e}_k - a_i b_k \delta_{ik} c_s \mathbf{e}_s$$
   $$= a_i b_k c_s (\delta_{is}\delta_{jk} - \delta_{ik}\delta_{js})\mathbf{e}_j$$

3. 两侧对比，由于 $a_i, b_k, c_s$ 为任意向量，得：
   $$e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}$$

> 💡 **理解关键**：这个证明是 HW1 Q1 的解题思路。核心手法是"对任意向量成立 → 系数相等"——这是张量推导的常见策略。证明中 $e_{ksp}e_{ipj}$ 的指标排列可能与讲义中的 $e_{ijk}e_{ist}$ 看起来不同，但本质一样——只是哑指标字母换了。关键是确认第一个指标（$i$ 和 $k$）相同，然后用恒等式。

---

> #### 检查你的理解：2.3
> 1. 证明正交矩阵的逆等于其转置：从正交基条件 $\mathbf{e}_i'\cdot\mathbf{e}_j' = \delta_{ij}$ 出发，利用基向量变换推导 $\mathbf{L}\mathbf{L}^T = \mathbf{I}$。
> 2. 已知在旧坐标系中某向量的分量为 $(a_1,a_2,a_3)$，新坐标系 $x_1'$ 轴沿 $(1,1,0)/\sqrt{2}$ 方向，$x_3'$ 轴与 $x_3$ 重合。求方向余弦矩阵 $\mathbf{L}$ 及各 $a_{i'}$。

---

### 2.4 Cartesian tensor（笛卡尔张量）

张量是描述物理状态或几何性质的量。"Tensor"一词源于其与应力（tensile）的历史关联。实际上，经典力学中刚体的转动惯量（moment of inertia）就是二阶张量——这一点在普通力学教材中常常未被明确指出。

> 💡 **理解关键**：张量的定义方式从之前的"会什么运算"升级为"满足什么变换律"。这是一个重要的思维转折——张量不是由其分量个数定义的，而是由其在坐标变换下的行为定义的。$3 \times 3$ 矩阵恰好满足变换律才是二阶张量。

我们将从坐标变换律的角度严格定义标量（零阶张量）、向量（一阶张量）和二阶张量。

#### 标量（零阶张量）

**定义**：由一个实数表达的物理量或几何量，在坐标变换下保持不变。

普通的标量（ordinary scalar）如密度、温度，在任何坐标系下值相同。这类标量**不依赖**坐标系选择。

**验证示例：两点间距离**

设空间两点 A 和 B，在 $\{x_i\}$ 系下坐标为 $a_i$ 和 $b_i$，距离平方：
$$\Delta s^2 = \Delta x_i \Delta x_i, \quad \Delta x_i = a_i - b_i$$

在新坐标系 $\{x_i'\}$ 中：
$$\Delta x_{i'} = a_{i'} - b_{i'} = l_{i'j}a_j - l_{i'j}b_j = l_{i'j}(a_j - b_j) = l_{i'j}\Delta x_j$$

则：
$$(\Delta s')^2 = \Delta x_{i'}\Delta x_{i'} = l_{i'k}l_{i'l}\Delta x_k\Delta x_l = \delta_{kl}\Delta x_k\Delta x_l = \Delta x_k\Delta x_k = (\Delta s)^2$$

故距离是标量。

#### 向量（一阶张量）

**严格定义**：由 3 个标量定义、依赖坐标系选择、且满足坐标分量变换律的量。

$$a_{i'} = l_{i'j}a_j$$

> 注意：并非任意三个标量的组合都是向量。例如（年龄，身高，体重）可组成 3 个数，尽管可以加法和数乘，但它不是向量——因为这 3 个数在坐标变换下不满足上述变换律。


**验证示例：位移是向量**

点 P 的坐标 $x_i = x_i(t)$，在时间间隔 $(t, t+\Delta t)$ 内位移：
$$u_i = x_i(t+\Delta t) - x_i(t)$$

在新坐标系下：
$$u_{i'} = x_{i'}(t+\Delta t) - x_{i'}(t) = l_{i'j}x_j(t+\Delta t) - l_{i'j}x_j(t) = l_{i'j}u_j$$

满足向量变换律。类似可证速度和加速度也都是向量。

#### 二阶张量

**严格定义**：由 9 个标量定义、依赖坐标系选择、且满足坐标分量变换律的量：

$$\boxed{a_{i'j'} = a_{mn}l_{i'm}l_{j'n}}$$

可写为矩阵形式：
$$a_{ij} = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}$$

> 与向量的情况类似，并非任意 $3\times 3$ 矩阵都是张量——必须满足上述坐标变换律。

> 🔗 **跨章连接**：二阶张量的变换律 $a_{i'j'} = a_{mn}l_{i'm}l_{j'n}$ 是整个弹性力学的核心公式。后面会反复用到：应力变换 $\sigma_{mn}' = \sigma_{ij}l_{mi}l_{nj}$、惯性矩变换、本构张量变换都是同一形式。记忆技巧——"每个新指标对应一个 $l$，旧指标缩并求和"。

**张量阶数与分量数**：$n$ 阶张量在 $d$ 维空间中具有 $d^n$ 个分量。
- 二维向量：$2^1 = 2$ 个分量
- 三维向量：$3^1 = 3$ 个分量
- 三维二阶张量：$3^2 = 9$ 个分量
- 三维三阶张量（如 $e_{ijk}$）：$3^3 = 27$ 个分量

**验证示例 1：Kronecker 符号是二阶张量**

$$
\delta_{i'j'} = \mathbf{e}_{i'}'\cdot\mathbf{e}_{j'}' = l_{i'k}\mathbf{e}_k \cdot l_{j's}\mathbf{e}_s = l_{i'k}l_{j's}\delta_{ks}
$$

满足二阶张量变换律。

> 💡 **理解关键**：$\delta_{ij}$ 的特殊性——它在任何直角坐标系下都是单位矩阵。因为 $\delta_{i'j'} = l_{i'k}l_{j'k} = \delta_{i'j'}$（正交性条件），所以 $\delta_{ij}$ 是"各向同性张量"——分量值不随坐标系旋转而改变。$\delta_{ij}$ 和 $e_{ijk}$ 是仅有的两个各向同性张量（前者是二阶，后者是三阶）。

**验证示例 2：转动惯量是二阶张量**

设刚体以角速度 $\boldsymbol{\omega}$ 旋转，任意点的速度 $\mathbf{v} = \boldsymbol{\omega} \times \mathbf{r}$（其中 $\mathbf{r}$ 是该点的位置向量）。

刚体的角动量：
$$\mathbf{I} = \int \mathbf{r} \times \mathbf{v}\, dm = \int \mathbf{r} \times (\boldsymbol{\omega} \times \mathbf{r})\, dm$$

用张量展开：$I_i = \int (x_p\mathbf{e}_p) \times [(\omega_j\mathbf{e}_j) \times (x_k\mathbf{e}_k)] dm$

$$= \int (x_p\mathbf{e}_p) \times (\omega_j x_k e_{jkq}\mathbf{e}_q) dm = \left(e_{pqi}e_{jkq}\int x_p x_k dm\right) \omega_j \mathbf{e}_i$$

转动惯量张量：
$$I_{ij} = e_{pqi}e_{jkq}\int x_p x_k dm$$

利用 $\delta$-$e$ 恒等式：$e_{qip}e_{qjk} = \delta_{ij}\delta_{pk} - \delta_{ik}\delta_{pj}$

$$I_{ij} = \int (\delta_{ij}\delta_{pk} - \delta_{ik}\delta_{pj})x_p x_k dm = \int (\delta_{ij}x_p x_p - x_i x_j) dm$$

在新坐标系下：
$$I_{i'j'} = \int (\delta_{i'j'}x_{p'}x_{p'} - x_{i'}x_{j'}) dm$$

利用 $\delta_{i'j'} = l_{i'm}l_{j'n}\delta_{mn}$ 和 $x_{i'}x_{j'} = l_{i'm}l_{j'n}x_m x_n$：
$$I_{i'j'} = l_{i'm}l_{j'n}\int (\delta_{mn}x_p x_p - x_m x_n) dm = l_{i'm}l_{j'n}I_{mn}$$

满足二阶张量变换律。Q.E.D.

> 💡 **理解关键**：转动惯量的张量证明是一个完整的"示范性推导"——它展示了如何用 $\delta$-$e$ 恒等式化简复杂表达式，以及如何验证一个物理量满足张量变换律。这是 HW1 Q2（证明 $\varepsilon_{ij}$ 是张量）的解题范本。

#### 张量运算

1. **加法与减法**：同阶张量对应分量相加
   $$C_{ij} = A_{ij} + B_{ij}$$

2. **张量与标量的乘积**：每个分量乘以该标量
   $$C_{ij} = \lambda A_{ij}$$

3. **张量积（外积）**：阶数相加
   $$C_{ijkl} = A_{ij}B_{kl} \quad \text{（二阶 $\times$ 二阶 = 四阶）}$$

4. **缩并（Contraction）**：令两个指标相同并求和，阶数降 2
   $$A_{ii} = A_{11} + A_{22} + A_{33} = \text{tr}(\mathbf{A})$$

5. **对称/反对称分解**：
   $$A_{ij} = \underbrace{\frac{1}{2}(A_{ij}+A_{ji})}_{\text{对称部分}} + \underbrace{\frac{1}{2}(A_{ij}-A_{ji})}_{\text{反对称部分}}$$

> 💡 **理解关键**：对称/反对称分解是后面应变分析中位移梯度分解 $u_{i,j} = \varepsilon_{ij} + \omega_{ij}$ 的数学基础。对称部分对应"伸长/缩"（→应力），反对称部分对应"刚体转动"（→不产生应力）。

6. **张量不变量**：与矩阵的特征值不变量类似
   $$I_1 = A_{ii}, \quad I_2 = \frac{1}{2}(A_{ii}A_{jj} - A_{ij}A_{ji}), \quad I_3 = \det(\mathbf{A})$$

> 🔗 **跨章连接**：张量不变量在后面应力分析中重新出现为"应力不变量 $I_1, I_2, I_3$"（§3.1.3），在应变分析中还有"应变不变量"。不变量的物理意义是——它们不随坐标系旋转而改变，是刻画应力/应变状态本质特征的特征值。

7. **张量的导数记号**：以下逗号记法广泛用于弹性力学
   $$\frac{\partial (\cdot)}{\partial x_i} \equiv (\cdot)_{,i} \qquad \frac{\partial^2 (\cdot)}{\partial x_i \partial x_j} \equiv (\cdot)_{,ij}$$

> ❌ **易错点**：逗号记法虽然简洁，但初学者容易漏看。$\sigma_{ij,j}$ 不是 $\sigma_{ij}$ 乘上什么——而是 $\frac{\partial \sigma_{ij}}{\partial x_j}$（对 $j$ 求偏导并求和）。$uv_{,i}$ 是 $\frac{\partial (uv)}{\partial x_i}$，不是 $u \cdot v_{i}$。读公式时养成"看到逗号就默读为对逗号后面指标求导"的习惯。

---

> #### 检查你的理解：2.4
> 1. "任意 $3 \times 3$ 矩阵都是二阶张量"这个说法对吗？为什么？
> 2. 证明置换符号 $e_{ijk}$ 在右手坐标系间的正交变换下是三阶张量（提示：行列式 $\det(\mathbf{L}) = 1$）。
> 3. 已知 $A_{ij}$ 和 $B_{kl}$ 分别是二阶张量，证明 $C_{ijkl} = A_{ij}B_{kl}$ 是四阶张量。

---

## 3. 弹性力学：微分表述

以下进入弹性力学的主体部分。基本框架由三部分组成：

1. **内力分析**——应力状态、Cauchy 公式、坐标变换、主应力与不变量、应力张量分解、平衡方程
2. **变形分析**——构形与位移、Green 应变至小变形应变、应变物理意义、刚体变形与旋转张量
3. **本构理论**——热力学框架下的应变能函数、线弹性本构推导（81 $\to$ 36 $\to$ 21 $\to$ 2）

> 💡 **理解关键**：弹性力学的逻辑链条非常清晰——内力（应力）→ 变形（应变）→ 应力-应变关系（本构）。三组方程各管一摊，合在一起刚好 15 个方程 15 个未知量，封闭可解。这一章的任务就是导出这三组方程。

---

### 3.1 内力分析（Internal Forces）

#### 3.1.1 一点的应力状态

物体内一点的应力状态由**九个应力分量**完全确定，前提是：任意斜截面上的应力可由这几个分量表示。这一结论在弹性力学中被严格证明。

考虑右图所示的无限小四面体微元，其三个面分别平行于坐标平面，第四个面为任意方向的斜面（外法向为 $\mathbf{n}$）。

**Cauchy 应力公式的推导**：

设斜面的外法向 $\mathbf{n}$ 的方向余弦为 $l_1 = \cos(\mathbf{n}, x_1), l_2 = \cos(\mathbf{n}, x_2), l_3 = \cos(\mathbf{n}, x_3)$。斜面面积记为 $dA$，则三个坐标面上的面积分别为 $l_1 dA, l_2 dA, l_3 dA$。

考虑四面体在三个坐标轴方向的力平衡。以 $x_1$ 方向为例：

$$T_1^n dA = \sigma_{11}(l_1 dA) + \sigma_{21}(l_2 dA) + \sigma_{31}(l_3 dA)$$

消去 $dA$，得：

$$\boxed{T_i^n = \sigma_{ji}l_j = \sigma_{ij}l_j} \quad (i = 1,2,3)$$

> ⚠️ **重难点**：Cauchy 公式的推导前提是"四面体微元无限小"→ 体力（体积力 $\propto$ 体积）远小于面力（面积力 $\propto$ 面积）→ 体力项在力平衡中被忽略。这是连续介质力学中"微元体"方法的典型手法。另一个陷阱——注意应力分量的指标顺序：$\sigma_{21}$ 表示作用在法向为 $x_2$ 的面上、沿 $x_1$ 方向的应力分量。$\sigma_{ij}$ 的第一个指标是面的法向，第二个指标是力的方向。

> **物理含义**：任意斜面上的应力向量分量 $T_i^n$ 等于应力张量 $\sigma_{ij}$ 与斜面法向方向余弦 $l_j$ 的缩并。这就是著名的 **Cauchy 应力公式**。

**展开形式**：
$$\begin{pmatrix} T_1^n \\ T_2^n \\ T_3^n \end{pmatrix} = \begin{pmatrix} \sigma_{11} & \sigma_{21} & \sigma_{31} \\ \sigma_{12} & \sigma_{22} & \sigma_{32} \\ \sigma_{13} & \sigma_{23} & \sigma_{33} \end{pmatrix} \begin{pmatrix} l_1 \\ l_2 \\ l_3 \end{pmatrix}$$

> ❌ **易错点**：Cauchy 公式展开矩阵中，矩阵的行是应力矩阵的**行**（即 $i$ 指标），列是 $l_j$ 的列。注意矩阵是 $\sigma_{ij}$（$i$ 行 $j$ 列），$T_i^n = \sigma_{ij}l_j$ 是 $i$ 行与 $j$ 列求和的缩并。如果用 Python/Matlab 写代码，`T = sigma @ l` 即可（不要写成 `l @ sigma`）。

**剪应力互等定理**：考虑微元体的力矩平衡可得
$$\sigma_{12} = \sigma_{21}, \quad \sigma_{23} = \sigma_{32}, \quad \sigma_{31} = \sigma_{13}$$

即应力张量是对称的：$\sigma_{ij} = \sigma_{ji}$。这使独立的应力分量从 9 个减少到 6 个。


#### 3.1.2 应力坐标变换

当坐标轴旋转时，应力分量在新坐标系下的表达式可利用 Cauchy 公式导出。

设新坐标系 $\{x_1', x_2', x_3'\}$，在旧坐标系中取一斜面，其法向为新系的 $x_n'$ 轴。斜面在新系的 $x_m'$ 方向上的投影应力 $\sigma_{nm}$ 为：

$$\sigma_{nm} = T_i^n l_{mi} = \sigma_{ji}l_{nj}l_{mi}$$

对任意指标 $m,n = 1,2,3$，即得新系下的 $9$ 个应力分量：

$$\boxed{\sigma_{mn}' = \sigma_{ij} l_{mi} l_{nj}}$$

**矩阵形式**：
$$\boxed{\boldsymbol{\sigma}' = \mathbf{L}\boldsymbol{\sigma}\mathbf{L}^T}$$

其中 $\mathbf{L}$ 为方向余弦矩阵。展开写：

$$\boldsymbol{\sigma}' = \begin{bmatrix} l_{11} & l_{12} & l_{13} \\ l_{21} & l_{22} & l_{23} \\ l_{31} & l_{32} & l_{33} \end{bmatrix} \begin{bmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\ \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{bmatrix} \begin{bmatrix} l_{11} & l_{21} & l_{31} \\ l_{12} & l_{22} & l_{32} \\ l_{13} & l_{23} & l_{33} \end{bmatrix}$$

> 应力分量按二阶张量的变换律进行坐标变换。这是应力是二阶张量的有力证据。


#### 3.1.3 主应力与应力不变量

**主应力的物理意义**：应力张量在某方向截面上的应力向量只有正应力分量而无剪应力分量时，该正应力称为主应力（principal stress），该截面称为主平面（principal plane），其法向称为主方向（principal direction）。

根据 Cauchy 公式，主应力 $\sigma$ 和主方向 $n_i$ 满足：
$$T_i^n = \sigma_{ij}n_j = \sigma n_i \quad \text{即} \quad \sigma_{ij}n_j - \sigma n_i = 0$$

利用 $\delta_{ij}$ 替换：$(\sigma_{ij} - \sigma\delta_{ij})n_j = 0$

这是矩阵特征值问题。非零解的条件是特征行列式为零：
$$\boxed{|\sigma_{ij} - \sigma\delta_{ij}| = 0}$$

展开得三次方程（特征多项式）：
$$\boxed{\sigma^3 - I_1\sigma^2 + I_2\sigma - I_3 = 0}$$

其中 $\sigma_1, \sigma_2, \sigma_3$ 是三个主应力（特征值），$I_1, I_2, I_3$ 分别为第一、第二、第三**应力不变量**（stress invariants）。

> 💡 **理解关键**：主应力 = 应力矩阵的特征值，主方向 = 应力矩阵的特征向量。这完全等同于线性代数中的矩阵对角化问题。$I_1, I_2, I_3$ 就是特征多项式的系数（只是符号有正负），验证很简单——展开 $\det(\lambda I - A) = \lambda^3 - \text{tr}(A)\lambda^2 + \ldots$ 对比即可。

**应力不变量的表达式**：

$$\boxed{I_1 = \sigma_{11} + \sigma_{22} + \sigma_{33} = \sigma_{ii} = \sigma_1 + \sigma_2 + \sigma_3}$$

$$\boxed{I_2 = \begin{vmatrix} \sigma_{11} & \sigma_{12} \\ \sigma_{21} & \sigma_{22} \end{vmatrix} + \begin{vmatrix} \sigma_{22} & \sigma_{23} \\ \sigma_{32} & \sigma_{33} \end{vmatrix} + \begin{vmatrix} \sigma_{33} & \sigma_{31} \\ \sigma_{13} & \sigma_{11} \end{vmatrix} = \sigma_1\sigma_2 + \sigma_2\sigma_3 + \sigma_3\sigma_1}$$

$$\boxed{I_3 = \begin{vmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\ \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{vmatrix} = \sigma_1\sigma_2\sigma_3}$$

**不变量验证示例（$I_1$）**：
$$\sigma_{i'i'} = \sigma_{mn}l_{i'm}l_{i'n} = \sigma_{mn}\delta_{mn} = \sigma_{mm} = \sigma_{ii}$$

> **物理意义**：应力不变量不依赖坐标系选择，反映了应力状态的本质特征。$I_1$ 与静水压力相关，$I_2$ 与剪应力水平相关，$I_3$ 与体积变形相关。

> ⚠️ **重难点**：第二不变量 $I_2$ 的表达式有两种常见的写法——用行列式和用特征值乘积。考试中要会根据需要选择形式。另外，特征方程 $\sigma^3 - I_1\sigma^2 + I_2\sigma - I_3 = 0$ 中 $I_2$ 前的符号是 +（不是 -），这和某些线性代数教材中的写法不同，注意区分。

**重要结论**：由于应力张量是实对称矩阵，一定存在三个实主应力和一组正交主方向。

#### 3.1.4 应力张量的分解（球 + 偏）

应力张量可以分解为两部分：

**球应力张量（spherical stress tensor）**：
$$\begin{bmatrix} \sigma_m & 0 & 0 \\ 0 & \sigma_m & 0 \\ 0 & 0 & \sigma_m \end{bmatrix}, \quad \sigma_m = \frac{1}{3}\sigma_{ii} = \frac{1}{3}(\sigma_{11}+\sigma_{22}+\sigma_{33})$$

$\sigma_m$ 称为平均应力。球应力张量只引起体积变化，不引起形状变化。

**偏应力张量（deviatoric stress tensor）**：
$$S_{ij} = \sigma_{ij} - \sigma_m\delta_{ij}$$

偏应力张量代表应力中偏离各向同性（静水）压力的部分，只引起形状变化（剪切变形）。在塑性力学中尤为重要。

> 💡 **理解关键**：球 + 偏分解的思想——将复杂的应力状态拆成"各方向均匀压/拉"（球，像在水中受压）+ "偏离均匀的剪/扭"（偏，引起形状改变）。这在本构方程中也会出现：$\lambda e \delta_{ij}$ 项对应球应力部分，$2\mu \varepsilon_{ij}$ 包含偏量效应。塑性力学中的 von Mises 屈服准则用的就是偏应力第二不变量 $J_2$。

#### 3.1.5 平衡方程（Equilibrium Equation）

从整体平衡出发导出微元平衡方程。

**预备：Gauss 公式**（散度定理）

$$\iint_S (P\cos\alpha + Q\cos\beta + R\cos\gamma)dS = \iiint_V \left(\frac{\partial P}{\partial x} + \frac{\partial Q}{\partial y} + \frac{\partial R}{\partial z}\right)dV$$
（这里 $V$ 是物体所占的**体积域** (volume domain)，$S$ 是包围该体积的**边界面** (bounding surface)。）

用应力记号写为：
$$\int_S \sigma_{ij}l_j dS = \int_V \sigma_{ij,j} dV \quad (i=1,2,3)$$

> ⚠️ **重难点**：Gauss 公式是弹性力学推导中出场率最高的积分定理——把面积分转为体积分。记忆方法：$\int_S \sigma_{ij}l_j dS = \int_V \sigma_{ij,j} dV$，即"面力积分 → 应力散度积分"。注意 $\sigma_{ij,j}$ 是对 $j$ 求偏导并求和（缩并），$j$ 既是积分面法向指标又是散度中的导数指标。

**平衡推导**：

取物体内由闭合曲面 $S$ 包围的体积 $V$。面力为 $T_i$，体力为 $F_i$。整体的力平衡条件为：
$$\int_S T_i dS + \int_V F_i dV = 0 \quad (i=1,2,3)$$
（$F_i$ 是**体力分量** (body force components)，即单位体积上作用的体力（如重力）在 $i$ 方向上的分量。）

利用 Cauchy 公式 $T_i = \sigma_{ij}l_j$：
$$\int_S \sigma_{ij}l_j dS + \int_V F_i dV = 0$$

运用 Gauss 公式将面积分转换为体积分：
$$\int_V \sigma_{ij,j} dV + \int_V F_i dV = 0$$
$$\int_V (\sigma_{ij,j} + F_i) dV = 0$$

由于 $V$ 是任意的，被积函数必须恒为零，得**平衡方程**：

$$\boxed{\sigma_{ij,j} + F_i = 0} \quad (i = 1,2,3)$$

展开为三个方程：
$$\begin{aligned}
\frac{\partial \sigma_{11}}{\partial x_1} + \frac{\partial \sigma_{12}}{\partial x_2} + \frac{\partial \sigma_{13}}{\partial x_3} + F_1 = 0 \\
\frac{\partial \sigma_{21}}{\partial x_1} + \frac{\partial \sigma_{22}}{\partial x_2} + \frac{\partial \sigma_{23}}{\partial x_3} + F_2 = 0 \\
\frac{\partial \sigma_{31}}{\partial x_1} + \frac{\partial \sigma_{32}}{\partial x_2} + \frac{\partial \sigma_{33}}{\partial x_3} + F_3 = 0
\end{aligned}$$

> **核心思想**：平衡方程 + Gauss 公式 = 体积分到面积分的转换 $\to$ 微元平衡。这种"面积分 + Gauss 公式 $\to$ 被积函数为零"的方法在弹性力学中反复出现。

> 🔗 **跨章连接**：平衡方程 $\sigma_{ij,j} + F_i = 0$ 是后续第 3 章"变分法中推导 Euler 方程"和第 4 章"FEM 刚度矩阵推导"中反复出现的三大方程之一。三个基本方程在有限元框架中各有分工：平衡方程→力的平衡（$[\partial]^T\{\sigma\} + \{F\} = 0$）、几何方程→位移-应变关系（$\{\varepsilon\} = [\partial]\{u\}$）、本构方程→应力-应变关系（$\{\sigma\} = [C]\{\varepsilon\}$）。

---

> #### 检查你的理解：3.1
> 1. 写出 Cauchy 公式的物理含义：斜面应力向量的分量与应力张量及法向方向余弦之间的关系。如果已知应力张量和法向，如何求斜面正应力和剪应力？
> 2. 证明应力张量的第二不变量 $I_2$ 在坐标变换下不变（可以只证 2D 情况）。
> 3. 平衡方程是否适用于各向异性材料？连续性假设对平衡方程意味着什么？

---

### 3.2 变形分析（Deformation Analysis）

#### 3.2.1 构形与位移（Configuration and Displacement）

**构形（configuration）**：物体瞬时占据的空间区域，描绘了物体的位置和形状。

- $t = 0$ 时：**初始构形**（initial configuration）$V_0$，点 P 用坐标 $(x_1,x_2,x_3)$ 表达，位置向量 $\mathbf{r}$。
- $t$ 时刻：**当前构形**（current configuration）$V$，点 P 移动到 Q，坐标 $(\xi_1,\xi_2,\xi_3)$，位置向量 $\mathbf{R}$。

从初始构形到当前构形的数学变换：
$$\xi_i = \xi_i(x_1, x_2, x_3, t) \quad (i=1,2,3)$$

当 $t=0$ 时：$\xi_i(x_1, x_2, x_3, 0) = x_i$

**位移（displacement）**：
$$u_i = \xi_i - x_i \quad (i=1,2,3)$$

位移分为两类：
- **刚体位移（rigid body displacement）**：物体内各点的相对位置不变，仅由于空间中的刚体运动。
- **变形位移（deformation displacement）**：物体内各点的相对位置发生变化——这是弹性力学关心的位移类型，因为它与应力直接相关。

> 💡 **理解关键**：区分 Lagrangian 描述（跟踪物质点 P 从初始位置 $x_i$ 到当前位置 $\xi_i$）和 Eulerian 描述（看空间中固定位置处流过什么物质）。弹性力学通常用 Lagrangian 描述，因为固体的变形是和"从哪儿来的"（初始构形）绑定的。注意本章中 $u_i$ 的定义用初始坐标 $(x_1,x_2,x_3)$ 表示，而不是当前坐标。

**变形分析的核心**：通过描述物体内任意两点的距离变化来刻画变形。

#### 3.2.2 Green 应变与小变形应变张量

考虑相邻两点 P $(x_1,x_2,x_3)$ 和 P' $(x_1+dx_1, x_2+dx_2, x_3+dx_3)$。

**变形前**线元 PP' 长度的平方：
$$ds_0^2 = dx_1^2 + dx_2^2 + dx_3^2 = dx_i dx_i$$

**变形后**对应线元 QQ' 长度的平方：
$$ds^2 = d\xi_1^2 + d\xi_2^2 + d\xi_3^2 = d\xi_i d\xi_i$$

由 $\xi_i = \xi_i(x_1,x_2,x_3)$：$d\xi_i = \frac{\partial \xi_i}{\partial x_j}dx_j = \xi_{i,j}dx_j$

线元长度平方之差（变形前后的差变化）：
$$ds^2 - ds_0^2 = d\xi_k d\xi_k - dx_i dx_i = \xi_{k,i}dx_i \xi_{k,j}dx_j - \delta_{ij}dx_i dx_j$$
$$= (\xi_{k,i}\xi_{k,j} - \delta_{ij})dx_i dx_j$$

引入 **Green 应变张量** $E_{ij}$：
$$\boxed{E_{ij} = \frac{1}{2}(\xi_{k,i}\xi_{k,j} - \delta_{ij})}$$

则：
$$ds^2 - ds_0^2 = 2E_{ij}dx_i dx_j$$

**将 Green 应变用位移表达**：

利用 $u_i = \xi_i - x_i$，即 $\xi_i = u_i + x_i$：
$$\frac{\partial \xi_i}{\partial x_j} = \frac{\partial (u_i + x_i)}{\partial x_j} = u_{i,j} + \delta_{ij}$$

代入 Green 应变表达式：
$$E_{ij} = \frac{1}{2}[(u_{k,i} + \delta_{ki})(u_{k,j} + \delta_{kj}) - \delta_{ij}]$$
$$= \frac{1}{2}(u_{i,j} + u_{j,i} + u_{k,i}u_{k,j})$$

其中 $u_{k,i}u_{k,j}$ 是位移梯度的二阶项（非线性项）。

> ⚠️ **重难点**：Green 应变展开中有一个容易跳过的细节——$(u_{k,i} + \delta_{ki})(u_{k,j} + \delta_{kj})$ 乘开得四项：$u_{k,i}u_{k,j} + u_{k,i}\delta_{kj} + \delta_{ki}u_{k,j} + \delta_{ki}\delta_{kj}$。第三项 $\delta_{ki}u_{k,j} = u_{i,j}$（$\delta$ 把 $k$ 换成 $i$），第二项 $u_{k,i}\delta_{kj} = u_{j,i}$，第四项 $\delta_{ki}\delta_{kj} = \delta_{ij}$。所以 $E_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i} + u_{k,i}u_{k,j})$。一定要亲手推导一遍。

**小变形假设下的简化**：在小变形条件下，位移梯度很小（$|u_{i,j}| \ll 1$），二阶项 $u_{k,i}u_{k,j}$ 远小于一阶项，可以略去：

$$\boxed{\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})}$$

这就是**小变形应变张量**（又称 Cauchy 应变张量或工程应变张量），是 Green 应变在小变形条件下的线性近似。


> **重要理解**：
> - **几何方程**：$\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$——将位移梯度与应变建立联系
> - Green 应变包含完整的非线性项，适用于大变形分析
> - 小变形应变只有线性项，是弹性力学中使用最广泛的几何关系
> - 可以证明 $\varepsilon_{ij}$ 是二阶张量（满足张量坐标变换律）

#### 3.2.3 应变张量的物理意义

**正应变的物理意义**：

考虑平行于 $x_1$ 轴方向的线元：$dx_1 = ds_0, dx_2 = dx_3 = 0$。变形后长度为 $ds$，相对伸长（relative elongation）：
$$L_1 = \frac{ds - ds_0}{ds_0}, \quad ds = (1+L_1)ds_0$$

由 $ds^2 - ds_0^2 = 2E_{ij}dx_i dx_j$，小变形下 $E_{ij} \approx \varepsilon_{ij}$：
$$ds^2 - ds_0^2 \approx 2\varepsilon_{11}dx_1^2 = 2\varepsilon_{11}ds_0^2$$

又 $ds^2 - ds_0^2 = [(1+L_1)^2 - 1]ds_0^2$，对比得：
$$(1+L_1)^2 - 1 = 2\varepsilon_{11} \;\Rightarrow\; L_1 = \sqrt{1+2\varepsilon_{11}} - 1 \approx \varepsilon_{11}$$

> **结论**：$\varepsilon_{11}$（即 $\varepsilon_{xx}$）在小变形下等于工程正应变（沿 $x_1$ 方向的相对伸长）。

**剪应变的物理意义**：

考虑两个互相垂直的线元：一个沿 $x_1$ 方向（$dx_1 = ds_0, dx_2 = dx_3 = 0$），另一个沿 $x_2$ 方向（$dx_2' = ds_0', dx_1' = dx_3' = 0$）。变形后两线元的内积：
$$d\mathbf{s} \cdot d\mathbf{s}' = ds\,ds'\cos\theta = d\xi_k d\xi_k' = \xi_{k,1}\xi_{k,2}\,ds_0 ds_0'$$

由 Green 应变 $E_{12} = \frac{1}{2}(\xi_{k,1}\xi_{k,2} - \delta_{12}) = \frac{1}{2}\xi_{k,1}\xi_{k,2}$，而在小变形下：
$$ds = \sqrt{1+2\varepsilon_{11}}ds_0 \approx ds_0, \quad ds' = \sqrt{1+2\varepsilon_{22}}ds_0' \approx ds_0'$$

$$\cos\theta = \frac{2E_{12}}{\sqrt{1+2E_{11}}\sqrt{1+2E_{22}}} \approx 2\varepsilon_{12}$$

设 $\alpha_{12} = \frac{\pi}{2} - \theta$ 表示直角的角度变化量（剪应变）。在小变形下：
$$\sin\alpha_{12} \approx \alpha_{12} \approx 2\varepsilon_{12} = u_{1,2} + u_{2,1}$$

> **结论**：$2\varepsilon_{12}$ 在小变形下等于工程剪应变 $\gamma_{xy}$（$x_1$ 与 $x_2$ 方向间直角的角度变化量）。

**总结**：
$$\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33} \quad\text{对应}\quad \text{正应变 } \varepsilon_x, \varepsilon_y, \varepsilon_z$$
$$2\varepsilon_{12} = \gamma_{xy},\; 2\varepsilon_{23} = \gamma_{yz},\; 2\varepsilon_{13} = \gamma_{xz} \quad\text{对应}\quad \text{工程剪应变}$$

> 注意：张量剪应变 $\varepsilon_{ij}(i \neq j)$ 是工程剪应变 $\gamma_{ij}$ 的一半！这是很多教材中容易混淆的地方。

> ⚠️ **重难点**：张量剪应变 vs 工程剪应变的区别是 Voigt 记法中"剪应变是否乘 2"的根源。张量 $\varepsilon_{12} = \frac{1}{2}(u_{1,2}+u_{2,1})$，工程 $\gamma_{12} = u_{1,2} + u_{2,1} = 2\varepsilon_{12}$。两者的物理意义不同：张量剪应变是"对称梯度的一半"，工程剪应变是"直角角度的变化量"。在 MATLAB/FEA 编程中要特别注意用的是哪个——商业 FEA 软件内部通常用张量应变的 Voigt 记法（剪应变乘 2 的那种）。

#### 3.2.4 刚体变形与旋转张量

物体在外力作用下产生的位移由刚体运动和变形共同引起。下面在小变形假设下分离刚体转动和伸长变形。

相邻点 B 相对点 A 的位移差（Taylor 展开，略去高阶项）：
$$du_i = \frac{\partial u_i}{\partial x_j}dx_j = u_{i,j}dx_j$$

**位移梯度张量** $u_{i,j}$ 通常不对称（即 $\partial u_i/\partial x_j \neq \partial u_j/\partial x_i$），因此可以分解为对称部分和反对称部分：

$$\boxed{u_{i,j} = \underbrace{\frac{1}{2}(u_{i,j} + u_{j,i})}_{\varepsilon_{ij}\text{（应变张量）}} + \underbrace{\frac{1}{2}(u_{i,j} - u_{j,i})}_{\omega_{ij}\text{（旋转张量）}}}$$

- **对称部分** $\varepsilon_{ij}$ = 小变形应变张量（表征纯变形）
- **反对称部分** $\omega_{ij}$ = 旋转张量（rotation tensor），满足 $\omega_{ij} = -\omega_{ji}$

> 💡 **理解关键**：位移梯度 $u_{i,j}$ 的物理含义——相邻点的位移变化。拆成对称 + 反对称后：对称部分 $\varepsilon_{ij}$ 描述"两点间距的变化"（→应力），反对称部分 $\omega_{ij}$ 描述"两点的相对旋转"（→不产生应力）。这就像一块橡皮——你可以捏它（对称→应变），也可以整体旋转它（反对称→刚体转动），只有前者产生内部力。

**旋转张量的物理意义**：

考虑 $x_1x_2$ 平面内，通过点 P 作两条分别平行于 $x_1$ 和 $x_2$ 轴的线元 PA 和 PB。在小变形假设下：
$$\alpha_1 = \frac{\partial u_2}{\partial x_1} \quad (\text{逆时针方向}), \qquad \alpha_2 = \frac{\partial u_1}{\partial x_2} \quad (\text{顺时针方向})$$

其中 $\alpha_1$ 是 PA 的转角，$\alpha_2$ 是 PB 的转角。

这个旋转运动可分解为：
- **刚体转动**：整体逆时针转过 $(\alpha_1 - \alpha_2)/2$（对应旋转张量 $\omega_{12}$）
- **剪切变形**：PA 再逆时针转过 $(\alpha_1 + \alpha_2)/2$，PB 再顺时针转过 $(\alpha_1 + \alpha_2)/2$（对应应变张量 $\varepsilon_{12}$）

剪切变形（工程剪应变）：
$$\gamma_{12} = \alpha_1 + \alpha_2 = \frac{\partial u_1}{\partial x_2} + \frac{\partial u_2}{\partial x_1}$$

旋转角（绕 $x_3$ 轴）：
$$\Omega_3 = \frac{\alpha_1 - \alpha_2}{2} = \frac{1}{2}\left(\frac{\partial u_1}{\partial x_2} - \frac{\partial u_2}{\partial x_1}\right) = \omega_{12}$$

**旋转张量和旋转向量的关系**：用置换符号表达
$$\boxed{\Omega_k = -\frac{1}{2}e_{ijk}\omega_{ij}} \quad\text{或}\quad \boxed{\omega_{pq} = -e_{pqk}\Omega_k}$$

> **物理理解**：位移梯度 = 应变（对称） + 旋转（反对称）。在小变形下，刚体转动不影响应力——只有应变张量才与应力直接相关。


---

> #### 检查你的理解：3.2
> 1. Green 应变中的 $u_{k,i}u_{k,j}$ 项何时可以忽略？如果位移梯度达到 $10^{-2}$ 量级，这项大约是一阶项的百分之几？
> 2. 如果在 $x_1x_2$ 平面内只有 $\partial u_1/\partial x_2 = 0.01$，$\partial u_2/\partial x_1 = 0.01$，这代表纯刚体转动还是纯剪切变形？
> 3. 证明 $\varepsilon_{ij}$ 是二阶张量：从 $\varepsilon_{ij} = \frac{1}{2}(u_{i,j}+u_{j,i})$ 出发，利用位移是向量（$u_{i'} = l_{i'j}u_j$）和导数链式法则推导 $\varepsilon_{i'j'} = l_{i'm}l_{j'n}\varepsilon_{mn}$。

---

### 3.3 本构关系（Constitutive Relations）

平衡方程和几何方程分别来源于力学原理和几何学，它们反映了变形过程中必须遵循的运动规律和连续性条件，但**不包含任何材料的物理特征**——因此对任何连续体都适用。

弹性力学的目的是研究物体在外部作用（载荷、温度等）下的力学响应（应力、应变等）。仅有平衡方程和几何方程是不够的，还需要建立**应力与应变之间**的关系——这就是**本构方程（constitutive equation）**。

> 💡 **理解关键**：为什么需要三组方程？平衡方程管"力必须平衡"（力学公理），几何方程管"位移和应变的关系"（几何必然性），但两者之间没有直接联系——同样的力作用在铁和橡皮上，变形完全不同。这就是本构方程的作用：架起应力到应变的桥，引入材料的"性格"。铁→$E$ 大→应变压应力小，橡皮→$E$ 小→应变压应力大。

#### 3.3.1 热力学框架

变形过程本质上是热力学过程。根据热力学基本定律建立理想弹性体的应变能函数框架。

**热力学第一定律**：封闭系统在无限小时间间隔内，内能的变化等于外界输入的热量加上外力做的功。

物体总能量 = 动能 $K$ + 内能 $U$。仅考虑静力学问题（$dK = 0$）：
$$dA + dQ = dU \quad\text{即}\quad dU - dQ = dA$$

外力功：体力 $F_i$ + 面力 $p_i$，对应位移增量 $du_i$：
$$dA = \int_V F_i du_i dV + \int_S p_i du_i dS$$

代入应力边界条件 $p_i = \sigma_{ij}l_j$，用 Gauss 公式转化面积分：
$$dA = \int_V F_i du_i dV + \int_S \sigma_{ij}l_j du_i dS = \int_V F_i du_i dV + \int_V (\sigma_{ij}du_i)_{,j}dV$$

展开括号并利用平衡方程 $\sigma_{ij,j} + F_i = 0$：
$$
\begin{aligned}
dA &= \int_V [F_i du_i + \sigma_{ij,j}du_i + \sigma_{ij}du_{i,j}] dV \\
&= \int_V \underbrace{(F_i + \sigma_{ij,j})}_{=0} du_i dV + \int_V \sigma_{ij} \cdot \frac{1}{2}(du_{i,j} + du_{j,i}) dV \\
&= \int_V \sigma_{ij} d\varepsilon_{ij} dV
\end{aligned}
$$

因此：
$$dU - dQ = \int_V \sigma_{ij} d\varepsilon_{ij} dV$$

> ⚠️ **重难点**：热力学推导的核心一步——从 $\sigma_{ij}du_{i,j}$ 到 $\sigma_{ij}d\varepsilon_{ij}$。$du_{i,j}$ 不是对称的（不可直接等于 $d\varepsilon_{ij}$），但因为 $\sigma_{ij}$ 是对称的，所以 $\sigma_{ij}du_{i,j} = \sigma_{ij}\frac{1}{2}(du_{i,j}+du_{j,i}) = \sigma_{ij}d\varepsilon_{ij}$。这里用到了一个恒等式：对称张量与任意张量的缩并等于对称张量与该张量对称部分的缩并。考试中可能考察这个推导逻辑。

**热力学第二定律**：对于理想弹性材料的可逆变形过程，在热平衡过程中始终存在统一的温度 $T$：
$$T dS = dQ$$

其中 $S$ 是熵（entropy），是与材料状态相关的状态量。

#### 3.3.2 应变能函数

引入 **Helmholtz 自由能** $F = U - ST$：
$$dF = dU - TdS - SdT$$

结合 $dU - TdS = \int_V \sigma_{ij}d\varepsilon_{ij} dV$，得：
$$dF = \int_V \sigma_{ij} d\varepsilon_{ij} dV - SdT$$

**绝热过程**（如冲击加载，$dQ = 0$）：
$$dU = \int_V \sigma_{ij} d\varepsilon_{ij} dV$$

**等温过程**（缓慢加载，$dT = 0$）：
$$dF = \int_V \sigma_{ij} d\varepsilon_{ij} dV$$

两种过程中，内能或自由能的变化仅由应变状态的改变引起，因此可以定义为应变能函数。引入**应变能密度** $W$（单位体积的应变能）：
$$d\int_V W dV = \int_V dW dV = \int_V \sigma_{ij} d\varepsilon_{ij} dV$$

由于积分区域任意，对任意点有：
$$dW = \sigma_{ij} d\varepsilon_{ij}$$

$W$ 是应变状态的单值函数（与变形路径无关），因此：

$$\boxed{\sigma_{ij} = \frac{\partial W}{\partial \varepsilon_{ij}}}$$

> 这就是 **Green 公式**——建立了应力分量等于应变能密度对对应应变分量偏导数的关系。这给出了建立应力-应变关系的一般规则。

> 🔗 **跨章连接**：Green 公式 $\sigma_{ij} = \partial W / \partial \varepsilon_{ij}$ 是建立本构关系的核心公式。在第 3 章变分法中，它也是推导 Euler 方程时的重要出发点——应变能函数 $W$ 出现在势能泛函中，最终通过变分导出平衡方程。

**关键认识**：在绝热和等温两种过程中，应变能密度函数 $W$ 的具体数值可能略有不同（因为弹性系数不同），但实验证明差异非常小，在小变形条件下可认为两种过程中的本构关系相同。

#### 3.3.3 线弹性本构推导——从 81 到 2 的简化

将应变能密度 $W$ 在 $\varepsilon_{ij} = 0$ 处 Taylor 展开。在小变形假设下（$\varepsilon_{ij} \ll 1$），略去三阶及以上项：
$$W(\varepsilon_{ij}) = a_{kl}\varepsilon_{kl} + b_{klrs}\varepsilon_{kl}\varepsilon_{rs}$$

由 Green 公式：
$$\sigma_{ij} = \frac{\partial W}{\partial \varepsilon_{ij}} = a_{kl}\frac{\partial\varepsilon_{kl}}{\partial\varepsilon_{ij}} + b_{klrs}\left(\frac{\partial\varepsilon_{kl}}{\partial\varepsilon_{ij}}\varepsilon_{rs} + \frac{\partial\varepsilon_{rs}}{\partial\varepsilon_{ij}}\varepsilon_{kl}\right)$$

利用 $\frac{\partial\varepsilon_{kl}}{\partial\varepsilon_{ij}} = \delta_{ki}\delta_{lj}$：
$$\sigma_{ij} = a_{ij} + (b_{ijkl} + b_{klij})\varepsilon_{kl}$$

当 $\varepsilon_{ij} = 0$ 时 $\sigma_{ij} = 0$，故 $a_{ij} = 0$。定义弹性系数张量：
$$c_{ijkl} = b_{klij} + b_{ijkl}$$

得**广义 Hooke 定律**：

$$\boxed{\sigma_{ij} = c_{ijkl}\varepsilon_{kl}}$$

> - 81→36：$\sigma_{ij}$ 的对称性（力矩平衡→剪应力互等）
> - 36→（保持36）：$\varepsilon_{kl}$ 的对称性（应变定义决定）
> - 36→21：应变能存在 → 混合偏导与次序无关 → $c_{ijkl} = c_{klij}$（"主对称性"）
> - 21→2：各向同性假设 → 材料性质与方向无关 → 只剩 2 个独立常数

**简化步骤**：

| 阶段 | 独立常数 | 原因 |
|------|---------|------|
| 初始 | **81** | $c_{ijkl}$ 是四阶张量，$3^4 = 81$ 个分量 |
| $\sigma_{ij}$ 对称（$i \leftrightarrow j$） | **36** | $\sigma_{ij} = \sigma_{ji} \Rightarrow c_{ijkl} = c_{jikl}$，$3 \times 3 \times 3 \times 3 \to 6 \times 9$ |
| $\varepsilon_{kl}$ 对称（$k \leftrightarrow l$） | **36** | $\varepsilon_{kl} = \varepsilon_{lk} \Rightarrow c_{ijkl} = c_{ijlk}$ |
| 由应变能存在：$c_{ijkl} = c_{klij}$ | **21** | 应变能是势函数，$\frac{\partial^2 W}{\partial\varepsilon_{ij}\partial\varepsilon_{kl}}$ 与求导次序无关 |
| 各向同性（材料性质与方向无关） | **2** | 应变能退化为仅依赖应变不变量的函数 |

> 💡 **理解关键**：36→21 这一步不那么直观。$c_{ijkl} = c_{klij}$ 源于 $\sigma_{ij} = \partial W/\partial\varepsilon_{ij}$ 和 $\partial^2 W/(\partial\varepsilon_{ij}\partial\varepsilon_{kl}) = \partial^2 W/(\partial\varepsilon_{kl}\partial\varepsilon_{ij})$（偏导次序无关）。这称为"主对称性"，与应力/应变各自的对称性（"次对称性"）不同。在 Voigt 记法中，主对称性意味着弹性矩阵 $[C]$ 是对称的。

> **各向同性线性弹性材料的最终结论**：只有 **2 个** 独立的弹性常数。

#### 3.3.4 Lame 常数与 Hooke 定律

对于各向同性线弹性材料，应变能密度是标量函数——应变不变量的函数。应变张量只有三个独立不变量 $I_1, I_2, I_3$（分别为 $\varepsilon_{ij}$ 的一次、二次、三次齐次式），更高次的不变量可由低次组合表示。

应变能密度取应变分量的二次齐次式：
$$W = AI_1^2 + BI_2$$

用 Lame 常数 $\lambda$ 和 $\mu$ 代替 A、B：
（$\lambda$ 和 $\mu$ 是描述各向同性线弹性材料力学行为的两个独立**Lame 弹性常数**；$\mu$ 即剪切模量 (shear modulus)，反映材料抵抗剪切变形的能力；$\lambda$ 则与材料的体积压缩特性有关。）

$$W = \frac{\lambda}{2}I_1^2 + \mu(\varepsilon_{11}^2 + \varepsilon_{22}^2 + \varepsilon_{33}^2) + \frac{\mu}{2}(\gamma_{12}^2 + \gamma_{23}^2 + \gamma_{31}^2)$$

或等价地：
$$W = \left(\frac{\lambda}{2}+\mu\right)I_1^2 - 2\mu I_2$$

由 Green 公式 $\sigma_{ij} = \partial W/\partial\varepsilon_{ij}$ 逐项求导，得各向同性线弹性本构方程：

$$\boxed{\sigma_{ij} = \lambda e\delta_{ij} + 2\mu\varepsilon_{ij}}$$

其中 $e = \varepsilon_{kk} = \varepsilon_{11} + \varepsilon_{22} + \varepsilon_{33}$ 是体积应变（volumetric strain）。


展开为工程分量形式：
$$\begin{cases}
\sigma_{11} = \lambda e + 2\mu\varepsilon_{11}, & \tau_{12} = \mu\gamma_{12} \\
\sigma_{22} = \lambda e + 2\mu\varepsilon_{22}, & \tau_{23} = \mu\gamma_{23} \\
\sigma_{33} = \lambda e + 2\mu\varepsilon_{33}, & \tau_{13} = \mu\gamma_{13}
\end{cases}$$

> $\lambda, \mu$ 称为 **Lame 弹性常数**。$\mu$ 就是剪切模量（shear modulus），也常记为 $G$。

#### 3.3.5 Lame 常数与工程常数的关系

工程中更常用的是杨氏模量 $E$ 和泊松比 $\nu$。
（$E$ 是**杨氏模量** (Young's modulus)，衡量材料抵抗单向拉伸或压缩变形的刚度；$\nu$ 是**泊松比** (Poisson's ratio)，定义为横向收缩应变与纵向拉伸应变之比的绝对值。）既然各向同性材料只有两个独立弹性常数，$\lambda,\mu$ 和 $E,\nu$ 之间必然存在确定关系。

**简单拉伸试验**：
$$\sigma_{11} = \sigma, \quad \sigma_{22} = \sigma_{33} = \tau_{12} = \tau_{23} = \tau_{13} = 0$$

工程关系：
$$\varepsilon_{11} = \frac{\sigma}{E}, \quad \varepsilon_{22} = \varepsilon_{33} = -\nu\frac{\sigma}{E}, \quad \gamma_{12} = \gamma_{23} = \gamma_{13} = 0$$

体积应变：$e = \varepsilon_{11} + \varepsilon_{22} + \varepsilon_{33} = \frac{\sigma}{E}(1-2\nu)$

代入 Lame 形式的本构方程：
$$\sigma_{11} = \lambda e + 2\mu\varepsilon_{11} = \frac{\sigma}{E}[\lambda(1-2\nu) + 2\mu] = \sigma$$
$$\sigma_{22} = \lambda e + 2\mu\varepsilon_{22} = \frac{\sigma}{E}[\lambda(1-2\nu) - 2\mu\nu] = 0$$

解此二元方程组得：

$$\boxed{\mu = \frac{E}{2(1+\nu)}, \qquad \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}}$$

**简单剪切试验**（薄壁圆筒扭转）：

$$\tau_{12} = \tau, \quad \sigma_{11} = \sigma_{22} = \sigma_{33} = \tau_{23} = \tau_{13} = 0$$

工程关系：
$$\gamma_{12} = \frac{\tau}{G}, \quad \varepsilon_{11} = \varepsilon_{22} = \varepsilon_{33} = \gamma_{23} = \gamma_{13} = 0, \quad e = 0$$

代入 Lame 本构：$\tau_{12} = \mu\gamma_{12} \;\Rightarrow\; \tau = \mu \cdot \frac{\tau}{G} \;\Rightarrow\; G = \mu$

$$\boxed{G = \mu = \frac{E}{2(1+\nu)}}$$

> **总结**：$\lambda$ 和 $\mu$（即 $G$）与工程常数 $E,\nu$ 之间的关系如上。在编程中常直接用 $E,\nu$ 计算 $\lambda,\mu$，然后在代码中使用 Lame 常数的形式——因为张量表达式 $\sigma_{ij} = \lambda e\delta_{ij} + 2\mu\varepsilon_{ij}$ 是最紧凑的形式。

> ❌ **易错点**：$\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}$ 的分母中有 $(1-2\nu)$。对于不可压缩材料（如橡胶，$\nu = 0.5$），$1-2\nu = 0$ 导致 $\lambda \to \infty$。很多学生在计算 Lame 常数时忘记检查 $\nu$ 的取值范围——如果 $\nu$ 接近 0.5，数值上会出现除零问题。工程中常用 $E$ 和 $\nu$，FEM 编程中通常预计算 $\lambda$ 和 $\mu$（或 $G$）。

**常用弹性常数关系速查**：

| Lame 常数 | 用 $E,\nu$ 表示 | 物理意义 |
|----------|----------------|---------|
| $\mu = G$ | $\frac{E}{2(1+\nu)}$ | 剪切模量 |
| $\lambda$ | $\frac{E\nu}{(1+\nu)(1-2\nu)}$ | Lame 第一常数 |
| $K$（体积模量） | $\frac{E}{3(1-2\nu)} = \lambda + \frac{2}{3}\mu$ | 静水压力与体积应变比 |

---

> #### 检查你的理解：3.3
> 1. 为什么一般各向异性线弹性材料最多有 21 个独立弹性常数？推导过程中哪一步将 36 降到 21？
> 2. 已知 $E = 200$ GPa, $\nu = 0.3$，求 Lame 常数 $\lambda$ 和 $\mu$。
> 3. 在 $\sigma_{ij} = \lambda e\delta_{ij} + 2\mu\varepsilon_{ij}$ 中取 $i=j$ 并对 $i$ 求和，推导体积应力与体积应变的关系 $\sigma_{kk} = (3\lambda + 2\mu)e$。

---

### 3.4 弹性力学微分方程汇总

将平衡方程、几何方程、本构方程联立，再加上边界条件，构成完整的弹性力学边值问题。

#### 3.4.1 15 个方程汇总

| 方程类型 | 张量形式 | 分量数 | 未知量 |
|---------|---------|-------|--------|
| 平衡方程 | $\sigma_{ij,j} + F_i = 0$ | 3 | $\sigma_{ij}$（6 个独立） |
| 几何方程 | $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$ | 6 | $u_i$（3 个）, $\varepsilon_{ij}$（6 个独立） |
| 本构方程 | $\sigma_{ij} = \lambda \delta_{ij} e + 2\mu\varepsilon_{ij}$ | 6 | — |

**未知量总计**：15 个（6 个应力 + 6 个应变 + 3 个位移）
**方程总计**：15 个（3 + 6 + 6 = 15）

> 15 个方程，15 个未知量——方程组封闭。加上适当的边界条件，理论上可解。


#### 3.4.2 基本假设与各方程的关联

分析各基本假设对各组方程的影响：

| 假设 | 平衡方程 | 几何方程 | 本构方程 |
|------|---------|---------|---------|
| 连续性 | $\checkmark$ | $\checkmark$ | $\checkmark$ |
| 均匀性 | | | $\checkmark$ |
| 各向同性 | | | $\checkmark$ |
| 线弹性 | | | $\checkmark$ |
| 小变形 | $\checkmark$ | $\checkmark$ | $\checkmark$ |

**关键理解**：
- **平衡方程和几何方程**与材料属性无关——对非均匀、各向异性、不满足 Hooke 定律的材料同样适用。
- **小变形假设**对平衡方程和几何方程至关重要：只有小变形才可用变形前位置建立平衡方程；只有小变形下几何方程才只含位移导数的线性项。
- **本构方程**需要所有五个假设，其中线弹性假设结合小变形假设才能得到线性应力-应变关系。


#### 3.4.3 边界条件

基本方程反映弹性体**内部**的应力、应变和位移规律。要解决具体的边值问题，必须给出**边界条件**。

**三种边界条件**：

1. **应力边界条件**（在 $S_\sigma$ 上）：
   $$\boxed{\sigma_{ij}l_j = \overline{p}_i}$$
   式中 $\overline{p}_i$ 为边界上给定的面力分量。

2. **位移边界条件**（在 $S_u$ 上）：
   $$\boxed{u_i = \overline{u}_i}$$
   式中 $\overline{u}_i$ 为边界上给定的位移分布。

3. **混合边界条件**：部分边界给应力条件 $S_\sigma$，其余部分给位移条件 $S_u$。注意：
   - $S_\sigma \cup S_u = S$（边界全覆盖）
   - $S_\sigma \cap S_u = \emptyset$（不重叠）
   - 边界上每点需要且仅需给定三个正交方向的条件
   - **同一点同一方向不能同时给应力条件和位移条件**

> 💡 **理解关键**：应力边界条件 $\sigma_{ij}l_j = \overline{p}_i$ 本质上就是 Cauchy 公式——它说的是边界处的斜面应力向量必须等于外加面力。位移边界条件 $u_i = \overline{u}_i$ 等价于材料力学中的"固定端"或"简支约束"。

#### 3.4.4 弹性力学边值问题的数学表述

完整而严格的弹性力学问题数学表述（张量形式）：

$$\begin{cases}
\text{在 V 内:} & \sigma_{ij,j} + F_i = 0 & \text{（平衡方程）} \\
\text{在 V 内:} & \varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i}) & \text{（几何方程）} \\
\text{在 V 内:} & \sigma_{ij} = D_{ijkl}\varepsilon_{kl} & \text{（本构方程）} \\
\text{在 } S_\sigma \text{ 上:} & \sigma_{ij}l_j = \overline{p}_i & \text{（应力边界条件）} \\
\text{在 } S_u \text{ 上:} & u_i = \overline{u}_i & \text{（位移边界条件）}
\end{cases}$$

其中 $D_{ijkl}$ 为一般线弹性本构张量（对各向同性材料退化为 $\lambda \delta_{ij}\delta_{kl} + \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})$）。

> 以上是弹性力学问题的微分方程表述。除此之外，还有等价的变分表述——通过建立泛函并利用泛函的极值条件求解。两种方法本质相同，泛函的 Euler 方程正好对应于弹性力学的基本方程和边界条件。

> 🔗 **跨章连接**：这组边值问题表述是第 3 章变分法的出发点。第 3 章的核心思想就是——把"求解微分方程+边界条件"转化为"寻找使某个泛函取极值的函数"，从而为 Ritz 法和有限元法建立理论基础。简而言之：微分表述（本章）→ 变分表述（第 3 章）→ FEM 代数化（第 4-6 章）。

---

> #### 检查你的理解：3.4
> 1. 为什么说"小变形假设对平衡方程和几何方程至关重要"？如果在大变形下，平衡方程应在哪种构形上建立？
> 2. 某弹性体表面某点外法向为 $(1,0,0)$，已知该点应力张量，写出该点的应力边界条件。如果该点位移被固定呢？
> 3. 一个纯应力边值问题和一个纯位移边值问题分别需要满足什么条件才有唯一解？

---

### 3.5 补充内容（Supplementary）

#### 3.5.1 Voigt 规则

张量记法在理论推导中极具优势，但在工程应用和有限元编程中常不如矩阵和列向量方便。**Voigt 规则**（Voigt notation）将高阶张量降阶为低阶矩阵/向量形式。

> 💡 **理解关键**：Voigt 规则是理论推导（张量）和实际编程（矩阵/向量）之间的桥梁。在 FEM 代码中，你永远不会看到一个 3×3×3×3 的四阶张量——全都被 Voigt 压缩成 6×6 矩阵了。

**动态量 Voigt 规则**（应力型，如 $\sigma_{ij}$）：

2D 情形（$\sigma_{ij}$ 只有 $i,j = 1,2$）：
$$\begin{bmatrix} \sigma_{11} & \sigma_{12} \\ \sigma_{21} & \sigma_{22} \end{bmatrix} \to \begin{Bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{Bmatrix} = \begin{Bmatrix} \sigma_{xx} \\ \sigma_{yy} \\ \tau_{xy} \end{Bmatrix}$$

索引映射：

| $ij$ | $a$（Voigt 索引） |
|------|------------------|
| 11 | 1 |
| 22 | 2 |
| 12 | 3 |

3D 情形：
$$\begin{bmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\ \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{bmatrix} \to \begin{Bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{Bmatrix}$$

索引映射（记忆方法：沿主对角线下行，在最后一列上转，回到第一行）：

| $ij$ | 11 | 22 | 33 | 23 | 13 | 12 |
|------|----|----|----|----|----|
| $a$ | 1 | 2 | 3 | 4 | 5 | 6 |

> ❌ **易错点**：Voigt 索引的 3D 排列不规范！11→1, 22→2, 33→3，但之后不是"12→4, 13→5, 23→6"，而是"23→4, 13→5, 12→6"。记忆口诀："对角线下行 123，最后一列回头 4，第一行最后 5，第二行最后 6"——或者说"11-22-33-23-13-12"。在 FEA 编程中索引映射错误是常见的 bug 来源。

**运动量 Voigt 规则**（应变型，如 $\varepsilon_{ij}$）：

与应力型规则相同，但**剪切分量乘以 2**：

2D：
$$\begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} \\ \varepsilon_{21} & \varepsilon_{22} \end{bmatrix} \to \begin{Bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{Bmatrix} = \begin{Bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \gamma_{xy} \end{Bmatrix}$$

3D：
$$\begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} & \varepsilon_{13} \\ \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23} \\ \varepsilon_{31} & \varepsilon_{32} & \varepsilon_{33} \end{bmatrix} \to \begin{Bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\ 2\varepsilon_{23} \\ 2\varepsilon_{13} \\ 2\varepsilon_{12} \end{Bmatrix} = \begin{Bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \\ \gamma_{yz} \\ \gamma_{zx} \\ \gamma_{xy} \end{Bmatrix}$$

> **为什么应变要乘 2？** 为了保证能量表达式的一致性：
> $$d\varepsilon_{ij}\sigma_{ij} = \{d\boldsymbol{\varepsilon}\}^T\{\boldsymbol{\sigma}\}$$
> 在 Voigt 记法下，张量的双重收缩等于 Voigt 向量的内积。

> ⚠️ **重难点**：动态量（应力）不乘 2，运动量（应变）乘 2——这是 Voigt 规则最大的坑。为什么不对称？因为能量 $dW = \sigma_{ij}d\varepsilon_{ij}$ 中应变的非对角项（$i \neq j$）在缩并时出现两次（$\sigma_{12}\varepsilon_{12} + \sigma_{21}\varepsilon_{21} = 2\sigma_{12}\varepsilon_{12}$），应变乘 2 后恰好抵消这个 2 倍，保证 $\{\varepsilon\}^T\{\sigma\} = \sigma_{ij}\varepsilon_{ij}$。如果搞反（应力乘 2 而应变不乘），总能量会出错。

**四阶张量的 Voigt 规则**（本构矩阵）：

线弹性定律的指标形式：
$$\sigma_{ij} = C_{ijkl}\varepsilon_{kl}$$

Voigt 矩阵形式：
$$\{\boldsymbol{\sigma}\} = [\mathbf{C}]\{\boldsymbol{\varepsilon}\}$$

其中 $a \leftrightarrow ij$，$b \leftrightarrow kl$。例如平面应变问题的 $[\mathbf{C}]$ 矩阵：
$$[\mathbf{C}] = \begin{bmatrix} C_{11} & C_{12} & C_{13} \\ C_{21} & C_{22} & C_{23} \\ C_{31} & C_{32} & C_{33} \end{bmatrix} = \begin{bmatrix} C_{1111} & C_{1122} & C_{1112} \\ C_{2211} & C_{2222} & C_{2212} \\ C_{1211} & C_{1222} & C_{1212} \end{bmatrix}$$

**注意要点**：四阶张量降阶时要考虑 Voigt 规则中剪应变的 2 倍因子。验证：
$$\sigma_{12} = C_{1211}\varepsilon_{11} + C_{1212}\varepsilon_{12} + C_{1221}\varepsilon_{21} + C_{1222}\varepsilon_{22}$$

在 Voigt 记法中（$\varepsilon_3^{\text{Voigt}} = \varepsilon_{12}+\varepsilon_{21} = 2\varepsilon_{12}$）：
$$\sigma_3 = C_{31}\varepsilon_1 + C_{33}\varepsilon_3 + C_{32}\varepsilon_2$$

两种表达一致性由 $C_{1212} = C_{1221}$（对称性）保证。

> 🔗 **跨章连接**：Voigt 规则在第 5 章"FEM 公式推导"和第 6 章"单元构造"中反复使用——每个单元的 $\mathbf{B}$ 矩阵（几何矩阵）和 $\mathbf{D}$ 矩阵（弹性矩阵）都是按 Voigt 规则构建的。特别是应变-位移矩阵 $\mathbf{B}$：$\{\varepsilon\} = \mathbf{B}\{u\}$ 中的 $\mathbf{B}$ 就是 $[\partial]$ 算子（见下一节）作用于形函数的矩阵形式。

#### 3.5.2 矩阵算子形式

为了便于有限元方程离散化，引入矩阵算子表达式。

**一阶微分算子矩阵** $[\partial]$：

$$[\partial] = \begin{bmatrix} \frac{\partial}{\partial x} & 0 & 0 \\ 0 & \frac{\partial}{\partial y} & 0 \\ 0 & 0 & \frac{\partial}{\partial z} \\ 0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\ \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \end{bmatrix}_{6 \times 3}$$

则三大基本方程的矩阵算子形式为：

**平衡方程**：$\boxed{[\partial]^T\{\boldsymbol{\sigma}\} + \{F\} = 0}$

**几何方程**：$\boxed{\{\boldsymbol{\varepsilon}\} = [\partial]\{\mathbf{u}\}}$

**本构方程**：$\boxed{\{\boldsymbol{\sigma}\} = [\mathbf{C}]\{\boldsymbol{\varepsilon}\}}$

> 这组关系在后续的有限元离散中反复使用。注意 $\{\boldsymbol{\sigma}\}$ 和 $\{\boldsymbol{\varepsilon}\}$ 已按 Voigt 规则写为 $6\times 1$ 列向量。


#### 3.5.3 柱坐标系下的控制方程

普通弹性力学教材通常采用"取微元体，列平衡"的方法直接在柱坐标系下推导方程。这里介绍更系统的方法——**通过坐标变换**从直角坐标系方程直接转换到柱坐标系。

**坐标变换关系**：
$$x = r\cos\theta, \quad y = r\sin\theta, \quad z = z$$

其逆变换：
$$r = \sqrt{x^2 + y^2}, \quad \theta = \tan^{-1}(y/x), \quad z = z$$

**偏导数算子变换**：
$$\frac{\partial}{\partial x} = \cos\theta\frac{\partial}{\partial r} - \frac{\sin\theta}{r}\frac{\partial}{\partial\theta}$$
$$\frac{\partial}{\partial y} = \sin\theta\frac{\partial}{\partial r} + \frac{\cos\theta}{r}\frac{\partial}{\partial\theta}$$
$$\frac{\partial}{\partial z} = \frac{\partial}{\partial z}$$

矩阵形式：
$$\begin{Bmatrix} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y} \\ \frac{\partial}{\partial z} \end{Bmatrix} = \begin{bmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{Bmatrix} \frac{\partial}{\partial r} \\ \frac{1}{r}\frac{\partial}{\partial\theta} \\ \frac{\partial}{\partial z} \end{Bmatrix}$$

**位移变换**：
$$\begin{Bmatrix} u \\ v \\ w \end{Bmatrix} = \begin{bmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{Bmatrix} u_r \\ u_\theta \\ u_z \end{Bmatrix}$$

其中 $u_r, u_\theta, u_z$ 分别为径向、周向和轴向位移分量。

**柱坐标系下的平衡方程**（最终结果）：
$$\begin{aligned}
\frac{\partial\sigma_r}{\partial r} + \frac{1}{r}\frac{\partial\tau_{r\theta}}{\partial\theta} + \frac{\partial\tau_{zr}}{\partial z} + \frac{\sigma_r - \sigma_\theta}{r} + F_r = 0 \\
\frac{\partial\tau_{r\theta}}{\partial r} + \frac{1}{r}\frac{\partial\sigma_\theta}{\partial\theta} + \frac{\partial\tau_{\theta z}}{\partial z} + \frac{2\tau_{r\theta}}{r} + F_\theta = 0 \\
\frac{\partial\tau_{zr}}{\partial r} + \frac{1}{r}\frac{\partial\tau_{z\theta}}{\partial\theta} + \frac{\partial\sigma_z}{\partial z} + \frac{\tau_{zr}}{r} + F_z = 0
\end{aligned}$$

**柱坐标系下的几何方程**（部分结果）：
$$\varepsilon_r = \frac{\partial u_r}{\partial r}, \quad \varepsilon_\theta = \frac{1}{r}\frac{\partial u_\theta}{\partial\theta} + \frac{u_r}{r}, \quad \varepsilon_z = \frac{\partial u_z}{\partial z}$$
$$\gamma_{r\theta} = \frac{1}{r}\frac{\partial u_r}{\partial\theta} + \frac{\partial u_\theta}{\partial r} - \frac{u_\theta}{r}$$

> **注意**：柱坐标下的平衡方程比直角坐标多了两个"附加项"：$\frac{\sigma_r - \sigma_\theta}{r}$ 和 $\frac{2\tau_{r\theta}}{r}$。这些附加项来自坐标系的曲率——在推导过程中不需要刻意保留任何微小量，坐标变换法自动处理。

> 💡 **理解关键**：柱坐标平衡方程中的"附加项"本质上是因为坐标系的基向量随位置变化。直角坐标系中 $\mathbf{e}_i$ 是常向量（导数=0），柱坐标中 $\mathbf{e}_r$ 和 $\mathbf{e}_\theta$ 随 $\theta$ 变化。所以"应力散度"项中出现了坐标曲率引起的附加项——这在物理上是合理的：一个纯环向应力 $\sigma_\theta$ 在径向产生净力（因为微小弧段的应力方向不完全相反）。

---

> **思考题**：球坐标系下的基本方程可以通过两次柱坐标变换获得。试推导球坐标系下的平衡方程。

---

### 3.6 结论

建立并求解弹性力学的边值问题构成了弹性理论的核心。边值问题将具体的物理问题转化为数学问题。

弹性力学边值问题的求解任务是：**求解 15 个未知量以满足 15 个方程和边界条件**。

**三类求解方法**：

| 方法 | 典型技术 | 特点 |
|------|---------|------|
| **解析法** | 偏微分方程、复变函数、积分方程、群论、微分几何、泛函分析 | 精确解，但适用范围有限 |
| **实验法** | 电测法、光弹性法、类比法、云纹法、全息法 | 直观可信，但成本高、通用性差 |
| **数值分析法** | 有限差分法、**有限元法**、边界元法、无网格法 | 通用性强，适合复杂问题 |

**有限元分析是最重要和实用的方法**，已被广泛应用于科学研究和工程计算。为了深入理解有限元分析，接下来进入经典变分法的基础理论。

> 🔗 **跨章连接**：三类求解方法的对比是考试选择题常客。"有限元法属于哪种方法？"→数值分析法。"弹性力学问题的主要数值方法有哪些？"→有限差分法、有限元法、边界元法、无网格法。

---

## 三大基本方程速查表

| 方程名称 | 张量形式 | 展开形式（直角坐标） | 方程数 |
|---------|---------|-------------------|--------|
| **平衡方程** | $\sigma_{ij,j} + F_i = 0$ | $\frac{\partial\sigma_{xx}}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{xz}}{\partial z} + F_x = 0$ 等 | 3 |
| **几何方程** | $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$ | $\varepsilon_{xx} = \frac{\partial u}{\partial x}, \gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}$ 等 | 6 |
| **本构方程** | $\sigma_{ij} = \lambda e\delta_{ij} + 2\mu\varepsilon_{ij}$ | $\sigma_{xx} = \lambda e + 2\mu\varepsilon_{xx}, \tau_{xy} = \mu\gamma_{xy}$ 等 | 6 |

---

## 检查你的理解（综合）

1. 弹性力学的求解框架中，15 个未知量如何通过 15 个方程确定？哪些已知量出现在方程中？
2. Cauchy 应力公式、应力变换公式和应力边界条件之间有什么关系？用指标记号写出这三者的异同。
3. Green 应变 $E_{ij} = \frac{1}{2}(\xi_{k,i}\xi_{k,j} - \delta_{ij})$ 如何退化为小变形应变 $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$？退化过程中哪些项被略去了？
4. 各向同性线弹性材料为什么最终只有 2 个独立弹性常数？说出从 81 $\to$ 36 $\to$ 21 $\to$ 2 每一步简化的物理/数学依据。
5. Voigt 规则中为什么运动量（应变）的剪切项要乘以 2？如果不乘 2 会有什么后果？

---

> **对应作业**：[HW1 Q1（恒等式证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q2（张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q3（$e_{ijk}$ 证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
>
> **下一章**：[变分法基础](../01-Lecture-Notes/1-3-Variational-Methods.md) — 弹性力学的变分表述是有限元方法的理论基础。
