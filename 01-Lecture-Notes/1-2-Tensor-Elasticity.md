# 第2章：张量分析与弹性力学

> **对应 PDF**：[`Chapter 2 Elastic theory.pdf`](../06-References/pdfs-originals/Chapter%202%20Elastic%20theory.pdf)
> **相关作业**：[HW1 Q1-Q4](../04-Homework-Solutions/2026w/HW1-Problem.md)（全为张量题）
> **前置知识**：线性代数（向量、矩阵乘法、坐标变换、行列式）、高等数学（偏导数、链式法则）

---

## 2.1 为什么学习张量？

弹性力学涉及大量的坐标变换、多分量物理量、以及在曲线坐标系下的微分运算。用分量形式逐项书写不仅冗长，而且容易出错。张量记号（tensor notation）提供了一种**简洁、统一**的语言：

| 不用张量 | 用张量后 |
|----------|----------|
| 几何方程：$\varepsilon_x = \frac{\partial u}{\partial x}, \cdots$（6个式子） | $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$（1个式子） |
| 平衡方程：$\frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{xz}}{\partial z} + f_x = 0, \cdots$ | $\sigma_{ij,j} + f_i = 0$ |
| 本构关系：需要写出 $6\times6$ 矩阵元素 | $\sigma_{ij} = D_{ijkl}\varepsilon_{kl}$ |
| 坐标变换：逐项展开 9 个分量 | $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$ |

> 张量分析的核心洞察：物理量本身（如应力 $\boldsymbol{\sigma}$）是与坐标系选取**无关**的客观实体，它的**分量**才随坐标系变化而变化。

---

## 2.2 指标记号约定

### 2.2.1 自由指标与哑指标

- **自由指标**（Free Index）：在表达式左右两侧各出现一次 → 表示分量的编号
- **哑指标**（Dummy Index）：在同一项中出现两次 → 表示对该指标从 1 到 3 求和（省略了求和符号 $\sum$）

**示例**：$a_{ij} = b_{ik}c_{kj}$
- $i, j$：自由指标 → 这代表 $3\times3 = 9$ 个方程
- $k$：哑指标 → 对 $k=1,2,3$ 求和

**展开验证**（以 $i=1, j=1$ 为例）：
$$a_{11} = \sum_{k=1}^3 b_{1k}c_{k1} = b_{11}c_{11} + b_{12}c_{21} + b_{13}c_{31}$$

### 2.2.2 Einstein 求和约定

同一项中重复出现的下标表示对该下标求和，省略求和符号 $\sum$：
$$a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^3 a_i b_i$$

> **重要规则**：哑指标不能在同一项中出现**三次及以上**。$a_i b_i c_i$ 是无意义的，因为无法确定对哪些下标求和。

---

## 2.3 两个重要的特殊符号

### 2.3.1 Kronecker $\delta$ 符号

$$ \delta_{ij} = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases} $$

几何意义：单位矩阵的元素 $I_{ij} = \delta_{ij}$。

**重要性质**（需要能熟练运用）：

| 性质 | 含义 |
|------|------|
| $\delta_{ii} = 3$ | 迹（trace）等于空间维数 |
| $\delta_{ij}a_j = a_i$ | 指标替换——把 $j$ 替换成 $i$ |
| $\delta_{ij}\delta_{jk} = \delta_{ik}$ | 缩并 |
| $\delta_{ij}\delta_{ij} = 3$ | 双缩并 |
| $\delta_{ij}A_{jk} = A_{ik}$ | 矩阵乘法中的单位元 |

### 2.3.2 置换符号 $e_{ijk}$（Levi-Civita 符号）

$$ e_{ijk} = \begin{cases}
0 & \text{任意两个指标相等} \\
+1 & (i,j,k) \text{ 是 } (1,2,3) \text{ 的偶排列} \\
-1 & (i,j,k) \text{ 是 } (1,2,3) \text{ 的奇排列}
\end{cases} $$

偶排列（$e_{ijk}=1$）：$(1,2,3), (2,3,1), (3,1,2)$ —— **循环置换**
奇排列（$e_{ijk}=-1$）：$(1,3,2), (2,1,3), (3,2,1)$ —— **交换一次**

**重要恒等式**（考试必考！）：

$$ e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks} $$

这是 $\varepsilon$-$\delta$ 恒等式的标准形式。推导思路：将相同指标放在最后（利用循环置换），然后应用标准恒等式。

$$ e_{ijk}e_{ijt} = 2\delta_{kt} $$

向量叉积的张量形式：
$$ (\mathbf{a}\times\mathbf{b})_i = e_{ijk}a_jb_k $$

**证明方案**：以标准恒等式 $e_{ijk}e_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}$ 为例。注意到两个 $e$ 在第一个指标上缩并，结果只与剩下的两个指标有关。通过枚举所有可能情况（共 $3^4=81$ 种，但多数为零）可证实。

---

## 2.4 坐标变换与张量定义

### 2.4.1 基向量变换

设新旧两个笛卡尔直角坐标系 $\{O; x_1,x_2,x_3\}$ 和 $\{O; x_1',x_2',x_3'\}$，它们共享原点。

方向余弦（direction cosine）：
$$ l_{i'j} = \cos(x_i', x_j) $$

基向量变换：
$$ \mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j \quad \text{或矩阵形式} \quad \mathbf{e}' = \mathbf{L}\mathbf{e} $$

由于 $\mathbf{L}$ 是正交矩阵（$\mathbf{L}^{-1} = \mathbf{L}^T$），逆变换为：
$$ \mathbf{e} = \mathbf{L}^T\mathbf{e}' $$

### 2.4.2 分量变换

一个向量 $\mathbf{a}$ 在新旧坐标系下可以表示为：
$$ \mathbf{a} = a_j\mathbf{e}_j = a_{j'}'\mathbf{e}_{j'}' $$

将基向量变换代入：
$$ a_j\mathbf{e}_j = a_{j'}'l_{j'j}\mathbf{e}_j \quad \Rightarrow \quad a_j = l_{j'j}a_{j'}' $$

即 $$ a_{i'} = l_{i'j}a_j $$

### 2.4.3 张量的定义（引入法）

一个变量 $T$ 是 $n$ 阶张量，当且仅当其分量在坐标变换下按以下规律变换：

| 阶数 | 名称 | 变换规律 | 例子 |
|------|------|----------|------|
| 0 | 标量 | $\phi' = \phi$ | 温度、密度 |
| 1 | 向量 | $a_{i'}' = l_{i'j}a_j$ | 位移、力 |
| 2 | 二阶张量 | $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$ | 应力、应变 |
| 3 | 三阶张量 | $a_{i'j'k'}' = l_{i'i}l_{j'j}l_{k'k}a_{ijk}$ | $e_{ijk}$（赝张量） |

### 2.4.4 张量运算

| 运算 | 公式 | 阶数变化 | 说明 |
|------|------|----------|------|
| 加法 | $C_{ij} = A_{ij} + B_{ij}$ | 同阶 | 对应分量相加 |
| 标量乘 | $B_{ij} = \lambda A_{ij}$ | 不变 | 每个分量乘常数 |
| 张量积（外积） | $C_{ijkl} = A_{ij}B_{kl}$ | 升 2 阶 | 从 $m$ 阶到 $m+n$ 阶 |
| 缩并 | $A_{ii}$ | 降 2 阶 | 令两个指标相等并求和 |
| 内积 | $A_{ij}b_j = c_i$ | 降 1 阶 | 对乘积缩并 |

### 2.4.5 张量运算的几何意义

- **缩并 $A_{ii}$** = 求矩阵的迹（trace）→ 标量不变量
- **张量对称分解**：任意二阶张量可分解为对称部分和反对称部分
  $$A_{ij} = \underbrace{\frac12(A_{ij}+A_{ji})}_{A_{(ij)}\text{ 对称}} + \underbrace{\frac12(A_{ij}-A_{ji})}_{A_{[ij]}\text{ 反对称}}$$
  应力张量 $\sigma_{ij}$ 和应变张量 $\varepsilon_{ij}$ 都**天然对称**（$\sigma_{ij} = \sigma_{ji}$）
- **张量的不变量**（坐标变换下保持不变）：
  - 第一不变量（迹）：$I_1 = A_{ii}$
  - 第二不变量：$I_2 = \frac12(A_{ii}A_{jj} - A_{ij}A_{ji})$
  - 第三不变量：$I_3 = \det(\mathbf{A})$

---

## 2.5 习题级推导：如何证明一个量是张量？

这是考试的标准题型。通用策略：

**Step 1**：写出坐标变换 $x_{i'} = l_{i'j}x_j$
**Step 2**：写出该量在新旧坐标系中的定义
**Step 3**：利用链式法则变换导数 $\partial/\partial x_{j'} = l_{mj'}\,\partial/\partial x_m$
**Step 4**：代入整理成 $l_{i'i}l_{j'j}\cdots$ 乘以原分量的形式
**Step 5**：与张量变换律对比 → 得证

### 例 1：证明 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 是二阶张量

**证明**：
$$\begin{aligned}
\varepsilon_{i'j'} &= \frac12(u_{i',j'} + u_{j',i'}) \\
&= \frac12\left[ \frac{\partial}{\partial x_{j'}}(l_{i'i}u_i) + \frac{\partial}{\partial x_{i'}}(l_{j'j}u_j) \right] \\
&= \frac12\left[ l_{i'i}l_{mj'}u_{i,m} + l_{j'j}l_{ni'}u_{j,n} \right]
\end{aligned}$$

将第二项的哑指标 $(j,n)\to(i,m)$：
$$ \varepsilon_{i'j'} = \frac12\left[ l_{i'i}l_{mj'} + l_{j'i}l_{mi'} \right] u_{i,m} $$

利用 $l_{mj'} = l_{j'm}$ 和 $l_{mi'} = l_{i'm}$：
$$ \varepsilon_{i'j'} = \frac12\left[ l_{i'i}l_{j'm} + l_{j'i}l_{i'm} \right] u_{i,m} $$

**关键步骤**：注意到 $u_i$ 是向量（一阶张量），$u_{i,m}$ 是二阶张量的分量。我们希望将其表达为：
$$ \varepsilon_{i'j'} = l_{i'i}l_{j'j}\varepsilon_{ij} $$

而 $\varepsilon_{ij} = \frac12(u_{i,j} + u_{j,i})$，代入右侧：
$$ l_{i'i}l_{j'j}\varepsilon_{ij} = \frac12\left(l_{i'i}l_{j'j} + l_{i'j}l_{j'i}\right)u_{i,j} $$

交换第二项中的哑指标名字 $(j,i)\to(m,i)$：
$$ = \frac12\left(l_{i'i}l_{j'm} + l_{i'm}l_{j'i}\right)u_{i,m} $$

与左侧式比较，完全相同。因此 $\varepsilon_{i'j'} = l_{i'i}l_{j'j}\varepsilon_{ij}$。■

### 例 2：证明 $e_{ijk}$ 是三阶张量（对正常转动）

需要验证：
$$ e_{i'j'k'} = l_{i'i}l_{j'j}l_{k'k}e_{ijk} $$

右侧 $l_{i'i}l_{j'j}l_{k'k}e_{ijk}$ 的几何意义是矩阵 $\mathbf{L}$ 三行向量的**标量三重积**，即 $\det(\mathbf{L})$。

对**正常转动**（$\det(\mathbf{L}) = 1$，不包含反射）：
$$ l_{i'i}l_{j'j}l_{k'k}e_{ijk} = \det(\mathbf{L}) \cdot e_{i'j'k'} = e_{i'j'k'} $$

因此满足三阶张量变换律。对包含反射的变换（$\det(\mathbf{L}) = -1$），$e_{ijk}$ 需额外乘以 $\det(\mathbf{L})$ 的符号，因此称为**赝张量（pseudotensor）**。

---

## 2.6 弹性力学的张量形式

使用张量记号后，弹性力学基本方程变得极其简洁：

| 方程 | 矩阵形式（第 1 章） | 张量形式 |
|------|-------------------|----------|
| 几何 | $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$ | $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ |
| 物理 | $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$ | $\sigma_{ij} = D_{ijkl}\varepsilon_{kl}$ 或 $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$ |
| 平衡 | $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$ | $\sigma_{ij,j} + f_i = 0$ |
| 位移 BC | $\mathbf{u}|_{S_u} = \bar{\mathbf{u}}$ | $u_i|_{S_u} = \bar{u}_i$ |
| 力 BC | $[\mathbf{n}]\boldsymbol{\sigma}|_{S_\sigma} = \bar{\mathbf{T}}$ | $\sigma_{ij}n_j|_{S_\sigma} = \bar{T}_i$ |

其中 $D_{ijkl}$ 是四阶弹性张量（$3^4 = 81$ 个分量），但对称性使其缩减为 21 个独立分量（各向异性）→ 2 个独立分量（各向同性）。

---

## 2.7 解题实用技能

### ε-δ 恒等式使用三步法

遇到 $e_{ijk}e_{pqr}$ 形式的乘积时：

1. **对齐指标**：利用循环置换 $e_{ijk} = e_{jki} = e_{kij}$，将相同指标放在最后一个位置
2. **应用恒等式**：$e_{ijk}e_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}$
3. **利用 $\delta$ 简化**：$\delta_{ij}a_j = a_i$，$\delta_{ii} = 3$

### 常见 ε-δ 恒等式速查

| 恒等式 | 记忆法 |
|--------|--------|
| $e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$ | 第一个指标相同，后面的"交叉减交换" |
| $e_{ijk}e_{ijt} = 2\delta_{kt}$ | 前两个指标缩并后还剩两个 |
| $e_{ijk}e_{ijk} = 6$ | 全部缩并 = 3 维空间的全排列数 |

---

## 检查你的理解

1. 什么是哑指标和自由指标？在 $a_i = B_{ij}c_j$ 中，哪些是哑指标，哪些是自由指标？这个式子代表多少个方程？
2. $\delta_{ij}\delta_{ij}$ 等于多少？用求和展开验证。
3. $e_{123}$ 和 $e_{213}$ 分别等于多少？
4. 为什么 $e_{ijk}$ 对 $\det(\mathbf{L}) = -1$ 的变换是赝张量而不是真正张量？
5. 将 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 展开成矩阵形式，验证与第 1 章的几何方程一致。
6. 张量运算中，缩并（contraction）的作用是什么？举例说明。

---

> **对应作业**：[HW1 Q1](../04-Homework-Solutions/2026w/HW1-Problem.md)（恒等式证明）· [Q2](../04-Homework-Solutions/2026w/HW1-Problem.md)（张量证明）· [Q3](../04-Homework-Solutions/2026w/HW1-Problem.md)（$e_{ijk}$ 证明）· [Q4](../04-Homework-Solutions/2026w/HW1-Problem.md)（指标运算）
> **往年参考**：[Homework1.1](../04-Homework-Solutions/past/comprehensive/Homework1.1.md) · [Homework1.2](../04-Homework-Solutions/past/comprehensive/Homework1.2.md)
