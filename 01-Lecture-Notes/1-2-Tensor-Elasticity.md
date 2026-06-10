# 张量分析与弹性力学

> **对应课件**：[`Chapter 3 Variation theory and applications-1.pdf`](../06-References/pdfs-originals/Chapter%203%20Variation%20theory%20and%20applications-1.pdf) 课程第 1 章
> **章节定位**：Solid mechanics: theory of elasticity → Introduction → Introduction to tensor → Index and summation → Common symbols → Revision of linear algebra → Cartesian tensor
> **相关作业**：[HW1 Q1-Q4](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：线性代数（向量、矩阵、行列式）、高等数学（偏导数）

---

## 2.1 弹性力学发展简介

弹性力学本身是一门具有 300 多年历史的学科。1667 年 Hooke 发现弹性体的位移与外力成正比。1821 年和 1823 年，Navier 和 Cauchy 分别推导出线弹性边值问题的平衡方程。1855 年 Saint-Venant 的扭转和弯曲理论，1933 年 Muskhelishvili 的复变函数方法等都是弹性力学的经典工作。

弹性力学需要强有力的数学工具，而**张量分析（tensor analysis）**正是描述弹性力学问题的最简洁、最系统的语言。

---

## 2.2 张量引论（Introduction to Tensor）

### 2.1 Index and summation（指标与求和）

#### Einstein 求和约定

在张量运算中，**在同一项中重复出现的下标表示对该下标求和**：
$$a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^3 a_i b_i$$

其中 $i$ 是**哑指标（dummy index）**。求和约定大大简化了书写。

**例**：
$$c_{ik} = a_{ij}b_{jk} = \sum_{j=1}^3 a_{ij}b_{jk} = a_{i1}b_{1k} + a_{i2}b_{2k} + a_{i3}b_{3k}$$

#### 自由指标（Free index）
不重复的指标表示分量的编号。方程两边的自由指标必须一致。

**例**：$a_{ij} = b_{ik}c_{kj}$
- $i, j$：自由指标 → 代表 9 个方程
- $k$：哑指标 → 对 $k=1,2,3$ 求和

### 2.2 Common symbols（常用符号）

#### Kronecker 符号 $\delta_{ij}$

$$\delta_{ij} = \begin{cases} 1 & i=j \\ 0 & i \neq j \end{cases}$$

**重要性质**：
- $\delta_{ii} = 3$（迹）
- $\delta_{ij}a_j = a_i$（指标替换）
- $\delta_{ij}\delta_{jk} = \delta_{ik}$

#### 置换符号（Permutation symbol）$e_{ijk}$

$$e_{ijk} = \begin{cases}
0 & \text{任意两指标相等} \\
+1 & (i,j,k)\text{ 为偶排列：(1,2,3),(2,3,1),(3,1,2)} \\
-1 & (i,j,k)\text{ 为奇排列：(1,3,2),(2,1,3),(3,2,1)}
\end{cases}$$

**核心恒等式**：
$$e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$$
$$e_{ijk}e_{ijt} = 2\delta_{kt}$$
$$e_{ijk}e_{ijk} = 6$$

向量叉积的张量形式：
$$(\mathbf{a}\times\mathbf{b})_i = e_{ijk}a_jb_k$$

### 2.3 Revision of linear algebra（线性代数复习）

#### 向量与坐标

在笛卡尔坐标系 $\{O; x_1,x_2,x_3\}$ 中，向量 $\mathbf{r}$ 可表示为：
$$\mathbf{r} = x_1\mathbf{e}_1 + x_2\mathbf{e}_2 + x_3\mathbf{e}_3$$

#### 基向量变换

设有新坐标系 $\{O; x_1',x_2',x_3'\}$，方向余弦 $l_{i'j} = \cos(x_i', x_j)$：
$$\mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j \quad\text{或}\quad \mathbf{e}' = \mathbf{L}\mathbf{e}$$

L 为正交矩阵（$\mathbf{L}^{-1} = \mathbf{L}^T$），逆变换为 $\mathbf{e} = \mathbf{L}^T\mathbf{e}'$。

#### 向量分量变换

旧到新：$a_{i'} = l_{i'j}a_j$
新到旧：$a_i = l_{ij'}a_{j'}$

#### 矩阵变换

线性变换 $\Psi$ 在基 $\{\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3\}$ 下的矩阵为 $\mathbf{A}$：
$$\Psi\mathbf{a} = \mathbf{e}^T\mathbf{A}\mathbf{a}$$

在新基下：$\Psi\mathbf{a} = \mathbf{e}'^T\mathbf{A}'\mathbf{a}'$，可得 $\mathbf{A}' = \mathbf{L}\mathbf{A}\mathbf{L}^T$。

#### 标量积与向量积

标量积（dot product）：$\mathbf{a}\cdot\mathbf{b} = a_ib_i$

向量积（cross product）：$(\mathbf{a}\times\mathbf{b})_i = e_{ijk}a_jb_k$

### 2.4 Cartesian tensor（笛卡尔张量）

#### 标量（零阶张量）

坐标变换下保持不变：$\phi' = \phi$

#### 向量（一阶张量）

按向量分量变换规则变换：$a_{i'}' = l_{i'j}a_j$

**例**：验证向量 $\mathbf{a}$ 确实是一阶张量——即 $a_{i'}' = l_{i'j}a_j$。

事实上 $\mathbf{a} = a_j\mathbf{e}_j = a_{i'}'\mathbf{e}_{i'}'$，代入基向量变换：
$$a_j\mathbf{e}_j = a_{i'}'l_{i'j}\mathbf{e}_j \Rightarrow a_j = l_{i'j}a_{i'}' \Rightarrow a_{i'}' = l_{i'j}a_j$$

#### 二阶张量

按以下规则变换：
$$a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$$

**例**：应力 $\sigma_{ij}$ 和应变 $\varepsilon_{ij}$ 都是二阶张量。

#### 张量运算

1. **加法**：同阶张量对应分量相加，$C_{ij} = A_{ij} + B_{ij}$
2. **张量积（外积）**：$C_{ijkl} = A_{ij}B_{kl}$，升 2 阶
3. **缩并（Contraction）**：令两个指标相同并求和，$A_{ij} \to A_{ii} = \text{tr}(\mathbf{A})$，降 2 阶
4. **张量的分解**：$A_{ij} = \frac12(A_{ij}+A_{ji}) + \frac12(A_{ij}-A_{ji})$（对称部分 + 反对称部分）
5. **张量的不变量**：$I_1 = A_{ii}$，$I_2 = \frac12(A_{ii}A_{jj} - A_{ij}A_{ji})$，$I_3 = \det(\mathbf{A})$

---

## 2.3 弹性力学的张量表示

| 方程 | 张量形式 | 分量个数 |
|------|----------|----------|
| 几何方程 | $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ | 6 |
| 本构方程 | $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$ | 6 |
| 平衡方程 | $\sigma_{ij,j} + f_i = 0$ | 3 |

---

## 检查你的理解

1. 哑指标和自由指标有什么区别？$a_{ij}b_{jk} = c_{ik}$ 中各包含哪些？
2. $\delta_{ij}\delta_{ij}$ 等于多少？
3. 证明 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 满足二阶张量变换律。
4. 为什么 $e_{ijk}$ 对 $\det(\mathbf{L})=-1$ 的变换被称为赝张量？

---

> **对应作业**：[HW1 Q1（恒等式证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q2（张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q3（$e_{ijk}$ 证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
