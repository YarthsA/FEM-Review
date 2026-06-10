# 第2章：张量分析与弹性力学

> **对应 PDF**：[`Chapter 2 Elastic theory.pdf`](../06-References/pdfs-originals/Chapter%202%20Elastic%20theory.pdf)
> **相关作业**：[HW1 Q1-Q4](../04-Homework-Solutions/2026w/HW1-Problem.md)（全为张量题）
> **前置知识**：线性代数（向量、矩阵、行列式、坐标变换）、高等数学（偏导数、链式法则）

---

## 2.1 为什么要用张量？

### 2.1.1 弹性力学需要张量的原因

弹性力学涉及大量多分量物理量（应力 6 个分量、应变 6 个分量），这些量在不同坐标系之间遵循特定的变换规律。如果我们用分量形式逐一书写，不仅冗长，而且容易迷失在符号的海洋中。

**张量记号（tensor notation）** 的优势：

| 操作 | 不用张量 | 用张量后 |
|------|----------|----------|
| 几何方程 | $\varepsilon_x = \frac{\partial u}{\partial x}, \varepsilon_y = \frac{\partial v}{\partial y}, \cdots$（6 个式子） | $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$（1 行） |
| 平衡方程 | $\frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{xz}}{\partial z} + f_x = 0$ 等 3 式 | $\sigma_{ij,j} + f_i = 0$（1 行） |
| 坐标变换 | 逐项展开 $3\times3=9$ 个分量 | $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$ |
| 本构关系 | 写出 $6\times6$ 矩阵 | $\sigma_{ij} = C_{ijkl}\varepsilon_{kl}$ |

### 2.1.2 物理量 vs 分量

张量分析的核心洞察：**物理量本身（如应力 $\boldsymbol{\sigma}$）是与坐标系选取无关的客观实体**，改变坐标系只改变它的**分量（components）**，不改变物理量本身。

例如，在 $xyz$ 坐标系中应力状态为 $\sigma_{ij}$，在旋转后的 $x'y'z'$ 坐标系中变成了 $\sigma_{i'j'}'$。$\sigma_{ij}$ 和 $\sigma_{i'j'}'$ 是同一个物理应力的不同面貌，它们通过一定的变换规律相互联系。

---

## 2.2 指标记号与求和约定

### 2.2.1 Einstein 求和约定

在张量运算中，**重复出现的下标表示对该下标从 1 到 3 求和**，求和符号 $\sum$ 省略不写：

$$a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^3 a_i b_i$$

$$a_{ij} b_j = a_{i1} b_1 + a_{i2} b_2 + a_{i3} b_3$$

> **规则**：同一项中同一个下标不能出现三次及以上。例如 $a_i b_i c_i$ 是无意义的。

### 2.2.2 哑指标与自由指标

- **哑指标（Dummy Index）**：同一项中重复出现的下标 → 对该指标求和
- **自由指标（Free Index）**：同一项中只出现一次的下标 → 表示分量的编号

**例**：$a_{ij} = b_{ik} c_{kj}$
- $i, j$：自由指标 → 这个式子代表 $3\times3 = 9$ 个方程
- $k$：哑指标 → 对 $k=1,2,3$ 求和

**展开**（取 $i=1, j=1$）：
$$a_{11} = \sum_{k=1}^3 b_{1k}c_{k1} = b_{11}c_{11} + b_{12}c_{21} + b_{13}c_{31}$$

### 2.2.3 练习

判断以下表达式的含义：

| 表达式 | 方程个数 | 求和指标 |
|--------|----------|----------|
| $a_i = b_{ij}c_j$ | 3 | $j$ |
| $A_{ij} = B_{ik}C_{kj}$ | 9 | $k$ |
| $a_i b_i = c$ | 1 | $i$ |
| $A_{ii} = B_{ij}C_{ji}$ | 1 | $i, j$ 都求和 |

---

## 2.3 Kronecker $\delta$ 符号

### 2.3.1 定义

$$\delta_{ij} = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}$$

即单位矩阵的元素：$I_{ij} = \delta_{ij}$。

### 2.3.2 重要性质

| 性质 | 含义 |
|------|------|
| $\delta_{ii} = 3$ | 迹 = 空间维数（1+1+1=3） |
| $\delta_{ij} a_j = a_i$ | **指标替换功能**——把 $j$ 换成 $i$ |
| $\delta_{ij} \delta_{jk} = \delta_{ik}$ | 缩并后仍是 $\delta$ |
| $\delta_{ij} \delta_{ij} = 3$ | 双缩并 |
| $\delta_{ij} A_{jk} = A_{ik}$ | 矩阵乘法中的单位元 |

**例**：验证 $\delta_{ij}a_j = a_i$
$$\delta_{1j}a_j = \delta_{11}a_1 + \delta_{12}a_2 + \delta_{13}a_3 = 1\cdot a_1 + 0 + 0 = a_1$$
$$\delta_{2j}a_j = \delta_{21}a_1 + \delta_{22}a_2 + \delta_{23}a_3 = 0 + 1\cdot a_2 + 0 = a_2$$
$$\delta_{3j}a_j = a_3$$

---

## 2.4 置换符号 $e_{ijk}$（Levi-Civita 符号）

### 2.4.1 定义

$$e_{ijk} = \begin{cases}
0 & \text{任意两指标相等} \\
+1 & (i,j,k) \text{ 是 } (1,2,3) \text{ 的偶排列} \\
-1 & (i,j,k) \text{ 是 } (1,2,3) \text{ 的奇排列}
\end{cases}$$

**偶排列**（循环置换）：$(1,2,3), (2,3,1), (3,1,2)$ → $e=1$
**奇排列**（交换一次）：$(1,3,2), (2,1,3), (3,2,1)$ → $e=-1$

### 2.4.2 向量运算中的用途

**叉积**的张量形式：
$$(\mathbf{a}\times\mathbf{b})_i = e_{ijk} a_j b_k$$

**展开**（$i=1$）：
$$(\mathbf{a}\times\mathbf{b})_1 = e_{123}a_2b_3 + e_{132}a_3b_2 = 1\cdot a_2b_3 + (-1)\cdot a_3b_2 = a_2b_3 - a_3b_2$$

这正是一阶行列式的形式。✅

### 2.4.3 核心恒等式（考试必考！）

**恒等式 1**（最常用）：
$$e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$$

**恒等式 2**（前两个指标缩并）：
$$e_{ijk}e_{ijt} = 2\delta_{kt}$$

**恒等式 3**（全部缩并）：
$$e_{ijk}e_{ijk} = 6$$

**记忆方法**：
- 第一个指标相同（$i$），后面的 "交叉相乘减交换"
- "交叉"：$j$→$s$，$k$→$t$
- "交换"：$j$→$t$，$k$→$s$

### 2.4.4 恒等式的应用步骤

遇到 $e_{ijk}e_{pqr}$ 形式的乘积：
1. **对齐指标**：利用 $e_{ijk} = e_{jki} = e_{kij}$（循环置换），将相同指标放在同一位置
2. **应用恒等式**：代入 $e_{ijk}e_{imn} = \delta_{jm}\delta_{kn} - \delta_{jn}\delta_{km}$
3. **利用 $\delta$ 简化**：$\delta_{ij}a_j = a_i$，$\delta_{ii} = 3$

---

## 2.5 坐标变换

### 2.5.1 基向量变换

设有新旧两个笛卡尔直角坐标系 $\{x_1,x_2,x_3\}$ 和 $\{x_1',x_2',x_3'\}$，共享原点。

定义方向余弦 $l_{i'j} = \cos(x_i', x_j)$，即 $x_i'$ 轴与 $x_j$ 轴夹角的余弦。

**基向量变换**：
$$\mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j \quad \text{或矩阵形式} \quad \mathbf{e}' = \mathbf{L}\mathbf{e}$$

其中 $\mathbf{L}$ 是方向余弦矩阵（$3\times3$），其 $i'$ 行 $j$ 列元素为 $l_{i'j}$。

### 2.5.2 正交性

方向余弦矩阵是**正交矩阵**：
$$\mathbf{L}^{-1} = \mathbf{L}^T \quad \Rightarrow \quad l_{ji'} = l_{i'j}$$

这一性质非常重要，意味着逆变换只需转置：
$$\mathbf{e} = \mathbf{L}^T\mathbf{e}'$$

### 2.5.3 向量分量的变换

同一个向量 $\mathbf{a}$ 在新旧坐标系中表达：
$$\mathbf{a} = a_j\mathbf{e}_j = a_{j'}'\mathbf{e}_{j'}'$$

代入基向量变换：
$$a_j\mathbf{e}_j = a_{j'}'l_{j'j}\mathbf{e}_j \quad\Rightarrow\quad a_j = l_{j'j}a_{j'}'$$

**旧→新**：$a_{i'} = l_{i'j}a_j$
**新→旧**：$a_i = l_{ij'}a_{j'} = l_{j'i}a_{j'}$

---

## 2.6 张量的定义

### 2.6.1 定义（变换律法）

一个量 $T$ 是 $n$ 阶张量，当且仅当其分量在坐标变换下按以下规律变换：

| 阶数 | 名称 | 变换规律 | 例子 |
|------|------|----------|------|
| 0 | 标量 | $\phi' = \phi$ | 温度、密度 |
| 1 | 向量 | $a_{i'}' = l_{i'j}a_j$ | 位移、力 |
| 2 | 二阶张量 | $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$ | 应力、应变 |
| 3 | 三阶张量 | $a_{i'j'k'}' = l_{i'i}l_{j'j}l_{k'k}a_{ijk}$ | $e_{ijk}$（赝张量）|

### 2.6.2 张量运算

| 运算 | 公式 | 阶数变化 |
|------|------|----------|
| 加法 | $C_{ij} = A_{ij} + B_{ij}$ | 不变 |
| 张量积（外积） | $C_{ijkl} = A_{ij}B_{kl}$ | 升 2 阶 |
| 缩并 | $A_{ii} = \text{tr}(A)$ | 降 2 阶 |
| 内积 | $A_{ij}b_j = c_i$ | 降 1 阶 |

### 2.6.3 对称与反对称分解

任意二阶张量可唯一分解为：
$$A_{ij} = \underbrace{\frac12(A_{ij}+A_{ji})}_{A_{(ij)}\text{ 对称}} + \underbrace{\frac12(A_{ij}-A_{ji})}_{A_{[ij]}\text{ 反对称}}$$

- 应力张量 $\sigma_{ij}$：**对称**（$\sigma_{ij}=\sigma_{ji}$）
- 应变张量 $\varepsilon_{ij}$：**对称**（$\varepsilon_{ij}=\varepsilon_{ji}$）

### 2.6.4 张量的不变量

坐标变换下保持不变：
- $I_1 = A_{ii}$（迹）
- $I_2 = \frac12(A_{ii}A_{jj} - A_{ij}A_{ji})$
- $I_3 = \det(\mathbf{A})$

---

## 2.7 张量证明通用模板

考试中常要求证明某个量是张量。标准策略：

**Step 1**：写出坐标变换 $x_{i'} = l_{i'j}x_j$
**Step 2**：写出该量在新旧坐标系中的定义
**Step 3**：利用链式法则 $\partial/\partial x_{j'} = l_{mj'}\,\partial/\partial x_m$
**Step 4**：整理成 $l_{i'i}l_{j'j}\cdots$ 乘以原分量的形式
**Step 5**：与张量变换律对比 → 得证

### 例 1：证明 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 是二阶张量

$$\begin{aligned}
\varepsilon_{i'j'} &= \frac12(u_{i',j'} + u_{j',i'}) \\
&= \frac12[l_{i'i}l_{mj'}u_{i,m} + l_{j'j}l_{ni'}u_{j,n}] \\
&= \frac12[l_{i'i}l_{j'm} + l_{j'i}l_{i'm}]u_{i,m}
\end{aligned}$$

其中利用了 $l_{mj'} = l_{j'm}$（正交性）。另一方面：
$$l_{i'i}l_{j'j}\varepsilon_{ij} = \frac12(l_{i'i}l_{j'j} + l_{i'j}l_{j'i})u_{i,j}$$

两式相等，故 $\varepsilon_{i'j'} = l_{i'i}l_{j'j}\varepsilon_{ij}$。■

### 例 2：证明 $e_{ijk}$ 是三阶张量

需要验证 $e_{i'j'k'} = l_{i'i}l_{j'j}l_{k'k}e_{ijk}$。

右侧 = $\det(\mathbf{L})\,e_{i'j'k'}$（行列式的定义）。$\det(\mathbf{L}) = 1$（正常转动）时成立；$\det(\mathbf{L}) = -1$（含反射）时多一个负号 → 称为**赝张量**。

---

## 2.8 弹性力学的张量形式

| 方程 | 矩阵形式 | 张量形式 |
|------|----------|----------|
| 几何 | $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$ | $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ |
| 物理 | $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$ | $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$ |
| 平衡 | $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$ | $\sigma_{ij,j} + f_i = 0$ |

**体积应变**：$\varepsilon_v = \varepsilon_{kk} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}$

---

## 2.9 恒等式速查表

| 恒等式 | 用途 |
|--------|------|
| $\delta_{ii} = 3$ | 计算迹 |
| $\delta_{ij}a_j = a_i$ | 指标替换 |
| $e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$ | **最常用！** |
| $e_{ijk}e_{ijt} = 2\delta_{kt}$ | 双重缩并 |
| $(\boldsymbol{a}\times\boldsymbol{b})\cdot(\boldsymbol{c}\times\boldsymbol{d}) = (\boldsymbol{a}\cdot\boldsymbol{c})(\boldsymbol{b}\cdot\boldsymbol{d}) - (\boldsymbol{a}\cdot\boldsymbol{d})(\boldsymbol{b}\cdot\boldsymbol{c})$ | Lagrange 恒等式 |

---

## 2.10 张量在 FEM 中的实际作用

### 2.10.1 从分量运算到抽象索引

在 FEM 编程中，很少直接使用张量符号。但张量符号是**推导公式**的利器。

**例**：弹性矩阵 $D_{ijkl}$ 的对称性

$$D_{ijkl} = \lambda\delta_{ij}\delta_{kl} + G(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})$$

利用这个表达式，我们可以直接推导出应力-应变关系：
$$\sigma_{ij} = D_{ijkl}\varepsilon_{kl} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$$

而不需要展开 $6\times6$ 的矩阵乘法。

### 2.10.2 FEM 编程中的索引映射

编程实现时，张量的多指标需要映射到一维存储：
- $\sigma_{11} \to \sigma_1,\; \sigma_{22} \to \sigma_2,\; \sigma_{33} \to \sigma_3$
- $\sigma_{12} \to \sigma_4,\; \sigma_{23} \to \sigma_5,\; \sigma_{31} \to \sigma_6$

这就是工程中 Voigt 记号的来源。

### 2.10.3 张量分析工具的价值

| 场景 | 不用张量 | 用张量 |
|------|----------|--------|
| 推导新本构关系 | 需要展开所有分量 | 张量代数一行搞定 |
| 检查坐标不变性 | 需要逐项验证 | 缩并结果自然满足 |
| 代码实现一致性 | 容易遗漏或重复 | 爱因斯坦求和 + 循环 |

---

## 检查你的理解

1. 什么是 Einstein 求和约定？在 $a_{ij}b_{jk} = c_{ik}$ 中，哪个是哑指标，哪个是自由指标？
2. $\delta_{ii}$ 等于多少？$\delta_{ij}\delta_{ij}$ 等于多少？
3. 证明 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 是二阶张量的关键步骤是什么？
4. 为什么 $e_{ijk}$ 对反射变换是赝张量？
5. 弹性力学三类方程的张量形式和矩阵形式有什么对应关系？
6. 什么是张量的缩并？它有什么几何意义？

---

> **对应作业**：[HW1 Q1（恒等式证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q2（张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q3（$e_{ijk}$ 证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
