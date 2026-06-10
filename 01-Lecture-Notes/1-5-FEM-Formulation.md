# 第5章：弹性力学有限元公式推导

> **对应 PDF**：[`5 FEM_formulation.pdf`](../06-References/pdfs-originals/5%20FEM_formulation.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 4 章（Ritz/Galerkin）、线性代数、弹性力学

---

## 5.1 从强形式到弱形式

### 5.1.1 强形式的困难

微分方程的**强形式**要求解在每个点上精确满足方程。以弹性力学为例：
$$\sigma_{ij,j} + f_i = 0 \quad \text{在 } \Omega \text{ 上每一点}$$

这对解的光滑度要求很高（$C^2$ 连续），且对试函数的要求也很高。FEM 不直接处理强形式。

### 5.1.2 弱形式（Weak Form）

通过**加权残量法** + **分部积分**，将 2 阶 PDE 降为 1 阶积分方程。以 Poisson 方程 $-\Delta u = f$ 为例：

取权函数 $\varphi$（$\varphi|_{\partial\Omega}=0$）：
$$-\iint_\Omega \Delta u\cdot\varphi\,d\Omega = \iint_\Omega f\varphi\,d\Omega$$

分部积分（Green 公式）：
$$-\iint_\Omega \nabla\cdot(\nabla u)\varphi = \iint_\Omega \nabla u\cdot\nabla\varphi\,d\Omega - \int_{\partial\Omega}(\nabla u\cdot\mathbf{n})\varphi\,d\Gamma$$

边界项因 $\varphi=0$ 消失：
$$\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,d\Omega = \iint_\Omega f\varphi\,d\Omega$$

这就是弱形式——只要求 $u \in C^0$（一阶导数可积），远低于强形式的 $C^2$。

| | 强形式 | 弱形式 |
|--|--------|--------|
| 方程 | $-\Delta u = f$ | $\int\nabla u\cdot\nabla\varphi = \int f\varphi$ |
| 光滑度 | $C^2$ | $C^0$ |
| 离散方法 | 有限差分 | **有限元** |

### 5.1.3 弹性力学的弱形式（虚位移原理）

从**虚位移原理**（$\delta W_{\text{int}} = \delta W_{\text{ext}}$）直接可得：
$$\int_\Omega \delta\boldsymbol{\varepsilon}^T\boldsymbol{\sigma}\,dV = \int_\Omega \delta\mathbf{u}^T\mathbf{f}\,dV + \int_{S_\sigma} \delta\mathbf{u}^T\bar{\mathbf{T}}\,dS$$

其中 $\delta\mathbf{u}$ 是虚位移，$\delta\boldsymbol{\varepsilon} = [B]\delta\mathbf{u}$ 是虚应变。这个方程是 FEM 离散化的起点。

---

## 5.2 3 节点三角形单元（CST）

### 5.2.1 单元描述

三角形单元 $e = \triangle P_i P_j P_m$，节点按**逆时针**编号。每个节点 2 个 DOF：$(u,v)$。

单元位移向量（6 个分量）：
$$\{\delta\}_e = \begin{pmatrix} u_i & v_i & u_j & v_j & u_m & v_m \end{pmatrix}^T$$

### 5.2.2 形函数

单元内线性插值：
$$u(x,y) = a_1 + a_2x + a_3y,\quad v(x,y) = a_4 + a_5x + a_6y$$

由 3 个节点的 $u$ 值：
$$\begin{cases}
u_i = a_1 + a_2x_i + a_3y_i \\
u_j = a_1 + a_2x_j + a_3y_j \\
u_m = a_1 + a_2x_m + a_3y_m
\end{cases}$$

解此方程组。利用 Cramer 法则：
$$u(x,y) = N_i u_i + N_j u_j + N_m u_m$$

其中：
$$N_i = \frac{1}{2\Delta_e}(a_i + b_i x + c_i y)$$

系数（循环置换 $i\to j\to m$）：
$$
\begin{aligned}
a_i &= x_j y_m - x_m y_j, & b_i &= y_j - y_m, & c_i &= x_m - x_j \\
a_j &= x_m y_i - x_i y_m, & b_j &= y_m - y_i, & c_j &= x_i - x_m \\
a_m &= x_i y_j - x_j y_i, & b_m &= y_i - y_j, & c_m &= x_j - x_i
\end{aligned}
$$

**单元面积**：
$$\Delta_e = \frac12\begin{vmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{vmatrix}$$

### 5.2.3 $[B]$ 矩阵（应变-位移）

应变场：
$$\boldsymbol{\varepsilon} = \begin{pmatrix}
\frac{\partial u}{\partial x} & \frac{\partial v}{\partial y} & \frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}
\end{pmatrix}^T = [B]\{\delta\}_e$$

$$[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}$$

$[B]$ 是常数矩阵 → 单元内应力和应变为常数。这就是**常应变三角形**（CST）名称的由来。

### 5.2.4 $[D]$ 矩阵（弹性矩阵）

**平面应力**：
$$[D] = \frac{E}{1-\nu^2}\begin{pmatrix}1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2}\end{pmatrix}$$

**平面应变**：
$$[D] = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}1 & \frac{\nu}{1-\nu} & 0 \\ \frac{\nu}{1-\nu} & 1 & 0 \\ 0 & 0 & \frac{1-2\nu}{2(1-\nu)}\end{pmatrix}$$

### 5.2.5 单元刚度矩阵

由最小势能原理：
$$[k]_e = \iint_e [B]^T[D][B]\,t\,dxdy = t\Delta_e [B]^T[D][B]$$

$[k]_e$ 是 $6\times6$ 对称矩阵。每个 $2\times2$ 子块 $[k_{rs}]$（$r,s=i,j,m$）对应两个节点间的耦合刚度。

**$K_{ij}$ 的物理意义**：当第 $j$ 个 DOF 为单位位移（其余为零）时，在第 $i$ 个 DOF 所需施加的节点力。

### 5.2.6 等效节点力

体力：
$$\{F_b\}_e = \iint_e [N]^T\mathbf{f}\,t\,dxdy$$

面力（在边界 $\Gamma_e$ 上）：
$$\{F_s\}_e = \int_{\Gamma_e} [N]^T\bar{\mathbf{T}}\,t\,d\Gamma$$

---

## 5.3 完整 CST 数值算例

**问题**：$P_1(0,0), P_2(4,0), P_3(2,3)$，$E=200$ GPa，$\nu=0.3$，$t=0.01$ m，平面应力，求 $[k]_e$。

**Step 1**：单元面积
$$\Delta = \frac12\begin{vmatrix}
1 & 0 & 0 \\ 1 & 4 & 0 \\ 1 & 2 & 3
\end{vmatrix} = 6$$

**Step 2**：$b_i, c_i$

| $i$ | $b_i$ | $c_i$ |
|-----|-------|-------|
| 1 | $y_2-y_3=0-3=-3$ | $x_3-x_2=2-4=-2$ |
| 2 | $y_3-y_1=3-0=3$ | $x_1-x_3=0-2=-2$ |
| 3 | $y_1-y_2=0-0=0$ | $x_2-x_1=4-0=4$ |

**Step 3**：$[B] = \frac1{12}\begin{pmatrix}
-3 & 0 & 3 & 0 & 0 & 0 \\
0 & -2 & 0 & -2 & 0 & 4 \\
-2 & -3 & -2 & 3 & 4 & 0
\end{pmatrix}$

**Step 4**：$[D] = \frac{200\times10^9}{0.91}\begin{pmatrix}1 & 0.3 & 0 \\ 0.3 & 1 & 0 \\ 0 & 0 & 0.35\end{pmatrix}$

**Step 5**：$[k]_e = 0.01 \times 6 \times [B]^T[D][B] = 0.06[B]^T[D][B]$

$[k]_e$ 展开后为 $6\times6$ 对称矩阵。

---

## 5.4 总体集成

### 5.4.1 直接刚度法

将各单元刚度矩阵按节点编号叠加：
$$[K] = \sum_{e=1}^{NE} [k]_e,\quad \{F\} = \sum_{e=1}^{NE} \{F\}_e$$

**局部→全局映射**：单元节点 $i,j,m$ 对应全局节点编号 $I,J,M$。$[k]_e$ 中的子块 $[k_{rs}]$ 叠加到 $[K]$ 的 $(R,S)$ 位置。

### 5.4.2 总刚性质

1. **对称性**：$[K]^T = [K]$
2. **稀疏性**：每个节点只与相邻节点耦合
3. **非负定性**：$\{\delta\}^T[K]\{\delta\} \geq 0$
4. **奇异性**（BC 处理前）：包含刚体位移模式

### 5.4.3 边界条件处理

**划行划列法**（精确）：删去已知位移对应的行和列。

**乘大数法**（编程方便）：令 $K_{ii} \leftarrow N\cdot K_{ii}$（$N$ 为大数如 $10^{15}$），$F_i \leftarrow N\cdot K_{ii}\cdot\bar{u}_i$。

---

## 5.5 求解后处理

1. **单元应变**：$\boldsymbol{\varepsilon}^{(e)} = [B]\{\delta\}_e$
2. **单元应力**：$\boldsymbol{\sigma}^{(e)} = [D]\boldsymbol{\varepsilon}^{(e)}$

CST 单元内应力为常数。后处理通常取节点周围单元应力的平均值作为节点应力。

---

## 5.6 1D 杆单元

作为最简单的 FEM 演示：

**形函数**（$\xi\in[-1,1]$）：$N_1 = \frac{1-\xi}{2},\; N_2 = \frac{1+\xi}{2}$

**应变**：$\varepsilon = \frac{du}{dx} = \left[-\frac1l\;\frac1l\right]\{u\}$

**刚度矩阵**：$[k]_e = \frac{EA}{l}\begin{pmatrix}1 & -1 \\ -1 & 1\end{pmatrix}$

**组装示例**（三杆串联）：
$$\begin{pmatrix}
k_1 & -k_1 & 0 & 0 \\
-k_1 & k_1+k_2 & -k_2 & 0 \\
0 & -k_2 & k_2+k_3 & -k_3 \\
0 & 0 & -k_3 & k_3
\end{pmatrix}
\begin{pmatrix}u_1\\u_2\\u_3\\u_4\end{pmatrix}
= \begin{pmatrix}F_1\\F_2\\F_3\\F_4\end{pmatrix}$$

---

## 检查你的理解

1. 强形式和弱形式的本质区别是什么？弱形式对试函数的光滑度要求为什么更低？
2. CST 单元中 $[B]$ 是常数矩阵，这对精度有何影响？
3. 总体刚度矩阵为什么在施加边界条件前是奇异的？
4. 总装过程中，单元局部节点编号与全局节点编号如何映射？
5. $k_{11}$ 的物理意义是什么？它为什么必须为正？

---

> **对应作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
