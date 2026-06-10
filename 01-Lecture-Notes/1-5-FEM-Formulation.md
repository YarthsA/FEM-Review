# 第5章：弹性力学有限元公式推导

> **对应 PDF**：[`5 FEM_formulation.pdf`](../06-References/pdfs-originals/5%20FEM_formulation.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 4 章（Ritz/Galerkin 法）、线性代数（矩阵乘法、转置）、材料力学（应力-应变关系）

---

## 5.1 从微分方程到等效积分弱形式

### 5.1.1 强形式 vs 弱形式

FEM 不直接处理微分方程的"强形式"，而是处理其"弱形式"。

| | 强形式（Strong Form） | 弱形式（Weak Form） |
|--|---------------------|-------------------|
| 对象 | 微分方程 + 边界条件 | 加权积分方程 |
| 光滑度要求 | $u \in C^2(\Omega)$ | $u \in C^0(\Omega)$（分部积分后降阶） |
| 物理意义 | 逐点平衡 | 加权平均平衡 |
| 数值方法 | 有限差分法 | **有限元法** |

### 5.1.2 弱形式的推导（以 Poisson 方程为例）

控制方程：$-\Delta u = f \text{ in } \Omega$，$u = u_0 \text{ on } \partial\Omega$

取权函数 $\varphi$（满足 $\varphi|_{\partial\Omega}=0$），乘方程后积分：
$$-\iint_\Omega \Delta u \cdot \varphi\,d\Omega = \iint_\Omega f\varphi\,d\Omega$$

利用分部积分（Green 公式）：
$$-\iint_\Omega \nabla\cdot(\nabla u)\varphi\,d\Omega = \iint_\Omega \nabla u \cdot \nabla\varphi\,d\Omega - \int_{\partial\Omega} (\nabla u\cdot\mathbf{n})\varphi\,d\Gamma$$

由于 $\varphi|_{\partial\Omega}=0$，边界项消失，得：
$$\boxed{\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,d\Omega = \iint_\Omega f\varphi\,d\Omega}$$

这就是 Poisson 方程的 **Galerkin 弱形式**。导数的阶数从 2 阶降到了 1 阶——对 $u$ 和 $\varphi$ 都只要求一阶导数连续（$C^0$），大幅降低了对试函数光滑度的要求。

### 5.1.3 弹性力学问题的弱形式

利用虚位移原理（$\delta W_{\text{int}} = \delta W_{\text{ext}}$），直接可得弹性力学 FEM 的弱形式：
$$\int_\Omega \delta\boldsymbol{\varepsilon}^T\boldsymbol{\sigma}\,dV = \int_\Omega \delta\mathbf{u}^T\mathbf{f}\,dV + \int_{S_\sigma} \delta\mathbf{u}^T\bar{\mathbf{T}}\,dS$$

其中 $\delta\mathbf{u}$ 是虚位移，$\delta\boldsymbol{\varepsilon} = [B]\delta\mathbf{u}$ 是虚应变。

---

## 5.2 3 节点三角形单元（CST）

这是最简单也是最基本的有限元单元——**常应变三角形**（Constant Strain Triangle, CST）。

### 5.2.1 单元划分

三角形单元 $e = \triangle P_i P_j P_m$，3 个节点按**逆时针**编号（保证面积为正）。

每个节点在弹性力学平面问题中有 2 个自由度 $(u, v)$ → 单元共 6 个自由度：
$$\{\delta\}_e = \begin{pmatrix} u_i & v_i & u_j & v_j & u_m & v_m \end{pmatrix}^T$$

### 5.2.2 形函数

单元内位移场采用线性插值：
$$u(x,y) = a_1 + a_2x + a_3y$$
$$v(x,y) = a_4 + a_5x + a_6y$$

对 $u$，由 3 个节点的 $u$ 值确定 $a_1, a_2, a_3$：
$$\begin{cases}
u_i = a_1 + a_2x_i + a_3y_i \\
u_j = a_1 + a_2x_j + a_3y_j \\
u_m = a_1 + a_2x_m + a_3y_m
\end{cases}$$

写成形函数形式：
$$u(x,y) = N_i u_i + N_j u_j + N_m u_m$$

其中：
$$N_i = \frac{1}{2\Delta_e}(a_i + b_ix + c_iy)$$

系数 $a_i, b_i, c_i$（$i,j,m$ 轮换）：
$$\begin{aligned}
a_i &= x_j y_m - x_m y_j, & b_i &= y_j - y_m, & c_i &= x_m - x_j \\
a_j &= x_m y_i - x_i y_m, & b_j &= y_m - y_i, & c_j &= x_i - x_m \\
a_m &= x_i y_j - x_j y_i, & b_m &= y_i - y_j, & c_m &= x_j - x_i
\end{aligned}$$

**单元面积**：
$$\Delta_e = \frac12\begin{vmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{vmatrix}$$

> **验证**：$N_i(x_i,y_i)=1$，$N_i(x_j,y_j)=N_i(x_m,y_m)=0$ ✅
> $\sum N_i = N_i + N_j + N_m = 1$ ✅（能表示刚体平移）

### 5.2.3 $[B]$ 矩阵（应变-位移矩阵）

应变场：
$$\boldsymbol{\varepsilon} = \begin{pmatrix}
\frac{\partial u}{\partial x} & \frac{\partial v}{\partial y} & \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}
\end{pmatrix}^T = [B]\{\delta\}_e$$

$$[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}$$

> **CST 的关键特征**：$[B]$ 是**常数矩阵**（$b_i, c_i$ 由节点坐标决定，与 $(x,y)$ 无关）→ **单元内应变和应力为常数**。这也是"常应变三角形"名称的来源。代价是需要较密的网格才能捕捉应力梯度。

### 5.2.4 弹性矩阵 $[D]$

**平面应力**（薄板）：
$$[D] = \frac{E}{1-\nu^2}\begin{pmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{pmatrix}$$

**平面应变**（长结构）：
$$[D] = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1 & \frac{\nu}{1-\nu} & 0 \\
\frac{\nu}{1-\nu} & 1 & 0 \\
0 & 0 & \frac{1-2\nu}{2(1-\nu)}
\end{pmatrix}$$

### 5.2.5 单元刚度矩阵

由最小势能原理，单元刚度矩阵为：
$$[k]_e = \iint_e [B]^T[D][B]\,t\,dxdy = t\Delta_e[B]^T[D][B]$$

其中 $t$ 为单元厚度。由于 $[B]$ 和 $[D]$ 在单元内都是常数，积分简化为乘积。

$[k]_e$ 是一个 $6\times6$ 的对称矩阵，每个子块 $[k_{rs}]$（$r,s = i,j,m$）是 $2\times2$ 的矩阵。

**$K_{ij}$ 的物理意义**：当第 $j$ 个节点位移为单位位移（其余节点位移为零）时，需要在第 $i$ 个节点施加的节点力大小。

### 5.2.6 等效节点力

体力等效：
$$\{F_b\}_e = \iint_e [N]^T\mathbf{f}\,t\,dxdy$$

面力等效（在边界 $\Gamma_e$ 上）：
$$\{F_s\}_e = \int_{\Gamma_e} [N]^T\bar{\mathbf{T}}\,t\,d\Gamma$$

---

## 5.3 CST 单元完整数值算例

**问题**：3 节点三角形单元，$P_1(0,0), P_2(4,0), P_3(2,3)$，$E=200$ GPa，$\nu=0.3$，$t=0.01$ m，平面应力，求 $[k]_e$。

### Step 1：单元面积
$$\Delta = \frac12\begin{vmatrix}
1 & 0 & 0 \\
1 & 4 & 0 \\
1 & 2 & 3
\end{vmatrix} = \frac12 \times 12 = 6$$

### Step 2：$b_i, c_i$ 系数

| 节点 | $b_i$ | $c_i$ |
|------|-------|-------|
| 1 ($i$) | $y_2-y_3 = 0-3 = -3$ | $x_3-x_2 = 2-4 = -2$ |
| 2 ($j$) | $y_3-y_1 = 3-0 = 3$ | $x_1-x_3 = 0-2 = -2$ |
| 3 ($m$) | $y_1-y_2 = 0-0 = 0$ | $x_2-x_1 = 4-0 = 4$ |

### Step 3：$[B]$ 矩阵
$$[B] = \frac{1}{12}\begin{pmatrix}
-3 & 0 & 3 & 0 & 0 & 0 \\
0 & -2 & 0 & -2 & 0 & 4 \\
-2 & -3 & -2 & 3 & 4 & 0
\end{pmatrix}$$

### Step 4：$[D]$ 矩阵（平面应力）
$$\frac{E}{1-\nu^2} = \frac{200\times10^9}{1-0.09} = 219.78\times10^9\text{ Pa}$$

$$[D] = 219.78\times10^9\begin{pmatrix}
1 & 0.3 & 0 \\
0.3 & 1 & 0 \\
0 & 0 & 0.35
\end{pmatrix}$$

### Step 5：$[k]_e = t\Delta[B]^T[D][B]$
$$[k]_e = 0.01 \times 6 \times [B]^T[D][B] = 0.06[B]^T[D][B]$$

展开后为 $6\times6$ 对称矩阵。检验：$k_{11} > 0$，$k_{12} = k_{21}$。

---

## 5.4 总体集成（Assembly）

### 5.4.1 直接刚度法

将各单元的刚度矩阵和载荷向量按节点编号叠加到总体矩阵中：

$$[K] = \sum_{e=1}^{NE} [k]_e,\quad \{F\} = \sum_{e=1}^{NE} \{F\}_e$$

**组装示例**（杆单元串联）：
$$\begin{pmatrix}
k_1 & -k_1 & 0 & 0 \\
-k_1 & k_1+k_2 & -k_2 & 0 \\
0 & -k_2 & k_2+k_3 & -k_3 \\
0 & 0 & -k_3 & k_3
\end{pmatrix}
\begin{pmatrix}u_1\\u_2\\u_3\\u_4\end{pmatrix}
= \begin{pmatrix}F_1\\F_2\\F_3\\F_4\end{pmatrix}$$

可以看到：总体刚度矩阵中 $K_{22} = k_1 + k_2$（在节点 2 上汇合了两个单元的贡献）。

### 5.4.2 总体方程
$$[K]\{\delta\} = \{F\}$$

### 5.4.3 总刚性质
1. **对称性**：$[K]^T = [K]$（由 Betti 互等定理保证）
2. **稀疏性**：每个单元只连接少数节点 → 每个节点只与相邻节点耦合
3. **非负定性**：$\{\delta\}^T[K]\{\delta\} \geq 0$
4. **奇异性**（引入边界条件前）：包含了刚体位移模式

### 5.4.4 边界条件处理方法

**划行划列法**（精确法）：已知 $\delta_i = \bar{u}_i$，删去第 $i$ 行和第 $i$ 列，右端项减去 $K_{ji}\bar{u}_i$。

**乘大数法**（罚函数法，编程实现方便）：将 $K_{ii}$ 乘以一个大数 $N$（如 $10^{15}$），同时 $F_i = K_{ii} \times N \times \bar{u}_i$。

---

## 5.5 求解后处理

求解 $[K]\{\delta\} = \{F\}$ 后：

1. **单元应变**：$\boldsymbol{\varepsilon}^{(e)} = [B]\{\delta\}_e$
2. **单元应力**：$\boldsymbol{\sigma}^{(e)} = [D]\boldsymbol{\varepsilon}^{(e)}$

对于 CST 单元，每个单元内的应力和应变为常数。通常把得到的应力值分配到节点上做平均处理。

---

## 5.6 1D 杆单元

作为最简单的 FEM 单元，杆单元有助于理解 FEM 的完整流程。

### 形函数
$$N_1 = \frac{1-\xi}{2},\quad N_2 = \frac{1+\xi}{2},\quad \xi = \frac{2x}{l}$$

### 应变-位移矩阵
$$\varepsilon = \frac{du}{dx} = \left[-\frac1l\;\frac1l\right]\{u\}$$

### 刚度矩阵
$$[k]_e = \int_0^l [B]^T EA [B]\,dx = \frac{EA}{l}\begin{pmatrix}1 & -1 \\ -1 & 1\end{pmatrix}$$

### 等效节点力（均布荷载 $q$）
$$\{f\}_e = \int_0^l [N]^T q\,dx = \frac{ql}{2}\begin{pmatrix}1 \\ 1\end{pmatrix}$$

---

## 检查你的理解

1. 为什么 FEM 使用弱形式而非强形式？弱形式的推导中利用了哪个数学工具？
2. 对于 CST 单元，为什么 $[B]$ 矩阵是常数？这对应力和应变分布有何影响？
3. 单元刚度的 $k_{ii}$ 和 $k_{ij}$ 物理意义分别是什么？
4. 总体刚度矩阵的稀疏性是由什么决定的？
5. 为什么在求解前必须引入足量的位移边界条件？

---

> **对应作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **往年相关**：[Homework3 (past)](../04-Homework-Solutions/past/HW3/Homework3.md)
