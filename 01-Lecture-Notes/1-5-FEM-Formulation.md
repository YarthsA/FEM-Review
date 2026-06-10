# 第5章：弹性力学有限元公式推导

> **对应 PDF**：[`5 FEM_formulation.pdf`](../06-References/pdfs-originals/5%20FEM_formulation.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 4 章（Ritz/Galerkin）、线性代数（矩阵乘法、转置）、弹性力学（三类方程、最小势能原理）

---

## 5.1 从强形式到弱形式

### 5.1.1 微分方程的强形式（Strong Form）

弹性力学问题的"强形式"是指**在每个点上精确满足**的微分方程和边界条件。以三维弹性力学为例：

$$\begin{cases}
\sigma_{ij,j} + f_i = 0 &\text{在域 }\Omega\text{ 内} \\
\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i}) &\text{在域 }\Omega\text{ 内} \\
\sigma_{ij} = D_{ijkl}\varepsilon_{kl} &\text{在域 }\Omega\text{ 内} \\
u_i = \bar{u}_i &\text{在位移边界 }S_u\text{ 上} \\
\sigma_{ij}n_j = \bar{T}_i &\text{在力边界 }S_\sigma\text{ 上}
\end{cases}$$

这要求解 $u_i$ 在 $\Omega$ 内**二阶连续可导**（$C^2$）。对于有限元使用的分片多项式试函数，这个要求太高了——两个单元交界处导数可能不连续。

### 5.1.2 等效积分形式

加权残量法的基本思想是：不要求逐点满足方程，只要求**残差的加权积分为零**。

设 $u$ 是试探解，残差 $R = L(u) - f$。令权函数 $\omega$，要求：
$$\int_\Omega \omega R\,d\Omega = 0$$

如果这对任意 $\omega$ 都成立，则等价于 $R \equiv 0$，即强形式。如果我们限制 $\omega$ 在有限维空间中，则得到近似解。

### 5.1.3 分部积分与弱形式

以 Poisson 方程 $-\Delta u = f$ 为例。加权残量法给出：
$$-\iint_\Omega \Delta u \cdot \varphi\,d\Omega = \iint_\Omega f\varphi\,d\Omega$$

直接离散这里会有问题：$\Delta u$ 要求 $u$ 有二阶导数。如果我们对 $\varphi$ 也有一阶导数的要求，可以用分部积分降阶。

利用 Green 公式（二维分部积分）：
$$-\iint_\Omega \nabla\cdot(\nabla u)\varphi\,d\Omega = \iint_\Omega \nabla u \cdot \nabla\varphi\,d\Omega - \int_{\partial\Omega} (\nabla u\cdot\mathbf{n})\varphi\,d\Gamma$$

如果 $\varphi$ 在边界上为零（即满足齐次本质边界条件），边界项消失：
$$\boxed{\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,d\Omega = \iint_\Omega f\varphi\,d\Omega}$$

这就是 **Galerkin 弱形式**。

**强形式 vs 弱形式的关键区别**：

| | 强形式 | 弱形式 |
|--|--------|--------|
| 方程类型 | 微分方程 | 积分方程 |
| 对 $u$ 的要求 | $C^2$（二阶导数连续） | $C^0$（一阶导数平方可积） |
| 对 $\varphi$ 的要求 | 无 | $C^0$（一阶导数平方可积） |
| 物理意义 | 逐点平衡 | 加权平均平衡 |
| 离散方法 | 有限差分 | **有限元** |

### 5.1.4 弹性力学弱形式（虚位移原理）

弹性力学 FEM 的出发点是**虚位移原理**：
$$\int_\Omega \delta\boldsymbol{\varepsilon}^T\boldsymbol{\sigma}\,dV = \int_\Omega \delta\mathbf{u}^T\mathbf{f}\,dV + \int_{S_\sigma} \delta\mathbf{u}^T\bar{\mathbf{T}}\,dS$$

其中：
- $\delta\mathbf{u}$：**虚位移**（满足 $\delta\mathbf{u}|_{S_u}=0$ 的任意微小位移）
- $\delta\boldsymbol{\varepsilon} = [B]\delta\mathbf{u}$：**虚应变**
- 左边：内力的虚功（应变能增量）
- 右边：外力的虚功

这个方程是 FEM 离散化的**起点**。它等价于平衡方程，但只要求 $C^0$ 连续。

---

## 5.2 3 节点三角形单元（CST）的完整推导

### 5.2.1 为什么选择三角形？

三角形单元是 FEM 中最基本的二维单元：
1. **网格划分简单**：任意形状的平面域都可以用三角形填充
2. **适应性强**：可以局部加密，过渡自然
3. **公式简单**：线性插值使 $[B]$ 矩阵为常数

### 5.2.2 单元几何与自由度

三角形单元 $e = \triangle P_i P_j P_m$，3 个节点按**逆时针**编号（这保证面积 $\Delta_e > 0$）。

每个节点在弹性力学平面问题中有 2 个自由度（$x$ 方向位移 $u$，$y$ 方向位移 $v$），所以单元共 6 个自由度：

$$\{\delta\}_e = \begin{pmatrix} u_i & v_i & u_j & v_j & u_m & v_m \end{pmatrix}^T$$

### 5.2.3 形函数的推导

单元内位移场取**线性插值**：
$$u(x,y) = a_1 + a_2 x + a_3 y$$
$$v(x,y) = a_4 + a_5 x + a_6 y$$

对 $u$ 和 $v$ 独立地由 3 个节点值确定 3 个系数。以 $u$ 为例：
$$\begin{cases}
u_i = a_1 + a_2 x_i + a_3 y_i \\
u_j = a_1 + a_2 x_j + a_3 y_j \\
u_m = a_1 + a_2 x_m + a_3 y_m
\end{cases}$$

写成矩阵形式：
$$\begin{pmatrix} u_i \\ u_j \\ u_m \end{pmatrix} = \begin{pmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{pmatrix} \begin{pmatrix} a_1 \\ a_2 \\ a_3 \end{pmatrix}$$

对系数矩阵求逆（利用 Cramer 法则）：
$$\begin{pmatrix} a_1 \\ a_2 \\ a_3 \end{pmatrix} = \frac{1}{2\Delta_e} \begin{pmatrix}
a_i & a_j & a_m \\
b_i & b_j & b_m \\
c_i & c_j & c_m
\end{pmatrix} \begin{pmatrix} u_i \\ u_j \\ u_m \end{pmatrix}$$

其中：
- $a_i = x_j y_m - x_m y_j,\quad b_i = y_j - y_m,\quad c_i = x_m - x_j$
- $a_j = x_m y_i - x_i y_m,\quad b_j = y_m - y_i,\quad c_j = x_i - x_m$
- $a_m = x_i y_j - x_j y_i,\quad b_m = y_i - y_j,\quad c_m = x_j - x_i$

代入 $u(x,y) = a_1 + a_2 x + a_3 y$ 得：
$$u(x,y) = N_i u_i + N_j u_j + N_m u_m$$

其中形函数：
$$\boxed{N_i = \frac{1}{2\Delta_e}(a_i + b_i x + c_i y),\quad i = i,j,m\text{（轮换）}}$$

**单元面积**：
$$\Delta_e = \frac12\begin{vmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{vmatrix}$$

**检验形函数性质**：
- $N_i(x_i,y_i) = 1$，$N_i(x_j,y_j) = N_i(x_m,y_m) = 0$ ✅（Kronecker 性质）
- $N_i + N_j + N_m = 1$ ✅（刚体平移，完备性要求）

类似地 $v(x,y) = N_i v_i + N_j v_j + N_m v_m$。

### 5.2.4 $[B]$ 矩阵（应变-位移矩阵）

对于二维弹性力学问题，应变向量有 3 个分量：
$$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy} \end{pmatrix} = \begin{pmatrix}
\frac{\partial u}{\partial x} \\
\frac{\partial v}{\partial y} \\
\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}
\end{pmatrix}$$

将 $u = N_i u_i + N_j u_j + N_m u_m$ 和 $v$ 的类似表达式代入：

$$\varepsilon_x = \frac{\partial u}{\partial x} = \frac{\partial N_i}{\partial x}u_i + \frac{\partial N_j}{\partial x}u_j + \frac{\partial N_m}{\partial x}u_m$$

其中 $\frac{\partial N_i}{\partial x} = \frac{b_i}{2\Delta_e}$（常数！），所以：
$$\varepsilon_x = \frac{1}{2\Delta_e}(b_i u_i + b_j u_j + b_m u_m)$$

类似地：
$$\varepsilon_y = \frac{1}{2\Delta_e}(c_i v_i + c_j v_j + c_m v_m)$$
$$\gamma_{xy} = \frac{1}{2\Delta_e}[(c_i u_i + b_i v_i) + (c_j u_j + b_j v_j) + (c_m u_m + b_m v_m)]$$

写成矩阵 $\boldsymbol{\varepsilon} = [B]\{\delta\}_e$：

$$\boxed{[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}}$$

**关键观察**：$[B]$ 中的所有元素都是常数（由节点坐标唯一决定，与 $(x,y)$ 无关）。这意味着**单元内应变处处相等**——这就是"常应变三角形"（Constant Strain Triangle，CST）名称的由来。

### 5.2.5 弹性矩阵 $[D]$

**平面应力问题**（薄板，$\sigma_z=0$）：
$$\boxed{[D] = \frac{E}{1-\nu^2}\begin{pmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{pmatrix}}$$

**平面应变问题**（长结构，$\varepsilon_z=0$）：
$$\boxed{[D] = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1 & \frac{\nu}{1-\nu} & 0 \\
\frac{\nu}{1-\nu} & 1 & 0 \\
0 & 0 & \frac{1-2\nu}{2(1-\nu)}
\end{pmatrix}}$$

两个矩阵的区别记住口诀：**"薄板无应力，厚板无应变"**——平面应力令 $\sigma_z=0$，平面应变令 $\varepsilon_z=0$。

### 5.2.6 单元刚度矩阵

由最小势能原理，单元刚度矩阵为：
$$[k]_e = \iint_e [B]^T[D][B]\,dV = \iint_e [B]^T[D][B]\,t\,dxdy$$

由于 $[B]$ 和 $[D]$ 都是常数矩阵，积分简化为矩阵乘积乘以单元面积和厚度：
$$\boxed{[k]_e = t\Delta_e [B]^T[D][B]}$$

$[k]_e$ 是一个 $6\times6$ 的对称矩阵。也可以写成 $2\times2$ 子块的形式：
$$[k]_e = \begin{pmatrix}
[k_{ii}] & [k_{ij}] & [k_{im}] \\
[k_{ji}] & [k_{jj}] & [k_{jm}] \\
[k_{mi}] & [k_{mj}] & [k_{mm}]
\end{pmatrix}$$

每个子块 $[k_{rs}]$（$r,s = i,j,m$）是 $2\times2$ 矩阵：
$$[k_{rs}] = \frac{tE}{4(1-\nu^2)\Delta_e} \begin{pmatrix}
b_r b_s + \frac{1-\nu}{2}c_r c_s & \nu b_r c_s + \frac{1-\nu}{2}c_r b_s \\
\nu c_r b_s + \frac{1-\nu}{2}b_r c_s & c_r c_s + \frac{1-\nu}{2}b_r b_s
\end{pmatrix}$$

（上面是平面应力的公式，平面应变只需替换 $E$ 和 $\nu$）

### 5.2.7 等效节点力

**体力等效节点力**：
$$\{F_b\}_e = \iint_e [N]^T \mathbf{f}\,t\,dxdy$$

其中 $\mathbf{f} = (f_x, f_y)^T$ 是单位体积体力。若体力为常数（如重力），积分可解析计算。对 CST 单元，等效体力节点力等于总体力按 $1/3$ 分配到 3 个节点上。

**面力等效节点力**：
$$\{F_s\}_e = \int_{\Gamma_e} [N]^T \bar{\mathbf{T}}\,t\,d\Gamma$$

其中 $\bar{\mathbf{T}}$ 是单元边界 $\Gamma_e$ 上的分布面力。

---

## 5.3 CST 完整数值算例

**问题**：如图所示的 3 节点三角形单元，节点坐标 $P_1(0,0)$，$P_2(4,0)$，$P_3(2,3)$（单位：m）。材料 $E=200$ GPa（$2\times10^{11}$ Pa），$\nu=0.3$，厚度 $t=0.01$ m。平面应力状态。求单元刚度矩阵 $[k]_e$。

### Step 1：计算单元面积

$$\Delta_e = \frac12\begin{vmatrix}
1 & 0 & 0 \\
1 & 4 & 0 \\
1 & 2 & 3
\end{vmatrix}
= \frac12\left[1\cdot4\cdot3 + 0\cdot0\cdot1 + 0\cdot1\cdot2 - 0\cdot4\cdot1 - 0\cdot0\cdot1 - 1\cdot0\cdot2\right]$$

$$= \frac12(12 + 0 + 0 - 0 - 0 - 0) = \frac12 \times 12 = 6$$

### Step 2：计算 $b_i, c_i$ 系数

按照 $i,j,m$ 对应节点 $1,2,3$（逆时针）：

| $i$ | $b_i$ | 计算 | $c_i$ | 计算 |
|-----|-------|------|-------|------|
| 1 | $b_1 = y_2 - y_3 = 0 - 3 = \mathbf{-3}$ | | $c_1 = x_3 - x_2 = 2 - 4 = \mathbf{-2}$ |
| 2 | $b_2 = y_3 - y_1 = 3 - 0 = \mathbf{3}$ | | $c_2 = x_1 - x_3 = 0 - 2 = \mathbf{-2}$ |
| 3 | $b_3 = y_1 - y_2 = 0 - 0 = \mathbf{0}$ | | $c_3 = x_2 - x_1 = 4 - 0 = \mathbf{4}$ |

**验算**：$b_1 + b_2 + b_3 = -3 + 3 + 0 = 0$，$c_1 + c_2 + c_3 = -2 + (-2) + 4 = 0$ ✅

### Step 3：构造 $[B]$ 矩阵

$$[B] = \frac{1}{2\times6}\begin{pmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{pmatrix}
= \frac1{12}\begin{pmatrix}
-3 & 0 & 3 & 0 & 0 & 0 \\
0 & -2 & 0 & -2 & 0 & 4 \\
-2 & -3 & -2 & 3 & 4 & 0
\end{pmatrix}$$

### Step 4：确定弹性矩阵 $[D]$

平面应力，$\frac{E}{1-\nu^2} = \frac{200\times10^9}{1-0.09} = 219.78\times10^9$ Pa

$$[D] = 219.78\times10^9\begin{pmatrix}
1 & 0.3 & 0 \\
0.3 & 1 & 0 \\
0 & 0 & 0.35
\end{pmatrix} \text{Pa}$$

### Step 5：计算 $[k]_e = t\Delta_e[B]^T[D][B]$

$$[k]_e = 0.01 \times 6 \times [B]^T[D][B] = 0.06[B]^T[D][B]$$

**完整展开**（部分元素）：

$[B]^T$ 是 $6\times3$，$[D]$ 是 $3\times3$，$[B]$ 是 $3\times6$，结果 $[k]_e$ 是 $6\times6$。

先计算 $[D][B]$：
$$[D][B] = 219.78\times10^9\times\frac1{12}\begin{pmatrix}
1 & 0.3 & 0 \\ 0.3 & 1 & 0 \\ 0 & 0 & 0.35
\end{pmatrix}
\begin{pmatrix}
-3 & 0 & 3 & 0 & 0 & 0 \\
0 & -2 & 0 & -2 & 0 & 4 \\
-2 & -3 & -2 & 3 & 4 & 0
\end{pmatrix}$$

$$= 18.315\times10^9\begin{pmatrix}
-3 & -0.6 & 3 & -0.6 & 0 & 1.2 \\
-0.9 & -2 & 0.9 & -2 & 0 & 4 \\
-0.7 & -1.05 & -0.7 & 1.05 & 1.4 & 0
\end{pmatrix}$$

然后 $[k]_e = 0.06[B]^T([D][B])$。最终结果是一个 $6\times6$ 对称矩阵（单位 N/m）。

**快速检查**：
- $k_{11} > 0$（对角线元素应为正）✅
- $k_{12} = k_{21}$（对称性）✅

---

## 5.4 总体集成（Assembly）

### 5.4.1 直接刚度法

FEM 的总体集成采用**直接刚度法**（Direct Stiffness Method）。过程极其直观：

$$[K] = \sum_{e=1}^{NE} [k]_e,\quad \{F\} = \sum_{e=1}^{NE} \{F\}_e$$

"叠加"是按**全局节点编号**进行的。每个单元计算出的局部刚度矩阵 $[k]_e$（$6\times6$，对应局部节点 $i,j,m$）需要映射到全局节点编号 $I,J,M$ 对应的位置。

**映射示例**：
单元局部节点 $i \to$ 全局节点 5，$j \to$ 全局节点 8，$m \to$ 全局节点 12，则：
- $k_{11}^{(e)}$（对应 $u_i$）→ 叠加到 $K$ 的 $(5\times2-1,\;5\times2-1)$ 位置
- $k_{12}^{(e)}$（对应 $u_i, v_i$）→ 叠加到 $K$ 的 $(5\times2-1,\;5\times2)$ 位置
- $k_{13}^{(e)}$（对应 $u_i, u_j$）→ 叠加到 $K$ 的 $(5\times2-1,\;8\times2-1)$ 位置
- 依此类推...

### 5.4.2 总体方程

$$[K]\{\delta\} = \{F\}$$

其中：
- $[K]$：总体刚度矩阵（$N\times N$，$N$ 为总自由度数）
- $\{\delta\}$：节点位移向量（全部未知量）
- $\{F\}$：节点力向量

### 5.4.3 总刚性质

1. **对称性**：$K_{ij} = K_{ji}$（由 Betti 互等定理保证）
2. **稀疏性**：每个节点只与相邻节点通过单元连接，大多数 $K_{ij} = 0$
3. **带状分布**：如果节点编号优化，非零元素集中在主对角线附近
4. **非负定性**：$\{\delta\}^T[K]\{\delta\} \geq 0$（物理上等于 $2\times$ 总应变能）
5. **奇异性**：引入边界条件前，$[K]$ 奇异（刚体位移模式）

### 5.4.4 边界条件处理

**划行划列法**：
已知 $\delta_i = \bar{u}_i$，直接删去第 $i$ 行和第 $i$ 列，并将右端项减去 $K_{ji}\bar{u}_i$。适用于手算。

**乘大数法（罚函数法）**：
将 $K_{ii}$ 乘以一个大数 $N$（如 $10^{15}$），同时将 $F_i$ 改为 $K_{ii}\times N\times\bar{u}_i$。这样第 $i$ 个方程强制给出 $\delta_i \approx \bar{u}_i$。适用于编程实现，不改变矩阵维数。

---

## 5.5 求解与后处理

### 5.5.1 求解流程

```
输入：节点坐标、单元连接、材料属性、边界条件、载荷
  ↓
对每个单元 e：
  ① 提取节点坐标 (x_i,y_i), (x_j,y_j), (x_m,y_m)
  ② 计算 b_i, c_i 系数
  ③ 组装 [B] 矩阵
  ④ 确定 [D] 矩阵
  ⑤ 计算 [k]_e = t·Δ_e·[B]^T[D][B]
  ⑥ 计算等效节点力 {f}_e
  ⑦ 按局部→全局映射组装到 [K] 和 {F}
  ↓
引入边界条件 → 修正 [K] 和 {F}
  ↓
求解线性方程组 [K]{δ} = {F}
  ↓
计算每个单元的应变 ε = [B]{δ}_e
  ↓
计算每个单元的应力 σ = [D]ε
  ↓
输出：节点位移、单元应力、后处理可视化
```

### 5.5.2 应变与应力计算

一旦求得节点位移 $\{\delta\}$，各单元内的应变和应力为：

$$\boldsymbol{\varepsilon}^{(e)} = [B]\{\delta\}_e$$
$$\boldsymbol{\sigma}^{(e)} = [D]\boldsymbol{\varepsilon}^{(e)}$$

对 CST 单元，每个单元内应力和应变是**常数**。这导致相邻单元间应力有跳跃。后处理中通常取节点周围单元应力的平均值作为该节点的应力。

---

## 5.6 1D 杆单元

### 5.6.1 单元公式

杆单元是最简单的 FEM 单元，有助于理解 FEM 的标准流程。

**形函数**（自然坐标 $\xi = 2x/l \in [-1,1]$）：
$$N_1 = \frac{1-\xi}{2},\quad N_2 = \frac{1+\xi}{2}$$

**应变-位移矩阵**：
$$\varepsilon = \frac{du}{dx} = \left[-\frac1l\;\frac1l\right]\{u\}$$

**单元刚度矩阵**：
$$[k]_e = \frac{EA}{l}\begin{pmatrix}1 & -1 \\ -1 & 1\end{pmatrix}$$

**等效节点力**（均布荷载 $q$）：
$$\{f\}_e = \frac{ql}{2}\begin{pmatrix}1 \\ 1\end{pmatrix}$$

### 5.6.2 杆系组装示例

**问题**：3 根杆串联，节点 $1\to2\to3\to4$，受节点力 $F_1=0, F_2, F_3, F_4$。

三根杆的刚度分别为 $k_1 = \frac{EA_1}{l_1}$，$k_2 = \frac{EA_2}{l_2}$，$k_3 = \frac{EA_3}{l_3}$。

各杆的单元方程：
$$\begin{pmatrix}k_i & -k_i \\ -k_i & k_i\end{pmatrix}\begin{pmatrix}u_i \\ u_{i+1}\end{pmatrix} = \begin{pmatrix}P_{i1} \\ P_{i2}\end{pmatrix}$$

**组装过程**（按 $4\times4$ 总刚叠加）：

杆 1 贡献 → 叠加到 $(1,1),(1,2),(2,1),(2,2)$ 位置
杆 2 贡献 → 叠加到 $(2,2),(2,3),(3,2),(3,3)$ 位置
杆 3 贡献 → 叠加到 $(3,3),(3,4),(4,3),(4,4)$ 位置

结果：
$$\begin{pmatrix}
k_1 & -k_1 & 0 & 0 \\
-k_1 & k_1+k_2 & -k_2 & 0 \\
0 & -k_2 & k_2+k_3 & -k_3 \\
0 & 0 & -k_3 & k_3
\end{pmatrix}
\begin{pmatrix}u_1 \\ u_2 \\ u_3 \\ u_4\end{pmatrix}
= \begin{pmatrix}P_{11} \\ P_{12}+P_{21} \\ P_{22}+P_{31} \\ P_{32}\end{pmatrix}
= \begin{pmatrix}F_1 \\ F_2 \\ F_3 \\ F_4\end{pmatrix}$$

这个例子清晰地展示了 FEM 总装的核心：**共享节点的单元刚度叠加**。

---

## 5.7 悬臂梁建模指导

对平面悬臂梁的 CST 分析：

1. **网格划分**：将梁纵向划分成若干段，再从每段矩形对角线分成 2 个三角形各
2. **节点编号**：优化编号以减小带宽（相邻节点编号之差最小化）
3. **计算步骤**：按 5.5.1 的流程逐步推进
4. **注意事项**：
   - 固定端施加 $u=v=0$
   - 自由端施加等效节点力
   - CST 单元内应力恒定 → 至少需要 2-3 层网格才能捕捉弯曲应力梯度

---

## 检查你的理解

1. 强形式和弱形式的主要区别是什么？为什么 FEM 使用弱形式？
2. 写出 CST 单元的 $[B]$ 矩阵，并解释为什么它是常数矩阵。
3. 平面应力和平面应变的 $[D]$ 矩阵有何不同？各自适用于什么场景？
4. 总体刚度矩阵 $[K]$ 的稀疏性是由什么决定的？怎样优化节点编号可减小带宽？
5. 为什么引入边界条件前 $[K]$ 是奇异的？这对求解有什么影响？
6. $k_{ij}$ 的物理意义是什么？利用物理意义判断刚度矩阵中元素的符号。
7. 单元刚度矩阵 $[k]_e$ 的计算过程中，为什么 $[B]^T[D][B]$ 可以提到积分号外？
8. 后处理中，CST 单元的应力是否直接可用？通常需要做什么处理？

---

> **对应作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
