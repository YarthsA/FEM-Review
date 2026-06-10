# 第5章：弹性力学有限元公式推导

## 1. 从微分方程到等效积分弱形式

### 1.1 加权残量法的等效积分形式

对微分方程 $L(u) - f = 0$，加权残量法给出：
$$\int_\Omega \omega_i [L(u) - f]\,d\Omega = 0$$

### 1.2 弱形式 (Weak Form)

通过分部积分降低对试探函数的光滑度要求。以 Poisson 方程 $-\Delta u = f$ 为例：

**强形式** (Strong form)：
$$-\Delta u = f \text{ in } \Omega,\quad u = u_0 \text{ on } \partial\Omega$$
要求：$u \in C^2(\Omega)$（二阶导数连续）

**弱形式** (Weak form)：
$$\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,dxdy = \iint_\Omega f\varphi\,dxdy$$
只要求：$u \in C^0(\Omega)$（函数值连续，一阶导数可积）

> **强形式 vs 弱形式**：两者等价（对充分光滑的解）。但弱形式的连续性要求更低，更适合有限元的分片多项式逼近。

### 1.3 Galerkin 弱形式推导（Poisson 方程）

控制方程：$-\Delta u = f,\quad u|_{\partial\Omega} = u_0$

Galerkin 法：取权函数 $\varphi$（满足 $\varphi|_{\partial\Omega}=0$），乘方程积分：
$$-\iint_\Omega \Delta u \cdot \varphi\,dxdy = \iint_\Omega f\varphi\,dxdy$$

分部积分（Green 公式）：
$$-\iint_\Omega \nabla\cdot(\nabla u)\varphi\,dxdy = \iint_\Omega \nabla u \cdot \nabla\varphi\,dxdy - \int_{\partial\Omega} (\nabla u\cdot\mathbf{n})\varphi\,d\Gamma$$

边界项因 $\varphi|_{\partial\Omega}=0$ 消失，得：
$$\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,dxdy = \iint_\Omega f\varphi\,dxdy$$

这就是 Poisson 方程的 Galerkin 弱形式。

---

## 2. 3节点三角形单元 (CST) — 详细推导

### 2.1 单元划分

三角形单元 $e = \Delta P_i P_j P_m$，节点按**逆时针**编号。

### 2.2 插值多项式（形函数）

对单元 $e$ 内的解 $u(x,y)$，采用线性插值：
$$
u(x,y) = a_1 + a_2x + a_3y
$$

3 个节点正好确定 3 个系数：

$$
\begin{cases}
u_i = a_1 + a_2x_i + a_3y_i \\
u_j = a_1 + a_2x_j + a_3y_j \\
u_m = a_1 + a_2x_m + a_3y_m
\end{cases}
$$

### 2.3 形函数的显式表达

$$
u(x,y) = N_i u_i + N_j u_j + N_m u_m = [N]\{\delta\}_e
$$

其中形函数：
$$
N_i = \frac{1}{2\Delta_e}(a_i + b_i x + c_i y)
$$

系数 $a_i, b_i, c_i$（循环置换 i→j→m）：

$$
a_i = x_j y_m - x_m y_j
$$
$$
b_i = y_j - y_m
$$
$$
c_i = x_m - x_j
$$

单元面积：
$$
\Delta_e = \frac{1}{2}\begin{vmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{vmatrix}
$$

### 2.4 $[B]$ 矩阵（应变-位移矩阵）

对 Poisson 方程，$\nabla u$ 的梯度矩阵：
$$
[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & b_j & b_m \\
c_i & c_j & c_m
\end{pmatrix}
$$

即：
$$
\begin{pmatrix} u_x \\ u_y \end{pmatrix} = [B]\{\delta\}_e
$$

### 2.5 单元刚度矩阵

$$
[k]_e = \iint_e [B]^T[B]\,dxdy = \Delta_e[B]^T[B]
$$

显式：
$$
[k]_e = \begin{pmatrix}
k_{ii}^e & k_{ij}^e & k_{im}^e \\
k_{ji}^e & k_{jj}^e & k_{jm}^e \\
k_{mi}^e & k_{mj}^e & k_{mm}^e
\end{pmatrix}
$$

刚度系数：
$$
k_{st}^e = \Delta_e\left[\frac{\partial N_s}{\partial x}\frac{\partial N_t}{\partial x} + \frac{\partial N_s}{\partial y}\frac{\partial N_t}{\partial y}\right]
$$
$$
= \frac{1}{4\Delta_e}(b_s b_t + c_s c_t)
$$

### 2.6 单元载荷向量

$$
\{F\}_e = \iint_e [N]^T f\,dxdy = \begin{pmatrix} F_i^e \\ F_j^e \\ F_m^e \end{pmatrix}
$$

### 2.7 CST 单元刚度矩阵计算示例（完整手算）

**题目**：3 节点三角形单元，节点坐标 $P_1(0,0), P_2(4,0), P_3(2,3)$，材料 $E=200$ GPa, $\nu=0.3$，厚度 $t=0.01$ m，平面应力状态，求 $[k]_e$。

**Step 1**：计算单元面积 $\Delta_e$

$$\Delta_e = \frac12\begin{vmatrix}
1 & 0 & 0 \\
1 & 4 & 0 \\
1 & 2 & 3
\end{vmatrix} = \frac12(1\cdot4\cdot3 + 0 + 0 - 0 - 0 - 1\cdot0\cdot2) = \frac12(12) = 6$$

**Step 2**：计算 $b_i, c_i$（$i=1,2,3$, 对应 $i,j,m$）

$$\begin{aligned}
b_1 &= y_2 - y_3 = 0 - 3 = -3, \quad &c_1 &= x_3 - x_2 = 2 - 4 = -2 \\
b_2 &= y_3 - y_1 = 3 - 0 = 3, \quad &c_2 &= x_1 - x_3 = 0 - 2 = -2 \\
b_3 &= y_1 - y_2 = 0 - 0 = 0, \quad &c_3 &= x_2 - x_1 = 4 - 0 = 4
\end{aligned}$$

**Step 3**：构造 $[B]$ 矩阵（$3 \times 6$）

$$[B] = \frac1{2\times6}\begin{pmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{pmatrix}
= \frac1{12}\begin{pmatrix}
-3 & 0 & 3 & 0 & 0 & 0 \\
0 & -2 & 0 & -2 & 0 & 4 \\
-2 & -3 & -2 & 3 & 4 & 0
\end{pmatrix}$$

**Step 4**：确定弹性矩阵 $[D]$（平面应力）

首先计算 $\frac{E}{1-\nu^2} = \frac{200}{1-0.09} = 219.78$ GPa

$$[D] = 219.78 \times 10^9 \times \begin{pmatrix}
1 & 0.3 & 0 \\
0.3 & 1 & 0 \\
0 & 0 & 0.35
\end{pmatrix} \text{Pa}$$

**Step 5**：计算 $[k]_e = t\Delta_e [B]^T[D][B]$

$$[k]_e = 0.01 \times 6 \times [B]^T[D][B] = 0.06 [B]^T[D][B]$$

$[B]^T[D][B]$ 展开计算得到一个 $6\times6$ 矩阵，每个元素是 $b_i,c_i$ 和 $D_{ij}$ 的组合。此处略去完整乘法，考试中写到 Step 5 的矩阵形式即可。

**关键检查**：
- $[k]_e$ 应对称 ✅
- 对角线元素 $k_{ii} > 0$ ✅
- 单位：N/m ✅

---

## 3. 总体集成 (Assembly)

### 总体刚度矩阵
将各单元刚度矩阵按节点编号叠加：
$$
[K] = \sum_{n=1}^{NE} [k]_{e_n}
$$

### 总体载荷向量
$$
\{F\} = \sum_{n=1}^{NE} \{F\}_{e_n}
$$

### 总体方程
$$
[K]\{\delta\} = \{F\}
$$

### 总体刚度矩阵的性质
1. **对称性**：$[K]^T = [K]$
2. **非负定性**：$\{\delta\}^T[K]\{\delta\} \geq 0$
3. **带状稀疏性**：因形函数为低阶分片多项式，大多数 $K_{ij} = 0$

### 边界条件处理
对已知位移节点 $\delta_i = \bar{u}_i$：
- 划行划列法
- 乘大数法（罚函数法）

---

## 4. 弹性力学平面问题 FEM 公式

### 4.1 矩阵表达式

**位移**：$\mathbf{u} = \begin{pmatrix} u \\ v \end{pmatrix}$

**应变**：$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy} \end{pmatrix}$

**应力**：$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x \\ \sigma_y \\ \tau_{xy} \end{pmatrix}$

### 4.2 弹性矩阵 D

**平面应力**：
$$
\mathbf{D} = \frac{E}{1-\nu^2}\begin{pmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2}
\end{pmatrix}
$$

**平面应变**：
$$
\mathbf{D} = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1 & \frac{\nu}{1-\nu} & 0 \\
\frac{\nu}{1-\nu} & 1 & 0 \\
0 & 0 & \frac{1-2\nu}{2(1-\nu)}
\end{pmatrix}
$$

### 4.3 单元位移场（3节点三角形）

$$
\mathbf{u} = \begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix}
N_i & 0 & N_j & 0 & N_m & 0 \\
0 & N_i & 0 & N_j & 0 & N_m
\end{pmatrix}
\begin{pmatrix} u_i & v_i & u_j & v_j & u_m & v_m \end{pmatrix}^T
$$

$$
= [N]\{\delta\}_e
$$

### 4.4 单元应变场

$$
\boldsymbol{\varepsilon} = \begin{pmatrix} \frac{\partial u}{\partial x} \\ \frac{\partial v}{\partial y} \\ \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \end{pmatrix}
= [B]\{\delta\}_e
$$

$$
[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}
$$

> **CST 单元特点**：$[B]$ 是常数矩阵 → 单元内应力和应变为常数（Constant Strain Triangle）

### 4.5 单元势能

$$
\Pi_e = \frac{1}{2}\{\delta\}_e^T[k]_e\{\delta\}_e - \{\delta\}_e^T\{F\}_e
$$

### 4.6 单元刚度矩阵

$$
[k]_e = \iint_e [B]^T[D][B]\,t\,dxdy = t\Delta_e[B]^T[D][B]
$$

其中 $t$ 为厚度（平面问题）。

### 4.7 等效节点力

体力等效：
$$
\{F_b\}_e = \iint_e [N]^T\mathbf{f}\,t\,dxdy
$$

面力等效（边界 $\Gamma_e$）：
$$
\{F_s\}_e = \int_{\Gamma_e} [N]^T\mathbf{T}\,t\,d\Gamma
$$

---

## 5. 求解过程总结

```
1. 划分网格（节点编号、单元编号、节点坐标）
2. 对每个单元 e:
   a. 计算形函数系数 b_i, c_i
   b. 计算 [B] 矩阵
   c. 计算 [k]_e = tΔ_e[B]^T[D][B]
   d. 计算 {F}_e
   e. 组装到总体 [K] 和 {F}
3. 引入位移边界条件
4. 求解 [K]{δ} = {F} → 得到节点位移
5. 计算单元应变: ε = [B]{δ}_e
6. 计算单元应力: σ = [D]ε
```

---

## 6. 杆单元（1D Bar Element）

### 6.1 基本公式

**单元描述**：2 节点直线杆，每节点 1 个 DOF（轴向位移 $u$）。

**形函数**（自然坐标 $\xi \in [-1,1]$）：
$$N_1 = \frac{1-\xi}{2},\quad N_2 = \frac{1+\xi}{2}$$

**坐标变换**：
$$x = N_1 x_1 + N_2 x_2,\quad \frac{dx}{d\xi} = \frac{l}{2}$$

**应变-位移矩阵**：
$$\varepsilon = \frac{du}{dx} = \frac{d}{dx}[N_1\; N_2]\{u\} = \left[-\frac1l\;\frac1l\right]\{u\}$$

**单元刚度矩阵**：
$$[k]_e = \int_0^l [B]^T EA [B]\,dx = \int_{-1}^1 [B]^T EA [B]\frac{l}{2}d\xi = \frac{EA}{l}\begin{pmatrix}1 & -1 \\ -1 & 1\end{pmatrix}$$

### 6.2 杆系组装示例

**问题**：三杆串联，受节点力 $F_1, F_2, F_3, F_4$，求位移。

每个杆的刚度 $k_i = \frac{EA_i}{l_i}$，单元方程：
$$\begin{pmatrix}k_i & -k_i \\ -k_i & k_i\end{pmatrix}\begin{pmatrix}u_i \\ u_{i+1}\end{pmatrix} = \begin{pmatrix}P_{i1} \\ P_{i2}\end{pmatrix}$$

**总体组装**（直接刚度法）：
$$\begin{pmatrix}
k_1 & -k_1 & 0 & 0 \\
-k_1 & k_1+k_2 & -k_2 & 0 \\
0 & -k_2 & k_2+k_3 & -k_3 \\
0 & 0 & -k_3 & k_3
\end{pmatrix}
\begin{pmatrix}u_1\\u_2\\u_3\\u_4\end{pmatrix}
= \begin{pmatrix}F_1\\F_2\\F_3\\F_4\end{pmatrix}$$

引入 $u_1=0$ 后删去第 1 行第 1 列，解得 $u_2,u_3,u_4$。

---

## 7. 平面悬臂梁例题（CST 单元）

### 建模步骤
1. **网格划分**：将梁划分为若干个三角形单元
2. **每个单元描述**：
   - 节点编号（逆时针）
   - 节点坐标 (x, y)
   - 节点位移 DOF
3. **计算单元刚度矩阵**
4. **组装总刚**
5. **施加约束**（固定端位移 = 0）和**载荷**
6. **求解**

### 注意事项
- 单元节点的编号顺序必须是**逆时针**（保证面积 $\Delta_e > 0$）
- 总体刚度矩阵是奇异的 → 必须施加足够的位移约束消除刚体位移
- CST 单元在每个单元内应力恒定 → 需要较密网格才能获得好的应力精度
