# 弹性力学有限元公式推导

> **对应课件**：[`5 FEM_formulation.pdf`](../06-References/pdfs-originals/5%20FEM_formulation.pdf) 课程第 3 章
> **章节定位**：Finite Element Method of Elastic Mechanics Problems → I. Summary → II. 2D Poisson equations → III. General format of FEM → IV. Plane problems → V. Further discussions
> **相关作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 1-4 章基础、线性代数、偏微分方程基本概念

---

## I. Summary（概述）

变分原理普遍适用于弹性力学，我们将在此范围内讨论有限元分析。

从变分问题的直接解法可知，如果试函数足够接近真实函数，结果将会相当精确。这可以称为**经典变分法**，它在解决一些工程问题时确实有效。

然而如前所述，一旦边界条件相对复杂，选择试函数将变得非常困难。此外试函数随问题不同而异，不便于计算机编码。**有限单元法的发展正是为了克服这些缺点。**

### FEM 的基本思想

有限单元法首先将连续体划分为有限个单元。在数学意义上，它将求解域划分为子域，并对每个子域采用统一的位移函数。

这种位移函数等价于经典变分法中的试函数，但它是针对**单元**设定的——单元之间可能不满足变形协调条件，这使得建立位移函数变得**容易得多**。

以三角形平面单元为例：它可以在单元内部满足变形协调条件，在边界上满足位移边界条件，但不在单元交界处满足变形协调性——因为位移连续但偏导数不连续。

因此问题将关联到**每个单元**，而不是整个连续的全局区域。一旦位移函数满足特定条件，可以证明通过 FEM 得到的结果将收敛到连续体的真实位移。

---

## II. Two-dimensional Poisson equations（二维 Poisson 方程）

### 2.1 The weak form of equivalent integral（等效积分的弱形式）

我们之前已用经典变分法（Galerkin 法）求解过二维 Poisson 方程。

现在来看平面有界域 $\Omega$ 中 Poisson 方程的一般形式：
$$\begin{cases}
-\Delta u = f \\
u|_{\partial\Omega} = u_0(x,y)
\end{cases}$$

其中 $\Delta = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ 是 Laplace 算子。

这是最简单也是最基本的"边值问题"——椭圆型方程，在物理或力学中很常见。例如面外外力作用下的薄膜平衡问题和稳态温度场问题。

无外力作用时 Poisson 方程变为 Laplace 方程。

我们可以直接用 Galerkin 法获得其**等效积分的弱形式**。如前所述，此情况下对试函数的连续性要求很高，分部积分后可以转化为相应的弱形式。

Poisson 方程的弱形式为：
$$\iint_\Omega (u_x\varphi_x + u_y\varphi_y)dxdy = \iint_\Omega f\varphi\,dxdy$$

### 2.2 Element Analysis（单元分析）

#### Element division（单元划分）

将区域 $\Omega$ 划分为若干个单元（例如三角形单元）。采用 **3 节点三角形单元**，单元节点按**逆时针**编号，以保证面积 $\Delta_e > 0$。

记单元 $e = \triangle P_iP_jP_m$，节点为 $P_s(s=i,j,m)$，函数值 $u_s$。

#### Interpolation polynomial（插值多项式）

在每个三角形单元内，假设函数 $u(x,y)$ 为线性插值：
$$u(x,y) = a_1 + a_2x + a_3y$$

由 3 个节点的函数值确定 3 个系数：
$$\begin{cases}
u_i = a_1 + a_2x_i + a_3y_i \\
u_j = a_1 + a_2x_j + a_3y_j \\
u_m = a_1 + a_2x_m + a_3y_m
\end{cases}$$

解出 $a_1,a_2,a_3$ 后，可写为形函数形式：
$$u(x,y) = N_i u_i + N_j u_j + N_m u_m = [N]\{\delta\}_e$$

形函数：
$$N_i = \frac{1}{2\Delta_e}(a_i + b_i x + c_i y),\quad i = i,j,m$$

其中系数：
$$a_i = x_jy_m - x_my_j,\quad b_i = y_j - y_m,\quad c_i = x_m - x_j$$

单元面积：
$$\Delta_e = \frac12\begin{vmatrix}
1 & x_i & y_i \\
1 & x_j & y_j \\
1 & x_m & y_m
\end{vmatrix}$$

#### 单元刚度矩阵与单元荷载向量

弱形式对每个单元求和：
$$\sum_{n=1}^{NE}\iint_{e_n}(u_x\varphi_x+u_y\varphi_y)dxdy = \sum_{n=1}^{NE}\iint_{e_n}f\varphi\,dxdy$$

代入 $u = [N]\{\delta\}_e$，$\nabla u = [B]\{\delta\}_e$，可得：
$$\iint_e (u_x\varphi_x+u_y\varphi_y)dxdy = \iint_e([B]\{\delta^*\}_e)^T([B]\{\delta\}_e)dxdy = \{\delta^*\}_e^T[k]_e\{\delta\}_e$$

其中单元刚度矩阵 $[k]_e = \iint_e [B]^T[B]dxdy = \Delta_e[B]^T[B]$

$$[k]_e = \begin{pmatrix}
k_{ii}^e & k_{ij}^e & k_{im}^e \\
k_{ji}^e & k_{jj}^e & k_{jm}^e \\
k_{mi}^e & k_{mj}^e & k_{mm}^e
\end{pmatrix}$$

刚度系数 $k_{st}^e = \Delta_e\left[\frac{\partial N_s}{\partial x}\frac{\partial N_t}{\partial x} + \frac{\partial N_s}{\partial y}\frac{\partial N_t}{\partial y}\right] = \frac{1}{4\Delta_e}(b_sb_t + c_sc_t)$

单元荷载向量 $\{F\}_e = \iint_e[N]^T f\,dxdy$

### 2.3 Global integration（总体集成）

#### 总体刚度矩阵与总体荷载向量

将各单元刚度矩阵和荷载向量按节点编号展开后叠加：
$$[K] = \sum_{n=1}^{NE}[k]_{e_n},\quad \{F\} = \sum_{n=1}^{NE}\{F\}_{e_n}$$

得到整体方程：
$$\{\delta^*\}^T([K]\{\delta\} - \{F\}) = 0 \quad\Rightarrow\quad [K]\{\delta\} = \{F\}$$

$[K]$ 是**对称**且**非负定**的。因为：
$$\{\delta\}^T[K]\{\delta\} = \sum_{n=1}^{NE}\iint_{e_n}|\nabla u|^2dxdy \geq 0$$

刚度矩阵通常是**带状稀疏**的——因为有限元基函数由低阶分片多项式函数构成，这为离散化和数值求解带来了很大便利。

从插值多项式可得近似解：
$$u(x,y) = [N]\{\delta\}_{e_n} = u_iN_i + u_jN_j + u_mN_m$$

---

## III. General format of finite element method（FEM 的一般形式）

### 3.1 Matrix expression of the basic equations（基本方程的矩阵表达）

现在讨论有限元法在弹性力学中的应用。FEM 起源于结构矩阵方法，但其真正的吸引力在于它成功地解决了**连续体（场）问题**。

弹性力学的常规 FEM 基于**虚位移原理**。所有类型问题的处理流程是一样的——只需将控制方程、坐标、位移和应变替换为相应的表达式即可。

采用最简单的线性四面体单元来学习三维空间弹性问题的有限元推导。

#### 基本变量的矩阵表达

三维弹性问题的位移、应变和应力一般用矩阵形式表达：

**位移**：$\mathbf{u} = \begin{pmatrix} u & v & w \end{pmatrix}^T$

**应变**：$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x & \varepsilon_y & \varepsilon_z & \gamma_{xy} & \gamma_{yz} & \gamma_{zx} \end{pmatrix}^T$

**应力**：$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x & \sigma_y & \sigma_z & \tau_{xy} & \tau_{yz} & \tau_{zx} \end{pmatrix}^T$

#### 微分算子矩阵

$$[\partial] = \begin{pmatrix}
\frac{\partial}{\partial x} & 0 & 0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial x} & \frac{\partial}{\partial z} & 0 \\
0 & 0 & \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
\end{pmatrix}^T$$

#### 三类基本方程

几何方程：$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$

平衡方程：$[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$

本构方程：$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$

其中弹性矩阵 $\mathbf{D}$ 用 Lame 常数表示为：
$$\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}$$

总势能泛函：
$$\Pi = \int_\Omega \frac12\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\mathbf{T}\,dS$$

### 3.2 Discretization and the establishment of the element（离散化与单元建立）

假设求解域是一个多面体，可以用多面体单元完全表示。四面体单元是三维单元中最简单的。

离散化后得到单元 $e_n$ 和节点 $P_i$。设每个节点 $P_i$ 的位移为：
$$\mathbf{u}_i = \begin{pmatrix} u_i & v_i & w_i \end{pmatrix}^T$$

单元 e 的四个顶点为 $P_i, P_j, P_m, P_l$，按**右手规则**排列。

根据四个顶点的位移值可以在单元内确定线性插值：
$$\mathbf{u} = \boldsymbol{\alpha}_1 x + \boldsymbol{\alpha}_2 y + \boldsymbol{\alpha}_3 z + \boldsymbol{\alpha}_4$$

求解 $\boldsymbol{\alpha}$ 并代入，得形函数形式：
$$\mathbf{u} = N_i\mathbf{u}_i + N_j\mathbf{u}_j + N_m\mathbf{u}_m + N_l\mathbf{u}_l = [N]\{\delta\}_e$$

其中 $N_i = \frac{1}{6V_e}\begin{vmatrix}1 & x & y & z \\ 1 & x_j & y_j & z_j \\ 1 & x_m & y_m & z_m \\ 1 & x_l & y_l & z_l\end{vmatrix}$，$V_e$ 为四面体体积。

### 3.3 Solving process（求解过程）

求解 $[K]\{\delta\} = \{F\}$ 得到节点位移后，可以回代计算每个单元的应变和应力。

---

## IV. Plane problems of elastic mechanics（弹性力学平面问题）

### 4.1 Element analysis（单元分析）

对平面弹性力学问题，采用 3 节点三角形单元。取节点 $P_i, P_j, P_m$，每个节点有 2 个自由度 $(u,v)$。

#### 单元位移场
$$\begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix}
N_i & 0 & N_j & 0 & N_m & 0 \\
0 & N_i & 0 & N_j & 0 & N_m
\end{pmatrix}
\begin{pmatrix} u_i & v_i & u_j & v_j & u_m & v_m \end{pmatrix}^T$$

#### 单元应变场
$$\boldsymbol{\varepsilon} = [B]\{\delta\}_e,\quad
[B] = \frac{1}{2\Delta_e}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}$$

#### 单元应力场
$$\boldsymbol{\sigma} = [D][B]\{\delta\}_e$$

#### 单元刚度矩阵
$$[k]_e = \iint_e [B]^T[D][B]\,t\,dxdy = t\Delta_e[B]^T[D][B]$$

平面应力与平面应变的 $[D]$ 矩阵参见第 1 章。

---

## V. Further discussions（进一步讨论）

有限元方法有效克服了经典变分法的局限——试函数选择困难和程序化不易。通过将求解域划分为单元，在单元级别上建立统一的位移函数，整个过程可以完全**标准化**和**自动化**，只需要根据不同的工程问题修改输入数据。

---

## 检查你的理解

1. 为什么 FEM 可以被看作是经典变分法（Ritz 法）的继承和发展？
2. Poisson 方程的弱形式是如何从强形式经过分部积分得到的？
3. 三角形单元刚度矩阵 $k_{st}^e$ 的公式中 $b_s, b_t, c_s, c_t$ 的含义是什么？
4. 总体刚度矩阵 $[K]$ 的对称性和非负定性是如何证明的？

---

> **对应作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
