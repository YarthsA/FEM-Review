# 弹性力学有限元公式推导

> **对应课件**：[`5 FEM_formulation.pdf`](../06-References/pdfs-originals/5%20FEM_formulation.pdf) 课程第 3 章 · [原文MD](../../md_output/5%20FEM_formulation.md)
> **章节定位**：Finite Element Method of Elastic Mechanics Problems → I. Summary → II. 2D Poisson equations → III. General format of FEM → IV. Plane problems → V. Further discussions
> **相关作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 1-4 章基础、线性代数、偏微分方程基本概念
> **记忆辅助**：本章核心是所有 FEM 推导的中枢——Poisson 方程的弱形式 → CST 单元 → 总刚集成 → 数值算例 → 收敛性与误差分析。完整掌握本章，后续的等参元和高阶单元才有坚实基础。

> **📋 考试范围覆盖**
>
> | 本讲义章节 | 考试大纲考点 |
> |-----------|-------------|
> | §5.1 概述 | [FEA] Concepts of elements, nodes, DOFs |
> | §5.2 Poisson 弱形式 | [FEA] Formulation of 1D/2D FEA (element analysis, assembling, BC, solving) |
> | §5.2.2 单元分析 | [FEA] Shape functions of 1D, triangular elements |
> | §5.3 一般形式 | [FEA] Characteristics of element stiffness matrix; [FEA] Physical representation of each element |
> | §5.4 平面问题 | [FEA] Formulation of 1D/2D FEA |
> | §5.5 进一步讨论 | [FEA] Convergence criteria of FEA |

---

## 5.1 概述（Summary）

**本节目标**：理解 FEM 为什么是经典变分法的继承者，以及它在工程应用中克服了哪些关键困难。

### 5.1.1 从经典变分法到有限元

> 🔗 **跨章连接**：这里直接承接第 3 章变分法和第 4 章 Ritz/Galerkin 法。经典变分法的痛点——"试函数难选、不通用、不便编程"——正是 FEM 的出发点。理解这个"痛点驱动"的逻辑链，整章结构就清晰了。

变分原理普遍适用于弹性力学，我们将在此范围内讨论有限元分析。

从变分问题的直接解法可知，如果试函数足够接近真实函数，结果将会相当精确。这可以称为**经典变分法**，它在解决一些工程问题时确实有效。

然而如前所述，一旦边界条件相对复杂，选择试函数将变得非常困难。此外试函数随问题不同而异，不便于计算机编码。**有限单元法的发展正是为了克服这些缺点。**

### 5.1.2 FEM 的基本思想

> 💡 **理解关键**：FEM 的核心智慧在于"分而治之"——与其找一个全局的复杂试函数，不如把区域切成小块，每块用最简单的函数（线性多项式），然后拼起来。这就像用直线段逼近曲线：每段很简单，拼起来却能逼近任意形状。

有限单元法首先将连续体划分为有限个单元。在数学意义上，它将求解域划分为子域，并对每个子域采用统一的位移函数。

这种位移函数等价于经典变分法中的试函数，但它是针对**单元**设定的——单元之间可能不满足变形协调条件，这使得建立位移函数变得**容易得多**。

以三角形平面单元为例：它可以在单元内部满足变形协调条件，在边界上满足位移边界条件，但不在单元交界处满足变形协调性——因为位移连续但偏导数不连续。

因此问题将关联到**每个单元**，而不是整个连续的全局区域。一旦位移函数满足特定条件，可以证明通过 FEM 得到的结果将收敛到连续体的真实位移。

### 5.1.3 本章路线图

本章将从纯数学问题（Poisson 方程）切入，逐步过渡到弹性力学的工程问题，循以下路径：

1. **II. 二维 Poisson 方程**（数学原型，无物理量纲）→ 弱形式 → Galerkin 离散 → CST 形函数 → 单元刚度 → 总体集成
2. **III. FEM 一般形式**（三维弹性力学的矩阵/算子化描述）→ 四面体单元推导 → 总体势能极值 → 乘大数法
3. **IV. 平面弹性问题**（与 II 呼应，但引入 $\mathbf{D}$ 矩阵和物理量纲）→ CST 完整推导 + **3 个完整数值算例**
4. **V. 进一步讨论** → 收敛准则、分片试验、位移元下限性、精度估计

> 💡 **理解关键**：注意这个编排逻辑——先用 Poisson 方程（纯数学、无物理量纲）讲清楚 FEM 的"骨架"，再用弹性力学（叠加 D 矩阵和物理量纲）填充"血肉"。如果你在学弹性力学部分时感到公式繁多，回到 5.2 节看对应关系——两套公式的结构是完全平行的：`[B]^T[B]`（Poisson）→ `[B]^T[D][B]`（弹性）。

---

## 检查你的理解（5.1）

1. 经典变分法的主要困难是什么？FEM 是如何克服的？
2. "单元之间可能不满足变形协调条件"——为什么这反而是 FEM 的优势？
3. 本章 II 节（Poisson 方程）与 IV 节（弹性力学平面问题）之间是什么关系？为什么按这个顺序编排？

---

## 5.2 二维 Poisson 方程（Two-dimensional Poisson Equations）

**本节目标**：以 Poisson 方程为数学原型，完成 FEM 的完整推导流程——从强形式到弱形式，从连续域到离散单元，从单元分析到总体集成。这是所有 FEM 推导的"模板"，后续的弹性力学问题只是在此模板上叠加物理方程。

### 5.2.1 等效积分的弱形式（The weak form of equivalent integral）

> ⚠️ **重难点**：弱形式是整章推导的基石。核心操作就是**分部积分（Green 公式，此处指高斯散度定理，非本构关系的 Green 公式）+ 齐次边界条件消去边界项**。很多人卡在这里是因为没理解"为什么边界项可以丢掉"——答案是试函数 $\varphi$ 在边界上取 0（这是 Galerkin 法的设定），不是数学上恒为 0。

我们之前已用经典变分法（Galerkin 法）求解过二维 Poisson 方程。

现在来看平面有界域 $\Omega$ 中 Poisson 方程的一般形式：

$$\begin{cases}
-\Delta u = f \\
u|_{\partial\Omega} = u_0(x,y)
\end{cases}$$

其中 $f = f(x,y)$ 是源项（Source term），代表单位面积内的热源强度或体力等已知函数；$u_0(x,y)$ 是 Dirichlet 边界上给定的边界函数值。

其中 $\Delta = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ 是 Laplace 算子。

这是最简单也是最基本的第一类"边值问题"——椭圆型方程，在物理或力学中很常见。例如面外外力作用下的薄膜平衡问题和稳态温度场问题。无外力作用时 Poisson 方程变为 Laplace 方程。

我们可以直接用 Galerkin 法获得其**等效积分的弱形式**。如前所述，此情况下对试函数的连续性要求很高（需要 $u \in C^2(\Omega) \cap C^1(\bar{\Omega})$），分部积分后可以转化为相应的弱形式。

**推导过程（Green 公式，此处指高斯散度定理）**：

设 $u(x,y) \in C^2(\Omega) \cap C^1(\bar{\Omega})$，对任意 $\varphi(x,y) \in C^1(\bar{\Omega})$：

$$\begin{aligned}
\int_\Omega (-\Delta u - f)\varphi\,dxdy &= \iint_\Omega (-\Delta u)\varphi\,dxdy + \iint_\Omega (-f)\varphi\,dxdy \\
&= \iint_\Omega \left[\frac{\partial}{\partial x}(u_x\varphi) + \frac{\partial}{\partial y}(u_y\varphi)\right] dxdy \\ 
&\quad + \iint_\Omega [u_x\varphi_x + u_y\varphi_y] dxdy - \iint_\Omega f\varphi\,dxdy \\
&= -\oint_{\partial\Omega} [u_x\varphi n_x + u_y\varphi n_y] ds + \iint_\Omega [u_x\varphi_x + u_y\varphi_y] dxdy - \iint_\Omega f\varphi\,dxdy \\
&= 0
\end{aligned}$$

由于 $\varphi$ 在边界 $\partial\Omega$ 上可取为零（齐次边界条件），边界积分项消失，得：

$$\boxed{\iint_\Omega (u_x\varphi_x + u_y\varphi_y)\,dxdy = \iint_\Omega f\varphi\,dxdy}$$

这就是 Poisson 方程的**等效积分弱形式**。重要特征是：导数阶次从二阶降为一阶，降低了对试函数的连续性要求。


---

### 5.2.2 单元分析（Element Analysis）

#### 单元划分（Element division）

将区域 $\Omega$ 划分为若干个三角形单元。三角形单元最简单、最灵活，能很好地逼近边界，因此被广泛采用。

> ❌ **易错点**：单元划分最容易犯三个错：(1) 节点取在边的内部而非顶点——这破坏了单元间的位移连续性；(2) 存在大的钝角（>120度）——会导致刚度矩阵条件数变差；(3) 相邻节点编号差过大——导致总刚带宽增大、求解效率骤降。

单元划分的注意事项：
1. **节点必须是相邻单元的顶点**，不能是边界上的内点
2. **避免大的钝角**（影响精度）
3. **网格疏密合理**：$u(x,y)$ 梯度变化激烈处加密，平缓处稀疏
4. **节点编号影响带宽**：相邻节点编号差尽量小

记节点为 $P_i(x_i, y_i)\;(i=1,\ldots,NP)$，单元为 $e_k\;(k=1,\ldots,NE)$。

#### 插值多项式（Interpolation polynomial）

在每个三角形单元内采用线性插值：

$$u(x,y) = a + bx + cy$$

有 3 个待定系数。取单元的 3 个顶点作为插值点。设单元 $e = \triangle P_iP_jP_m$，按**逆时针**顺序排列，以保证面积 $\Delta_e > 0$。

> ❌ **易错点**：节点顺序必须逆时针！如果用了顺时针排列，面积行列式 $\Delta_e$ 为负，后续的 $b_i$、$c_i$ 符号会出错，整个刚度矩阵都会算错。验证方法：代三个坐标到面积行列式，正数才是对的。

由节点条件：

$$\boxed{\begin{cases}
a + bx_i + cy_i = u_i \\
a + bx_j + cy_j = u_j \\
a + bx_m + cy_m = u_m
\end{cases}}$$

**用 Cramer 法则求解系数**（这是理解形函数来源的关键步骤）：

> 💡 **理解关键**：Cramer 法则在这里虽然写出了复杂的行列式形式，但本质上就是解 3 元一次方程组。记不住行列式形式没关系——考试时直接解方程组也完全可以。关键是要理解出的结果 $N_s$ 是 $x, y$ 的线性函数。

$$a = \frac{1}{2\Delta_e}\left[\begin{vmatrix} y_j & 1 \\ y_m & 1 \end{vmatrix} u_i + \begin{vmatrix} y_m & 1 \\ y_i & 1 \end{vmatrix} u_j + \begin{vmatrix} y_i & 1 \\ y_j & 1 \end{vmatrix} u_m\right]$$

$$b = \frac{1}{2\Delta_e}\left[-\begin{vmatrix} x_j & 1 \\ x_m & 1 \end{vmatrix} u_i - \begin{vmatrix} x_m & 1 \\ x_i & 1 \end{vmatrix} u_j - \begin{vmatrix} x_i & 1 \\ x_j & 1 \end{vmatrix} u_m\right]$$

$$c = \frac{1}{2\Delta_e}\left[\begin{vmatrix} x_j & y_j \\ x_m & y_m \end{vmatrix} u_i + \begin{vmatrix} x_m & y_m \\ x_i & y_i \end{vmatrix} u_j + \begin{vmatrix} x_i & y_i \\ x_j & y_j \end{vmatrix} u_m\right]$$

其中单元面积：

$$\boxed{2\Delta_e = \begin{vmatrix} 1 & x_i & y_i \\ 1 & x_j & y_j \\ 1 & x_m & y_m \end{vmatrix}}$$

$P_i, P_j, P_m$ 为逆时针排列，故 $\Delta_e > 0$，即为三角形单元的面积。

代入线性函数的一般形式得单元 e 的插值函数：

$$\boxed{u = N_i(x,y)u_i + N_j(x,y)u_j + N_m(x,y)u_m = [N]\{\delta\}_e}$$

**形函数的显式表达式**：

$$N_i = \frac{1}{2\Delta_e}\left[\begin{vmatrix} y_j & 1 \\ y_m & 1 \end{vmatrix} x - \begin{vmatrix} x_j & 1 \\ x_m & 1 \end{vmatrix} y + \begin{vmatrix} x_j & y_j \\ x_m & y_m \end{vmatrix}\right] \equiv \frac{1}{2\Delta_e}(a_i x + b_i y + c_i)$$

$N_j, N_m$ 的表达式通过**轮换下标**获得。引入记号：

$$\boxed{\begin{cases}
a_i = x_j y_m - x_m y_j \\
b_i = y_j - y_m \\
c_i = x_m - x_j
\end{cases}}\quad\text{（}i,j,m\text{ 轮换）}$$

则 $\boxed{N_s = \dfrac{1}{2\Delta_e}(a_s x + b_s y + c_s)},\;(s = i,j,m)$

#### 形函数的六条性质（必须牢记）


**性质 1**：$N_s(x,y)$ 在单元 e 内是一次多项式（$s = i,j,m$）

**性质 2**：$N_s(x_t, y_t) = \delta_{st}$（$s,t = i,j,m$），即形函数在自身节点为 1，其他节点为 0

**性质 3**（几何意义）：$u = N_i(x,y)$ 在空间 $(x,y,u)$ 中表示一张通过 $(x_i,y_i,1)$、$(x_j,y_j,0)$、$(x_m,y_m,0)$ 的平面；$N_j, N_m$ 同理

**性质 4**（单位分解 / 线性再生性）：由于线性插值的唯一性，线性函数插值等于其自身，因此恒有：
$$\begin{aligned}
1 &= N_i + N_j + N_m \\
x &= x_i N_i + x_j N_j + x_m N_m \\
y &= y_i N_i + y_j N_j + y_m N_m
\end{aligned}$$

**性质 5**（面积坐标变换）：考虑线性变换
$$\lambda_1 = N_i(x,y),\quad \lambda_2 = N_j(x,y)$$
则 $N_m = 1 - \lambda_1 - \lambda_2$。逆变换为：
$$\begin{aligned}
x &= (x_i - x_m)\lambda_1 + (x_j - x_m)\lambda_2 + x_m \\
y &= (y_i - y_m)\lambda_1 + (y_j - y_m)\lambda_2 + y_m
\end{aligned}$$
Jacobi 行列式：
$$\left|\frac{\partial(x,y)}{\partial(\lambda_1,\lambda_2)}\right| = \begin{vmatrix} x_i - x_m & x_j - x_m \\ y_i - y_m & y_j - y_m \end{vmatrix} = 2\Delta_e$$

此变换将任意三角形 e 映射为 $\lambda_1$—$\lambda_2$ 平面上的标准三角形 OAB：$P_i \to A(1,0)$，$P_j \to B(0,1)$，$P_m \to O(0,0)$。这是等参元理论的基础。

> 🔗 **跨章连接**：性质 5（面积坐标）直接通向第 6 章的等参元理论。等参元的本质就是用一个标准单元（母单元）通过坐标变换映射到实际单元（子单元），而面积坐标就是这个变换的一维版本。Jacobi 行列式 $2\Delta_e$ 在等参元中对应 $|J|$——提前理解这个对应关系，学第 6 章会轻松很多。

**性质 6**（边界上的线性变化）：在单元的每一条边上，$N_s$ 是该边弧长参数 $t$ 的一次函数。以边 $\overline{P_i P_j}$ 为例（设 $P_i: t=0$，$P_j: t=l = |\overline{P_i P_j}|$）：
$$N_i|_{\overline{P_i P_j}} = 1 - \frac{t}{l},\quad N_j|_{\overline{P_i P_j}} = \frac{t}{l},\quad N_m|_{\overline{P_i P_j}} = 0$$

其中 $N_m|_{\overline{P_i P_j}} = 0$ 是因为 $N_i + N_j + N_m \equiv 1$。

> 💡 **理解关键**：性质 6 是性质 2 和性质 4 的直接推论——在一条边上，第三个节点的形函数必为 0（因为它在边的两端均为 0 且线性变化），而前两个节点的形函数必在端点取 1 和 0，所以只能是 $(1-t/l)$ 和 $t/l$。这条性质保证了相邻单元在边界上的位移连续（C^0 协调）。

#### 梯度向量与 $[B]$ 矩阵

$u$ 的梯度向量：

$$\nabla u = \begin{pmatrix} \frac{\partial u}{\partial x} \\ \frac{\partial u}{\partial y} \end{pmatrix} = \begin{pmatrix} \frac{\partial N_i}{\partial x} & \frac{\partial N_j}{\partial x} & \frac{\partial N_m}{\partial x} \\ \frac{\partial N_i}{\partial y} & \frac{\partial N_j}{\partial y} & \frac{\partial N_m}{\partial y} \end{pmatrix} \{\delta\}_e = \boxed{\frac{1}{2\Delta_e}\begin{pmatrix} a_i & a_j & a_m \\ b_i & b_j & b_m \end{pmatrix}}\{\delta\}_e \equiv [B]\{\delta\}_e$$

注意：$[B]$ 是 $2 \times 3$ 的**常数矩阵**（对三角形线性元而言）。这是 CST（Constant Strain Triangle）名称的由来。

> 💡 **理解关键**："$[B]$ 是常数矩阵"这个看似简单的结论是整个 CST 单元的命门——它意味着单元内的应变处处相同，应力也处处相同。这既是优点（单元刚度矩阵可以直接写成 $[B]^T[B]\Delta_e$，无需数值积分），也是缺点（在应变梯度大的地方精度差）。当你在算例中看到单元刚度矩阵直接等于 $\Delta_e[B]^T[B]$ 时，原因就在这里。

#### 单元刚度矩阵和单元荷载向量

将弱形式按单元分片表达：

$$\sum_{n=1}^{NE}\iint_{e_n} (u_x\varphi_x + u_y\varphi_y)dxdy = \sum_{n=1}^{NE}\iint_{e_n} f\varphi\,dxdy$$

设单元 $e = \triangle P_iP_jP_m$，记 $\varphi$ 的节点值为 $\{\delta^*\}_e$。

代入 $\nabla u = [B]\{\delta\}_e$：

$$\begin{aligned}
\iint_e (u_x\varphi_x + u_y\varphi_y)dxdy &= \iint_e \{\nabla\varphi\}^T\{\nabla u\}dxdy \\
&= \iint_e ([B]\{\delta^*\}_e)^T([B]\{\delta\}_e)dxdy \\
&= \{\delta^*\}_e^T \iint_e ([B]^T[B])dxdy \;\{\delta\}_e \\
&= \{\delta^*\}_e^T [k]_e \{\delta\}_e
\end{aligned}$$

其中**单元刚度矩阵** $[k]_e$ 为 $3\times 3$ 矩阵：

$$\boxed{[k]_e = \iint_e [B]^T[B]dxdy = \Delta_e [B]^T[B]}$$


$$[k]_e = \begin{pmatrix}
k_{ii}^e & k_{ij}^e & k_{im}^e \\
k_{ji}^e & k_{jj}^e & k_{jm}^e \\
k_{mi}^e & k_{mj}^e & k_{mm}^e
\end{pmatrix}$$

**刚度系数的显式公式**：

$$\boxed{k_{st}^e = \Delta_e\left[\frac{\partial N_s}{\partial x}\frac{\partial N_t}{\partial x} + \frac{\partial N_s}{\partial y}\frac{\partial N_t}{\partial y}\right] = \frac{1}{4\Delta_e}(a_s a_t + b_s b_t),\quad (s,t = i,j,m)}$$

右端项：

$$\iint_e f\varphi\,dxdy = \iint_e ([N]\{\delta^*\}_e)^T f\,dxdy = \{\delta^*\}_e^T \{F\}_e$$

其中**单元荷载向量**：

$$\boxed{\{F\}_e = \iint_e [N]^T f\,dxdy = \begin{pmatrix} F_i^e \\ F_j^e \\ F_m^e \end{pmatrix}}$$

---

### 5.2.3 总体集成（Global integration）

#### 总体刚度矩阵与总体荷载向量

> ⚠️ **重难点**：总体集成（assembly）常被忽略，但它是最容易在编程/作业中出错的地方。核心操作是"对号入座"——单元刚度矩阵的局部下标 (i, j, m) 要映射到全局节点编号，然后叠加到总刚对应位置。理解了这个操作，才能理解为什么总刚是稀疏带状的、为什么非对角项只在相邻节点间非零。

将各单元的 $\{\delta\}_e$、$\{F\}_e$、$[k]_e$ 扩展到 NP 维（NP 为总节点数）以便叠加。方法是将单元刚度系数和荷载分量按全局节点编号"对号入座"写入总刚和总荷的对应位置。

**需要注意的是**：三角形单元三顶点必须按逆时针排列，全局编号可能不是有序的（例如 $j<i<m$）。

叠加后得到：

$$\{\delta^*\}^T\left(\sum_{n=1}^{NE}[k]_{e_n}\right)\{\delta\} = \{\delta^*\}^T\left(\sum_{n=1}^{NE}\{F\}_{e_n}\right)$$

即：

$$\boxed{[K] = \sum_{n=1}^{NE}[k]_{e_n},\qquad \{F\} = \sum_{n=1}^{NE}\{F\}_{e_n}}$$

由于 $\{\delta^*\}$ 是任意 NP 维向量（$\varphi$ 在边界上无约束），得到线性方程组：

$$\boxed{[K]\{\delta\} = \{F\}}$$

#### $[K]$ 的对称性和非负定性

$[K]$ 是对称且非负定的。证明：

> ⚠️ **重难点**：为什么总刚是"非负定"而不是"正定"？非负定意味着 $\{\delta\}^T[K]\{\delta\} \geq 0$，但可能等于 0。等于 0 的情况对应刚体位移——所有节点的 $|\nabla u| = 0$，即 $u$ 为常数。引入位移边界约束后，刚体位移被消除，总刚才变为正定。这个区别分析题目中经常被问到。

$$\begin{aligned}
\{\delta\}^T[K]\{\delta\} &= \sum_{n=1}^{NE} \{\delta\}_{e_n}^T [K]_{e_n} \{\delta\}_{e_n} \\
&= \sum_{n=1}^{NE} \iint_{e_n} ([B]\{\delta\}_{e_n})^T([B]\{\delta\}_{e_n})dxdy \\
&= \sum_{n=1}^{NE} \iint_{e_n} |\nabla u|^2 dxdy \geq 0
\end{aligned}$$

刚度矩阵通常是**带状稀疏**的——因为有限元基函数由低阶分片多项式函数构成。只有当两节点属于同一单元时，对应的刚度子块才非零。这为离散化和数值求解带来了很大便利。

#### 求解和近似解

处理约束条件后求解 $[K]\{\delta\} = \{F\}$ 得到各节点函数值 $u_i\;(i=1,\ldots,NP)$。任一单元内的近似解为：

$$u(x,y) = [N]\{\delta\}_{e_n} = u_i N_i(x,y) + u_j N_j(x,y) + u_m N_m(x,y)$$

---

## 检查你的理解（5.2）

1. 弱形式的本质优势是什么？为什么 FEM 中普遍采用弱形式而非强形式？
2. 形函数六条性质中，哪两条直接保证了单元间的位移连续性（协调性）？
3. 单元刚度矩阵 $k_{st}^e = \frac{1}{4\Delta_e}(a_s a_t + b_s b_t)$ 是如何从 $[B]^T[B]$ 推导出来的？
4. 总刚矩阵 $[K]$ 为什么是"非负定"而非"正定"？什么条件下才变成正定？

---

## 5.3 FEM 的一般形式（General Format of FEM）

**本节目标**：推导三维弹性力学问题中基于最小势能原理的 FEM 一般形式。这是所有位移元 FEM 的统一框架——平面问题和三维问题的区别仅在于 $[B]$、$[D]$ 矩阵的维度和内容。

### 5.3.1 基本方程的矩阵表达（Matrix expression of the governing equations）

现在讨论有限元法在弹性力学中的应用。FEM 起源于结构矩阵方法，但其真正的吸引力在于它成功地解决了**连续体（场）问题**。

> 💡 **理解关键**：5.3 节用矩阵/算子语言重写了弹性力学的三大方程，看起来很抽象，但目的很明确——把所有东西都写成统一形式，使得推导可以机械化。关键对应：几何方程 $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$、本构方程 $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$、平衡方程 $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$。合并后就是 $[\partial]^T\mathbf{D}[\partial]\mathbf{u} + \mathbf{f} = \mathbf{0}$——这就是用位移表达的平衡方程，是位移元的出发点。

弹性力学的常规 FEM 基于**虚位移原理**（等价于最小势能原理）。所有类型问题的处理流程是一样的——只需将控制方程、坐标、位移和应变替换为相应的表达式即可。

#### 基本变量的矩阵形式

三维弹性问题：
- **位移**：$\mathbf{u} = \begin{pmatrix} u & v & w \end{pmatrix}^T$
- **应变**：$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x & \varepsilon_y & \varepsilon_z & \gamma_{xy} & \gamma_{yz} & \gamma_{zx} \end{pmatrix}^T$
- **应力**：$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x & \sigma_y & \sigma_z & \tau_{xy} & \tau_{yz} & \tau_{zx} \end{pmatrix}^T$

#### 微分算子矩阵

$$[\partial] = \begin{pmatrix}
\frac{\partial}{\partial x} & 0 & 0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial x} & \frac{\partial}{\partial z} & 0 \\
0 & 0 & \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
\end{pmatrix}^T \quad (6\times 3)$$

#### 三类基本方程

**几何方程**（$6\times 1 = (6\times 3)(3\times 1)$）：$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$

**平衡方程**（$3\times 1 = (3\times 6)(6\times 1) + (3\times 1)$）：$[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$，其中 $\mathbf{f} = (f_x, f_y, f_z)^T$ 是体力（Body force）向量，表示单位体积内所受的外力。

**本构方程**：$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$

其中弹性矩阵 $\mathbf{D}$ 用 Lamé 常数表示（$6 \times 6$），$\lambda$ 和 $G$ 称为 Lamé 常数，与弹性模量 $E$ 和泊松比 $\nu$ 的关系为 $\lambda = \dfrac{E\nu}{(1+\nu)(1-2\nu)}$、$G = \dfrac{E}{2(1+\nu)}$：

$$\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}$$

代入得平衡方程用位移表示的形式：

$$[\partial]^T\mathbf{D}[\partial]\mathbf{u} + \mathbf{f} = \mathbf{0}$$

约定边界记号：弹性体的边界 $S$ 分为位移边界 $S_u$ 和外力边界 $S_\sigma$（$S = S_u \cup S_\sigma$，$S_u \cap S_\sigma = \emptyset$）。$\bar{\mathbf{u}}$ 和 $\mathbf{T}$ 分别表示 $S_u$ 和 $S_\sigma$ 上给定的已知位移和已知面力。

**位移边界条件**：$\mathbf{u}|_{S_u} = \bar{\mathbf{u}}$

**力边界条件**（Cauchy 公式的矩阵形式）：

$$\begin{pmatrix} \cos(n,x) & 0 & 0 & \cos(n,y) & 0 & \cos(n,z) \\ 0 & \cos(n,y) & 0 & \cos(n,x) & \cos(n,z) & 0 \\ 0 & 0 & \cos(n,z) & 0 & \cos(n,y) & \cos(n,x) \end{pmatrix}\boldsymbol{\sigma} = \mathbf{T}$$

**总势能泛函**的矩阵表达：

$$\Pi = \int_\Omega \frac{1}{2}\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\mathbf{T}\,dS$$

---

### 5.3.2 离散化与单元建立（Discretization and the establishment of the element）

#### 四面体单元

假设求解域是一个多面体，可用四面体单元完全表示。四面体单元是三维单元中最简单的。

离散化后得单元 $e_n\;(n=1,\ldots,NE)$ 和节点 $P_i(x_i, y_i, z_i)\;(i=1,\ldots,NP)$。

每个节点的位移：$\mathbf{u}_i = \begin{pmatrix} u_i & v_i & w_i \end{pmatrix}^T$

单元 $e$ 的四个顶点按**右手规则**排列为 $P_i, P_j, P_m, P_l$——即 $P_i, P_j, P_m$ 右手旋转时，$P_l$ 在拇指方向。

> ❌ **易错点**：四面体单元的节点排列必须满足右手规则，这和三角形单元的逆时针规则类似但更复杂。如果排列不对，体积 $V_e$ 可能为负。而且三维轮换不能简单照搬二维的规则——如果从偶数下标开始轮换，需要修正符号（因为会从右手规则变成左手规则）。

#### 插值位移场

在四面体单元内采用线性插值（每个方向 4 个待定系数，共 12 个）：

$$\mathbf{u} = \boldsymbol{\alpha}_1 x + \boldsymbol{\alpha}_2 y + \boldsymbol{\alpha}_3 z + \boldsymbol{\alpha}_4$$

满足节点条件 $\mathbf{u}_s = \boldsymbol{\alpha}_1 x_s + \boldsymbol{\alpha}_2 y_s + \boldsymbol{\alpha}_3 z_s + \boldsymbol{\alpha}_4\;(s = i,j,m,l)$。用 Cramer 法则求解 $\boldsymbol{\alpha}_i$，得形函数形式：

$$\mathbf{u} = N_i\mathbf{u}_i + N_j\mathbf{u}_j + N_m\mathbf{u}_m + N_l\mathbf{u}_l = [N]\{\delta\}_e$$

其中：

$$\{\delta\}_e = \begin{pmatrix} u_i & v_i & w_i & u_j & v_j & w_j & u_m & v_m & w_m & u_l & v_l & w_l \end{pmatrix}^T$$

$$[N] = \begin{bmatrix} N_i\mathbf{I}_3 & N_j\mathbf{I}_3 & N_m\mathbf{I}_3 & N_l\mathbf{I}_3 \end{bmatrix}_{3\times 12}$$，其中 $\mathbf{I}_3$ 是 $3\times 3$ 单位矩阵（Identity matrix）。

形函数 $N_s$ 通过行列式表达：

$$N_i = \frac{1}{6V_e}\begin{vmatrix} 1 & x & y & z \\ 1 & x_j & y_j & z_j \\ 1 & x_m & y_m & z_m \\ 1 & x_l & y_l & z_l \end{vmatrix},\quad\ldots$$

$$V_e = \frac{1}{6}\begin{vmatrix} 1 & x_i & y_i & z_i \\ 1 & x_j & y_j & z_j \\ 1 & x_m & y_m & z_m \\ 1 & x_l & y_l & z_l \end{vmatrix}$$

$V_e$ 为四面体体积。形函数也可写成标准形式：

$$N_s = \frac{1}{6V_e}(a_s x + b_s y + c_s z + d_s),\quad (s = i,j,m,l)$$

其中 $a_s, b_s, c_s, d_s$ 是与节点几何位置相关的系数。例如对节点 $i$：

$$a_i = -\begin{vmatrix} 1 & y_j & z_j \\ 1 & y_m & z_m \\ 1 & y_l & z_l \end{vmatrix},\; b_i = \begin{vmatrix} 1 & x_j & z_j \\ 1 & x_m & z_m \\ 1 & x_l & z_l \end{vmatrix},\; c_i = -\begin{vmatrix} 1 & x_j & y_j \\ 1 & x_m & y_m \\ 1 & x_l & y_l \end{vmatrix},\; d_i = -\begin{vmatrix} x_j & y_j & z_j \\ x_m & y_m & z_m \\ x_l & y_l & z_l \end{vmatrix}$$

**注意**：三维情形下的轮换与平面问题不同——因为右手规则不能自动保证。如果轮换起始于偶序数（如从 j=2 开始），则启用左手规则，需修正符号。

---

### 5.3.3 求解过程（Solving process）

#### 应变和应力的离散表达

将位移插值代入几何方程：

$$\boldsymbol{\varepsilon} = [\partial]\mathbf{u} = [\partial][N]\{\delta\}_e = [B]\{\delta\}_e$$

$$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon} = \mathbf{D}[B]\{\delta\}_e$$

其中 $[B]$ 为 $(6 \times 12)$ 应变矩阵：

$$[B] = [\partial][N] = \begin{bmatrix} B_i & B_j & B_m & B_l \end{bmatrix}$$

经过分块矩阵运算：

$$B_s = [\partial]N_s\mathbf{I}_3 = \frac{1}{6V_e}\begin{pmatrix}
a_s & 0 & 0 \\
0 & b_s & 0 \\
0 & 0 & c_s \\
b_s & a_s & 0 \\
0 & c_s & b_s \\
c_s & 0 & a_s
\end{pmatrix},\quad (s = i,j,m,l)$$

#### 从势能泛函到 FEM 方程

离散模型的总势能等于各单元势能之和：

> 💡 **理解关键**：这一段的推导流程是 FEM 最核心的"生产线"——(1) 势能泛函 → (2) 单元离散 → (3) 代入 $\mathbf{u}=[N]\{\delta\}_e$、$\boldsymbol{\varepsilon}=[B]\{\delta\}_e$ → (4) 提取 $\{\delta\}_e^T(\cdots)\{\delta\}_e$ 得到 $\mathbf{K}^e$ → (5) 组装总刚 → (6) $\delta\Pi=0$ 得 $\mathbf{K}\mathbf{a}=\mathbf{P}$。所有不同问题（平面应力、平面应变、三维、板壳……）的区别只在 $[B]$ 和 $\mathbf{D}$ 的具体内容，上述步骤完全一样。

$$\begin{aligned}
\Pi &= \sum_e \Pi_p^e = \sum_e \left(\int_{\Omega_e}\frac{1}{2}\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}dV - \int_{\Omega_e}\mathbf{u}^T\mathbf{f}dV - \int_{S_\sigma^e}\mathbf{u}^T\mathbf{T}dS\right) \\
&= \sum_e \left(\{\delta\}_e^T \int_{\Omega_e}\frac{1}{2}\mathbf{B}^T\mathbf{D}\mathbf{B}dV \{\delta\}_e\right) - \sum_e \left(\{\delta\}_e^T \int_{\Omega_e}\mathbf{N}^T\mathbf{f}dV\right) - \sum_e \left(\{\delta\}_e^T \int_{S_\sigma^e}\mathbf{N}^T\mathbf{T}dS\right)
\end{aligned}$$

定义：
- **单元刚度矩阵**：$\boxed{\mathbf{K}_{(12\times 12)}^e = \int_{\Omega_e} \mathbf{B}^T\mathbf{D}\mathbf{B}\,dV}$
- **等效节点荷载列阵**：$\mathbf{P}^e = \int_{\Omega_e} \mathbf{N}^T\mathbf{f}\,dV + \int_{S_\sigma^e} \mathbf{N}^T\mathbf{T}\,dS$

则势能化为：

$$\Pi = \frac{1}{2}\sum_e (\{\delta\}_e^T \mathbf{K}^e \{\delta\}_e) - \sum_e (\{\delta\}_e^T \mathbf{P}^e)$$

组装得总刚和总荷载：

$$\Pi = \frac{1}{2}\mathbf{a}^T\mathbf{K}\mathbf{a} - \mathbf{a}^T\mathbf{P}$$

此处记号从单元节点位移 $\{\delta\}_e$ 切换为总体节点位移 $\mathbf{a}$（即 $\mathbf{a}$ 是全体节点位移组成的总体位移向量）。

由 $\delta\Pi = 0$（势能极值条件），得 FEM 的基本方程：

$$\boxed{\mathbf{K}\mathbf{a} = \mathbf{P}}$$

#### 四面体单元刚度矩阵的分块形式

对于线性四面体单元，$[B]$ 是常数矩阵，所以：

$$\mathbf{K}_{(3\times 3)_{st}}^e = \mathbf{B}_s^T\mathbf{D}\mathbf{B}_t V_e,\quad (s,t = i,j,m,l)$$

显式表达式：

$$\mathbf{K}_{st}^e = \frac{1}{36V_e}\begin{bmatrix}
a_s a_t(\lambda+2G)+(b_s b_t+c_s c_t)G & a_s b_t\lambda + b_s a_t G & a_s c_t\lambda + c_s a_t G \\
b_s a_t\lambda + a_s b_t G & b_s b_t(\lambda+2G)+(a_s a_t+c_s c_t)G & b_s c_t\lambda + c_s b_t G \\
c_s a_t\lambda + a_s c_t G & c_s b_t\lambda + b_s c_t G & c_s c_t(\lambda+2G)+(b_s b_t+a_s a_t)G
\end{bmatrix}$$

总刚度矩阵为 $(12 \times 12)$：

$$\mathbf{K}^e = \begin{bmatrix}
K_{ii}^e & K_{ij}^e & K_{im}^e & K_{il}^e \\
K_{ji}^e & K_{jj}^e & K_{jm}^e & K_{jl}^e \\
K_{mi}^e & K_{mj}^e & K_{mm}^e & K_{ml}^e \\
K_{li}^e & K_{lj}^e & K_{lm}^e & K_{ll}^e
\end{bmatrix}_{12\times 12}$$

#### 边界条件的引入——乘大数法

当节点位移被指定 $a_j = \bar{a}_j$ 时，修改第 j 个方程：

- 将 $K_{jj}$ 乘以一个大数 $\alpha$（如 $\alpha \approx 10^{10}$）
- 将右端项 $P_j$ 替换为 $\alpha K_{jj}\bar{a}_j$

即：
$$\alpha K_{jj} a_j \approx \alpha K_{jj}\bar{a}_j \;\Rightarrow\; a_j = \bar{a}_j$$

**乘大数法的优点**：不改变方程阶数和节点编号顺序，极便于编程实现，因此被广泛采用于 FEM 程序。


---

## 检查你的理解（5.3）

1. 三维弹性力学 FEM 的 $[B]$ 矩阵维度是多少？与二维情况有何本质区别？
2. $\mathbf{K}^e = \int_{\Omega_e} \mathbf{B}^T\mathbf{D}\mathbf{B}\,dV$ 中，如果 $[B]$ 不是常数矩阵（如高阶单元），这个积分如何处理？
3. "乘大数法"处理的物理本质是什么？为什么 $\alpha$ 取 $10^{10}$ 就足够？

---

## 5.4 弹性力学平面问题（Plane Problems of Elastic Mechanics）

**本节目标**：完成 CST（Constant Strain Triangle）单元在弹性力学平面问题中的完整推导，并通过三个完整数值算例展示 FEM 的全流程计算能力。

弹性平面问题包括**平面应力问题**和**平面应变问题**。两者的区别仅在于弹性矩阵 $[D]$ 的取法，推导过程完全相同。

### 5.4.1 平面问题的矩阵形式

$$\text{位移}\;\mathbf{u} \equiv \begin{pmatrix} u \\ v \end{pmatrix},\quad \text{应变}\;\boldsymbol{\varepsilon} \equiv \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \gamma_{xy} \end{pmatrix},\quad \text{应力}\;\boldsymbol{\sigma} \equiv \begin{pmatrix} \sigma_x \\ \sigma_y \\ \tau_{xy} \end{pmatrix}$$

微分算子矩阵：$[\partial] = \begin{pmatrix} \frac{\partial}{\partial x} & 0 \\ 0 & \frac{\partial}{\partial y} \\ \frac{\partial}{\partial y} & \frac{\partial}{\partial x} \end{pmatrix}$

### 5.4.2 单元分析（Element Analysis）

仍以 3 节点三角形单元为例。该单元有 6 个自由度。

> 🔗 **跨章连接**：5.4.2 与 5.2.2 的结构完全平行——都是 (A) 位移场构造 → (B) B 矩阵 → (C) 应力/D 矩阵 → (D) 单元刚度 → (E) 等效荷载 → (F) 单元方程。区别在于弹性力学多了一个 D 矩阵，且每个节点的自由度从 1 个（标量 u）变成了 2 个（向量 u, v）。对比阅读这两个小节，你会清晰看到 FEM 标准流程是如何从数学问题"升级"到物理问题的。

节点位移列阵和节点力列阵：

$$\mathbf{q}_{(6\times 1)}^e = \begin{pmatrix} u_1 & v_1 & u_2 & v_2 & u_3 & v_3 \end{pmatrix}^T$$

$$\mathbf{P}_{(6\times 1)}^e = \begin{pmatrix} P_{1x} & P_{1y} & P_{2x} & P_{2y} & P_{3x} & P_{3y} \end{pmatrix}^T$$

#### (A) 单元位移场的构造

假设每个方向的位移由 3 个节点位移确定，因此每个方向有 3 个待定系数：

$$\begin{cases}
u(x,y) = \beta_1 + \beta_2 x + \beta_3 y \\
v(x,y) = \beta_4 + \beta_5 x + \beta_6 y
\end{cases}$$

这是最简单的满足**完备性**（含刚体位移和常应变模态）和**唯一性**的位移模式。

**由节点条件确定系数**（以 u 方向为例）：

节点条件：$u(x_i, y_i) = u_i\;(i=1,2,3)$

用 Cramer 法则求解：

$$\begin{aligned}
\beta_1 &= \frac{1}{2A}\begin{vmatrix} u_1 & x_1 & y_1 \\ u_2 & x_2 & y_2 \\ u_3 & x_3 & y_3 \end{vmatrix} = \frac{1}{2A}(a_1 u_1 + a_2 u_2 + a_3 u_3) \\
\beta_2 &= \frac{1}{2A}\begin{vmatrix} 1 & u_1 & y_1 \\ 1 & u_2 & y_2 \\ 1 & u_3 & y_3 \end{vmatrix} = \frac{1}{2A}(b_1 u_1 + b_2 u_2 + b_3 u_3) \\
\beta_3 &= \frac{1}{2A}\begin{vmatrix} 1 & x_1 & u_1 \\ 1 & x_2 & u_2 \\ 1 & x_3 & u_3 \end{vmatrix} = \frac{1}{2A}(c_1 u_1 + c_2 u_2 + c_3 u_3)
\end{aligned}$$

其中：

$$\boxed{\begin{aligned}
a_1 &= \begin{vmatrix} x_2 & y_2 \\ x_3 & y_3 \end{vmatrix} = x_2 y_3 - x_3 y_2 \\
b_1 &= -\begin{vmatrix} 1 & y_2 \\ 1 & y_3 \end{vmatrix} = y_2 - y_3 \\
c_1 &= \begin{vmatrix} 1 & x_2 \\ 1 & x_3 \end{vmatrix} = x_3 - x_2
\end{aligned}}\quad\text{（1,2,3 轮换）}$$

$$A = \frac{1}{2}\begin{vmatrix} 1 & x_1 & y_1 \\ 1 & x_2 & y_2 \\ 1 & x_3 & y_3 \end{vmatrix}$$

此处 $A = \Delta_e$（与 5.2 节记号一致），表示三角形单元的面积。

同理可解出 v 方向的 $\beta_4, \beta_5, \beta_6$（将 $u$ 替换为 $v$）。

**整理为形函数形式**：

$$\begin{cases}
u(x,y) = N_1 u_1 + N_2 u_2 + N_3 u_3 \\
v(x,y) = N_1 v_1 + N_2 v_2 + N_3 v_3
\end{cases}$$

其中 $N_s = \dfrac{1}{2A}(a_s + b_s x + c_s y),\;(s=1,2,3)$

#### (B) 单元应变场——$[B]$ 矩阵的构造

由几何方程：

$$\boldsymbol{\varepsilon}_{(3\times 1)} = [\partial]_{(3\times 2)}\mathbf{u}_{(2\times 1)} = [\partial][N]\mathbf{q}^e = [B]\mathbf{q}^e$$

代入 $[N]$ 的表达式，得**$[B]$ 矩阵**（$3 \times 6$ 常数矩阵）：

$$\boxed{[B]_{(3\times 6)} = \frac{1}{2A}\begin{pmatrix}
b_1 & 0 & b_2 & 0 & b_3 & 0 \\
0 & c_1 & 0 & c_2 & 0 & c_3 \\
c_1 & b_1 & c_2 & b_2 & c_3 & b_3
\end{pmatrix}}$$

可写为分块形式：$[B] = \begin{bmatrix} B_1 & B_2 & B_3 \end{bmatrix}$，其中：

$$B_s = \frac{1}{2A}\begin{pmatrix} b_s & 0 \\ 0 & c_s \\ c_s & b_s \end{pmatrix},\quad (s=1,2,3)$$

#### (C) 单元应力场

对于**平面应力问题**，本构矩阵为：

> ❌ **易错点**：平面应力和平面应变的 D 矩阵极易混淆。记住：平面应力是对薄板（$\sigma_z = 0$），用 $E/(1-\mu^2)$；平面应变是对长柱体（$\varepsilon_z = 0$），需要替换 $E \to E/(1-\mu^2)$、$\mu \to \mu/(1-\mu)$。考试时如果不确定，就写清楚"按平面应力计算"或"按平面应变计算"。

$$\boxed{[D]_{(3\times 3)} = \frac{E}{1-\mu^2}\begin{pmatrix}
1 & \mu & 0 \\
\mu & 1 & 0 \\
0 & 0 & \frac{1-\mu}{2}
\end{pmatrix}}$$

对于**平面应变问题**，只需将 $E \to \frac{E}{1-\mu^2}$、$\mu \to \frac{\mu}{1-\mu}$。

$$\boxed{[D] = \frac{E(1-\mu)}{(1+\mu)(1-2\mu)}\begin{pmatrix}
1 & \frac{\mu}{1-\mu} & 0 \\
\frac{\mu}{1-\mu} & 1 & 0 \\
0 & 0 & \frac{1-2\mu}{2(1-\mu)}
\end{pmatrix}}$$

应力场：

$$\boxed{\boldsymbol{\sigma} = [D]\boldsymbol{\varepsilon} = [D][B]\mathbf{q}^e = [S]\mathbf{q}^e}$$

其中 $[S]$ 称为**应力矩阵**。对 CST 单元，$[B]$ 为常数，因此 $[S]$ 也是常数——单元内应力和应变处处相同（这就是 CST 名称的来源）。

**注**：因此三角形 3 节点单元也称**常应变（常应力）单元**。在应变梯度（应力梯度）变化剧烈的区域需要加密网格，否则会产生较大误差。

#### (D) 单元刚度矩阵的显式公式

$$\begin{aligned}
\mathbf{K}_{(6\times 6)}^e &= \int_{\Omega^e} \mathbf{B}^T\mathbf{D}\mathbf{B}\,dV = t\int_{A^e} \mathbf{B}^T\mathbf{D}\mathbf{B}\,dA \quad (\text{其中 } t \text{ 是平面单元厚度}) \\
&= \mathbf{B}^T\mathbf{D}\mathbf{B}\,tA \quad (\text{因 } [B] \text{ 为常数})
\end{aligned}$$

分块形式：

$$\mathbf{K}_{(6\times 6)}^e = \begin{pmatrix}
k_{11} & k_{12} & k_{13} \\
k_{21} & k_{22} & k_{23} \\
k_{31} & k_{32} & k_{33}
\end{pmatrix}$$

其中每个子块 $k_{rs}$ 为 $2 \times 2$ 矩阵：

$$\boxed{k_{rs} = B_r^T D B_s\,tA = \frac{Et}{4(1-\mu^2)A}\begin{pmatrix}
b_r b_s + \frac{1-\mu}{2}c_r c_s & \mu b_r c_s + \frac{1-\mu}{2}c_r b_s \\
\mu c_r b_s + \frac{1-\mu}{2}b_r c_s & c_r c_s + \frac{1-\mu}{2}b_r b_s
\end{pmatrix}}$$

展开记为：

$$k_{rs} = \frac{Et}{4(1-\mu^2)A}\begin{pmatrix} k_1 & k_3 \\ k_2 & k_4 \end{pmatrix}$$

其中：

$$\begin{aligned}
k_1 &= b_r b_s + \frac{1-\mu}{2}c_r c_s \\
k_2 &= \mu c_r b_s + \frac{1-\mu}{2}b_r c_s \\
k_3 &= \mu b_r c_s + \frac{1-\mu}{2}c_r b_s \\
k_4 &= c_r c_s + \frac{1-\mu}{2}b_r b_s
\end{aligned}$$

#### (E) 等效节点荷载

$$\boxed{\mathbf{P}_{(6\times 1)}^e = \int_{\Omega^e} \mathbf{N}^T\mathbf{f}\,dV + \int_{S_\sigma^e} \mathbf{N}^T\mathbf{T}\,dS}$$

对均布体力（如重力），可用面积坐标简化积分。对复杂分布，需数值积分。

#### (F) 单元刚度方程

对节点位移取势能的一次极值，得单元刚度方程：

$$\boxed{\mathbf{K}_{(6\times 6)}^e \mathbf{q}_{(6\times 1)}^e = \mathbf{P}_{(6\times 1)}^e}$$

---

### 5.4.3 完整数值算例


#### 算例一：悬臂梁受集中力（2 单元 4 节点）

**问题描述**：悬臂梁长度 2 m, 高度 1 m，右端受集中力 F。材料弹性模量 E，泊松比 $\mu = \frac{1}{3}$，厚度 t。按平面应力问题求解节点位移和支座反力。

**第 1 步：结构离散化与编号**

划分为 2 个三角形单元。力 F 按静力等效原则分配到节点 1 和 2（各 $F/2$）。

节点位移列阵：

$$\mathbf{q} = \begin{pmatrix} u_1 & v_1 & u_2 & v_2 & u_3 & v_3 & u_4 & v_4 \end{pmatrix}^T$$

节点外力列阵：

$$\mathbf{F} = \begin{pmatrix} 0 & -\frac{F}{2} & 0 & -\frac{F}{2} & 0 & 0 & 0 & 0 \end{pmatrix}^T$$

支座反力列阵：

$$\mathbf{R} = \begin{pmatrix} 0 & 0 & 0 & 0 & R_{3x} & R_{3y} & R_{4x} & R_{4y} \end{pmatrix}^T$$

**第 2 步：各单元的描述**

两单元取相同局部编号时具有相同的刚度矩阵：

$$\mathbf{K}^{(1)} = \mathbf{K}^{(2)} = \begin{pmatrix}
k_{ii} & k_{ij} & k_{im} \\
k_{ji} & k_{jj} & k_{jm} \\
k_{mi} & k_{mj} & k_{mm}
\end{pmatrix}$$

代入几何参数计算后（过程略），得数值矩阵：

$$\mathbf{K}^{(1)} = \mathbf{K}^{(2)} = \frac{9Et}{32}\begin{pmatrix}
1      & 0        & 0       & \frac{2}{3} & -1      & -\frac{2}{3} \\
0      & \frac{1}{3} & \frac{2}{3} & 0        & -\frac{2}{3} & -\frac{1}{3} \\
0      & \frac{2}{3} & \frac{4}{3} & 0        & -\frac{4}{3} & -\frac{2}{3} \\
\frac{2}{3} & 0        & 0       & 4            & -\frac{2}{3} & -4 \\
-1     & -\frac{2}{3} & -\frac{4}{3} & -\frac{2}{3} & \frac{7}{3} & \frac{4}{3} \\
-\frac{2}{3} & -\frac{1}{3} & -\frac{2}{3} & -4            & \frac{4}{3} & \frac{13}{3}
\end{pmatrix}$$

**第 3 步：组装总体刚度矩阵**

按位移自由度"对号入座"组装：

$$\mathbf{K}_{(8\times 8)} = \mathbf{K}^{(1)} + \mathbf{K}^{(2)}$$

$$\mathbf{K} = \frac{9Et}{32}\begin{pmatrix}
\frac{7}{3} & \frac{4}{3} & -\frac{4}{3} & -\frac{2}{3} & -1 & -\frac{2}{3} & 0 & 0 \\
\frac{4}{3} & \frac{13}{3} & -\frac{2}{3} & -4 & -\frac{2}{3} & -\frac{1}{3} & 0 & 0 \\
-\frac{4}{3} & -\frac{2}{3} & \frac{7}{3} & 0 & 0 & \frac{4}{3} & -1 & -\frac{2}{3} \\
-\frac{2}{3} & -4 & 0 & \frac{13}{3} & \frac{4}{3} & 0 & -\frac{2}{3} & -\frac{1}{3} \\
-1 & -\frac{2}{3} & 0 & \frac{4}{3} & \frac{7}{3} & 0 & -\frac{4}{3} & -\frac{2}{3} \\
-\frac{2}{3} & -\frac{1}{3} & \frac{4}{3} & 0 & 0 & \frac{13}{3} & -\frac{2}{3} & -4 \\
0 & 0 & -1 & -\frac{2}{3} & -\frac{4}{3} & -\frac{2}{3} & \frac{7}{3} & \frac{4}{3} \\
0 & 0 & -\frac{2}{3} & -\frac{1}{3} & -\frac{2}{3} & -4 & \frac{4}{3} & \frac{13}{3}
\end{pmatrix}$$

**第 4 步：边界条件处理与求解**

位移边界条件：$u_3 = v_3 = u_4 = v_4 = 0$。（删行法：删除第 5-8 行/列）

> ❌ **易错点**：删行删列法操作时，要注意原方程的行列编号。删除约束自由度对应的行和列后，剩余矩阵的维度变为"总自由度 - 约束数"。此处 8 个自由度删去 4 个，剩余 4×4 矩阵。注意缩聚后刚度矩阵的对称性仍然保持。

缩聚后的方程：

$$\frac{9Et}{32}\begin{pmatrix}
\frac{7}{3} & \frac{4}{3} & -\frac{4}{3} & -\frac{2}{3} \\
\frac{4}{3} & \frac{13}{3} & -\frac{2}{3} & -4 \\
-\frac{4}{3} & -\frac{2}{3} & \frac{7}{3} & 0 \\
-\frac{2}{3} & -4 & 0 & \frac{13}{3}
\end{pmatrix}
\begin{pmatrix} u_1 \\ v_1 \\ u_2 \\ v_2 \end{pmatrix}
= \begin{pmatrix} 0 \\ -F/2 \\ 0 \\ -F/2 \end{pmatrix}$$

**求解得节点位移（数值结果）**：

$$\boxed{\begin{pmatrix} u_1 & v_1 & u_2 & v_2 \end{pmatrix}^T = \frac{F}{Et}\begin{pmatrix} 1.88 & -8.99 & -1.50 & -8.42 \end{pmatrix}^T}$$

**第 5 步：支座反力计算**

将位移代入总刚方程的回代部分：

$$\begin{aligned}
R_{3x} &= \frac{9Et}{32}(-u_1 - \frac{2}{3}v_1 + \frac{4}{3}v_2) = -2F \\
R_{3y} &= \frac{9Et}{32}(-\frac{2}{3}u_1 - \frac{1}{3}v_1 + \frac{4}{3}u_2) = -0.07F \\
R_{4x} &= \frac{9Et}{32}(-u_2 - \frac{2}{3}v_2) = 2F \\
R_{4y} &= \frac{9Et}{32}(-\frac{2}{3}u_2 - \frac{1}{3}v_2) = 1.07F
\end{aligned}$$

**验证**：$\sum F_x = R_{3x} + R_{4x} = -2F + 2F = 0$；$\sum F_y = -F/2 - F/2 + R_{3y} + R_{4y} \approx -F + F = 0$。支座反力与外力构成平衡力系。

> 💡 **理解关键**：最后的平衡验证不是多余的——它是检验计算结果正确性的重要手段。如果支座反力不能平衡外力，说明组装或求解过程中有错误。FEM 计算中，力平衡检查是最基本的"合理性检验"。$R_{3y} = -0.07F$ 非常小是因为仅用 2 个 CST 单元模拟弯曲时，单元过刚导致剪应力被低估——这是位移元下限性的一个直观体现。

---

#### 算例二：正方形薄板均匀拉伸（对称性利用）

**问题描述**：均质正方形薄板，厚度 $t = 1\text{ cm}$，上下两边受均匀拉伸 $q = 10^6\text{ N/m}$。$E = 1000$，$\mu = 1/3$，忽略自重。用 FEM 求解应力分量。

**第 1 步：力学建模**

板的长宽远大于厚度，荷载在板面内且沿厚度均匀分布——这是**平面应力问题**。利用对称性取**四分之一结构**分析。

> 💡 **理解关键**：对称性利用是简化 FEM 计算的重要技巧。本算例中结构双向对称、荷载对称，因此可以只取 1/4 分析——对称面上的边界条件是法向位移为零（$u = 0$ 或 $v = 0$），切向自由。考试时如果题目给了对称结构，不利用对称性直接算全模型，虽然正确但计算量大得多。但如果荷载不对称，就不能直接用对称性了。

**第 2 步：结构离散化**

四分之一结构离散为 2 个三角形单元，4 个节点。节点坐标：

| 坐标/节点 | 1     | 2     | 3     | 4     |
|-----------|-------|-------|-------|-------|
| x         | 0     | 1     | 1     | 0     |
| y         | 0     | 0     | 1     | 1     |

**第 3 步：计算单元刚度矩阵**

**单元 1**（节点编号 1→2→3，对应局部 i,j,m）：

计算几何参数：
$$b_1 = y_2 - y_3 = -1,\; b_2 = y_3 - y_1 = 1,\; b_3 = y_1 - y_2 = 0$$
$$c_1 = -(x_2 - x_3) = 0,\; c_2 = -(x_3 - x_1) = -1,\; c_3 = -(x_1 - x_2) = 1$$
$$\Delta^{(1)} = \frac{1}{2}(b_2 c_3 - b_3 c_2) = \frac{1}{2}$$

预备常数：$\frac{1-\mu}{2} = \frac{1}{3}$，$\frac{Et}{4(1-\mu^2)\Delta} = \frac{9E}{16}$

计算各子块（$k_{rs} = \frac{9E}{16}\begin{pmatrix} k_1 & k_3 \\ k_2 & k_4 \end{pmatrix}$）：

$$k_{11} = \frac{9E}{16}\begin{pmatrix} (-1)(-1)+\frac{1}{3}\cdot0\cdot0 & \frac{1}{3}\cdot(-1)\cdot0+\frac{1}{3}\cdot0\cdot(-1) \\ \frac{1}{3}\cdot0\cdot(-1)+\frac{1}{3}\cdot(-1)\cdot0 & 0\cdot0+\frac{1}{3}\cdot(-1)\cdot(-1) \end{pmatrix} = \frac{9E}{16}\begin{pmatrix} 1 & 0 \\ 0 & \frac{1}{3} \end{pmatrix}$$

$$k_{12} = \frac{9E}{16}\begin{pmatrix} -1 & \frac{1}{3} \\ \frac{1}{3} & -\frac{1}{3} \end{pmatrix}, \quad k_{13} = \frac{9E}{16}\begin{pmatrix} 0 & -\frac{1}{3} \\ -\frac{1}{3} & 0 \end{pmatrix}$$

$$k_{22} = \frac{9E}{16}\begin{pmatrix} \frac{4}{3} & -\frac{2}{3} \\ -\frac{2}{3} & \frac{4}{3} \end{pmatrix}, \quad k_{23} = \frac{9E}{16}\begin{pmatrix} -\frac{1}{3} & \frac{1}{3} \\ \frac{1}{3} & -1 \end{pmatrix}, \quad k_{33} = \frac{9E}{16}\begin{pmatrix} \frac{1}{3} & 0 \\ 0 & 1 \end{pmatrix}$$

**单元 1 的刚度矩阵**（$6\times 6$ 对称）：

$$[K]_{6\times 6}^{(1)} = \frac{9E}{16}\begin{pmatrix}
1 & 0 & -1 & \frac{1}{3} & 0 & -\frac{1}{3} \\
0 & \frac{1}{3} & \frac{1}{3} & -\frac{1}{3} & -\frac{1}{3} & 0 \\
-1 & \frac{1}{3} & \frac{4}{3} & -\frac{2}{3} & -\frac{1}{3} & \frac{1}{3} \\
\frac{1}{3} & -\frac{1}{3} & -\frac{2}{3} & \frac{4}{3} & \frac{1}{3} & -1 \\
0 & -\frac{1}{3} & -\frac{1}{3} & \frac{1}{3} & \frac{1}{3} & 0 \\
-\frac{1}{3} & 0 & \frac{1}{3} & -1 & 0 & 1
\end{pmatrix}$$

**单元 2**（节点编号 3→4→1 对应按 1→2→3 顺序）：因坐标关系相同，刚度矩阵数值与单元 1 完全相同，只是下标对应关系为 $(3,4,1) \leftrightarrow (1,2,3)$。

**第 4 步：组装总体刚度矩阵**

$$[K]_{8\times 8} = \begin{pmatrix}
[K_{11}]^{(1)+(2)} & [K_{12}]^{(1)} & [K_{13}]^{(1)+(2)} & [K_{14}]^{(2)} \\
[K_{21}]^{(1)} & [K_{22}]^{(1)} & [K_{23}]^{(1)} & 0 \\
[K_{31}]^{(1)+(2)} & [K_{32}]^{(1)} & [K_{33}]^{(1)+(2)} & [K_{34}]^{(2)} \\
[K_{41}]^{(2)} & 0 & [K_{43}]^{(2)} & [K_{44}]^{(2)}
\end{pmatrix}$$

利用对称性 $[K_{rs}] = [K_{sr}]^T$ 和编号对应关系，计算得到：

$$[K]_{8\times 8} = \frac{3E}{16}\begin{pmatrix}
4 & 0 & -3 & 1 & 0 & -2 & -1 & 1 \\
0 & 4 & 1 & -1 & -2 & 0 & 1 & -3 \\
-3 & 1 & 4 & -2 & -1 & 1 & 0 & 0 \\
1 & -1 & -2 & 4 & 0 & -3 & 0 & 0 \\
0 & -2 & -1 & 0 & 4 & 0 & -3 & 1 \\
-2 & 0 & 1 & -3 & 0 & 4 & 1 & -1 \\
-1 & 1 & 0 & 0 & -3 & 1 & 4 & -2 \\
1 & -3 & 0 & 0 & 1 & -1 & -2 & 4
\end{pmatrix}$$

**第 5 步：引入边界条件**

约束条件（对称面上）：$u_1 = v_1 = 0$；$v_2 = 0$；$u_4 = 0$

等效节点力：$\{F\} = \{0, 0, 0, 0, 0, q/2, 0, q/2\}^T$（上边界均布荷载 $q$ 分配到节点 3 和 4）

删除位移为零对应的行列（1,2,4,7），缩聚方程：

$$\frac{3E}{16}\begin{pmatrix}
4 & -1 & 1 & 0 \\
-1 & 4 & 0 & 1 \\
1 & 0 & 4 & -1 \\
0 & 1 & -1 & 4
\end{pmatrix}
\begin{pmatrix} u_2 \\ u_3 \\ v_3 \\ v_4 \end{pmatrix}
= \begin{pmatrix} 0 \\ 0 \\ q/2 \\ q/2 \end{pmatrix}$$

**求解得**：

$$\begin{pmatrix} u_2 & u_3 & v_3 & v_4 \end{pmatrix}^T = \frac{q}{E}\begin{pmatrix} -1/3 & -1/3 & 1 & 1 \end{pmatrix}^T$$

即全局位移向量：$\{\delta\} = \dfrac{q}{E}\begin{pmatrix} 0 & 0 & -1/3 & 0 & -1/3 & 1 & 0 & 1 \end{pmatrix}^T$

**第 6 步：应力回算**

单元 1 的应力矩阵及应力计算：

$$\{\sigma\}^{(1)} = [S]^{(1)}\{\delta\}^{(1)} = \frac{3E}{8}\begin{pmatrix}
-3 & 0 & 3 & -1 & 0 & 1 \\
-1 & 0 & 1 & -3 & 0 & 3 \\
0 & -1 & -1 & 1 & 1 & 0
\end{pmatrix}
\begin{pmatrix} 0 \\ 0 \\ -q/3E \\ 0 \\ -q/3E \\ q/E \end{pmatrix}
= q\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$$

$$\{\sigma\}^{(2)} = q\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\quad\text{（同理）}$$

**结果解读**：$\sigma_x = 0$，$\sigma_y = q$，$\tau_{xy} = 0$。全部单元应力相同，与解析解完全一致（该问题本身就是常应力状态，CST 单元精确捕捉）。

> 💡 **理解关键**：这个算例有一个重要启示——当真实解恰好是常应力状态时，即使只用 2 个 CST 单元也能得到**精确解**。这完美验证了 CST 单元的"常应变完备性"：任何常应变模态都能被 CST 精确表示。但如果应力不均匀（如弯曲），同样的 2 个单元就不够用了——这也是算例一中 $R_{3y}$ 不精确的原因。

---

#### 算例三：4 单元平面应力问题——乘大数法

**问题描述**：平面应力问题，$\mu = 0$，单元厚度 $t=1$。结构离散化如图，共 4 个单元、6 个节点。求节点位移和单元应力、应变。

> ❌ **易错点**：$\mu = 0$ 是为了简化计算而设的假设，真实工程材料（钢材 $\mu \approx 0.3$，混凝土 $\mu \approx 0.2$）泊松比不为 0。考试时注意题目给定的材料参数——如果 $\mu = 0$，D 矩阵简化为对角阵 $\text{diag}(E, E, E/2)$，计算量大大简化，但求解思路完全不变。

**第 1 步：计算几何参数**

对于单元 1, 2, 4（形状相同）：
$$b_i = 0,\; b_j = -a,\; b_m = a$$ 其中 $a$ 是单元的特征尺寸（此处为等腰直角三角形单元的直角边长）。
$$c_i = a,\; c_j = -a,\; c_m = 0$$
$$A = a^2/2$$

对于单元 3：
$$b_i = -a,\; b_j = 0,\; b_m = a$$
$$c_i = 0,\; c_j = -a,\; c_m = a$$
$$A = a^2/2$$

**第 2 步：计算单元刚度矩阵**

**单元 1, 2, 4 的刚度矩阵**：

$$k^{(1,2,4)} = \frac{E}{4}\begin{pmatrix}
1 & 0 & -1 & -1 & 0 & 1 \\
0 & 2 & 0 & -2 & 0 & 0 \\
-1 & 0 & 3 & 1 & -2 & -1 \\
-1 & -2 & 1 & 3 & 0 & -1 \\
0 & 0 & -2 & 0 & 2 & 0 \\
1 & 0 & -1 & -1 & 0 & 1
\end{pmatrix}\begin{matrix} \leftarrow i \\ \leftarrow j \\ \leftarrow m \end{matrix}$$

**单元 3 的刚度矩阵**：

$$k^{(3)} = \frac{E}{4}\begin{pmatrix}
2 & 0 & 0 & 0 & -2 & 0 \\
0 & 1 & 1 & 0 & -1 & -1 \\
0 & 1 & 1 & 0 & -1 & -1 \\
0 & 0 & 0 & 2 & 0 & -2 \\
-2 & -1 & -1 & 0 & 3 & 1 \\
0 & -1 & -1 & -2 & 1 & 3
\end{pmatrix}\begin{matrix} \leftarrow i \\ \leftarrow j \\ \leftarrow m \end{matrix}$$

**单元节点编号对应表**：

| 单元号 | i   | j   | m   |
|--------|-----|-----|-----|
| 1      | 1   | 2   | 3   |
| 2      | 2   | 4   | 5   |
| 3      | 2   | 5   | 3   |
| 4      | 3   | 5   | 6   |

**第 3 步：扩阶并组装总体刚度矩阵**

将各单元刚度矩阵扩展为 $12\times 12$（6 个节点 $\times$ 2 个自由度），按对应关系"对号入座"叠加，得总体刚度矩阵（$12\times 12$，对称阵）：

> ⚠️ **重难点**："扩阶再叠加"是理解总刚组装的最佳视角。每个单元刚度矩阵先扩展为总自由度规模的零矩阵（仅在与该单元节点对应的位置填入局部刚度值），然后将所有扩展矩阵相加。这个过程在编程中是两重循环：遍历单元 → 遍历该单元的局部节点对 → 将 $k_{rs}^e$ 加到 $K_{global(r), global(s)}$。

$$K = \sum_{e=1}^4 k_{12\times 12}^{(e)} = \frac{E}{4}\begin{pmatrix}
1 & 0 & -1 & -1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
-1 & 0 & 6 & 1 & -4 & -1 & -1 & -1 & 0 & 1 & 0 & 0 \\
-1 & -2 & 1 & 6 & -1 & -2 & 0 & -2 & 1 & 0 & 0 & 0 \\
0 & 0 & -4 & -1 & 6 & 1 & 0 & 0 & -2 & -1 & 0 & 1 \\
1 & 0 & -1 & -2 & 1 & 6 & 0 & 0 & -1 & -4 & 0 & 0 \\
0 & 0 & -1 & 0 & 0 & 0 & 3 & 1 & -2 & -1 & 0 & 0 \\
0 & 0 & -1 & -2 & 0 & 0 & 1 & 3 & 0 & -1 & 0 & 0 \\
0 & 0 & 0 & 1 & -2 & -1 & -2 & 0 & 6 & 1 & -2 & -1 \\
0 & 0 & 1 & 0 & -1 & -4 & -1 & -1 & 1 & 6 & 0 & -1 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & 2 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & -1 & -1 & 0 & 1
\end{pmatrix}$$

**第 4 步：边界条件处理——乘大数法**

边界条件：$u_1 = u_2 = u_4 = v_4 = v_5 = v_6 = 0$

荷载向量：$\mathbf{R} = \begin{pmatrix} 0 & -P & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{pmatrix}^T$

将 $K_{11}, K_{33}, K_{77}, K_{88}, K_{10,10}, K_{12,12}$ 乘以 $10^{15}$，相应右端项置零：

$$\frac{E}{4}\begin{pmatrix}
1\times 10^{15} & 0 & -1 & -1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
-1 & 0 & 6\times10^{15} & 1 & -4 & -1 & -1 & -1 & 0 & 1 & 0 & 0 \\
\cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots & \cdots
\end{pmatrix}
\begin{pmatrix} u_1 \\ v_1 \\ u_2 \\ v_2 \\ u_3 \\ v_3 \\ u_4 \\ v_4 \\ u_5 \\ v_5 \\ u_6 \\ v_6 \end{pmatrix}
= \begin{pmatrix} 0 \\ -P \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$$

**第 5 步：求解节点位移**

$$\boxed{\boldsymbol{\delta} = \frac{P}{E}\begin{pmatrix} 0 \\ -3.252 \\ 0 \\ -1.252 \\ -0.088 \\ -0.374 \\ 0 \\ 0 \\ 0.176 \\ 0 \\ 0.176 \\ 0 \end{pmatrix} \begin{matrix} \leftarrow u_1,v_1 \\ \leftarrow u_2,v_2 \\ \leftarrow u_3,v_3 \\ \leftarrow u_4,v_4 \\ \leftarrow u_5,v_5 \\ \leftarrow u_6,v_6 \end{matrix}}$$

**第 6 步：应力与应变回算**

**单元 1**（节点 i=1, j=2, m=3）：

$$\boldsymbol{\varepsilon}^{(1)} = \frac{1}{a^2}\begin{pmatrix}
0 & 0 & -a & 0 & a & 0 \\
0 & a & 0 & -a & 0 & 0 \\
a & 0 & -a & -a & 0 & a
\end{pmatrix}
\begin{pmatrix} 0 \\ -3.252 \\ 0 \\ -1.252 \\ -0.088 \\ -0.374 \end{pmatrix}\frac{P}{E} = \frac{P}{Ea}\begin{pmatrix} -0.088 \\ -2.000 \\ 0.880 \end{pmatrix}$$

$$\boldsymbol{\sigma}^{(1)} = E\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0.5 \end{pmatrix}\boldsymbol{\varepsilon}^{(1)} = \frac{P}{a}\begin{pmatrix} -0.088 \\ -2.000 \\ 0.440 \end{pmatrix}$$

**单元 2**（i=2, j=4, m=5）：

$$\boldsymbol{\varepsilon}^{(2)} = \frac{P}{Ea}\begin{pmatrix} 0.176 \\ -1.252 \\ 0 \end{pmatrix}, \quad \boldsymbol{\sigma}^{(2)} = \frac{P}{a}\begin{pmatrix} 0.176 \\ -1.252 \\ 0 \end{pmatrix}$$

**单元 3**（i=2, j=5, m=3）：

$$\boldsymbol{\varepsilon}^{(3)} = \frac{P}{Ea}\begin{pmatrix} -0.088 \\ -0.374 \\ 0.614 \end{pmatrix}, \quad \boldsymbol{\sigma}^{(3)} = \frac{P}{a}\begin{pmatrix} -0.088 \\ -0.374 \\ 0.307 \end{pmatrix}$$

**单元 4**（i=3, j=5, m=6）：

$$\boldsymbol{\varepsilon}^{(4)} = \frac{P}{Ea}\begin{pmatrix} 0 \\ -0.374 \\ -0.264 \end{pmatrix}, \quad \boldsymbol{\sigma}^{(4)} = \frac{P}{a}\begin{pmatrix} 0 \\ -0.374 \\ -0.132 \end{pmatrix}$$

---

### 5.4.4 求解流程总结

以具有常体力（如重力）和线性面力的平面问题为例，完整求解步骤如下：

**1. 输入原始信息**
- 节点坐标 $(x_i, y_i),\; i = 1,\ldots,NP$
- 单元信息：每个单元的 $(n_1, n_2, n_3)$ 逆时针排列 + 材料类型号
- 材料特性：弹性常数、密度（重力用）
- 集中力信息：节点号及其量值
- 面力信息：受面力线元两端的编号及荷载密度
- 约束信息：约束节点号及 $x, y$ 方向的约束情况

**2. 建立求解方程**
- 清零：$K_{ij} = 0,\; F_i = 0$
- 单元扫描：逐一计算 $[k]_e$ 和 $\{F\}_e$，对号组装入总刚（可利用对称性只算一半，带宽外元素可忽略）
- 线单元分析：计算面力的等效节点荷载
- 集中力直接累加

**3. 边界条件处理**（"变 1 法"或"乘大数法"）

**4. 求解线性方程组** $[K]\{\delta\} = \{F\}$（利用对阵性、稀疏性、带状分布）

**5. 应力回算**：对每个单元 $\boldsymbol{\sigma} = [D][B]\{\delta\}_e$（CST 单元内应力为常数）

---

## 检查你的理解（5.4）

1. CST 中 $[B]$ 矩阵的每个元素 $b_i, c_i$ 的意义是什么？为什么 $[B]$ 是常数矩阵？
2. 算例二中利用了什么对称性简化计算？如果在非对称情况下需要如何处理？
3. 算例三中 $\mu = 0$ 这个条件使 $[D]$ 矩阵简化为什么形式？这与实际工程材料是否吻合？

---

## 5.5 进一步讨论（Further Discussions）

**本节目标**：掌握单元刚度矩阵的四大特性、FEM 收敛的充要条件、分片试验的原理、位移解的精度估计和下限性证明。这些是判断 FEM 解是否可靠的理论基础。

### 5.5.1 单元刚度矩阵的四大特性


#### 1. 物理意义

单元刚度方程 $\mathbf{K}^e\mathbf{a}^e = \mathbf{P}^e$ 本质上是**单元节点的平衡方程**——每个节点在 x 和 y 方向各有一个平衡方程（对 CST 共 6 个）。

左边为用节点位移表示的单元节点内力，右边为外荷载的等效节点力。

刚度系数 $K_{ij}$ 的物理含义是：**当单元的第 j 个自由度产生单位位移而其他自由度位移为零时，需要在第 i 个自由度上施加的节点力大小**。

验证：令 $a_1 = 1(u_i = 1),\; a_2 = a_3 = \cdots = a_6 = 0$，由平衡方程得：

$$\begin{pmatrix} K_{11} \\ K_{21} \\ \vdots \\ K_{61} \end{pmatrix}_{a_1=1} = \begin{pmatrix} P_1 \\ P_2 \\ \vdots \\ P_6 \end{pmatrix}$$

即第一列即为使节点 i 产生单位 x 向位移所需的各节点力。

由于单元在这组节点力作用下处于平衡状态，各方向合力为零：

$$x: K_{11} + K_{31} + K_{51} = 0,\quad y: K_{21} + K_{41} + K_{61} = 0$$

#### 2. 对称性

单元刚度矩阵是对称的：$K_{ij} = K_{ji}$。这在前述的加权残量法和 Galerkin 变分法中均已讨论。基于最小势能原理的 FEM 其本质相同——$\mathbf{K}^e = \int \mathbf{B}^T\mathbf{D}\mathbf{B}\,dV$ 自动对称。

#### 3. 奇异性

对 $6 \times 6$ 的线性三角形单元刚度矩阵，只有 3 行（或 3 列）是独立的。证明：

$$\begin{aligned}
K_{1j} + K_{3j} + K_{5j} &= K_{j1} + K_{j3} + K_{j5} = 0 \\
K_{2j} + K_{4j} + K_{6j} &= K_{j2} + K_{j4} + K_{j6} = 0
\end{aligned}$$

奇异的物理原因：未引入边界条件前，单元有刚体位移自由度，位移解不唯一。引入边界约束后总刚才变为正定。

> 💡 **理解关键**：奇异性 = 无约束 = 刚体位移未消除。CST 单元有 3 个刚体位移模态（x 向平动、y 向平动、面内转动），每个模态对应刚度矩阵的一个零特征值。6×6 矩阵秩为 3（6 - 3 个刚体模态 = 3），所以只有 3 行独立。这个秩分析推广到三维四面体单元：12×12 矩阵秩为 6（12 - 6 个刚体模态 = 6）。

#### 4. 稀疏性与带状分布

一个节点仅通过相邻单元与周围少数节点相关（该节点形函数的"支集区域"）。因此尽管总刚阶数很高，非零元素很少——**稀疏性**。

若节点编号合理，非零元素集中在带状区域内——**带状分布**。这对大型方程求解效率至关重要。

---

### 5.5.2 收敛准则（Requirements of Convergence）


FEM 可以看作是 Ritz 法的一种特殊形式——试函数定义在单元（子域）而非全局域。因此 FEM 的收敛性要求与 Ritz 法一致。

Ritz 法的收敛要求：试函数必须**完备**且**连续**。当 $n\to\infty$ 时近似解收敛到精确解。

#### (1) 完备性要求（Completeness）

设泛函 $\Pi$ 中包含场函数 $\varphi$ 及其最高 $m$ 阶导数。则收敛的**必要条件**：单元内的试函数（插值函数）必须是至少 $m$ **次完备多项式**。

即试函数必须包含：自身以及直到 $m$ 阶导数为常数的所有项。

当单元尺寸趋于零时，$m$ 阶以内导数趋于常数。若试函数为 $p$ 阶完备多项式（$p \geq m$），各单元的泛函才能逼近精确值。

**对弹性力学平面问题**：$\Pi_p^e = \int_{\Omega_e} \frac{1}{2}\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \cdots$ 含位移一阶导数，$m=1$。

完备性要求位移函数至少是 $x, y$ 的**一次完备多项式**——即至少包含：
$$u = a_1 + a_2 x + a_3 y,\quad v = a_4 + a_5 x + a_6 y$$

**物理意义**：上式中常数项代表**刚体位移**模态，一次项代表**常应变**模态。完备性 = FEM 解必须能反映单元的刚体位移和常应变状态。

#### (2) 协调性要求（Compatibility）

如果泛函中最高导数为 $m$ 阶，则单元间界面上的试函数必须具有 $C^{m-1}$ 阶连续性——即界面上函数直到 $m-1$ 阶导数连续。

**对弹性力学平面问题**（$m=1$）：要求 $C^0$ 连续——位移函数本身在单元界面上连续，但不要求导数连续。

满足上述完备性和协调性要求的单元称为**协调元**。当单元同时满足完备性和协调性时，FEM 解收敛。

线性三角形单元（CST）显然同时满足完备性和协调性，因此其解是收敛的。

**重要补充**：
- 对于二维/三维弹性力学，泛函一般含一阶导数（$m=1$），仅需 $C^0$ 连续——即函数本身在界面连续
- 对于梁、板、壳问题，泛函含二阶导数（$m=2$），需 $C^1$ 连续——即函数及其一阶导数在界面连续

> 🔗 **跨章连接**：$C^0$ vs $C^1$ 连续性的区别在这里首次明确提出，但在第 6 章构造梁/板单元形函数时才会真正体会到它的"杀伤力"——$C^1$ 连续性要求意味着需要 Hermite 插值（不仅节点值连续，导数也连续），比 Lagrange 插值复杂得多。这也是为什么板壳单元的构造比平面/三维单元难一个数量级。

---

### 5.5.3 分片试验（Patch Test）

B.M. Irons 提出的分片试验是判断单元是否满足收敛性要求的简易方法——试验通过是收敛的**充分条件**。

> ⚠️ **重难点**：分片试验是 FEM 中较抽象的概念。理解它的关键是"必要 vs 充分"的区分：满足完备性+协调性 → 一定收敛；但即使不满足协调性，只要通过分片试验，也一定收敛。分片试验的核心思想是：无需逐一检查协调条件，直接用"常应变加载 → 检查内节点是否给出常应变"来一次性验证。不通过 = 肯定不收敛；通过 = 一定收敛。

**试验步骤**：

1. 用若干单元拼成一小块弹性体，确保至少有一个内节点被单元完全包围（内节点可改变以获得任意形状的单元）
2. 在边界节点上施加与**常应变状态**对应的位移或外力
3. 内节点既不施加外力，也不施加位移约束

**判定准则**：若计算得到的单元位移、应力和应变与给定的常应变状态一致，则该单元通过分片试验。

**意义**：通过分片试验的单元，即使不严格满足协调性条件（称为**非协调元**），其 FEM 解仍然收敛。

**非协调元有时精度更好**——原因：
- 基于假设位移场的 FEM 结构通常比实际结构**更刚**
- 非协调元允许单元间分离或重叠，使 FEM 结构**变柔**
- 这两种相反的效应可以相互抵消，有时反而得到更好的结果

> 💡 **理解关键**：非协调元的"反直觉优势"是考试简答题的经典素材。记住一句话：**位移元偏刚，非协调元引入"软"成分，两者抵消可能更好**。但这不是理所当然的——非协调元必须通过分片试验才能保证收敛，否则加密网格也不收敛。

---

### 5.5.4 位移解的精度估计

#### 误差阶

将单元位移场作 Taylor 展开：
$$\mathbf{u} = \mathbf{u}_i + \left(\frac{\partial\mathbf{u}}{\partial x}\right)_i \Delta x + \left(\frac{\partial\mathbf{u}}{\partial y}\right)_i \Delta y + \cdots$$

设单元尺寸为 $h$，则 $\Delta x, \Delta y$ 的量级为 $O(h)$。

若位移函数采用 $p$ 阶完备多项式（能精确逼近 Taylor 级数前 $p$ 阶），则：

$$\boxed{\|\mathbf{u} - \mathbf{u}^h\| \leq C h^{p+1}}$$

即位移误差为 $O(h^{p+1})$ 阶。


对 CST（$p=1$）：位移误差 $\sim O(h^2)$，收敛速率也为 $O(h^2)$。

#### 应力和应变精度

应变基于位移的 $m$ 阶导数（弹性力学 $m=1$），误差为：

$$\boxed{\|\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^h\| \leq C h^{p-m+1}}$$

对 CST（$p=m=1$）：**应变误差** $\sim O(h)$。

这一结论揭示了有限元解的"不平衡"特性：**位移精度高于应变（应力）精度**，两者差一个量级。

> ❌ **易错点**：很多同学会误以为"位移算准了，应力就也准了"——大错特错。应力是通过对位移求导得到的（$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$），求导操作会放大误差，令精度减一阶。这意味着即使位移看起来收敛得很好，应力可能还很粗糙。工程上通常以应力收敛作为网格加密的判断标准，而不是位移。

#### 应变能精度

应变能由应变的二次项表达，误差为：

$$U \text{ 误差} \sim O(h^{2(p-m+1)})$$

对 CST：应变能误差 $\sim O(h^2)$。

#### 外推提高精度

协调元在 $h\to 0$ 时单调收敛。利用两次不同网格的结果可外推：

$$\frac{u_i^{(1)} - u_i}{u_i^{(2)} - u_i} \approx \frac{O(h^2)}{O((h/2)^2)} = 4$$

得外推公式：

$$\boxed{u_i \approx \frac{1}{3}(4u_i^{(2)} - u_i^{(1)})}$$

同规格外推可提高应力精度。

---

### 5.5.5 位移元解的下限性（Lower Bound Property）

**这是位移元（基于最小势能原理的有限元）最重要的性质之一。**


#### 推理过程

离散系统的总势能：
$$\tilde{\Pi} = \frac{1}{2}\tilde{\mathbf{a}}^T\tilde{\mathbf{K}}\tilde{\mathbf{a}} - \tilde{\mathbf{a}}^T\tilde{\mathbf{P}}$$

由 $\delta\tilde{\Pi} = 0$ 得 $\tilde{\mathbf{K}}\tilde{\mathbf{a}} = \tilde{\mathbf{P}}$，代入：

$$\tilde{\Pi} = \frac{1}{2}\tilde{\mathbf{a}}^T\tilde{\mathbf{K}}\tilde{\mathbf{a}} - \tilde{\mathbf{a}}^T\tilde{\mathbf{K}}\tilde{\mathbf{a}} = -\frac{1}{2}\tilde{\mathbf{a}}^T\tilde{\mathbf{K}}\tilde{\mathbf{a}} = -\tilde{U}$$

总势能 = 负的应变能（平衡态）。

由于最小值原理：$\tilde{\Pi} \geq \Pi_{\min}$（近似解的势能大于等于精确解的势能）

因此：$\tilde{U} \leq U$（近似解的应变能**小于**精确解的应变能）

由于 $\tilde{U} = \frac{1}{2}\tilde{\mathbf{a}}^T\tilde{\mathbf{P}} = \frac{1}{2}\tilde{\mathbf{P}}^T\tilde{\mathbf{a}}$，$U = \frac{1}{2}\mathbf{P}^T\mathbf{a}$

得：

$$\boxed{\tilde{\mathbf{P}}^T\tilde{\mathbf{a}} \leq \mathbf{P}^T\mathbf{a}}$$

即近似位移解**整体上不大于**精确解——位移有限元解具有**下限性**。

#### 物理解释

连续体原本具有**无限自由度**。假设单元位移函数后，自由度退化为有限个（仅由节点位移表达）。位移函数**约束和限制了单元的变形**，使得：

- 单元刚度 > 实际连续体在对应区域的刚度
- 离散总刚 > 实际连续体的刚度
- 因此在同样荷载下，计算位移 < 实际位移（整体上偏小）

这就是"位移元偏刚、位移偏小"的直观理解。

> 💡 **理解关键**：下限性的价值不只在理论分析，更有工程实践意义——如果你用位移元算出来的位移已经满足刚度设计要求，那真实位移更大，说明真实结构更柔、你的设计是**偏保守**的。这在安全关键的结构设计中非常重要。（但如果需要的是位移的上界估计，位移元就无能为力了——此时需要应力元或混合元。）

---

### 5.5.6 FEM 求解步骤总结

弹性静力问题的主要环节：

1. **建立 FEM 理论基础、选择单元类型**——最常用的是基于最小势能原理的**位移元**
2. **建立 FEM 模型**——选择单元形式、划分网格、设定边界条件
3. **建立特征矩阵和求解方程**：
   - 构造插值函数矩阵 $[N]$
   - 建立应变矩阵 $[B] = [\partial][N]$
   - 建立单元刚度矩阵 $\mathbf{K}^e = \int_{\Omega_e} \mathbf{B}^T\mathbf{D}\mathbf{B}\,dV$
   - 计算等效节点荷载
   - 组装总刚和总荷得 $\mathbf{K}\mathbf{a} = \mathbf{P}$
4. **选择数值方法求解**——利用对称性、稀疏性、带状分布提高效率
5. **收敛性检验**——保证完备性（含刚体位移 + 常应变模态）和协调性

后续章节将讨论如何建立更高精度的单元（高阶插值、等参元），以及相应的数值积分技术。

---

## 检查你的理解（5.5）

1. 为什么单元刚度矩阵在无约束时是奇异的？从物理角度解释。
2. "完备性"和"协调性"各自的物理内涵是什么？CST 单元同时满足两者吗？
3. 分片试验为什么可以作为收敛的充分条件？非协调元为何有时比协调元精度更高？
4. CST 单元的位移误差是 $O(h^2)$ 而应力误差是 $O(h)$——这说明了什么？实际工程中应如何应对？
5. 位移元下限性的本质是什么？如果某种情况我们需要的正好是位移的下界估计，这是优势还是劣势？

---

## 核心公式速查

| 公式 | 表达式 | 用途 |
|------|--------|------|
| 形函数 | $N_s = \frac{1}{2\Delta_e}(a_s + b_s x + c_s y)$ | CST 位移插值 |
| 系数 | $b_i = y_j - y_m,\; c_i = x_m - x_j$ | 几何描述 |
| 应变矩阵 | $[B]_{(3\times 6)} = \frac{1}{2A}\begin{pmatrix} b_1 & 0 & b_2 & 0 & b_3 & 0 \\ 0 & c_1 & 0 & c_2 & 0 & c_3 \\ c_1 & b_1 & c_2 & b_2 & c_3 & b_3 \end{pmatrix}$ | 应变-位移关系 |
| 单元刚度（CST 平面应力） | $k_{rs} = \frac{Et}{4(1-\mu^2)A}\begin{pmatrix} b_r b_s + \frac{1-\mu}{2}c_r c_s & \mu b_r c_s + \frac{1-\mu}{2}c_r b_s \\ \mu c_r b_s + \frac{1-\mu}{2}b_r c_s & c_r c_s + \frac{1-\mu}{2}b_r b_s \end{pmatrix}$ | 刚度系数 |
| 总刚组装 | $[K] = \sum_{e}[k]_e$ | 总体集成 |
| FEM 方程 | $[K]\{\delta\} = \{F\}$ | 求解 |
| 位移误差 | $\|\mathbf{u}-\mathbf{u}^h\| \leq C h^{p+1}$ | 精度估计 |
| 应力误差 | $\|\boldsymbol{\sigma}-\boldsymbol{\sigma}^h\| \leq C h^p$ | 精度估计 |
| 位移下限性 | $\tilde{\mathbf{P}}^T\tilde{\mathbf{a}} \leq \mathbf{P}^T\mathbf{a}$ | 解的性质 |

---

> **对应作业**：[HW3 Q4（梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
