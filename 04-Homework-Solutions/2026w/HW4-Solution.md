# Homework 4 参考答案

> **题目**：[HW4-Problem](./HW4-Problem.md) · **提交截止**：June 22, 2026
> **对应知识**：[1-5 FEM 公式推导 §5.1-5.2](../../01-Lecture-Notes/1-5-FEM-Formulation.md)、[1-6 单元构造 §6.4](../../01-Lecture-Notes/1-6-Element-Construction.md)

---

## Problem 1 — 杆单元刚度矩阵

### 问题重述

等截面弹性直杆，两端固定，面积 $A$ 常数，总长 $3L$，受均匀体力 $f$。划分为 3 个线性杆单元（每段长 $L$）。

**节点编号**：$1-2-3-4$，**单元编号**：$(1)-(2)-(3)$

---

### 1.1 单元分析（单元 e）

对一维线性杆单元，取局部坐标 $x \in [0, L]$，节点 $i$（左）、$j$（右）。

#### 形函数

线性插值：

$$u^{(e)}(x) = N_i(x) u_i + N_j(x) u_j$$

其中：

$$N_i(x) = 1 - \frac{x}{L},\qquad N_j(x) = \frac{x}{L}$$

或写为矩阵形式：

$$\{u^{(e)}\} = [N]\{\delta\}_e = \begin{pmatrix} N_i & N_j \end{pmatrix} \begin{pmatrix} u_i \\ u_j \end{pmatrix}$$

#### 应变-位移矩阵 [B]

$$\varepsilon = \frac{du}{dx} = \begin{pmatrix} \frac{dN_i}{dx} & \frac{dN_j}{dx} \end{pmatrix} \begin{pmatrix} u_i \\ u_j \end{pmatrix} = [B]\{\delta\}_e$$

$$[B] = \begin{pmatrix} -\frac{1}{L} & \frac{1}{L} \end{pmatrix}$$

#### 单元刚度矩阵

对于一维杆单元，弹性矩阵 $[D] = E$（杨氏模量），因此：

$$[k]_e = \int_0^L [B]^T E [B]\,A\,dx = A \int_0^L \begin{pmatrix} -\frac{1}{L} \\ \frac{1}{L} \end{pmatrix} E \begin{pmatrix} -\frac{1}{L} & \frac{1}{L} \end{pmatrix} dx$$

$$[k]_e = \frac{AE}{L} \int_0^L \begin{pmatrix} \frac{1}{L} & -\frac{1}{L} \\ -\frac{1}{L} & \frac{1}{L} \end{pmatrix} dx = \frac{AE}{L} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$$

所以每个单元的刚度矩阵（局部坐标）为：

$$\boxed{[k]_e = \frac{AE}{L} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}}$$

#### 单元荷载向量（体力 $f$）

$$[F]_e = \int_0^L [N]^T f\,A\,dx = Af \int_0^L \begin{pmatrix} 1-\frac{x}{L} \\ \frac{x}{L} \end{pmatrix} dx = Af \begin{pmatrix} \frac{L}{2} \\ \frac{L}{2} \end{pmatrix}$$

即每个单元的体力等效节点力各为 $\frac{AfL}{2}$。

---

### 1.2 总体刚度矩阵组装

3 个单元，4 个节点，总体刚度矩阵 $[K]$ 为 $4\times4$。

按直接刚度法组装：

- 单元 (1)：节点 $1,2$ → 叠加到 $[K]$ 的 $(1,1),(1,2),(2,1),(2,2)$
- 单元 (2)：节点 $2,3$ → 叠加到 $(2,2),(2,3),(3,2),(3,3)$
- 单元 (3)：节点 $3,4$ → 叠加到 $(3,3),(3,4),(4,3),(4,4)$

$$[K] = \frac{AE}{L} \left[ \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} + \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} + \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & -1 & 1 \end{pmatrix} \right]$$

$$\boxed{[K] = \frac{AE}{L} \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 1 \end{pmatrix}}$$

### 1.3 引入边界条件并求解

两端固定：$u_1 = 0,\; u_4 = 0$。约束后方程（自由度 $u_2, u_3$）：

$$\frac{AE}{L} \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix} \begin{pmatrix} u_2 \\ u_3 \end{pmatrix} = AfL \begin{pmatrix} 1 \\ 1 \end{pmatrix}$$

解之：

$$\begin{pmatrix} u_2 \\ u_3 \end{pmatrix} = \frac{fL^2}{AE} \cdot \frac{1}{3} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix} \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{fL^2}{AE} \begin{pmatrix} 1 \\ 1 \end{pmatrix}$$

$$\boxed{u_2 = u_3 = \frac{fL^2}{AE}}$$

由对称性可知中点位移相等，结果合理。

---

## Problem 2 — FEM 求解 ODE

### 问题重述

求解：

$$-u''(x) + u(x) = -x,\quad x\in[0,1]$$
$$u(0) = u(1) = 0$$

区间等分为 4 个单元（5 个节点），每个单元长度 $h = \frac14$。

弱形式（已给出）：

$$Q[u(x)] = \frac12 \int_0^1 \bigl[ u'^2 + u^2 + 2x u \bigr] dx$$

---

### 2.1 单元分析

对单元 $e$（节点 $i, j$，长度 $h$），取线性形函数（同 Problem 1）：

$$N_i(\xi) = 1-\xi,\quad N_j(\xi) = \xi,\quad \xi = \frac{x-x_i}{h}\in[0,1]$$

#### 单元刚度矩阵

$$[k]_e = \int_0^h \left( [B]^T[B] + [N]^T[N] \right) dx$$

第一项（对应 $-u''$ 项）：

$$[B] = \frac{1}{h}\begin{pmatrix} -1 & 1 \end{pmatrix}$$

$$\int_0^h [B]^T[B]\,dx = \int_0^h \frac{1}{h^2} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} dx = \frac{1}{h} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$$

第二项（对应 $+u$ 项）：

$$\int_0^h [N]^T[N]\,dx = \int_0^h \begin{pmatrix} (1-\xi)^2 & \xi(1-\xi) \\ \xi(1-\xi) & \xi^2 \end{pmatrix} h\,d\xi = \frac{h}{6} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

因此单元刚度矩阵：

$$[k]_e = \frac{1}{h} \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} + \frac{h}{6} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

代入 $h = \frac14$：

$$[k]_e = 4 \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} + \frac{1}{24} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 4 + \frac{2}{24} & -4 + \frac{1}{24} \\ -4 + \frac{1}{24} & 4 + \frac{2}{24} \end{pmatrix}$$

$$\boxed{[k]_e = \begin{pmatrix} \frac{49}{12} & -\frac{95}{24} \\ -\frac{95}{24} & \frac{49}{12} \end{pmatrix}}$$

#### 单元荷载向量

$$\{F\}_e = \int_0^h [N]^T (-x)\,dx$$

在单元 $e$ 上的 $x$ 坐标：$x = x_i + \xi h$

$$\{F\}_e = -\int_0^1 \begin{pmatrix} 1-\xi \\ \xi \end{pmatrix} (x_i + \xi h)\, h\,d\xi$$

$$= -h \begin{pmatrix} x_i\int_0^1(1-\xi)d\xi + h\int_0^1\xi(1-\xi)d\xi \\ x_i\int_0^1\xi d\xi + h\int_0^1\xi^2 d\xi \end{pmatrix}$$

$$= -h \begin{pmatrix} \frac{x_i}{2} + \frac{h}{6} \\ \frac{x_i}{2} + \frac{h}{3} \end{pmatrix}$$

代入 $h = \frac14$：

$$\{F\}_e = -\frac14 \begin{pmatrix} \frac{x_i}{2} + \frac{1}{24} \\ \frac{x_i}{2} + \frac{1}{12} \end{pmatrix}$$

---

### 2.2 总体组装

各单元节点范围：

| 单元 | 节点 $i$ | 节点 $j$ | $x_i$ |
|------|---------|---------|-------|
| 1 | 1 | 2 | 0 |
| 2 | 2 | 3 | 1/4 |
| 3 | 3 | 4 | 1/2 |
| 4 | 4 | 5 | 3/4 |

**总体刚度矩阵** $[K]$ （$5\times5$，边界条件 $u_1=u_5=0$ 待引入）：

$$[K] = \frac{1}{h} \begin{pmatrix}
1 & -1 & 0 & 0 & 0 \\
-1 & 2 & -1 & 0 & 0 \\
0 & -1 & 2 & -1 & 0 \\
0 & 0 & -1 & 2 & -1 \\
0 & 0 & 0 & -1 & 1
\end{pmatrix} + \frac{h}{6} \begin{pmatrix}
2 & 1 & 0 & 0 & 0 \\
1 & 4 & 1 & 0 & 0 \\
0 & 1 & 4 & 1 & 0 \\
0 & 0 & 1 & 4 & 1 \\
0 & 0 & 0 & 1 & 2
\end{pmatrix}$$

代入 $h = \frac14$：

$$[K] = \begin{pmatrix}
4 & -4 & 0 & 0 & 0 \\
-4 & 8 & -4 & 0 & 0 \\
0 & -4 & 8 & -4 & 0 \\
0 & 0 & -4 & 8 & -4 \\
0 & 0 & 0 & -4 & 4
\end{pmatrix} + \frac{1}{24} \begin{pmatrix}
2 & 1 & 0 & 0 & 0 \\
1 & 4 & 1 & 0 & 0 \\
0 & 1 & 4 & 1 & 0 \\
0 & 0 & 1 & 4 & 1 \\
0 & 0 & 0 & 1 & 2
\end{pmatrix}$$

$$[K] = \begin{pmatrix}
\frac{49}{12} & -\frac{95}{24} & 0 & 0 & 0 \\[2pt]
-\frac{95}{24} & \frac{97}{12} & -\frac{95}{24} & 0 & 0 \\[2pt]
0 & -\frac{95}{24} & \frac{97}{12} & -\frac{95}{24} & 0 \\[2pt]
0 & 0 & -\frac{95}{24} & \frac{97}{12} & -\frac{95}{24} \\[2pt]
0 & 0 & 0 & -\frac{95}{24} & \frac{49}{12}
\end{pmatrix}$$

**总体荷载向量** $\{F\}$：

$$\{F\} = \begin{pmatrix}
F_1^{(1)} \\
F_2^{(1)} + F_2^{(2)} \\
F_3^{(2)} + F_3^{(3)} \\
F_4^{(3)} + F_4^{(4)} \\
F_5^{(4)}
\end{pmatrix}$$

各单元节点力分量：

| 单元 | $x_i$ | $F_i^{(e)}$ | $F_j^{(e)}$ |
|------|-------|-------------|-------------|
| 1 | 0 | $-1/96$ | $-1/48$ |
| 2 | 1/4 | $-13/384$ | $-7/192$ |
| 3 | 1/2 | $-25/384$ | $-11/192$ |
| 4 | 3/4 | $-37/384$ | $-5/48$ |

计算：$F_1 = -1/96 \approx -0.0104$，$F_2 = -1/48-13/384 = -21/384 \approx -0.0547$，$F_3 = -7/192-25/384 = -39/384 \approx -0.1016$，$F_4 = -11/192-37/384 = -59/384 \approx -0.1536$，$F_5 = -5/48 \approx -0.1042$

### 2.3 引入边界条件

$u_1 = u_5 = 0$，删除第 1 行第 5 行和第 1 列第 5 列，得缩减方程：

$$\begin{pmatrix}
\frac{97}{12} & -\frac{95}{24} & 0 \\[2pt]
-\frac{95}{24} & \frac{97}{12} & -\frac{95}{24} \\[2pt]
0 & -\frac{95}{24} & \frac{97}{12}
\end{pmatrix}
\begin{pmatrix} u_2 \\ u_3 \\ u_4 \end{pmatrix}
= -\begin{pmatrix} 21/384 \\ 39/384 \\ 59/384 \end{pmatrix}$$

### 2.4 数值求解

首先将矩阵归一化。令：

$$A = \frac{97}{12} \approx 8.0833,\quad B = \frac{95}{24} \approx 3.9583$$

方程组：

$$\begin{cases}
A u_2 - B u_3 = -21/384 \\
-B u_2 + A u_3 - B u_4 = -39/384 \\
- B u_3 + A u_4 = -59/384
\end{cases}$$

由对称性 $u(0)=u(1)$，此处右端不对称 → 不对称荷载。用数值求解：

$$\begin{pmatrix}
8.0833 & -3.9583 & 0 \\
-3.9583 & 8.0833 & -3.9583 \\
0 & -3.9583 & 8.0833
\end{pmatrix}
\begin{pmatrix} u_2 \\ u_3 \\ u_4 \end{pmatrix}
= \begin{pmatrix} -0.05469 \\ -0.10156 \\ -0.15365 \end{pmatrix}$$

解得：

$$u_2 \approx -0.0187,\quad u_3 \approx -0.0276,\quad u_4 \approx -0.0298$$

即各节点 FEM 近似解：

$$\boxed{u(0)=0,\; u(0.25)\approx -0.0187,\; u(0.5)\approx -0.0276,\; u(0.75)\approx -0.0298,\; u(1)=0}$$

---

## Problem 3 — 划线法构造形函数

### 问题重述

使用**划线法**构造 4 节点四边形单元的形函数。节点在 $(\xi,\eta)$ 局部坐标系中的坐标：

$$1(1,1),\quad 2(-1,1),\quad 3(-1,-1),\quad 4(1,-1)$$

所有边长相等（正方形），这是标准 Q4 等参单元。

### 3.1 划线法原理

划线法：对节点 $i$，画出**通过除 $i$ 外所有其他节点**的直线方程，将这些直线方程的左端表达式相乘，再归一化（使 $N_i=1$ 在节点 $i$ 处成立）。

### 3.2 各节点形函数构造

#### 节点 1 $(1,1)$

除节点 1 外的节点为 $2(-1,1)$, $3(-1,-1)$, $4(1,-1)$。

构造步骤：

1. **过点 2 和 3 的直线**：$\xi = -1$ → 直线方程 $\boxed{\xi+1=0}$
2. **过点 3 和 4 的直线**：$\eta = -1$ → 直线方程 $\boxed{\eta+1=0}$

> 注意：划线法要求通过**除目标节点外所有节点**的**最小直线集合**。这里 $\xi+1=0$ 经过节点 2 和 3，$\eta+1=0$ 经过节点 3 和 4，两条直线一起覆盖了节点 2,3,4。

乘积极形式：$(\xi+1)(\eta+1)$

归一化：在节点 1 $(1,1)$ 处，$(\xi+1)(\eta+1) = (2)(2) = 4$

$$N_1 = \frac{(\xi+1)(\eta+1)}{4} = \frac14 (1+\xi)(1+\eta)$$

#### 节点 2 $(-1,1)$

除节点 2 外的节点为 $1(1,1)$, $3(-1,-1)$, $4(1,-1)$。

1. **过点 1 和 4 的直线**：$\eta = -1$ → $\boxed{\eta+1=0}$  
   但这不过点 3，还需另一条
2. **过点 1 的另一种选取**：过点 1 和 4 的直线为 $\xi$ 方向 → 不对
   
重新审视：划线法的标准做法是找出过**除目标节点外所有节点**的直线。

- 过节点 $1(1,1)$ 和 $3(-1,-1)$ 的直线：$\xi-\eta=0$。但还需要覆盖节点 $4(1,-1)$。
- 过节点 $3(-1,-1)$ 和 $4(1,-1)$ 的直线：$\eta=-1$ → $\eta+1=0$
- 过节点 $1(1,1)$ 和 $4(1,-1)$ 的直线：$\xi=1$ → $\xi-1=0$

取 $\xi-1=0$（通过节点 1,4）和 $\eta-1=0$（通过节点 1,3）→ 不对。

实际上，划线法的核心是：对每个节点 $i$，找出所有不通过节点 $i$ 的边的方程，它们的乘积构成形函数的分子。

**更直观的理解**：对 Q4 单元，形函数就是两个方向线性插值的乘积：

$$N_i = \frac14 (1+\xi_i\xi)(1+\eta_i\eta)$$

不需要每步划线。但题目要求用划线法推导。

划线法步骤：

对节点 2 $(-1,1)$：
- 边 3-4：$\eta=-1$，方程 $\eta+1=0$ → 不通过节点 2 ✓
- 边 1-4：$\xi=1$，方程 $\xi-1=0$ → 不通过节点 2 ✓

乘积：$(\eta+1)(\xi-1)$

归一化：在节点 2 $(-1,1)$ 处，$(\eta+1)(\xi-1) = (2)(-2) = -4$

$$N_2 = \frac{(\eta+1)(\xi-1)}{-4} = \frac14 (1-\xi)(1+\eta)$$

#### 节点 3 $(-1,-1)$

- 边 1-2：$\eta=1$，方程 $\eta-1=0$ → 不通过节点 3 ✓
- 边 1-4：$\xi=1$，方程 $\xi-1=0$ → 不通过节点 3 ✓

乘积：$(\eta-1)(\xi-1)$

归一化：在节点 3 $(-1,-1)$ 处，$(\eta-1)(\xi-1) = (-2)(-2) = 4$

$$N_3 = \frac{(\eta-1)(\xi-1)}{4} = \frac14 (1-\xi)(1-\eta)$$

#### 节点 4 $(1,-1)$

- 边 1-2：$\eta=1$，方程 $\eta-1=0$ → 不通过节点 4 ✓
- 边 2-3：$\xi=-1$，方程 $\xi+1=0$ → 不通过节点 4 ✓

乘积：$(\eta-1)(\xi+1)$

归一化：在节点 4 $(1,-1)$ 处，$(\eta-1)(\xi+1) = (-2)(2) = -4$

$$N_4 = \frac{(\eta-1)(\xi+1)}{-4} = \frac14 (1+\xi)(1-\eta)$$

### 3.3 结果汇总

$$\boxed{\begin{aligned}
N_1 &= \frac14 (1+\xi)(1+\eta) \\
N_2 &= \frac14 (1-\xi)(1+\eta) \\
N_3 &= \frac14 (1-\xi)(1-\eta) \\
N_4 &= \frac14 (1+\xi)(1-\eta)
\end{aligned}}$$

### 3.4 性质验证

**性质 1**：$\sum_{i=1}^4 N_i = 1$（刚体位移）

$$\frac14[(1+\xi)(1+\eta) + (1-\xi)(1+\eta) + (1-\xi)(1-\eta) + (1+\xi)(1-\eta)]$$
$$= \frac14[2(1+\eta) + 2(1-\eta)] = \frac14[4] = 1 \quad \checkmark$$

**性质 2**：$N_i(\xi_j,\eta_j) = \delta_{ij}$（插值性）

以 $N_1$ 为例：
- 在节点 1：$N_1(1,1)= \frac14(2)(2)=1$
- 在节点 2：$N_1(-1,1)= \frac14(0)(2)=0$
- 在节点 3：$N_1(-1,-1)= \frac14(0)(0)=0$
- 在节点 4：$N_1(1,-1)= \frac14(2)(0)=0$

其余节点同理。 $\checkmark$

### 3.5 等参变换说明

这 4 个形函数就是标准 Q4 等参单元的形函数。坐标变换同样使用它们：

$$x = \sum_{i=1}^4 N_i(\xi,\eta)\, x_i,\qquad y = \sum_{i=1}^4 N_i(\xi,\eta)\, y_i$$

从而将标准正方形 $(\xi,\eta)\in[-1,1]^2$ 映射到任意四边形。
