# 第6章：单元构造与形函数

> **对应 PDF**：[`6 FEM_Element construction.pdf`](../06-References/pdfs-originals/6%20FEM_Element%20construction.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §5
> **相关作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 5 章、线性代数（插值理论）、数值分析

---

## 6.1 形函数的作用

### 6.1.1 形函数是什么？

形函数（Shape Function）是 FEM 中连接**节点值**和**单元内任意点值**的桥梁：
$$u^{(e)}(x) = \sum_{i=1}^n N_i(x)\,u_i^{(e)}$$

其中 $u_i^{(e)}$ 是待求的节点位移，$N_i(x)$ 是已知的插值函数。

### 6.1.2 为什么形函数如此重要？

1. **标准化**：形函数 → 单元刚度矩阵 → 总装 → 求解，整个过程可以完全程序化
2. **收敛性**：解的精度和收敛行为直接由形函数的性质决定
3. **计算效率**：形函数的阶数决定了单元的自由度数，直接影响计算量

### 6.1.3 构造单元前需确定的三个因素

| 因素 | 选择 | 示例 |
|------|------|------|
| **几何形状** | 1D/2D/3D | 线、三角形、四边形、四面体 |
| **节点分布** | 顶点 / 边中点 / 内部点 | 3 节点三角、6 节点三角 |
| **节点自由度** | Lagrange / Hermite | 仅位移 / 位移+转角 |

---

## 6.2 插值基础

### 6.2.1 Lagrange 插值

已知 $n$ 个点的函数值，求 $n-1$ 次多项式：
$$P_{n-1}(x) = \sum_{i=1}^n f(x_i)L_i(x),\quad L_i(x) = \prod_{j\neq i}\frac{x-x_j}{x_i-x_j}$$

$L_i(x_j) = \delta_{ij}$ ✅

### 6.2.2 Hermite 插值

不仅函数值，导数值也已知。梁弯曲问题（4 个条件 → 三次多项式）是典型应用。

### 6.2.3 广义 Lagrange 公式

通用形函数构造：
$$N_i = \prod_{j\neq i}\frac{f_j(\xi)}{f_j(\xi_i)},\quad f_j(\xi) = \xi - \xi_j$$

自动满足 $N_i(\xi_j) = \delta_{ij}$。

---

## 6.3 一维单元与长度坐标

### 6.3.1 长度坐标

单元 $e_i = [x_i, x_{i+1}]$，长度 $L_i$。定义：
$$\lambda_1 = \frac{x_{i+1} - x}{L_i},\quad \lambda_2 = \frac{x - x_i}{L_i},\quad \lambda_1 + \lambda_2 = 1$$

**长度坐标下的积分公式**（核心工具）：
$$\int_{x_i}^{x_{i+1}} \lambda_1^{\alpha_1}\lambda_2^{\alpha_2}\,dx = L_i\int_0^1 \lambda_1^{\alpha_1}(1-\lambda_1)^{\alpha_2}d\lambda_1 = \frac{\alpha_1!\,\alpha_2!}{(\alpha_1+\alpha_2+1)!}L_i$$

这个公式使单元刚度矩阵的积分计算变成简单的代数运算——代指数即可。

### 6.3.2 线性 Lagrange 单元（2 节点）

$$N_1 = \frac{1-\xi}{2},\quad N_2 = \frac{1+\xi}{2},\quad \xi = \frac{2x}{l}\in[-1,1]$$

验证：$N_1(-1)=1, N_1(1)=0$ ✅，$N_1+N_2=1$ ✅

### 6.3.3 二次 Lagrange 单元（3 节点）

$$N_1 = \frac{\xi(\xi-1)}{2},\; N_2 = \frac{\xi(\xi+1)}{2},\; N_3 = 1-\xi^2$$

### 6.3.4 Hermite 三次单元（梁单元）

每节点 2 个自由度（挠度 $w$ + 转角 $\theta = dw/dx$），共 4 个条件 → 三次多项式。

设 $w(\xi) = c_0 + c_1\xi + c_2\xi^2 + c_3\xi^3$，$\xi = x/L$：

节点条件：
$$w(0)=w_1,\; \frac{1}{L}\frac{dw}{d\xi}\big|_0 = \theta_1,\; w(1)=w_2,\; \frac{1}{L}\frac{dw}{d\xi}\big|_1 = \theta_2$$

求解得标准 Hermite 形函数：
$$\boxed{\begin{aligned}
N_1 &= 1 - 3\xi^2 + 2\xi^3 \quad &\text{(对应 } w_1 \text{)}\\
N_2 &= L(\xi - 2\xi^2 + \xi^3) \quad &\text{(对应 } \theta_1 \text{)}\\
N_3 &= 3\xi^2 - 2\xi^3 \quad &\text{(对应 } w_2 \text{)}\\
N_4 &= L(-\xi^2 + \xi^3) \quad &\text{(对应 } \theta_2 \text{)}
\end{aligned}}$$

**性质验证**：$N_1+N_3=1$（刚体平移），$N_2+N_4 = L(1-\xi)$（刚体转动）。

由最小势能原理得梁单元刚度矩阵：
$$[k]_e = \frac{EI}{L^3}\begin{pmatrix}
12 & 6L & -12 & 6L \\
6L & 4L^2 & -6L & 2L^2 \\
-12 & -6L & 12 & -6L \\
6L & 2L^2 & -6L & 4L^2
\end{pmatrix}$$

---

## 6.4 二维单元

### 6.4.1 面积坐标（三角形单元的"自然坐标"）

对三角形 $\triangle Q_1Q_2Q_3$，定义：
$$L_1 = \frac{A_1}{A},\quad L_2 = \frac{A_2}{A},\quad L_3 = \frac{A_3}{A}$$

其中 $A$ 是总面积，$A_i$ 是点 $P$ 与对边 $Q_jQ_k$ 围成的子三角形面积。

**性质**：
1. $L_1+L_2+L_3=1$（只有 2 个独立坐标）
2. 顶点坐标：$(1,0,0), (0,1,0), (0,0,1)$；形心 $(\frac13,\frac13,\frac13)$
3. 三边方程：$L_1=0$（$Q_2Q_3$），$L_2=0$（$Q_1Q_3$），$L_3=0$（$Q_1Q_2$）
4. 平行于对边的直线上 $L_i$ 为常数

**与直角坐标的互化**：
$$\begin{pmatrix}x\\y\\1\end{pmatrix} = \begin{pmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
1 & 1 & 1
\end{pmatrix}\begin{pmatrix}L_1\\L_2\\L_3\end{pmatrix}$$

**面积坐标下的积分公式（核心！）**：
$$\boxed{\iint_{\Delta_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}\,dxdy = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!}\,2\Delta_e}$$

**例**：$\iint L_1L_2\,dxdy = \frac{1!\,1!\,0!}{(1+1+0+2)!}2\Delta = \frac{2\Delta}{24} = \frac{\Delta}{12}$

**导数变换**（面积坐标→直角坐标）：
$$\begin{cases}
\frac{\partial}{\partial x} = \frac{1}{2\Delta_e}\left(b_1\frac{\partial}{\partial L_1} + b_2\frac{\partial}{\partial L_2} + b_3\frac{\partial}{\partial L_3}\right) \\[4pt]
\frac{\partial}{\partial y} = \frac{1}{2\Delta_e}\left(c_1\frac{\partial}{\partial L_1} + c_2\frac{\partial}{\partial L_2} + c_3\frac{\partial}{\partial L_3}\right)
\end{cases}$$

### 6.4.2 CST（3 节点三角形）

$$N_i = L_i = \frac{1}{2\Delta}(a_i + b_i x + c_i y)$$

线性位移场 → 常应变。$[B]$ 常数矩阵。

### 6.4.3 划线法构造 LST（6 节点三角形）

"划线法"：对节点 $i$，找出除 $i$ 外所有经过其他节点的直线，将每条直线的方程左端相乘，除以在节点 $i$ 处的值归一化。

**例**：对节点 1（顶点），经过节点 2,5,3 的直线为 $L_1=0$，经过节点 4,6 的直线为 $L_1 - \frac12 = 0$。
$$N_1 = \frac{L_1 - \frac12}{\frac12}\cdot\frac{L_1}{1} = L_1(2L_1-1)$$

对节点 4（边中点），经过节点 1,6 的直线为 $L_2=0$，经过节点 2,5 的直线为 $L_3=0$。
$$N_4 = \frac{L_2}{1/2}\cdot\frac{L_3}{1/2} = 4L_1L_2$$

完整组：
$$\begin{aligned}
N_1 &= L_1(2L_1-1), & N_2 &= L_2(2L_2-1), & N_3 &= L_3(2L_3-1) \\
N_4 &= 4L_1L_2, & N_5 &= 4L_2L_3, & N_6 &= 4L_3L_1
\end{aligned}$$

### 6.4.4 Q4（4 节点矩形）

自然坐标 $(\xi,\eta)\in[-1,1]$：
$$N_i = \frac14(1+\xi_i\xi)(1+\eta_i\eta)$$

双线性插值。$[B]$ 线性变化 → 精度高于 CST。

### 6.4.5 Serendipity 单元（Q8）

节点集中在边界上：4 角点 + 4 边中点。比 9 节点 Lagrange 单元少 1 个内部节点，精度相近。

形函数（角节点 1）：
$$N_1 = \frac14(1+\xi_1\xi)(1+\eta_1\eta)(\xi_1\xi+\eta_1\eta-1)$$

边中点 5：$N_5 = \frac12(1-\xi^2)(1+\eta_5\eta)$

修正：添加边节点后，原角节点的形函数需要减去修正项（如 $N_1 \leftarrow N_1 - N_5/2 - N_8/2$）。

### 6.4.6 Pascal 三角形

决定单元多项式的项次选择：
```
1
x    y
x²  xy  y²
x³  x²y xy² y³
```
三角形单元取对称的三角形部分，矩形单元取矩形部分。

---

## 6.5 三维单元

### 6.5.1 4 节点四面体

体积坐标 $L_i = V_i/V$。积分公式：
$$\iiint_{V_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}L_4^{\alpha_4}\,dV = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!\,\alpha_4!}{(\sum\alpha_i+3)!}\,6V_e$$

### 6.5.2 8 节点六面体

$$N_i = \frac18(1+\xi\xi_i)(1+\eta\eta_i)(1+\zeta\zeta_i)$$

---

## 6.6 等参元

### 6.6.1 核心思想

坐标和位移使用**相同的形函数**：
$$x = \sum N_i(\xi,\eta)x_i,\quad u = \sum N_i(\xi,\eta)u_i$$

### 6.6.2 Jacobian 矩阵

$$\begin{pmatrix} \partial N_i/\partial x \\ \partial N_i/\partial y \end{pmatrix} = \mathbf{J}^{-1}\begin{pmatrix} \partial N_i/\partial \xi \\ \partial N_i/\partial \eta \end{pmatrix}$$

$$\mathbf{J} = \begin{pmatrix}
\partial x/\partial\xi & \partial y/\partial\xi \\
\partial x/\partial\eta & \partial y/\partial\eta
\end{pmatrix}$$

**关键条件**：$|\mathbf{J}| \neq 0$ 在整个单元内成立。

### 6.6.3 等参元刚度矩阵

$$[k]_e = \int_{-1}^1\int_{-1}^1 [B]^T[D][B]\,t\,|\mathbf{J}|\,d\xi d\eta$$

---

## 6.7 Gauss 数值积分

$$\int_{-1}^1 f(\xi)d\xi \approx \sum_{i=1}^n w_i f(\xi_i)$$

| $n$ | $\xi_i$ | $w_i$ | 精度 |
|-----|---------|-------|------|
| 1 | 0 | 2 | 1 次 |
| 2 | $\pm 1/\sqrt{3}$ | 1 | 3 次 |
| 3 | $\pm\sqrt{0.6}, 0$ | $5/9, 8/9$ | 5 次 |

**2D**：$\iint f(\xi,\eta)d\xi d\eta \approx \sum\sum w_i w_j f(\xi_i,\eta_j)$

**减缩积分**：用低于精确所需的阶数 → 避免剪切自锁，但可能引入零能模式。

---

## 6.8 收敛准则

1. **完备性**：能表示刚体位移 + 常应变 → $\sum N_i = 1$
2. **协调性**：单元间 $C^{m-1}$ 连续 → 协调元解单调收敛
3. **分片试验**：判断新单元是否满足收敛性的数值测试

---

## 6.9 等参/超参/次参

| 类型 | 坐标节点 vs 插值节点 | 说明 |
|------|---------------------|------|
| **等参元** | 相等 | 最常用 |
| **超参元** | 坐标节点 > 插值节点 | 几何更精确 |
| **次参元** | 坐标节点 < 插值节点 | 函数更精确 |

---

## 检查你的理解

1. 形函数必须满足哪两个基本条件？分别对应什么物理意义？
2. CST、LST 和 Q4 三种单元的精度有什么差别？各自适用于什么场景？
3. 什么是面积坐标？它相对于直角坐标有什么优势？
4. 等参元的核心思想是什么？等参、超参、次参三者的区别是什么？
5. Gauss 积分 $n=2$ 能精确积分几次多项式？对于 $4\times4$ 矩形的刚度矩阵，需要多少 Gauss 点？
6. 什么是协调元和非协调元？为什么非协调元有时能得到更好的结果？

### 形函数验证练习

**练习**：验证 Hermite 梁单元形函数 $N_1 = 1-3\xi^2+2\xi^3$ 满足：
- $N_1(0)=1$，$N_1(1)=0$
- $N_1'(0)=0$，$N_1'(1)=0$

**练习**：对 LST 单元验证 $N_1 + N_2 + N_3 + N_4 + N_5 + N_6 = 1$

**练习**：用面积坐标积分公式计算 $\iint L_1^2 L_2\,dxdy$。

---

> **对应作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
