# 第6章：单元构造与形函数

> **对应 PDF**：[`6 FEM_Element construction.pdf`](../06-References/pdfs-originals/6%20FEM_Element%20construction.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §5
> **相关作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 5 章（FEM 公式）、线性代数（多项式插值）、数值分析（数值积分）

---

## 6.1 形函数的作用与重要性

### 6.1.1 形函数做什么？

在 FEM 中，形函数（Shape Function）是在**单元级别**上建立节点值与单元内任意点值之间关系的桥梁：

$$u^{(e)}(x,y) = \sum_{i=1}^{n} N_i(x,y) \, u_i^{(e)}$$

其中：
- $u_i^{(e)}$：节点 $i$ 的位移值（需要求解的未知数）
- $N_i(x,y)$：节点 $i$ 的形函数（已知函数，由单元类型决定）

### 6.1.2 为什么形函数重要？

1. **标准化**：形函数 → 单元刚度矩阵 → 总装 → 求解，整个过程可以完全程序化
2. **收敛性**：FEM 解的精度和收敛性直接由形函数的性质决定（完备性、协调性）
3. **求解效率**：形函数的阶数决定了单元的自由度数和计算量

### 6.1.3 单元构造三要素

在构造任何单元之前，需要确定三个因素：

| 因素 | 选择 | 例子 |
|------|------|------|
| **几何形状** | 1D/2D/3D 的形状 | 线、三角形、四边形、四面体 |
| **节点分布** | 边界节点 + 内部节点 | 顶点、边中点、形心 |
| **节点自由度** | Lagrange 型 / Hermite 型 | 仅位移 / 位移+转角 |

---

## 6.2 插值基本理论

### 6.2.1 Lagrange 插值

已知函数在 $n$ 个点上的取值，求一个 $n-1$ 次多项式满足这些点：

$$P_{n-1}(x) = \sum_{i=1}^n f(x_i) L_i(x)$$

其中 Lagrange 基函数：
$$L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$$

性质：$L_i(x_j) = \delta_{ij}$（在自身节点为 1，在其他节点为 0）

### 6.2.2 Hermite 插值

不仅知道函数值，还知道导数值（如梁单元的挠度 $w$ 和转角 $\theta = dw/dx$）。

### 6.2.3 广义 Lagrange 插值公式

对于任意 $n$ 节点单元的形函数，通用构造公式：
$$N_i = \prod_{j=1,\,j\neq i}^{n} \frac{f_j(\xi)}{f_j(\xi_i)}$$

其中 $f_j(\xi) = \xi - \xi_j$ 表示点到节点 $j$ 的距离。该公式自动保证 $N_i(\xi_j) = \delta_{ij}$。

---

## 6.3 一维单元

### 6.3.1 长度坐标

在 1D 单元中，引入**长度坐标** $\lambda_1, \lambda_2$ 大大简化了计算：
$$\lambda_1 = \frac{x_{i+1} - x}{L_i},\quad \lambda_2 = \frac{x - x_i}{L_i},\quad \lambda_1 + \lambda_2 = 1$$

长度坐标下的积分公式（核心！）：
$$\int_{x_i}^{x_{i+1}} \lambda_1^{\alpha_1}\lambda_2^{\alpha_2}\,dx = L_i \int_0^1 \lambda_1^{\alpha_1}(1-\lambda_1)^{\alpha_2}d\lambda_1 = \frac{\alpha_1!\,\alpha_2!}{(\alpha_1+\alpha_2+1)!}L_i$$

### 6.3.2 线性 Lagrange 单元（2 节点）

$$N_1 = \frac{1-\xi}{2},\quad N_2 = \frac{1+\xi}{2},\quad \xi \in [-1,1]$$

性质验证：
- $N_1(-1)=1, N_1(1)=0$ ✅
- $N_1 + N_2 = 1$ ✅（刚体平移）
- $-1\cdot N_1 + 1\cdot N_2 = \xi$ ✅（线性场）

### 6.3.3 二次 Lagrange 单元（3 节点）

$$N_1 = \frac{\xi(\xi-1)}{2},\quad N_2 = \frac{\xi(\xi+1)}{2},\quad N_3 = 1-\xi^2$$

### 6.3.4 Hermite 三次单元（梁单元）

每节点 2 个自由度：$w_i$（挠度）和 $\theta_i = dw/dx$（转角）。

形函数（自然坐标 $\xi = x/L$）：
$$w(\xi) = N_1 w_1 + N_2\theta_1 + N_3 w_2 + N_4\theta_2$$

其中：
$$\boxed{\begin{aligned}
N_1 &= 1 - 3\xi^2 + 2\xi^3 \\
N_2 &= L(\xi - 2\xi^2 + \xi^3) \\
N_3 &= 3\xi^2 - 2\xi^3 \\
N_4 &= L(-\xi^2 + \xi^3)
\end{aligned}}$$

**性质验证**：

| 形函数 | $\xi=0$ | $\xi=1$ | $d/d\xi$ at 0 | $d/d\xi$ at 1 |
|--------|---------|---------|---------------|---------------|
| $N_1$（$w_1$ 的系数） | 1 | 0 | 0 | 0 |
| $N_2$（$\theta_1$ 的系数） | 0 | 0 | $L$ | 0 |
| $N_3$（$w_2$ 的系数） | 0 | 1 | 0 | 0 |
| $N_4$（$\theta_2$ 的系数） | 0 | 0 | 0 | $L$ |

$$N_1 + N_3 = 1 \quad\text{（刚体平移）}$$
$$-N_2 + N_4 = L(1-\xi) \quad\text{（系数匹配保证刚体转动）}$$

---

## 6.4 二维单元

### 6.4.1 面积坐标

对于三角形单元，**面积坐标** $L_1, L_2, L_3$ 是最自然的坐标系统：
$$L_1 = \frac{A_1}{A},\quad L_2 = \frac{A_2}{A},\quad L_3 = \frac{A_3}{A}$$

其中 $A$ 为三角形总面积，$A_i$ 为对边子三角形面积。

性质：
1. $L_1 + L_2 + L_3 = 1$（只有 2 个独立坐标）
2. 顶点坐标：$(1,0,0), (0,1,0), (0,0,1)$
3. 三边方程：$L_1=0$（对面），$L_2=0$（对面），$L_3=0$（对面）

**与直角坐标互化**：
$$\begin{pmatrix} x \\ y \\ 1 \end{pmatrix} = \begin{pmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
1 & 1 & 1
\end{pmatrix} \begin{pmatrix} L_1 \\ L_2 \\ L_3 \end{pmatrix}$$

**面积坐标下的积分公式**（极其有用！）：
$$\boxed{\iint_{\Delta_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}\,dxdy = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!}\,2\Delta_e}$$

**例**：计算 $\iint L_1 L_2\,dxdy$（$\alpha_1=\alpha_2=1,\alpha_3=0$）
$$= \frac{1!\,1!\,0!}{(1+1+0+2)!}2\Delta = \frac{1}{4!}2\Delta = \frac{2\Delta}{24} = \frac{\Delta}{12}$$

这在计算单元刚度矩阵中频繁出现。

### 6.4.2 CST（3 节点三角形）

形函数：
$$N_i = L_i = \frac{1}{2\Delta}(a_i + b_i x + c_i y)\quad(i=1,2,3)$$

- 线性位移场 → 常应变
- $[B]$ 常数 → 单元内应力为常数

### 6.4.3 LST（6 节点三角形）

6 个节点：3 个顶点 + 3 条边中点。二次位移场。

形函数（用"划线法"构造）：

节点 1（顶点）：$N_1 = L_1(2L_1-1)$
节点 4（边中点）：$N_4 = 4L_1L_2$

完整组：
$$\begin{aligned}
N_1 &= L_1(2L_1-1),\quad N_2 = L_2(2L_2-1),\quad N_3 = L_3(2L_3-1) \\
N_4 &= 4L_1L_2,\quad N_5 = 4L_2L_3,\quad N_6 = 4L_3L_1
\end{aligned}$$

验证 $N_1(L_1=1) = 1\cdot(2-1)=1$ ✅，$N_4(L_1=L_2=1/2)=4\cdot\frac12\cdot\frac12=1$ ✅

### 6.4.4 Q4（4 节点矩形单元）

自然坐标 $\xi,\eta \in [-1,1]$：
$$N_i = \frac14(1+\xi_i\xi)(1+\eta_i\eta)$$

- 双线性位移场（固定一个方向时，另一个方向线性变化）
- 单元内应变线性变化 → 精度高于 CST

### 6.4.5 Pascal 三角形

决定单元多项式应包含哪些项：
```
1
x    y
x²  xy  y²
x³  x²y xy² y³
```
- 三角形单元：取对称的三角形部分
- 矩形单元：取矩形部分（如 Q4：$1, \xi, \eta, \xi\eta$）

---

## 6.5 三维单元

### 6.5.1 4 节点四面体（线性）
- 每个节点 3 个 DOF $(u,v,w)$
- 体积坐标 $L_i = V_i/V$
- 体积坐标下积分公式：
$$\iiint_{V_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}L_4^{\alpha_4}\,dV = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!\,\alpha_4!}{(\sum\alpha_i +3)!}\,6V_e$$

### 6.5.2 8 节点六面体
$$N_i = \frac18(1+\xi\xi_i)(1+\eta\eta_i)(1+\zeta\zeta_i)$$

---

## 6.6 等参元

### 6.6.1 核心思想

**坐标和位移用相同的形函数**：
$$x = \sum N_i(\xi,\eta)x_i,\quad y = \sum N_i(\xi,\eta)y_i$$
$$u = \sum N_i(\xi,\eta)u_i,\quad v = \sum N_i(\xi,\eta)v_i$$

**优势**：任意四边形/曲面单元可以映射到标准正方形/正方体 → 程序统一处理。

### 6.6.2 Jacobian 矩阵

链式法则给出坐标变换的导数关系：
$$\begin{pmatrix} \partial N_i/\partial x \\ \partial N_i/\partial y \end{pmatrix} = \mathbf{J}^{-1} \begin{pmatrix} \partial N_i/\partial \xi \\ \partial N_i/\partial \eta \end{pmatrix}$$

Jacobian 矩阵：
$$\mathbf{J} = \begin{pmatrix}
\partial x/\partial\xi & \partial y/\partial\xi \\
\partial x/\partial\eta & \partial y/\partial\eta
\end{pmatrix}$$

**重要条件**：$|\mathbf{J}| \neq 0$ 必须在整个单元内成立，否则坐标变换不是一对一的。这要求单元不能太畸变（内角不能接近 $180^\circ$，边不能过度弯曲）。

### 6.6.3 等参元刚度矩阵

在局部坐标系下数值积分：
$$[k]_e = \int_{-1}^1\int_{-1}^1 [B(\xi,\eta)]^T[D][B(\xi,\eta)]\,t\,|\mathbf{J}|\,d\xi d\eta$$

---

## 6.7 数值积分

等参元的刚度矩阵无法解析积分 → 必须使用数值积分。

### Gauss 积分

$$ \int_{-1}^1 f(\xi)\,d\xi \approx \sum_{i=1}^n w_i f(\xi_i) $$

| $n$ | $\xi_i$ | $w_i$ | 精确度 |
|-----|---------|-------|--------|
| 1 | 0 | 2 | 线性 |
| 2 | $\pm 1/\sqrt{3}$ | 1 | 三次 |
| 3 | $\pm\sqrt{0.6},\;0$ | $5/9,\;8/9$ | 五次 |

**二维**：
$$\int_{-1}^1\int_{-1}^1 f(\xi,\eta)\,d\xi d\eta \approx \sum_{i=1}^n\sum_{j=1}^n w_i w_j f(\xi_i,\eta_j)$$

**减缩积分**：使用比精确积分更少的 Gauss 点 → 避免剪切自锁，但可能引入零能模态。

---

## 6.8 Serendipity 单元

Serendipity 单元将节点集中在边界上，减少内部节点。以 8 节点四边形（Q8）为例——4 个角点 + 4 个边中点，比 9 节点 Lagrange 单元少 1 个内部节点，精度相近。

形函数（角节点 1）：
$$N_1 = \frac14(1+\xi_1\xi)(1+\eta_1\eta)(\xi_1\xi + \eta_1\eta - 1)$$

边中点 5：$$N_5 = \frac12(1-\xi^2)(1+\eta_5\eta)$$

---

## 6.9 收敛准则（必考！）

### 完备性（Completeness）
形函数必须能表示**刚体位移**和**常应变**状态。数学上：$\sum N_i = 1$（刚体平移），$\sum N_i x_i = x$（刚体转动/线性场）。

### 协调性（Compatibility）
单元间边界上位移连续（$C^0$ 连续对仅含一阶导数的泛函）。
- **协调元**（同时满足完备性和协调性）：解单调收敛
- **非协调元**（仅满足完备性）：可能更准（刚柔抵消），但不保证单调收敛

### 分片试验（Patch Test）
取几个单元拼装成小块，对边界节点施加常应变对应的位移。若内部节点计算出的应变为同一常数 → 通过分片试验。

---

## 6.10 等参/超参/次参单元对比

| 类型 | 坐标节点数 $m$ vs 插值节点数 $n$ | 说明 |
|------|--------------------------------|------|
| **等参元** | $m = n$ | 最常用，坐标和位移同阶近似 |
| **超参元** | $m > n$ | 几何描述更精确（如曲线边界） |
| **次参元** | $m < n$ | 函数逼近更精确 |

---

## 检查你的理解

1. 形函数必须满足哪两个基本性质？为什么它们重要？
2. CST 和 LST 的主要区别是什么？什么时候需要选择 LST？
3. 等参元的核心思想是什么？为什么要保证 $|\mathbf{J}| \neq 0$？
4. Gauss 积分 $n=2$ 能精确积分几次多项式？如果被积函数是 4 次多项式，应该用几个点？
5. 什么情况下非协调单元可能比协调单元更精确？

---

> **对应作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **往年相关**：[Homework3 (past)](../04-Homework-Solutions/past/HW3/Homework3.md) · [LIU Sai 答案](../04-Homework-Solutions/past/HW3/Ans%20to%20HM3_LIU%20Sai_handed%20in.md)
