# 单元构造与形函数

> **对应课件**：[`Chapter 2 Elastic theory.pdf`](../06-References/pdfs-originals/Chapter%202%20Elastic%20theory.pdf) 课程第 4 章
> **章节定位**：Construction of element and shape functions → I. Introduction → II. Two-dimensional situation → III. Three-dimensional → IV. Isoparametric element and numerical integration → V. Conclusion
> **相关作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **前置知识**：第 5 章（FEM 公式）、插值理论

---

## I. Introduction（引言）

### 1.1 形函数的重要性

单元形函数的构造是 FEM 中最重要也最技巧性的环节之一，因为：

1. **标准化与自动化**：单元刚度矩阵的计算、总体刚度矩阵的合成到有限元方程组的求解，全过程可以标准化和自动化。
2. **收敛性**：有限元空间中的基函数由形函数生成，因此形函数关系到解的**收敛性和收敛精度**。

一般来说，单元的选择取决于结构或求解域的几何特征、方程类型和预期的解决精度。

有限元的插值函数/形函数取决于：
- 单元的**形状**（1D/2D/3D）
- **节点类型**（Lagrange/Hermite）
- **节点数目**

### 构造单元需确定的因素

**① 单元的几何形状**
- 一维单元：直线或曲线
- 二维单元：三角形、矩形或任意四边形
- 三维单元：四面体、五面体、三棱柱或六面体

**② 节点个数和分布**
- 内节点：位于单元内部；外节点：位于边界上
- 为便于构造形函数，节点通常布置在特殊位置：一维单元的两端点、三角形的三个顶点、三边中点或形心

**③ 节点的自由度（DOF）**
- **Lagrange 型**：仅包含位移场，$n$ 维问题结点有 $n$ 个自由度
- **Hermite 型**：包含位移和导数（转角），$n$ 维问题结点有 $2n$ 个自由度

对应两类插值：
- **Lagrange 插值**：仅要求插值多项式在插值点处的函数值已知
- **Hermite 插值**：还要求多项式的导数（包括一阶、二阶等）已知

---

## II. Two-dimensional situation（二维情况）

二维插值是一维插值的自然推广。但随着区域维数的增加，新的特点和困难也随之出现：
1. 二维插值与区域划分方式密切相关（三角形、矩形、任意四边形）
2. 连接两个单元的不仅是节点，而是共有边——节点处的连续性不能保证整体公共边上的连续性
3. 插值点可以是顶点、边界点或内点

### 2.1 Triangular element（三角形单元）

#### 线性插值与面积坐标

三角形单元在二维问题中应用广泛——不仅因为形状简单灵活易于适应区域形状，还因为采用**面积坐标**后单元形函数的建立简便且标准化。

**面积坐标（Area Coordinate）** $L_1, L_2, L_3$：
$$L_1 = \frac{\triangle QQ_2Q_3}{\triangle Q_1Q_2Q_3},\quad L_2 = \frac{\triangle Q_1QQ_3}{\triangle Q_1Q_2Q_3},\quad L_3 = \frac{\triangle Q_1Q_2Q}{\triangle Q_1Q_2Q_3}$$

其中 $\triangle Q_1Q_2Q_3$ 为三角形总面积。

**面积坐标的性质**：
1. 三个顶点的面积为坐标：$(1,0,0), (0,1,0), (0,0,1)$；形心 $(\frac13,\frac13,\frac13)$
2. 三边方程：$L_1=0$（$Q_2Q_3$ 边），$L_2=0$（$Q_1Q_3$ 边），$L_3=0$（$Q_1Q_2$ 边）
3. 平行于对边的直线上任一点的 $L_i$ 值相同
4. 任一 $x,y$ 的 $k$ 次多项式可唯一转化为 $L_1,L_2,L_3$ 的 $k$ 次齐次多项式

**面积坐标与直角坐标的互化**：
$$\begin{pmatrix}x\\y\\1\end{pmatrix} = \begin{pmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
1 & 1 & 1
\end{pmatrix}\begin{pmatrix}L_1\\L_2\\L_3\end{pmatrix}$$

**导数的链式变换**：
$$\begin{cases}
\frac{\partial}{\partial x} = \frac{1}{2\Delta_e}(a_1\frac{\partial}{\partial L_1} + a_2\frac{\partial}{\partial L_2} + a_3\frac{\partial}{\partial L_3}) \\
\frac{\partial}{\partial y} = \frac{1}{2\Delta_e}(b_1\frac{\partial}{\partial L_1} + b_2\frac{\partial}{\partial L_2} + b_3\frac{\partial}{\partial L_3})
\end{cases}$$

**面积坐标下的积分公式**（极为有用）：
$$\iint_{\Omega_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}dxdy = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!}\,2\Delta_e$$

#### 高阶 Lagrange 插值

以 6 节点二次三角形单元为例，各节点的面积坐标：
- 节点1(1,0,0)、节点2(0,1,0)、节点3(0,0,1)
- 节点4($\frac12,\frac12,0$)、节点5($0,\frac12,\frac12$)、节点6($\frac12,0,\frac12$)

**划线法**：建立形函数时，对节点 $i$，找出除 $i$ 外所有经过其他节点的直线方程，将左侧表达式相乘再归一化。

例如节点 1：经过节点 4 和 6 的直线为 $L_1=0$，经过节点 2,5,3 的直线为 $L_1-\frac12=0$：
$$N_1 = \frac{L_1 - 1/2}{1-1/2}\cdot\frac{L_1}{1} = (2L_1-1)L_1$$

节点 4：经过节点 1,6 的直线为 $L_2=0$，经过节点 2,5 的直线为 $L_3=0$：
$$N_4 = \frac{L_2}{1/2}\cdot\frac{L_3}{1/2} = 4L_1L_2$$

### 2.2 Rectangular element（矩形单元）

矩形单元对边界形状的适应性不如三角形单元，但精度较高，与其他单元配合时可以发挥优势。

采用局部坐标 $(\xi,\eta)\in[-1,1]$。标准矩形为 $D: \{-1\leq\xi\leq1, -1\leq\eta\leq1\}$。

#### Lagrange 单元

双线性插值：
$$N_i = \frac14(1+\xi_i\xi)(1+\eta_i\eta),\quad i=1,2,3,4$$

这是两个一维线性插值基函数的乘积。

#### Serendipity 单元

为减少单元内部的节点而提出。以 8 节点二次 Serendipity 四边形为例——4 个角点 + 4 个边中点，无内部节点。角节点形函数：
$$N_1 = \frac14(1+\xi)(1+\eta)(\xi+\eta-1)$$

---

## III. Three-dimensional situation（三维情况）

### 3.1 Tetrahedron element（四面体单元）

采用**体积坐标** $L_1,L_2,L_3,L_4$，积分公式：
$$\iiint_{V_e} L_1^{\alpha_1}L_2^{\alpha_2}L_3^{\alpha_3}L_4^{\alpha_4}dV = \frac{\alpha_1!\,\alpha_2!\,\alpha_3!\,\alpha_4!}{(\sum\alpha_i+3)!}\,6V_e$$

---

## IV. Isoparametric element and numerical integration（等参元与数值积分）

### 4.1 任意四边形单元

结合三角形和矩形单元的优点——内部精度高，边界逼近好。

通过等参变换将标准矩形映射到任意四边形。坐标变换公式：
$$x = \sum_{i=1}^4 x_i N_i(\xi,\eta),\quad y = \sum_{i=1}^4 y_i N_i(\xi,\eta)$$

即坐标变换和插值函数**采用相同的形函数和相同的节点**——称为**等参元**。

### 4.2 等参变换与 Jacobian 矩阵

$$\mathbf{J} = \begin{pmatrix}
\partial x/\partial\xi & \partial y/\partial\xi \\
\partial x/\partial\eta & \partial y/\partial\eta
\end{pmatrix}$$

导数关系：
$$\begin{pmatrix}\partial N_i/\partial x \\ \partial N_i/\partial y\end{pmatrix} = \mathbf{J}^{-1}\begin{pmatrix}\partial N_i/\partial\xi \\ \partial N_i/\partial\eta\end{pmatrix}$$

**坐标变换一对一的充要条件**：$|\mathbf{J}| \neq 0$。

### 4.3 Gauss 数值积分

对等参元，单元刚度矩阵需要数值积分：
$$\iint_{\Omega_e} f(x,y)dxdy = \iint_{-1}^1 f(\xi,\eta)|\mathbf{J}|d\xi d\eta \approx \sum_{i=1}^n\sum_{j=1}^n w_i w_j f(\xi_i,\eta_j)$$

一维 Gauss 积分点：
| $n$ | $\xi_i$ | $w_i$ |
|-----|---------|-------|
| 1 | 0 | 2 |
| 2 | $\pm1/\sqrt{3}$ | 1 |
| 3 | $\pm\sqrt{0.6}, 0$ | $5/9, 8/9$ |

---

## V. Conclusion（结论）

长度坐标、面积坐标和体积坐标均与单元形式无直接关系，统称为**自然坐标**。

采用多项式建立插值函数的原因：
1. 便于计算
2. 易于证明收敛性

实际计算中通常采用线性或二次单元以减少计算量。提高 FEM 精度的策略：
- **$h$ 方法**：不改变插值函数形式，逐步细化网格
- **$p$ 方法**：不改变网格划分，提高单元插值函数次数

---

## 检查你的理解

1. 形函数的重要性体现在哪些方面？单元构造前需要确定哪些因素？
2. 面积坐标与直角坐标的转换关系是什么？面积坐标下的积分公式是什么？
3. 等参元、超参元和次参元有何区别？
4. 在等参元中为什么要保证 $|\mathbf{J}| \neq 0$？
5. $h$ 方法和 $p$ 方法的区别是什么？

---

> **对应作业**：[HW3 Q4（Hermite 梁单元形函数）](../04-Homework-Solutions/2026w/HW3-Problem.md)
