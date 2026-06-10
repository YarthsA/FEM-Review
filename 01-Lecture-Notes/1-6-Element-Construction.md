# 第6章：单元构造与形函数

> **对应 PDF**：[`6 FEM_Element construction.pdf`](../06-References/pdfs-originals/6%20FEM_Element%20construction.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §5
> **相关作业**：[HW3 Q4](../04-Homework-Solutions/2026w/HW3-Problem.md)（Hermite 梁单元）
> **前置知识**：第 5 章（FEM 公式）、线性代数（多项式插值）

形函数的构造是 FEM 中最重要也最技巧性的环节，因为：
1. 单元刚度矩阵计算、总体集成和求解可标准化、自动化
2. 基函数由形函数生成，直接关系到收敛性和精度

### 单元构造前需确定的因素

**① 单元几何形状**
- 1D：直线段 / 曲线
- 2D：三角形 / 矩形 / 任意四边形
- 3D：四面体 / 五面体 / 三棱柱 / 六面体

**② 节点数量和分布**
- 边界节点 / 内部节点
- 通常布置在特殊位置（端点、顶点、边中点、形心）

**③ 节点自由度 (DOF)**
- **Lagrange 型**：仅位移，$n$ 维问题每节点 $n$ 个 DOF
- **Hermite 型**：位移+导数（转角），$n$ 维问题每节点 $2n$ 个 DOF

---

## 2. 插值基本理论

### 2.1 Lagrange 插值
仅要求插值多项式在插值点上的函数值已知：
$$
P_n(x) = \sum_{i=0}^n f(x_i)L_i(x)
$$

其中 Lagrange 基函数：
$$
L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}
$$

### 2.2 Hermite 插值
除函数值外，还要求导数已知（包括一阶、二阶等）。

---

## 3. 一维单元

### 3.1 线性 Lagrange 单元（2 节点）

形函数（自然坐标 $\xi \in [-1, 1]$）：
$$
N_1 = \frac{1-\xi}{2}, \quad N_2 = \frac{1+\xi}{2}
$$

### 3.2 二次 Lagrange 单元（3 节点）

形函数：
$$
N_1 = \frac{(\xi-1)\xi}{2}, \quad N_2 = \frac{(\xi+1)\xi}{2}, \quad N_3 = 1-\xi^2
$$

### 3.3 Hermite 单元（梁单元）

每节点 2 个 DOF：$u_i$ 和 $\theta_i = \frac{du}{dx}\big|_i$

形函数（Hermite 三次插值）：
$$
H_1^0 = \frac{1}{4}(\xi-1)^2(\xi+2), \quad H_2^0 = \frac{1}{4}(\xi+1)^2(2-\xi)
$$
$$
H_1^1 = \frac{L}{8}(\xi-1)^2(\xi+1), \quad H_2^1 = \frac{L}{8}(\xi+1)^2(\xi-1)
$$

---

## 4. 二维单元

### 4.1 3节点三角形（CST / 线性三角形）

形函数（面积坐标 $L_1, L_2, L_3$）：
$$
L_1 = \frac{A_1}{A}, \quad L_2 = \frac{A_2}{A}, \quad L_3 = \frac{A_3}{A}
$$

其中 $A$ 为三角形总面积，$A_i$ 为对边子三角形面积。

性质：$L_1 + L_2 + L_3 = 1$，$L_i$ 只有 2 个独立变量。
线性插值：$u = L_1 u_1 + L_2 u_2 + L_3 u_3$

**面积坐标与直角坐标的互化**：
$$
\begin{pmatrix} x \\ y \\ 1 \end{pmatrix}
= \begin{pmatrix} x_1 & x_2 & x_3 \\ y_1 & y_2 & y_3 \\ 1 & 1 & 1 \end{pmatrix}
\begin{pmatrix} L_1 \\ L_2 \\ L_3 \end{pmatrix}
$$

**面积坐标的重要性质**：
1. 三顶点坐标：$(1,0,0), (0,1,0), (0,0,1)$；形心 $(\frac13,\frac13,\frac13)$
2. 三边方程：$L_1=0$（$Q_2Q_3$ 边），$L_2=0$（$Q_1Q_3$ 边），$L_3=0$（$Q_1Q_2$ 边）
3. 平行于对边的直线上的点，$L_i$ 值相同

**面积坐标下的积分公式（核心！）**：
$$
\iint_{\Delta_e} L_1^{\alpha_1} L_2^{\alpha_2} L_3^{\alpha_3}\,dxdy
= \frac{\alpha_1! \,\alpha_2! \,\alpha_3!}{(\alpha_1 + \alpha_2 + \alpha_3 + 2)!}\,2\Delta_e
$$

> 这个公式在计算单元刚度矩阵和载荷向量时非常有用——直接用 $L_i$ 的指数确定积分值，无需显式积分。

**导数变换**（面积坐标 → 直角坐标）：
$$
\begin{cases}
\frac{\partial}{\partial x} = \frac{1}{2\Delta_e}\left(b_1\frac{\partial}{\partial L_1} + b_2\frac{\partial}{\partial L_2} + b_3\frac{\partial}{\partial L_3}\right) \\[4pt]
\frac{\partial}{\partial y} = \frac{1}{2\Delta_e}\left(c_1\frac{\partial}{\partial L_1} + c_2\frac{\partial}{\partial L_2} + c_3\frac{\partial}{\partial L_3}\right)
\end{cases}
$$

### 4.2 广义 Lagrange 插值公式（节点形函数通用构造法）

对任意 $n$ 节点单元，形函数的通用表达式：
$$
N_i = \prod_{j=1,\,j\neq i}^{n} \frac{f_j(\xi)}{f_j(\xi_i)}
$$

其中 $f_j(\xi) = \xi - \xi_j$ 表示点到节点 $j$ 的距离。
这个公式保证了 $N_i(\xi_j) = \delta_{ij}$。

### 4.3 6节点三角形（LST / 二次三角形）

节点：3 个顶点 + 3 条边中点

形函数（面积坐标，用"划线法"构造）：
$$
\begin{aligned}
N_1 &= \frac{L_1 - 1/2}{1/2}\cdot\frac{L_1}{1} = L_1(2L_1-1) \\
N_4 &= 4L_1L_2
\end{aligned}
$$

完整组：
$$
\begin{aligned}
N_1 &= L_1(2L_1-1),\quad N_2 = L_2(2L_2-1),\quad N_3 = L_3(2L_3-1) \\
N_4 &= 4L_1L_2,\quad N_5 = 4L_2L_3,\quad N_6 = 4L_3L_1
\end{aligned}
$$

**划线法原理**：对节点 $i$，找出所有除 $i$ 外经过其他节点的直线，将每条直线的左端表达式相乘，除以在节点 $i$ 处的取值归一化。

### 4.4 4节点矩形单元（Q4 / 双线性）

节点：3 个顶点 + 3 条边中点

形函数（面积坐标）：
$$
\begin{aligned}
N_1 &= L_1(2L_1 - 1), \quad N_2 = L_2(2L_2 - 1), \quad N_3 = L_3(2L_3 - 1) \\
N_4 &= 4L_1L_2, \quad N_5 = 4L_2L_3, \quad N_6 = 4L_3L_1
\end{aligned}
$$

特点：二次位移场 → 应变线性变化，精度高于 CST

### 4.3 4节点矩形单元（Q4 / 双线性）

形函数（自然坐标 $\xi, \eta \in [-1,1]$）：
$$
N_1 = \frac{(1-\xi)(1-\eta)}{4}, \quad N_2 = \frac{(1+\xi)(1-\eta)}{4}
$$
$$
N_3 = \frac{(1+\xi)(1+\eta)}{4}, \quad N_4 = \frac{(1-\xi)(1+\eta)}{4}
$$

### 4.4 Pascal 三角形（多项式完备性选择）

```
1
x    y
x²  xy  y²
x³  x²y xy² y³
x⁴  x³y x²y² xy³ y⁴
```

- 三角形单元：选取三角形内的项（对称）
- 矩形单元：选取矩形边界上的项

---

## 5. 三维单元

### 5.1 4节点四面体（线性四面体）

每个节点 3 个 DOF $(u,v,w)$

形函数（体积坐标 $L_1, L_2, L_3, L_4$）：
$$
N_i = L_i = \frac{V_i}{V} \quad (i=1,2,3,4)
$$

其中 $V$ 为四面体体积，$V_i$ 为对面子四面体体积。

线性位移场：
$$
\mathbf{u} = \sum_{i=1}^4 N_i \mathbf{u}_i = [N]\{\delta\}_e
$$

体积计算公式：
$$
V_e = \frac{1}{6}\begin{vmatrix}
1 & x_i & y_i & z_i \\
1 & x_j & y_j & z_j \\
1 & x_m & y_m & z_m \\
1 & x_l & y_l & z_l
\end{vmatrix}
$$

### 5.2 8节点六面体单元

类似 2D Q4 扩展到 3D：
$$
N_i = \frac{1}{8}(1+\xi\xi_i)(1+\eta\eta_i)(1+\zeta\zeta_i)
$$

---

## 6. 等参元 (Isoparametric Element)

### 核心思想
**坐标和位移使用相同的形函数**：
$$
x = \sum N_i(\xi,\eta)x_i, \quad y = \sum N_i(\xi,\eta)y_i
$$
$$
u = \sum N_i(\xi,\eta)u_i, \quad v = \sum N_i(\xi,\eta)v_i
$$

### 优点
- 可处理曲线边界
- 同一程序可处理不同形状单元
- 精度高

### Jacobian 矩阵

坐标变换的导数关系：
$$
\begin{pmatrix} \frac{\partial N_i}{\partial x} \\ \frac{\partial N_i}{\partial y} \end{pmatrix}
= \mathbf{J}^{-1} \begin{pmatrix} \frac{\partial N_i}{\partial \xi} \\ \frac{\partial N_i}{\partial \eta} \end{pmatrix}
$$

其中 Jacobian 矩阵：
$$
\mathbf{J} = \begin{pmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\
\frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
\end{pmatrix}
$$

### 6.1 等参元刚度矩阵数值计算

等参元中，单元刚度矩阵必须在局部坐标系下计算：
$$
[k]_e = \iint_{\Omega_e} [B]^T[D][B]\,t\,dxdy
        = \int_{-1}^1\int_{-1}^1 [B(\xi,\eta)]^T[D][B(\xi,\eta)]\,t\,|\mathbf{J}|\,d\xi d\eta
$$

其中 $[B]$ 矩阵中的 $\partial N_i/\partial x$ 通过 Jacobian 逆变换得到。

### 6.2 数值积分（Gauss 积分）

$$
\int_{-1}^{1} \int_{-1}^{1} f(\xi,\eta)\,d\xi d\eta \approx \sum_{i=1}^{n_g} \sum_{j=1}^{n_g} w_i w_j f(\xi_i,\eta_j)
$$

一维 Gauss 点：
| $n_g$ | 积分点 $\xi_i$ | 权系数 $w_i$ |
|-------|---------------|-------------|
| 1 | 0 | 2 |
| 2 | $\pm 1/\sqrt{3}$ | 1 |
| 3 | $\pm\sqrt{0.6}, 0$ | $5/9, 8/9$ | 5 次 |

**二维 Gauss 积分**（张量积形式）：
$$
\int_{-1}^1\int_{-1}^1 f(\xi,\eta)\,d\xi d\eta \approx \sum_{i=1}^n\sum_{j=1}^n w_i w_j f(\xi_i,\eta_j)
$$

**积分阶数选择原则**：
- Q4 单元精确积分双线性项需要 $2\times2$ Gauss 点
- Q8 单元需要 $3\times3$ Gauss 点
- **减缩积分 (Reduced Integration)**：使用比精确积分所需更低的阶数 → 可避免剪切自锁，但可能引入零能模式（hourglass mode）

---

## 7. Serendipity 单元

**动机**：Lagrange 矩形单元（如 9 节点 Q9）有内部节点，包含大量高次项但对精度提升贡献有限。Serendipity 单元将节点集中在边界上。

**8 节点 Serendipity 四边形**（Q8）：4 个角点 + 4 个边中点，精度与 Q9 接近但少 1 个内部节点。

形函数构造示例（角节点 1 和边中点 5）：
$$
N_1 = \frac14(1+\xi_1\xi)(1+\eta_1\eta)(\xi_1\xi + \eta_1\eta - 1),\quad
N_5 = \frac12(1-\xi^2)(1+\eta_5\eta)
$$

修正方法：当添加边节点后，原角节点的形函数会受新节点影响，需要减去修正项。例如 Q8 中角节点 1 需要减去 $N_5/2$ 和 $N_8/2$。

---

## 8. 等参/超参/次参单元

| 类型 | 坐标变换节点数 $m$ vs 插值节点数 $n$ | 说明 |
|------|-------------------------------------|------|
| **等参元** (Isoparametric) | $m = n$，形函数相同 | 最常用 |
| **超参元** (Superparametric) | $m > n$ | 坐标描述比函数更精确 |
| **次参元** (Subparametric) | $m < n$ | 函数描述比坐标更精确 |

等参元的收敛性验证：
- **协调性**：相邻单元在共享边上节点相同、形函数相同 → 位移连续 ✅
- **完备性**：验证 $\sum N_i = 1$ → 若成立，则能表示刚体位移 $u = a + bx + cy + dz$ ✅

**Jacobian 矩阵非奇异的条件**：
- $|\mathbf{J}| \neq 0$ 在整个单元内
- 避免单元畸形（边退化、内角过大、内角变号）

---

## 9. 收敛准则

1. **完备性 (Completeness)**：形函数必须能表示**刚体位移**和**常应变**状态
2. **协调性 (Compatibility)**：单元间边界上位移连续（$C^0$ 连续）
3. 满足完备性和协调性的单元称为**协调单元 (conforming element)** → 解单调收敛
4. 仅满足完备性的单元称为**非协调单元 (non-conforming element)** → 工程中有时使用

### 单元选择指南

| 问题类型 | 推荐单元 | 备注 |
|----------|----------|------|
| 杆/梁 | 1D Lagrange/Hermite | 梁建议用 Hermite |
| 2D 简单几何 | 3节点三角形 (CST) | 最简单，应力为常数 |
| 2D 复杂几何 | 6节点三角形 (LST) | 精度高 |
| 2D 规则区域 | 4节点矩形 (Q4) | 双线性 |
| 3D 任意几何 | 4节点四面体 | 网格划分容易 |
| 3D 规则区域 | 8节点六面体 | 精度更高 |
| 弯曲边界 | 等参元 | 精确描述边界 |

---

## 10. 形函数性质验证方法

考试中常要求验证形函数是否合法，标准检查项：

1. **Kronecker 性质**：$N_i(\text{节点}_j) = \delta_{ij}$
   - 在自身节点为 1，在其他节点为 0
2. **单位分解**：$\sum N_i = 1$
   - 能表示刚体平移（$u = \text{常数}$）
3. **线性完备性**：$\sum N_i x_i = x$，$\sum N_i y_i = y$
   - 能表示刚体转动和常应变状态

**例**：验证梁单元 Hermite 形函数 $N_1 = 1-3\xi^2+2\xi^3$
- $N_1(0) = 1$，$N_1(1) = 0$ ✅
- $N_1'(0) = 0$，$N_1'(1) = 0$ ✅（对转角无贡献）

---

> **对应作业**：[HW3 Q4: Hermite 梁单元形函数](../04-Homework-Solutions/2026w/HW3-Problem.md)
> **往年相关**：[Homework3 (past)](../04-Homework-Solutions/past/HW3/Homework3.md) · [答案(LIU Sai)](../04-Homework-Solutions/past/HW3/Ans%20to%20HM3_LIU%20Sai_handed%20in.md)
