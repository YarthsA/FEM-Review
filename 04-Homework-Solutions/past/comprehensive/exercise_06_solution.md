# 习题六 解答

---

## 1. Euler 积分公式推导

**题目**：利用 Euler 积分公式，推导

$$\iint_{OAB} \lambda_1^{\alpha_1} \lambda_2^{\alpha_2} \lambda_3^{\alpha_3} dxdy = \frac{\alpha_1!\alpha_2!\alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!} 2\Delta_e$$

**解答**：

**一维形式**：

$$\int_0^1 t^m (1-t)^n dt = \frac{m! \, n!}{(m+n+1)!}$$

### 二维面积坐标积分

**面积坐标定义**：

$$L_1 = \frac{\text{面积}(P, 2, 3)}{\Delta_e}, \quad L_2 = \frac{\text{面积}(P, 3, 1)}{\Delta_e}, \quad L_3 = \frac{\text{面积}(P, 1, 2)}{\Delta_e}$$

满足：$L_1 + L_2 + L_3 = 1$

**积分变换**：

$$\iint_{\Delta_e} f(L_1, L_2, L_3) \, dxdy = 2\Delta_e \int_0^1 \int_0^{1-L_1} f(L_1, L_2, 1-L_1-L_2) \, dL_2 \, dL_1$$

### 推导

设 $I = \iint_{\Delta_e} L_1^{\alpha_1} L_2^{\alpha_2} L_3^{\alpha_3} \, dxdy$

$$I = 2\Delta_e \int_0^1 \int_0^{1-L_1} L_1^{\alpha_1} L_2^{\alpha_2} (1-L_1-L_2)^{\alpha_3} \, dL_2 \, dL_1$$

对内层积分（$L_2$）：

令 $u = L_2/(1-L_1)$，则 $L_2 = u(1-L_1)$，$dL_2 = (1-L_1)du$

$$\int_0^{1-L_1} L_2^{\alpha_2} (1-L_1-L_2)^{\alpha_3} dL_2 = (1-L_1)^{\alpha_2+\alpha_3+1} \int_0^1 u^{\alpha_2}(1-u)^{\alpha_3} du$$

$$= (1-L_1)^{\alpha_2+\alpha_3+1} \cdot \frac{\alpha_2! \alpha_3!}{(\alpha_2+\alpha_3+1)!}$$

代入外层积分：

$$I = 2\Delta_e \cdot \frac{\alpha_2! \alpha_3!}{(\alpha_2+\alpha_3+1)!} \int_0^1 L_1^{\alpha_1} (1-L_1)^{\alpha_2+\alpha_3+1} dL_1$$

$$= 2\Delta_e \cdot \frac{\alpha_2! \alpha_3!}{(\alpha_2+\alpha_3+1)!} \cdot \frac{\alpha_1! (\alpha_2+\alpha_3+1)!}{(\alpha_1+\alpha_2+\alpha_3+2)!}$$

$$= 2\Delta_e \cdot \frac{\alpha_1! \alpha_2! \alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!}$$

**结论**：

$$\boxed{\iint_{OAB} \lambda_1^{\alpha_1} \lambda_2^{\alpha_2} \lambda_3^{\alpha_3} \, dxdy = \frac{\alpha_1! \alpha_2! \alpha_3!}{(\alpha_1+\alpha_2+\alpha_3+2)!} 2\Delta_e}$$

$\blacksquare$

---

## 2. 划线法构造二阶三角形单元插值函数

**题目**：利用划线法构造二阶三角形单元的插值函数。

**解答**：

**基本思想**：形函数 $N_i$ 在节点 $i$ 处为 1，在其他节点处为 0。用划线法构造：

1. 过除节点 $i$ 外的其他节点作曲面（直线），使每条曲面方程为 0
2. 将这些曲面方程相乘 → 乘积在所有其他节点处为 0，仅在节点 $i$ 处非零
3. 归一化使得 $N_i$ 在节点 $i$ 处为 1

### 二阶三角形单元

6 个节点：3 个顶点（1, 2, 3）+ 3 个边中点（4, 5, 6）

**顶点 $Q_1$ 的形函数 $N_1$** 需在 $Q_2, Q_3, Q_4, Q_5, Q_6$ 处为零。作两条直线：

- 直线 $Q_2Q_6$（即边 $Q_2Q_3$ 上的 $L_1=0$ 线）：方程 $L_1 = 0$
- 直线 $Q_3Q_5$（过 $Q_3$ 和 $Q_5$，即 $L_1=1/2$ 的等值线）：方程 $L_1 - \frac12 = 0$

$$N_1 = c \cdot L_1 \cdot \left(L_1 - \frac12\right)$$

归一化：在 $Q_1$ 处 $L_1=1$，要求 $N_1=1$：
$$c \cdot 1 \cdot \left(1 - \frac12\right) = 1 \quad\Rightarrow\quad c = 2$$
$$\boxed{N_1 = 2L_1\left(L_1 - \frac12\right) = L_1(2L_1 - 1)}$$

验证：
- $Q_1$($L_1=1$)：$N_1 = 1\cdot(2-1)=1$ ✓
- $Q_2$($L_1=0$)：$N_1=0$ ✓
- $Q_4$($L_1=1/2$)：$N_1=\frac12\cdot(1-1)=0$ ✓

**类似地**（顶点形函数通式 $N_i = L_i(2L_i-1)$）：

$$N_2 = L_2(2L_2 - 1)$$
$$N_3 = L_3(2L_3 - 1)$$

**边中点形函数**：以 $Q_4$（边 $Q_1Q_2$ 中点）为例，需在 $Q_1, Q_2, Q_3, Q_5, Q_6$ 处为零。作两条直线：

- 直线 $Q_2Q_3$（即 $L_1=0$）：方程 $L_1 = 0$
- 直线 $Q_1Q_3$（即 $L_2=0$）：方程 $L_2 = 0$

$$N_4 = c \cdot L_1 \cdot L_2$$

归一化：在 $Q_4$（边中点）处 $L_1=L_2=1/2$，要求 $N_4=1$：
$$c \cdot \frac12 \cdot \frac12 = 1 \quad\Rightarrow\quad c = 4$$
$$\boxed{N_4 = 4L_1L_2}$$

验证：
- $Q_1$($L_1=1,L_2=0$)：$N_4=0$ ✓
- $Q_4$($L_1=1/2,L_2=1/2$)：$N_4 = 4\cdot\frac12\cdot\frac12=1$ ✓

同理（边中点形函数通式 $N_k = 4L_iL_j$，其中 $i,j$ 为该边两端顶点编号）：

$$N_5 = 4L_2L_3$$
$$N_6 = 4L_3L_1$$

---

## 3. 二次四边形单元的形函数偏导数

**题目**：图所示为二次四边形单元，试计算 $\frac{\partial N_1}{\partial x}$、$\frac{\partial N_2}{\partial y}$ 在自然坐标为 $(1/2, 1/2)$ 的点 Q 的值。节点坐标：1(40,50), 2(5,40), 3(10,10), 4(30,20)。

![四边形单元](exercise_06/images/f5b2e7f7723a3386a61ce674aa57ab43106ebc9379532dac0077fcbd5ae6150d.jpg)

**解答**：

**节点坐标**：

- 节点 1：$(x_1, y_1) = (40, 50)$　节点 2：$(x_2, y_2) = (5, 40)$
- 节点 3：$(x_3, y_3) = (10, 10)$　节点 4：$(x_4, y_4) = (30, 20)$

**形函数**：

$$N_i = \frac{1}{4}(1 + \xi\xi_i)(1 + \eta\eta_i)$$

节点自然坐标：$(\xi_1,\eta_1)=(1,1)$, $(\xi_2,\eta_2)=(-1,1)$, $(\xi_3,\eta_3)=(-1,-1)$, $(\xi_4,\eta_4)=(1,-1)$

**Step 1：形函数在 $Q(\xi,\eta)=(1/2,1/2)$ 处的偏导数**

$$\frac{\partial N_i}{\partial \xi} = \frac{\xi_i}{4}(1 + \eta\eta_i),\quad \frac{\partial N_i}{\partial \eta} = \frac{\eta_i}{4}(1 + \xi\xi_i)$$

| $i$ | $\xi_i$ | $\eta_i$ | $\partial N_i/\partial\xi$ | $\partial N_i/\partial\eta$ |
|-----|---------|---------|---------------------------|---------------------------|
| 1 | 1 | 1 | $\frac14(1+\frac12)=\frac38$ | $\frac14(1+\frac12)=\frac38$ |
| 2 | -1 | 1 | $-\frac14(1+\frac12)=-\frac38$ | $\frac14(1-\frac12)=\frac18$ |
| 3 | -1 | -1 | $-\frac14(1-\frac12)=-\frac18$ | $-\frac14(1-\frac12)=-\frac18$ |
| 4 | 1 | -1 | $\frac14(1-\frac12)=\frac18$ | $-\frac14(1+\frac12)=-\frac38$ |

**Step 2：Jacobi 矩阵在 Q 点的值**

$$J = \begin{bmatrix} \frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\ \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta} \end{bmatrix}$$

$$\frac{\partial x}{\partial \xi} = \sum_{i=1}^4 \frac{\partial N_i}{\partial \xi} x_i = \frac38(40) + \left(-\frac38\right)(5) + \left(-\frac18\right)(10) + \frac18(30) = \frac{125}{8}$$

$$\frac{\partial x}{\partial \eta} = \sum_{i=1}^4 \frac{\partial N_i}{\partial \eta} x_i = \frac38(40) + \frac18(5) + \left(-\frac18\right)(10) + \left(-\frac38\right)(30) = \frac{25}{8}$$

$$\frac{\partial y}{\partial \xi} = \sum_{i=1}^4 \frac{\partial N_i}{\partial \xi} y_i = \frac38(50) + \left(-\frac38\right)(40) + \left(-\frac18\right)(10) + \frac18(20) = 5$$

$$\frac{\partial y}{\partial \eta} = \sum_{i=1}^4 \frac{\partial N_i}{\partial \eta} y_i = \frac38(50) + \frac18(40) + \left(-\frac18\right)(10) + \left(-\frac38\right)(20) = 15$$

$$J = \begin{bmatrix} \frac{125}{8} & 5 \\ \frac{25}{8} & 15 \end{bmatrix},\quad \det(J) = \frac{125}{8}\cdot15 - 5\cdot\frac{25}{8} = \frac{1750}{8} = \frac{875}{4}$$

**Step 3：求 $J^{-1}$**

$$J^{-1} = \frac{1}{\det(J)}\begin{bmatrix} 15 & -5 \\ -\frac{25}{8} & \frac{125}{8} \end{bmatrix} = \frac{8}{1750}\begin{bmatrix} 15 & -5 \\ -\frac{25}{8} & \frac{125}{8} \end{bmatrix} = \begin{bmatrix} \frac{12}{175} & -\frac{4}{175} \\ -\frac{1}{70} & \frac{1}{14} \end{bmatrix}$$

**Step 4：计算 $\partial N_1/\partial x$ 和 $\partial N_2/\partial y$（Q 点处的值）**

$$\begin{aligned}
\frac{\partial N_1}{\partial x} &= J^{-1}_{11}\frac{\partial N_1}{\partial \xi} + J^{-1}_{12}\frac{\partial N_1}{\partial \eta} \\
&= \frac{12}{175}\cdot\frac38 + \left(-\frac{4}{175}\right)\cdot\frac38 \\
&= \frac{36-12}{1400} = \frac{24}{1400} = \boxed{\frac{3}{175}}
\end{aligned}$$

$$\begin{aligned}
\frac{\partial N_2}{\partial y} &= J^{-1}_{21}\frac{\partial N_2}{\partial \xi} + J^{-1}_{22}\frac{\partial N_2}{\partial \eta} \\
&= \left(-\frac{1}{70}\right)\left(-\frac38\right) + \frac{1}{14}\cdot\frac18 \\
&= \frac{3}{560} + \frac{5}{560} = \frac{8}{560} = \boxed{\frac{1}{70}}
\end{aligned}$$

---

## 4. 正六面体的变形状态和线应变

**题目**：正六面体弹性体，位移分量 $U = a_1xyz$，$V = a_2xyz$，$W = a_3xyz$。变形前 E 点坐标 (1.5, 1.0, 2.0)，变形后移至 (1.053, 1.001, 1.997)。求 E 点的变形状态和沿 EA 方向的线应变。

![正六面体](exercise_06/images/0a05f71afa537b654be1dee1c4e76959cba0890b36edcc213ff14459eda82222.jpg)

**解答**：

由图可知：正六面体尺寸为 1cm × 1.5cm × 2cm（x × y × z），O 为原点，A 点坐标 (1,0,0)。

变形前 E 点坐标：$(x, y, z) = (1.5, 1.0, 2.0)$

变形后 E 点坐标：$(x', y', z') = (1.053, 1.001, 1.997)$

**变形关系**：

$$x' = x + U = x + a_1 xyz$$

$$y' = y + V = y + a_2 xyz$$

$$z' = z + W = z + a_3 xyz$$

代入 E 点坐标：

$$1.053 = 1.5 + a_1(1.5)(1.0)(2.0) = 1.5 + 3a_1$$

$$a_1 = \frac{1.053 - 1.5}{3} = \frac{-0.447}{3} = -0.149$$

$$1.001 = 1.0 + a_2(1.5)(1.0)(2.0) = 1.0 + 3a_2$$

$$a_2 = \frac{1.001 - 1.0}{3} = \frac{0.001}{3} = 0.000333$$

$$1.997 = 2.0 + a_3(1.5)(1.0)(2.0) = 2.0 + 3a_3$$

$$a_3 = \frac{1.997 - 2.0}{3} = \frac{-0.003}{3} = -0.001$$

### 应变分量

**小变形假设下的应变**（$\varepsilon_{ij}$ 为张量应变，工程剪应变 $\gamma_{ij}=2\varepsilon_{ij}$）：

$$\varepsilon_{xx} = \frac{\partial U}{\partial x} = a_1 yz$$

$$\varepsilon_{yy} = \frac{\partial V}{\partial y} = a_2 xz$$

$$\varepsilon_{zz} = \frac{\partial W}{\partial z} = a_3 xy$$

$$\varepsilon_{xy} = \frac{1}{2}\left(\frac{\partial U}{\partial y} + \frac{\partial V}{\partial x}\right) = \frac{1}{2}(a_1 xz + a_2 yz)$$

$$\varepsilon_{yz} = \frac{1}{2}\left(\frac{\partial V}{\partial z} + \frac{\partial W}{\partial y}\right) = \frac{1}{2}(a_2 xy + a_3 xz)$$

$$\varepsilon_{zx} = \frac{1}{2}\left(\frac{\partial W}{\partial x} + \frac{\partial U}{\partial z}\right) = \frac{1}{2}(a_3 yz + a_1 xy)$$

### E 点的应变

代入 $(x, y, z) = (1.5, 1.0, 2.0)$：

$$\varepsilon_{xx} = -0.149 \times 1.0 \times 2.0 = -0.298$$

$$\varepsilon_{yy} = 0.000333 \times 1.5 \times 2.0 = 0.001$$

$$\varepsilon_{zz} = -0.001 \times 1.5 \times 1.0 = -0.0015$$

### EA 方向的线应变

**EA 方向**：从 E 到 A 的方向向量。

假设 A 点坐标为 $(1, 0, 0)$（根据图示），则：

方向向量 $\mathbf{n} = \frac{A - E}{|A - E|}$

$$\mathbf{n} = \frac{(1-1.5, 0-1.0, 0-2.0)}{|(1-1.5, 0-1.0, 0-2.0)|} = \frac{(-0.5, -1.0, -2.0)}{\sqrt{0.25 + 1 + 4}} = \frac{(-0.5, -1.0, -2.0)}{\sqrt{5.25}}$$

**线应变**：

$$\varepsilon_{EA} = \mathbf{n}^T [\varepsilon] \mathbf{n}$$

展开计算（$n_i n_j \varepsilon_{ij}$ 求和）：

$$\begin{aligned}
\varepsilon_{EA} &= n_1^2\varepsilon_{xx} + n_2^2\varepsilon_{yy} + n_3^2\varepsilon_{zz} + 2n_1n_2\varepsilon_{xy} + 2n_2n_3\varepsilon_{yz} + 2n_3n_1\varepsilon_{zx}
\end{aligned}$$

代入 $n_1=-0.5/\sqrt{5.25}$, $n_2=-1/\sqrt{5.25}$, $n_3=-2/\sqrt{5.25}$ 及各应变分量：

$$n_1^2=\frac{0.25}{5.25},\ n_2^2=\frac{1}{5.25},\ n_3^2=\frac{4}{5.25},\ n_1n_2=\frac{0.5}{5.25},\ n_2n_3=\frac{2}{5.25},\ n_3n_1=\frac{1}{5.25}$$

代入 E 点应变值 $\varepsilon_{xx}=-0.298$, $\varepsilon_{yy}=0.001$, $\varepsilon_{zz}=-0.0015$, $\varepsilon_{xy}=$（需补充 $\varepsilon_{yz},\varepsilon_{zx}$ 数值）：

补充计算各剪应变在 E 点的值：
$$\varepsilon_{xy} = \frac12(-0.149\times1.5\times2.0 + 0.000333\times1.0\times2.0) = \frac12(-0.447+0.000666) = -0.223$$
$$\varepsilon_{yz} = \frac12(0.000333\times1.5\times1.0 + (-0.001)\times1.5\times2.0) = \frac12(0.0005-0.003) = -0.00125$$
$$\varepsilon_{zx} = \frac12((-0.001)\times1.0\times2.0 + (-0.149)\times1.5\times1.0) = \frac12(-0.002-0.2235) = -0.113$$

$$\begin{aligned}
\varepsilon_{EA} &= \frac{0.25}{5.25}(-0.298) + \frac{1}{5.25}(0.001) + \frac{4}{5.25}(-0.0015) \\
&\quad + 2\cdot\frac{0.5}{5.25}(-0.223) + 2\cdot\frac{2}{5.25}(-0.00125) + 2\cdot\frac{1}{5.25}(-0.113) \\[4pt]
&= \frac{1}{5.25}\big[-0.0745 + 0.001 - 0.006 - 0.223 - 0.005 - 0.226\big] \\[4pt]
&= \frac{-0.5335}{5.25} = \boxed{-0.102}
\end{aligned}$$

---

## 5. 二杆结构的有限元分析

**题目**：两根杆组成的结构，$E_1 = E_2 = 2\times10^6$ kg/cm²，$A_1 = 2A_2 = 2$ cm²。试完成：(1) 各单元刚度矩阵；(2) 总刚度矩阵；(3) 节点 2 位移；(4) 各单元应力；(5) 支反力。

![二杆结构](exercise_06/images/cd76c0075cbba5b53ba20bdab9f713fc53773a1f9f07bb9bd8a757560022f5fa.jpg)

**解答**：

由图可知：
- 节点 1：左侧固定端，节点 2：中间自由端，节点 3：上方固定端
- 单元 1（水平）：连接节点 1→2，长度 $L_1=10$cm，$E_1A_1=2\times10^6\times2=4\times10^6$
- 单元 2（垂直）：连接节点 2→3，长度 $L_2=10$cm，$E_2A_2=2\times10^6\times1=2\times10^6$
- 力 $P=\sqrt{2}$ kg 于节点 2，方向 45° → $P_x = P_y = 1$ kg


### (1) 单元刚度矩阵（从形函数推导）

**A. 局部坐标系下的 1D 杆单元**

先将杆放在沿杆轴的局部坐标 $\bar{x}$（$0\le\bar{x}\le L$）下推导。

**Step A1 — 形函数**（$\xi = \bar{x}/L \in [0,1]$）：

$$N_1(\xi) = 1 - \xi,\quad N_2(\xi) = \xi$$

**Step A2 — 位移场**：

$$u'(\xi) = \begin{bmatrix} N_1 & N_2 \end{bmatrix} \begin{Bmatrix} u'_1 \\ u'_2 \end{Bmatrix}$$

**Step A3 — 形函数对 $\bar{x}$ 求偏导**：

$$\frac{dN_1}{d\bar{x}} = -\frac{1}{L},\quad \frac{dN_2}{d\bar{x}} = \frac{1}{L}$$

**Step A4 — 几何方程（轴向应变）**：

$$\varepsilon = \frac{du'}{d\bar{x}} = \underbrace{\begin{bmatrix} -\dfrac{1}{L} & \dfrac{1}{L} \end{bmatrix}}_{[B']_{1\times2}} \begin{Bmatrix} u'_1 \\ u'_2 \end{Bmatrix}$$

**Step A5 — 局部单元刚度矩阵**（$[D]=E$）：

$$[k'] = \int_0^L [B']^T E A [B']\,d\bar{x} = \frac{EA}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

**B. 坐标变换到全局**

**Step B1 — 变换矩阵 $[T]$**：方向余弦 $l=\cos\theta,\ m=\sin\theta$，将局部轴向位移投影到全局 DOF：

$$\begin{Bmatrix} u'_1 \\ u'_2 \end{Bmatrix}
= \underbrace{\begin{bmatrix} l & m & 0 & 0 \\ 0 & 0 & l & m \end{bmatrix}}_{[T]_{2\times4}}
\begin{Bmatrix} u_{1x} \\ u_{1y} \\ u_{2x} \\ u_{2y} \end{Bmatrix}$$

**Step B2 — 全局刚度矩阵**：$[k]^e = [T]^T [k'] [T]$

$$\boxed{[k]^e = \frac{EA}{L}
\begin{bmatrix}
l^2 & lm & -l^2 & -lm\\
lm & m^2 & -lm & -m^2\\
-l^2 & -lm & l^2 & lm\\
-lm & -m^2 & lm & m^2
\end{bmatrix}}$$

**C. 代入本题数据**

**单元 1**（水平，1→2）：$L_1=10$ cm，$E_1A_1=2\times10^6\times2=4\times10^6$，$(l,m) = (1,0)$

$$k^{(1)} = 4\times10^5\begin{bmatrix}
1 & 0 & -1 & 0\\
0 & 0 & 0 & 0\\
-1 & 0 & 1 & 0\\
0 & 0 & 0 & 0
\end{bmatrix}
\begin{array}{l}
\text{DOF 1}(u_{1x})\\
\text{DOF 2}(u_{1y})\\
\text{DOF 3}(u_{2x})\\
\text{DOF 4}(u_{2y})
\end{array}$$

**单元 2**（垂直，2→3）：$L_2=10$ cm，$E_2A_2=2\times10^6\times1=2\times10^6$，$(l,m) = (0,1)$

$$k^{(2)} = 2\times10^5\begin{bmatrix}
0 & 0 & 0 & 0\\
0 & 1 & 0 & -1\\
0 & 0 & 0 & 0\\
0 & -1 & 0 & 1
\end{bmatrix}
\begin{array}{l}
\text{DOF 3}(u_{2x})\\
\text{DOF 4}(u_{2y})\\
\text{DOF 5}(u_{3x})\\
\text{DOF 6}(u_{3y})
\end{array}$$

### (2) 总体刚度矩阵（$6\times6$）

组装得：

$$K = \begin{bmatrix}
4e5 & 0   & -4e5 & 0   & 0   & 0\\
0   & 0   & 0   & 0   & 0   & 0\\
-4e5& 0   & 4e5 & 0   & 0   & 0\\
0   & 0   & 0   & 2e5 & 0   & -2e5\\
0   & 0   & 0   & 0   & 0   & 0\\
0   & 0   & 0   & -2e5& 0   & 2e5
\end{bmatrix}
\begin{array}{l}
u_{1x}\\
u_{1y}\\
u_{2x}\\
u_{2y}\\
u_{3x}\\
u_{3y}
\end{array}$$

引入边界条件 $u_{1x}=u_{1y}=u_{3x}=u_{3y}=0$，缩减为 $2\times2$：

$$K_{\text{red}} = \begin{bmatrix}
4\times10^5 & 0\\
0 & 2\times10^5
\end{bmatrix}$$

### (3) 节点 2 位移

$$\begin{bmatrix}
4\times10^5 & 0\\
0 & 2\times10^5
\end{bmatrix}
\begin{Bmatrix}u_{2x}\\u_{2y}\end{Bmatrix}
= \begin{Bmatrix}1\\1\end{Bmatrix}$$

$$\boxed{u_{2x} = \frac{1}{4\times10^5} = 2.5\times10^{-6}\ \text{cm}},\quad
\boxed{u_{2y} = \frac{1}{2\times10^5} = 5\times10^{-6}\ \text{cm}}$$

### (4) 各单元应力

$$\sigma_1 = E_1\frac{u_{2x}-u_{1x}}{L_1} = 2\times10^6\times\frac{2.5\times10^{-6}}{10} = \boxed{0.5\ \text{kg/cm}^2}$$

$$\sigma_2 = E_2\frac{u_{2y}-u_{3y}}{L_2} = 2\times10^6\times\frac{5\times10^{-6}}{10} = \boxed{1.0\ \text{kg/cm}^2}$$

### (5) 支反力

$$R = K\{u\} - \{F\}$$

$$R_{1x} = -4\times10^5\cdot u_{2x} = -4\times10^5\cdot2.5\times10^{-6} = \boxed{-1\ \text{kg}}$$
$$R_{3y} = -2\times10^5\cdot u_{2y} = -2\times10^5\cdot5\times10^{-6} = \boxed{-1\ \text{kg}}$$

其余支反力分量为零。水平方向由节点 1 提供支反力，垂直方向由节点 3 提供支反力，与 $P_x=P_y=1$ 平衡。$\blacksquare$

---

## 6. 平面桁架的节点位移和单元内力

**题目**：求平面桁架的节点位移和单元内力。$E = 2\times10^6$ MPa，$A = 1$ cm²。

![平面桁架](exercise_06/images/4624a345311b2e9a0a85262b494caa250adf94c850d1c6c28cb2bbc20e2110d5.jpg)

**解答**：

由图及坐标标注可知：
- 节点 1 $(40,0)$：右侧自由端，受 10 N 竖直向下荷载
- 节点 2 $(0,0)$：左下角固定端（$u_2=v_2=0$）
- 节点 3 $(0,30)$：左上角固定端（$u_3=v_3=0$）

| 单元 | 连接 | 方向 | $L$ (cm) | $\frac{EA}{L}$ (N/cm) |
|------|------|------|----------|----------------------|
| 1 | 1→2 | 水平 | 40 | $\frac{2\times10^8}{40}=5\times10^6$ |
| 2 | 2→3 | 垂直 | 30 | $\frac{2\times10^8}{30}=6.667\times10^6$ |
| 3 | 1→3 | 斜杆 | 50 | $\frac{2\times10^8}{50}=4\times10^6$ |

**单位说明**：$E=2\times10^6$ MPa $=2\times10^{12}$ Pa $=2\times10^{12}$ N/m²，$A=1$ cm² $=10^{-4}$ m²，故 $EA=2\times10^8$ N。

DOF 编号：节点 1→(1,2)、节点 2→(3,4)、节点 3→(5,6)。

### (1) 单元刚度矩阵（全局坐标）

平面杆单元在全局坐标下的刚度矩阵：

$$k^e = \frac{EA}{L}\begin{bmatrix}
l^2 & lm & -l^2 & -lm\\
lm & m^2 & -lm & -m^2\\
-l^2 & -lm & l^2 & lm\\
-lm & -m^2 & lm & m^2
\end{bmatrix}$$

其中 $(l,m)$ 为从单元起点到终点的方向余弦。

**单元 1（1→2，水平）**：$l=(0-40)/40=-1,\ m=0$，DOF {1,2,3,4}

$$k^{(1)} = 5\times10^6\begin{bmatrix}
1 & 0 & -1 & 0\\
0 & 0 & 0 & 0\\
-1 & 0 & 1 & 0\\
0 & 0 & 0 & 0
\end{bmatrix}$$

**单元 2（2→3，垂直）**：$l=0,\ m=(30-0)/30=1$，DOF {3,4,5,6}

$$k^{(2)} = 6.667\times10^6\begin{bmatrix}
0 & 0 & 0 & 0\\
0 & 1 & 0 & -1\\
0 & 0 & 0 & 0\\
0 & -1 & 0 & 1
\end{bmatrix}$$

**单元 3（1→3，斜杆）**：$l=(0-40)/50=-\frac45,\ m=(30-0)/50=\frac35$，DOF {1,2,5,6}

$$l^2=\frac{16}{25},\ m^2=\frac{9}{25},\ lm=-\frac{12}{25}$$

$$k^{(3)} = 4\times10^4\begin{bmatrix}
\frac{16}{25} & -\frac{12}{25} & -\frac{16}{25} & \frac{12}{25}\\[4pt]
-\frac{12}{25} & \frac{9}{25} & \frac{12}{25} & -\frac{9}{25}\\[4pt]
-\frac{16}{25} & \frac{12}{25} & \frac{16}{25} & -\frac{12}{25}\\[4pt]
\frac{12}{25} & -\frac{9}{25} & -\frac{12}{25} & \frac{9}{25}
\end{bmatrix}$$

### (2) 总体刚度矩阵（$6\times6$）

组装得：

$$
K = \begin{bmatrix}
7.56 & -1.92 & -5 & 0 & -2.56 & 1.92\\
-1.92 & 1.44 & 0 & 0 & 1.92 & -1.44\\
-5 & 0 & 5 & 0 & 0 & 0\\
0 & 0 & 0 & 6.667 & 0 & -6.667\\
-2.56 & 1.92 & 0 & 0 & 2.56 & -1.92\\
1.92 & -1.44 & 0 & -6.667 & -1.92 & 8.107
\end{bmatrix}\times10^6
\begin{array}{l}
u_{1x}\\
u_{1y}\\
u_{2x}\\
u_{2y}\\
u_{3x}\\
u_{3y}
\end{array}
$$

### (3) 边界条件与缩减

$u_{2x}=u_{2y}=u_{3x}=u_{3y}=0$（节点 2、3 固定），缩减为节点 1 的 $2\times2$ 方程：

$$\begin{bmatrix}
7.56\times10^6 & -1.92\times10^6\\
-1.92\times10^6 & 1.44\times10^6
\end{bmatrix}
\begin{Bmatrix}u_{1x}\\u_{1y}\end{Bmatrix}
= \begin{Bmatrix}0\\-10\end{Bmatrix}$$

### (4) 节点位移

$$\det = (7.56\times10^6)(1.44\times10^6) - (-1.92\times10^6)^2 = 7.2\times10^{12}$$

$$u_{1x} = \frac{0\times1.44\times10^6 - (-1.92\times10^6)\times(-10)}{7.2\times10^{12}}
= \boxed{-2.67\times10^{-6}\ \text{cm}}$$

$$u_{1y} = \frac{7.56\times10^6\times(-10) - 0}{7.2\times10^{12}}
= \boxed{-1.05\times10^{-5}\ \text{cm}}$$

### (5) 各单元内力

杆单元轴力：$F = \frac{EA}{L}\big[(u_j-u_i)l + (v_j-v_i)m\big]$

**单元 1**（1→2, $l=-1,m=0$）：

$$F_1 = 5\times10^6\big[(0-(-2.67\times10^{-6}))(-1) + 0\big] = 5\times10^6\times(-2.67\times10^{-6}) = \boxed{-13.3\ \text{N}}$$

（负号表示受压）

**单元 2**（2→3, $l=0,m=1$）：

$$F_2 = 6.667\times10^6\big[0 + 0\big] = \boxed{0}$$

（两端固定无相对位移）

**单元 3**（1→3, $l=-\frac45,m=\frac35$）：

$$\begin{aligned}
F_3 &= 4\times10^6\left[(0-(-2.67\times10^{-6}))\left(-\frac45\right) + (0-(-1.05\times10^{-5}))\left(\frac35\right)\right] \\
&= 4\times10^6\big[(-2.136\times10^{-6}) + (6.30\times10^{-6})\big] \\
&= 4\times10^6\times4.164\times10^{-6} = \boxed{16.7\ \text{N}}
\end{aligned}$$

（正号表示受拉）

### (6) 支反力

$$R_{2x} = -5\times10^6\times u_{1x} = -5\times10^6\times(-2.67\times10^{-6}) = \boxed{13.3\ \text{N}}$$
$$R_{3x} = -2.56\times10^6\times u_{1x} + 1.92\times10^6\times u_{1y} = 6.84 - 20.16 = \boxed{-13.3\ \text{N}}$$
$$R_{3y} = 1.92\times10^6\times u_{1x} - 1.44\times10^6\times u_{1y} - 6.667\times10^6\times0 = -5.13 + 15.12 = \boxed{10\ \text{N}}$$

验算：$R_{2x}+R_{3x}=0$（水平平衡 ✓），$R_{3y}=10$ N（竖直平衡 ✓）。$\blacksquare$
