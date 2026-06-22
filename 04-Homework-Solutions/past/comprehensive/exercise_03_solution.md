# 习题三 解答

---

## 1. 矩形截面杆的杆长变化

**题目**：考虑图中所示的矩形截面杆，中部受一对夹持力 P 作用，求由 P 引起的杆长变化 $\Delta$。

![矩形截面杆](../../../../md_output/exercises/exercise_03/images/b5f5e5b7954a8c6ffc14f215ef2a345661cba28da506ed74dda5a42a0cf721c5.jpg)

**解答**：

由图可知：矩形截面杆（宽 $b$，高 $h$，长 $L$），力 $P$ 沿 $h$ 方向施加（横向压缩），而非沿杆轴向。这是一个**横向加载 + Poisson 效应**问题。

**Step 1：应力状态**

力 $P$ 作用于杆的顶面和底面（面积 $b \times L$），沿 $h$ 方向压缩。

$$\sigma_h = -\frac{P}{bL} \quad \text{（压应力）}$$

杆两端自由，轴向无外力：$\sigma_L = 0$，$\sigma_b = 0$

**Step 2：应变分析（广义 Hooke 定律）**

仅 $\sigma_h \neq 0$，由 Poisson 效应：

$$\varepsilon_L = -\frac{\nu \sigma_h}{E} = \frac{\nu P}{EbL} \quad \text{（$L$ 方向伸长）}$$

**Step 3：长度变化**

$$\boxed{\Delta = \varepsilon_L \cdot L = \frac{\nu P}{Eb}}$$

---

## 2. 求泛函 $Q[y] = \int_0^1 (y'^2 - 4y^2)dx$ 的极值函数

**题目**：求泛函的极值函数。分别用直接法（Euler 方程法）和间接法（Ritz 法，取三项）。$Q[y] = \int_{0}^{1}(y'^{2} - 4y^{2})dx$，$y(0) = 0, y(1) = 1$

**解答**：

$$F = y'^2 - 4y^2$$

$$\frac{\partial F}{\partial y} = -8y$$

$$\frac{\partial F}{\partial y'} = 2y'$$

Euler 方程：

$$-8y - \frac{d}{dx}(2y') = 0$$

$$-8y - 2y'' = 0$$

$$y'' + 4y = 0$$

特征方程：$r^2 + 4 = 0$，$r = \pm 2i$

通解：

$$y = C_1\cos 2x + C_2\sin 2x$$

**代入边界条件**：

$y(0) = 0$：$C_1 = 0$

$y(1) = 1$：$C_2\sin 2 = 1 \Rightarrow C_2 = \frac{1}{\sin 2}$

**精确解**：

$$\boxed{y_{\text{exact}} = \frac{\sin 2x}{\sin 2}}$$

### 解法 B：Ritz 法（三项近似）

> 参考公式 3.16：非齐次 BC 时令 $u_n = u_0 + \sum a_i\varphi_i$，其中 $u_0$ 满足非齐次 BC，$\varphi_i$ 满足齐次 BC。

**Step 1：处理非齐次边界条件**

边界条件 $y(0)=0$（齐次），$y(1)=1$（非齐次）。

取 $u_0 = x$，满足 $u_0(0)=0$，$u_0(1)=1$。

令 $y = x + w$，其中 $w(0)=w(1)=0$（齐次 BC）。

**Step 2：选取基函数**

取三项近似，基函数需满足 $w(0)=w(1)=0$：

$$w \approx a_1\varphi_1 + a_2\varphi_2 + a_3\varphi_3$$

$$\varphi_1 = x(1-x),\quad \varphi_2 = x^2(1-x),\quad \varphi_3 = x^3(1-x)$$

**Step 3：代入泛函**

$$Q = \int_0^1 (y'^2 - 4y^2)\,dx$$

代入 $y = x + w$：

$$y' = 1 + w'$$

$$Q = \int_0^1 \left[(1+w')^2 - 4(x+w)^2\right]dx$$

展开：

$$Q = \int_0^1 \left[1 + 2w' + w'^2 - 4x^2 - 8xw - 4w^2\right]dx$$

**Step 4：代入 $w = a_1\varphi_1 + a_2\varphi_2 + a_3\varphi_3$，计算 $Q(a_1,a_2,a_3)$**

$$w' = a_1(1-2x) + a_2(2x-3x^2) + a_3(3x^2-4x^3)$$

计算各积分项（利用 $\int_0^1 x^m(1-x)^n dx = \frac{m!n!}{(m+n+1)!}$）：

$$\int_0^1 w'^2 dx = \sum_{i,j} a_i a_j \int_0^1 \varphi_i'\varphi_j' dx$$

$$\int_0^1 w^2 dx = \sum_{i,j} a_i a_j \int_0^1 \varphi_i\varphi_j dx$$

$$\int_0^1 xw\,dx = \sum_i a_i \int_0^1 x\varphi_i\,dx$$

$$\int_0^1 w'\,dx = \sum_i a_i \int_0^1 \varphi_i'\,dx = \sum_i a_i [\varphi_i]_0^1 = 0$$

**Step 5：列方程 $\frac{\partial Q}{\partial a_i} = 0$**

$$\frac{\partial Q}{\partial a_i} = 2\sum_j a_j \int_0^1 \varphi_j'\varphi_i' dx - 8\sum_j a_j \int_0^1 \varphi_j\varphi_i dx - 8\int_0^1 x\varphi_i\,dx = 0$$

写成矩阵形式：

$$\mathbf{K}\mathbf{a} = \mathbf{b}$$

其中 $K_{ij} = \int_0^1 (\varphi_i'\varphi_j' - 4\varphi_i\varphi_j)\,dx$，$b_i = 4\int_0^1 x\varphi_i\,dx$

**Step 6：计算系数并求解**

利用积分公式计算各元素：

$$K_{11} = \int_0^1 [(1-2x)^2 - 4x^2(1-x)^2]dx = \frac{1}{3} - \frac{4}{30} = \frac{3}{10}$$

$$K_{12} = \int_0^1 [(1-2x)(2x-3x^2) - 4x^3(1-x)^2]dx = \frac{1}{10} - \frac{4}{60} = \frac{1}{15}$$

类似计算其他元素，解 $3\times3$ 方程组得 $a_1, a_2, a_3$。

**Step 7：回代**

$$y = x + a_1 x(1-x) + a_2 x^2(1-x) + a_3 x^3(1-x)$$

---

## 3. 求泛函的欧拉方程和自然边界条件

**题目**：求下列泛函的欧拉方程和自然边界条件：$Q[y] = \frac{1}{2}\int_{x_0}^{x_1}[p(x)y''^2 + q(x)y'^2 + r(x)y^2 - 2s(x)y]dx$

**解答**：

$$F = \frac{1}{2}[p y''^2 + q y'^2 + r y^2 - 2sy]$$

计算偏导：

$$\frac{\partial F}{\partial y} = ry - s$$

$$\frac{\partial F}{\partial y'} = qy'$$

$$\frac{\partial F}{\partial y''} = py''$$

Euler 方程（含四阶导数）：

$$\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) + \frac{d^2}{dx^2}\left(\frac{\partial F}{\partial y''}\right) = 0$$

$$ry - s - \frac{d}{dx}(qy') + \frac{d^2}{dx^2}(py'') = 0$$

展开：

$$\boxed{\frac{d^2}{dx^2}(py'') - \frac{d}{dx}(qy') + ry = s}$$

**自然边界条件**：

在 $x = x_0$ 和 $x = x_1$ 处：

- $\frac{\partial F}{\partial y''} = py'' = 0$（弯矩为零）或 $y$ 给定
- $\frac{d}{dx}\left(\frac{\partial F}{\partial y''}\right) - \frac{\partial F}{\partial y'} = (py'')' - qy' = 0$（剪力为零）或 $y'$ 给定

---

## 4. 悬臂梁问题

**题目**：受均布荷载的悬臂梁。(1) 用挠度方程求出精确解；(2) 写出两种以上的许可位移场（试函数）；(3) 分别用最小势能原理（Rayleigh-Ritz 法）、Galerkin 加权残值法、残值最小二乘法求挠度曲线。

![悬臂梁](../../../../md_output/exercises/exercise_03/images/bbfeeb9f78e116822479211dbf49ea83f3864faeeaf88e9812250ac4ed55cbe1.jpg)

**解答**：

### (1) 精确解

控制方程：$EI\frac{d^4 w}{dx^4} = p$

边界条件（固定端 $x=0$，自由端 $x=l$）：

- $w(0) = 0$（位移为零）
- $w'(0) = 0$（转角为零）
- $w''(l) = 0$（弯矩为零）
- $w'''(l) = 0$（剪力为零）

积分四次：

$$w = \frac{p}{24EI}x^4 + C_1 x^3 + C_2 x^2 + C_3 x + C_4$$

由 $w(0)=0$：$C_4=0$

由 $w'(0)=0$：$C_3=0$

由 $w'''(l)=0$：$\frac{p}{EI}l + 6C_1 = 0 \Rightarrow C_1 = -\frac{pl}{6EI}$

由 $w''(l)=0$：$\frac{p}{2EI}l^2 + 6C_1 l + 2C_2 = 0 \Rightarrow C_2 = \frac{pl^2}{4EI}$

**精确解**：

$$\boxed{w(x) = \frac{p}{24EI}(x^4 - 4lx^3 + 6l^2x^2)}$$

---

### (2) 许可位移场（试函数）

> **核心区别**：不同方法对试函数的要求不同，需准备两套。

**用于 Ritz 法**：只需满足**本质边界条件**（$w(0)=0$，$w'(0)=0$）

$$\varphi_1 = x^2, \quad \varphi_2 = x^3, \quad \varphi_3 = x^4$$

验证：$\varphi_i(0)=0$ ✓，$\varphi_i'(0)=0$ ✓

**用于 Galerkin / 最小二乘法**：需满足**全部边界条件**（本质 + 自然 BC）

检查各单项：

| 基函数 | $w''(l)$ | $w'''(l)$ | 满足全部 BC？ |
|--------|----------|-----------|--------------|
| $x^2$ | $2$ | $0$ | ✗ |
| $x^3$ | $6l$ | $6$ | ✗ |
| $x^4$ | $12l^2$ | $24l$ | ✗ |

需要构造满足全部 BC 的试函数：

$$\varphi_1 = x^4 - 4lx^3 + 6l^2x^2$$

验证：$\varphi_1''(l) = 12l^2 - 24l^2 + 12l^2 = 0$ ✓，$\varphi_1'''(l) = 24l - 24l = 0$ ✓

---

### (3) 三种近似解法

#### 方法 A：最小势能原理（Rayleigh-Ritz 法）

Ritz 法只需满足**本质边界条件**（$w(0)=0, w'(0)=0$），取一项近似 $w = a_1 x^2$。

**势能泛函**：

$$\Pi = \int_0^l \left[\frac{EI}{2}(w'')^2 - pw\right] dx$$

**代入** $w = a_1 x^2$，$w'' = 2a_1$：

$$\Pi = \int_0^l \left[2EI a_1^2 - pa_1 x^2\right] dx = 2EI a_1^2 l - pa_1 \cdot \frac{l^3}{3}$$

**极值条件**：

$$\frac{\partial \Pi}{\partial a_1} = 4EI a_1 l - \frac{pl^3}{3} = 0 \quad \Rightarrow \quad \boxed{a_1 = \frac{pl^2}{12EI}}$$

**近似解**：$w_A = \frac{pl^2}{12EI}x^2$，$w_A(l) = \frac{pl^4}{12EI}$（误差 33%，偏小）

---

#### 方法 B：Galerkin 加权残值法

Galerkin 法要求试函数满足**全部边界条件**（包括力边界 $w''(l)=0, w'''(l)=0$）。

取第 (2) 问构造的试函数：$\varphi_1 = x^4 - 4lx^3 + 6l^2x^2$

验证：$\varphi_1''(l) = 12l^2 - 24l^2 + 12l^2 = 0$ ✓，$\varphi_1'''(l) = 24l - 24l = 0$ ✓

**控制方程**：$EI w'''' - p = 0$

代入 $w = a_1\varphi_1$：$w'''' = 24a_1$，残量 $R = 24EI a_1 - p$

**Galerkin 条件**：

$$\int_0^l R \cdot \varphi_1\,dx = (24EI a_1 - p)\int_0^l (x^4 - 4lx^3 + 6l^2x^2)dx = 0$$

计算积分：$\int_0^l (x^4 - 4lx^3 + 6l^2x^2)dx = \frac{l^5}{5} - l^4\cdot l + 2l^2\cdot l^3 = \frac{6l^5}{5} \neq 0$

$$\boxed{a_1 = \frac{p}{24EI}}$$

**近似解**：$w_B = \frac{p}{24EI}(x^4 - 4lx^3 + 6l^2x^2)$，$w_B(l) = \frac{p}{24EI}\cdot 3l^4 = \frac{pl^4}{8EI}$

**与精确解完全一致！** 因为试函数 $\varphi_1$ 恰好就是精确解的形式。

---

#### 方法 C：残值最小二乘法

同样取满足全部 BC 的试函数：$w = a_1\varphi_1 = a_1(x^4 - 4lx^3 + 6l^2x^2)$

**残差**：$R = EIw'''' - p = 24EI a_1 - p$

**最小二乘目标**：

$$J = \int_0^l R^2\,dx = (24EI a_1 - p)^2 l$$

$$\frac{\partial J}{\partial a_1} = 2(24EI a_1 - p)\cdot 24EI\cdot l = 0$$

$$\boxed{a_1 = \frac{p}{24EI}}$$

**近似解**：$w_C = \frac{p}{24EI}(x^4 - 4lx^3 + 6l^2x^2) = w_{\text{exact}}$

---

### 精确解与近似解对比

| 方法 | 试函数 | 近似解 | $w(l)$ | 误差 |
|------|--------|--------|--------|------|
| 精确解 | — | $\frac{p}{24EI}(x^4-4lx^3+6l^2x^2)$ | $\frac{pl^4}{8EI}$ | — |
| Rayleigh-Ritz | $x^2$（仅本质 BC） | $\frac{pl^2}{12EI}x^2$ | $\frac{pl^4}{12EI}$ | 33% 偏小 |
| Galerkin | $\varphi_1$（全部 BC） | $\frac{p}{24EI}\varphi_1$ | $\frac{pl^4}{8EI}$ | **0%** |
| 最小二乘 | $\varphi_1$（全部 BC） | $\frac{p}{24EI}\varphi_1$ | $\frac{pl^4}{8EI}$ | **0%** |

> **关键结论**：
> - Galerkin 法和最小二乘法用满足全部 BC 的试函数后，精度远超 Ritz 法
> - 试函数的选择比方法本身更重要——$\varphi_1$ 恰好是精确解的形式，所以一阶近似就得到精确解
> - Rayleigh-Ritz 法只需本质 BC，试函数容易选但精度低；Galerkin 法需全部 BC，试函数难选但精度高

---

## 5. 加权残值法求解

**题目**：设 1D 问题：$\frac{d^2\varphi}{dx^2} + \varphi + x = 0$（$0 \leq x \leq 1$），边界条件 $\varphi(0) = \varphi(1) = 0$。试函数 $\varphi = c_1\varphi_1 + c_2\varphi_2$，其中 $\varphi_1 = x(1-x)$，$\varphi_2 = x^2(1-x)$。试用 Galerkin 方法、最小二乘法、Rayleigh-Ritz 方法求解。

**解答**：

**残差**：

$$R = \frac{d^2\varphi}{dx^2} + \varphi + x$$

计算导数：

$$\varphi_1'' = -2, \quad \varphi_2'' = 2 - 6x$$

$$R = c_1(-2 + x - x^2) + c_2(2 - 6x + x^2 - x^3) + x$$

**Galerkin 条件**：

$$\int_0^1 R \cdot \varphi_1 \, dx = 0$$

$$\int_0^1 R \cdot \varphi_2 \, dx = 0$$

代入计算得两个方程，解方程组得 $c_1, c_2$。

### (2) 最小二乘法

**目标**：最小化 $J = \int_0^1 R^2 dx$

$$\frac{\partial J}{\partial c_1} = 2\int_0^1 R \frac{\partial R}{\partial c_1} dx = 0$$

$$\frac{\partial J}{\partial c_2} = 2\int_0^1 R \frac{\partial R}{\partial c_2} dx = 0$$

计算得方程组，解得 $c_1, c_2$。

### (3) Rayleigh-Ritz 方法

**等效泛函**：

$$Q[\varphi] = \int_0^1 \left[\frac{1}{2}(\varphi')^2 - \frac{1}{2}\varphi^2 - x\varphi\right] dx$$

代入 $\varphi = c_1\varphi_1 + c_2\varphi_2$，对 $c_1, c_2$ 求偏导并令为零。

---

## 6. Galerkin 方法推导二维热传导有限元单元

**题目**：使用 Galerkin 方法推导二维稳定状态热传导问题的一个四节点有限元单元。

**解答**：

$$\int_\Omega \nabla T \cdot \nabla N_i \, dxdy = \int_\Omega f N_i \, dxdy + \int_{\Gamma_q} \bar{q} N_i \, ds$$

**四节点四边形单元**：

形函数（等参元）：

$$N_i = \frac{1}{4}(1 + \xi\xi_i)(1 + \eta\eta_i), \quad i = 1,2,3,4$$

**单元刚度矩阵**：

$$K_{ij}^e = \int_{-1}^{1}\int_{-1}^{1} k \left(\frac{\partial N_i}{\partial x}\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}\frac{\partial N_j}{\partial y}\right) |J| \, d\xi d\eta$$

其中 $J$ 是 Jacobi 矩阵：

$$J = \begin{bmatrix} \frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\ \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta} \end{bmatrix}$$

**单元载荷向量**：

$$f_i^e = \int_{-1}^{1}\int_{-1}^{1} f N_i |J| \, d\xi d\eta$$

组装总体刚度矩阵和载荷向量，引入边界条件后求解。
