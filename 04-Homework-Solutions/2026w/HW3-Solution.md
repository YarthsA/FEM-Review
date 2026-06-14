# Variation Method & FEA Course — Homework 3 参考答案

> 交大 FEM 课程 · Galerkin 法与 FEM 应用专题

---

## 1. 矩形截面杆在夹持力下的长度变化

**题目**：矩形截面杆（宽 $b$，高 $h$，长 $L$）中部受一对夹持力 $P$ 作用，求杆的长度变化 $\Delta$。

> ⚠️ 注意：由图可知，力 $P$ 沿 $h$ 方向施加（横向压缩），而非沿杆轴向。因此这是一个**横向加载 + Poisson 效应**的问题。

### 解答

**Step 1**：物理模型与应力状态

杆的截面面积：$A = b \cdot h$

力 $P$ 作用于杆的顶面和底面（面积 $b \times L$），沿 $h$ 方向压缩。

$h$ 方向的正应力：
$$\sigma_h = -\frac{P}{bL} \quad \text{（压应力）}$$

杆两端自由，轴向（$L$ 方向）无外力：
$$\sigma_L = 0, \quad \sigma_b = 0$$

**Step 2**：应变分析（广义 Hooke 定律）

对于**单轴应力状态**（仅 $\sigma_h \neq 0$），三个方向的应变为：

$$\varepsilon_h = \frac{\sigma_h}{E} = -\frac{P}{EbL} \quad \text{（$h$ 方向压缩）}$$

$$\varepsilon_L = -\frac{\nu\,\sigma_h}{E} = \frac{\nu P}{EbL} \quad \text{（$L$ 方向因 Poisson 效应伸长）}$$

$$\varepsilon_b = -\frac{\nu\,\sigma_h}{E} = \frac{\nu P}{EbL} \quad \text{（$b$ 方向因 Poisson 效应伸长）}$$

其中 $\nu$ 为材料的 **Poisson 比**。

**Step 3**：长度变化

杆的长度 $L$ 方向的伸长量：

$$\boxed{\Delta = \varepsilon_L \cdot L = \frac{\nu P }{Eb}}$$


---

## 2. 悬臂梁 Galerkin 法试探函数检验

**题目**：证明 $w = a\left(1 - \cos\frac{\pi x}{2l}\right)$ 不是悬臂梁 Galerkin 法的合法试探函数。

### 解答

**Step 1**：悬臂梁的边界条件

悬臂梁在 $x=0$ 处固定，$x=l$ 处自由。边界条件：

固定端（$x=0$）：
$$w(0) = 0 \quad\text{(位移为零)}$$
$$w'(0) = 0 \quad\text{(转角为零)}$$

自由端（$x=l$）：
$$M(l) = 0 \Rightarrow w''(l) = 0 \quad\text{(弯矩为零)}$$
$$V(l) = 0 \Rightarrow w'''(l) = 0 \quad\text{(剪力为零)}$$

**Step 2**：检验边界条件

检验函数 $w(x) = a\left(1 - \cos\frac{\pi x}{2l}\right)$ 的各阶导数：

$$w'(x) = a \cdot \frac{\pi}{2l} \cdot \sin\frac{\pi x}{2l}$$
$$w''(x) = a \cdot \left(\frac{\pi}{2l}\right)^2 \cdot \cos\frac{\pi x}{2l}$$
$$w'''(x) = -a \cdot \left(\frac{\pi}{2l}\right)^3 \cdot \sin\frac{\pi x}{2l}$$

在固定端 $x=0$：
$$w(0) = a(1 - \cos 0) = a(1-1) = 0 \quad \checkmark$$
$$w'(0) = a \cdot \frac{\pi}{2l} \cdot 0 = 0 \quad \checkmark$$

**位移和转角边界条件均满足。**

在自由端 $x=l$：
$$w''(l) = a\left(\frac{\pi}{2l}\right)^2 \cos\frac{\pi}{2} = 0 \quad \checkmark$$
$$w'''(l) = -a\left(\frac{\pi}{2l}\right)^3 \sin\frac{\pi}{2} = -a\left(\frac{\pi}{2l}\right)^3 \neq 0 \quad \times$$

**Step 3**：Galekin 法对试探函数的要求

在经典 Galerkin 法中（强形式），试探函数必须满足**所有边界条件**（包括本质边界条件和自然边界条件）。这是因为 Galerkin 法要求：
- 试探函数 $\bar{w}$ 属于容许函数空间，其中所有函数均满足问题的全部边界条件

对于 Euler-Bernoulli 梁控制方程 $EI\,w'''' = q$ 的 Galerkin 弱形式：
$$\int_0^l EI\,w''\varphi''\,dx = \int_0^l q\varphi\,dx$$

试探函数可只满足**本质边界条件**（$w(0)=0, w'(0)=0$），自然边界条件（$w''(l)=0, w'''(l)=0$）在弱形式中自动包含。

对于**强形式的经典 Galerkin 法**（直接处理微分方程余量），试探函数需满足**全部**边界条件。

**Step 4**：结论

$$w'''(l) = -a\left(\frac{\pi}{2l}\right)^3 \neq 0$$

自由端剪力不为零，违背了悬臂梁自由端 $V(l) = 0$ 的自然边界条件。

因此，$w = a\left(1 - \cos\frac{\pi x}{2l}\right)$ **不是**悬臂梁 Galerkin 法的合法试探函数。■

---

## 3. 弹性地基梁 — 变分法与 Galerkin 法求解

**题目**：用变分法导出弹性地基梁的微分方程和边界条件，并用 Galerkin 法的一阶近似求解。

### 3.1 变分法推导控制方程

**Step 1**：体系总势能

弹性地基梁（Winkler 地基）的总势能包括三部分：
- 梁的弯曲应变能：$U_b = \int_0^l \frac{1}{2}EI (w'')^2 dx$
- 弹性地基应变能：$U_f = \int_0^l \frac{1}{2} k w^2 dx$（$k$ 为地基刚度系数）
- 外荷载势能：$V = -\int_0^l p\,w\,dx$

总势能泛函：
$$\Pi[w] = \int_0^l \left[\frac{1}{2}EI(w'')^2 + \frac{1}{2}kw^2 - pw\right]dx$$

**Step 2**：变分

令 $\delta\Pi = 0$：
$$\delta\Pi = \int_0^l \left[EI\,w''\delta w'' + kw\,\delta w - p\,\delta w\right]dx = 0$$

**Step 3**：分部积分

对 $EI\,w''\delta w''$ 项分部积分：
$$\int_0^l EI\,w''\delta w''\,dx = \left[EI\,w''\delta w'\right]_0^l - \int_0^l EI\,w'''\delta w'\,dx$$

再分部积分一次：
$$= \left[EI\,w''\delta w'\right]_0^l - \left[EI\,w'''\delta w\right]_0^l + \int_0^l EI\,w''''\delta w\,dx$$

**Step 4**：代入变分方程

$$[EI\,w''\delta w']_0^l - [EI\,w'''\delta w]_0^l + \int_0^l (EI\,w'''' + kw - p)\delta w\,dx = 0$$

由 $\delta w$ 的任意性，得**控制微分方程**：
$$\boxed{EI\frac{d^4w}{dx^4} + kw = p(x)} \quad (0 < x < l)$$

**Step 5**：边界条件

由分部积分产生的边界项：
$$[EI\,w''\delta w']_0^l - [EI\,w'''\delta w]_0^l = 0$$

在每一端点，以下组合之一须成立：
- 若 $w$ 给定 → $\delta w = 0$（本质边界条件）
- 若 $w$ 未给定 → $EI\,w''' = 0$（自然边界条件，剪力为零）
- 若 $w'$ 给定 → $\delta w' = 0$（本质边界条件）
- 若 $w'$ 未给定 → $EI\,w'' = 0$（自然边界条件，弯矩为零）

### 3.2 一阶 Galerkin 近似

**Step 1**：设定试函数

对于简支梁（$w(0)=w(l)=0, w''(0)=w''(l)=0$），取一阶试函数：
$$w_1(x) = a_1 \sin\frac{\pi x}{l}$$

验证边界条件：$w_1(0) = a_1 \sin 0 = 0 \checkmark$, $w_1(l) = a_1 \sin\pi = 0 \checkmark$, $w_1''(0) = -a_1(\pi/l)^2\sin 0 = 0 \checkmark$, $w_1''(l) = 0 \checkmark$。

**Step 2**：Galerkin 加权残量方程

控制方程的余量：
$$R(x) = EI\frac{d^4w_1}{dx^4} + kw_1 - p(x)$$

Galerkin 条件：余量与权函数（取基函数本身）正交
$$\int_0^l R(x)\cdot\sin\frac{\pi x}{l}\,dx = 0$$

**Step 3**：代入计算

$w_1(x) = a_1\sin\frac{\pi x}{l}$ 的各阶导数：
$$w_1' = a_1\frac{\pi}{l}\cos\frac{\pi x}{l}, \quad w_1'' = -a_1\left(\frac{\pi}{l}\right)^2\sin\frac{\pi x}{l}$$
$$w_1''' = -a_1\left(\frac{\pi}{l}\right)^3\cos\frac{\pi x}{l}, \quad w_1'''' = a_1\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l}$$

代入 Galerkin 方程：
$$\int_0^l \left[EI\,a_1\left(\frac{\pi}{l}\right)^4\sin\frac{\pi x}{l} + k\,a_1\sin\frac{\pi x}{l} - p(x)\right]\sin\frac{\pi x}{l}\,dx = 0$$

**Step 4**：求解 $a_1$

$$a_1\left[EI\left(\frac{\pi}{l}\right)^4 + k\right]\int_0^l \sin^2\frac{\pi x}{l}\,dx = \int_0^l p(x)\sin\frac{\pi x}{l}\,dx$$

利用 $\int_0^l \sin^2\frac{\pi x}{l}dx = \frac{l}{2}$：

$$a_1 = \frac{2}{l}\cdot\frac{\int_0^l p(x)\sin\frac{\pi x}{l}\,dx}{EI\left(\frac{\pi}{l}\right)^4 + k}$$

**Step 5**：若 $p(x) = p_0$（均布荷载）

则：
$$\int_0^l p_0\sin\frac{\pi x}{l}dx = p_0\cdot\frac{2l}{\pi}$$

$$a_1 = \frac{2}{l}\cdot\frac{p_0\cdot 2l/\pi}{EI(\pi/l)^4 + k} = \frac{4p_0}{\pi\left[EI(\pi/l)^4 + k\right]}$$

**最终解答**：
$$\boxed{w(x) = \frac{4p_0}{\pi\left[EI(\pi/l)^4 + k\right]}\sin\frac{\pi x}{l}}$$

---

## 4. 两节点 Euler-Bernoulli 梁单元形函数

**题目**：建立两节点 Euler-Bernoulli 梁单元的形函数（Hermite 插值）。

### 解答

**Step 1**：单元自由度

两节点 Euler-Bernoulli 梁单元，每个节点有 2 个自由度：
- $u_1$：节点 1 的挠度 $w_1$
- $\theta_1$：节点 1 的转角 $\frac{dw}{dx}\big|_1$
- $u_2$：节点 2 的挠度 $w_2$
- $\theta_2$：节点 2 的转角 $\frac{dw}{dx}\big|_2$

共 4 个自由度，需要 4 个条件确定三次多项式。

**Step 2**：位移场假设

每节点 2 个自由度 → 采用**三次 Hermite 插值**（$C^1$ 连续）：
$$w(x) = a_0 + a_1x + a_2x^2 + a_3x^3$$

转角：
$$\theta(x) = \frac{dw}{dx} = a_1 + 2a_2x + 3a_3x^2$$

**Step 3**：引入自然坐标

设单元长度 $L = x_2 - x_1$，定义自然坐标 $\xi = \frac{x - x_1}{L}$，$\xi \in [0, 1]$：

$$w(\xi) = c_0 + c_1\xi + c_2\xi^2 + c_3\xi^3$$

节点条件（$\xi=0$ 和 $\xi=1$）：
$$w(0) = w_1, \quad \theta(0) = \frac{1}{L}\frac{dw}{d\xi}\bigg|_0 = \theta_1$$
$$w(1) = w_2, \quad \theta(1) = \frac{1}{L}\frac{dw}{d\xi}\bigg|_1 = \theta_2$$

**Step 4**：系数求解

由节点条件：
$$\begin{pmatrix}
w_1 \\ L\theta_1 \\ w_2 \\ L\theta_2
\end{pmatrix}
= \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
1 & 1 & 1 & 1 \\
0 & 1 & 2 & 3
\end{pmatrix}
\begin{pmatrix}
c_0 \\ c_1 \\ c_2 \\ c_3
\end{pmatrix}$$

求逆得系数：
$$\begin{pmatrix}
c_0 \\ c_1 \\ c_2 \\ c_3
\end{pmatrix}
= \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
-3 & -2 & 3 & -1 \\
2 & 1 & -2 & 1
\end{pmatrix}
\begin{pmatrix}
w_1 \\ L\theta_1 \\ w_2 \\ L\theta_2
\end{pmatrix}$$

**Step 5**：形函数表达式

将 $w(\xi) = \sum_{i=1}^4 N_i(\xi) d_i$ 展开，其中 $d = (w_1, \theta_1, w_2, \theta_2)^T$。

代入各系数得 Hermite 形函数：
$$\boxed{N_1(\xi) = 1 - 3\xi^2 + 2\xi^3 \quad \text{(对应 $w_1$)}}$$
$$\boxed{N_2(\xi) = L(\xi - 2\xi^2 + \xi^3) \quad \text{(对应 $\theta_1$)}}$$
$$\boxed{N_3(\xi) = 3\xi^2 - 2\xi^3 \quad \text{(对应 $w_2$)}}$$
$$\boxed{N_4(\xi) = L(-\xi^2 + \xi^3) \quad \text{(对应 $\theta_2$)}}$$

其中 $\xi = \frac{x-x_1}{L} \in [0, 1]$。

**Step 6**：形函数性质验证

| 性质 | 验证 |
|------|------|
| $N_1(0)=1, N_1(1)=0$ | 节点 1 的挠度在自身节点为 1，在节点 2 为 0 ✓ |
| $N_1'(0)=0, N_1'(1)=0$ | 对转角无贡献 ✓ |
| $N_2(0)=0, N_2(1)=0$ | 对挠度无贡献 ✓ |
| $N_2'(0)=L, N_2'(1)=0$ | 节点 1 的转角在自身为 $L$（因 $\frac{dw}{dx}= \frac{1}{L}\frac{dw}{d\xi}$）✓ |
| 完备性：$\sum N_i = 1$ | 刚体位移（平移）模式 ✓ |
| $-\sum N_i x_i + \sum N_i'L = 0$... | 刚体转动模式 ✓ |

**Step 7**：位移场插值表示

$$w(\xi) = N_1(\xi)w_1 + N_2(\xi)\theta_1 + N_3(\xi)w_2 + N_4(\xi)\theta_2$$

或写成矩阵形式：
$$w = \begin{pmatrix} N_1 & N_2 & N_3 & N_4 \end{pmatrix} \begin{pmatrix} w_1 & \theta_1 & w_2 & \theta_2 \end{pmatrix}^T = [N]\{\delta\}_e$$

**Step 8**：单元刚度矩阵

基于 Euler-Bernoulli 梁理论（$M = EI\frac{d^2w}{dx^2}$），单元刚度矩阵：
$$[k]_e = \int_0^L EI [B]^T[B]\,dx$$

其中 $B = \frac{d^2N}{dx^2} = \frac{1}{L^2}\frac{d^2N}{d\xi^2}$。

标准结果：
$$[k]_e = \frac{EI}{L^3}\begin{pmatrix}
12 & 6L & -12 & 6L \\
6L & 4L^2 & -6L & 2L^2 \\
-12 & -6L & 12 & -6L \\
6L & 2L^2 & -6L & 4L^2
\end{pmatrix}$$

■
