# 第3章：变分法基础（Variational Calculus）

> **对应 PDF**：[`Chapter 3 Variation theory and applications-1.pdf`](../06-References/pdfs-originals/Chapter%203%20Variation%20theory%20and%20applications-1.pdf) · [`有限元复习.pdf`](../06-References/pdfs-originals/有限元复习.pdf) §4
> **相关作业**：[HW2 Q1-Q4](../04-Homework-Solutions/2026w/HW2-Problem.md)（全部变分法题）· [HW3 Q3](../04-Homework-Solutions/2026w/HW3-Problem.md)（弹性地基梁）
> **前置知识**：微积分（微分定义、分部积分、常微分方程）、线性代数

**变分法 (Variational Method)** — 泛函分析的一个分支，研究**泛函极值**问题。

### 核心问题

从合适的函数集中选择某个函数 $y(x)$，使某个积分型泛函 $Q[y]$ 达到极值：
$$
Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx \to \text{extremum}
$$

### 变分法发展简史

- **1696** — Johann Bernoulli 提出**最速降线问题 (Brachistochrone)**：质点在重力作用下从 A 到 B 沿哪条曲线时间最短？
- **1736** — Euler 给出一般解法 → **Euler 方程**
- **1755** — Lagrange 改进简化，变分法正式建立
- **1908** — Ritz 提出直接法（近似解法）
- **1960** — FEM 诞生，变分法成为其理论基石

---

## 2. 泛函 (Functional)

### 定义

泛函是"函数的函数"——从函数空间到实数域的映射：
$$
Q = Q[f(x)]
$$

与函数的本质区别：
- **函数** $y = f(x)$：依赖于自变量 $x$ 在某个点的取值
- **泛函** $Q[f]$：依赖于函数 $f$ 的**整体形状**（整个函数关系）

**例**：$Q[f] = \int_0^1 f(x)\,dx$ 输入一个函数 $f$，输出一个标量（曲线下的面积）。

### 线性泛函

严格定义：$Q$ 是线性泛函当且仅当同时满足：
1. **齐次性** (Homogeneity)：$Q[cy(x)] = c\,Q[y(x)]$
2. **可加性** (Additivity)：$Q[y_1(x) + y_2(x)] = Q[y_1(x)] + Q[y_2(x)]$

统一表达：
$$
Q[c_1y_1 + c_2y_2] = c_1Q[y_1] + c_2Q[y_2]
$$

**例**：$Q[y] = \int_a^b y(x)dx$ — ✅ **是**线性泛函
**例**：$Q[y] = \int_a^b y^2(x)dx$ — ❌ **不是**线性泛函

> 证明后者不是线性：
> $$Q[c_1y_1 + c_2y_2] = \int (c_1y_1 + c_2y_2)^2dx = \int (c_1^2y_1^2 + c_2^2y_2^2 + 2c_1c_2y_1y_2)dx$$
> $$c_1Q[y_1] + c_2Q[y_2] = c_1\int y_1^2dx + c_2\int y_2^2dx$$
> 两者不等（交叉项 $2c_1c_2y_1y_2$ 不消失）→ 不是线性泛函。

---

## 3. 变分 (Variation)

### 3.1 自变函数的变分

两个函数之间的差：
$$
\delta y = y(x) - y_1(x),\quad \delta y' = y'(x) - y_1'(x)
$$

**K 阶接近度 (K-th order proximity)**：$\delta y, \delta y', \ldots, \delta y^{(k)}$ 都很小 ⟺ 两条曲线有 K 阶接近度。

### 3.2 函数微分回顾（两种定义）

**定义1（常规）**：增量分解为线性主部 + 高阶小量
$$
\Delta y = A(x)\Delta x + o(\Delta x),\quad dy = A(x)dx = y'(x)dx
$$

**定义2（Lagrange）**：引入参数 $\varepsilon$ 求导
$$
\frac{\partial}{\partial\varepsilon}y(x + \varepsilon\Delta x)\bigg|_{\varepsilon=0} = y'(x)\Delta x = dy(x)
$$

### 3.3 泛函的变分

模仿函数的微分，泛函的变分也有两种定义。

**定义1（常规）**：将增量分解为线性泛函 + 高阶小量
$$
\Delta Q = Q[y + \delta y] - Q[y] = T[y, \delta y] + o(\delta y)
$$
其中 $T$ 对 $\delta y$ 是线性的，则 $\delta Q = T[y, \delta y]$。

**例**：$Q[y] = \int_a^b y^2dx$
$$
\Delta Q = \int_a^b (y+\delta y)^2dx - \int_a^b y^2dx = \int_a^b 2y\,\delta y\,dx + \int_a^b (\delta y)^2dx
$$
第一项是 $\delta y$ 的线性泛函 → $\delta Q = \int_a^b 2y\,\delta y\,dx$

**定义2（Lagrange）**：计算利器
$$
\delta Q = \frac{\partial}{\partial\varepsilon} Q[y + \varepsilon\delta y]\bigg|_{\varepsilon=0}
$$

> **等价性证明**：$Q[y+\varepsilon\delta y]$ 是 $\varepsilon$ 的函数，在 $\varepsilon=0$ 处 Taylor 展开的一次项系数就是 $\delta Q$。

### 3.4 变分记号运算规则（重要！）

微分算子和变分算子可交换顺序：
$$
\frac{d}{dx}(\delta y) = \delta(y')
$$

这一性质在推导 Euler 方程时分部积分时非常方便。

---

## 4. 泛函极值条件

类比函数极值（$dy=0$）：
$$
\delta Q = 0 \quad \text{(泛函极值的必要条件)}
$$

### 强极值与弱极值
- **强极值**：在零阶邻域（函数值接近）内成立
- **弱极值**：在一阶邻域（函数值和一阶导数都接近）内成立

> 考试中通常只用到一阶条件 $\delta Q=0$，充分性由物理背景保证。

---

## 5. Euler 方程

### 5.1 变分法预备定理（Fundamental Lemma）

若 $f(x)$ 在 $[a,b]$ 上连续，且对任意满足 $\eta(a)=\eta(b)=0$ 的 $\eta(x)$ 都有：
$$
\int_a^b f(x)\eta(x)\,dx = 0
$$
则 $f(x) \equiv 0$ 在 $[a,b]$ 上。

> 反证法思路：若在某点 $f(x_0)>0$，构造函数 $\eta(x)$ 在该点附近为正、其余为零 → 积分 > 0 → 矛盾。

### 5.2 简单 Euler 方程（完整推导）

考虑泛函 $Q[y] = \int_{x_1}^{x_2} F(x, y, y')\,dx$，边界条件 $y(x_1)=y_1, y(x_2)=y_2$。

**Step 1**：构造邻域函数

令 $y_1(x) = y(x) + \alpha\eta(x)$，其中 $\eta(x_1)=\eta(x_2)=0$，$\alpha$ 是小参数。则：
$$
\varphi(\alpha) = Q[y + \alpha\eta] = \int_{x_1}^{x_2} F(x, y+\alpha\eta, y'+\alpha\eta')\,dx
$$

极值条件 $\varphi'(0) = 0$。

**Step 2**：求导（在积分号内求导）
$$
\varphi'(0) = \int_{x_1}^{x_2} \left[\frac{\partial F}{\partial y}\eta + \frac{\partial F}{\partial y'}\eta'\right]dx = 0
$$

**Step 3**：分部积分处理第二项
$$
\int_{x_1}^{x_2} \frac{\partial F}{\partial y'}\eta'\,dx = \left[\frac{\partial F}{\partial y'}\eta\right]_{x_1}^{x_2} - \int_{x_1}^{x_2} \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right)\eta\,dx
$$

由于 $\eta(x_1)=\eta(x_2)=0$，边界项为零，得：
$$
\int_{x_1}^{x_2} \left[\frac{\partial F}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right)\right]\eta(x)\,dx = 0
$$

**Step 4**：应用预备定理

$\eta(x)$ 是任意满足端点为零的函数，被积函数连续 → 括号内必须恒为零：
$$
\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right) = 0}
$$

这就是 **Euler-Lagrange 方程**（简称为 Euler 方程）。

### 5.3 $F$ 不显含 $x$ 时的首次积分

当 $F = F(y, y')$ 时，有首次积分：
$$
F - y'\frac{\partial F}{\partial y'} = C\quad(\text{常数})
$$

**推导**：
$$
\frac{d}{dx}(F - y'F_{y'}) = F_yy' + F_{y'}y'' - y''F_{y'} - y'(F_{y'y}y' + F_{y'y'}y'') = y'(F_y - \frac{d}{dx}F_{y'}) = 0
$$

这个公式在解最速降线等问题时非常有用。

### 5.4 简单算例

**例**：求 $Q[y]=\int_0^1 (y'^2 - 4y^2)dx$，$y(0)=0, y(1)=1$ 的极值函数。

解：$F = y'^2 - 4y^2$，$F_y = -8y$，$F_{y'} = 2y'$
Euler 方程：$-8y - \frac{d}{dx}(2y') = 0$ → $y'' + 4y = 0$
通解：$y = C_1\cos 2x + C_2\sin 2x$
代入 BC：$y(0)=0 \Rightarrow C_1=0$，$y(1)=1 \Rightarrow C_2 = 1/\sin 2$
**解**：$y(x) = \frac{\sin 2x}{\sin 2}$

**例**：求 $Q[y]=\int_0^1 (y'^2 + 4y^2 - 8xy)dx$，$y(0)=1, y(1)=2$ 的极值函数。

解：$F = y'^2 + 4y^2 - 8xy$，$F_y = 8y - 8x$，$F_{y'} = 2y'$
Euler 方程：$8y - 8x - 2y'' = 0$ → $y'' - 4y = -4x$
齐次解 $y_h = C_1e^{2x} + C_2e^{-2x}$，特解 $y_p = x$
通解 $y = C_1e^{2x} + C_2e^{-2x} + x$
代入 BC 得 $C_1 = \frac{1}{e^2+1}, C_2 = \frac{e^2}{e^2+1}$

---

## 6. Euler 方程推广

### 6.1 含二阶导数的泛函

$Q[y] = \int_a^b F(x, y, y', y'')dx$

$$\delta Q = \int (F_y\delta y + F_{y'}\delta y' + F_{y''}\delta y'')dx = 0$$

分部积分（对 $F_{y''}\delta y''$ 项需要两次分部积分）：
$$F_y - \frac{d}{dx}F_{y'} + \frac{d^2}{dx^2}F_{y''} = 0$$

### 6.2 含 $n$ 阶导数的泛函（通用形式）

$$Q[y] = \int_a^b F(x, y, y', \ldots, y^{(n)})dx$$

$$\boxed{\sum_{k=0}^n (-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right) = 0}$$

其中 $y^{(0)} = y$。规律：符号正负交替 $+ - + - \cdots$。

### 6.3 多个独立函数的泛函

$$Q[y_1, \ldots, y_n] = \int_a^b F(x, y_1, \ldots, y_n, y_1', \ldots, y_n')dx$$

对每个 $y_i$ 独立地有：
$$F_{y_i} - \frac{d}{dx}F_{y_i'} = 0 \quad (i=1,2,\ldots,n)$$

### 6.4 多元函数的泛函（多变量）

$$Q[z(x,y)] = \iint_D F(x, y, z, p, q)\,dxdy,\quad p = \frac{\partial z}{\partial x}, q = \frac{\partial z}{\partial y}$$

Euler 方程：
$$F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0$$

**例**：$Q[z] = \iint_D (z_x^2 + z_y^2)dxdy$ → $F = p^2+q^2$ → $F_z=0, F_p=2p, F_q=2q$ → $\frac{\partial^2 z}{\partial x^2} + \frac{\partial^2 z}{\partial y^2} = 0$（Laplace 方程！）

---

## 7. 边界条件分类

### 本质边界条件 (Essential BC)
边界值**预先给定**，与变分无关。变分时 $\delta y$ 在边界上为零。
- 例：固定端 $y(a)=y_a$，简支端 $w(0)=0$

### 自然边界条件 (Natural BC)
由泛函极值条件**自动导出**的边界条件。
- 当边界值未给定时，分部积分产生的边界项必须单独为零：
  $$\left.\frac{\partial F}{\partial y'}\right|_{x=a} = 0 \quad \text{或} \quad \left.\frac{\partial F}{\partial y'}\right|_{x=b} = 0$$
- 例：梁的自由端 $EIw''(l)=0$（弯矩为零）

### 无条件极值 vs 条件极值
- **无条件极值**：函数满足边界条件外无其他约束
- **条件极值**：附加等周条件 $\int \varphi\,dx = \alpha$ 等 → 用 **Lagrange 乘子法**

---

## 8. 条件极值与 Lagrange 乘子法

求泛函 $Q[y] = \int_a^b F(x,y,y')dx$ 在约束 $\int_a^b \varphi(x,y,y')dx = \alpha$ 下的极值。

构造新泛函：
$$Q^*[y] = \int_a^b (F + \lambda\varphi)dx - \lambda\alpha = \int_a^b F^*dx - \lambda\alpha$$
其中 $F^* = F + \lambda\varphi$，$\lambda$ 为 Lagrange 乘子（额外未知数）。

对新泛函使用 Euler 方程：
$$F^*_y - \frac{d}{dx}F^*_{y'} = 0$$

通解含 2 个积分常数 + $\lambda$，由 2 个边界条件 + 1 个等周条件确定。

---

## 9. 最速降线问题的完整求解

**物理模型**：质点在重力作用下从 A 到 B 沿某曲线滑下，求耗时最短的路径。

**Step 1**：建立泛函

能量守恒：$\frac12 mv^2 = mgy \Rightarrow v = \sqrt{2gy}$
弧长微分：$dS = \sqrt{1+y'^2}dx$
时间元：$dt = \frac{dS}{v} = \frac{\sqrt{1+y'^2}}{\sqrt{2gy}}dx$

总时间泛函（略去常数 $\frac{1}{\sqrt{2g}}$）：
$$Q[y] = \int_{x_A}^{x_B} \sqrt{\frac{1+y'^2}{y}}\,dx,\quad F(y,y') = \sqrt{\frac{1+y'^2}{y}}$$

**Step 2**：利用首次积分（$F$ 不显含 $x$）
$$F - y'F_{y'} = C$$

计算 $F_{y'} = \frac{y'}{\sqrt{y(1+y'^2)}}$，代入得：
$$\sqrt{\frac{1+y'^2}{y}} - \frac{y'^2}{\sqrt{y(1+y'^2)}} = C \Rightarrow \frac{1}{\sqrt{y(1+y'^2)}} = C$$

即：$y(1+y'^2) = D$（$D = 1/C^2$ 为常数）

**Step 3**：参数化解

令 $y' = \tan\theta$，代入得：
$$y = \frac{D}{1+\tan^2\theta} = D\cos^2\theta = \frac{D}{2}(1+\cos 2\theta)$$
$$dx = \frac{dy}{y'} = \frac{-D\sin 2\theta\,d\theta}{\tan\theta} = -D(1+\cos 2\theta)d\theta$$

积分得参数方程：
$$\begin{cases}
x = -\frac{D}{2}(2\theta + \sin 2\theta) + E \\
y = \frac{D}{2}(1 + \cos 2\theta)
\end{cases}$$

令 $2\theta = \pi - \phi$，得旋轮线（cycloid）的标准形式：
$$\boxed{\begin{cases}
x = r(\phi - \sin\phi) \\
y = r(1 - \cos\phi)
\end{cases}},\quad r = \frac{D}{2}$$

这就是最速降线的解——一条旋轮线！

---

## 10. 变分法记号与运算法则

变分法中常用运算法则可极大简化推导：

1. **变分算子与微分算子可交换**：
   $$\frac{d}{dx}(\delta y) = \delta(y')$$

2. **变分算子与积分算子可交换**：
   $$\delta \int_a^b F\,dx = \int_a^b \delta F\,dx$$

3. **复合函数的变分**（与微分完全类似）：
   $$\delta F(x,y,y') = \frac{\partial F}{\partial y}\delta y + \frac{\partial F}{\partial y'}\delta y'$$

利用这些法则，Euler 方程的推导可以写为：
$$\delta Q = \int (\frac{\partial F}{\partial y}\delta y + \frac{\partial F}{\partial y'}\delta y')dx = \int (\frac{\partial F}{\partial y} - \frac{d}{dx}\frac{\partial F}{\partial y'})\delta y\,dx = 0$$

---

## 11. 经典变分法与 FEM 的关系

**Ritz 法**：设 $u \approx \sum c_i\varphi_i(x)$，使 $\partial\Pi/\partial c_i = 0$
- 局限：$\varphi_i$ 必须满足全部边界条件
- 复杂几何下很难构造

**FEM**：将求解域划分为单元，每个单元独立建立形函数
- 形函数不必满足单元间的协调条件（可有间断）
- 程序标准化、自动化
- 可以看作"在单元层面上应用 Ritz/Galerkin 法"

```
经典变分法 (Euler-Lagrange) → 精确但求解困难
   ↓
Ritz 法 (整体 trial function) → 边界条件难满足
   ↓
Courant 分片近似 (1943) → 克服边界困难
   ↓
FEM (单元级 trial function, 1960) → 程序化、通用化
```

---

## 12. 解题常见陷阱

### 陷阱 1：错把全导当偏导
Euler 方程中 $\frac{\partial F}{\partial y'}$ 是偏导，但 $\frac{d}{dx}$ 是全导：
$$\frac{d}{dx}F_{y'} = \frac{\partial F_{y'}}{\partial x} + \frac{\partial F_{y'}}{\partial y}y' + \frac{\partial F_{y'}}{\partial y'}y''$$

### 陷阱 2：忘记边界项
分部积分时，边界项 $[F_{y'}\delta y]_{x_1}^{x_2}$ 必须处理。如果 $\delta y$ 在端点为 0，该项消失。如果端点自由，必须令 $F_{y'}=0$（自然边界条件）。

### 陷阱 3：线性泛函判断不全
只检查了齐次性而漏了可加性，或反之。**两者必须同时满足**。

### 陷阱 4：条件极值中 $\lambda$ 的处理
引入 $\lambda$ 后，$F^* = F + \lambda\varphi$，$\lambda$ 和 $y(x)$ 都是未知函数，对两者都要用 Euler 方程。

---

> **对应作业**：[HW2 Q1: 最短路径](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q2: 三阶导数 Euler 方程](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q3: Lagrange 乘子](../04-Homework-Solutions/2026w/HW2-Problem.md) · [Q4: 泛函极值函数](../04-Homework-Solutions/2026w/HW2-Problem.md)
> **往年相关**：[past/HW2/homework 2](../04-Homework-Solutions/past/HW2/homework%202.md) · [答案(LIU Sai)](../04-Homework-Solutions/past/HW2/Ans%20to%20HM2_LIU%20Sai_handed%20in.md)
