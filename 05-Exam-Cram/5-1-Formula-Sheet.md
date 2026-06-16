# 考前冲刺速查 — 一页纸公式 + 证明模板

> 考前 30 分钟快速翻阅，配合 [概念速查](../02-Concepts-Formulas/2-1-Concepts-Glossary.md) 使用

---

## 1. 张量与指标记法

| 项目 | 公式 |
|------|------|
| Einstein 求和 | $a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3$ |
| Kronecker $\delta$ | $\delta_{ij} = \begin{cases}1 & i=j\\0 & i\neq j\end{cases}$，$\delta_{ii}=3$，$\delta_{ij}a_j = a_i$ |
| 置换符号 | $e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$，$e_{ijk}e_{ijt} = 2\delta_{kt}$ |
| 向量叉积 | $(\boldsymbol{a}\times\boldsymbol{b})_i = e_{ijk}a_j b_k$ |
| 张量变换 | $a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}$（$l_{i'j'}$ 为方向余弦矩阵元素），$e_{i'j'k'}' = \det(\boldsymbol{L}) e_{ijk}$（三阶赝张量） |

**恒等式证明模板**：$\varepsilon$-$\delta$ 恒等式 ← 利用 Lagrange 公式 $(\boldsymbol{a}\times\boldsymbol{b})\cdot(\boldsymbol{c}\times\boldsymbol{d}) = (\boldsymbol{a}\cdot\boldsymbol{c})(\boldsymbol{b}\cdot\boldsymbol{d}) - (\boldsymbol{a}\cdot\boldsymbol{d})(\boldsymbol{b}\cdot\boldsymbol{c})$

---

## 2. 变分法

| 项目 | 公式 |
|------|------|
| 泛函定义 | $Q[f] : D \to \mathbb{R}$，$D$ 为函数集 |
| 线性泛函 | $Q[c_1 y_1 + c_2 y_2] = c_1 Q[y_1] + c_2 Q[y_2]$ |
| 变分（Lagrange定义） | $\delta Q = \frac{\partial}{\partial\alpha} Q[y+\alpha\delta y]\big|_{\alpha=0}$ |
| **Euler 方程** | $\boxed{\frac{\partial F}{\partial y} - \frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right) = 0}$ |
| 含 $n$ 阶导数的 Euler | $\sum_{k=0}^n (-1)^k\frac{d^k}{dx^k}\left(\frac{\partial F}{\partial y^{(k)}}\right) = 0$ |
| 多元函数 Euler | $F_{y_i} - \frac{d}{dx}F_{y_i'} = 0\quad(i=1,\ldots,n)$ |
| 多变量 Euler | $F_z - \frac{\partial}{\partial x}F_p - \frac{\partial}{\partial y}F_q = 0$（$p=z_x, q=z_y$） |
| $F$ 不显含 $x$ 时 | $F - y'F_{y'} = C$（首次积分） |

---

## 3. FEM — 三角形单元刚度矩阵

### 3.1 形函数（CST）

$$N_i = \frac{1}{2\Delta}(a_i + b_i x + c_i y),\quad 
b_i = y_j - y_m,\quad c_i = x_m - x_j$$

$$\Delta = \frac{1}{2}\begin{vmatrix}1&x_i&y_i\\1&x_j&y_j\\1&x_m&y_m\end{vmatrix}$$

### 3.2 $[B]$ 矩阵

$$[B] = \frac{1}{2\Delta}\begin{pmatrix}
b_i & 0 & b_j & 0 & b_m & 0 \\
0 & c_i & 0 & c_j & 0 & c_m \\
c_i & b_i & c_j & b_j & c_m & b_m
\end{pmatrix}$$

### 3.3 弹性矩阵 $[D]$

**平面应力**（$E$ 为杨氏模量，$\nu$ 为泊松比）：
$$[D] = \frac{E}{1-\nu^2}\begin{pmatrix}1&\nu&0\\\nu&1&0\\0&0&\frac{1-\nu}{2}\end{pmatrix}$$

**平面应变**（$E$ 为杨氏模量，$\nu$ 为泊松比）：
$$[D] = \frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{pmatrix}1&\frac{\nu}{1-\nu}&0\\\frac{\nu}{1-\nu}&1&0\\0&0&\frac{1-2\nu}{2(1-\nu)}\end{pmatrix}$$

### 3.4 单元刚度矩阵

$$\boxed{[k]_e = t\Delta_e [B]^T[D][B]}\quad (\text{$t$ 为单元厚度，$\Delta_e$ 为单元面积})$$

### 3.5 总体集成

$$[K] = \sum_{n=1}^{NE} [k]_{e_n},\quad \{F\} = \sum_{n=1}^{NE} \{F\}_{e_n},\quad [K]\{\delta\} = \{F\}\quad (\text{$NE$ 为单元总数})$$

---

## 4. Euler-Bernoulli 梁单元（Hermite 形函数）

$$w(\xi) = N_1 w_1 + N_2 \theta_1 + N_3 w_2 + N_4 \theta_2,\quad \xi = \frac{x-x_1}{L}\quad (\text{$w_i$ 为挠度，$\theta_i$ 为转角})$$

$$N_1 = 1-3\xi^2+2\xi^3,\quad N_2 = L(\xi-2\xi^2+\xi^3)$$
$$N_3 = 3\xi^2-2\xi^3,\quad N_4 = L(-\xi^2+\xi^3)$$

$$[k]_e = \frac{EI}{L^3}\begin{pmatrix}
12 & 6L & -12 & 6L\\
6L & 4L^2 & -6L & 2L^2\\
-12 & -6L & 12 & -6L\\
6L & 2L^2 & -6L & 4L^2
\end{pmatrix}$$

---

## 5. Gauss 积分

$$\int_{-1}^1 f(\xi)d\xi \approx \sum_{i=1}^n w_i f(\xi_i)$$

| $n$ | $\xi_i$ | $w_i$ | 精确度 |
|-----|---------|-------|--------|
| 1 | 0 | 2 | 线性 |
| 2 | $\pm 1/\sqrt{3}$ | 1 | 三次 |
| 3 | $\pm\sqrt{0.6},\;0$ | $5/9,\;8/9$ | 五次 |

**2D**: $\int_{-1}^1\int_{-1}^1 f(\xi,\eta)d\xi d\eta \approx \sum_{i=1}^n\sum_{j=1}^n w_i w_j f(\xi_i,\eta_j)$

**减缩积分**：用低阶 Gauss 代替精确积分 → 避免剪切自锁

---

## 6. 弹性力学三大方程

| 方程 | 矩阵形式 | 张量形式 |
|------|----------|----------|
| 几何 | $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$（$[\partial]$ 为微分算子矩阵） | $\varepsilon_{ij} = \frac{1}{2}(u_{i,j}+u_{j,i})$ |
| 物理 | $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$ | $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$ |
| 平衡 | $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$（$\mathbf{f}$ 为体力向量） | $\sigma_{ij,j} + f_i = 0$ |

**Lame 常数**：$\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)},\; G = \frac{E}{2(1+\nu)}$

---

## 7. 常见试题证明模板

### 7.1 证明某量是张量
1. 写出在旧坐标系中的定义
2. 写出坐标变换 $x_{i'} = l_{i'j}x_j$
3. 利用链式法则变换导数 $\frac{\partial}{\partial x_{j'}} = l_{mj'}\frac{\partial}{\partial x_m}$
4. 代入新坐标，整理成 $l_{i'i}l_{j'j}\cdots\times\text{原分量}$ 形式
5. 与张量变换律对比 → 得证

### 7.2 推导 Euler 方程
1. 写出泛函 $Q[y] = \int F dx$
2. 一阶变分 $\delta Q = \int(F_y\delta y + F_{y'}\delta y')dx = 0$
3. 分部积分 $\int F_{y'}\delta y' dx = [F_{y'}\delta y] - \int\frac{d}{dx}F_{y'}\delta y\,dx$
4. 利用端点条件 $\delta y(a)=\delta y(b)=0$ 消去边界项
5. 得 $\int(F_y - \frac{d}{dx}F_{y'})\delta y\,dx = 0$
6. 由 $\delta y$ 任意性 → $F_y - \frac{d}{dx}F_{y'} = 0$

### 7.3 证明总刚度矩阵对称
$$\{\delta\}^T[K]\{\delta\} = \sum_e \{\delta\}_e^T[k]_e\{\delta\}_e \quad\text{且 $[k]_e$ 对称}$$

### 7.4 证明位移元下限性
位移元假设位移场 → 结构比实际刚硬 → 总势能 $\Pi \geq \Pi_{\text{true}}$ → 位移解 $\leq$ 真实解

---

## 8. 判断题 / 选择题常考结论

- ❌ 有限元解总是比精确解小 → ✓ 位移元满足下限性，但非协调单元可能不成立
- ✅ Galerkin 法的试探函数必须满足全部边界条件（位移+力）
- ✅ Ritz 法的试探函数只需满足位移（本质）边界条件
- ✅ 线性三角形单元（CST）内应力和应变为**常数**
- ❌ $Q[y] = \int_a^b y^2 dx$ 是线性泛函 → ✓ 不是！（不满足可加性）
- ✅ 总刚度矩阵 $[K]$ 具有对称性、稀疏性、非负定性
- ✅ $K_{ij}$ 的物理意义：第 $j$ 个 DOF 为单位位移时在第 $i$ 个 DOF 所需施加的力
- ✅ 等参元满足收敛条件：坐标与位移用相同形函数 + 形函数满足 $\sum N_i = 1$

---

> 祝考试顺利！🎯
