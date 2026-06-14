# Variation Method & FEA Course — Homework 1 参考答案

> 交大 FEM 课程 · 张量分析专题
>
> 本套题考察张量分析的基本功：指标运算、Kronecker δ 和置换符号的性质、张量的判定与证明。掌握这些是后续学习弹性力学张量形式的基础。

---

## 1. 利用 Lagrange 恒等式证明置换符号恒等式

> 核心思路：把 $e_{ijk}$ 的乘积转化为向量叉积，再用 Lagrange 恒等式展开成点积组合，最后用 δ 改写点积。

**题目**：利用 $(\boldsymbol{a} \times \boldsymbol{b}) \cdot (\boldsymbol{c} \times \boldsymbol{d}) = (\boldsymbol{a} \cdot \boldsymbol{c})(\boldsymbol{b} \cdot \boldsymbol{d}) - (\boldsymbol{a} \cdot \boldsymbol{d})(\boldsymbol{b} \cdot \boldsymbol{c})$ 证明：
$$e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}$$

**证明**：

**Step 1**：将 Lagrange 恒等式左右两侧写成指标形式。

左侧用置换符号展开：
$$
(\boldsymbol{a} \times \boldsymbol{b}) \cdot (\boldsymbol{c} \times \boldsymbol{d}) = (e_{ijk}a_j b_k)(e_{ipq}c_p d_q) = e_{ijk}e_{ipq}\,a_j b_k c_p d_q
$$

右侧展开：
$$
(\boldsymbol{a} \cdot \boldsymbol{c})(\boldsymbol{b} \cdot \boldsymbol{d}) - (\boldsymbol{a} \cdot \boldsymbol{d})(\boldsymbol{b} \cdot \boldsymbol{c}) = (a_j c_j)(b_k d_k) - (a_j d_j)(b_k c_k) = a_j b_k(c_j d_k - c_k d_j)
$$

**Step 2**：将右侧写成与左侧一致的形式。

引入 Kronecker $\delta$：
$$
a_j b_k(c_j d_k - c_k d_j) = a_j b_k(\delta_{jp}\delta_{kq} - \delta_{jq}\delta_{kp})c_p d_q
$$

**Step 3**：比较两侧。

由于 $\boldsymbol{a}, \boldsymbol{b}, \boldsymbol{c}, \boldsymbol{d}$ 为任意向量，比较得：
$$
e_{ijk}e_{ipq} = \delta_{jp}\delta_{kq} - \delta_{jq}\delta_{kp} \tag{★}
$$

**Step 4**：通过指标重命名将 (★) 转化为待证形式。

在 (★) 式中，分别将指标 $(i,j,k,p,q)$ 置换为 $(k,s,p,i,j)$：
$$
e_{ksp}e_{kij} = \delta_{si}\delta_{pj} - \delta_{sj}\delta_{pi}
$$

注意 $e_{kij} = e_{ijk}$（偶排列），且 $e_{ipj} = -e_{ijp}$ ... 

更简洁的做法：利用置换符号的循环对称性。

由 (★) 式 $e_{ijk}e_{ipq} = \delta_{jp}\delta_{kq} - \delta_{jq}\delta_{kp}$，对第二项做指标替换 $(i,j,k,p,q) \to (k,s,p,\Box,\blacksquare)$：

令 $i \to k, j \to s, k \to p$，则 $e_{ijk} \to e_{ksp}$；
对于 $e_{ipq}$，令 $i \to \color{blue}{k}$（第一个指标与 $e_{ksp}$ 的第一个指标...）

**更好方法**：使用标准的 $\varepsilon$-$\delta$ 恒等式 $e_{ijm}e_{ikn} = \delta_{jk}\delta_{mn} - \delta_{jn}\delta_{mk}$。

取 $e_{ksp}e_{ipj}$，先将 $e_{ipj}$ 循环置换：
$$
e_{ipj} = e_{jip} \quad (\text{偶排列，循环：} i \to p \to j \to i)
$$

则：
$$
e_{ksp}e_{ipj} = e_{ksp}e_{jip}
$$

使用标准恒等式 $e_{ijk}e_{mnk} = \delta_{im}\delta_{jn} - \delta_{in}\delta_{jm}$（在第三指标上缩并），令：
- 第一组：$i=k, j=s, k=p$
- 第二组：$m=j, n=i, k=p$

得：
$$
e_{ksp}e_{jip} = \delta_{kj}\delta_{si} - \delta_{ki}\delta_{sj}
$$

由于 Kronecker $\delta$ 对称：
$$
\delta_{kj}\delta_{si} = \delta_{is}\delta_{jk}, \quad \delta_{ki}\delta_{sj} = \delta_{ik}\delta_{js}
$$

因此：
$$
e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}
$$

**证毕**。■

---

## 2. 证明 $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$ 是二阶张量

**证明**：

需证明 $\varepsilon_{ij}$ 在坐标变换下服从二阶张量的变换律。

设坐标变换：$x_{i'} = l_{i'j}x_j$，其中 $l_{i'j} = \cos(x_{i'}, x_j)$ 为方向余弦矩阵 $\mathbf{L}$ 的元素。

**Step 1**：位移向量 $\mathbf{u}$ 是一阶张量（向量），变换规律：
$$
u_{i'} = l_{i'i}u_i
$$

**Step 2**：导数算子的变换（链式法则）：
$$
\frac{\partial}{\partial x_{j'}} = \frac{\partial x_m}{\partial x_{j'}}\frac{\partial}{\partial x_m} = l_{mj'}\frac{\partial}{\partial x_m}
$$

其中 $l_{mj'} = \cos(x_m, x_{j'})$ 是逆变换矩阵的元素。由于 $\mathbf{L}$ 正交，$\mathbf{L}^{-1} = \mathbf{L}^T$，故 $l_{mj'} = l_{j'm}$。

**Step 3**：在新坐标系中计算 $\varepsilon_{i'j'}$：

$$
\begin{aligned}
\varepsilon_{i'j'} &= \frac{1}{2}(u_{i',j'} + u_{j',i'}) \\
&= \frac{1}{2}\left(\frac{\partial u_{i'}}{\partial x_{j'}} + \frac{\partial u_{j'}}{\partial x_{i'}}\right) \\
&= \frac{1}{2}\left[\frac{\partial}{\partial x_{j'}}(l_{i'i}u_i) + \frac{\partial}{\partial x_{i'}}(l_{j'j}u_j)\right]
\end{aligned}
$$

**Step 4**：代入导数变换：

$$
\begin{aligned}
\varepsilon_{i'j'} &= \frac{1}{2}\left[ l_{i'i}l_{mj'}\frac{\partial u_i}{\partial x_m} + l_{j'j}l_{ni'}\frac{\partial u_j}{\partial x_n} \right] \\
&= \frac{1}{2}\left[ l_{i'i}l_{mj'}u_{i,m} + l_{j'j}l_{ni'}u_{j,n} \right]
\end{aligned}
$$

**Step 5**：在第二项中交换哑指标 $(j,n) \to (i,m)$：

$$
\varepsilon_{i'j'} = \frac{1}{2}\left[ l_{i'i}l_{mj'}u_{i,m} + l_{j'i}l_{mi'}u_{i,m} \right]
$$

利用 $l_{mj'} = l_{j'm}$（正交性）：
$$
l_{j'm} = l_{mj'}, \quad l_{mi'} = l_{i'm}
$$

所以：
$$
\varepsilon_{i'j'} = \frac{1}{2}\left[ l_{i'i}l_{j'm} + l_{j'i}l_{i'm} \right] u_{i,m}
$$

**Step 6**：而二阶张量变换的目标是：
$$
l_{i'i}l_{j'j}\varepsilon_{ij} = l_{i'i}l_{j'j} \cdot \frac{1}{2}(u_{i,j}+u_{j,i}) = \frac{1}{2}\left(l_{i'i}l_{j'j} + l_{i'j}l_{j'i}\right)u_{i,j}
$$

比较得两边均等于 $\frac{1}{2}\left(l_{i'i}l_{j'j} + l_{i'j}l_{j'i}\right)u_{i,j}$，即：
$$
\varepsilon_{i'j'} = l_{i'i}l_{j'j}\varepsilon_{ij}
$$

**Step 7**：这正是二阶张量（协变）的变换规律。因此 $\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$ 是二阶张量。

**证毕**。■

---

## 3. 证明置换符号 $e_{ijk}$ 是三阶张量

**证明**：

**Step 1**：张量变换律。若 $e_{ijk}$ 是三阶张量，则在坐标变换下应满足：
$$
e_{i'j'k'} = l_{i'i}l_{j'j}l_{k'k}e_{ijk}
$$

**Step 2**：考虑 $l_{i'i}l_{j'j}l_{k'k}e_{ijk}$。这是三阶交错张量的分量。根据行列式的定义：
$$
l_{i'i}l_{j'j}l_{k'k}e_{ijk} = \det(\mathbf{L})\, e_{i'j'k'}
$$

其中方向余弦矩阵 $\mathbf{L}$ 的行列式 $\det(\mathbf{L}) = \pm 1$。

**Step 3**：证明确实满足张量变换律。

- 若 $(i',j',k')$ 是 $(1,2,3)$ 的偶排列，则 $e_{i'j'k'} = 1$：
  $$
  l_{i'i}l_{j'j}l_{k'k}e_{ijk} = \det(\mathbf{L}) \cdot 1 = \det(\mathbf{L}) \cdot e_{i'j'k'}
  $$

- 若 $(i',j',k')$ 是 $(1,2,3)$ 的奇排列，则 $e_{i'j'k'} = -1$：
  $$
  l_{i'i}l_{j'j}l_{k'k}e_{ijk} = \det(\mathbf{L}) \cdot (-1) = \det(\mathbf{L}) \cdot e_{i'j'k'}
  $$

- 若 $(i',j',k')$ 有重复指标，$e_{i'j'k'} = 0$：
  $$
  l_{i'i}l_{j'j}l_{k'k}e_{ijk} = 0 = \det(\mathbf{L}) \cdot e_{i'j'k'}
  $$

**Step 4**：对**正常转动**（$\det(\mathbf{L}) = 1$）：
$$
e_{i'j'k'} = l_{i'i}l_{j'j}l_{k'k}e_{ijk}
$$

这满足三阶张量的变换规律。因此 $e_{ijk}$ 是三阶张量（或称**三阶交错张量**，在正常转动下是真正张量）。

> **Note**: 当包含反射变换（$\det(\mathbf{L}) = -1$）时，$e_{ijk}$ 是**赝张量 (pseudotensor)**，因为引入了符号 $\det(\mathbf{L})$。

**证毕**。■

---

## 4. 给定张量的指标运算

> 缩并就是对一对重复指标求和。$T_{ii}$ 是迹（对角线之和），$T_{ij}T_{ij}$ 是所有元素平方和，$T_{ij}T_{ji}$ 先转置再做内积。

**题目**：
$$T = \begin{bmatrix} 4 & -1 & 2 \\ -1 & 3 & 0 \\ 2 & 0 & 5 \end{bmatrix}$$

### (1) $T_{ii}$（迹 / trace）

$$
T_{ii} = T_{11} + T_{22} + T_{33} = 4 + 3 + 5 = \boxed{12}
$$

### (2) $T_{ij}T_{ij}$

$$
\begin{aligned}
T_{ij}T_{ij} &= \sum_{i=1}^3\sum_{j=1}^3 T_{ij}^2 \\
&= T_{11}^2 + T_{12}^2 + T_{13}^2 + T_{21}^2 + T_{22}^2 + T_{23}^2 + T_{31}^2 + T_{32}^2 + T_{33}^2 \\
&= 4^2 + (-1)^2 + 2^2 + (-1)^2 + 3^2 + 0^2 + 2^2 + 0^2 + 5^2 \\
&= 16 + 1 + 4 + 1 + 9 + 0 + 4 + 0 + 25 \\
&= \boxed{60}
\end{aligned}
$$

### (3) $T_{ij}T_{ji}$

由于 $T$ 是对称矩阵（$T_{ij}=T_{ji}$），$T_{ij}T_{ji} = T_{ij}T_{ij}$。

计算 $T^2$：
$$
T^2 = \begin{bmatrix} 4 & -1 & 2 \\ -1 & 3 & 0 \\ 2 & 0 & 5 \end{bmatrix}
\begin{bmatrix} 4 & -1 & 2 \\ -1 & 3 & 0 \\ 2 & 0 & 5 \end{bmatrix}
= \begin{bmatrix} 21 & -7 & 18 \\ -7 & 10 & -2 \\ 18 & -2 & 29 \end{bmatrix}
$$

$$
T_{ij}T_{ji} = \text{tr}(T^2) = 21 + 10 + 29 = \boxed{60}
$$

> **验证**：$T_{ij}T_{ji} = \sum_i\sum_j T_{ij}T_{ji}$
> $$
> \begin{aligned}
> &= T_{11}T_{11} + T_{12}T_{21} + T_{13}T_{31} \\
> &\quad + T_{21}T_{12} + T_{22}T_{22} + T_{23}T_{32} \\
> &\quad + T_{31}T_{13} + T_{32}T_{23} + T_{33}T_{33} \\
> &= 4\cdot4 + (-1)(-1) + 2\cdot2 \\
> &\quad + (-1)(-1) + 3\cdot3 + 0\cdot0 \\
> &\quad + 2\cdot2 + 0\cdot0 + 5\cdot5 \\
> &= 16+1+4+1+9+0+4+0+25 = 60
> \end{aligned}
> $$
> 结果与 (2) 相同，因为 $T$ 是对称张量。

---

### 关键公式备忘

| 运算 | 含义 | 本例结果 |
|------|------|----------|
| $T_{ii}$ | 迹 (trace) | 12 |
| $T_{ij}T_{ij}$ | Frobenius 范数平方 | 60 |
| $T_{ij}T_{ji}$ | $\text{tr}(T^2)$ | 60 |
