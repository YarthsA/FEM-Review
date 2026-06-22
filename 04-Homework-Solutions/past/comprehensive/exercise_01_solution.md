# 习题一 解答

---

## 1. 推导线性变换及其逆变换的矩阵表达与指标表达

**题目**：推导线性变换及其逆变换的矩阵表达与指标表达。

**解答**：

**矩阵表达**：

设线性变换 $y_i = a_{ij}x_j$，写成矩阵形式：

$$\mathbf{y} = \mathbf{A}\mathbf{x}$$

其中 $\mathbf{A} = [a_{ij}]$ 是变换矩阵。

**逆变换**：

若 $\mathbf{A}$ 可逆，则逆变换为：

$$\mathbf{x} = \mathbf{A}^{-1}\mathbf{y}$$

**指标表达**：

正变换：$y_i = a_{ij}x_j$（$j$ 为哑指标，对 $j$ 求和）

逆变换：$x_i = b_{ij}y_j$，其中 $b_{ij}$ 是 $\mathbf{A}^{-1}$ 的元素，满足 $a_{ik}b_{kj} = \delta_{ij}$

---

## 2. 证明 $e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}$

**题目**：利用 Lagrange 公式 $(a \times b) \bullet (c \times d) = (a \bullet c)(b \bullet d) - (a \bullet d)(b \bullet c)$，证明 $e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}$。

**解答**：

将 $e_{ksp}e_{ipj}$ 写成向量形式：

$$e_{ksp}e_{ipj} = (\mathbf{e}_k \times \mathbf{e}_s) \cdot (\mathbf{e}_i \times \mathbf{e}_p) \cdot \delta_{pj}$$

利用 Lagrange 恒等式展开：

$$e_{ksp}e_{ipj} = (\mathbf{e}_k \cdot \mathbf{e}_i)(\mathbf{e}_s \cdot \mathbf{e}_p) - (\mathbf{e}_k \cdot \mathbf{e}_p)(\mathbf{e}_s \cdot \mathbf{e}_i)$$

由于 $\mathbf{e}_i \cdot \mathbf{e}_j = \delta_{ij}$（标准正交基）：

$$e_{ksp}e_{ipj} = \delta_{ki}\delta_{sp} - \delta_{kp}\delta_{si}$$

对 $p$ 求和（$p=j$）：

$$\boxed{e_{ksp}e_{ipj} = \delta_{is}\delta_{jk} - \delta_{ik}\delta_{js}}$$

---

## 3. 证明速度、加速度为矢量

**题目**：证明速度、加速度为矢量。

**解答**：

**速度**：$\mathbf{v} = \frac{d\mathbf{r}}{dt}$，**加速度**：$\mathbf{a} = \frac{d^2\mathbf{r}}{dt^2}$

矢量的定义是满足坐标变换律的量。在坐标系 $K$ 中，$\mathbf{r} = x_i\mathbf{e}_i$；在坐标系 $K'$ 中，$\mathbf{r} = x'_i\mathbf{e}'_i$。

坐标变换关系：$x'_i = Q_{ij}x_j$，其中 $Q_{ij} = \mathbf{e}'_i \cdot \mathbf{e}_j$ 是正交变换矩阵。

**速度的变换**：

$$v'_i = \frac{dx'_i}{dt} = Q_{ij}\frac{dx_j}{dt} = Q_{ij}v_j$$

满足矢量变换律，故速度是矢量。

**加速度的变换**：

$$a'_i = \frac{dv'_i}{dt} = Q_{ij}\frac{dv_j}{dt} = Q_{ij}a_j$$

满足矢量变换律，故加速度是矢量。 $\blacksquare$

---

## 4. 排列符号 $e_{ijk}$ 是三阶张量

**题目**：证明排列符号 $e_{ijk}$ 是三阶张量。

**解答**：

需要证明 $e_{ijk}$ 在坐标变换下满足张量变换律。

设坐标变换 $x'_i = Q_{ij}x_j$，其中 $Q_{ij}$ 是正交变换矩阵。

三阶张量的变换律为：$T'_{ijk} = Q_{ip}Q_{jq}Q_{kr}T_{pqr}$

对 $e_{ijk}$：

$$Q_{ip}Q_{jq}Q_{kr}e_{pqr} = \det(\mathbf{Q}) \cdot e_{ijk}$$

由于 $\mathbf{Q}$ 是正交矩阵，$\det(\mathbf{Q}) = \pm 1$。

- 当 $\det(\mathbf{Q}) = +1$（旋转）：$e'_{ijk} = e_{ijk}$
- 当 $\det(\mathbf{Q}) = -1$（反射）：$e'_{ijk} = -e_{ijk}$

因此 $e_{ijk}$ 是**三阶张量**（更准确地说是**三阶伪张量**）。 $\blacksquare$

---

## 5. 证明应变为二阶张量

**题目**：证明应变为二阶张量。

**解答**：

应变张量定义为：$\varepsilon_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)$

在坐标变换 $x'_i = Q_{ij}x_j$ 下，位移矢量变换为 $u'_i = Q_{ij}u_j$。

$$\varepsilon'_{ij} = \frac{1}{2}\left(\frac{\partial u'_i}{\partial x'_j} + \frac{\partial u'_j}{\partial x'_i}\right)$$

利用链式法则：

$$\frac{\partial u'_i}{\partial x'_j} = \frac{\partial (Q_{ik}u_k)}{\partial (Q_{jl}x_l)} = Q_{ik}Q_{jl}\frac{\partial u_k}{\partial x_l}$$

因此：

$$\varepsilon'_{ij} = Q_{ik}Q_{jl}\frac{1}{2}\left(\frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k}\right) = Q_{ik}Q_{jl}\varepsilon_{kl}$$

满足二阶张量变换律，故应变是二阶张量。 $\blacksquare$

---

## 6. 张量运算

**题目**：已知 $T_{ij} = \begin{bmatrix} 1 & 1 & 0\\ 1 & 2 & 2\\ 0 & 2 & 3 \end{bmatrix}$，试求：(1) $T_{ii}$ (2) $T_{ij}T_{ij}$

**解答**：

**(1) 求 $T_{ii}$（迹）**

$$T_{ii} = T_{11} + T_{22} + T_{33} = 1 + 2 + 3 = \boxed{6}$$

**(2) 求 $T_{ij}T_{ij}$（Frobenius 范数的平方）**

$$T_{ij}T_{ij} = \sum_{i=1}^3\sum_{j=1}^3 T_{ij}^2$$

$$= 1^2 + 1^2 + 0^2 + 1^2 + 2^2 + 2^2 + 0^2 + 2^2 + 3^2 = \boxed{24}$$

---

## 7. 验证 $\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = (\mathbf{a} \cdot \mathbf{c})\mathbf{b} - (\mathbf{a} \cdot \mathbf{b})\mathbf{c}$

**题目**：试验证 $\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = (\mathbf{a} \cdot \mathbf{c})\mathbf{b} - (\mathbf{a} \cdot \mathbf{b})\mathbf{c}$。

**解答**：

**指标形式证明**：

左边：$[\mathbf{a} \times (\mathbf{b} \times \mathbf{c})]_i = e_{ijk}a_j(e_{klm}b_l c_m)$

利用 $e_{ijk}e_{klm} = e_{kij}e_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}$：

$$= (\delta_{il}\delta_{jm} - \delta_{im}\delta_{jl})a_j b_l c_m = a_j b_i c_j - a_j b_j c_i$$

$$= (a_j c_j)b_i - (a_j b_j)c_i = (\mathbf{a} \cdot \mathbf{c})b_i - (\mathbf{a} \cdot \mathbf{b})c_i$$

因此：$\mathbf{a} \times (\mathbf{b} \times \mathbf{c}) = (\mathbf{a} \cdot \mathbf{c})\mathbf{b} - (\mathbf{a} \cdot \mathbf{b})\mathbf{c}$ $\blacksquare$

---

## 8. 证明 $\nabla(\varphi\gamma) = \varphi\nabla\gamma + \gamma\nabla\varphi$

**题目**：证明 $\nabla(\varphi\gamma) = \varphi\nabla\gamma + \gamma\nabla\varphi$（$\varphi, \gamma$ 为标量函数）。

**解答**：

梯度算子在直角坐标系中定义为：$\nabla = \mathbf{e}_i\frac{\partial}{\partial x_i}$

对标量函数乘积 $\varphi\gamma$ 取梯度：

$$[\nabla(\varphi\gamma)]_i = \frac{\partial(\varphi\gamma)}{\partial x_i}$$

由偏导数乘积法则：

$$= \frac{\partial\varphi}{\partial x_i}\gamma + \varphi\frac{\partial\gamma}{\partial x_i}$$

写回向量形式：

$$\nabla(\varphi\gamma) = \gamma\nabla\varphi + \varphi\nabla\gamma = \varphi\nabla\gamma + \gamma\nabla\varphi$$

$\blacksquare$

---

## 9. 选作：球坐标系下的平衡方程

**题目**：如何采用坐标变换法推导球坐标系下的平衡方程？（提示：可以把球坐标系看作直角坐标系经过两次柱坐标系变换得到）

**解答**：

**步骤**：

1. **第一次变换**：直角坐标 $(x,y,z)$ → 柱坐标 $(r,\theta,z)$
   - $x = r\cos\theta$, $y = r\sin\theta$, $z = z$

2. **第二次变换**：柱坐标 $(r,\theta,z)$ → 球坐标 $(R,\phi,\theta)$
   - $r = R\sin\phi$, $z = R\cos\phi$, $\theta = \theta$

3. **坐标变换矩阵**：计算 Jacobi 矩阵 $\frac{\partial(x,y,z)}{\partial(R,\phi,\theta)}$

4. **平衡方程变换**：利用张量变换律将直角坐标系下的平衡方程 $\sigma_{ij,j} + f_i = 0$ 变换到球坐标系。

最终得到球坐标系下的平衡方程（略去详细推导）：
- $R$ 方向：$\frac{\partial\sigma_{RR}}{\partial R} + \frac{1}{R}\frac{\partial\sigma_{R\phi}}{\partial\phi} + \frac{1}{R\sin\phi}\frac{\partial\sigma_{R\theta}}{\partial\theta} + \cdots + f_R = 0$
- $\phi$ 方向和 $\theta$ 方向类似
