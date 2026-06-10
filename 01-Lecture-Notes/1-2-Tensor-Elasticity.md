# 第2章：张量分析与弹性力学

## 1. 指标记号与求和约定

### Einstein 求和约定
在同一项中，重复出现的下标表示对该下标求和（哑指标 dummy index）：
$$
a_i b_i = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^3 a_i b_i
$$

### Kronecker 符号
$$
\delta_{ij} = \begin{cases}
1 & (i = j) \\
0 & (i \neq j)
\end{cases}
$$

性质：
$$
\delta_{ii} = 3,\quad \delta_{ij}a_j = a_i,\quad \delta_{ij}\delta_{jk} = \delta_{ik}
$$

### 置换符号（Levi-Civita 符号）
$$
e_{ijk} = \begin{cases}
0 & \text{任意两指标相等} \\
+1 & (i,j,k)\text{为偶排列 }(123,231,312) \\
-1 & (i,j,k)\text{为奇排列 }(132,213,321)
\end{cases}
$$

性质：
$$
e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}
$$
$$
e_{ijk}e_{ijt} = 2\delta_{kt}
$$

向量叉积：
$$
(\mathbf{a} \times \mathbf{b})_i = e_{ijk}a_jb_k
$$

---

## 2. 向量与坐标变换

### 基向量变换
设旧坐标系 $\{O; x_1, x_2, x_3\}$，基向量 $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$
新坐标系 $\{O; x_1', x_2', x_3'\}$，基向量 $\{\mathbf{e}_1', \mathbf{e}_2', \mathbf{e}_3'\}$

方向余弦矩阵（正交矩阵）：
$$
l_{i'j} = \cos(x_i', x_j)
$$

基向量变换：
$$
\mathbf{e}_{i'}' = l_{i'j}\mathbf{e}_j \quad\text{或}\quad \mathbf{e}' = \mathbf{L}\mathbf{e}
$$

逆变换（L 为正交矩阵，$\mathbf{L}^{-1} = \mathbf{L}^T$）：
$$
\mathbf{e} = \mathbf{L}^T\mathbf{e}'
$$

### 向量分量变换
$$
a_{i'} = l_{i'j}a_j \quad\text{(旧→新)}
$$
$$
a_i = l_{ij'}a_{j'} \quad\text{(新→旧)}
$$

---

## 3. 笛卡尔张量

### 标量（零阶张量）
- 在坐标变换下保持不变
  $$
  \phi' = \phi
  $$

### 向量（一阶张量）
- 按向量分量变换规则变换
  $$
  a_{i'}' = l_{i'j}a_j
  $$

### 张量（二阶张量）
- 按下列规则变换：
  $$
  a_{i'j'}' = l_{i'i}l_{j'j}a_{ij}
  $$

### 张量的基本运算

1. **加法**：同阶张量对应分量相加
2. **张量积（外积）**：$C_{ijkl} = A_{ij}B_{kl}$，升阶
3. **缩并（Contraction）**：令两个指标相同，降二阶
   - $A_{ij} \to A_{ii} = \text{tr}(\mathbf{A})$ (标量不变量)
4. **张量商定理（Quotient law）**

---

## 4. 矩阵变换

线性变换 $\Psi$ 在基 $\{\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3\}$ 下的矩阵表示为 $\mathbf{A}$：
$$
\Psi\mathbf{a} = \mathbf{e}^T\mathbf{A}\mathbf{a}
$$

坐标变换下矩阵的变换规则：
$$
\mathbf{A}' = \mathbf{L}\mathbf{A}\mathbf{L}^T
$$

### 张量的分解
任意二阶张量可分解为：
- **对称部分**：$A_{(ij)} = \frac{1}{2}(A_{ij} + A_{ji})$
- **反对称部分**：$A_{[ij]} = \frac{1}{2}(A_{ij} - A_{ji})$

### 张量的不变量
- 第一不变量（迹）：$I_1 = A_{ii}$
- 第二不变量：$I_2 = \frac{1}{2}(A_{ii}A_{jj} - A_{ij}A_{ji})$
- 第三不变量：$I_3 = \det(\mathbf{A})$

---

## 5. 弹性力学中的常用张量

### 应力张量 $\sigma_{ij}$
- $\sigma_{11}, \sigma_{22}, \sigma_{33}$ — 正应力
- $\sigma_{12}, \sigma_{23}, \sigma_{31}$ — 剪应力
- 对称张量：$\sigma_{ij} = \sigma_{ji}$
- 平衡方程张量形式：$\sigma_{ij,j} + f_i = 0$

### 应变张量 $\varepsilon_{ij}$
- 几何方程张量形式：$\varepsilon_{ij} = \frac{1}{2}(u_{i,j} + u_{j,i})$
- 对称张量：$\varepsilon_{ij} = \varepsilon_{ji}$

### 本构方程张量形式
$$
\sigma_{ij} = C_{ijkl}\varepsilon_{kl}
$$

对各向同性线弹性材料：
$$
\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}
$$

其中 $\lambda$ 和 $G$ 为 Lame 常数，$\varepsilon_{kk} = \varepsilon_{11}+\varepsilon_{22}+\varepsilon_{33} = \text{体积应变}$

---

## 6. 解题重要公式速查

| 概念 | 向量/矩阵形式 | 张量/指标形式 |
|------|-------------|-------------|
| 位移 | $\mathbf{u} = (u,v,w)^T$ | $u_i$ |
| 应变-位移 | $\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$ | $\varepsilon_{ij} = \frac{1}{2}(u_{i,j}+u_{j,i})$ |
| 本构关系 | $\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$ | $\sigma_{ij} = \lambda\varepsilon_{kk}\delta_{ij} + 2G\varepsilon_{ij}$ |
| 平衡方程 | $[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$ | $\sigma_{ij,j} + f_i = 0$ |
| 体积应变 | — | $\varepsilon_v = \varepsilon_{kk} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z}$ |

---

## 7. 张量证明题通用模板

考试中常出现"证明某量是张量"的题目，通用步骤：

**Step 1**：写出坐标变换关系
$$x_{i'} = l_{i'j}x_j,\quad l_{i'j} = \cos(x_{i'}, x_j)$$

**Step 2**：写出该量的定义（在旧坐标系中）

**Step 3**：在新坐标系中表达该量，利用链式法则变换导数
$$\frac{\partial}{\partial x_{j'}} = l_{km}\frac{\partial}{\partial x_m}$$

**Step 4**：整理成 $l_{i'i}l_{j'j}\cdots$ 乘以原分量的形式

**Step 5**：与张量变换律对比 → 得证

### 例：证明 $\varepsilon_{ij} = \frac12(u_{i,j}+u_{j,i})$ 是二阶张量

**证明**：
$$\begin{aligned}
\varepsilon_{i'j'} &= \frac12(u_{i',j'} + u_{j',i'}) \\
&= \frac12[l_{i'i}l_{mj'}u_{i,m} + l_{j'j}l_{ni'}u_{j,n}] \\
&= \frac12[l_{i'i}l_{j'm} + l_{j'i}l_{i'm}]u_{i,m} \\
&= l_{i'i}l_{j'j}\cdot\frac12(u_{i,j}+u_{j,i}) = l_{i'i}l_{j'j}\varepsilon_{ij}
\end{aligned}$$

符合二阶张量变换律，得证。■

### 例：证明 $e_{ijk}$ 是三阶张量（对正常转动）

$$\det(\mathbf{L}) = 1\text{ 时：}\quad e_{i'j'k'} = l_{i'i}l_{j'j}l_{k'k}e_{ijk}$$

即满足三阶张量变换律。对反射变换 $\det(\mathbf{L}) = -1$，$e_{ijk}$ 是**赝张量**。

---

## 8. 常见恒等式速查

| 恒等式 | 用途 |
|--------|------|
| $\delta_{ii} = 3$ | 计算迹 |
| $\delta_{ij}a_j = a_i$ | 指标替换 |
| $\delta_{ij}\delta_{jk} = \delta_{ik}$ | 缩并消元 |
| $e_{ijk}e_{ist} = \delta_{js}\delta_{kt} - \delta_{jt}\delta_{ks}$ | **最常用！** |
| $e_{ijk}e_{ijt} = 2\delta_{kt}$ | 双重缩并 |
| $(\boldsymbol{a}\times\boldsymbol{b})\cdot(\boldsymbol{c}\times\boldsymbol{d}) = (\boldsymbol{a}\cdot\boldsymbol{c})(\boldsymbol{b}\cdot\boldsymbol{d}) - (\boldsymbol{a}\cdot\boldsymbol{d})(\boldsymbol{b}\cdot\boldsymbol{c})$ | Lagrange 恒等式 |
