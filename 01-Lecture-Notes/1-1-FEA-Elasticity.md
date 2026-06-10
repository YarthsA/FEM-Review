# 第1章：引言与弹性力学基础

## 1. FEA/FEM 基本概念

**FEM (Finite Element Method)** / **FEA (Finite Element Analysis)** — 一种数值求解场问题（field problems）的方法。

### 核心思想
1. **离散化** — 将连续体（continuum）划分为有限个单元（elements）
2. **单元分析** — 每个单元采用统一的位移函数（形函数 interpolation）
3. **总体集成** — 将所有单元组装回整体结构进行分析

### 连续体 → 有限元的关键步骤
1. 将求解域划分为子域（subdomains）
2. 每个子域建立位移函数（类似于 classical variation method 中的 trial function）
3. 单元间不强制满足变形协调条件 → 大大简化了位移函数的建立
4. 当位移函数满足一定条件时，FEM 结果收敛到真实解

### FEM 与其他方法的对比

| 方法 | 特点 | 局限 |
|------|------|------|
| 解析法 | 精确解 | 仅适用于简单几何和边界条件 |
| 经典变分法（Ritz） | 近似解，trial function 需满足边界条件 | 复杂边界下 trial function 选择困难 |
| **FEM** | 逐单元建立函数，程序自动化 | 需足够网格密度 |

### FEM 的基本流程
```
连续体 → 离散化（网格划分）
       → 选择位移函数（形函数）
       → 建立单元刚度矩阵 [k]ₑ
       → 组装总体刚度矩阵 [K]
       → 引入边界条件
       → 求解线性方程组 [K]{δ} = {F}
       → 计算应变、应力
```

---

## 2. 弹性力学（Theory of Elasticity）基础

### 定义
- **弹性 (Elasticity)** — 外力去除后，弹性材料制成的物体恢复到原始状态
- **弹性体 (Elastic body)** — 以弹性为唯一材料性质的物理对象
- **弹性力学** — 研究弹性体在外载荷下应力和变形分布的学科

### 弹性力学 vs 材料力学

| 方向 | 材料力学 | 弹性力学 |
|------|----------|----------|
| 研究对象 | 杆、梁、柱等结构构件 | 板、壳、实体等一般结构 |
| 假设 | 需要较多简化假设（如平截面假设） | 假设较少，更接近真实 |
| 精度 | 工程近似 | 更精确 |
| 适用性 | 规则构件 | 任意几何形状和载荷条件 |

### 弹性力学的基本假设
1. 连续性 (Continuity)
2. 均匀性 (Homogeneity)
3. 各向同性 (Isotropy)
4. 完全弹性 (Perfect elasticity — Hooke's law)
5. 小变形 (Small deformation)
6. 无初始应力 (No initial stress)

---

## 3. 弹性力学的基本变量（三维）

### 位移向量
$$
\mathbf{u} = \begin{pmatrix} u \\ v \\ w \end{pmatrix}
$$

### 应变向量（6 个分量）
$$
\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \varepsilon_z \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{pmatrix}
$$

### 应力向量（6 个分量）
$$
\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x \\ \sigma_y \\ \sigma_z \\ \tau_{xy} \\ \tau_{yz} \\ \tau_{zx} \end{pmatrix}
$$

---

## 4. 三类基本方程

### (1) 几何方程（应变-位移关系）
$$
\boldsymbol{\varepsilon} = [\partial]\mathbf{u}
$$

其中微分算子矩阵：
$$
[\partial] = \begin{pmatrix} 
\frac{\partial}{\partial x} & 0 & 0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial x} & \frac{\partial}{\partial z} & 0 \\
0 & 0 & \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
\end{pmatrix}^T
$$

展开形式：
$$
\varepsilon_x = \frac{\partial u}{\partial x}, \quad
\varepsilon_y = \frac{\partial v}{\partial y}, \quad
\varepsilon_z = \frac{\partial w}{\partial z}
$$

$$
\gamma_{xy} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}, \quad
\gamma_{yz} = \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y}, \quad
\gamma_{zx} = \frac{\partial w}{\partial x} + \frac{\partial u}{\partial z}
$$

### (2) 物理方程（本构关系-各向同性线弹性）
$$
\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}
$$

弹性矩阵 D（Lame 常数形式）：
$$
\mathbf{D} = \begin{pmatrix}
\lambda + 2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda + 2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda + 2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}
$$

其中 Lame 常数：
$$
\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}, \quad
G = \frac{E}{2(1+\nu)}
$$

### (3) 平衡方程（应力-外力关系）
$$
[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}
$$

展开：
$$
\frac{\partial\sigma_x}{\partial x} + \frac{\partial\tau_{xy}}{\partial y} + \frac{\partial\tau_{zx}}{\partial z} + f_x = 0
$$
$$
\frac{\partial\tau_{xy}}{\partial x} + \frac{\partial\sigma_y}{\partial y} + \frac{\partial\tau_{yz}}{\partial z} + f_y = 0
$$
$$
\frac{\partial\tau_{zx}}{\partial x} + \frac{\partial\tau_{yz}}{\partial y} + \frac{\partial\sigma_z}{\partial z} + f_z = 0
$$

### (4) 边界条件

位移边界：
$$
\mathbf{u}|_{S_u} = \bar{\mathbf{u}}
$$

力边界：
$$
[\mathbf{n}]\boldsymbol{\sigma}|_{S_\sigma} = \mathbf{T}
$$

---

## 5. 最小势能原理（矩阵形式）

总势能泛函：
$$
\Pi = \int_\Omega \frac{1}{2}\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\mathbf{T}\,dS
$$

其中：
- 第一项：应变能
- 第二项：体力势能
- 第三项：面力势能

> **重要**: 虚功方程、虚位移原理和最小势能原理均可用于推导有限元的离散形式

---

## 6. 弹性力学发展简史

- **1678** — Hooke 发现弹性体的位移与外力成正比
- **1821/1823** — Navier / Cauchy 推导出线弹性边值问题的平衡方程
- **1855** — Saint-Venant 扭转和弯曲理论
- **1933** — Muskhelishvili 复变函数方法
- **1908** — Ritz 提出近似方法
- **1943** — Courant 扩展 Ritz 法，分片假设函数
- **1956** — Turner, Clough, Martin, Topp 引入矩阵位移法解平面应力
- **1960** — Clough 首次提出 "Finite Element Method"

---

## 7. FEA 常用术语

| 术语 | 含义 |
|------|------|
| **节点 (Node)** | 离散模型中连接单元的点，存储位移自由度 |
| **单元 (Element)** | 离散模型的基本构建块，有特定几何形状和插值函数 |
| **自由度 (DOF)** | 节点上独立未知量的个数（位移元每节点 $n$ 个，$n$ 为问题维数） |
| **网格 (Mesh)** | 所有节点和单元的集合 |
| **插值函数** | 用节点值表达单元内部任意点值的函数 |
| **形函数** | 插值函数中对应各节点的那部分基函数 |

### 单元类型

```
1D 单元：杆单元 (truss)、梁单元 (beam) → 线/曲线
2D 单元：三角形 (tri)、四边形 (quad) → 平面/曲面
3D 单元：四面体 (tet)、六面体 (hex) → 体
```

- **低阶单元**：线性插值（如 3 节点三角形 CST）
- **高阶单元**：二次及以上插值（如 6 节点三角形 LST）

### FEM 求解三阶段

```
① 前处理 (Preprocessing)
   └── 几何建模 → 网格划分 → 材料属性 → 边界条件 → 载荷
② 求解 (Solution)
   └── 形成单元矩阵 → 总装 → 施加BC → 求解线性方程组
③ 后处理 (Postprocessing)
   └── 位移 → 应变 → 应力 → 可视化
```

### 解的收敛方法
- **h 方法**：不改变单元阶数，细化网格 → 更密的网格
- **p 方法**：不改变网格，提高单元阶数 → 更高次多项式
- **hp 方法**：两者结合

### FEM 软件
ANSYS、ABAQUS、NASTRAN、MARC、ADINA、COMSOL、LS-DYNA
