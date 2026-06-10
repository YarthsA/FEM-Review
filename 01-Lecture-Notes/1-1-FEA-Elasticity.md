# FEA 概述与弹性力学基础

> **对应课件**：[`Chapter 1 Introduction to FEA.pdf`](../06-References/pdfs-originals/Chapter%201%20Introduction%20to%20FEA.pdf) 绪论部分 · [原文MD](../06-References/../06-References/../md_output/Chapter%201%20Introduction%20to%20FEA.md)
> **相关作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
> **前置知识**：高等数学、线性代数、大学物理

---

## 1.1 什么是有限元法？

> The finite element method (FEM), also called as finite element analysis (FEA), is a method for numerical solution of field problems.
> — 课程原始定义

### 1.1.1 FEM 的基本原理

FEM 将连续体（continuum）离散为有限个单元（elements），通过在每个单元上假设位移函数（形函数），利用变分原理建立代数方程组来求解全场问题。将连续体离散化后，整个问题的分析就转化为对每个单元的分析和所有单元的重新组合。

### 1.1.2 FEM 的三种理解途径

1. **结构矩阵法（Structural Matrix Method）**：从结构力学的矩阵位移法发展而来
2. **变分法（Variational Method）**：基于最小势能原理或虚位移原理
3. **加权残量法（Method of Weighted Residuals）**：直接处理微分方程

FEM 的通用性极强——整个计算过程由计算机自动完成，只需根据不同的工程问题改变输入即可。这种方法彻底改变了分析解的限制。

### 1.1.3 FEM 的发展简史

- **1908**：Ritz 提出用带未知量的试探函数近似能量泛函，得到求解未知量的方程组——但试探函数必须满足边界条件，复杂几何下极为困难
- **1943**：Courant 在解扭转问题时将截面划分为三角形区域，在每个三角形内假设翘曲函数均匀分布——克服了 Ritz 法要求整体函数满足边界条件的困难，预示了 FEM 思想
- **1956**：Turner, Clough, Martin, Topp 在航空学会年会上提出从矩阵位移法发展而来的一种新计算方法——将结构划分为三角形和矩形单元求解平面应力
- **1955**：Argyris 发表能量理论和结构分析的多篇论文，统一了弹性结构的基本能量原理
- **1960**：**Clough** 在论文"The finite element method in plane stress analysis"中首次引入"FEM"术语——学科诞生标志

### 1.1.4 商业 FEM 软件

NASTRAN, MARC, SAP—NONSAP, ADINA, ANSYS, ABAQUS, COMSOL, LS-DYNA, HyperMesh 等。

---

## 1.2 弹性力学（Theory of Elasticity）

### 1.2.1 基本定义

- **弹性（Elasticity）**：外力去除后，物体恢复到原始状态的性质
- **弹性体（Elastic body）**：以弹性为唯一材料特性的物理对象
- **弹性力学**：研究弹性体在外载荷下应力和变形分布规律的学科

### 1.2.2 弹性力学 vs 材料力学

弹性力学是固体力学的重要分支，而材料力学可能是固体力学中历史最悠久的分支。材料力学关注结构构件的强度、刚度和稳定性问题，但材料力学只能处理杆、梁、柱等结构构件，板、壳和实体结构则很难处理。即使对杆梁构件，仍有一些问题留待解决。

弹性力学在原则上对弹性体的几何形状和载荷形式没有限制，允许比材料力学更真实的假设。例如梁弯曲理论中的"平截面假设"在某些条件下并不适用。

### 1.2.3 弹性力学的基本假设

1. **连续性（Continuity）**
2. **均匀性（Homogeneity）**
3. **各向同性（Isotropy）**
4. **完全弹性（Perfect elasticity — Hooke's law）**
5. **小变形（Small deformation）**
6. **无初始应力（No initial stress）**

### 1.2.4 弹性力学的三类基本变量

三维弹性力学的基本变量：

**① 位移（3个分量）**：
$$\mathbf{u} = \begin{pmatrix} u \\ v \\ w \end{pmatrix}$$

**② 应变（6个分量）**：
$$\boldsymbol{\varepsilon} = \begin{pmatrix} \varepsilon_x \\ \varepsilon_y \\ \varepsilon_z \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{pmatrix}$$

**③ 应力（6个分量）**：
$$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_x \\ \sigma_y \\ \sigma_z \\ \tau_{xy} \\ \tau_{yz} \\ \tau_{zx} \end{pmatrix}$$

---

## 1.3 弹性力学的三类基本方程

### 1.3.1 几何方程
$$\boldsymbol{\varepsilon} = [\partial]\mathbf{u}$$

其中微分算子矩阵：
$$[\partial] = \begin{pmatrix} 
\frac{\partial}{\partial x} & 0 & 0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial z} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{\partial}{\partial x} & \frac{\partial}{\partial z} & 0 \\
0 & 0 & \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
\end{pmatrix}^T$$

### 1.3.2 物理方程
$$\boldsymbol{\sigma} = \mathbf{D}\boldsymbol{\varepsilon}$$

弹性矩阵 $\mathbf{D}$（Lame 常数形式）：
$$\mathbf{D} = \begin{pmatrix}
\lambda+2G & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2G & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2G & 0 & 0 & 0 \\
0 & 0 & 0 & G & 0 & 0 \\
0 & 0 & 0 & 0 & G & 0 \\
0 & 0 & 0 & 0 & 0 & G
\end{pmatrix}$$

Lame 常数：$\lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}$，$G = \frac{E}{2(1+\nu)}$

### 1.3.3 平衡方程
$$[\partial]^T\boldsymbol{\sigma} + \mathbf{f} = \mathbf{0}$$

### 1.3.4 边界条件
- 位移边界：$\mathbf{u}|_{S_u} = \bar{\mathbf{u}}$
- 外力边界：$[\mathbf{n}]\boldsymbol{\sigma}|_{S_\sigma} = \mathbf{T}$

### 1.3.5 最小势能原理

总势能泛函：
$$\Pi = \int_\Omega \frac12\boldsymbol{\varepsilon}^T\mathbf{D}\boldsymbol{\varepsilon}\,dV - \int_\Omega \mathbf{u}^T\mathbf{f}\,dV - \int_{S_\sigma} \mathbf{u}^T\mathbf{T}\,dS$$

- 第一项：应变能
- 第二项：体力势能
- 第三项：面力势能

> 在一切可能位移场中，真实位移场使总势能取最小值，即 $\delta\Pi = 0$。

---

## 1.4 工程设计的任务

工程设计的主要任务是运用固体力学理论（包括结构力学、弹性力学和塑性力学）对结构进行**强度、刚度和稳定性**分析。

---

## 1.5 FEM 的求解步骤

```
结构离散 → 按几何特点和精度要求划分单元和网格
  ↓
形成单元刚度矩阵和等效节点荷载列阵
  ↓
集成结构的总体刚度矩阵和荷载列阵
  ↓
引入强制边界条件
  ↓
求解方程 → 得到节点位移
  ↓
计算单元应变和应力
  ↓
必要的后处理
```

---

> **对应作业**：[HW1 Q2（应变张量证明）](../04-Homework-Solutions/2026w/HW1-Problem.md) · [HW1 Q4（指标运算）](../04-Homework-Solutions/2026w/HW1-Problem.md)
