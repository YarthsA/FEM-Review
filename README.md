# 变分法与有限元分析 — 课程复习指南

> **上海交通大学 · 有限元方法课程（土木工程方向）**
> 基于课堂笔记、作业及往年参考资料系统整理的复习资料库。

---

## 仓库结构

```
FEM-Review/
├── 00-Cross-Reference/       ← 交叉引用索引（笔记→PDF→作业双向链接）
│   └── index.md
│
├── 01-Lecture-Notes/          ← 六章讲义笔记（教材级深度）
│   ├── 1-1-FEA-Elasticity.md        # FEM 概述 + 弹性力学基础
│   ├── 1-2-Tensor-Elasticity.md     # 张量分析 + 弹性力学张量形式
│   ├── 1-3-Variational-Methods.md   # 变分法 + Euler 方程
│   ├── 1-4-FEM-Theory.md            # Ritz 法 + 加权残量法 + FEM 起源
│   ├── 1-5-FEM-Formulation.md       # CST 单元刚度矩阵 + 总刚组装
│   └── 1-6-Element-Construction.md  # 形函数 + 等参元 + Gauss 积分
│
├── 02-Concepts-Formulas/      ← 概念术语表 + 公式参考
│   └── 2-1-Concepts-Glossary.md     # 70+ 条定义、对比表
│
├── 03-Problem-Templates/      ← 题型分类 + 解题策略
│   └── 3-1-Solving-Guide.md        # 5 类题型、分步解题模板
│
├── 04-Homework-Solutions/     ← 作业题目与详解
│   ├── 2026w/                       # 本学年（题目+详解配对）
│   │   ├── HW1-Problem.md          # 作业题目
│   │   ├── HW1-Solution.md         # 完整解答
│   │   ├── HW2-Problem.md          # 作业题目
│   │   ├── HW2-Solution.md         # 完整解答
│   │   ├── HW3-Problem.md          # 作业题目
│   │   ├── HW3-Solution.md         # 完整解答
│   │   ├── HW4-Problem.md          # 作业题目
│   │   └── HW4-Solution.md         # 完整解答
│   └── past/                        # 往年作业 — 参考与补充练习
│       ├── HW2/                     # 往年 HW2 题目 + LIU Sai 答案
│       ├── HW3/                     # 往年 HW3（两版）+ LIU Sai 答案
│       └── comprehensive/           # HowardWoolley 完整答案本 + HW1.1/1.2
│
├── 05-Exam-Cram/              ← 考前公式速查
│   └── 5-1-Formula-Sheet.md       # 一页纸公式速查 + 证明模板
│
├── 06-References/            ← 参考资料
│   ├── 往年参考答案-LIU-Sai.md       # 往年作业答案 + 概念定义
│   └── 变分原理名词简答-宫婷.md      # 术语简答题集
│
└── README.md                  ← 本文件
```

---

## 复习策略

资料按三个阶段组织，由浅入深、循序渐进。

### 第一阶段：系统理解

按顺序通读六章讲义笔记。每章笔记从课件中提炼核心概念、关键公式和重要推导过程。

| 步骤 | 主题 | 核心目标 |
|------|------|---------|
| ① | [FEM 与弹性力学](01-Lecture-Notes/1-1-FEA-Elasticity.md) | 理解 FEM 的本质与必要性；掌握弹性力学三大基本方程（几何、本构、平衡）的矩阵与张量形式；理解最小势能原理 |
| ② | [张量分析](01-Lecture-Notes/1-2-Tensor-Elasticity.md) | 掌握 Einstein 求和约定与指标记法；理解 Kronecker δ 和 Levi-Civita 置换符号的恒等式；掌握张量在坐标旋转下的变换律 |
| ③ | [变分法](01-Lecture-Notes/1-3-Variational-Methods.md) | 定义泛函、变分、线性泛函；从 δQ=0 推导 Euler-Lagrange 方程；推广到高阶导数和多变量情形 |
| ④ | [FEM 理论](01-Lecture-Notes/1-4-FEM-Theory.md) | 理解 Ritz 法：试函数→能量极小→线性方程组；理解 Galerkin 法：以基函数为权函数的加权残量法；区分 FEM 的三种理解途径 |
| ⑤ | [FEM 公式推导](01-Lecture-Notes/1-5-FEM-Formulation.md) | 逐步推导 CST（常应变三角形）单元刚度矩阵；从单元贡献组装总体刚度矩阵；引入边界条件并求解 [K]{δ}={F} |
| ⑥ | [单元构造](01-Lecture-Notes/1-6-Element-Construction.md) | 构造 Lagrange（仅位移）和 Hermite（位移+转角）形函数；理解等参元概念与 Jacobian 变换；应用 Gauss 数值积分 |

**第一阶段目标**：能用自己的语言解释每个核心概念，并能从记忆中写出关键公式。

### 第二阶段：解题训练

将理论应用于具体题目。以作业详解为主线，往年参考答案提供补充练习。

| 步骤 | 资料 | 重点 |
|------|------|------|
| ① | [解题指南](03-Problem-Templates/3-1-Solving-Guide.md) | 熟悉 5 类题型及各自的解题策略 |
| ② | [HW1 解答](04-Homework-Solutions/2026w/HW1-Solution.md) | 张量恒等式证明（ε-δ 恒等式）、指标缩并练习 |
| ③ | [HW2 解答](04-Homework-Solutions/2026w/HW2-Solution.md) | Euler 方程推导与 ODE 求解；Lagrange 乘子法求条件极值 |
| ④ | [HW3 解答](04-Homework-Solutions/2026w/HW3-Solution.md) | Galerkin 法试函数合法性验证；弹性地基梁；Hermite 梁单元形函数 |
| ⑤ | [HW4 解答](04-Homework-Solutions/2026w/HW4-Solution.md) | CST 单元刚度矩阵计算；总刚组装；边界条件引入与求解 |
| ⑥ | [往年参考答案](06-References/往年参考答案-LIU-Sai.md) | 补充练习：Ritz vs Galerkin 数值算例、概念定义及评分要点 |

**练习建议**：先独立尝试每道题，再对照解答。重点关注：
- 分部积分的应用位置及边界项消失的条件
- 本质边界条件与自然边界条件的处理方式
- 非齐次边界条件（如 HW2 Q4 中 y(1)=1）的处理

### 第三阶段：巩固与考前冲刺

将知识整合为紧凑的心智模型。

| 步骤 | 资料 | 重点 |
|------|------|------|
| ① | [概念速查表](02-Concepts-Formulas/2-1-Concepts-Glossary.md) | 系统回顾所有定义，能为每个条目写出简洁定义和一个示例 |
| ② | [公式速查表](05-Exam-Cram/5-1-Formula-Sheet.md) | 背诵一页纸公式速查，能从基本原理重建每个公式 |
| ③ | 证明模板（速查表内） | 练习四个经典推导，做到不看笔记也能写出来 |

---

## 考试范围

> 基于课件截图整理，详见 [exam-scope.md](../exam/exam-scope.md)

### 一、Tensor（张量）

- 标量、向量、张量的定义
- 证明某量是张量
- 特殊符号：Kronecker δ、Levi-Civita 置换符号
- 代数与张量的转换
- 坐标旋转
- 向量在不同坐标系下的变换
- 矩阵在不同坐标系下的变换
- 点积与叉积

### 二、Elasticity（弹性力学）

- Cauchy 公式
- 坐标旋转下的应力变换
- 主应力/主应变
- 应力/应变张量不变量
- 平衡方程
- 几何方程
- Green 应变
- 本构关系
- 各向同性线弹性材料的本构关系
- 边界条件
- 三组方程推导中的假设

### 三、Variational Principles（变分原理）

- 泛函
- 泛函的变分
- 泛函极值
- Euler 方程
- 本质边界条件与自然边界条件
- 泛函的条件极值
- Euler 方程的推广形式（高阶导数、多元函数、多变量）
- 变分法在力学中的应用
- 虚功原理
- 可能位移与可能应力、虚位移与虚应力
- 功的互等定理
- 虚位移原理
- 虚应力原理
- 最小势能原理
- 弹性力学中的 Euler 方程
- 变分问题的直接法与间接法
- 微分方程的等效变分方程
- 有限差分法
- Ritz 法
- Galerkin 法
- 加权残量法

### 四、Basics of FEA（有限元分析基础）

- 一维、三角形、矩形单元的形函数
- 1D 和 2D FEA 公式推导（单元分析、总刚组装、边界条件引入、求解流程）
- 单元、节点、DOF 的概念
- 单元刚度矩阵的特性
- 单元刚度矩阵各元素的物理意义
- FEM 收敛准则
- 一维单元长度坐标
- 一维 Lagrange 线性与高次插值
- Hermite 三次插值
- 三角形单元面积坐标
- 三角形 Lagrange 线性与高次插值
- 划线法
- 矩形单元形函数
- Serendipity 单元形函数
- 等参元
- Jacobian 矩阵

---

## 补充资料

- **课件全文转录**：`courses/FEM/md_output/` — 6 份课件 PDF 完整转换为 Markdown，含所有公式和图片引用
- **补充习题**：`courses/FEM/others/md_output/` — 往年作业题目与答案
- **原始知识库**：`courses/FEM/knowledge_base/` — 本复习资料的原始素材来源

---

## 使用说明

- 文件采用标准 Markdown 格式，数学公式使用 LaTeX 语法（`$$` 显示公式，`$` 行内公式）。可用 VS Code、Typora、GitHub 等任何 Markdown 查看器渲染。
- `§1-5` 等引用指向对应讲义笔记中的章节。
- 本仓库基于 2026 年春季学期课程整理，不同年份的考试范围和侧重点可能有所不同，请以当年课堂内容为准。
- 欢迎纠错和补充，欢迎提 Issue 或 PR。

---

> 整理于 2026 年 6 月 · *上海交通大学 · 变分法与有限元分析*
