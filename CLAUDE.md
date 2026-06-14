# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

上海交通大学有限元方法（FEM）课程的结构化复习资料库，2026 年夏季学期。Git 仓库，远程地址 `github.com:YarthsA/FEM-Review.git`。

**无代码文件**——纯 Markdown + LaTeX 数学公式 + PDF 参考资料。

## 目录结构与用途

```
00-Cross-Reference/    双向索引（笔记↔PDF↔作业）
01-Lecture-Notes/      6 章教材级笔记（基于课件 PDF 重写）
02-Concepts-Formulas/  概念术语表 + 易混概念对比表
03-Problem-Templates/  5 类题型解题模板 + 考试重点预测
04-Homework-Solutions/ 2026 HW1-HW4 题目与详解 + 往年作业
05-Exam-Cram/          考前公式速查表（1 页）
06-References/         20 份原始 PDF + 补充参考材料
```

## 内容修改约定

- 所有笔记为 Markdown + `$...$`/`$$...$$` LaTeX 数学公式，编辑时保持格式一致
- 文件命名规则：`序号-标题.md`（如 `1-3-Variational-Methods.md`、`3-1-Solving-Guide.md`）
- 作业文件成对出现：`HW*-Problem.md` + `HW*-Solution.md`，修改时保持同步
- 图片引用使用相对路径，存放在各目录的 `images/` 子目录
- `00-Cross-Reference/index.md` 是全局索引，新增章节或作业后需同步更新

## 课程知识体系（6 章）

1. **FEA 概述与弹性力学** — 应力/应变/位移、三大基本方程
2. **张量与弹性力学张量形式** — 指标记号、Kronecker delta、置换符号
3. **变分法基础** — 泛函、Euler 方程推导、Lagrange 乘子
4. **FEM 理论** — Ritz 法、Galerkin 法、加权残量法
5. **FEM 公式推导** — 杆单元刚度矩阵、总体集成、CST 单元
6. **单元构造** — Lagrange/Hermite 形函数、等参元、Gauss 积分

## 辅助 FEM 解题

当用户发送 FEM 作业/考试题目时，查阅顺序：
1. `03-Problem-Templates/3-1-Solving-Guide.md` — 题型分类与解题模板
2. `02-Concepts-Formulas/2-1-Concepts-Glossary.md` — 概念速查
3. `05-Exam-Cram/5-1-Formula-Sheet.md` — 公式速查
4. 对应章节的 `01-Lecture-Notes/` 笔记 — 详细推导
5. `04-Homework-Solutions/` — 类似题目的解答参考
