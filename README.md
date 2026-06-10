# Variation Method & Finite Element Analysis — Course Review Guide

> **Shanghai Jiao Tong University · FEM Course (Civil Engineering)**
> A systematically organized review repository based on lecture notes, homework assignments, and past reference solutions.

---

## Repository Structure

```
FEM-Review/
├── 01-Lecture-Notes/          ← Chapter summaries distilled from course slides
│   ├── 1-1-FEA-Elasticity.md        # FEA overview + fundamentals of elasticity
│   ├── 1-2-Tensor-Elasticity.md     # Tensor notation + elasticity in tensor form
│   ├── 1-3-Variational-Methods.md   # Variational calculus + Euler equations
│   ├── 1-4-FEM-Theory.md            # Ritz method + weighted residuals + FEM origins
│   ├── 1-5-FEM-Formulation.md       # CST element stiffness + global assembly
│   └── 1-6-Element-Construction.md  # Shape functions + isoparametric + Gauss quadrature
│
├── 02-Concepts-Formulas/      ← Glossary + formula reference
│   └── 2-1-Concepts-Glossary.md     # 70+ definitions, comparison tables
│
├── 03-Problem-Templates/      ← Problem taxonomy + solution strategies
│   └── 3-1-Solving-Guide.md        # 5 problem types, step-by-step templates
│
├── 04-Homework-Solutions/     ← Problems + solutions (current & past)
│   ├── 2026w/                       # Current year (problem + solution pairs)
│   │   ├── HW1-Problem.md          # Problem statement
│   │   ├── HW1-Solution.md         # Full solution
│   │   ├── HW2-Problem/Solution.md
│   │   └── HW3-Problem/Solution.md
│   └── past/                        # Past years — for reference & extra practice
│       ├── HW2/                     # Old HW2 problem + LIU Sai's solution
│       ├── HW3/                     # Old HW3 (two versions) + LIU Sai's solution
│       └── comprehensive/           # HowardWoolley's full answer book + HW1.1/1.2
│
├── 05-Exam-Cram/              ← Condensed reference for last-minute review
│   └── 5-1-Formula-Sheet.md       # One-page formula summary + proof templates
│
├── 06-References/            ← Historical reference materials
│   ├── 往年参考答案-LIU-Sai.md       # Past homework solutions with concept definitions
│   └── 变分原理名词简答-宫婷.md      # Terminology Q&A collection
│
└── README.md                  ← This file
```

---

## Review Strategy

The material is organized into three sequential phases. Each phase builds on the previous one; mastery of earlier topics is assumed in later sections.

### Phase I: Systematic Understanding

Work through the six lecture-note chapters in order. Each note is a structured summary that extracts the core concepts, key formulas, and important procedures from the original course slides.

| Step | Topic | Core Objectives |
|------|-------|-----------------|
| ① | [FEA & Elasticity](01-Lecture-Notes/1-1-FEA-Elasticity.md) | — Understand what FEM is and why it is needed<br/>— Know the three fundamental equations of elasticity (geometric, constitutive, equilibrium) in both matrix and tensor forms<br/>— Understand the principle of minimum potential energy |
| ② | [Tensor Analysis](01-Lecture-Notes/1-2-Tensor-Elasticity.md) | — Master Einstein summation convention and index notation<br/>— Understand Kronecker delta ($\delta_{ij}$) and permutation symbol ($e_{ijk}$) identities<br/>— Know tensor transformation laws under coordinate rotation |
| ③ | [Variational Calculus](01-Lecture-Notes/1-3-Variational-Methods.md) | — Define functionals, variation, linear functionals<br/>— Derive the Euler-Lagrange equation from $\delta Q=0$<br/>— Generalize to higher-order derivatives and multiple independent functions |
| ④ | [FEM Theory](01-Lecture-Notes/1-4-FEM-Theory.md) | — Understand Ritz method: trial functions → energy minimization → linear system<br/>— Understand Galerkin method: weighted residuals with basis functions as weights<br/>— Distinguish the three perspectives on FEM (matrix method, variational, weighted residuals) |
| ⑤ | [FEM Formulation](01-Lecture-Notes/1-5-FEM-Formulation.md) | — Derive the CST (constant strain triangle) element stiffness matrix step by step<br/>— Assemble the global stiffness matrix from element contributions<br/>— Handle boundary conditions and solve $[K]\{\delta\}=\{F\}$ |
| ⑥ | [Element Construction](01-Lecture-Notes/1-6-Element-Construction.md) | — Construct Lagrange (nodal displacement only) and Hermite (displacement + rotation) shape functions<br/>— Understand the isoparametric concept and Jacobian transformation<br/>— Apply Gauss quadrature for numerical integration |

**Expected outcome after Phase I**: You can explain each core concept in your own words and write down the key formulas from memory.

### Phase II: Problem Solving

Apply the theory to concrete problems. The homework solutions are the primary vehicle; past reference solutions provide additional practice and expose alternative solution styles.

| Step | Materials | Focus |
|------|-----------|-------|
| ① | [Problem-solving guide](03-Problem-Templates/3-1-Solving-Guide.md) | Familiarize yourself with the 5 problem categories and the recommended strategy for each |
| ② | [HW1 solution](04-Homework-Solutions/HW1-Solution.md) | Tensor identity proofs ($\varepsilon$-$\delta$ identities), index contraction exercises |
| ③ | [HW2 solution](04-Homework-Solutions/HW2-Solution.md) | Euler equation derivation and ODE solving; Lagrange multiplier method for constrained extrema |
| ④ | [HW3 solution](04-Homework-Solutions/HW3-Solution.md) | Galerkin method trial-function validity; elastic foundation beam; Hermite beam element shape functions |
| ⑤ | [Past reference solutions](06-References/往年参考答案-LIU-Sai.md) | Additional practice: Ritz vs Galerkin numerical examples, concept definitions with scoring hints |

**Recommended practice routine**: Attempt each problem independently before consulting the solution. Compare your derivation step-by-step with the solution, paying special attention to:
- Where integration by parts is applied and how boundary terms vanish
- How boundary conditions are incorporated (essential vs natural)
- The handling of non-homogeneous boundary conditions (e.g., $y(1)=1$ in HW2 Q4)

### Phase III: Consolidation & Exam Preparation

Synthesize the material into a compact mental model.

| Step | Materials | Focus |
|------|-----------|-------|
| ① | [Concepts glossary](02-Concepts-Formulas/2-1-Concepts-Glossary.md) | Review all definitions systematically. Aim to produce a concise definition and one illustrative example for each entry |
| ② | [Formula sheet](05-Exam-Cram/5-1-Formula-Sheet.md) | Commit the one-page summary to memory. Verify you can reconstruct each formula from first principles |
| ③ | Proof templates (within Formula Sheet) | Practice the four canonical derivations until they can be produced without notes |

---

## Exam Topic Analysis

Based on lecture content and past exam patterns, the following hierarchy of importance emerges:

### Tier 1 — Core (appear near-universally)
- **CST element stiffness matrix calculation**: full pipeline $[N] \to [B] \to [k]_e = t\Delta_e[B]^T[D][B]$
- **Euler-Lagrange equation derivation**: $\delta Q=0 \to$ integration by parts $\to$ fundamental lemma $\to$ $F_y - \frac{d}{dx}F_{y'} = 0$
- **Ritz method for functional extremum**: trial function selection, energy minimization, linear system solution
- **Conceptual comparisons**: Ritz vs Galerkin, displacement element lower-bound property, isoparametric/superparametric/subparametric elements

### Tier 2 — High frequency
- Tensor identities ($\varepsilon$-$\delta$ relations, proving tensor character of a quantity)
- Weighted residual method family (Galerkin, least squares, collocation, subdomain, moment — accuracy comparison)
- Shape function construction (Lagrange vs Hermite interpolation)
- Global stiffness matrix assembly procedure
- Principle of minimum potential energy and principle of virtual work

### Tier 3 — Moderate frequency
- Isoparametric element concept and Jacobian matrix
- Gauss quadrature: point locations, weights, and order of accuracy
- Strain energy, complementary energy, total potential energy
- Essential vs natural boundary conditions

---

## Common Pitfalls

| Pitfall | Correction | Reference |
|---------|-----------|-----------|
| Element node numbering order | Must be **counter-clockwise**; otherwise element area $\Delta_e$ becomes negative | §1-5 |
| Solving $[K]\{\delta\}=\{F\}$ without boundary conditions | The global stiffness matrix is singular until rigid-body modes are eliminated | §1-5 |
| Confusing plane stress with plane strain | Thin plate → plane stress ($\sigma_z=0$); thick/long structure → plane strain ($\varepsilon_z=0$) | §1-5 |
| Ritz vs Galerkin trial-function requirements | Ritz: displacement BC only; Galerkin: **all** BCs (displacement + traction) | §1-4 |
| Jacobian determinant changing sign | Element is excessively distorted; remesh | §1-6 |
| Testing linearity of a functional | Must satisfy **both** homogeneity $Q[cy]=cQ[y]$ **and** additivity $Q[y_1+y_2]=Q[y_1]+Q[y_2]$ | §1-3 |
| Confusing $w'''(l)=0$ (shear-free) condition | Verify all boundary conditions, not just the displacement ones, when assessing trial-function admissibility | §1-4 |

---

## Supplementary Materials

- **Full transcripts from slides**: `courses/FEM/md_output/` — complete conversion of 6 lecture PDFs to Markdown, including all formulas and image references
- **Additional exercises**: `courses/FEM/others/md_output/` — past homework papers and solutions
- **Original knowledge base**: `courses/FEM/knowledge_base/` — the working knowledge base from which this review was distilled

---

## Notes on Use

- Files are formatted in standard Markdown with LaTeX math delimiters (`$$` for display, `$` for inline). Render with any Markdown-compatible viewer (VS Code, Typora, GitHub).
- References to `§1-5` etc. refer to sections within the corresponding lecture-note file.
- This repository reflects the Spring 2026 offering of the course. Syllabi and emphasis may vary between years; consult the current lecture for the definitive topic list.
- Corrections and additions are welcome. Please open an issue or submit a PR.

---

> Compiled June 2026 · *Variation Method & FEA — SJTU*
