## 习题五

1 利用有限元法求解方程

$$
\left\{ \begin{array}{c} - u ^ {\prime \prime} (x) + u (x) = - x \\ u (0) = u (1) = 0 \end{array} \right.
$$

要求：用结点0，1/4，1/2，3/4，1将区间分为四个单元。

提示:

该方程的等效若积分形式如下:

$$
Q [ u (x) ] = \frac {1}{2} \int_ {0} ^ {1} \left[ u ^ {\prime 2} (x) + u ^ {2} (x) + 2 x u (x) \right] d x
$$

2 针对垂直外力作用下的薄膜平衡问题，用变分的方法导出其对应的Poisson方程。  
3 推导线性三角形单元的单元刚度矩阵。

4 推导图示矩形单元的单元刚度矩阵。

![](exercise_05/images/10d85f05df99fed56a8c5604feb788b8a30d9efb4930a27c965c6017785a0ac3.jpg)

<details>
<summary>line chart</summary>

| Point | x | y |
|---|---|---|
| 1 | 0 | 1 |
| 2 | 2 | 1 |
| 3 | 2 | 1 |
The chart displays a rectangular region defined by the line y = x, with labels indicating the coordinates of x and y for each point. The x-axis ranges from 0 to 2 and the y-axis ranges from 0 to 1. The label 'y' appears in the top left corner.
</details>

5 试证明面积坐标与直角坐标满足下列转换关系。

$$
\begin{array}{l} X = x _ {i} L _ {i} + x _ {j} L _ {j} + x _ {m} L _ {m} \\ y = y _ {i} L _ {i} + y _ {j} L _ {j} + y _ {m} L _ {m} \\ \end{array}
$$

6 证明二维平行四边形单元的Jacobi矩阵是常数。

7 对于承受轴对称载荷的回转体，若取3节点三角形环形单元，试求:

(1) 以转速 $\omega$ 旋转时节点的等效载荷。  
(2) 若回转轴方向有 $\alpha_{z}$ 的加速度时，如何计算节点的等效荷载。