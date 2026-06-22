## 习题四

1 对于以下方程: $\frac{d^2\varphi}{dx^2} + \varphi = X$

边界条件为 $\varphi (x = 0) = 0,\varphi (x = 1) = 1$

试推导出与它等效的泛函。若采用近似函数 $\varphi (X) = a_0 + a_1X + a_2X^2$

求解时，试用泛函极值的方法求解待定参数 $a_{0}$ ， $a_{1}$ ， $a_{2}$ 。

2 若函数 $Y(X)$ 的二次泛函为

$$
\Pi [ y ] = \int_ {x _ {1}} ^ {x _ {2}} \Big [ p _ {(x)} y ^ {2} + 2 q _ {(x)} y y ^ {\prime} + r _ {(x)} (y ^ {\prime}) ^ {2} + 2 f _ {(x)} y + 2 g (x) y ^ {\prime} \Big ] d x
$$

试证明所对应的控制方程（即欧拉方程）为

$$
(r y ^ {\prime}) ^ {\prime} + (q ^ {\prime} - p) y + g ^ {\prime} - f = 0
$$

3 证明 $w = a(1 - \cos \frac{\pi x}{2l})$ 不可用作下图悬臂梁问题的Galerkin法试函数。

![](exercise_04/images/4d50149193b4256bcb454398c7cfc278b2983f6e4cf8adef67d9631c7608f122.jpg)

<details>
<summary>text_image</summary>

l
x
y
</details>

4 利用变分法推导图示弹性地基梁的微分方程及边界条件。并利用Galerkin法取一阶近似进行求解。

![](exercise_04/images/5b04ea4003513cb3f499c4b9cc533955e1b309564f2bcb94cd03f75f9647b36b.jpg)

<details>
<summary>text_image</summary>

p
x
y
l/2
l/2
</details>

5 采用课件中给出二阶试函数，任选两种加权残值法计算下图所示的简支梁问题。

![](exercise_04/images/21326fbfe2dc8f5fe8ae1a0d548baff31085b23277b50bd851fd201d49240f3d.jpg)

<details>
<summary>text_image</summary>

y
l
x
</details>

6 设横截面面积为常数的弹性杆两端固定，杆长为3L，弹性杆各处受相同的体积力f作用。采用3个长为L的线性元，试给出单元形函数矩阵、单元刚度阵、总体刚度矩阵。

![](exercise_04/images/79106b9876ddbeadb7f3bbd21cb009adcc46e4d8b7faa804190dd9711c4fd793.jpg)

<details>
<summary>text_image</summary>

1 2 3 4
1 2 3
L L L
</details>