## 习题六

1 利用Euler积分公式，推导

$$
\iint_ {O A B} \lambda_ {1} ^ {\alpha_ {1}} \lambda_ {2} ^ {\alpha_ {2}} \lambda_ {3} ^ {\alpha_ {3}} d x d y = \frac {\alpha_ {1} ! \alpha_ {2} ! \alpha_ {3} !}{(\alpha_ {1} + \alpha_ {2} + \alpha_ {3} + 2) !} 2 \Delta_ {e}
$$

2 利用划线法构造二阶三角形单元的插值函数。

3 图所示为二次四边形单元，试计算 $\frac{\partial N_{1}}{\partial x}, \frac{\partial N_{2}}{\partial y}$ 在自然坐标为 $(1/2, 1/2)$ 的点 Q 的值。

提示，因为单元的边是直线，可用四个结点定义单元几何形状的等参变换。

![](exercise_06/images/f5b2e7f7723a3386a61ce674aa57ab43106ebc9379532dac0077fcbd5ae6150d.jpg)

<details>
<summary>text_image</summary>

y
η
1 (40,50)
2(5,40)
Q
ξ
3 (10,10)
4 (30,20)
x
</details>

4 如图所示为一正六面体的弹性体，其位移分量为 $U = a_{1}xyz, V = a_{2}xyz$ ， $W = a_{3}xyz$ ，其中 $a_{1}, a_{2}, a_{3}$ 为常数。若变形前E点坐标为（1.5, 1.0, 2.0）变后移至（1.05.3, 1.001, 1.997），求此时E点的变形状态和E点沿EA方向的线应变。

![](exercise_06/images/0a05f71afa537b654be1dee1c4e76959cba0890b36edcc213ff14459eda82222.jpg)

<details>
<summary>text_image</summary>

z
E
0
2cm
y
A
5cm
1
x
1cm
</details>

5 如图所示，为一个由两根杆组成的结构（二杆分别沿x, y方向）。

结构参数为： $E_{1}=E_{2}=2\times10^{6}Kg/cm^{2},A_{1}=2A_{2}=2cm^{2}$ ，试完成下列有限元分析。

（1）写出各单元的刚度矩阵。  
(2) 写出总刚度矩阵。  
（3）求节点2的位移 $u_{x2}, v_{y2}$ 。  
（4）求各单元的应力。  
(5) 求支反力。

![](exercise_06/images/cd76c0075cbba5b53ba20bdab9f713fc53773a1f9f07bb9bd8a757560022f5fa.jpg)

<details>
<summary>text_image</summary>

y
o x
3
2
10cm
1
2 45°
10cm
p = √2kg
</details>

6 求如图所示平面桁架的节点位移和单元内力。设 $E=2\times10^{6}MPa, A=1cm^{2}$ 。

![](exercise_06/images/4624a345311b2e9a0a85262b494caa250adf94c850d1c6c28cb2bbc20e2110d5.jpg)

<details>
<summary>text_image</summary>

y
V₃
u₃
30cm
2
3
10N
V₁
u₁
x
40cm
2
V₂
u₂
</details>