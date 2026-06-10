## 一、绪论

有限元法的要点:

(1) 将一个表示结构或连续体的求解域离散为若干个子域 (单元), 并通过它们的边界上的节点相互联结成为组合体;  
(2) 用每个单元内所假设的近似函数来分片地表示全求解域内待求的未知场变量;  
(3) 通过和原问题数学模型（基本方程、边界条件）等效的变分原理或加权余量法，建立求解基本未知量（场函数的节点值）的代数方程组或常微分方程组。

## 二、有限元法的基本理论基础——加权余量法和变分原理

## 1、等效积分形式

对于微分方程和边界条件:

$$
A (u) = \left( \begin{array}{c} A _ {1} (u) \\ A _ {2} (u) \\ \dots \end{array} \right) = 0 (\text {在} \Omega \text {内}), B (u) = \left( \begin{array}{c} B _ {1} (u) \\ B _ {2} (u) \\ \dots \end{array} \right) = 0 (\text {在} \Gamma \text {上})
$$

有

$$
\begin{array}{l} \int_ {\Omega} v ^ {T} A (u) d \Omega \equiv \int_ {\Omega} \left(v _ {1} A _ {1} (u) + v _ {2} A _ {2} (u) + \dots\right) d \Omega \equiv 0 \\ \int_ {\Gamma} \bar {v} ^ {\mathrm{T}} B (u) \mathrm{d} \Gamma \equiv \int_ {\Gamma} \left(\bar {v} _ {1} B _ {1} (u) + \bar {v} _ {2} B _ {2} (u) + \dots\right) \mathrm{d} \Gamma \equiv 0 \\ \end{array}
$$

其中， $v=\begin{pmatrix}v_{1}\\v_{2}\\...\end{pmatrix},\quad\overline{v}=\begin{pmatrix}\overline{v}_{1}\\ \overline{v}_{2}\\...\end{pmatrix}$ 为任意函数。

方程 $\int_{\Omega} v^{\mathrm{T}} A(u) \, \mathrm{d}\Omega + \int_{\Gamma} \bar{v}^{\mathrm{T}} B(u) \, \mathrm{d}\Gamma = 0$ 称为该微分方程边值问题的等效积分形式。

## 2、等效积分弱形式

通过适当提高对任意函数 v 及 $\bar{v}$ 的连续性要求，以降低对微分方程场函数 u 的连续性要求所建立的等效积分形式称为微分方程的等效积分弱形式。

$$
\int_ {\Omega} C ^ {T} (v) D (u) d \Omega + \int_ {\Gamma} E ^ {T} (\bar {v}) F (u) d \Gamma = 0
$$

其中，C、D、E、F 是微分算子，所包含的导数阶数低于 A。

## 3、 $C_{n - 1}$ 连续性

在等效积分形式中，如果微分算子 A 出现的最高阶导数是 n 阶，则要求函数 u 必须有连续的 n-1 阶导数，即 $C_{n-1}$ 连续性。

## 4、加权余量法

采用使余量的加权积分为零来求得微分方程近似解得方法称为加权余量法（weighted residual method，WRM）。余量：

$$
R = A (N a), \quad \overline {{{{R}}}} = B (N a)
$$

$$
\int_ {\Omega} C ^ {T} \left(W _ {j}\right) D (N a) d \Omega + \int_ {\Gamma} E ^ {T} \left(\bar {W} _ {j}\right) F (N a) d \Gamma = 0 (j = 1, \dots , n)
$$

其中， $u \approx \widetilde{u} = \sum_{i=1}^{n} N_i a_i = Na$ ， $W_j$ 和 $\overline{W_j}$ 为权函数。近似函数所取的试探函数项数 n 越多，
近似解的精度越高。n无穷大时，近似解收敛于精确解。

常用权函数:

## (1) 配点法

简单地强迫余量在域内 n 个点上等于零， $W_{j}=\delta(x-x_{j})$ ， $\int_{\Omega}W_{j}d\Omega=I$ 。

## (2) 子域法

强迫余量在 $\mathbf{n}$ 个子域 $\Omega_{j}$ 的积分为零。

## (3) 最小二乘法

$$
W _ {j} = \frac {\partial}{\partial a _ {j}} A \left(\sum_ {i = 1} ^ {n} N _ {i} a _ {i}\right)
$$

$$
I \left(a _ {i}\right) = \int_ {\Omega} A ^ {\mathrm{T}} \left(\sum_ {i = 1} ^ {n} N _ {i} a _ {i}\right) A \left(\sum_ {i = 1} ^ {n} N _ {i} a _ {i}\right) \mathrm{d} \Omega
$$

$$
\frac {\partial I}{\partial a _ {i}} = 0
$$

## (4) 力矩法

强迫余量的各次矩等于零， $W_{j}=1,x,x^{2},\ldots$ ，与 $A(Na)$ 内积为零。

## 5、伽辽金法

简单地利用近似解得试探函数序列作为权函数， $W_{j}=N_{j}$ ， $\overline{W}_{j}=-W_{j}=-N_{j}$

$$
\begin{array}{l} \int_ {\Omega} N _ {j} ^ {\mathrm{T}} A \left(\sum_ {i = 1} ^ {n} N _ {i} a _ {i}\right) \mathrm{d} \Omega - \int_ {\Gamma} N _ {j} ^ {\mathrm{T}} B \left(\sum_ {i = 1} ^ {n} N _ {i} a _ {i}\right) \mathrm{d} \Gamma = 0 \\ \delta \widetilde {u} = N _ {1} \delta a _ {1} + N _ {2} \delta a _ {2} + \dots + N _ {n} \delta a _ {n} \\ \end{array}
$$

等效积分形式： $\int_{\Omega}\delta \widetilde{u}^{\mathrm{T}}A(\widetilde{u})\mathrm{d}\Omega -\int_{\Gamma}\delta \widetilde{u}^{\mathrm{T}}B(\widetilde{u})\mathrm{d}\Gamma = 0$

等效积分弱形式： $\int_{\Omega}C^{\mathrm{T}}\bigl (\delta \widetilde{u}\bigr)D(\widetilde{u})\mathrm{d}\Omega -\int_{\Gamma}E^{\mathrm{T}}\bigl (\delta \widetilde{u}\bigr)F(\widetilde{u})\mathrm{d}\Gamma = 0$

如果算子 A 是 2m 阶线性自伴随的，采用伽辽金法求解的方程组系数矩阵是对称的。

## 6、线性、自伴随算子

若微分方程 $L(u)+b=0$ （在 $\Omega$ 内），有 $L(\alpha u_{1}+\beta u_{2})=\alpha L(u_{1})+\beta L(u_{2})$ 称为线性；

与任意函数 v 内积进行分部积分至 u 的导数消失， $\int_{\Omega} L(u) v \, \mathrm{d}\Omega = \int_{\Omega} u L^{*}(v) \, \mathrm{d}\Omega + b.t.(u,v)$ ，且 $L = L^{*}$ 称为自伴随。

## 7、泛函和变分原理

对微分方程边值问题:

$$
A (u) = L (u) + f = 0 \text {(在} \Omega \text {内)}, B (u) = 0 \text {(在} \Gamma \text {上})
$$

伽辽金法： $\int_{\Omega}\delta u^{\mathrm{T}}[L(u) + f]\mathrm{d}\Omega -\int_{\Gamma}\delta u^{\mathrm{T}}B(u)\mathrm{d}\Gamma = 0$

算子线性、自伴随： $\int_{\Omega}\delta u^{\mathrm{T}}L(u)\mathrm{d}\Omega = \delta \left[\int_{\Omega}\frac{1}{2} u^{\mathrm{T}}L(u)\mathrm{d}\Omega \right] + b.t.\bigl (\delta u,u\bigr)$

原问题的变分原理： $\delta \Pi (u) = 0$ ， $\Pi (u) = \int_{\Omega}\left[\frac{1}{2} u^{\mathrm{T}}L(u) + u^{\mathrm{T}}f\right]\mathrm{d}\Omega +b.t.(u)$

微分方程边值问题等效于泛函变分等于零，即泛函取驻值。

泛函可通过等效积分的伽辽金法得到，则称为自然变分原理。

## 8、二次泛函

$\Pi (u) = \int_{\Omega}\left[\frac{1}{2} u^{\mathrm{T}}L(u) + u^{\mathrm{T}}f\right]\mathrm{d}\Omega +b.t.(u)$ 中，u的最高次为二次，称为二次泛函。

## 9、强制边界条件

在选择函数 u 时已满足的边界条件称为强制边界条件。

## 10、自然边界条件

在形成等效积分弱形式过程中（分部积分法），分部积分产生的边值项为零所形成的边界条件称为自然边界条件。通常为导数形式。

## 11、泛函的驻值和极值

对于 2m 阶线性、自伴随的微分方程，通过伽辽金法弱形式建立的变分原理，只有在近似场函数事先满足强制边界条件的情况下，才可能使泛函具有极值性。

## 12、里兹法

$u \approx \widetilde{u} = \sum_{i=1}^{n} N_i a_i = Na$ 代入泛函 $\Pi$ :

$$
\delta \Pi = \frac {\partial \Pi}{\partial a _ {1}} \delta a _ {1} + \frac {\partial \Pi}{\partial a _ {2}} \delta a _ {2} +... + \frac {\partial \Pi}{\partial a _ {n}} \delta a _ {n}, \text {由于} \delta a _ {i} \text {的任意性，} \frac {\partial \Pi}{\partial a} = \left[ \begin{array}{l} \frac {\partial \Pi}{\partial a _ {1}} \\ \frac {\partial \Pi}{\partial a _ {2}} \\ ... \\ \frac {\partial \Pi}{\partial a _ {n}} \end{array} \right] = 0
$$

通过待定参数 a 的方程组求解 a 的方法，称为里兹法。

当问题存在自然变分原理时，里兹法与伽辽金法所得的结果相同。

里兹法的收敛性:

（1）试探函数 N 具备完备性；  
（2）试探函数 N 满足 $C_{n-1}$ 连续性，即协调性。

## 13、虚功原理

虚功原理是平衡方程与几何方程的等效积分弱形式，其表述为：

变形体中任意满足平衡的利息在任意满足协调条件的变形状态上做的虚功等于零, 即体系外力虚功与内力虚功之和等于零。

虚功原理是虚位移原理和虚应力原理的总称。

虚功原理不包含物理方程，所以可用于非线性弹性问题和弹塑性问题。但几何方程和平衡方程都是基于小变形理论，不能直接用于大变形问题。

## 14、虚位移原理

平衡方程： $\sigma_{ij,j} + f_i = 0$ （在 V 内）

力边界条件： $\sigma_{ij}n_{j}-\overline{T}_{i}=0$ （在 $S_{\sigma}$ 上）

等效积分形式： $\int_{V}\delta u_{i}\bigl (\sigma_{ij,j} - f_{i}\bigr)\mathrm{dV} - \int_{S_{\sigma}}\delta u_{i}\bigl (\sigma_{ij}n_{j} - \overline{T}_{i}\bigr)\mathrm{d}S = 0$

等效积分弱形式： $\int_{V}\left(-\delta \varepsilon_{ij}\sigma_{ij} + \delta u_i f\right)\mathrm{dV} + \int_{S_\sigma}\delta u_i\overline{T}_i\mathrm{d}S = 0$

力学意义：如果力系是平衡的，则它们在虚位移和虚应变上做的功总和为零。反之，如果力系在虚位移和虚应变上做的功总和为零，则它们一定是平衡的。

## 15、虚应力原理

几何方程： $\varepsilon_{ij}=\frac{1}{2}\left(u_{i,j}+u_{j,i}\right)$ （在 V 内）

位移边界条件： $u_{i}=\bar{u}_{i}$ （在 $S_{u}$ 上）

等效积分形式： $\int_{V}\delta \sigma_{ij}\left(\varepsilon_{ij} - \frac{1}{2}\big(u_{i,j} + u_{j,i}\big)\right)\mathrm{dV} + \int_{S_u}\delta T_i\big(u_i - \overline{u}_i\big)\mathrm{d}S = 0$

简化为： $\int_{V}\delta \sigma_{ij}\varepsilon_{ij}\mathrm{dV} - \int_{S_u}\delta T_i\overline{u}_i\mathrm{d}S = 0$

力学意义：如果位移是协调的，则虚应力和虚边界约束力在它们上面所做的功总和为零。反之，如果虚力系在它们上面做的功之和为零，则它们一定是协调的。

## 16、最小位能原理

虚位移原理： $\int_V\Bigl (\delta \varepsilon_{ij}D_{ijkl}\varepsilon_{kl} - \delta u_if_i\bigr)\mathrm{dV} - \int_{S_\sigma}\delta u_i\overline{T}_i\mathrm{d}S = 0$

体积应变能的变分： $\left(\delta \varepsilon_{ij}\right)D_{ijkl}\varepsilon_{kl} = \delta \left(\frac{1}{2} D_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\right) = \delta U\left(\varepsilon_{mn}\right)$

系统总位能： $\Pi_{\mathrm{p}} = \Pi_{\mathrm{p}}(u_i) = \int_V\left(\frac{1}{2} D_{ijkl}\varepsilon_{ij}\varepsilon_{kl} - u_if_i\right)\mathrm{dV} - \int_{S_\sigma}u_i\overline{T}_i\mathrm{d}S = 0$

真实位移使系统总位能取最小值： $\delta\Pi_{p}=0$

最小位能原理求得近似解的弹性变性能是精确解变形能的下界, 即近似的位移场在总体上偏小, 结构的计算模型偏于刚硬。

## 17、最小余能原理

虚应力原理： $\int_{V}\delta \sigma_{ij}C_{ijkl}\sigma_{kl}\mathrm{dV} - \int_{S_u}\delta T_i\overline{u}_i\mathrm{d}S = 0$

体积应变余能的变分： $\left(\delta \sigma_{ij}\right)C_{ijkl}\sigma_{kl} = \delta \left(\frac{1}{2} C_{ijkl}\sigma_{ij}\sigma_{kl}\right) = \delta V\left(\sigma_{mn}\right)$

系统总余能： $\Pi_{\mathrm{c}} = \Pi_{\mathrm{c}}\left(\sigma_{ij}\right) = \int_{V}\frac{1}{2} C_{ijkl}\sigma_{ij}\sigma_{kl}\mathrm{dV} - \int_{S_u}u_i\overline{T}_i\mathrm{d}S = 0$

真实应力使系统总余能取最小值： $\delta\Pi_{c}=0$

最小余能原理得到的应力近似解得弹性与能是精确解余能的上界, 即近似的应力解在总体上偏大, 结构的计算模型偏于柔软。

## 三、弹性力学问题有限元方法的一般原理和表达格式

## 1、广义坐标

单元位移函数采用多项式近似，待定系数称为广义坐标。

## 2、位移模式

位移模式即单元位移函数的基函数。

## 3、位移插值函数

位移函数表示为节点位移的函数，称为位移插值函数或形函数。

位移插值函数的性质:

（1）节点上具有 Kronecker delta 性质，即某一节点的形函数在自身节点上为 1，其余节点上为零；

（2）在单元中任意点各插值函数之和应等于1，反映刚体位移；

（3）插值函数是线性的，则单元内部及单元边界上的位移也是线性的。

## 4、单元刚度矩阵及刚度系数

$$
K ^ {e} = \int_ {\Omega^ {e}} B ^ {\mathrm{T}} D B t \mathrm{d} x \mathrm{d} y
$$

单元刚度矩阵中任一元素 $K_{ij}$ 的物理意义：当单元的第j个节点位移为单位位移而其他节点位移为零时，需在单元第i个节点位移方向上施加的节点力的大小。

## 5、单元刚度矩阵的对称性和奇异性

（1）对称性

(2) 奇异性: 刚度矩阵中不是所有行都是独立的。

物理意义: 通过节点位移可以计算出单元的节点力, 但不能通过节点荷载求得单元节点位移, 可能有刚体运动。

（3）主元恒正：在节点上施加的节点力与节点位移方向一致。

## 6、单元节点荷载列阵

$$
P ^ {e} = P _ {f} ^ {e} + P _ {S} ^ {e} = \int_ {\Omega^ {e}} N ^ {T} f t \mathrm{d} x \mathrm{d} y + \int_ {S _ {\sigma} ^ {e}} N ^ {T} T t \mathrm{d} S
$$

## 7、结构刚度矩阵的集成

$$
K = \sum_ {e} G ^ {\mathrm{T}} K ^ {e} G
$$

结构刚度矩阵的性质:

(1) 对称性

(2) 奇异性: 引入位移边界条件后奇异性消失

（3）稀疏性：每个节点的相关单元只是围绕在该节点周围为数较少的几个。

（4）非零元素呈带状分布：只要节点编号合理，稀疏的非零元素将集中在以主对角线为中心的带状区域。

引入边界条件的方法：直接代入法、对角元素该1法（对0位移）、对角元素乘大数法

## 8、有限元解的收敛准则

有限元可看做是里兹法的一种特殊形式，试探函数定义于单元而不是全域。

两类误差:

（1）离散误差：结构离散过程中，单元试探函数近似整体域的场函数所引起的误差；  
（2）舍入误差：计算机有限的有效位数引起的。

## 9、位移元的完备性和协调性

（1）完备性要求：如果泛函中场函数的最高阶导数为 m 阶，则试探函数至少是 m 次完全多项式；  
（2）协调性要求：如果泛函中场函数的最高阶导数为 $\mathrm{m}$ 阶，则试探函数在单元交界面上必须具有 $C_{n - 1}$ 连续性。

## 10、位移元解得下限性质

位移元得到的位移解总体上不大于精确解。

单元原是连续体的一部分，具有无限多个自由度。在假定了单元的位移函数后，自由度限制为只有以节点位移表示的有限自由度，即位移函数对单元的变形形成了约束和限制，使单元刚度较连续体加强了。

## 11、弹性力学问题有限元分析的步骤

（1）对结构进行离散，按问题的几何特点和精度要求等进行单元划分并形成网格  
(2) 形成单元刚度矩阵和等效节点荷载列阵  
（3）集成结构的刚度矩阵和等效节点荷载列阵  
（4）引入强制边界条件（给定位移）  
（5）求解有限元方程（线性代数方程），得到节点位移  
(6) 计算单元应变和应力  
(7) 进行必要的后处理

## 四、单元和插值函数的构造

## 1、自然坐标

无量纲的局部坐标，只与单元形状相关。

## 2、面积坐标

三角形内任意点到三边的距离与三个顶点到对边的距离之比。

特点:

（1）三角形三个顶点坐标为： $i(1,0,0)$ ， $j(0,1,0)$ ， $m(0,0,1)$  
（2）三个面积坐标并不独立，必然满足 $L_{i} + L_{j} + L_{m} = 1$

## 3、体积坐标

三棱锥内任意点到四个面的距离与四个顶点到对面的距离之比。

## 4、拉格朗日单元

采用拉格朗日差值多项式作为插值函数的单元。

$$
\phi = \sum_ {i = 1} ^ {n} N _ {i} \phi_ {i}, \quad N _ {i} (x _ {j}) = \delta_ {i j}, \quad \sum_ {i = 1} ^ {n} N _ {i} (x) = 1
$$

拉格朗日矩形单元随着插值函数方次的增加，内节点增加，从而增加了单元的自由度，而这些自由度通常不能提高单元的精度。

## 5、Serendipity 单元

对于高次单元，减少单元内部的节点数，而保留边界上的节点数。

四次及以上的 Serendipity 单元需要增加内部节点才能保证相应次多项式的完备性。

## 五、等参元

## 1、等（超、次）参变化

如果坐标变换和函数产值采用相同的节点，并且采用相同的插值函数，则称这种变换为等参变换。

如果坐标变换节点数多于函数插值的节点数，则称为超参变换。

如果坐标变换节点数少于函数插值的节点数，则称为次（亚）参变换。

## 2、雅可比矩阵和行列式

$$
\left[ \begin{array}{c} \frac {\partial N _ {i}}{\partial \xi} \\ \frac {\partial N _ {i}}{\partial \eta} \\ \frac {\partial N _ {i}}{\partial \zeta} \end{array} \right] = \left[ \begin{array}{c c c} \frac {\partial x}{\partial \xi} & \frac {\partial y}{\partial \xi} & \frac {\partial z}{\partial \xi} \\ \frac {\partial x}{\partial \eta} & \frac {\partial y}{\partial \eta} & \frac {\partial z}{\partial \eta} \\ \frac {\partial x}{\partial \zeta} & \frac {\partial y}{\partial \zeta} & \frac {\partial z}{\partial \zeta} \end{array} \right] \left[ \begin{array}{c} \frac {\partial N _ {i}}{\partial x} \\ \frac {\partial N _ {i}}{\partial y} \\ \frac {\partial N _ {i}}{\partial z} \end{array} \right] = J \left[ \begin{array}{c} \frac {\partial N _ {i}}{\partial x} \\ \frac {\partial N _ {i}}{\partial y} \\ \frac {\partial N _ {i}}{\partial z} \end{array} \right]
$$

$$
J \equiv \frac {\partial (x , y , z)}{\partial (\xi , \eta , \zeta)} = \left[ \begin{array}{c c c} \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \xi} x _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \xi} y _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \xi} z _ {i} \\ \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \eta} x _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \eta} y _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \eta} z _ {i} \\ \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \zeta} x _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \zeta} y _ {i} & \sum_ {i = 1} ^ {n} \frac {\partial N _ {i} ^ {\prime}}{\partial \zeta} z _ {i} \end{array} \right]
$$

## 3、等参变换的条件

$$
\left| J \right| \neq 0
$$

$$
\left| J \right| = \frac {\left| \mathrm{d} \xi \right| \left| \mathrm{d} \eta \right| \sin (\mathrm{d} \xi , \mathrm{d} \eta)}{\mathrm{d} \xi \mathrm{d} \eta}
$$

（1）防止任意的两个节点退化为一个节点： $|d\xi|,|d\eta|$ 不为零；

(2) 防止单元过分歪曲导致 $d\xi, d\eta$ 共线。

## 4、等参元的收敛性

协调性、完备性

## 六、有限元法应用中的若干实际考虑

## 1、单元选择

三角形单元：形状简单、随意性大、适应区域形状能力强。

矩形单元：插值函数构造差异性大，但精度较高。

## 2、网格划分

一般情况下，为了使结果达到必要的精度，可以采用以下措施：

（1）对于应力变化激烈的区域局部加密网格进行重分析；  
（2）采用自适应分析方法。

## 3、应力解的近似性

应力解的近似性表现在:

（1）单元内部一般不满足平衡方程；  
（2）单元与单元的交界面应力一般不连续；  
（3）在力边界上一般也不满足力的边界条件。

## 七、约束变分原理

## 1、约束变分原理

引入附加条件构造修正泛函的变分原理。

## 2、拉格朗日乘子法

$\Pi^{*} = \Pi +\int_{\Omega}\lambda^{\mathrm{T}}C(u)\mathrm{d}\Omega ,\lambda$ 为域内一组独立坐标的函数向量

$$
\delta \Pi^ {*} = \delta \Pi + \int_ {\Omega} \delta \lambda^ {T} C (u) d \Omega + \int_ {\Omega} \lambda^ {T} \delta C (u) d \Omega
$$

直接使用的两个问题:

（1）方程组的阶数随附加条件的增加而增加，从而增加了计算量；  
（2）对 $\lambda$ 求偏导得到的方程中必然不包含 $\lambda$ ，方程组系数矩阵必定存在零对角元素。

## 3、罚函数法

$\Pi^{**}=\Pi+\alpha\int_{\Omega}C^{\mathrm{T}}(u)C(u)\mathrm{d}\Omega,\quad\alpha$ 称为罚参数。

利用罚函数求解条件驻值问题不增加未知参量的个数，并且不改变驻值的性质。