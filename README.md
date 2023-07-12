# ACSwithPatternsearch
以求下面这个公式的最小值为例子
$$\sum_{i=1}^{3}[(x_{i+1}-{x_i}^2)+(x_i-1)^2]  x_i \in [-30,30]$$
算法描述
1. 蚂蚁随机分布在四维空间，根据位置赋予初始信息素 ${T_0}$
2. 对每一只蚂蚁 $ant_k$，先运用模式搜索法局部搜索。然后依据概率 $p_{ks}$选择是否向较优蚂蚁s移动，C是上一次迭代向蚂蚁S移动的蚂蚁数量，如果过多，会让这个概率减小，起到防止算法早熟的作用。2. 对每一只蚂蚁 $ant_k$，先运用模式搜索法局部搜索。然后依据概率 $p_{ks}$选择是否向较优蚂蚁s移动，C是上一次迭代向蚂蚁S移动的蚂蚁数量，如果过多，会让这个概率减小，起到防止算法早熟的作用。
   ![image]([pic/pks.png](https://github.com/GodOdimm/ACSwithPatternsearch/blob/main/pic/pks.png)https://github.com/GodOdimm/ACSwithPatternsearch/blob/main/pic/pks.png)
4. 更新信息素
5. k=k+1,若k<最大迭代次数，转2，否则转5
6. 算法结束
