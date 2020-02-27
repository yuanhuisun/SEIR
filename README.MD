# n-Cov 2019 Virus Infectiion Simulation (C++)

## SIER模型参考

基于Harrison Ho在GIT上《2019-nCoV感染动力学分析及防护的必要性》[https://github.com/Harrison-Ho/Analysis-of-2019-nCoV-propagation-based-on-Model/tree/master/markdown],
我用C++实现的nCov19的病毒传染模拟。

## SEIR Model is used in the simulation

通常，把病毒传染空间范围内的人群分为如下4类：  
	S（Susceptible），易感者，指未被传染的健康者但缺乏免疫力，与感染者接触后容易被感染；  
	E（Exposed），暴露者，指接触过感染者，但暂时无能力传染他人，适用于潜伏期长的病毒；  
	I（Infectious），感染者，指已经确诊的感染人群，可以传染S类人群将其变为E类或者I类群体；  
	R（Recovered），康复者，指因隔离或痊愈后有免疫力的人群，包括死亡人数。  
	D (Died)，死亡者  (= R * (1-g1)) ，g1为康复者的康复率，1-g1为康复者的死亡率，参数g1对SEIR模型没有影响 
	
易感者S一开始先成为潜伏者E经历一段潜伏期，一段时间之后才会成为感染者I，感染者I被发现后经治疗或依靠自身免疫成为康复者R。

假设一个地区总人数为N，I个感染者如不隔离，每人每天会碰到r个人，以概率$\beta$传染这r个人；潜伏者E会有概率$\alpha$成为感染者；同时感染者因治疗或免疫有$\gamma$的概率成为康复者R，可以建立如下方程组:
$$
\left\{ 
\begin{array}{c}
\frac{dS}{dt}=-\frac{r \beta I S}{N} \\
\frac{dE}{dt}=\frac{r \beta I S}{N}-aE \\
\frac{dI}{dt}=aE-\gamma I \\
\frac{dR}{dt}=\gamma I \\
N = S+E+I+R
\end{array}
\right.  
$$
上述方程组没有解析解，但是可以利用计算机软件得到其数值解。	
 
