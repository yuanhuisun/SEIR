#include <bits/stdc++.h>

using namespace std;

/**
	 通常，把病毒传染空间范围内的人群分为如下4类：
	S（Susceptible），易感者，指未被传染的健康者但缺乏免疫力，与感染者接触后容易被感染；
	E（Exposed），暴露者，指接触过感染者，但暂时无能力传染他人，适用于潜伏期长的病毒；
	I（Infectious），感染者，指已经确诊的感染人群，可以传染S类人群将其变为E类或者I类群体；
	R（Recovered），康复者，指因隔离或痊愈后有免疫力的人群，包括死亡人数。
	D (Died)，死亡者  (= R * (1-g1))
*/ 
struct SEIR_POINT {
	double Susceptible;
	double Exposed;
	double Infectious;
	double Recovered;    //
	double Died1; //缺乏治疗的死亡者
	double Died2; //有效治疗的死亡者
};

/**
seri0 = [9999, #S,  
			0, #E
			1, #I
			0] #R

seri = seri0

para = [	20, 	# I每天接触的人 r
			0.03, 	# S和E被感染的概率 beta
			0.1, 	# 潜伏者变为感染者的概率 alpha
			0.1, 	# 康复/治愈概率 gamma
			10000] 	# 总人数 N

paraII = [	20, 	# 感染者每天接触的人 r
			20, 	# 潜伏者每天接触的人 r2
			0.03, 	# 感染者传染的概率 beta
			0.03, 	# 潜伏者传染的概率 beta2
			0.1, 	# 潜伏者变为感染者的概率 a
			0.199, 	# 康复/治愈概率 gamma
			0.815,  # 自然康复率 g1
			0.983,  # 治疗康复率 g2
			10000]	# 总人数 N
*/

struct PARA {
	int r;
	int r2;
	double beta;
	double beta2;
	double alpha;
	double gamma;
	double g1; // 
	double g2; // 
	int N;
};

vector <SEIR_POINT> trend, trend2; 

SEIR_POINT SEIR(SEIR_POINT seir,PARA para,int steps)
{
	SEIR_POINT next;
	
	double dS = -(para.beta*para.r*seir.Infectious*seir.Susceptible)/para.N;
	double dE = (para.beta*para.r*seir.Infectious*seir.Susceptible)/para.N - para.alpha*seir.Exposed;
	double dI = para.alpha*seir.Exposed - para.gamma*seir.Infectious;
	double dR = para.gamma*seir.Infectious;
	
	next.Susceptible = seir.Susceptible+dS*steps;
	next.Exposed     = seir.Exposed+dE*steps;
	next.Infectious  = seir.Infectious+dI*steps;
	next.Recovered   = seir.Recovered+dR*steps;
	next.Died1        = next.Recovered * (1-para.g1);
	next.Died2        = next.Recovered * (1-para.g2);


//	cout << "sier (S, E, I, R) = (" << seir.Susceptible << ","<< seir.Exposed << "," << seir.Infectious << "," << seir.Recovered << ")" << endl;
//	cout << "     (dS, dE, dI, dR) = (" << dS << ","<< dE << "," << dI << "," << dR << ")" << endl;
//	cout << "next (S, E, I, R) = (" << next.Susceptible << ","<< next.Exposed << "," << next.Infectious << "," << next.Recovered << ")" << endl;

	return next;
}

SEIR_POINT SEIR_II(SEIR_POINT seir,PARA para,int steps)
{
	SEIR_POINT next;
	
	double dS = -(para.beta*para.r*seir.Infectious*seir.Susceptible)/para.N 
			- (para.beta2*para.r2*seir.Exposed*seir.Susceptible)/para.N;
	double dE = (para.beta*para.r*seir.Infectious*seir.Susceptible)/para.N 
			- para.alpha*seir.Exposed 
			+ (para.beta2*para.r2*seir.Exposed*seir.Susceptible)/para.N;
	double dI = para.alpha*seir.Exposed - para.gamma*seir.Infectious;
	double dR = para.gamma*seir.Infectious;

//	cout << "dS = " << dS <<", dE = " << dE <<", dI = " << dI << ", dR = " << dR << endl;

	next.Susceptible = seir.Susceptible+dS*steps;
	next.Exposed     = seir.Exposed+dE*steps;
	next.Infectious  = seir.Infectious+dI*steps;
	next.Recovered   = seir.Recovered+dR*steps;
	next.Died1        = next.Recovered * (1-para.g1);
	next.Died2        = next.Recovered * (1-para.g2);

	return next;
}

/*
	func:ģ��
	seri:ģ�͵�SERIȺ��
	para:ģ�Ͳ���
	intervene_N����N����и�Ԥ����
	"""
*/
void simulate(SEIR_POINT seir0, PARA para,int steps)
{
	SEIR_POINT seir, seir2, seir_next, seir2_next;
	seir = seir0;
	seir2 = seir0;

	trend.push_back(seir);
	trend2.push_back(seir2);

	for (int i=0;i<=steps;i++)
	{
		seir_next=SEIR(seir,para,1);
		trend.push_back(seir_next);
		seir=seir_next;
		seir2_next=SEIR(seir2,para,1);
		trend2.push_back(seir2_next);
		seir2=seir2_next;
	}

	return;
}

int main()
{
	SEIR_POINT seir_init;
	PARA para;
	freopen("nCoV2019.in","r",stdin);
	freopen("nCov2019.csv","w",stdout);
	ios::sync_with_stdio(false);

	// 读入传染模型初始SEIR
	cin >> seir_init.Susceptible;
	cin >> seir_init.Exposed;
	cin >> seir_init.Infectious;
	cin >> seir_init.Recovered;	
	cin >> seir_init.Died1;
	cin >> seir_init.Died2; 


	// 读入传染模型参数
	cin >> para.r >> para.r2;
	cin >> para.beta >> para.beta2;
	cin >> para.alpha >> para.gamma;
	cin >> para.g1 >> para.g2;
	cin >> para.N;

	cout << "SEIR Polulation: " << endl;
	cout << "# 易感者 S        : " << seir_init.Susceptible << endl;
	cout << "# 暴露者 E        : " << seir_init.Exposed << endl;
	cout << "# 感染者 I        : " << seir_init.Infectious << endl;
	cout << "# 康复者 R        : " << seir_init.Exposed << endl;
	cout << "# 死亡者 D1       : " << seir_init.Died1 << endl;
	cout << "# 死亡者 D2       : " << seir_init.Died2 << endl;

	// output parameters
	cout << "nCOV parameters: " << endl;
	cout << "# 感染者每天接触的人 r        : " << para.r << endl;
	cout << "# 潜伏者每天接触的人 r2       : " << para.r2 << endl;
	cout << "# 感染者传染的概率 beta       : " << para.beta << endl;
	cout << "# 潜伏者传染的概率 beta2      : " << para.beta2 << endl;
	cout << "# 潜伏者变为感染者的概率 alpha : " << para.alpha << endl;
	cout << "# 康复/治愈概率 gamma         :" << para.gamma << endl;
	cout << "# 自然康复率 g1               :" << para.g1 << endl;
	cout << "# 治疗康复率 g2               :" << para.g2 << endl;
	cout << "# 总人数 N                   :" << para.N << endl;

	simulate(seir_init, para, 500);

	cout << "易感者S 暴露者E 感染者I 康复者R 死亡者D1 死亡者D2" << endl;


	for (int i=0;i < trend.size();i++ )
	{
		cout << trend[i].Susceptible << " ";
		cout << trend[i].Exposed << " ";
		cout << trend[i].Infectious << " ";
		cout << trend[i].Recovered << " ";
		cout << trend[i].Died1 << " ";
		cout << trend[i].Died2 << endl;
	}

	cout << "易感者S 暴露者E 感染者I 康复者R 死亡者D1 死亡者D2" << endl;

	for (int i=0;i < trend.size();i++ )
	{
		cout << trend2[i].Susceptible << " ";
		cout << trend2[i].Exposed << " ";
		cout << trend2[i].Infectious << " ";
		cout << trend2[i].Recovered << " ";
		cout << trend2[i].Died1 << " ";
		cout << trend2[i].Died2 << endl;
	}

	fclose(stdin);
	fclose(stdout);
	return 0;
}
 
