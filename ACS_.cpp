#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include<ctime>
#include<vector>
#include<string.h>
#include <fstream>
#include <sstream>

#define D_LEN 4
#define G_NUM 30
#define P 120
#define lb -30
#define rb 30
#define ro 0.9
using namespace std;

unsigned seed = (unsigned)time(0);

double randN(int l, int r) {//产生[l,r)之间的随机数
	if (l > r) swap(l, r);
	double d = (double)rand() / RAND_MAX;
	return (r - l) * d + l;
}

typedef class Vec {//向量类，重载完成加减乘运算
public:
	double a[D_LEN];
public:
	Vec(int l = 0, int r = 0);

	Vec operator+(const Vec& c);

	Vec(double a[D_LEN]);

	Vec(vector<double>& v) {
		for (int i = 0; i < D_LEN; i++)
			this->a[i] = v[i];
	}

	double& operator[](int index) {
		return a[index];
	}

	Vec operator-(const Vec & c);

	Vec& operator=(const Vec & c);

	friend Vec operator*(double num, Vec c);

} Ant;

Vec::Vec(double a[D_LEN]) {
	for (int i = 0; i < D_LEN; i++)
		this->a[i] = a[i];
}

Vec::Vec(int l, int r) {
	for (int i = 0; i < D_LEN; i++)
		a[i] = randN(l, r);
}

class Vec Vec::operator+(const class Vec& c) {
	Vec r(0, 0);
	int i;
	for (i = 0; i < D_LEN; i++)
		r.a[i] = this->a[i] + c.a[i];
	return r;
}

class Vec Vec::operator-(const class Vec& c) {
	Vec r(0, 0);
	int i;
	for (i = 0; i < D_LEN; i++)
		r.a[i] = this->a[i] - c.a[i];
	return r;
}

class Vec& Vec::operator=(const class Vec& c) {
	if (this == &c)
		return *this;
	for (int i = 0; i < D_LEN; i++)
		this->a[i] = c.a[i];
	return *this;
}

Vec operator*(double num, Vec c) {
	for (int i = 0; i < D_LEN; i++)
		c.a[i] *= num;
	return c;
}

double f(Vec & v) {
	int i;
	double res = 0;
	double* x = v.a;
	for (i = 0; i < 3; i++)
	{
		res += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(x[i] - 1, 2);
	}
	return res;
}

Vec e[D_LEN];

void init_e() {//初始化维度，方便模式搜索法
	for (int i = 0; i < D_LEN; i++)
	{
		for (int j = 0; j < D_LEN; j++)
		{
			if (i == j)
				e[i][j] = 1;
			else
				e[i][j] = 0;
		}
	}
}

Vec& mssy(Vec & y, int a, double b, double buchang) {//模式搜索法
	Vec x = y;
	int k;
	int ite_max = 200;
	for (k = 0; k < ite_max; k++)
	{
		int j = 0;
		for (j = 0; j < D_LEN; j++)
		{
			Vec tmp = y + buchang * e[j];
			double n1 = f(tmp);
			double n2 = f(y);
			if (n1 < n2)
			{
				y = tmp;
				continue;
			}
			tmp = y - buchang * e[j];
			n1 = f(tmp);
			if (n1 < n2)
			{
				y = tmp;
				continue;
			}

		}
		if (f(y) < f(x))
		{
			Vec tmp = x;
			x = y;
			y = x + a * (x - tmp);
		}
		else
		{
			if (buchang < 0.0000001)
				break;
			buchang = buchang * b;
			y = x;
		}
	}
	//if (k < ite_max) cout << k << endl;
	return y;

}

class ACS {//蚁群系统
public:
	vector<Ant> x;//蚂蚁
	double f_ave;//平均值
	vector<double> f_value;//蚂蚁函数值
	vector<double> T;//信息素
	int c[P];//向蚂蚁i移动的的蚂蚁数量
	double pro[P][P];
public:
	ACS();//构造函数

	double get_fave();//得到f(x)的平均值

	void cacul_f();//计算所有f(x)

	void cacul_pro();//计算pro[P][P]，即蚂蚁移动到其他更优蚂蚁的概率

	int best_index();//f值最小的蚂蚁下标

	int Kernel();//

	int move_toward(int i);//轮盘赌选择蚂蚁i向谁移动，返回下标

	void update_T();//更新信息素，最优解增加，其他减少

	//Vec &mssy(Vec &y, int a, double b, double buchang);

};

int ACS::best_index() {
	int best = 0;
	for (int i = 1; i < P; i++)
	{
		if (f_value[i] < f_value[best])
		{
			best = i;
		}
	}
	return best;
}

ACS::ACS() {
	f_value = vector<double>(P, 0);
	T = vector<double>(P, 6.0);
	memset(this->c, 0, sizeof(this->c));
	int i, j;
	for (i = 0; i < P; i++)
	{
		x.push_back(Ant(lb, rb));
	}
	f_ave = get_fave();
}

double ACS::get_fave() {
	double res = 0.0;
	int i;
	for (i = 0; i < f_value.size(); i++)
		res += f_value[i];
	this->f_ave = res / f_value.size();
	return this->f_ave;
}

int ACS::Kernel() {
	int k = 0;
	double a = 0.08;//蚂蚁全局搜索
	for (k = 0; k < 10; k++)
	{
		int i, j;
		for (i = 0; i < P; i++)
			mssy(x[i], 1, 0.2, 0.25);

		cacul_f();

		cacul_pro();

		//每只蚂蚁向谁移动
		int best = best_index();
		for (i = 0; i < P; i++)
		{
			if (i == best)
				continue;
			int s = move_toward(i);
			c[s]++;
			x[i] = x[i] + a * (x[s] - x[i]);
		}

		//更新信息素
		update_T();


		memset(this->c, 0, sizeof(this->c));

	}
	return k;
}

void ACS::update_T() {
	int i, s;
	double t0 = 6.0;
	double ave = get_fave();

	int best = 0;
	for (i = 1; i < P; i++)
	{
		if (f_value[i] < f_value[best])
			best = i;
	}

	for (s = 0; s < P; s++)
	{
		if (s == best)
		{
			T[s] += t0 * (1.0 * c[s] / P) + t0;
		}
		else
			T[s] = ro * T[s];
		//T[s] = ro * T[s] + t0 * (1.0 * c[s] / P) + t0 * exp((f_value[s] - f_value[best]) / (ave - f_value[best]));
	}
}

int ACS::move_toward(int i) {
	int j;
	double r1 = randN(0, 1);
	for (j = 0; j < P; j++)
	{
		r1 -= pro[i][j];
		if (r1 <= 0)
			return j;
	}
	return P - 1;
}

void ACS::cacul_f() {
	int i;
	for (i = 0; i < x.size(); i++)
		this->f_value[i] = f(x[i]);
}

void ACS::cacul_pro() {
	int i, s;

	double ave = get_fave();

	vector<double> all(P, 0.0);
	for (i = 0; i < f_value.size(); i++)
	{
		double n = (f_value[i] - ave) / 10;
		for (s = 0; s < f_value.size(); s++)
		{
			double ds = f_value[i] - f_value[s];
			if (ds > 0)// i向小的s转移
			{
				this->pro[i][s] = T[s] * exp(-n / ds) * exp(-1.0 * c[s] / P);
				all[i] += this->pro[i][s];
			}
			else
				this->pro[i][s] = 0.0;
		}
	}
	for (i = 0; i < f_value.size(); i++)
	{
		if (all[i] == 0) continue;
		for (s = 0; s < f_value.size(); s++)
		{
			pro[i][s] /= all[i];
		}
	}

}

int main() {
    srand(seed);
    init_e();

    int i;
    ACS acs;
    int k = acs.Kernel();
    //cout << "k=" << k << endl;
    int best = acs.best_index();
    Vec &v = acs.x[best];
    cout << setprecision(10);
    for (i = 0; i < D_LEN; i++) cout << v[i] << ' ';
    cout << endl << f(v) << endl;

    getchar();
    return 0;

}
