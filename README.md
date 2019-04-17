# PTA-代码以及部分STL函数的运用
//lower_bound函数的运用
/*
#include <iostream>
#include <algorithm>
#include <cstdlib>
using namespace std;

int main()
{
	int a[122], n;
	cin >> n;
	for (int i = 0; i < n; i++)
		cin >> a[i];
	sort(a, a + n);
	int m;
	cin >> m;
	//lower_bound函数如果找到则返回该数的地址，否则返回end。所以用返回值减去查找的首地址就可以得到该数在数组中的位置
	int t = lower_bound(a, a + n, m) - a;
	if (a[t] == m)
		cout << "Yes" << endl;
	else
		cout << "No" << endl;
	system("pause");
	return 0;
}
*/


//结构体运用
/*
#include <iostream>

using namespace std;

//struct Point {
//	int x, y;
//	Point(int x = 0, int y = 0) : x(x), y(y) {}
//};
struct Point {
	int x, y;
	Point(int x = 0, int y = 0)
	{
		this->x = x;
		this->y = y;
	}
};

Point operator + (const Point& A, const Point& B) {
	return Point(A.x + B.x, A.y + B.y);
}

ostream& operator << (ostream &out, const Point& p) {
	out << "(" << p.x << "," << p.y << ")";
	return out;
}

int main()
{
	Point a, b(1, 2);
	a.x = 3;
	cout << a + b << endl;
	system("pause");
	return 0;
}
*/



//string函数简单运用
/*
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <sstream>
using namespace std;

int main()
{
	string line;
	int cut  = 1;
	while (getline(cin, line))
	{
		int sum = 0, x;
		stringstream ss(line);
		while (ss >> x)
		{
			sum += x;
			cout << cut << ' ' << x << endl;
			cut++;
		}
		cout << sum << endl;
	}
	system("pause");
	return 0;
}
*/



//救济金发放
/*
#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

const int maxn = 25;

int n, k, m, a[maxn];

int go(int p, int d, int t)
{
	while (t--){
		do {
			p = (p + d + n - 1) % n + 1;
		} while (a[p] == 0);
	}
	return p;
}

int main()
{
	while (scanf_s("%d %d %d", &n, &k, &m) == 3 && n)
	{
		for (int i = 1; i <= n; i++)	a[i] = i;
		int left = n;
		int p1 = n, p2 = 1;
		while (left)
		{
			p1 = go(p1, 1, k);
			p2 = go(p2, -1, m);
			printf("%3d", p1); left--;
			if (p2 != p1)
			{
				printf("%3d", p2);
				left--;
			}
			a[p1] = a[p2] = 0;
			if (left)	printf(",");
		}
		printf("\n");
	}
	system("pause");
	return 0;
}
*/



//数位DP
/*
#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

typedef long long ll;

int main()
{
	char num1[20],num2[20];
	cin >> num1;
	ll sum = 0;
	int n = strlen(num1);
	int digit;
	if (n < 4)
		cout << "No" << endl;
	else
	{
		//digit = n * (n - 1) * (n - 2) * (n - 3) / 12;
		for (int i = 0; i <= n - 4; i++)//2前的数
		{
			digit = (n - i) * (n - i - 1) * (n - i - 2) * (n - i - 3) / 12;
			sum += digit * pow(10, (n - i - 4));
		}
	}
	ll T,T1;
	T = pow(10, n + 1);
	T1 = 0;
	for (int i = 0; i < n; i++)
		T1 += (num1[i] - 48) * pow(10, n - i - 1);
	for (int i = T1; i < T; i++)
	{
		for(int i)
	}
	cout << sum << endl;
	system("pause");
	return 0;
}
*/



//归并排序
/*
#include <iostream>
#include <cstdlib>

using namespace std;

const int N = 1e5 + 9;
int a[N], t[N];
int i, j, k;

void merge(int l, int r)//合并
{
	int mid = (l + r) >> 1;
	i = l; j = mid + 1;
	k = 0;
	while (i <= mid && j <= r)
	{
		if (a[i] < a[j])
			t[k++] = a[i++];
		else
			t[k++] = a[j++];
	}
	while (i <= mid)
		t[k++] = a[i++];
	while (j <= r)
		t[k++] = a[j++];
	for (i = 0; i < k; i++)//将t数组序列存入原数组
		a[l++] = t[i];
}
void merge_s(int l, int r)//分割
{
	if (l == r) return;
	int mid = (l + r) >> 1;
	merge_s(l, mid);
	merge_s(mid + 1, r);
	merge(l, r);
}
int main()
{
	int n;
	cin >> n;
	for (i = 0; i < n; i++)
	{
		cin >> a[i];
	}
	merge_s(0, n - 1);
	for (i = 0; i < n - 1; i++)
	{
		cout << a[i] << " ";
	}
	cout << a[n - 1] << endl;
	system("pause");
	return 0;
}
*/


//实验二
/*
#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

class Person
{
public:
	Person();
	Person(Person &p);
	void Register(string n, int y, char g);
	void showme();
	~Person()
	{
		cout << "Now destroying the instance of Person" << endl;
	}
private:
	string name;
	int year;
	char gender;
};
//构造函数
Person::Person()
{
	name = "peter";
	year = 18;
	gender = 'm';
}
//复制构造函数
Person::Person(Person &p)
{
	name = p.name;
	year = p.year;
	gender = p.gender;
}
//赋值函数
void Person::Register(string n, int y, char g)
{
	name = n;
	gender = g;
	year = y;
}
//输出函数
inline void Person::showme()
{
	cout << name << ' ' << year << ' ' << gender << endl;
}
int main()
{
	Person *p1, *p2;
	Person G, M;
	int y;
	string n;
	char g;
	p1 = &G;
	p2 = &M;
	p1->showme();
	p2->showme();
	cin >> n >> y >> g;
	p1->Register(n, y, g);
	p2->Register(n, y, g);
	(*p1).showme();
	(*p2).showme();
	system("pause");
	return 0;
}
*/



//实验一
/*
#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

class AccountType
{
public:
	AccountType(double x,double y,string n);
	void Deposit(double x);
	void Withdraw(double x);
	void WriteAccount();
	void Writename();
private:
	double Id, balance;
	string name;
};

//构造函数
AccountType::AccountType(double x, double y, string n)
{
	Id = x;
	balance = y;
	name = n;
}

//存钱函数
inline void AccountType::Deposit(double x)
{
	balance += x;
}

//取钱函数
inline void AccountType::Withdraw(double x)
{
	balance -= x;
}

//打印账户余额信息函数
inline void  AccountType::WriteAccount()
{
	cout << balance << endl;
}

//打印客户信息函数
inline void AccountType::Writename()
{
	cout << Id << ' ' << name << ' ' << balance << endl;
}

int main()
{
	double a, b;
	string c;
	cin >> a >> c >> b;
	AccountType My(a, b, c);
	double d, e;
	cin >> d;
	My.Deposit(d);
	cout << "存入钱后的账户余额：";
	My.WriteAccount();
	cin >> e;
	My.Withdraw(e);
	cout << "取出钱后账户余额为：";
	My.WriteAccount();
	My.Writename();
	system("pause");
	return 0;
}
*/



//链式前向星
/*
#include <iostream>
#include <cstdlib>

using namespace std;

#define MAXN 100501

struct NODE {
	int w;//边权值
	int e;//边终点
	int next; //next[i]表示与第i条边同起点的下一条边的储存位置
}edge[MAXN];

int cnt;
int head[MAXN];

void add(int u, int v, int w) {
	edge[cnt].w = w;
	edge[cnt].e = v;    //edge[i]表示第i条边的终点 
	edge[cnt].next = head[u]; //head[i]表示以i为起点的最后一条边的储存位置 
	head[u] = cnt++;
}
int main() {
	memset(head, 0, sizeof(head));
	cnt = 1;
	int n;
	cin >> n;
	int a, b, c;//a为起点，b为终点，c为权值
	while (n--) {
		cin >> a >> b >> c;
		add(a, b, c);
	}
	int start;
	cin >> start;
	for (int i = head[start]; i != 0; i = edge[i].next)
		cout << start << "->" << edge[i].e << " " << edge[i].w << endl;
	system("pause");
	return 0;
}
*/



//高精度算法
/*
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
using namespace std;

const int N = 1e5 + 9;

char poly[N], R[N], sum[N] = {0};

char symbol;

int numl, numr, cut, flag;

void Add(char x[], char y[]);
//void Redu(char x[], char y[]);

int main()
{
	cin >> poly;
	cut = 0;
	int n = strlen(poly);
	for (int i = 0; i < n; i++)
	{
		cut++;
		if (poly[i] == '+' || poly[i] == '-' || poly[i] == '=')
		{
			symbol = poly[i];
			break;
		}
		sum[i + 1] = poly[i];
	}
	for (int i = cut,j = 0; i < n - 1; i++,j++)
	{
		if (poly[i] == '+')
		{
			if (symbol = '+')
				Add(sum, R);
			//else
			//	Redu(sum, R);
			symbol = poly[i];
			memset(R, 0, sizeof(R));
		}
		else if (poly[i] == '-')
		{
			if (symbol = '+')
				Add(sum, R);
			//else
				//Redu(sum, R);
			symbol = poly[i];
			memset(R, 0, sizeof(R));
		}
		else
			R[j] = poly[i];
	}
	for (int i = 0; i < strlen(sum); i++)
		cout << sum[i];
	cout << endl;
	system("pause");
	return 0;
}

//加法
void Add(char x[], char y[])
{
	int n, m, t;
	n = strlen(x);
	m = strlen(y);
	if (n > m)
	{
		for (int i = m - 1, j = n; i > 0; i--, j--)
		{
			t = x[j] + y[i] - 96;
			if (t > 9)
			{
				x[j] = '9';
				x[j - 1] += 1;
			}
			else
				x[j] = t + 48;
		}
		if (sum[0] != '0')
			for (int i = 0; i <= n; i++)
				sum[i + 1] = sum[i];
	}
	else
	{
		for (int i = m - 1, j = n; i > 0; i--, j--)
		{
			t = x[j] + y[i] - 96;
			if (t > 9)
			{
				x[j] = '9';
				x[j - 1] += 1;
			}
			else
				x[j] = t + 48;
		}
		if (sum[0] != '0')
			for (int i = 0; i <= n; i++)
				sum[i + 1] = sum[i];
	}
}

//减法
void Redu(char x[], char y[])
{

}
*/



//多项式之积与和
/*
#include <string>
#include <iostream>
#include <cstdlib>
#include <algorithm>
using namespace std;

const int N = 1e5 + 9;

int Add[N] = {0},
prod[N] = {0}, 
poly1[N] = {0}, 
poly2[N] = {0};

int main()
{
	int n, m, a, b, max1, max2, flag;
	flag = 0;
	max1 = max2 = -1;
	cin >> n;
	while (n--)
	{
		cin >> a >> b;
		Add[b] = Add[b] + a;
		poly1[b] = poly1[b] + a;
		if (max1 < b)
			max1 = b;
	}
	cin >> m;
	while (m--)
	{
		cin >> a >> b;
		Add[b] = Add[b] + a;
		poly2[b] = poly2[b] + a;
		if (max2 < b)
			max2 = b;
	}
	for (int i = 0; i <= max1; i++)
		for (int j = 0; j <= max2; j++)
			if (poly1[i] != 0 && poly2[j] != 0)
				prod[i + j] = prod[i + j] + poly1[i] * poly2[j];
	for (int i = max1 + max2; i >= 0; i--)
		if (prod[i] != 0)
		{
			cout << prod[i] << ' ' << i << ' ';
			flag = 1;
		}
	if (flag == 1)
		cout << "\b" << endl;
	else
		cout << 0 << ' ' << 0 << endl;
	flag = 0;
	for (int i = max(max1, max2); i >= 0; i--)
		if (Add[i] != 0)
		{
			cout << Add[i] << ' ' << i << ' ';
			flag = 1;
		}
	if (flag == 1)
		cout << "\b" << endl;
	else
		cout << 0 << ' ' << 0 << endl;
	system("pause");
	return 0;
}
*/


//多项式之和
/*
#include <iostream>
#include <cstdlib>
using namespace std;

const int N = 1e4 + 9;

int num[N] = {0};

int maxl(int x, int y);
int main()
{
	int n, m, a, b,maxl = -1;
	cin >> n;
	for(int i = 0;i < n;i++)
	{
		cin >> a >> b;
		if (maxl < b)
			maxl = b;
		num[b] = num[b] + a;
	}
	cin >> m;
	for(int i = 0;i < m;i++)
	{
		cin >> a >> b;
		if (maxl < b)
			maxl = b;
		num[b] = num[b] + a;
	}
	for (int i = maxl; i >= 0; i--)
		if (num[i] != 0)
			cout << num[i] << ' ' << i << ' ';
	cout << "\b" << endl;
	system("pause");
	return 0;
}
*/


//PTA——01-复杂度2 Maximum Subsequence Sum （25 分)
/*
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

const int N = 1e5 + 9;

int num[N];

int main()
{
	int n, L, R, cut, flag;
	long long Thissum, Maxsum;
		Thissum = Maxsum = cut = flag = 0;
		cin >> n;
		for (int i = 0; i < n; i++)
		{
			cin >> num[i];
			if (num[i] >= 0)
				flag = 1;
		}
		if (flag == 1)
		{
			for (int i = 0; i < n; i++)
			{
				Thissum += num[i];
				cut++;
				if (Thissum > Maxsum)
				{
					Maxsum = Thissum;
					R = i;
					L = i - cut + 1;
				}
				else if (Thissum < 0)
				{
					Thissum = 0;
					cut = 0;
				}
			}
			if(Maxsum != 0)
				cout << Maxsum << ' ' << num[L] << ' ' << num[R] << endl;
			else 
				cout << 0 << ' ' << 0 << ' ' << 0 << endl;
		}
		else
			cout << 0 << ' ' << num[0] << ' ' << num[n - 1] << endl;
		system("pause");
		return 0;
}
*/



//最大连续子列和问题——多组输入
/*
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

const int N = 1e5 + 9;

int num[N];

int main()
{
	int n,L,R,cut,flag;
	long long Thissum, Maxsum;
	while (scanf_s("%d", &n) != EOF)
	{
		Thissum = Maxsum = cut = flag = 0;
		for (int i = 0; i < n; i++)
		{
			cin >> num[i];
			if (num[i] >= 0)
				flag = 1;
		}
		if (flag == 1)
		{
			for (int i = 0; i < n; i++)
			{
				Thissum += num[i];
				cut++;
				if (Thissum >= Maxsum)
				{
					Maxsum = Thissum;
					R = i;
					L = i - cut + 1;
				}
			    else if (Thissum < 0 || num[i + 1] - num[i] != 1)
				{
					Thissum = 0;
					cut = 0;
				}
			}
			if (Maxsum != 0)
				cout << Maxsum << ' ' << num[L] << ' ' << num[R] << endl;
			else
				cout << 0 << ' ' << 0 << ' ' << 0 << endl;
		}
		else
			cout << 0 << ' ' << num[0] << ' ' << num[n - 1] << endl;
	}
	return 0;
}
*/



//最大子列和问题
/*
#include <iostream>
#include <cstdlib>
using namespace std;

const int N = 1e5 + 9;

int num[N];

int main()
{
	int k;
	long long Thissum, Maxsum;
	Thissum = Maxsum = 0;
	cin >> k;
	for (int i = 0; i < k; i++)
	{
		cin >> num[i];
		Thissum += num[i];
		if (Thissum > Maxsum)
			Maxsum = Thissum;
		else if (Thissum < 0)
			Thissum = 0;
	}
	cout << Maxsum << endl;
	system("pause");
	return 0;
}
*/



//纸对折问题
/*
#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

int main()
{
	int n;
	cin >> n;
	for (long long i = 1; i < pow(2, n) - 1; i++)
	{
		if (i % 3 == 0)
			cout << 'U';
		else
			cout << 'D';
	}
	if (n != 1)
		cout << 'U';
	else
		cout << 'D';
	cout << endl;
	cout << pow(2, 25) << endl;
	system("pause");
	return 0;
}
*/



//这个糖不甜
/*
#include <iostream>
#include <cstdlib>

using namespace std;

int num[10009];

int main()
{
	int n,cut = 0,t = -1,a,b,s = 0,flag = 0;
	cin >> n;
	for (int i = 0; i < n; i++)
		cin >> num[i];
	for (int i = 0; i < n; i++)
	{
		if (num[i] >= 0)
			flag = 1;
	}
	if (flag == 0)
		cout << 0 << ' ' << num[0] << ' ' << num[n - 1] << endl;
	else
	{
		for (int i = 0; i < n; i++)
		{
			cut += num[i];
			s++;
			if (num[i] < 0)
			{
				if (t < cut)
				{
					t = cut;
					b = i - 1;
					a = i - s;
				}
				cut = 0;
				s = 0;
			}
			if (num[i + 1] - num[i] != 1)
			{
				if (t < cut)
				{
					t = cut;
					b = i;
					a = i - s + 1;
				}
				cut = 0;
				s = 0;
			}
		}
		cout << t << ' ' << num[a] << ' ' << num[b] << endl;
	}
	system("pause");
	return 0;
}
*/



//迷雾森林
/*
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
using namespace std;

const int N = 1e3 + 9;

string str;

int main()
{
	int n, m, x, y, flag = 1;
	cin >> n >> m >> x >> y;;
	for (int i = 0; i < n; i++)
	{
		getchar();
		getline(cin, str);
			for (int j = 0; j < str.size(); j++)
			{
				if (str[j] == 'E')
				{
					y = y + 1;
					//cout << 1 << endl;
					if (y > m)
						flag = 0;
					else
						flag = 1;
				}
				else if (str[j] == 'W')
				{
					y = y - 1;
					//cout << 2 << endl;
					if (y < 0)
						flag = 0;
					else
						flag = 1;
				}
				else if (str[j] == 'N')
				{
					x = x - 1;
					//cout << 3 << endl;
					if (x < 0)
						flag = 0;
					else
						flag = 1;
				}
				else
				{
					x = x + 1;
					//cout << 4 << endl;
					if (x > n)
						flag = 0;
					else
						flag = 1;
				}
			}
	}
	if (flag == 1)
		cout << "Never!" << endl;
	else
		cout << x << ' ' << y << endl;
	system("pause");
	return 0;
}
*/



//ACM招新了
/*
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
using namespace std;

const int N = 1e5 + 9;

int high[N];

int minx(int x, int y);
int maxx(int x, int y);

int main()
{
	int n, H, m, a, b;
	while (scanf_s("%d %d %d", &n, &H, &m) == 3)
	{
		for (int i = 0; i <= n; i++)
			high[i] = H;
		while (m--)
		{
			cin >> a >> b;
			for (int i = minx(a, b) + 1; i < maxx(a, b); i++)
			{
				if (high[i] >= minx(high[a], high[b]))
				{
					high[i] = minx(high[a], high[b]) - 1;
				}
			}
		}
		for (int i = 1; i <= n; i++)
			cout << high[i] << endl;
	}
	system("pause");
	return 0;
}

int minx(int x, int y)
{
	if (x > y)
		return y;
	return x;
}

int maxx(int x, int y)
{
	if (x > y)
		return x;
	return y;
}
*/



//数据结构——二分查找
/*
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

#define MAXSIZE 10
#define NotFound 0
typedef int ElementType;

typedef int Position;
typedef struct LNode *List;
struct LNode {
	ElementType Data[MAXSIZE];
	Position Last; // 保存线性表中最后一个元素的位置 
};

List ReadInput(); // 裁判实现，细节不表。元素从下标1开始存储
Position BinarySearch(List L, ElementType X);

int main()
{
	List L;
	ElementType X,N;
	Position P;
	struct LNode A;

	L = &A;
	cin >> N;
	for (int i = 1; i <= N; i++)
	{
		cin >> L->Data[i];
		//cout << L->Data[i] << endl;
		if (i == N)
			L->Last = N;
	}
	cin >> X;
	P = BinarySearch(L, X);
	cout << P << endl;
	system("pause");
	return 0;
}

Position BinarySearch(List L, ElementType X)
{
	Position l, r, mid;
	l = 1;	r = L->Last;
	//cout << l << ' ' << r << endl;
	while (l <= r)
	{
		mid = (l + r) / 2;
		if (L->Data[mid] > X)
			r = mid - 1;
		else if (L->Data[mid] < X)
			l = mid + 1;
		else
			return mid;
	}
	return NotFound;
}
*/



//币值转换
/*
#include <iostream>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <cstdlib>
using namespace std;

int main()
{
	char num[10];
	char number[10] = { 'a','b','c','d','e','f','g','h','i','j' };
	char numb[10][5] = {
		{},
		{"S"},
		{"B"},
		{"Q"},
		{"W"},
		{"S"},
		{"B"},
		{"Q"},
		{"Y"}
	};
	cin >> num;
	int n = strlen(num);
	for (int i = 0,j = n - 1; i < n; i++,j--)
	{
		if (num[i] != '0')
			cout << number[num[i] - 48] << numb[j];
		else if(num[i] == '0' && num[i - 1] != '0' && num[i + 1] != '0')
			cout << number[num[i] - 48];
	}
	cout << endl;
	system("pause");
	return 0;
}
*/



//sort()结构体排序应用
/*
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdlib>
using namespace std;

struct Man {
	char name[15], num[15];
	int fen;
}a[100009];

bool cmp(Man a, Man b)
{
	return a.fen < b.fen;
}


int main()
{
	int n;
	cin >> n;
	for (int i = 0; i < n; i++)
		cin >> a[i].name >> a[i].num >> a[i].fen;
	sort(a, a + n, cmp);
	cout << a[n - 1].name << ' ' << a[n - 1].num << endl;
	cout << a[0].name << ' ' << a[0].num << endl;
	system("pause");
	return 0;
}
*/



//sort()的用法 结构体数组，递增，递减
/*
#include <iostream>
#include <algorithm>
#include <cstdlib>
using namespace std;

struct node {
	int x, y;
}a[10];

bool cmpj(node a, node b)
{
	return a.x > b.x;
}

bool cmp(int a, int b)
{
	return a < b;
}

int main()
{
	for (int i = 0; i < 10; i++)
		a[i].x = i;
	int b[10] = {9,8,7,6,5,4,3,2,1,0};
	for (int i = 0; i < 10; i++)
		cout << b[i] << endl;
	sort(b, b + 10, cmp);
	for (int i = 0; i < 10; i++)
		cout << b[i] << endl;
	for (int i = 0; i < 10; i++)
		cout << a[i].x << endl;
	sort(a, a + 10, cmpj);
	for (int i = 0; i < 10; i++)
		cout << a[i].x << endl;
	system("pause");
	return 0;
}
*/



//string的用法 
/*
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

int main()
{
	string str1, str2, str3;
	//getline(cin,str3);
	cin >> str1;
	cin >> str2;
	getchar();
	getline(cin, str3);
	cout << str1 << ' ' << str2 << endl;
	//str1 += str2;
	cout << str3 + str1 + str2 << endl;
	str1 = str2;
	cout << str1 << endl;
	cout << str1.size() << endl;
	cout << str3 << endl;
	system("pause");
	return 0;
}
*/



//PAT违禁品
/*
#include <iostream>
#include <vector>
#include <map>

using namespace std;

int main() {
	int n, k, t1, t2;
	map<int, vector<int>> m;
	cin >> n >> k;
	for (int i = 0; i < n; i++) {
		cin >> t1 >> t2;
		m[t1].push_back(t2);
		m[t2].push_back(t1);
	}
	while (k--) {
		int cnt, flag = 0, a[100000] = { 0 };
		cin >> cnt;
		vector<int> v(cnt);
		for (int i = 0; i < cnt; i++) {
			cin >> v[i];
			a[v[i]] = 1;
		}
		for (int i = 0; i < v.size(); i++)
			for (int j = 0; j < m[v[i]].size(); j++)
				if (a[m[v[i]][j]] == 1) flag = 1;
		printf("%s\n", flag ? "No" : "Yes");
	}
	return 0;
}
*/



//月饼销量
/*
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
using namespace std;

const int N = 1e3 + 9;

long moon[N] = {0};

int main()
{
	int n, m, t, maxl = -1, k, cut = 0;
	cin >> n >> m;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			cin >> t;
			moon[j] = moon[j] + t;
			if (maxl < moon[j])
			{
				maxl = moon[j];
				k = j;
			}
		}
	cout << moon[k] << endl;
	for (int i = 0; i < n; i++)
	{
		if (moon[i] == moon[k] && cut == 0)
		{
			cout << i + 1;
			cut++;
		}
		else if(moon[i] == moon[k] && cut != 0)
			cout << ' ' << i + 1;
	}
	cout << endl;
	system("pause");
	return 0;
}
*/



//牛客--图书管理员
/*
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
using namespace std;

const int N = 1e3 + 9;

long long num[N];

typedef struct Man{
	int L;
	long long x_num;
}man;

man m[N];

int main()
{
	int n, q, flag;
	cin >> n >> q;//书 读者
	for (int i = 0; i < n; i++)
		cin >> num[i];
	sort(num, num + n);
	for (int i = 0; i < q; i++)
		cin >> m[i].L >> m[i].x_num;
	for (int i = 0; i < q; i++)
	{
		flag = 0;
		for (int j = 0; j < n; j++)
		{
			if (num[j] % (long long)pow(10, m[i].L) == m[i].x_num)
			{
				cout << num[j] << endl;
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			cout << "-1" << endl;
	}
	system("pause");
	return 0;
}
*/



//queue函数具体实现以及运用
/*
#include <iostream>
#include <queue>
#include <cstdlib>
#include <string>

using namespace std;
//push：入队插入队尾  
//pop:队首元素输出 
//size ： 返回队列中元素的个数  
//back:返回队列中的最后一个元素
//front :返回队列中的第一个元素
//empty :判断队列是否为空
int main()
{
	queue<int>q;
	for (int i = 0; i < 5; i++)
	{
		q.push(i);
		cout << q.front() << endl;
		q.pop();
	}
	system("pause");
	return 0;
}
*/



//中序遍历和后序遍历构建二叉树，并且层序遍历输出
/*
#include<iostream>
#include<queue>
using namespace std;
//结点权值作为结点编号
int postOrder[31];     //后序遍历结点
int inOrder[31];       //中序遍历结点
int leftNodes[31];              //保存某结点的左子树编号
int rightNodes[31];           //保存某结点的右子树编号

							  //根据inOrder[L1]到inOrder[R1]  和postOrder[L1]到postOrder[R1]的结点编号 来构建树
							  //返回根节点
int buildTree(int L1, int R1, int L2, int R2) 
{
	if (R1 < L1)  //空树
		return -1;
	int root = postOrder[R2];      //后序遍历序列最后一个结点一定是根结点
	int p = 0;
	while (inOrder[p] != root)    //找到中序遍历序列中对应哪个根结点的结点
		p++;
	int count = p - L1;          //左子树结点总数

								 //p是中序序列的根，从L1到p-1为左子树，对应的后续序列的从L2到L2+count-1
	leftNodes[root] = buildTree(L1, p - 1, L2, L2 + count - 1);
	//中序序列从p+1到R1为右子树，对应的后续序列从L2+count到R2 - 1 ！！因为根节点已经去掉了！！
	rightNodes[root] = buildTree(p + 1, R1, L2 + count, R2 - 1);
	return root;
}

//层序遍历
//传了个N进去是因为输出格式控制 = = 
void printVex(int root, int N) 
{
	queue<int> q;
	q.push(root);
	while (q.size()) 
	{
		int vex = q.front();
		if (N == 1)
			cout << vex;
		else
			cout << vex << " ";
		q.pop();
		if (leftNodes[vex] != -1)
			q.push(leftNodes[vex]);
		if (rightNodes[vex] != -1)
			q.push(rightNodes[vex]);
		N--;
	}
}
int main()
{
	int N;
	cin >> N;
	for (int i = 0; i < N; i++)
		cin >> postOrder[i];
	for (int j = 0; j < N; j++) 
		cin >> inOrder[j];
	int root = buildTree(0, N - 1, 0, N - 1);
	printVex(root, N);
	return 0;
}
*/



//有几个PAT
/*
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstdlib>
using namespace std;

const int N = 1e5 + 9;

char ch[N];

int num[5];

int main()
{
	cin >> ch;
	int n = strlen(ch);
	memset(num, 0, 5);
	for (int i = 0; i < n; i++)
	{
		if (ch[i] == 'P')
			num[0]++;
		else if (ch[i] == 'A')
			num[1] = (num[0] + num[1]) % 1000000007;
		else
			num[2] = (num[1] + num[2]) % 1000000007;
	}
	cout << num[2] << endl;
	system("pause");
	return 0;
}
*/



//二叉树
/*
#include <iostream>
#include <cstring>
#include <queue>
#include <algorithm>
using namespace std;

typedef struct TreeNode* BinTree;

const int N = 50;

typedef struct TreeNode {
	int data;
	BinTree L;
	BinTree R;
};

int pre[N], in[N], post[N];

int main()
{
	int n;
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
	{
		scanf("%d", &post[i]);
	}
	for (int i = 0; i < n; i++)
	{
		scanf("%d", &in[i]);
	}
	BinTree root = create(0, n - 1, 0, n - 1);
	bfs(root);
	return 0;
	return 0;
}

BinTree create(int postl, int postr, int inl, int inr) {
	if (postl > postr) {
		return NULL;
	}
	BinTree root = new TreeNode;
	root->data = post[postr];
	int k;
	for (k = inl; k <= inr; k++) {
		if (in[k] == post[postr]) {
			break;
		}
	}
	int numLeft = k - inl;
	root->L = create(postl, postl + numLeft - 1, inl, k - 1);
	root->R = create(postl + numLeft, postr - 1, k + 1, inr);
	return root;
}
int num = 0;

void bfs(BinTree root) {
	queue < BinTree > q;
	q.push(root);
	while (!q.empty()) {
		BinTree now = q.front();
		q.pop();
		printf("%d", now->data);
		num++;
		if (num < n) printf(" ");
		if (now->L != NULL) q.push(now->L);
		if(now->R!= NULL)q.push(now->R)
	}

}
*/



//进制转换2 - 10
/*
#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;

const int N = 1e5 +9;

char ch[N];

int main()
{
cin >> ch;
int t = strlen(ch);
int sum = ch[0] - 48;
for (int i = 0; i < t - 1; i++)
{
sum = (ch[i + 1] - 48) + 2 * sum;
}
cout << sum << endl;
system("pause");
return 0;
}
*/



//牛客五星——借教室
/*
#include <iostream>
#include <cstdio>
#include <cstdlib>
using namespace std;

typedef long long ll;

const int N = 1e6 + 9;

typedef struct Day {
	ll Num;//数量
	ll time_s;//开始时间
	ll time_e;//结束时间
}Room;

ll room[N];
Room room_j[N];

int main()
{
	int n, m, cut, flag = 1;
	cin >> n >> m;
	for (ll  i = 1; i <= n; i++)
		cin >> room[i];
	for (ll i = 1; i <= m; i++)
	{
		cin >> room_j[i].Num >> room_j[i].time_e >> room_j[i].time_s;
		if(flag != 0)
		for (ll j = room_j[i].time_e; j <= room_j[i].time_s; j++)
		{
			if (room_j[i].Num <= room[j])
			{
				flag = 1;
				room[j] = room[j] - room_j[i].Num;
				//cout << room[j] << ' ';
			}
			else
			{
				cut = i;
				flag = 0;
				break;
			}
		}
	}
	if (flag == 0)
		cout << "-1" << "\n" << cut << endl;
	else
		cout << '0' << endl;
	system("pause");
	return 0;
}
*/



//二分查找
/*
#include <iostream>
#include <algorithm>
#include <cstdio>
using namespace std;

typedef long long ll;

const int N = 1e6 + 9;

ll Tbl[N];

int ErFen(ll x[], int n, int m);

int main()
{
int n,m;
cin >> n >> m;
for (int i = 0; i < n; i++)
cin >> Tbl[i];
sort(Tbl, Tbl + n);
int t = ErFen(Tbl, n, m);
cout << t << ' ' << Tbl[t] << endl;
system("pause");
return 0;
}

int ErFen(ll x[],int n,int m)
{
int L, R, mid;
L = 0;
R = n - 1;
while (L <= R)
{
mid = (L + R) / 2;
if (x[mid] > m)
R = mid - 1;
else if (x[mid] < m)
L = mid + 1;
else
return mid;
}
return (-1);
}
*/



//L2-004月饼
/*
#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;

typedef struct Moon {
double inven;//库存
double price;//总价
double profits;//单价
}Cake;

Cake y[1009];//声明Cake数组

void PaiXu(Cake x[], int n);//声明排序

int main(){
int n;//种类
double m,sum;//总价，最后价格
sum = 0;
cin >> n >> m;//输入
for (int i = 0; i < n; i++) {
cin >> y[i].inven;//输入库存
}
for (int i = 0; i < n; i++){
cin >> y[i].price;//输入总价
y[i].profits = y[i].price / y[i].inven;//计算单价
}
PaiXu(y,n);//排序
for (int i = 0; i < n; i++){
if (m <= y[i].inven){
sum += m * y[i].profits;
break;
}
else{
sum += y[i].price;
m = m - y[i].inven;
}
}
if (sum != 0)
printf("%.2f\n", sum);
else
printf("0\n");
system("pause");
return 0;
}

void PaiXu(Cake x[],int n) {
int maxl; Cake t;
for (int i = 0; i < n; i++) {
maxl = i;
for (int j = i; j < n; j++) {
if (x[maxl].profits < x[j].profits)
maxl = j;
}
if (maxl != i) {
t = x[i];
x[i] = x[maxl];
x[maxl] = t;
}
}
}
*/
