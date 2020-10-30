#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <queue>
#include <set>
#include <unordered_map>
#include <cassert>
#include <fstream>
#include <functional>
#include <vector>
#include <time.h>
#include "pruned_landmark_labeling.h"
#include <windows.h>
#include <ctime>
#include <stdlib.h>
using namespace std;
string str = "D:\\dataset_comp\\Epinions.txt";//CA-HepPh,ca-AstroPh.mtx,in500_3810.txt,in19_23.txt
//graph1000_10000,Epinions,CA-CondMat,Slashdot
const int k = 64;
const int inf = 1 << 29;//1<<29 表示1左移29位，每左移一位乘以2，所以就是1*2^29=
int cnt;
int n = 0, m = 0;
int *edg_size, *r_node, *pecc;
int c = 0;//number of unequal node
int **dis,*its;
struct Vertex *Vertex_array;
struct eccbound_ *eccbound, *peccbound;


vector<vector<pair<int, int>>> edg(1); // storing graph : to, cost
std::vector<std::vector<std::pair<int, int> > > Lz;
std::vector<int> id(1); // get vertex according to position
std::vector<std::pair<int, int>> node(1);
std::vector<std::pair<int, int>> el;
std::vector<std::pair<int, int>> R_id;

PrunedLandmarkLabeling<> pll;
typedef struct eccbound_ {
	int upper;
	int lower;
};
struct ListNode {
	int val;
	int dist;
	ListNode* next;
};
typedef struct Vertex {
	struct edg_ *next;
};
typedef struct edg_ {
	int first;
	int second;
	struct edg_ *next;
};
int BFS(int x)
{
	std::queue<pair<int, int>> q_bfs;
	//std::vector<int> flag(n + 1);
	struct edg_ *r;
	int *flag;
	flag = new int[n + 1];
	memset(flag, -1, sizeof(int) * (n + 1));
	int u, dist = 0;
	q_bfs.push(make_pair(0, x));
	flag[x] = 1;
	while (!q_bfs.empty())
	{
		u = q_bfs.front().second;
		dist = q_bfs.front().first;
		q_bfs.pop();
		r = Vertex_array[u].next;
		while (r != NULL)
			//for (int i = 0; i < edg[u].size(); i++)
		{
			if (flag[r->first] != 1)
			{

				q_bfs.push(make_pair(dist + 1, r->first));
				flag[r->first] = 1;
			}
			r = r->next;
		}
	}
	return dist;
}
bool cmp(const pair<int, int> a, const pair<int, int> b) {
	return a.first > b.first; //降序排列，如果改为return a<b，则为升序
}
int LocateVex(int num) {
	for (int i = 1; i < id.size(); i++) {
		if (id[i] == num)
			return i;
	}
	id.emplace_back(num);
	edg.push_back(std::vector<std::pair<int, int> >(0));
	return (id.size() - 1);
}

bool delete_edge(int uId, int vId) {
	for (int i = 0; i < edg[uId].size(); i++)
	{
		if (edg[uId][i].first == vId)
			return false;
	}
	return true;
}
void graph_init(ifstream &ifs) {
	assert(ifs);
	int u, v, uId, vId;
	bool flag;
	struct edg_ *e, *r;

	Lz.resize(k + 1, vector<pair<int, int>>(0));

	while (!ifs.eof()) {
		//cout << " m=" << m << endl;
		ifs >> u >> v;
		uId = LocateVex(u);
		vId = LocateVex(v);

		flag = delete_edge(uId, vId);
		if (flag == true) {
			m++;
			//cout << "u=" << u << " " << "v=" << v << endl;
			edg[uId].emplace_back(vId, 1); edg[vId].emplace_back(uId, 1);
			el.push_back(make_pair(vId, uId));
		}
	}
	n = id.size() - 1;
	cout << "n=" << n << endl;
	cout << "m=" << m << endl;
	/*dis = new int*[k + 1];
	for (int i = 1; i <= k; ++i) {
		dis[i] = new int[n + 1];
		memset(dis[i], -1, sizeof(int) * (n + 1));
	}
	*/
	pecc = new int[k + 1];
	Vertex_array = new Vertex[n + 1];
	//edg_size = new int[n + 1];1260357639
	for (int i = 1; i <= n; i++) {
		Vertex_array[i].next = NULL;
		e = (edg_*)malloc(sizeof(edg_));
		Vertex_array[i].next = e;
		e->first = edg[i][0].first;
		e->second = edg[i][0].second;
		e->next = NULL;
		r = e;
		for (int j = 1; j < edg[i].size(); j++) {
			e = (edg_*)malloc(sizeof(edg_));
			e->first = edg[i][j].first;
			e->second = edg[i][j].second;
			e->next = NULL;

			r->next = e;
			r = e;//r->next
		}
	}
	peccbound = new eccbound_[k + 1];
	eccbound = new eccbound_[n + 1];
	for (int i = 0; i < n + 1; i++) {
		eccbound[i].lower = inf;
		eccbound[i].upper = 0;
		if (i >= 1 && i <= k) {
			peccbound[i].lower = inf;
			peccbound[i].upper = 0;
		}
	}

	for (int i = 1; i <= n; i++)
		node.push_back(make_pair(edg[i].size(), i));

	sort(node.begin() + 1, node.end(), cmp);
	r_node = new int[n + 1];
	for (int i = 1; i <= n; ++i)
		r_node[node[i].second] = i;
	/*for (int i = 1; i <= n; i++)
		cout << "node[" << i << "].first=" <<
				node[i].first << "," << "node[" << i << "].second=" << id[node[i].second] << endl;
	*/
}
bool update_lb(int select, int lb) {
	if (lb > peccbound[select].lower) { peccbound[select].lower = lb; return true; }
	return false;
}
bool update_ub(int select, int ub) {
	if (ub < peccbound[select].upper) { peccbound[select].upper = ub; return true; }
	return false;
}
void Refpool() {
	int z, dist = 0;
	int *flag, count = 0;
	flag = new int[n + 1]();
	its = new int[k + 1];
	int i = 0;
	edg_ *r, *e;

	for (int j = 1; j <= k; ++j) {
		while (i <= n) {
			if (flag[node[i].second]) { i++; continue; }
			if (!node[i].first)
				i++;
			else {//用作参考节点
				flag[node[i].second] = j;
				count++;

				Lz[j].emplace_back(0, node[i].second);
				e = Vertex_array[node[i].second].next;//e指向参考节点的邻居
				while (e != NULL) {
					//cout << "e->first=" << id[e->first] << endl;

					if (!flag[e->first]) {
						flag[e->first] = j;//邻居设为Lz[k]的元素
						Lz[j].emplace_back(1, e->first);
						count++;
						node[r_node[e->first]].first--;
					}
					else { e = e->next; continue; }
					r = Vertex_array[e->first].next;//邻居的邻居度数减1
					while (r != NULL) {
						if (node[r_node[r->first]].first)
							node[r_node[r->first]].first--;
						r = r->next;
						/*for (int i = 1; i <= n; i++)
						{
							cout << "node[" << i << "].first=" << node[i].first << "," << "node[" << i << "].second=" << id[node[i].second] << endl;
						}
						cout << endl;*/
					}
					e = e->next;
				}
				pecc[j] = dist + 1;
				i++;
				break;
			}
		}
	}

	int *add = new int[k + 1]();
	bool f = true;
	while (count < n&&f) {
		dist++;
		f = false;
		for (int j = 1; j <= k; ++j) {
			//cout << "第" << dist << "次BFS" << endl;
			if (add[j] == Lz[j].size() - 1) {
				//pecc[j] = dist;
				continue;
			}
			add[j]++;

			while (Lz[j][add[j]].first == dist && count < n) {
				e = Vertex_array[Lz[j][add[j]].second].next;
				while (e != NULL && count < n) {
					if (!flag[e->first]) {//未访问
						//cout << "访问" << id[e->first] << endl;
						Lz[j].emplace_back(dist + 1, e->first);
						count++;
						flag[e->first] = j;
						f = true;
					}
					e = e->next;
				}
				if (add[j] == Lz[j].size() - 1 || Lz[j][add[j] + 1].first != dist) break;
				add[j]++;
			}
			if (add[j] != Lz[j].size() - 1)
				pecc[j] = dist + 1;
			if (count == n) break;
		}
	}



	cout << "final ref " << count << endl;
	/*for (int i = 1; i <= k; ++i){
		cout << i<<",pecc=" << pecc[i] << endl;
		for (int j = 0; j < Lz[i].size(); ++j){
			cout << "Lz[" << i << "][" << j << "]=" << Lz[i][j].first << "," << id[Lz[i][j].second] << " ";
			if (!(j % 9) && j)
				cout << endl;
		}
		cout << endl;
	}*/
}
int query_refpool(int max_lb) {
	//查询部分偏心率上下界最大的部分
	int select = -1, min_b = -1, max_b = 0;//记录最大部分的参考节点位置
	for (int j = 1; j <= k; ++j) {
		if (max_lb >= peccbound[j].upper) {
			peccbound[j].upper = -1;
		}
		if (max_b <= peccbound[j].upper)
			if (max_b == peccbound[j].upper) {
				if (min_b < peccbound[j].lower) {
					min_b = peccbound[j].lower;
					select = j;
				}
			}
			else {
				max_b = peccbound[j].upper;
				min_b = peccbound[j].lower;
				select = j;
			}
	}
	return select;
}

void Ecc_Comp(int x) {
	//找出最大范围的部分偏心率上下界
	int *dist_x_z = new int[k + 1];
	int dist, ecc_temp = -1, oldDist = INT_MAX;
	int select = -1, max_lb = -1, max_b = -1;//记录最大部分的参考节点位置
	int dist_z_i = -1;
	bool f;
	for (int j = 1; j <= k; ++j) {
		dist_x_z[j] = pll.QueryDistance(Lz[j][0].second, x);
		if (dist_x_z[j] != INT_MAX) {
			peccbound[j].upper = dist_x_z[j] + pecc[j];
			//cout << "dist_x_z[j]=" << dist_x_z[j] << endl;
			peccbound[j].lower = max(dist_x_z[j], pecc[j] - dist_x_z[j]);
			if (max_lb < peccbound[j].lower)
				max_lb = peccbound[j].lower;
			//cout << j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
		}
		else {
			peccbound[j].upper = -1;
		}
	}
	//cout << "max_lb=" << max_lb << endl;
	select = query_refpool(max_lb);
	if (select == -1) return;
	//cout << "secect=" << select << endl;
	//for (int j = 1; j <= k; ++j)
		//cout << j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
	for (int i = Lz[select].size() - 1; i >= 0; --i) {
		//for (int i = 0 ; i < Lz[select].size() ; ++i) {
		dist_z_i = Lz[select][i].first;
		//cout << i <<","<<peccbound[select].lower << "," << dist_x_z[select] << ","<< dist_z_i << endl;
		//cout << "oldDist=" << oldDist << "," << dist_z_i << endl;//cout << "Lz[select][i].second=" << id[Lz[select][i].second] << endl;

		if (oldDist > dist_z_i) {

			oldDist = dist_z_i;
			update_ub(select, max(peccbound[select].lower, dist_x_z[select] + dist_z_i));
		}//list_land[l][i].val==λ
		dist = pll.QueryDistance(x, Lz[select][i].second);

		//cout << "dist=" << dist << endl;

		update_lb(select, dist);
		//if (f&&peccbound[select].lower >= Lz[select][i].first)
			//peccbound[select].upper == peccbound[select].lower;
		//cout << select << "-> " << peccbound[select].lower << "," << peccbound[select].upper << endl;

		if (peccbound[select].lower == peccbound[select].upper) {
			//cout << "ecc[" << select << "]=" << peccbound[select].lower << "," << peccbound[select].upper << endl;
			if (ecc_temp < peccbound[select].lower)
				ecc_temp = peccbound[select].lower;
			//cout << "        ecc_temp=" << ecc_temp << endl;
			for (int j = 1; j <= k; ++j) {
				if (peccbound[j].upper == -1)continue;
				if (ecc_temp < peccbound[j].upper) {
					select = j;
					i = Lz[select].size();//-1
					oldDist = INT_MAX;
					//cout << "secect_again=" << select << endl;
					break;
				}
			}
			//for (int j = 1; j <= k; ++j) 
				//cout << "  "<<j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
			if (i != Lz[select].size()) {
				//cout << "break!\n";
				break;
			}
		}
	}
	eccbound[x].lower = ecc_temp;
	eccbound[x].upper = ecc_temp;
	//return;
}
void Ecc_Comp_d(int x) {
	//找出最大范围的部分偏心率上下界
	int *dist_x_z = new int[k + 1];
	int dist, ecc_temp = -1, oldDist = INT_MAX;
	int select = -1, max_lb = -1, max_b = -1;//记录最大部分的参考节点位置
	int dist_z_i = -1;
	bool f;
	for (int j = 1; j <= k; ++j) {
		dist_x_z[j] = pll.QueryDistance(Lz[j][0].second, x);
		if (dist_x_z[j] != INT_MAX) {
			peccbound[j].upper = dist_x_z[j] + pecc[j];
			//cout << "dist_x_z[j]=" << dist_x_z[j] << endl;
			peccbound[j].lower = max(dist_x_z[j], pecc[j] - dist_x_z[j]);
			if (max_lb < peccbound[j].lower)
				max_lb = peccbound[j].lower;
			//cout << j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
		}
		else {
			peccbound[j].upper = -1;
		}
	}
	//cout << "max_lb=" << max_lb << endl;
	select = query_refpool(max_lb);
	if (select == -1) return;
	//cout << "secect=" << select << endl;
	//for (int j = 1; j <= k; ++j)
		//cout << j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
	//cout << endl;
	//得出pecc(x|W)
	int pecc_lw, pecc_rw, t;
	while (1) {
		pecc_lw = INT_MAX; pecc_rw = -1; t = 0; f = false;
		while (pecc_lw > pecc_rw) {
			dist = pll.QueryDistance(x, Lz[select][t].second);
			//cout << "  dist=" << dist << endl;
			if (Lz[select][t].first + dist_x_z[select] <= dist) {
				pecc_rw = dist;
				pecc_lw = pecc_rw;
			}
			else
				pecc_rw = max(pecc_rw, dist);
			//cout << "  pecc_rw=" << pecc_rw << endl;
			t++;
		}
		//cout << "pecc_lw=" << pecc_lw << ",pecc_rw" << pecc_rw << endl;
		//得出pecc(x|V')
		for (int i = Lz[select].size() - 1; i >= t; --i) {
			dist_z_i = Lz[select][i].first;
			update_ub(select, max(peccbound[select].lower, dist_x_z[select] + dist_z_i));
			dist = pll.QueryDistance(x, Lz[select][i].second);
			update_lb(select, dist);
			if (peccbound[select].lower == peccbound[select].upper) break;
			if (i == t && peccbound[select].lower != peccbound[select].upper) {
				peccbound[select].lower = peccbound[select].upper = pecc_rw;
				break;
			}
		}

		if (peccbound[select].lower == peccbound[select].upper) {

			ecc_temp = max(ecc_temp, peccbound[select].upper);
				for (int j = 1; j <= k; ++j) {
					if (peccbound[j].upper == -1)continue;
					if (ecc_temp < peccbound[j].upper) {
						select = j;
						f = true;
						//cout << "secect_again=" << select << endl;
						break;
					}
				}
			//for (int j = 1; j <= k; ++j) 
				//cout << "  "<<j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
			if (!f) {
				//cout << "break!\n";
				break;
			}
		}
	}
	eccbound[x].lower = ecc_temp;
	eccbound[x].upper = ecc_temp;
}
/*void Ecc_Comp_d(int x) {
	//找出最大范围的部分偏心率上下界
	int *dist_x_z = new int[k + 1];
	int dist, ecc_temp = -1, oldDist = INT_MAX;
	int select = -1, max_lb = -1, max_b = -1;//记录最大部分的参考节点位置
	int dist_z_i = -1;
	bool f_ub, f_lb;
	for (int j = 1; j <= k; ++j) {
		dist_x_z[j] = pll.QueryDistance(Lz[j][0].second, x);
		if (dist_x_z[j] != INT_MAX) {
			peccbound[j].upper = dist_x_z[j] + pecc[j];
			//cout << "dist_x_z[j]=" << dist_x_z[j] << endl;
			peccbound[j].lower = max(dist_x_z[j], pecc[j] - dist_x_z[j]);
			if (max_lb < peccbound[j].lower)
				max_lb = peccbound[j].lower;
			cout << j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
		}
		else {
			peccbound[j].upper = -1;
		}
		its[j] = Lz[j].size() - 1;
	}
	
	cout << "max_lb=" << max_lb << endl;
	select = query_refpool(max_lb);
	cout << "secect=" << select << endl;
	for (int j = 1; j <= k; ++j)
		cout <<"1---"<<j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
	for (int i = its[select]; i >= 0; --i) {//Lz[select].size() - 1
		//for (int i = 0 ; i < Lz[select].size() ; ++i) {
		dist_z_i = Lz[select][i].first;
		cout << i <<","<<peccbound[select].lower << "," << dist_x_z[select] << ","<< dist_z_i << endl;
		//cout << "oldDist=" << oldDist << "," << dist_z_i << endl;//cout << "Lz[select][i].second=" << id[Lz[select][i].second] << endl;

		//if (oldDist > dist_z_i) {

			//oldDist = dist_z_i;
			f_ub = update_ub(select, max(peccbound[select].lower, dist_x_z[select] + dist_z_i));
		//}//list_land[l][i].val==λ
		dist = pll.QueryDistance(x, Lz[select][i].second);

		//cout << "dist=" << dist << endl;

		f_lb = update_lb(select, dist);

		if (max_lb < peccbound[select].lower)
			max_lb = peccbound[select].lower;
		cout << select << "-> " << peccbound[select].lower << "," << peccbound[select].upper << endl;

		if (peccbound[select].lower == peccbound[select].upper) {
			//cout << "ecc[" << select << "]=" << peccbound[select].lower << "," << peccbound[select].upper << endl;
			if (ecc_temp < peccbound[select].lower)
				ecc_temp = peccbound[select].lower;
			//cout << "        ecc_temp=" << ecc_temp << endl;
			//its[select] = i;
			select = query_refpool(max_lb);
			if (select == -1)break;
			cout << "secect_again=" << select << endl;
			i = its[select];
			oldDist = INT_MAX;
			for (int j = 1; j <= k; ++j) 
				cout <<"2---"<<j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
			if (i != its[select]) {
				//cout << "break!\n";
				break;
			}
		}
		else if (f_ub || f_lb) {
			its[select] = i;
			select = query_refpool(max_lb);
			if (select == -1)break;
			i = its[select];
			oldDist = INT_MAX;
			for (int j = 1; j <= k; ++j)
				cout << "3---"<<j << "-> " << peccbound[j].lower << "," << peccbound[j].upper << endl;
		}
	}
	eccbound[x].lower = ecc_temp;
	eccbound[x].upper = ecc_temp;
}
*/
void ECC() {
	int lb, ub;
	for (int i = 1; i <= n; ++i)
	{
		lb = eccbound[i].lower;
		ub = eccbound[i].upper;
		if (lb != ub)
		{
			//cout << "计算" << i <<","<<id[i]<< endl;
			//Ecc_Comp(i);
			Ecc_Comp_d(i);
			//if (eccbound[i].second > lb) update_lb_dfs(i, eccbound[i].second);//lowbound增大
			//if (eccbound[i].first < ub) update_ub_dfs(i, eccbound[i].first);//upperbound减小
		}
		else
		{
			//cout << "不用计算：" << id[i] << endl;
			c++;
		}
	}
}
void test(int &n_e) {
	for (int i = 1; i <= n; i++) {
		int dist_bfs = BFS(i);
		cout << id[i] << ",bfs_ECC[" << i << "]=" << dist_bfs << endl;

		if (eccbound[i].lower == eccbound[i].upper) {
			//c++;
			//cout << "不用计算 " << i << endl;
			cout << id[i] << ",ecc[" << i << "]=" << eccbound[i].lower << "," << eccbound[i].upper << endl;

			if (dist_bfs != eccbound[i].lower)
			{
				cout << "不相等：" << i << "," << id[i] << endl;
				n_e++;
			}
			continue;
		}
		cout << id[i] << ",ecc[" << i << "]=" << eccbound[i].lower << "," << eccbound[i].upper << endl;
		if (dist_bfs != eccbound[i].lower)
		{
			cout << "不相等：" << i << "," << id[i] << endl;
			n_e++;
		}
	}
}
int main()
{
	ifstream ifs, ifs_main;
	clock_t start, end, start_1, end_1, start_2, end_2, start_3, end_3;
	int n_e = 0;
	ifs.open(str);
	cout << "Graph preparing..." << endl;
	graph_init(ifs);

	cout << "PLL preparing..." << endl;
	start_1 = clock();
	pll.ConstructIndex(el);
	end_1 = clock();

	cout << "Refpool preparing..." << endl;
	start_2 = clock();
	Refpool();
	end_2 = clock();

	cout << "ECC preparing..." << endl;
	start = clock();
	ECC();
	end = clock();
	cout << "Finished ECC" << endl;

	//test(n_e);
	/*for (int i = 1; i <= n; i++) {
		int dist_bfs = BFS(i);
		cout << id[i] << ",bfs_ECC[" << i << "]=" << dist_bfs << endl;

		if (eccbound[i].lower == eccbound[i].upper) {
			//c++;
			//cout << "不用计算 " << i << endl;
			cout << id[i] << ",ecc[" << i << "]=" << eccbound[i].lower << "," << eccbound[i].upper << endl;
			
			if (dist_bfs != eccbound[i].lower)
			{
				cout << "不相等：" << i << "," << id[i] << endl;
				n_e++;
			}
			continue;
		}
		cout << id[i] << ",ecc[" << i << "]=" << eccbound[i].lower << "," << eccbound[i].upper << endl;
		if (dist_bfs != eccbound[i].lower)
		{
			cout << "不相等：" << i << "," << id[i] << endl;
			n_e++;
		}
	}
	*/

	double endtime_1 = (double)(end_1 - start_1) / CLOCKS_PER_SEC;
	double endtime_2 = (double)(end_2 - start_2) / CLOCKS_PER_SEC;
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;

	cout << "PLL time:" << endtime_1 * 1000 << "ms" << endl;	//ms为单位
	cout << "Refpool time:" << endtime_2 * 1000 << "ms" << endl;	//ms为单位
	cout << "EccOneNode time:" << endtime * 1000 << "ms" << endl;	//ms为单位
	cout << "不用计算的个数等于 ： " << c << endl;
	cout << "不相等的个数为：" << n_e << endl;
	system("pause");
	return 0;
}
