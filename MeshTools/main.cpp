#include <iostream>
#include <vector>
#include <array>
#include <map>
#include<string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cmath>

#include <Eigen/Dense>
#include "OBJ_Loader.h"	

const float pi = 3.14;

struct Mesh 
{
	std::vector<Eigen::Vector3f> positions;//焊接后的顶点
	std::vector<Eigen::Vector3i> faces;//三角面的索引
};

// 一条边用"无序顶点对"表示
using EdgeKey = std::pair<int, int>;

//边和顶点的邻接顶点
struct Adjacency
{
	std::map<EdgeKey, std::vector<int>> edges;//边及其对顶点
	std::vector<std::vector<int>> neighbors;//顶点及其邻居顶点
};

Mesh LoadAndWeld(const std::string& path)
{
	Mesh mesh;
	objl::Loader loader;
	if (!loader.LoadFile(path)) return mesh;
	std::map<std::array<float, 3>, int> lut;//位置和对应的新编号

	std::vector<int> remap(loader.LoadedVertices.size());//旧编号索引->新编号索引

	for (size_t i = 0; i < loader.LoadedVertices.size(); ++i)//遍历每个顶点
	{
		const auto& p = loader.LoadedVertices[i].Position;//记录每个顶点的信息
		std::array<float, 3> key{ p.X, p.Y, p.Z };//存储为位置

		auto it = lut.find(key);//查找是否存在位置
		if (it == lut.end()) {  //不存在
			int newindex = static_cast<int>(mesh.positions.size());//记录当前顶点表长度
			mesh.positions.push_back(Eigen::Vector3f(p.X, p.Y, p.Z));//顶点表添加新的位置
			it = lut.emplace(key, newindex).first;//添加到索引表中并返回添加的迭代器(key,value)
		}
		remap[i] = it->second;//把迭代器记录的value存进新编号序列
	}

	for (size_t i = 0; i < loader.LoadedIndices.size(); i += 3)
	{
		mesh.faces.push_back(Eigen::Vector3i(
			remap[loader.LoadedIndices[i]],
			remap[loader.LoadedIndices[i + 1]],
			remap[loader.LoadedIndices[i + 2]]
		)
		);
	}

	return mesh;
}

Adjacency buildAdjacency(const Mesh& mesh) 
{
	Adjacency adj;
	adj.neighbors.resize(mesh.positions.size());//开邻居表

	for (size_t f = 0; f < mesh.faces.size(); ++f)//遍历面，登记边
	{
		int a = mesh.faces[f][0];
		int b = mesh.faces[f][1];
		int c = mesh.faces[f][2];

		EdgeKey key_ab = { std::min(a,b),std::max(a,b) };
		adj.edges[key_ab].push_back(c);
		EdgeKey key_bc = { std::min(c,b),std::max(c,b) };
		adj.edges[key_bc].push_back(a);
		EdgeKey key_ac = { std::min(a,c),std::max(a,c) };
		adj.edges[key_ac].push_back(b);
	}

	for (const auto& kv : adj.edges)//从边推出邻接表
	{
		int i = kv.first.first;
		int j = kv.first.second;
		adj.neighbors[i].push_back(j);//一条边，其邻居顶点就是另一个端点
		adj.neighbors[j].push_back(i);
	}

	return adj;
}

Mesh subdivideLoop(const Mesh& mesh, const Adjacency& adj)//输入：原网格 + 邻接表；输出：细分一层的网格
{
	Mesh out;
	const int V = static_cast<int>(mesh.positions.size());
	const int E = static_cast<int>(adj.edges.size());

	out.positions.resize(V + E);
	out.faces.reserve(mesh.faces.size() * 4);

	//给边编号至[V,..]并计算边点位置
	std::map<EdgeKey, int> egdeIndex;
	int k = 0;
	for (const auto& kv : adj.edges)
	{
		int idx = V + k++;
		egdeIndex[kv.first] = idx;

		Eigen::Vector3f A = mesh.positions[kv.first.first];
		Eigen::Vector3f B = mesh.positions[kv.first.second];
		Eigen::Vector3f C = mesh.positions[kv.second[0]];
		Eigen::Vector3f D = mesh.positions[kv.second[1]];
		Eigen::Vector3f ep = 0.375f * (A + B) + 0.125f * (C + D);


		out.positions[idx] = ep;
	}

	//给旧顶点编号为[0,V-1]
	for (int v = 0; v < V; v++)
	{
		float n = static_cast<float> (adj.neighbors[v].size());
		float c = 0.375f + 0.25f * std::cos(6.28318530718f / n);
		float beta = (0.625f - c * c) / n;
		Eigen::Vector3f sum = Eigen::Vector3f::Zero();
		for (int nb : adj.neighbors[v])
		{
			sum += mesh.positions[nb];
		}
		out.positions[v] = (1.f - n * beta )* mesh.positions[v] + beta * sum;
	}

	for (const auto& f : mesh.faces)
	{
		int a = f[0], b = f[1], c = f[2];
		int ab = egdeIndex.at({ std::min(a,b),std::max(a,b) });
		int bc = egdeIndex.at({ std::min(c,b),std::max(c,b) });
		int ca = egdeIndex.at({ std::min(a,c),std::max(a,c) });

		out.faces.push_back(Eigen::Vector3i(a, ab, ca));   // 角 a
		out.faces.push_back(Eigen::Vector3i(b, bc, ab));   // 角 b
		out.faces.push_back(Eigen::Vector3i(c, ca, bc));   // 角 c
		out.faces.push_back(Eigen::Vector3i(ab, bc, ca));  // 中心
	}

	return out;
}

Mesh QEM(const Mesh& mesh, const Adjacency& adj, int targetFaces)
{
	Mesh work = mesh;                       // 工作副本，后面一直改它
	const int V = static_cast<int>(work.positions.size());
	std::vector<bool> alive(V, true);       // 顶点是否还活着

	// ---------- ① 顶点 Q 矩阵 ----------
	std::vector<Eigen::Matrix4d> Q(V, Eigen::Matrix4d::Zero());
	for (const auto& f : work.faces)
	{
		// TODO(1)：
		//   p0,p1,p2 取坐标（double）
		//   n = (p1-p0).cross(p2-p0)，归一化；d = -n.dot(p0)
		//   p4 = (n.x, n.y, n.z, d)
		//   K = p4 * p4.transpose()
		//   Q[f[0]] += K; Q[f[1]] += K; Q[f[2]] += K;
		Eigen::Vector3d p0 = work.positions[f[0]].cast<double>();
		Eigen::Vector3d p1 = work.positions[f[1]].cast<double>();
		Eigen::Vector3d p2 = work.positions[f[2]].cast<double>();

		Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
		double len = n.norm();
		if (len < 1e-12) continue;      // 退化面，跳过
		n /= len;
		double  d = -n.dot(p0);

		Eigen::Vector4d p4(n.x(), n.y(), n.z(), d);
		Eigen::Matrix4d K = p4 * p4.transpose();

		Q[f[0]] += K;
		Q[f[1]] += K;
		Q[f[2]] += K;
	}

	int iter = 0;
	Eigen::Vector3d bestX = Eigen::Vector3d::Zero();

	while (static_cast<int>(work.faces.size()) > targetFaces && iter < 20000)
	{
		Adjacency a = buildAdjacency(work);

		// ---------- ② 找代价最小的可塌边 ----------
		double bestCost = std::numeric_limits<double>::max();
		EdgeKey best{ -1, -1 };
		for (const auto& kv : a.edges)
		{
			int u = kv.first.first, v = kv.first.second;
			if (u == v) continue;

			// TODO(2)：
			//   n 个公共邻居 != 2 就跳过（link condition）
			//   x = 中点
			//   cost = x^T (Q[u]+Q[v]) x
			//   比 bestCost 小就记下
			int common = 0;
			for (int nb : a.neighbors[u]) {
				for (int nb2 : a.neighbors[v]) {
					if (nb == nb2) ++common;
				}
			}
			if (common != 2) continue;

			Eigen::Vector3d x = 0.5 * (work.positions[u].cast<double>()
				+ work.positions[v].cast<double>());
			Eigen::Vector4d x4(x.x(), x.y(), x.z(), 1.0);
			Eigen::Matrix4d Qbar = Q[u] + Q[v];
			double cost = x4.dot(Qbar * x4);     // = x^T · Qbar · x

			if (cost < bestCost)
			{
				bestCost = cost;
				best = { u, v };
				bestX = x;          // ★ 关键
			}

		}
		if (best.first < 0) break;          // 没有可塌的了

		// ---------- ③ 执行塌缩 ----------
		int u = best.first, v = best.second;
		// TODO(3)：
		//   work.positions[u] = x
		//   Q[u] += Q[v]
		//   遍历 faces：先 v->u 替换，再删除含重复顶点的面
		//   alive[v] = false
		work.positions[u] = bestX.cast<float>();
		Q[u] += Q[v];
		for (auto& f : work.faces) {
			for (int i = 0; i < 3; ++i) {
				if (f[i] == v) f[i] = u;
			}
		}
		work.faces.erase(
			std::remove_if(work.faces.begin(), work.faces.end(),
				[](const Eigen::Vector3i& f)
				{
					return f[0] == f[1] || f[1] == f[2] || f[2] == f[0];
				}),
			work.faces.end());
		alive[v] = false;

		++iter;
		if (iter % 500 == 0)
			std::cout << "iter " << iter << "  faces " << work.faces.size() << "\n";
	}
	return work;
}


// 去掉没有被任何面引用的顶点，并重映射索引
Mesh compactMesh(const Mesh& m)
{
	Mesh out;
	std::vector<int> remap(m.positions.size(), -1);
	for (const auto& f : m.faces)
		for (int i = 0; i < 3; ++i)
			if (remap[f[i]] < 0)
			{
				remap[f[i]] = static_cast<int>(out.positions.size());
				out.positions.push_back(m.positions[f[i]]);
			}
	for (const auto& f : m.faces)
		out.faces.push_back(Eigen::Vector3i(remap[f[0]], remap[f[1]], remap[f[2]]));
	return out;
}

void writeOBJ(const std::string& path, const Mesh& mesh)
{
	std::ofstream o(path);
	if (!o) { std::cout << "!! write failed: " << path << "\n"; return; }
	o << "# vertices " << mesh.positions.size()
	  << " faces " << mesh.faces.size() << "\n";
	for (const auto& p : mesh.positions)
		o << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
	for (const auto& f : mesh.faces)
		o << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << "\n";
}

void printBBox(const char* tag, const Mesh& mesh)
{
	if (mesh.positions.empty()) return;
	Eigen::Vector3f mn = mesh.positions[0], mx = mesh.positions[0];
	for (const auto& p : mesh.positions)
	{
		mn = mn.cwiseMin(p);
		mx = mx.cwiseMax(p);
	}
	std::cout << "  " << tag << " bbox min " << mn.transpose()
		<< "   max " << mx.transpose() << "\n";
}

void printMeshStats(const Mesh& mesh)
{
	std::cout << "顶点 = " << mesh.positions.size()
		<< "   三角形 = " << mesh.faces.size() << "\n";
}

void printAdjStats(const Mesh& mesh, const Adjacency& adj)
{
	// TODO：打印
	//   ① 边数            （期望 8784）
	//   ② 对顶点数 != 2 的边有多少条（期望 0）
	//   ③ 平均度数        （期望 ≈ 6.0）
	std::cout << "边数 =" << adj.edges.size() << std::endl;
	int i = 0;
	for (const auto& v : adj.edges) 
	{
		if (v.second.size() != 2) { ++i; }
	}
	std::cout << "对顶点数 != 2 的边数 =" << i << std::endl;
	const auto& pnum = mesh.positions.size();
	int total = 0;
	for (const auto& v : adj.neighbors)
	{
		total += v.size();
	}
	std::cout << "平均度数 =" << static_cast<double>(total) / pnum << std::endl;
}

void printSubdivide(const Mesh& subout)
{
	std::cout << "细分之后的顶点数 = " << subout.positions.size()<<std::endl;
	std::cout << "细分之后的面数 = " << subout.faces.size() << std::endl;
}

void printQEM(const Mesh& QEMout)
{
	std::cout << "简化之后的顶点数 = " << QEMout.positions.size() << std::endl;
	std::cout << "简化之后的面数 = " << QEMout.faces.size() << std::endl;
}

int main() 
{
	const std::string path =
		"E:/Code/01_Study/games101/Assignment3/models/spot/spot_triangulated_good.obj";
	Mesh mesh = LoadAndWeld(path);
	Adjacency adj = buildAdjacency(mesh);

	printMeshStats(mesh);
	printAdjStats(mesh, adj);

	// L0 = 焊接后的原模型（2930 顶点 / 5856 面）
	writeOBJ("E:/Temp/spot_L0.obj", mesh);
	printBBox("L0", mesh);

	// L1 = 细分一级（11714 / 23424）
	Mesh L1 = subdivideLoop(mesh, adj);
	printSubdivide(L1);
	writeOBJ("E:/Temp/spot_L1.obj", L1);
	printBBox("L1", L1);

	// L2 = 再细分一级（46850 / 93696）
	Adjacency adj1 = buildAdjacency(L1);
	Mesh L2 = subdivideLoop(L1, adj1);
	printSubdivide(L2);
	writeOBJ("E:/Temp/spot_L2.obj", L2);
	printBBox("L2", L2);

	// QEM 减面：5856 -> 1000
	Mesh qem = QEM(mesh, adj, 1000);
	Mesh qemClean = compactMesh(qem);
	printQEM(qemClean);
	writeOBJ("E:/Temp/spot_qem1000.obj", qemClean);
	printBBox("QEM", qemClean);

	// QEM 减面：5856 -> 2000
	Mesh qem2k = compactMesh(QEM(mesh, adj, 2000));
	printQEM(qem2k);
	writeOBJ("E:/Temp/spot_qem2000.obj", qem2k);
	printBBox("QEM2K", qem2k);

	return 0;
}
