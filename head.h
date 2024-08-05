/*
 * head.h
 *
 *  Created on: 24 August 2022
 *      Author: Xinjie ZHOU
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include "Heap.h"
//#include "labeling.hpp"
//#include <omp.h>

#define INF 99999999
#define NO_Boundary 0
#define Post_Boundary 1

//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;
extern vector<int> _DD2_;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderComp{// maximum-first, Higher-order first
    int ID;
    OrderComp(){ID=0;}
    OrderComp(int _ID){
        ID=_ID;
    }
    bool operator< (const OrderComp d) const{
        if(NodeOrder_[ID]!=NodeOrder_[d.ID])
            return NodeOrder_[ID]>NodeOrder_[d.ID];
        return ID>d.ID;
    }
};

struct OrderCompp{//prior to return the vertex with smaller order
    int x;
    OrderCompp(int _x){
        x=_x;
    }
    bool operator< (const OrderCompp& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};

struct OrderComp3{
    int ID1, ID2;
    OrderComp3(){
        ID1=0, ID2=0;
    }
    OrderComp3(int _ID1, int _ID2){
        ID1=_ID1;
        ID2=_ID2;
    }
    bool operator< (const OrderComp3 d) const{//return the larger order vertex
        return (NodeOrder_[ID1] > NodeOrder_[d.ID1]) || ((NodeOrder_[ID1] <= NodeOrder_[d.ID1]) && (NodeOrder_[ID2] > NodeOrder_[d.ID2]) );
    }
};

//return smallest Dis with largest-order vertex
struct MinComp{
    int s, t, Dis;//s is the source vertex while t is the target vertex
    MinComp(){
        s=0, t=0, Dis=0;
    }
    MinComp(int _ID1, int _ID2, int _Dis){
        s=_ID1; t=_ID2; Dis=_Dis;
    }
    bool operator< (const MinComp d) const{
        if(Dis != d.Dis){
            return Dis<d.Dis;
        }else{
            if(NodeOrder_[s]!=NodeOrder_[d.s]){
                return NodeOrder_[s]>NodeOrder_[d.s];
            }
            return NodeOrder_[t]>NodeOrder_[d.t];
        }
    }
};


struct DegComp{//min-first
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator < (const DegComp d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct DegComp2{//min-first
    int x;
    DegComp2(int _x){
        x=_x;
    }
    bool operator < (const DegComp2 d) const{
        if(_DD2_[x]!=_DD2_[d.x])
            return _DD2_[x]<_DD2_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
//    vector<pair<int,pair<int,int>>> vertNo;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<int> pos;
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
//    vector<int> disNo;//the distance value of post-boundary strategy
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
	//vector<set<int>> FromNode;
//	set<int> changedPos;
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
	set<int> DisRe;//record the vertex id that the distance label should be updated
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0
	int uniqueVertex;//?vertex id of this tree node?
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
//		neighInf.clear();
		pos.clear();
		dis.clear();
		cnt.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
//		changedPos.clear();
		FN.clear();
		DisRe.clear();
//		piv.clear();
//        treeroot=-1;
	}
};

class Graph{
public:
    string graphfile;
    string dataset;
	int node_num=0;    //vertex number
	unsigned long long edge_num=0;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
//    vector<unordered_map<vertex,int>> NeighborMap;//original graph, map version
    vector<vector<pair<vertex,int>>> NeighborsParti;//<node_number,<in-partition adjacency lists>>
    vector<unordered_map<vertex,int>> NeighborsPartiPost;//adjacency lists of post-boundary partitions
    vector<unordered_map<vertex,int>> NeighborsOverlay;//<node_number,<adjacency lists of overlay graph>>
//    vector<vector<pair<vertex,int>>> NeighborsOverlay;//<node_number,<adjacency lists of overlay graph>>
//    vector<unordered_map<int,int>> BoundaryShortcuts;
    vector<pair<int,bool>> PartiTag;//<node_number,<partition_id,if_boundary>>
    vector<vector<vertex>> PartiVertex;//<partition_number,<in-partition vertices>>, in increasing vertex order, higher-rank vertex first
    vector<vector<vertex>> BoundVertex;//boundary vertices of each partition
    vector<vertex> OverlayVertex;//overlay vertex in increasing vertex order
//    vector<unordered_map<vertex,pair<int,int>>> repairShortcuts;//<ID1,ID2,<overlay weight, in-partition weight>>
    int partiNum;   //partition number
    bool ifParallel = true;

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    int algoQuery=0;//algorithm for core update, (0: PDPLL; 1: SDPLL), default: 0
    int algoUpdate=0;//algorithm for core construction, (0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read)
    string algoParti="NC";

	//vertex order
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID


    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);;// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    //for overlay graph
    vector<Node> Tree;
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
//    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    //    vector<map<int, vector<int>>> SCconNodesMT;//supportive vertex, multiple thread of SCconNodes
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for overlay

    //for no-boundary partitions
    vector<vertex> IDMap;//map the old id to new id in partitions
    vector<vector<Node>> Trees;//Trees for no-boundary
    vector<vector<int>> toRMQs;
    vector<vector<vector<int>>> RMQIndexs;
    vector<vector<int>> ranks;
    vector<int> heightMaxs;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTP;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidP;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    //for post-boundary partitions
//    vector<vertex> IDMapPost;//map the old id to new id in partitions
    vector<vector<Node>> TreesPost;//Trees for post-boundary
    vector<vector<int>> toRMQsPost;
    vector<vector<vector<int>>> RMQIndexsPost;
    vector<vector<int>> ranksPost;
    vector<int> heightMaxsPost;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTPost;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidPost;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    vector<bool> ifRepaired;

    string sourcePath;

    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }


    /// For PH2H
    void PH2HIndexConstruct(); //PH2H index construction
    void ConstructBoundaryShortcutV(vector<int> & p);
    void ConstructBoundaryShortcut(int pid);
    void Construct_PartiIndex(bool ifParallel);
    void Construct_OverlayGraph(bool ifParallel);
    void Construct_OverlayIndex();

    void ConstructPartitionPost(bool ifParallel);
    void ConstructPostParti(int pid);
    void ConstructPostPartiV(vector<int>& p);
    void ConstructPartitionPostIndex(bool ifParallel);

    void Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);


    void ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID);
    void ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs);//Create tree for partition
    void H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP);//Create labels for partition

    void H2HCreateTree_Overlay();
    void H2HCreateIndex_Overlay();
    int ShortcutDisCheck(int ID1, int ID2);

    void makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP);

    void makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees);
    void makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees);
    void makeRMQCore();
    void makeRMQDFSCore(int p, int height, vector<int>& EulerSeq);
    void makeIndexDFS(int p, vector<int> &list);
//    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
//    void insertEMTorder(int u,int v,int w);
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void insertECoreMT(int u,int v,int w);
    int matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    int matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    void IndexSizePH2H();  //Core-tree index size computation
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);

    void PH2HVertexOrdering(int type);
    void SketchGraphBuild();
    void OverlayOrderingBuild();
    void OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor);
    void PartitionOrderingBuildMDE(bool ifParallel);
    void OrderingAssemblyMDEBoundaryFirst(int pNum);
    void OrderingAssemblyMDE(int pNum);
    void OrderingAssemblyBoundaryFirst(int pNum);
    void SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch);

    void PartitionOrderingV(vector<int>& p);
    void PartitionOrdering(int pid);
    void deleteEOrderGenerate(int u,int v);
    void insertEOrderGenerate(int u,int v,int w);

//    vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
    vector<vector<int>> NeighborSketch;
    vector<set<int>> NeighborSketchS;
    vector<map<int,int>> vNodeOrderParti;
    vector<int> vNodeOrderOverlay;

	///Query processing
    int QueryH2HPartition(int ID1, int ID2, int PID);
    int QueryH2HPartitionPost(int ID1, int ID2, int PID);
    void EffiCheck(string filename,int runtimes);
	int Query(int ID1, int ID2);
    int QueryDebug(int ID1, int ID2);
	int QueryPartiCore(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
    int QuerySamePartiPost(int ID1, int ID2);
	int QueryPartiParti(int ID1, int ID2);
	int QueryCore(int ID1, int ID2);
    int LCAQuery(int _p, int _q);
    int LCAQueryPartition(int _p, int _q, int PID);// query within partition
    int LCAQueryPartitionPost(int _p, int _q, int PID);// query within partition
    int LCAQueryOverlay(int _p, int _q);

    //Correctness Check
    void CorrectnessCheck(int runtimes);
    void CorrectnessCheckCore(int runtimes);
    void DFSTree(vector<int>& tNodes, int id);


    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    void RetrievePath(int ID1, int ID2, vector<int> & prece);
    int DijkstraCore(int ID1, int ID2);


	/// Index update
    void IndexMaintenance(int updateType, bool ifBatch, int batchNumber, int batchSize);
    void DecreaseSingle(int a, int b, int oldW, int newW);//process one update
    void IncreaseSingle(int a, int b, int oldW, int newW);//process one update
    void DecreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void IncreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
    void EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank);
    void IncreaseOverlay(int a, int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void IncreaseParti(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax);//batch decrease for overlay graph
    void DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void DecreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void IncreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);


	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
    void ODGeneParti(int num, string filename);
    void ODGeneSameParti(int num, string filename);
    void ODGeneCrossParti(int num, string filename);
	void UpdateGene(int num, string filename);
    void QueryGenerationParti(bool ifSame, string filePath);


    void WriteTreeIndexOverlay(string filename);
    void ReadTreeIndex(string file);
    void WriteTreeIndexParti(string filename);
    void WriteGraph(string graphfile);
    void WriteOrder(string filename);
    void ReadOrder(string filename);
    void CompareOrder(string filename1, string filename2);
    void GraphPartitionRead(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEAD_H_
