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
#include <omp.h>

#define INF 99999999
//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;

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

struct OrderCompp{//prior to reture the vertex with smaller order
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
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert), i.e., one-hop. there still be a case that the dis is also supported by neighbor's label even it is true.
	set<int> DisRe;//record the vertex id that the distance label should be updated
	vector<int> ch;
	int height=0;//tree height of a tree node
    int hdepth=0;//hdepty is the deepest node that a vertex still exists
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
    int partiNum=16;   //partition number

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    string algoParti="NC";

	//vertex order
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID

    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    vector<map<int, vector<int>>> SCconNodesMT;// <ID1,ID2,<supportive vertices>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
    vector<vector<int>> VidtoTNid;// (nodeID,vector<tree node rank>), record the tree node id whose unique vertex involves this vertex as neighbor
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
    vector<Node> Tree;
    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//For LCA query


    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }

    //// For H2H index constrction
    void H2HIndexConstruct(); //H2H index construction
    /// MDE contraction
    void MDEContraction(string orderfile);
    //for order generation
    void deleteEOrderGenerate(int u,int v);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    void insertEMTorder(int u,int v,int w);
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    /// Tree building
    void makeTree();
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    /// Index building
    void makeIndex();
    void makeRMQ();
    void makeRMQDFS(int p, int height);
    void makeIndexDFS(int p, vector<int>& list);

    void IndexsizeH2H();  //Core-tree index size computation

    //// For H2H query processing
    void CorrectnessCheckH2H(int runtimes);
    void EffiCheckH2H(string filename,int runtimes);
    int QueryH2H(int ID1,int ID2);
    int LCAQuery(int _p, int _q);
    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    void RetrievePath(int ID1, int ID2, vector<int> & prece);

    //// For H2H Index Maintenance
    void IndexMaintenanceH2H(int updateType, int updateBatch, bool ifBatch, int batchSize);
    void H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
    void H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel);



	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
	void UpdateGene(int num, string filename);

    void WriteOrder(string filename);
    void ReadOrder(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEAD_H_
