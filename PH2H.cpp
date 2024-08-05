/*
 * Construction.cpp
 *
 *  Created on: 24 August 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"
#include "PH2H.hpp"

/// Index Construction
void Graph::PH2HIndexConstruct(){
    double runT1, runT2, runT3, runT4;
    runT1=0, runT2=0, runT3=0, runT4=0;

    /// Read order and partitions
    string orderfile=graphfile+".orderP";
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_order";
//    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE";
//    orderfile=graphfile+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
//    Construct_PartiIndex(false);
    Construct_PartiIndex(true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


//    g.WriteOrder(graphfile+".order");
//    WriteCoreGraph(graphfile+"C");
//    exit(0);
    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraph(true);
//    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
//    Construct_core(algoCoreC);
//    WriteCoreIndex(graphfile);
    Construct_OverlayIndex();

    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair

    if(algoQuery==1){

        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)



    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}


//function for computing the index size
void Graph::IndexSizePH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].dis.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }


    //Periphery index
    if(algoQuery==NO_Boundary){
        for(int pid=0;pid<Trees.size();++pid){
            for(int i=0;i<Trees[pid].size();i++){
                m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
                m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
                m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
            }
        }

        for(int i=0;i< SCconNodesMTP.size();i++){
            for(auto it=SCconNodesMTP[i].begin(); it!=SCconNodesMTP[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }
    else if(algoQuery==Post_Boundary){
        for(int pid=0;pid<TreesPost.size();++pid){
            for(int i=0;i<TreesPost[pid].size();i++){
                m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
                m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
                m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
                m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
            }
        }

        for(int i=0;i< SCconNodesMTPost.size();i++){
            for(auto it=SCconNodesMTPost[i].begin(); it!=SCconNodesMTPost[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }



    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overlay graph index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Partition graphs index size "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PH2HVertexOrdering(int type){
    ReadGraph(graphfile);
    int pNum=partiNum;
    switch (type) {
        case 0:{//MDE partition + distant MDE overlay
            cout<<"MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
//            PartitionOrderingBuildMDE(false);
            OrderingAssemblyMDE(pNum);
            break;
        }
        case 1:{
            cout<<"Boundary-first ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuildBoundaryFirst(node_num, NeighborSketch);
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyBoundaryFirst(pNum);
            break;
        }
        case 2:{
            cout<<"Boundary-first MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyMDEBoundaryFirst(pNum);
            break;
        }
        default:{
            cout<<"Wrong ordering type! "<<type<<endl; exit(1);
        }
    }
    exit(0);
}
void Graph::OrderingAssemblyMDEBoundaryFirst(int pNum){
    string filename=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE2";
    filename=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE2";
    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    vNodeOrder.clear();
    int order_i=0;
    set<int> vertices;
    /// For non-boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
            ID=it->second;
            if(!PartiTag[ID].second){// if not boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }
    /// For boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
            ID=it->second;
            if(PartiTag[ID].second){// if boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }

//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num || vertices.size()!=node_num || vNodeOrder.size()!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<vertices.size()<<" "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }

    for(int i=0;i<node_num;++i){
        NodeOrder[vNodeOrder[i]]=i;
    }

    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of MDE ordering assemblying
void Graph::OrderingAssemblyMDE(int pNum){
    string filename=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    int order_i=0;
    set<int> vertices;
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
            ID=it->second;
            vertices.insert(ID);
            NodeOrder[ID]=order_i;
            ++order_i;
        }
    }
//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<node_num<<endl; exit(1);
    }
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of boundary-first assemblying
void Graph::OrderingAssemblyBoundaryFirst(int pNum){
    string orderfile=graphfile+"_"+algoParti+"_"+to_string(pNum)+"/vertex_order2";
    set<int> vcheck;//use to check the redundant ordered vertex
    vcheck.clear();
    vNodeOrder.clear();
    int pid,ID;
    //for vertex within partition
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];

        for(auto it=vNodeOrderParti[k].begin();it!=vNodeOrderParti[k].end();++it){
            ID=it->second;
            if(!PartiTag[ID].second){//if not boundary vertex
                vNodeOrder.push_back(ID);
                if(vcheck.find(ID)!=vcheck.end())
                    cout<<"wrong: redundant vertex ordered"<<endl;
                vcheck.insert(ID);
            }
        }
    }

    //for boundary vertex
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID=BoundVertex[pid][i];
            vNodeOrder.push_back(ID);
            if(vcheck.find(ID)!=vcheck.end())
                cout<<"wrong: redundant vertex ordered"<<endl;
            vcheck.insert(ID);
        }
    }

    //cout<<"total number of ordered vertex "<<vNodeOrder.size()<<endl;
    if(vNodeOrder.size()!=node_num)
        cout<<"Something wrong happened: some vertices do not have the vertex order!"<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    ofstream OF(orderfile);
    if(!OF){
        cout<<"Cannot open file "<<orderfile<<endl;
        exit(1);
    }
    OF<<NodeOrder.size()<<endl;
    for(int i=0;i<NodeOrder.size();i++){
        OF<<i<<" "<<NodeOrder[i]<<endl;
    }
    OF.close();
    cout<<"Finished."<<endl;
}

void Graph::SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch){
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(1,1)));
    }

    _DD_.assign(node_num,0);
    DD.assign(node_num,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;

    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

        vNodeOrderSketch.push_back(x);
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}

void Graph::SketchGraphBuild(){
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));

    bool flag_minus = false;

    string filename=graphfile+"_"+algoParti+"_"+to_string(partiNum);
    filename=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<filename+"/subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoParti == "NC"){
//        flag_minus = true;
        partiNum = pnum2;
    }else if(algoParti == "SC" || algoParti == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;

    PartiVertex.assign(partiNum,vector<vertex>());
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    PartiVertex[k].emplace_back(ID);
                }else{
                    cout<<"vertex already in one partition!"<<ID<<" "<<PartiTag[ID].first<<" "<<k<<endl;
                }
            }else{
                cout<<"Wrong vertex ID! "<<ID<<endl; exit(1);
            }

        }
    }
    //further check that each vertex is in one and only one partition
    for(int vid=0;vid<node_num;vid++){
        if(PartiTag[vid].first==-1){
            cout<<"vertex "<<vid<<" not within any partition"<<endl; exit(1);
        }
    }
    int nNum=0;
    for(int pid=0;pid<partiNum;++pid){
        nNum+=PartiVertex[pid].size();
    }
    if(nNum!=node_num){
        cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
    for(int k=0;k<pnum1;k++){
        int edgenum0,ID1,ID2,weight;
        IF>>edgenum0;
        for(int i=0;i<edgenum0;i++){
            IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
                NeighborsParti[ID1].emplace_back(ID2,weight);
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    vector<int> vecint;
    vecint.clear();
    NeighborSketch.assign(partiNum, vecint);
//    vector<set<int>> NeighborSketchS;
    NeighborSketchS.assign(partiNum, set<int>());

    BoundVertex.assign(partiNum,vector<vertex>());
    //read the cut edges
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    int PID1, PID2;

    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            PID1=PartiTag[ID1].first, PID2=PartiTag[ID2].first;
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second){
                PartiTag[ID1].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);

            if(NeighborSketchS[PID1].find(PID2)==NeighborSketchS[PID1].end()){//if not found, i.e., PID2 is not in the NeighborSketchS of PID1
                NeighborSketch[PID1].push_back(PID2);
                NeighborSketch[PID2].push_back(PID1);

                NeighborSketchS[PID1].insert(PID2);
                NeighborSketchS[PID2].insert(PID1);
            }
        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }




    /*for(int k=0;k<pnum;k++){
        cout<<k<<" "<<NeighborSketch[k].size()<<endl;
    }*/
}

void Graph::OverlayOrderingBuild(){

    int lastParti = -1;
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<NeighborSketch.size();i++){
        for(auto it=NeighborSketch[i].begin();it!=NeighborSketch[i].end();++it){
            E[i].insert(make_pair(*it,make_pair(1,1)));
        }
    }
    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);
    _DD2_.assign(partiNum,0);

    set<DegComp> Deg;
    set<DegComp2> Deg2;
    int ID,degree;
    for(ID=0;ID<partiNum;ID++){
        degree=NeighborSketch[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<set<int>> NeighborCon(partiNum,set<int>());
    vector<int> neix;
    neix.clear();

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int x;

    set<int> pSet;
    while(!Deg.empty()){

        x=(*Deg.begin()).x;
        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }
        Deg.erase(Deg.begin());



        if(lastParti!=-1 && NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if x is lastParti's neighbor
            neix.clear();
            while(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found
                _DD2_[x]++;
                neix.emplace_back(x);

                if(Deg.empty()){
                    break;
                }else{
                    x=Deg.begin()->x;
                    Deg.erase(Deg.begin());
                }
            }

            if(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found, i.e., x is the neighbor of lastParti
                if(neix.size()>1){//if neix has more than one element
                    if(Deg.empty()){
                        Deg2.clear();
                        for(int i=0;i<neix.size();++i){
                            Deg2.insert(DegComp2(neix[i]));
                        }
                        x=Deg2.begin()->x;
                        Deg2.erase(Deg2.begin());
                        if(!Deg2.empty()){
                            for(auto it=Deg2.begin();it!=Deg2.end();++it){
                                ID=it->x;
                                Deg.insert(DegComp(ID));
                            }
                            Deg2.clear();
                        }
                    }else{
                        cout<<"Wrong! "<<endl; exit(1);
                    }
                }
            }//if not the neighbor
            else{
                if(!neix.empty()){
                    for(int i=0;i<neix.size();++i){
                        Deg.insert(DegComp(neix[i]));
                    }
                }
            }
        }

//        cout<<x<<" "<<Deg.size()<<endl;
        vNodeOrderOverlay.emplace_back(x);
        pSet.insert(x);
        lastParti = x;
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
                NeighborCon[x].insert(it->first);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
    }

//    exit(0);
}

//based on SketchOrder function, to guarantee the order of neighboring vertices are not contiguous
void Graph::OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor){
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j],make_pair(1,1)));
    }

    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<Neighbor.size();i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(partiNum,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;
    int x;
    int lastx;
    vector<int> neix;
    neix.clear();
    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;

        while(true){
            if(Deg.empty()){//in case that all the remaining vertices are the neighbors of current vertex x
                x=neix[0];
                for(int j=1;j<neix.size();j++){
                    Deg.insert(DegComp(neix[j]));
//					cout<<"insert back/// "<<neix[j]<<endl;
                }
                neix.clear();
                break;///
            }
            else
                x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            if(count==1){
                lastx=x;
                break;
            }else if(NeighborSketchS[lastx].find(x)==NeighborSketchS[lastx].end()){//if not found
                lastx=x;
                break;
            }else{
                Deg.erase(DegComp(x));
                neix.push_back(x);
//				cout<<"erase "<<x<<endl;
            }

//            if(Deg.empty())////
//                break;
        }

        if(neix.size()!=0){
            for(int j=0;j<neix.size();j++){
                Deg.insert(DegComp(neix[j]));
//				cout<<"insert back "<<neix[j]<<endl;
            }
        }
        neix.clear();

        vNodeOrderOverlay.push_back(x);
        Deg.erase(DegComp(x));
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}

void Graph::PartitionOrderingBuildMDE(bool ifParallel){
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    _DD_.assign(node_num,0);
    DD.assign(node_num,0);
    vNodeOrderParti.assign(partiNum,map<int,int>());

    if(ifParallel){
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingV, this, boost::ref(processID[j]) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrdering, this, j));
            }
            thread.join_all();
        }
    }
    else{
        for(int i=0;i<partiNum;++i){
            PartitionOrdering(i);
        }
    }

}

void Graph::PartitionOrderingV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        PartitionOrdering(p[i]);
    }
}

void Graph::PartitionOrdering(int pid){

    set<DegComp> Deg;
    int ID,degree;
    for(int i=0;i<PartiVertex[pid].size();i++){
        ID = PartiVertex[pid][i];
        degree=NeighborsParti[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    while(!Deg.empty()){
//        if(count%10000==0)
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

        vNodeOrderParti[pid].insert({order_i,x});
        order_i++;
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderParti[pid].size() != PartiVertex[pid].size()){
        cout<<"Inconsistent size! "<< pid <<" "<<vNodeOrderParti[pid].size() <<" "<< PartiVertex[pid].size()<<endl; exit(1);
    }
}



/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){//delete u from v's neighbor
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::insertEOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//        DD2[u]++;
    }
}


/// Query Processing
//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... ";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
        int pid=rand()%partiNum;
        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
//        s=115891,t=114656;//NY
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
        }
        d1=Dijkstra(s,t,Neighbor);

        tt.start();
        d2=Query(s,t);
        tt.stop();
        runT+=tt.GetRuntime();
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
            QueryDebug(s,t);
            exit(1);
        }
    }
    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

//function for Query processing, debug version
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in core
        cout<<"Core-Core"<<endl;
//        dis=QueryCoreDebug(ID1, ID2);

    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
//        dis=QueryPartiCoreDebug(ID2, ID1);

    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
//        dis=QueryPartiCoreDebug(ID1, ID2);

    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            cout<<"Parti-Parti"<<endl;
            int d=INF;
            int b1,b2,d1,d2;//final results
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            vector<int> B1=BoundVertex[pid1];
            vector<int> B2=BoundVertex[pid2];

            map<int,int> m1,m2;
            m1.clear();
            m2.clear();
            int bID1, bID2, tempdis;
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];

//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2, QueryH2HPartition(ID2,bID2,pid2)));
            }

            for(int k=0;k<B1.size();k++){
                bID1=B1[k];

                if(m1[bID1]>d)
                    continue;

                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];

                    if(m2[bID2]>d)
                        continue;

                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d){
                        d=tempdis;
                        b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
                    }
                }
            }
            dis=d;
            int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
            cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<") "<<b2<<"("<<NodeOrder[b2]<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

        }
        else{//Case 4: in the same periphery
            cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
            int d=INF;
            int b1,b2,df1,df2;
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            int temp_dis = QueryH2HPartition(ID1, ID2, pid1);/// d2 may be wrong sometimes
            if (temp_dis < d){
                d = temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
                b1=b2=-1;
                df1=df2=-1;
            }

            vector<int> B = BoundVertex[pid1];
            map<int, int> m1, m2;
            m1.clear();
            m2.clear();
            vector<int> B1, B2;
            B1.clear();
            B2.clear();
            int bID, d1, d2;
            for (int i = 0; i < B.size(); i++) {
                bID = B[i];

                d1 = QueryH2HPartition(ID1,bID,pid1);
                d2 = QueryH2HPartition(ID2,bID,pid1);

                if (d1 < d) {
                    B1.push_back(bID);
                    m1.insert(make_pair(bID, d1));
                }
                if (d2 < d) {
                    B2.push_back(bID);
                    m2.insert(make_pair(bID, d2));
                }
            }

            int bID1, bID2, tempdis;
            if (!B1.empty() && !B2.empty()) {
                for (int k = 0; k < B1.size(); k++) {
                    bID1 = B1[k];
                    if (m1[bID1] > d)
                        continue;
                    for (int z = 0; z < B2.size(); z++) {
                        bID2 = B2[z];
                        if (m2[bID2] > d)
                            continue;
                        tempdis = m1[bID1] + QueryCore(bID1, bID2) + m2[bID2];
                        if (tempdis < d){
                            d = tempdis;
                            b1=bID1;b2=bID2;
                            df1=m1[bID1];df2=m2[bID2];
                        }
                    }
                }
            }

            if(b1!=-1){
                cout<<"d4: "<<ID1<<" "<<b1<<" "<<b2<<" "<<ID2<<" : "<<df1<<" "<<QueryCore(b1,b2)<<" "<<df2<<" ; "<<Dijkstra(ID1,b1,Neighbor)<<" "<<Dijkstra(b1,b2,Neighbor)<<" "<<Dijkstra(b2,ID2,Neighbor)<<endl;
            }else{
                int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
                cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
                if(d!=dDijk2){
//                        DijkstraPath(ID1,ID2);
                }
            }

            dis = d;

        }

    }
    return dis;
}

//function for core index correctness check
void Graph::CorrectnessCheckCore(int runtimes){
    srand (time(NULL));
    int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<node_num;++i){
        if(PartiTag[i].second){
            coreVertex.emplace_back(i);
        }
    }
    int corenum=coreVertex.size();
    cout<<"Core graph correctness check ("<<runtimes<<" rounds)..."<<endl;
    for(int i=0;i<runtimes;i++){
        s=coreVertex[rand()%corenum];
        t=coreVertex[rand()%corenum];
        if(PartiTag[s].second && PartiTag[t].second){//for core vertex
            d1=QueryCore(s,t);
            d2=DijkstraCore(s,t);

            if(d1!=d2){
                cout<<"InCorrect! "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<"): "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
                exit(1);
            }
        }else
            i--;
    }
}

//function for efficiency test
void Graph::EffiCheck(string filename,int runtimes){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.push_back(make_pair(ID1, ID2));
    }
    IF.close();
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0;
    int d1,d2;
    clock_t start = clock();
    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
//        if(PartiTag[ID1].first!=PartiTag[ID2].first){
//            cout<<"Different Partition: "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<endl;
//        }
        tt.start();
        d1=Query(ID1,ID2);
        tt.stop();
        runT+=tt.GetRuntime();
        results[i]=d1;
//        d2= Dijkstra(ID1,ID2,Neighbor);
//        if(d1!=d2){
//            cout<<"Incorrect! "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<"): "<<d1<<" "<<d2<<endl; exit(1);
//        }
    }


    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}

void Graph::DFSTree(vector<int>& tNodes, int id){
    tNodes.push_back(Tree[id].uniqueVertex);
    for(int i=0;i<Tree[id].ch.size();++i){
        DFSTree(tNodes,Tree[id].ch[i]);
    }
}





/// Index Maintenance

void Graph::IndexMaintenance(int updateType, bool ifBatch, int batchNumber, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = graphfile + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);
    cout<<"Update batch number: "<<batchNumber<<" ; Batch size: "<<batchSize<<endl;
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);
    Timer tt;
    double runT1=0, runT2 = 0;
    map<pair<int,int>,int> uEdges;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"\nUpdate type: Decrease"<<endl;
            Graph g2=*this;
            if(ifBatch){//for batch update
                if(batchNumber*batchSize>updateData.size()){
                    batchNumber=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNumber;u++){
                    wBatch.clear();
                    uEdges.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;

                        if(ID1>ID2){
                            int temp=ID1;
                            ID1=ID2, ID2=temp;
                        }
                        bool ifFind=false;
                        for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                            if(it->first==ID2){
                                ifFind=true;
                                oldW=it->second;
                                newW=0.5*oldW;
                                break;
                            }
                        }

                        if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                            uEdges.insert({make_pair(ID1,ID2),newW});
                        }else{
                            cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<newW<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                            exit(1);
                        }
                        if(!ifFind){
                            cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
                        }

//                        oldW = updateData[update_i].second;
//                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
//                        if(ifDebug){
//                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
//                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                        ++update_i;
                    }
                    cout<<"Batch "<<u<<". "<<wBatch.size()<<endl;
                    tt.start();
                    g2.DecreaseBatch(wBatch);
//                    DecreaseBatch(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    cout<<"Update time: "<<tt.GetRuntime()<<" s."<<endl;
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
//                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease update Time: "<<runT1/batchNumber<<" s; "<<runT1/(batchNumber*batchSize)<<" s."<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<batchNumber;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    tt.start();
                    g2.DecreaseSingle(ID1,ID2,oldW,newW);
//                    DecreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }

                }

                cout<<"Average Decrease update Time: "<<runT1/batchNumber<<" s."<<endl;
            }

//            break;
        }
        case 2:{
            //Increase update
            cout<<"\nUpdate type: Increase"<<endl;
            if(ifBatch){//for batch update
                if(batchNumber*batchSize>updateData.size()){
                    batchNumber=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNumber;u++){
                    wBatch.clear();
                    uEdges.clear();
                    for(int i=0;i<batchSize;++i){
                        update_i=u*batchSize+i;
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        if(ID1>ID2){
                            int temp=ID1;
                            ID1=ID2, ID2=temp;
                        }
                        bool ifFind=false;
                        for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                            if(it->first==ID2){
                                ifFind=true;
                                oldW=it->second;
                                newW=2*oldW;
                                break;
                            }
                        }

                        if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                            uEdges.insert({make_pair(ID1,ID2),newW});
                        }else{
                            cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<newW<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                            exit(1);
                        }
                        if(!ifFind){
                            cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
                        }

//                        oldW = updateData[update_i].second;
//                        newW=oldW*1.5;
//                        if(ifDebug){
//                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
//                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                        ++update_i;
                    }
                    cout<<"Batch "<<u<<". "<<wBatch.size()<<endl;
                    tt.start();
                    IncreaseBatch(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    cout<<"Update time: "<<tt.GetRuntime()<<" s."<<endl;
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase update Time: "<<runT2/batchNumber<<" s; "<<runT2/(batchNumber*batchSize)<<" s."<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<batchNumber;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<") "<<oldW<<" "<<newW<<endl;
                    }

                    tt.start();
                    IncreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }

                cout<<"Average Increase update Time: "<<runT2/batchNumber<<" s.\n"<<endl;
            }

            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::DecreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    if(!partiBatch.empty()){
//        if(partiBatch.size()>threadnum){
//            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
//        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
//            DecreasePartiBatchUpdateCheck(pid,it->second,overlayBatch);
        }
        thread.join_all();
    }
//    cout<<"update overlay."<<endl;
    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax);

//     repair the partition index
    if(algoQuery==1){
//        cout<<"Repair post index"<<endl;
        Repair_PartiIndex(true, false, partiBatch);
    }
}

void Graph::DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch decrease update
    DecreasePartiBatch(wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                sm->notify();
            }
        }
    }
}

void Graph::IncreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    if(!partiBatch.empty()){
//        if(partiBatch.size()>threadnum){
//            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
//        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        if(algoQuery==1){
//            Trees=TreesNo;
//            for(int pi=0;pi<ifRepaired.size();++pi){
//                if(ifRepaired[pi]){
//                    Trees[pi]=TreesNo[pi];
//                }
//            }
        }
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();

    }

    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);

    // repair the partition index
    if(algoQuery==1){
//        Trees=TreesNo;
//        cout<<"Post-boundary update."<<endl;
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
    }
}

void Graph::IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch Increase update
    IncreasePartiBatch(wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis>olddis){//if '=', not problem; if '<', problem
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
    if(algoQuery==1){
//        TreesNo[pid]=Trees[pid];
    }
}

//Function for single-edge decrease update
void Graph::DecreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbor[b][i].second=newW;
            break;
        }
    }
    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();
    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    if(pid1!=pid2){//for cut edge
        DecreaseOverlay(a, b, newW,NeighborsOverlay,Tree,rank,heightMax);
    }else{//for intra edge
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();
        DecreaseParti(a,b,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1]);

        //weightOverlay collect the changed edges on overlay graph
        vector<int> Bid=BoundVertex[pid1];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }

                newdis=QueryH2HPartition(bid1,bid2,pid1);
                if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                }
            }
        }

        DecreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax);
        //cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        //update the overlay graph index, after partition index update
        /*for(int l=0;l<weightOverlay.size();l++){
            Decrease(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
        }*/
    }
    // repair the partition index
    if(algoQuery==1){
        Repair_PartiIndex(true, false, partiBatch);
    }

}

//Function for single-edge increase update
void Graph::IncreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            Neighbor[b][i].second=newW;
            break;
        }
    }

    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();


//    cout<<"Before update "<<endl;
//    int b1=133362, b2=137086;
//    if(NeighborsOverlay[b1].find(b2) != NeighborsOverlay[b1].end()){//if found
//        cout<<"Overlay dis: "<<NeighborsOverlay[b1][b2]<<endl;
//    }
//    cout<<"Partition dis: "<<QueryH2HPartition(b1,b2,PartiTag[b1].first)<<endl;

    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    //cout<<"increase edge in partition "<<pid<<endl;
    if(pid1!=pid2) {//for cut edge
        cout<<"Inter edge update"<<endl;
        IncreaseOverlay(a, b, oldW, newW,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);
    }
    else {
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
//        cout<<"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"<<endl;
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();
        if(algoQuery==1){
//            Trees=TreesNo;
//            for(int pi=0;pi<ifRepaired.size();++pi){
//                if(ifRepaired[pi]){
//                    Trees[pi]=TreesNo[pi];
//                }
//            }
        }
        IncreaseParti(a,b,oldW,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1],SCconNodesMTP,VidtoTNidP);

        //cout<<"/////////////////////////////////////////"<<endl;

        //cout<<"boundary edge checkkkkkkkkkkkkkkkkkkkkkkkkkkk"<<endl;
        //boundary edges check
        vector<int> Bid=BoundVertex[pid1];
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];//only works for no-boundary
//                    olddis= QueryCore(bid1,bid2);
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }
                newdis=QueryH2HPartition(bid1,bid2,pid1);
//                newdis=ShortcutDisCheck(bid1,bid2);
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                int overlaydis=QueryCore(bid1,bid2);
                if(newdis>olddis)//if '=', not problem; if '<', problem
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                else if(newdis<olddis){
                    cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                    exit(1);
                }


            }
        }

        IncreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);
//        CorrectnessCheckCore(100);
        //update the overlay graph index, after partition index update
        /*cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        for(int l=0;l<weightOverlay.size();l++){
            Increase(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
        }*/
        //cout<<"''''''''''''''''''''''''''''''''''''''''''"<<endl;
    }
    if(algoQuery==1){
//        Trees=TreesNo;
        Repair_PartiIndex(true, true, partiBatch);
    }
}



