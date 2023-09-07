/*
 * TreeIndex.cpp
 *
 *  Created on: 16 June 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"

vector<int> NodeOrder_;//nodeID order
vector<int> _DD_;//true degree, temporal degree ,_DD2_
vector<int> _DD2_;//count

//// Index Construction
//Function of constructing tree index for partitions
void Graph::Construct_PartiIndex(bool ifParallel){
    //for H2H update
    SCconNodesMTP.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidP.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    Trees.assign(partiNum,vector<Node>());
    toRMQs.assign(partiNum,vector<int>());
    RMQIndexs.assign(partiNum,vector<vector<int>>());
    ranks.assign(partiNum,vector<int>());
    heightMaxs.assign(partiNum,0);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranks[pid].assign(PartiVertex[pid].size(),-1);
    }

    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
                thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            ConstructPH2H_Parti(pid, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
        }
    }
    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(Trees[i].size());
        aveHeight+=heightMaxs[i];
    }
    cout<<"Partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxs.begin(),heightMaxs.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}

//function of vertex allocation
void Graph::ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID){
//        processID.assign(threadnum, vector<NodeId>());
    int pid=0;
    for(int i=0;i<vertices.size();++i){
        pid=i%processID.size();
        processID[pid].emplace_back(vertices[i]);
    }
}
void Graph::ConstructBoundaryShortcut(int pid){
    //boundary edges
    int ID1,ID2,weight;
    for(int i=0;i<BoundVertex[pid].size();i++){
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID1=BoundVertex[pid][i];
            ID2=BoundVertex[pid][j];
            weight=QueryH2HPartition(ID1,ID2,pid);
            NeighborsOverlay[ID1][ID2]=weight;
            NeighborsOverlay[ID2][ID1]=weight;
        }
    }
}
void Graph::ConstructBoundaryShortcutV(vector<int> & p){
    for(int i=0;i<p.size();++i){
        ConstructBoundaryShortcut(p[i]);
    }
}
void Graph::Construct_OverlayGraph(bool ifParallel){
    if(ifParallel){
        //multiple threads
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
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutV, this, boost::ref(processID[j]) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcut, this, j));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcut(k);
        }
    }


}
//Function of constructing tree index for overlay grpah
void Graph::Construct_OverlayIndex(){
    //Create tree for partition
    H2HCreateTree_Overlay();
    //Create labels for partition
    H2HCreateIndex_Overlay();
}

void Graph::ConstructPartitionPost(bool ifParallel){
    NeighborsPartiPost.assign(node_num,unordered_map<int,int>());

    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostPartiV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostParti, this, j));
            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            ConstructPostParti(k);
        }
    }
}
void Graph::ConstructPostParti(int pid){
    int ID;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }
    int ID1,ID2,weight=-1;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
            assert(NeighborsOverlay[ID1].find(ID2)!=NeighborsOverlay[ID1].end());
//            weight=NeighborsOverlay[ID1][ID2];
//            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
//                NeighborsPartiPost[ID1].insert({ID2,weight});
//                NeighborsPartiPost[ID2].insert({ID1,weight});
//            }
            weight= QueryCore(ID1,ID2);
            NeighborsPartiPost[ID1][ID2]=weight;
            NeighborsPartiPost[ID2][ID1]=weight;
        }
    }
}
void Graph::ConstructPostPartiV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        ConstructPostParti(p[i]);
    }
}
void Graph::ConstructPartitionPostIndex(bool ifParallel){
    SCconNodesMTPost.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidPost.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    TreesPost.assign(partiNum,vector<Node>());
    toRMQsPost.assign(partiNum,vector<int>());
    RMQIndexsPost.assign(partiNum,vector<vector<int>>());
    ranksPost.assign(partiNum,vector<int>());
    heightMaxsPost.assign(partiNum,0);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsPartiPost.size();i++){
        for(auto it=NeighborsPartiPost[i].begin();it!=NeighborsPartiPost[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranksPost[pid].assign(PartiVertex[pid].size(),-1);
    }

    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
                thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost)));
            }
            thread.join_all();
        }
    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            ConstructPH2H_Parti(pid,TreesPost, ranksPost, SCconNodesMTPost, VidtoTNidPost, heightMaxsPost, toRMQsPost, RMQIndexsPost);
        }
    }
    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(TreesPost[i].size());
        aveHeight+=heightMaxsPost[i];
    }
    cout<<"Post-boundary partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxsPost.begin(),heightMaxsPost.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}


//Function of repair the partition index
void Graph::Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    repairShortcuts.assign(node_num, unordered_map<vertex,pair<int,int>>());
    ifRepaired.assign(partiNum, false);
    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch)));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncrease, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecrease, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }

        }
    }
    else{
        // single thread
        if(ifIncrease){
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexIncrease(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecrease(k,partiBatch);
            }
        }

    }

//    int pNum=0;
//    for(auto it=ifRepaired.begin();it!=ifRepaired.end();++it){
//        if(*it){
//            ++pNum;
//        }
//    }
//    cout<<"Repaired partition number: "<<pNum<<endl;

}

//Function of repair the partition index

void Graph::RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    if(ifIncrease){
        for(int i=0;i<p.size();++i){
            RepairPartitionIndexIncrease(p[i], partiBatch);
        }
    }
    else{
        for(int i=0;i<p.size();++i){
            RepairPartitionIndexDecrease(p[i], partiBatch);
        }
    }

}
void Graph::RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
        weightsParti=partiBatch[pid];
    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];

            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(woverlay<wlocal){
                if(ID1<ID2){
                    weightsParti.emplace_back(make_pair(ID1,ID2),make_pair(wlocal,woverlay));
                }else{
                    weightsParti.emplace_back(make_pair(ID2,ID1),make_pair(wlocal,woverlay));
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay>wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        DecreasePartiBatch(weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid]);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    map<pair<int,int>,pair<int,int>> updatesSet;

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
//        cout<<"Find for "<<pid<<" "<<partiBatch[pid].size()<<endl;
//        weightsParti=partiBatch[pid];
        for(auto it=partiBatch[pid].begin();it!=partiBatch[pid].end();++it){
            ID1=it->first.first, ID2=it->first.second;
            if(ID1<ID2){
                if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
            }else{
                if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
            }
        }
    }
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(woverlay>wlocal){
                // update partition index
                if(ID1<ID2){
                    if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID1,ID2),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl; exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
                }else{
                    if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID2,ID1),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl; exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay<wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    for(auto it=updatesSet.begin();it!=updatesSet.end();++it){
        weightsParti.emplace_back(*it);
    }


    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        IncreasePartiBatch(weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid],SCconNodesMTPost,VidtoTNidPost);
        ifRepaired[pid]=true;
    }

}

void Graph::ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    int PID;
    for(int i=0;i<P.size();++i){
        PID=P[i];
        ConstructPH2H_Parti(PID, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
    }
}
void Graph::ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    //Create tree for partition
    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}

//Function of Creating tree for partition
void Graph::H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs) {
//    cout<<"Partition "<<pid<<": boundary number is "<<BoundVertex[pid].size()<<endl;
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    unordered_map<vertex,bool> existCore; existCore.clear();
    for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
        existCore.insert({*it,true});
    }

//    unordered_map<vertex,int> rankP; rankP.clear();
    for(int id=PartiVertex[pid].size()-1;id>=0;--id){
        ID = PartiVertex[pid][id];
//        rankP.insert({ID,-1});
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore.find(it->first)==existCore.end()){// not found
                cout<<"Wrong neighbor! "<<ID<<"("<<PartiTag[ID].first<<") "<< it->first<<"("<<PartiTag[it->first].first<<")"<<endl; exit(1);
            }
            else{// if found
                if(existCore[it->first]){
                    Neigh.emplace_back(*it);
                }else{
                    cout<<"Not in core!"<<it->first<<endl; exit(1);
                }
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());
//        cout<<ID<<" "<<NeighborCon[ID].size()<<" "<<PartiTag[ID].second<<endl;

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();j++){
                ID2=Neigh[j].first;
                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                /// For TD update
                if(ID1<ID2){
                    SCconNodesMTP[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                }
                else{
                    SCconNodesMTP[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                }

            }
        }


    }

    /// Create Tree
    ID=PartiVertex[pid][0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    TreeP.push_back(root);
//    cout<<"0 "<<ID<<" "<<IDMap[ID]<<endl;
    rankP[IDMap[ID]] = 0;
//    rankP[ID] = 0;

    for(int id=1;id<PartiVertex[pid].size();++id){
        ID = PartiVertex[pid][id];
//        cout<<id<<" "<<ID<<" "<<IDMap[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCoreParti(ID,NeighborCon[ID], rankP);

//        cout<<"pa "<<pa<<" "<<TreeP[pa].height<<endl;

        TreeP[pa].ch.push_back(TreeP.size());
        nod.pa=pa;
        nod.height=TreeP[pa].height+1;
        /// for update
        nod.hdepth=TreeP[pa].height+1;
        for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[ID][i].first;
            if(PartiTag[nn].first != pid){
                cout<<"Wrong nn! "<<PartiTag[nn].first <<" "<< pid<<endl; exit(1);
            }
            VidtoTNidP[nn].emplace_back(TreeP.size());
//            if(rankP.find(nn)==rankP.end()){
//                cout<<"Not found "<<nn<<" in rank!"<<endl; exit(1);
//            }
            if(TreeP[rankP[IDMap[nn]]].hdepth<TreeP[pa].height+1){
                TreeP[rankP[IDMap[nn]]].hdepth=TreeP[pa].height+1;
            }

        }
        if(nod.height>heightMaxs[pid]){
            heightMaxs[pid]=nod.height;
        }

//        if(rankP.find(ID)==rankP.end()){
//            cout<<"Not found "<<ID<<" in rank!"<<endl; exit(1);
//        }
        rankP[IDMap[ID]]=TreeP.size();//the position of tree, higher-order vertex has lower rank
        TreeP.push_back(nod);

    }

//    cout<<pid<<"'s tree node number: "<<TreeP.size()<<" ; tree height: "<< heightMaxs[pid]<<endl;
}
//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP){
//    cout<<"Computing Tree Label for partition "<<pid<<endl;
    //initialize
    vector<int> list; //list.clear();
    list.push_back(TreeP[0].uniqueVertex);
    TreeP[0].pos.clear();
    TreeP[0].pos.push_back(0);

    for(int i=0;i<TreeP[0].ch.size();i++){
        makeTreeIndexDFSP(TreeP[0].ch[i],list,TreeP, rankP);
    }

}

//Function of Creating tree for partition
void Graph::H2HCreateTree_Overlay() {
    //for H2H update
    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);
    VidtoTNid.assign(node_num,vector<int>());
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsOverlay.size();i++){
        if(!NeighborsOverlay[i].empty()){
            for(auto it=NeighborsOverlay[i].begin();it!=NeighborsOverlay[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    vector<bool> existCore(node_num,true);
    rank.assign(node_num,-1);

    bool flagAdd = true;
    for(int id=OverlayVertex.size()-1;id>=0;--id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore[it->first]){
                Neigh.emplace_back(*it);
            }else{
                cout<<"Not in core!"<<it->first<<endl; exit(1);
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }

        if(Neigh.size()<=100){
            //single thread
            for(int i=0;i<Neigh.size();i++){
                ID1=Neigh[i].first;
                for(int j=i+1;j<Neigh.size();j++){
                    ID2=Neigh[j].first;
                    insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                    /// For TD update
                    if(ID1<ID2){
                        SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                    }else{
                        SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                    }
                }
            }
        }else{
//            cout<<"Multiple thread for contraction. "<<ID<<" "<<Neigh.size()<<endl;
            //multiple thread
            if(Neigh.size()>threadnum){
                int step=Neigh.size()/threadnum;
                boost::thread_group thread;
                for(int i=0;i<threadnum;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==threadnum-1)
                        p.second=Neigh.size();
                    else
                        p.second=(i+1)*step;
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(int i=0;i<Neigh.size();i++){
                    pair<int,int> p;
                    p.first=i; p.second=(i+1);
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }
        }

        //add all-pair neighbors
//        for(int i=0;i<Neigh.size();i++){
//            ID1=Neigh[i].first;
//            for(int j=i+1;j<Neigh.size();j++){
//                ID2=Neigh[j].first;
//                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
//                /// For TD update
//                SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//                SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
//            }
//        }


    }
//    cout<<"Flag 1"<<endl;
    /// Create Tree
    ID=OverlayVertex[0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    Tree.push_back(root);
//    rankP[IDMap[ID]] = 0;
    rank[ID] = 0;
    heightMax=0;

    for(int id=1;id<OverlayVertex.size();++id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCore(ID,NeighborCon[ID], rank);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[ID][i].first;
            VidtoTNid[nn].emplace_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1){
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
            }

        }
        if(nod.height>heightMax){
            heightMax=nod.height;
        }

        rank[ID]=Tree.size();//the position of tree, higher-order vertex has lower rank
        Tree.push_back(nod);

    }

    /// LCA index
    makeRMQCore();//build LCA index

    cout<<"Overlay graph! Tree node number: "<<Tree.size()<<" ; Tree height: "<< heightMax<<endl;
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertECoreMT(ID1,ID2,w1+w2);
            /// For TD update
            if(ID1<ID2){
                SCconNodesMT[ID1][ID2].emplace_back(x,w1+w2);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
            }
        }
    }
//    sm->notify();
}

//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Overlay(){
    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}

//function of computing the H2H label of peripheries: original version
void Graph::makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP){
//initialize
//    cout<<"Map "<<p<<" "<<IDMap[p]<<endl;
//    p=IDMap[p];
    int NeiNum=TreeP[p].vert.size();
    TreeP[p].pos.assign(NeiNum+1,0);
    TreeP[p].dis.assign(list.size(),INF);
    TreeP[p].cnt.assign(list.size(),0);
    TreeP[p].FN.assign(list.size(),true);

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(TreeP[p].vert[i].first==list[j]){
                TreeP[p].pos[i]=j;//record the position of neighbors
                TreeP[p].dis[j]=TreeP[p].vert[i].second.first;
                TreeP[p].cnt[j]=1;
                break;
            }
        }
    }
    TreeP[p].pos[NeiNum]=list.size();


    //dis
    for(int i=0;i<NeiNum;i++){
        int x=TreeP[p].vert[i].first;
        int disvb=TreeP[p].vert[i].second.first;
        int k=TreeP[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=TreeP[rankP[IDMap[y]]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=TreeP[rankP[IDMap[x]]].dis[j];

                if(TreeP[p].dis[j]>z+disvb){
                    TreeP[p].dis[j]=z+disvb;
                    TreeP[p].FN[j]=false;
                    TreeP[p].cnt[j]=1;
                }else if(TreeP[p].dis[j]==z+disvb){
                    TreeP[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(TreeP[p].uniqueVertex);
    for(int i=0;i<TreeP[p].ch.size();i++){
        makeTreeIndexDFSP(TreeP[p].ch[i],list, TreeP, rankP);
    }
    list.pop_back();
}

void Graph::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFS(Tree[p].ch[i],list);
    }
    list.pop_back();
}

//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteECore(int u,int v){
//	if(Emap[u].find(v)!=Emap[u].end()){
//		Emap[u].erase(Emap[u].find(v));
//		DD[u]--;
//	}

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}
//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//		DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}
void Graph::insertECoreMT(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

}

//compute the father tree node
int Graph::matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
    return rank[nearest];
}
//compute the father tree node
int Graph::matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
//    if(rank.find(nearest)==rank.end()){//if not found
//        cout<<"Not found "<<nearest<<" in rank!"<<endl; exit(1);
//    }
    return rank[IDMap[nearest]];
}

//construct RMQ index
void Graph::makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees){
    //EulerSeq.clear();
    toRMQs[pid].assign(PartiVertex[pid].size(),0);
    vector<int> EulerSeqP;
    //RMQIndex.clear();
    makeRMQDFSCoreP(pid, 0, 1, EulerSeqP, toRMQs, Trees);
    RMQIndexs[pid].push_back(EulerSeqP);

    int m = EulerSeqP.size();
//    cout<<"m: "<<m<<endl;
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndexs[pid][k - 1][j], y = RMQIndexs[pid][k - 1][j + i / 2];
//            cout<<"x and y: "<<x<<" "<<y<<endl;
            if (Trees[pid][x].height < Trees[pid][y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndexs[pid].push_back(tmp);
    }
}

void Graph::makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees){
    toRMQs[pid][p] = EulerSeqP.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeqP.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Trees[pid][p].ch.size(); i++){
        makeRMQDFSCoreP(pid,Trees[pid][p].ch[i], height + 1, EulerSeqP, toRMQs, Trees);
        EulerSeqP.push_back(p);
    }
}

//construct RMQ index
void Graph::makeRMQCore(){
    vector<int> EulerSeq;
    EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFSCore(0, 1, EulerSeq);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFSCore(int p, int height, vector<int>& EulerSeq){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFSCore(Tree[p].ch[i], height + 1, EulerSeq);
        EulerSeq.push_back(p);
    }
}


/// Query Processing
//function for Query processing
int Graph::Query(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//        cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCore(ID2, ID1);
    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCore(ID1, ID2);
    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            dis= QueryPartiParti(ID1,ID2);
        }else{//Case 4: in the same periphery
//            cout<<"Same partition!"<<endl;
            if(algoQuery==NO_Boundary){//no-boundary
                dis= QuerySameParti(ID1,ID2);
            }else if(algoQuery==Post_Boundary){
                dis= QuerySamePartiPost(ID1,ID2);
            }


        }
    }
    return dis;
}
//Case 1: query on overlay graph
int Graph::QueryCore(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);

    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}

//Case 2: one core, one tree
int Graph::QueryPartiCore(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    if(algoQuery==NO_Boundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartition(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }
    else if(algoQuery==Post_Boundary){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartitionPost(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }

    return d;
}

//Case 3: Different trees
int Graph::QueryPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;

        if(algoQuery==NO_Boundary){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2,QueryH2HPartition(ID2,bID2,pid2)));
            }
        }else if(algoQuery==Post_Boundary){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartitionPost(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2,QueryH2HPartitionPost(ID2,bID2,pid2)));
            }
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
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}

//Case 4: Same tree, for original version
int Graph::QuerySameParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryH2HPartition(ID1,ID2,pid1);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<B.size();i++){
            bID=B[i];
//            d1=Tree[rank[ID1]].disInf[i];
//            d2=Tree[rank[ID2]].disInf[i];
            d1=QueryH2HPartition(ID1,bID,pid1);
            d2=QueryH2HPartition(ID2,bID,pid1);

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//Case 4: Same tree, for query-orient version
int Graph::QuerySamePartiPost(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartitionPost(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int r1=ranks[PID][IDMap[ID1]], r2=ranks[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartition(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1)
        return Trees[PID][r2].dis[Trees[PID][r1].pos.back()];
    else if(LCA==r2)
        return Trees[PID][r1].dis[Trees[PID][r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Trees[PID][LCA].pos.size();i++){
            if(tmp>Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]])
                tmp=Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]];
        }
        return tmp;
    }
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartitionPost(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int r1=ranksPost[PID][IDMap[ID1]], r2=ranksPost[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartitionPost(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1)
        return TreesPost[PID][r2].dis[TreesPost[PID][r1].pos.back()];
    else if(LCA==r2)
        return TreesPost[PID][r1].dis[TreesPost[PID][r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<TreesPost[PID][LCA].pos.size();i++){
            if(tmp>TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]])
                tmp=TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]];
        }
        return tmp;
    }
}

int Graph::LCAQueryPartition(int _p, int _q, int PID){
    int p = toRMQs[PID][_p], q = toRMQs[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Trees[PID][RMQIndexs[PID][k][p]].height < Trees[PID][RMQIndexs[PID][k][q]].height)
        return RMQIndexs[PID][k][p];
    else return RMQIndexs[PID][k][q];
}

int Graph::LCAQueryPartitionPost(int _p, int _q, int PID){
    int p = toRMQsPost[PID][_p], q = toRMQsPost[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (TreesPost[PID][RMQIndexsPost[PID][k][p]].height < TreesPost[PID][RMQIndexsPost[PID][k][q]].height)
        return RMQIndexsPost[PID][k][p];
    else return RMQIndexsPost[PID][k][q];
}

int Graph::LCAQueryOverlay(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}


/// Index maintenance
//H2H index update
void Graph::DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;//map<tree node ID, distance index>

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first>newW){
                Tree[rank[lid]].vert[i].second.first=newW;
                Tree[rank[lid]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lid]].vert[i].second.first==newW){
                Tree[rank[lid]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
        //cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }else{
                        Cw=Vert[j].second.first;
                    }
                }

                if(Tree[rank[ProID]].dis[cidH]>=Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Tree[rank[Cid]].height].insert(hid);
                            if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;

                        }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lid]].vert[k].second.first>wsum){
                                Tree[rank[lid]].vert[k].second.first=wsum;
                                Tree[rank[lid]].vert[k].second.second=1;
                                SCre[Tree[rank[lid]].height].insert(Cid);
                                if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;

                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL, checkedDis, Tree, rank);
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis,Tree, rank);
    }
    line.pop_back();

}

void Graph::DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;//map<tree node ID, distance index>

    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first>newW){
                Tree[rank[lidM]].vert[i].second.first=newW;
                Tree[rank[lidM]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                Tree[rank[lidM]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
//        cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int CidM=IDMap[Cid];
                int cidH=Tree[rank[CidM]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }else{
                        Cw=Vert[j].second.first;
                    }
                }

                if(Tree[rank[ProIDM]].dis[cidH]>=Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                    hid=Tree[rank[CidM]].vert[j].first;hidHeight=Tree[rank[IDMap[hid]]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.first=wsum;
                            Tree[rank[CidM]].vert[j].second.second=1;
                            SCre[Tree[rank[CidM]].height].insert(hid);
                            if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;

                        }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    int lidM=IDMap[lid];
                    lidHeight=Tree[rank[lidM]].height-1;
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lidM]].vert[k].second.first>wsum){
                                Tree[rank[lidM]].vert[k].second.first=wsum;
                                Tree[rank[lidM]].vert[k].second.second=1;
                                SCre[Tree[rank[lidM]].height].insert(Cid);
                                if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;

                            }else if(Tree[rank[lidM]].vert[k].second.first==wsum){
                                Tree[rank[lidM]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5Parti(rank[IDMap[ProBeginID]], linee, vertexIDChL, checkedDis, Tree, rank);
    }
    else{
//        cout<<"Not trigger update! "<<lid<<" "<<hid<<endl;
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5Parti(Tree[child].ch[i], line, vertexIDChL,checkedDis,Tree, rank);
    }
    line.pop_back();

}

//batch update for overlay graph
void Graph::DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        if(Neighbor[a].find(b)!=Neighbor[a].end()){
            Neighbor[a][b]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }
        if(Neighbor[b].find(a)!=Neighbor[b].end()){
            Neighbor[b][a]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProID]].dis[cidH]>Cw){
                Tree[rank[ProID]].dis[cidH]=Cw;
                Tree[rank[ProID]].FN[cidH]=true;
                ProIDdisCha=true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                Tree[rank[ProID]].FN[cidH]=true;
            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis,Tree,rank);
    }
    //return checkedDis.size();
}

//batch update for partition graph
void Graph::DecreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
                Neighbors[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
                Neighbors[b][i].second=newW;
                break;
            }
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                Tree[rank[ProIDM]].dis[cidH]=Cw;
                Tree[rank[ProIDM]].FN[cidH]=true;
                ProIDdisCha=true;
                Tree[rank[ProIDM]].DisRe.insert(Cid);
            }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                Tree[rank[ProIDM]].FN[cidH]=true;
            }

            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompp(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompp(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            Tree[rank[lid2M]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[IDMap[ProBeginVertexSet[i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChL,checkedDis,Tree,rank);
    }
    //return checkedDis.size();
}

//batch update for partition graph
void Graph::DecreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        if(Neighbors[a].find(b)!=Neighbors[a].end()){
            Neighbors[a][b]=newW;
        }else{
            cout<<"Not found edge! "<<endl; exit(1);
        }

        if(Neighbors[b].find(a)!=Neighbors[b].end()){
            Neighbors[b][a]=newW;
        }else{
            cout<<"Not found edge! "<<endl; exit(1);
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                Tree[rank[ProIDM]].dis[cidH]=Cw;
                Tree[rank[ProIDM]].FN[cidH]=true;
                ProIDdisCha=true;
                Tree[rank[ProIDM]].DisRe.insert(Cid);
            }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                Tree[rank[ProIDM]].FN[cidH]=true;
            }

            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompp(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompp(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            Tree[rank[lid2M]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[IDMap[ProBeginVertexSet[i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChL,checkedDis,Tree,rank);
    }
    //return checkedDis.size();
}

void Graph::IncreaseOverlay(int a,int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[pachid]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first==oldW){
                Tree[rank[lid]].vert[i].second.second-=1;
                if(Tree[rank[lid]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second-=1;
                            if(Tree[rank[Cid]].vert[j].second.second<1){
                                SCre[Tree[rank[Cid]].height].insert(hid);
                                if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lid]].vert[k].second.second-=1;
                                if(Tree[rank[lid]].vert[k].second.second<1){
                                    SCre[Tree[rank[lid]].height].insert(Cid);
                                    if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                    if(it2->first==Cid){
                        Cw=it2->second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                    if(Tree[rank[ProID]].vert[i].first==Cid){
                        Tree[rank[ProID]].vert[i].second.first=Cw;
                        Tree[rank[ProID]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
    while(Tree[rank[pachid]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]==0){
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
                if(Tree[PID].height>Tree[children].height){///
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::IncreaseParti(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first==oldW){
                Tree[rank[lidM]].vert[i].second.second-=1;
                if(Tree[rank[lidM]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[IDMap[Cid]]].height-1;
                int CidM=IDMap[Cid];

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){// for higher-order vertices
                    hid=Tree[rank[CidM]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second-=1;
                            if(Tree[rank[CidM]].vert[j].second.second<1){
                                SCre[Tree[rank[CidM]].height].insert(hid);
                                if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){//for lower-order vertices
                    lid=Lnei[j].first; lidM=IDMap[lid];
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            if(Tree[rank[lidM]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lidM]].vert[k].second.second-=1;
                                if(Tree[rank[lidM]].vert[k].second.second<1){
                                    SCre[Tree[rank[lidM]].height].insert(Cid);
                                    if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProIDM]].FN[cidH]){//if the label is obtained from shortcut
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[Cid]]].dis[i]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProIDM]].FN[cidH]=false;
                    Tree[rank[ProIDM]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(int i=0;i<Neighbors[ProID].size();i++){
                    if(Neighbors[ProID][i].first==Cid){
                        Cw=Neighbors[ProID][i].second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){//for each supportive vertex of this shortcut
                    wid=Wnodes[i].first;
                    int widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                    if(Tree[rank[ProIDM]].vert[i].first==Cid){
                        Tree[rank[ProIDM]].vert[i].second.first=Cw;
                        Tree[rank[ProIDM]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginID]],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
//    int childM=IDMap[children];
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]==0){//the distance from child to line[i] should be updated
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int tid;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//check the tree nodes containing child
                tid=VidtoTNid[childID][k];
//                cout<<k<<" "<<tid<<" "<<IDMap[tid]<<" "<<childH<<" "<<Tree[tid].dis.size()<<" "<<Tree[tid].FN.size()<<endl;
                if(Tree[tid].FN[childH] && Tree[tid].dis[i]==disBF+Tree[tid].dis[childH]){//if the distance from tid to line[i] sources from child, valley path
                    Tree[tid].cnt[i]-=1;
                }
            }
//            cout<<"Flag 2"<<endl;
            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                tid=VidtoTNid[line[i]][k];
                if(Tree[tid].height>Tree[children].height){//children is the ancestor of tid
                    if(Tree[tid].FN[i] && Tree[tid].dis[childH]==disBF+Tree[tid].dis[i]){//if the distance from tid to child sources from line[i], peak path, out of scope
                        Tree[tid].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[IDMap[b]]].height-1;
                if(bH<i){//if b is the ancestor of line[i]
                    if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[IDMap[line[i]]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;//shortcut
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{//if line[i] is the ancestor of b
                    if(Dvb+Tree[rank[IDMap[b]]].dis[i]<dis){
                        dis=Dvb+Tree[rank[IDMap[b]]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[b]]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1Parti(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbor[a].find(b)!=Neighbor[a].end()){
                if(Neighbor[a][b]!=oldW){//only works for no-boundary
                    cout<<"Inconsistent! "<<Neighbor[a][b]<<" "<<oldW<<endl; exit(1);
                }
                Neighbor[a][b]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }
            if(Neighbor[b].find(a)!=Neighbor[b].end()){
                if(Neighbor[b][a]!=oldW){
                    cout<<"Inconsistent! "<<Neighbor[b][a]<<" "<<oldW<<endl; exit(1);
                }
                Neighbor[b][a]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompp(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid]].vert[k].second.second-=1;
                            if(Tree[rank[lid]].vert[k].second.second<1){
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompp(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProID]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                        Tree[rank[ProID]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProID]].FN[cidH]=false;
                Tree[rank[ProID]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                    if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                        Tree[rank[ProID]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbor[ProID].begin();it2!=Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=Cw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }

    }

    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
    }
    //return checknum;
}

void Graph::IncreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    Neighbors[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    Neighbors[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompp(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid2M]].vert[k].second.second-=1;
                            if(Tree[rank[lid2M]].vert[k].second.second<1){
                                SCre[lid2].insert(Cid);
                                OC.insert(OrderCompp(lid2));
                                OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProIDM]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProIDM]].FN[cidH]=false;
                Tree[rank[ProIDM]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=Cw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[IDMap[ProBeginVertexSet[i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }

    }

    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
    }
    //return checknum;
}

void Graph::IncreasePartiBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbors[a].find(b)!=Neighbors[a].end()){
                Neighbors[a][b]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            if(Neighbors[b].find(a)!=Neighbors[b].end()){
                Neighbors[b][a]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompp(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid2M]].vert[k].second.second-=1;
                            if(Tree[rank[lid2M]].vert[k].second.second<1){
                                SCre[lid2].insert(Cid);
                                OC.insert(OrderCompp(lid2));
                                OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProIDM]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProIDM]].FN[cidH]=false;
                Tree[rank[ProIDM]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=Cw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[IDMap[ProBeginVertexSet[i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }

    }

    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
    }
    //return checknum;
}

int Graph::ShortcutDisCheck(int ID1, int ID2){
    int d=INF;
    int ID;
    if(PartiTag[ID1].first!=PartiTag[ID2].first){
        cout<<ID1<<" and "<<ID2<<" are not in the same partition!"<<endl; exit(1);
    }
    int PID=PartiTag[ID1].first;
    int lid,hid;
    if(NodeOrder[ID1]>NodeOrder[ID2]){
        hid=ID1, lid=ID2;
    }else{
        hid=ID2, lid=ID1;
    }
    bool flag=false;
/*    for(int i=0;i<Trees[PID][ranks[PID][IDMap[lid]]].vert.size();i++){
        if(Trees[PID][ranks[PID][IDMap[lid]]].vert[i].first==hid){
            d=Trees[PID][ranks[PID][IDMap[lid]]].vert[i].second.first;
            flag=true;
            break;
        }
    }*/

    int Cw=INF;
    for(auto it2=NeighborsParti[lid].begin();it2!=NeighborsParti[lid].end();++it2){
        if(it2->first==hid){
            Cw=it2->second;//the weight value in the original graph
//            cout<<lid<<" "<<hid<<", Cw1: "<<Cw<<endl;
            flag=true;
            break;
        }
    }

    if(SCconNodesMTP[lid].find(hid) != SCconNodesMTP[lid].end()){//if found
        int widM;
        for(auto it=SCconNodesMTP[lid][hid].begin();it!=SCconNodesMTP[lid][hid].end();++it){
            ID=it->first;
            if(PartiTag[ID].second){//if boundary vertex
//                cout<<"Continue "<<ID<<endl;
                continue;
            }
            widM=IDMap[ID];
            int ssw=-1,wtt=-1;
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++){
                if(Trees[PID][ranks[PID][widM]].vert[j].first==lid){
                    ssw=Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }

            }
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++) {
                if (Trees[PID][ranks[PID][widM]].vert[j].first == hid) {
                    wtt = Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }
            }
            if(ssw==-1 || wtt==-1){
                cout<<"Wrong! "<<ssw<<" "<<wtt<<endl; exit(1);
            }

            if(ssw+wtt<Cw){
                Cw=ssw+wtt;
//                cout<<lid<<" "<<hid<<", Cw2: "<<Cw<<endl;
            }
        }
        //cout<<Cw<<endl;
        flag=true;
    }
//    else{
//        cout<<"Shortcut sc("<<ID1<<","<<ID2<<") does not exist!"<<endl;
//    }
    d=Cw;
    if(!flag){
        cout<<"Not found! "<<endl; exit(1);
    }
    return d;
}
