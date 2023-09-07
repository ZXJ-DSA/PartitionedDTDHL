/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 11){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> partition number, e.g. 64\n");
        printf("<arg4> (optional) partition method, (NC: PUNCH; MT: METIS), default: NC\n");
        printf("<arg5> (optional) query strategy, (0: No-boundary; 1: Post-boundary), default: 0\n");
        printf("<arg6> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg7> (optional) whether batch update, (0: No; 1: Yes), default: 0\n");
        printf("<arg8> (optional) batch size, default: 10\n");
        printf("<arg9> (optional) thread number, default: 15\n");
        printf("<arg10> (optional) preprocessing task (1: Partitioned MDE Ordering; 2: Partitioned Query Generation)\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int partitionNum = 20;
    int algoQuery = 0;
    int algoUpdate = 0;
    string algoParti = "NC";
    int updateType = 0;
    int runtimes = 10000;
    int updateBatch = 10;
    updateBatch = 1000;
    updateBatch = 100;
    int batchSize = 10;
    bool ifBatch = false;
    int threadNum = 15;
    int preTask=0;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];

        cout << "argv[3] (Partition Number): " << argv[3] << endl;//treewidth
        partitionNum = stoi(argv[3]);

        if(argc > 4){
            cout << "argv[4] (Partition Method): " << argv[4] << endl;//partition method
            algoParti = argv[4];
            if(algoParti != "NC" && algoParti != "MT"){
                cout<<"Wrong partition method! "<<algoParti<<endl; exit(1);
            }
        }
        if(argc > 5){
            cout << "argv[5] (Query Strategy): " << argv[5] << endl;//algorithm for query
            algoQuery = stoi(argv[5]);
        }

        if(argc > 6){
            cout << "argv[6] (Update Type): " << argv[6] << endl;//update type
            updateType = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (Whether Batch Update): " << argv[7] << endl;//algorithm for update
            ifBatch = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8] (Batch Size): " << argv[8] << endl;//algorithm for update
            batchSize = stoi(argv[8]);
        }
        if(argc > 9){
            cout << "argv[9] (Thread Number): " << argv[9] << endl;//thread number
            threadNum = stoi(argv[9]);
        }
        if(argc > 10){
            cout << "argv[10] (Preprocessing Task): " << argv[10] << endl;//preprocessing task
            preTask = stoi(argv[10]);
        }

    }

	//used for running time calculation
    Timer tt0;
    tt0.start();

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";

    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    g.graphfile=graphfile;
    g.ifParallel = true;
    g.dataset=dataset;
    g.algoQuery=algoQuery;
    g.algoUpdate=algoUpdate;
    g.algoParti=algoParti;
    g.partiNum=partitionNum;
    cout<<"Dataset: "<<dataset<<endl;
    if(g.algoQuery==0){
        cout<<"This is test for PH2H !!!!!!!"<<endl;
        cout<<"Query strategy: No-boundary"<<endl;
    }else if(g.algoQuery==1){
        cout<<"This is test for PH2H !!!!!!!"<<endl;
        cout<<"Query strategy: Post-boundary"<<endl;
    }else{
        cout<<"Wrong query strategy! "<<g.algoQuery<<endl; exit(1);
    }
    cout<<"Partition method: "<<g.algoParti<<endl;
    cout<<"Partition number: "<<g.partiNum<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    if(ifBatch){
        cout<<"Test for batch update! Batch size: "<<batchSize<<endl;
        updateBatch=updateBatch/batchSize;
    }else{
        cout<<"Test for single-edge update!"<<endl;
        batchSize=1;
    }

    if(preTask==1){
        g.PH2HVertexOrdering(0);//MDE ordering
//        g.PH2HVertexOrdering(1);//Boundary-first ordering
//        g.PH2HVertexOrdering(2);//Boundary-first MDE ordering
    }
    else if(preTask==2){
        g.QueryGenerationParti(true);//same partition and real-world simulation

    }

//    g.ReadGraph(graphfile);//
//    g.StainingMethod(0);

    ///Task 1: Index construction
    g.PH2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
//    g.CorrectnessCheckCore(100);
    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile+"Parti",runtimes);//query efficiency test
    g.EffiCheck(ODfile+"SameParti",runtimes);
    g.EffiCheck(ODfile+"CrossParti",runtimes);
//    exit(0);
    ///Task 3: Index update
    g.IndexMaintenance(updateType,updateBatch, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
