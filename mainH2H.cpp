/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "headH2H.h"

int main(int argc, char** argv){

    if( argc < 3 || argc > 7){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg4> (optional) whether batch update, (0: No; 1: Yes), default: 0\n");
        printf("<arg5> (optional) batch size, default: 10\n");
        printf("<arg6> (optional) thread number, default: 15\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";

    int updateType = 0;
    bool ifBatch = false;
    int runtimes = 10000;
    int updateBatch = 10;
    updateBatch = 100;
    int batchSize = 10;
    int threadNum = 15;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];


        if(argc > 3){
            cout << "argv[3] (Update Type): " << argv[3] << endl;//update type
            updateType = stoi(argv[3]);
        }
        if(argc > 4){
            cout << "argv[4] (Batch Update): " << argv[4] << endl;//batch update
            ifBatch = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5] (Batch Size): " << argv[5] << endl;//batch size
            batchSize = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (Thread Number): " << argv[6] << endl;//thread number
            threadNum = stoi(argv[6]);
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
    g.dataset=dataset;
    g.threadnum=threadNum;
    g.graphfile=graphfile;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Thread number: "<<threadNum<<endl;
    if(ifBatch){
        cout<<"Test for batch update! Batch size: "<<batchSize<<endl;
        updateBatch=updateBatch/batchSize;
    }else{
        cout<<"Test for single edge update!"<<endl;
        batchSize=1;
    }

    cout<<"This is test for H2H !!!!!!!"<<endl;

    g.ReadGraph(graphfile);//
//    g.StainingMethod(0);

    ///Task 1: Index construction
    g.H2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
//    g.CorrectnessCheckCore(100);
    g.CorrectnessCheckH2H(100);
    g.EffiCheckH2H(ODfile+"Parti",runtimes);//query efficiency test
    g.EffiCheckH2H(ODfile+"SameParti",runtimes);//query efficiency test
    g.EffiCheckH2H(ODfile+"CrossParti",runtimes);
//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);

    ///Task 3: Index update
    g.IndexMaintenanceH2H(updateType, updateBatch, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
