/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "headH2H.h"

int main(int argc, char** argv){

    if( argc < 3 || argc > 9){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg4> (optional) whether batch update, (0: No; 1: Yes), default: 0\n");
        printf("<arg5> (optional) batch number, default: 10\n");
        printf("<arg6> (optional) batch size, default: 100\n");
        printf("<arg7> (optional) thread number, default: 15\n");
        printf("<arg8> (optional) query path\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";

    int updateType = 0;
    bool ifBatch = false;
    int runtimes = 10000;
    int batchNumber = 10;
    int batchSize = 100;
    int threadNum = 15;
    string queryPath;

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
            cout << "argv[5] (Batch Number): " << argv[5] << endl;//batch number
            batchNumber = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (Batch Size): " << argv[6] << endl;//batch size
            batchSize = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (Thread Number): " << argv[7] << endl;//thread number
            threadNum = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8] (Query Path): " << argv[8] << endl;//Query path
            queryPath = argv[8];
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
    g.EffiCheckH2H(ODfile,runtimes);//query efficiency test
    g.EffiCheckH2H(queryPath+"/sameParti.query",runtimes);//query efficiency test
    g.EffiCheckH2H(queryPath+"/crossParti.query",runtimes);//query efficiency test
    g.EffiCheckH2H(queryPath+"/mixParti.query",runtimes);
//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);

    ///Task 3: Index update
    g.IndexMaintenanceH2H(updateType, ifBatch, batchNumber, batchSize);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
