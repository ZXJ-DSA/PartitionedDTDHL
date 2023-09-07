## Introduction
This is the source code of the TKDE paper "*Partitioned Dynamic Hub Labeling for Large Road Networks*" (submitted). Please refer to the paper for the algorithm details.

## Algorithms

The implementation code includes the index construction, query processing, and index update of our *DTDHL* algorithm. The runnable components include:

1. Dynamic Hub Labeling *DTDHL*: the entry is `mainH2H.cpp`
1. Partitioined Dynamic Hub Labeling *DTDHL-P* and *DTDHL-PR*: then entry is `main.cpp`





## Data
The datasets of this paper are sourced from [http://www.diag.uniroma1.it/challenge9/download.shtml](http://www.diag.uniroma1.it/challenge9/download.shtml). Please refer to the paper for details.

An example graph *NY* is provided in the directory *data* for your reference. You can run *DPH2H* and *H2H* on the example graph by using the source path `./data`. 

## Dependency

1. `g++` and `boost`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
