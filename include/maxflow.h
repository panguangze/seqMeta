//
// Created by wangxuedong on 11/6/23.
//

#ifndef SEQGRAPH_MAXFLOW_H
#define SEQGRAPH_MAXFLOW_H
#include "Graph.h"
#include "algorithm"
#include <deque>
#include "util.h"
extern int VERBOSE;
extern bool BREAK_C;
class maxFlow {
private:
    seqGraph::Graph* graph;
    uint V;
    std::vector<std::vector<int>> augmentationPaths;
    float ** ConjugateMatrix;
    int sId;
    int tId;
    std::map<uint,std::string> matrixIdx2nodeName;
public:
    seqGraph::Graph *getGraph() const;

private:
    float** currentMatrix;
    float** originalMatrix;
    std::vector<seqGraph::Vertex*>* originalVertices;
    seqGraph::Graph* originalGraph;
    int N;


public:
    explicit maxFlow(seqGraph::Graph* graph1);
    ~maxFlow();

    void printM(int i);

    static bool cmpVertex (int i,int j);
    inline int getN() const {
        return N;
    }

    inline float** getMatrix() const {
        return this->currentMatrix;
    };
//    inline seqGraph::SparseMatrix& getMatrix() const {
//        return this->graph->getConjugateMatrix();
//    }

    std::map<int, std::vector<int>*>* resolvePath(std::map<int, std::vector<int>*>* prevPaths);





    bool vertexLookup(int i, int j);
    int  inDegree(int idx);
    int  outDegree(int idx);
    bool bfs(float **rGraph, int s, int t, int parent[]);
    void dfs(float **rGraph, int s, int t, bool visited[], std::vector<int>& path);
    void fordFulkerson();
};


#endif SEQGRAPH_MAXFLOW_H

