//
// Created by wangxuedong on 11/6/23.
//

//
// Created by caronkey on 28/12/2021.
//

#include "../include/maxflow.h"
#include <cstring>
#include <algorithm>
#include <queue>
#include <iostream>
#include "limits"
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <limits.h>
#include <cfloat>

const float ZERO = 0.000001;
//






const float INF = std::numeric_limits<float>::max();
maxFlow::maxFlow(seqGraph::Graph* graph1) {
    this->graph = graph1;
    this->graph->initMatrix();
    ConjugateMatrix = this->graph->getConjugateMatrix2();

    /// print matrix
    int N = 4 * graph1->getVCount()-4;
    for (int i=0; i< N; i++) {
        if (i == N-4) {
            matrixIdx2nodeName[i] = "s5+";
            std::cout<<"s\t";
        } else if (i == N-3) {
            matrixIdx2nodeName[i] = "s3+";
            std::cout<<"t\t";
        } else if (i == N-2) {
            matrixIdx2nodeName[i] = "t5+";
            std::cout<<"s'\t";
        } else if (i == N-1) {
            matrixIdx2nodeName[i] = "t3+";
            std::cout<<"t'\t";
        } else {
            int node_idx = i/4;
            int node_type = i%4;
            std::string node_name = (*this->graph->getVertices())[node_idx]->getOriginId();
            std::string node_type_str;
            if (node_type == 0) {
                node_type_str = " 5+";
            } else if (node_type == 1) {
                node_type_str = " 3+";
            } else if (node_type == 2) {
                node_type_str = " 5-";
            } else if (node_type == 3) {
                node_type_str = " 3-";
            }
            matrixIdx2nodeName[i] = node_name+node_type_str;
            std::cout<<node_name<<node_type_str<<"\t";
        }
    }
    V=N;
    sId=N-3;
    tId=N-2;
    std::cout<<" \t";
    for (int i = 1; i < N + 1; i++) {
        std::cout<<i<<"\t";
    }
    std::cout<<std::endl;
    for (int i=0; i< N; i++) {
        std::cout<<i+1<<"\t";
        for (int j=0; j< N; j++) {
            std::cout<<graph1->getConjugateMatrix2()[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }



//    this->graph->removeByGeneAndScore();
//    N = 2 * graph1->getVCount();
//    this->matched = new int[N + 1];
//    std::fill_n(this->matched,N+1, -1);
//    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
//    currentMatrix = this->graph->getConjugateMatrix();
//    this->originalGraph = graph1;
//    this->originalMatrix = new float*[N+1];
//    for(int i = 0; i < N+1; ++i)
//        originalMatrix[i] = new float [N+1];
//    std::memcpy(originalMatrix, currentMatrix, sizeof(int)*(N+1)*(N+1));
//    this->originalVertices = new std::vector<seqGraph::Vertex*>();
//    this->originalGraph = graph1;
//    for (auto item : *this->graph->getVertices()) {
//        this->originalVertices->push_back(new seqGraph::Vertex(*item));
//    }
//    this->originalJunctions = new std::vector<seqGraph::Junction*>();
//    for (auto item : *this->graph->getJunctions()) {
//        this->originalJunctions->push_back(new seqGraph::Junction(*item));
//    }
//    this->originalVertices = this->graph->getVertices();
//    if (VERBOSE == 2) {
//        this->graph->getConjugateMatrix().debugPrint();
//    }
}

maxFlow::~maxFlow() {
//    delete this->matched;
    if (graph != nullptr)
        free(graph);


}



bool  maxFlow::bfs(float **rGraph, int s, int t, int parent[])
{
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    std::queue<int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph[u][v] > 0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    if (visited[t] == true) {
        std::vector<int> path;
        for (int v = t; v != s; v = parent[v]) {
            path.push_back(v);
        }
        path.push_back(s);
        reverse(path.begin(), path.end());
        augmentationPaths.push_back(path);
    }

    return (visited[t] == true);
}

void maxFlow::fordFulkerson()
{
    int u, v;

//    float rGraph[V][V];
    int s = sId;
    int t = tId;
    float** rGraph = new float*[V];
    for(int i = 0; i < V; ++i) {
        rGraph[i] = new float[V];
        for(int j = 0; j < V; ++j) {
            rGraph[i][j] = ConjugateMatrix[i][j];
        }
    }

    int parent[V];

    float max_flow = 0.0f;

    float max_single_path_flow = 0.0f;
    std::vector<int> max_flow_path;

    while (bfs(rGraph, s, t, parent)) {
        float path_flow = FLT_MAX;
        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            path_flow = std::min(path_flow, rGraph[u][v]);
        }

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }

        max_flow += path_flow;

        if (path_flow > max_single_path_flow) {
            max_single_path_flow = path_flow;

            max_flow_path.clear();
            for (v = t; v != s; v = parent[v]) {
                max_flow_path.push_back(v);
            }
            max_flow_path.push_back(s);
            reverse(max_flow_path.begin(), max_flow_path.end());
        }
    }

    std::cout << "\nThe maximum possible flow is " << max_flow;

    for(int i = 0; i < V; ++i)
        delete [] rGraph[i];
    delete [] rGraph;

    std::cout << "\nAugmentation paths:\n";
    for (const auto& path : augmentationPaths) {
        for (int i = 0; i < path.size(); ++i) {
            std::cout << matrixIdx2nodeName[path[i]];
            if (i != path.size() - 1) {
                std::cout << " --> ";
            }
        }
        std::cout << '\n';
    }
}
