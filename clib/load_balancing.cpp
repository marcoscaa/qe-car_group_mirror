#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
using namespace std;

class directed_graph
{
public:
    vector< list<int> > adj_list;
public:
    directed_graph(int nbsp);
    directed_graph(int nRows, int nCols, int* matrix);
    void addEdge(int node, int from);
    bool rmEdge(int node, int from);
    list<int> bfsToEmpty(int initNode, double aveNEdge);
    void adjustWithPath(list<int> path);
};
directed_graph::directed_graph(int nbsp)
{
    list<int> emptyEdgeList;
    for (int i = 0; i < nbsp; i++)
        adj_list.push_back(emptyEdgeList);
}
directed_graph::directed_graph(int nRows, int nCols, int* matrix)
{
    list<int> emptyEdgeList;
    for (int i = 0; i < nRows; i++)
    {
        adj_list.push_back(emptyEdgeList);
        for (int j = 0; j < nCols; j++)
        {
            if (*(matrix + i*nCols + j) != 0)
            {
                addEdge(i, *(matrix + i*nCols + j) - 1);
            }
        }
    }
}
void directed_graph::addEdge(int node, int from)
{
    adj_list[node].push_back(from);
}
bool directed_graph::rmEdge(int node, int from)
{
    list<int>::iterator itr = find(adj_list[node].begin(), adj_list[node].end(), from);
    if (itr != adj_list[node].end())
    {
        // successfully remove
        adj_list[node].erase(itr);
        return true;
    }
    else
    {
        // fail to find
        return false;
    }
}
list<int> directed_graph::bfsToEmpty(int initNode, double aveNEdge)
{
    list<int> path;
    list<int> nodeQueue;
    map<int,int> checkMap;
    nodeQueue.push_back(initNode);
    checkMap.insert(pair<int,int>(initNode, initNode));

    if (adj_list[initNode].size() >= aveNEdge)
    {
        int emptyNode = -1;
        while (!nodeQueue.empty())
        {
            int currentNode = nodeQueue.front();
            nodeQueue.pop_front();

            for (list<int>::iterator itr = adj_list[currentNode].begin(); itr != adj_list[currentNode].end(); itr++)
            {
                if (checkMap.find(*itr) == checkMap.end())
                {
                    checkMap.insert(pair<int,int>(*itr, currentNode));
                    nodeQueue.push_back(*itr);

                    if (adj_list[*itr].size() < aveNEdge)
                    {
                        emptyNode = *itr;
                        break;
                    }
                }
            }

            if (emptyNode != -1)
                break;
        }

        path.push_front(emptyNode);
        map<int,int>::iterator itr = checkMap.find(emptyNode);
        while (itr->first != itr->second)
        {
            path.push_front(itr->second);
            itr = checkMap.find(itr->second);
        }
    }
    else
    {
        path.push_front(initNode);
    }

    return path;
}
void directed_graph::adjustWithPath(list<int> path)
{
    list<int>::iterator itr = path.begin();
    int prevNode = *itr;
    itr++;
    rmEdge(prevNode, *itr);
    while (*itr != path.back())
    {
        addEdge(*itr, prevNode);
        prevNode = *itr;
        itr++;
        rmEdge(prevNode, *itr);
    }
    addEdge(*itr, prevNode);
}

// remove non existing edges from prev_load, remove existing edges from curr_overlap
void rmPrevCalc(directed_graph& curr_load, directed_graph& curr_overlap)
{
    for (unsigned int i = 0; i < curr_load.adj_list.size(); i++)
    {
        list<int>::iterator itr = curr_load.adj_list[i].begin();
        while (itr != curr_load.adj_list[i].end())
        {
            bool successfulRm = false;
            successfulRm = curr_overlap.rmEdge(i, *itr);
            successfulRm = curr_overlap.rmEdge(*itr, i);

            if (!successfulRm)
                itr = curr_load.adj_list[i].erase(itr);
            else
                itr++;
        }
    }
}

void addNewCalc(directed_graph& curr_load, directed_graph& curr_overlap, double aveNEdge)
{
    for (unsigned int i = 0; i < curr_overlap.adj_list.size(); i++)
    {
        for (list<int>::iterator itr = curr_overlap.adj_list[i].begin(); itr != curr_overlap.adj_list[i].end(); itr++)
        {
            curr_load.addEdge(*itr, i);
            list<int> path = curr_load.bfsToEmpty(*itr, aveNEdge);

            curr_load.adjustWithPath(path);

            curr_overlap.rmEdge(*itr, i);
        }
    }
}

double aveNumEdges(directed_graph& graph)
{
    double edgeSum = 0;
    for (unsigned int i = 0; i < graph.adj_list.size(); i++)
        edgeSum += graph.adj_list[i].size();

    return edgeSum/graph.adj_list.size();
}

void writeMatrix(directed_graph curr_load, int nRows, int nCols, int* matrix)
{
    for (unsigned int i = 0; i < curr_load.adj_list.size(); i++)
    {
        list<int>::iterator itr = curr_load.adj_list[i].begin();
        for (int j = 0; j < nCols; j++)
        {
            if (itr != curr_load.adj_list[i].end())
            {
                *(matrix + i*nCols + j) = *itr + 1;
                itr++;
            }
            else
            {
                *(matrix + i*nCols + j) = 0;
            }
        }
    }
}

void finalAdjustment(directed_graph curr_load, double aveNEdge)
{
    for (unsigned int i = 0; i < curr_load.adj_list.size(); i++)
    {
        if (curr_load.adj_list[i].size() > aveNEdge + 1)
        {
            list<int> path = curr_load.bfsToEmpty(i, aveNEdge);
            curr_load.adjustWithPath(path);
        }
    }
}

extern "C" void load_balancing(int* neigh_ptr, int* nbsp_ptr, int* prev_load_ptr, int* overlap_ptr, int* load_ptr)
{
    int neigh = *neigh_ptr;
    int nbsp  = *nbsp_ptr;
    directed_graph curr_load(nbsp, neigh/2, prev_load_ptr);
    directed_graph curr_overlap(nbsp, neigh, overlap_ptr);
    double aveNumJobs = aveNumEdges(curr_overlap)/2;
    rmPrevCalc(curr_load, curr_overlap);
    addNewCalc(curr_load, curr_overlap, aveNumJobs);
    finalAdjustment(curr_load, aveNumJobs);
    writeMatrix(curr_load, nbsp, neigh/2, load_ptr);
}

/* HK: comment out otherwise get multiple definition of main in compilation
 *int main()
 *{
 *    int neigh = 6;
 *    int nbsp  = 6;
 *    int prev_load[6][3] = {{6, 3, 0}, {1, 3, 0}, {4, 5, 0}, {2, 5, 0}, {1, 6, 0}, {2, 4, 0}};
 *    int overlap[6][6]   = {{2, 3, 4, 5, 0, 0}, {1, 4, 3, 5, 6, 0}, {2, 1, 5, 0, 0, 0}, {2, 1, 5, 6, 0, 0}, {1, 3, 2, 6, 4, 0}, {2, 4, 5, 0, 0, 0}};
 *    int load[6][3];
 *
 *    load_balancing(&neigh, &nbsp, &prev_load[0][0], &overlap[0][0], &load[0][0]);
 *
 *    for (int i = 0; i < 6; i++)
 *    {
 *        for (int j = 0; j < 3; j++)
 *            cout << load[i][j] << " ";
 *        cout << endl;
 *    }
 *
 *    cout << "hello, world" << endl;
 *
 *    return 0;
 *}
 */
