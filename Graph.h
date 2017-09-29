/*
 * Graph.h
 *
 *  Created on: 20 Aug 2014
 *      Author: dhf
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include "Vertex.h"
#include <string>
#include <stdio.h>


class Graph
    {
public:
    Graph();
    Graph(vector<Vertex*> vert, vector<Edge*>edge);
    virtual ~Graph();

    int getNVertices(){return nVertices;}

    vector<Vertex*> getVertices(){return vertices;}
    vector<Edge*> getEdges() {return edges;}

    Vertex* getVertex(int iVertex){return vertices[iVertex];}
    Edge* getEdge(int iEdge){return edges[iEdge];}

    vector<int> getNEdges(){return nEdges;}
    int getTotalEdges(){return nTotalEdges;}
    int getNConnectedComponents(){return nConnectedComponents;}

    void addVertex(Vertex* vert);
    void addVertex(Vertex* vert, int index);
    void addEdge(Edge* edge);

    void removeVertex(int iVertex);
    void removeEdge(int iEdge);

    double calculateTotalEdgeWeight();

    Vertex* selectRandomVertex();

    int removeIsolatedVertices();
    void createEdgeCatalogues();

    Vertex* findNearestVertex(Vector3D location);

    bool edgeInGraph(Edge* e);
    bool vertexInGraph(Vertex* v);

    void addGraph(Graph other);

    void resetVertexWeights();
    void calculateEdgeDistanceWeights();
    void clearComponentData();

    void findConnectedComponents(Vertex* start);
    void findConnectedComponents();
    void findConnectedComponents(int iVertex);

    vector<int> getComponentCounts(){return nInComponent;}
    Graph minimumSpanningTree(Vertex* start);

    Graph minimumSpanningForest();
    Graph minimumSpanningTree();

    Graph shortestPath(Vertex* start, Vertex* destination, double &distance);

    void readFromFile(string &inputFileString);
    void writeToFile(string &outputFileString);


    void clearAllEdges();
    void createNeighbourNetwork(double range);
    void initialiseLKSystems(double t,double dt);
    void updateLKSystems(double t);


protected:

    vector<Vertex*> vertices;
    vector<int> vertexIndex;
    vector<Edge*> edges;

    vector<Vertex*> componentCentres;
    vector<int> nInComponent;

    int nVertices;
    int nTotalEdges;
    int nConnectedComponents;

    vector<int> nEdges; // No. of edges for each vertex


    };

#endif /* GRAPH_H_ */
