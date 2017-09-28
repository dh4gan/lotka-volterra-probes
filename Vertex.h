/*
 * Vertex.h
 *
 *  Created on: 20 Aug 2014
 *      Author: dhf
 */

#ifndef VERTEX_H_
#define VERTEX_H_

#include "Vector3D.h"

using namespace std;

// Note: have to forward declare Edge
class Edge;

#include "Edge.h"

class Vertex
    {
public:
    Vertex();
    Vertex(Vector3D pos);
    virtual ~Vertex();

    int getID(){return ident;}
    void setID(int i){ident = i;}

    int getEdgeNumber(){return nEdge;}
    void setEdgeNumber(int e){nEdge = e;}

    void setWeight(double w){weight = w;}
    double getWeight(){return weight;}

    void setPosition(Vector3D pos){position = pos;}
    Vector3D getPosition(){return position;}

    void setConnected(){connected = true;}
    bool isConnected(){return connected;}

    void setVisited(bool choice){visited = choice;}
    bool isVisited(){return visited;}

    void setComponentID(int ID){componentID = ID;}
    int getComponentID(){return componentID;}

    double calcVertexSeparation(Vertex* other){return getPosition().subtractVector(other->getPosition()).magVector();}

    bool inVector(vector<Vertex*> vertices);

    void setEdges(vector<Edge*> catalogue){edges=catalogue;}
    vector<Edge*> getEdges(){return edges;}

    vector<Vertex*> getConnectedVertices() {return connectedVertices;}

    int getNEdges(){return nEdge;}

    void addConnectedEdge(Edge* e, Vertex* other);


protected:
    int ident;
    int componentID;
    int nEdge;
    Vector3D position;

    double weight;

    vector<Vertex*> connectedVertices;
    vector<Edge*> edges;

    bool connected;
    bool visited;


    };


#endif /* VERTEX_H_ */
