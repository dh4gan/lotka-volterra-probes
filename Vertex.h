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
    Vertex(int ID);
    Vertex(int ID, Vector3D pos);
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

    void clearEdges();
    void setEdges(vector<Edge*> catalogue){edges=catalogue;}
    vector<Edge*> getEdges(){return edges;}

    vector<Vertex*> getConnectedVertices() {return connectedVertices;}

    int getNEdges(){return nEdge;}

    void addConnectedEdge(Edge* e, Vertex* other);


    // virtual methods for LKVertex

    // set and get methods

        virtual void setNPrey(double n){}
        virtual double getNPrey(){return -1;}

        virtual void setNPredator(double n){}
        virtual double getNPredator(){return -1;}

        virtual void setPreyGrowth(double a){}
        virtual double getPreyGrowth(){return -1;}

        virtual void setPreyDeath(double a){}
        virtual double getPreyDeath(){return -1;}

        virtual void setPredatorGrowth(double a){}
        virtual double getPredatorGrowth(){return -1;}

        virtual void setPredatorDeath(double a){}
        virtual double getPredatorDeath(){return -1;}

        virtual void setMutationRate(double m){}
        virtual double getMutationRate(){return -1;}

        virtual void setOutflowRate(double o){}
        virtual double getOutflowRate(){return -1;}

        virtual void setProbeVelocity(double v){}
        virtual double getProbeVelocity(){return -1;}

        virtual void setTZero(double t){}
        virtual double getTZero(double t){return -1;}

        // Increment methods for inward/outward fluxes of predators/prey
        virtual void addOutwardPrey(double increment){}
        virtual void addInwardPrey(double increment){}

        virtual void addOutwardPredator(double increment){}
        virtual void addInwardPredator(double increment){}

        virtual void initialiseSystem(double time, double dt, double initialPrey, double initialPredator){};
        virtual void computeOutwardFlux(double t){};
        virtual void updateLKSystem(double t){};
        virtual void writeToFile(double time){};


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
