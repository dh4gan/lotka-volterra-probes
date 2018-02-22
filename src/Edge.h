/*
 * Edge.h
 *
 *  Created on: 20 Aug 2014
 *      Author: dhf
 */

#ifndef EDGE_H_
#define EDGE_H_

class Vertex;
#include "Vertex.h"

class Edge
    {
public:
    // Constructors and Destructor
    Edge();
    Edge(Vertex* b, Vertex* e);
    Edge(Vertex* b, Vertex* e, double w);

    virtual ~Edge();


    // Set and get Methods
    void setBeginning(Vertex* b){begin = b;}
    void setEnd(Vertex* e){end=e;}
    void setWeight(double w){weight = w;}

    Vertex* getBeginning(){return begin;}
    Vertex* getEnd(){return end;}
    double getWeight(){return weight;}


    // Other Methods

    void calcDistanceWeight();
    bool inVector(vector<Edge*> edges);

protected:
    double weight;
    Vertex* begin;
    Vertex* end;

    };

#endif /* EDGE_H_ */
