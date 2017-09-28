/*
 * Edge.cpp
 *
 *  Created on: 20 Aug 2014
 *      Author: dhf
 */

#include "Edge.h"
#include "Vector3D.h"
#include <algorithm>

Edge::Edge()
    {
    begin=0;
    end =0;
    weight = 0.0;
    }

Edge::Edge(Vertex* b, Vertex* e)
    {

    /*
     * NOTE: default weight for edges is the distance between vertices
     */
    begin = b;
    end = e;
    calcDistanceWeight();
    }


Edge::Edge(Vertex* b, Vertex* e, double w)
    {
    begin = b;
    end = e;
    weight = w;
    }


Edge::~Edge()
    {
    }

void Edge::calcDistanceWeight()
    {
    weight = begin->calcVertexSeparation(end);
    }

bool Edge::inVector(vector<Edge*> edges)
    {
    Edge *e = this;
    bool inVector = true;
     vector<Edge*>::iterator itEdge = find(edges.begin(),
     	    edges.end(), e);

     if(itEdge==edges.end())
 	{
 	inVector =false;
 	}

     return inVector;
    }
