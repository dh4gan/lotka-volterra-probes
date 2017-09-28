/*
 * Vertex.cpp
 *
 *  Created on: 20 Aug 2014
 *      Author: dhf
 */

#include "Vertex.h"

Vertex::Vertex()
    {
    Vector3D zerovector;

    ident=0;
    componentID = 0;
    position = zerovector;
    nEdge = 0;
    weight = 0.0;
    connected = false;
    visited = false;

    }

Vertex::Vertex(Vector3D pos)
    {
    ident = 0;
    componentID = 0;
    position = pos;
    nEdge = 0;
    weight = 0.0;
    connected = false;
    visited = false;
    }

Vertex::~Vertex()
    {
    }

void Vertex::addConnectedEdge(Edge* e, Vertex* other)
    {

    edges.push_back(e);
    connectedVertices.push_back(other);
    nEdge = nEdge+1;
    }

bool Vertex::inVector(vector<Vertex*> vertices)
    {

    Vertex* v = this;
   bool inVector = true;
    vector<Vertex*>::iterator itVertex = find(vertices.begin(),
    	    vertices.end(), v);

    if(itVertex==vertices.end())
	{
	inVector =false;
	}

    return inVector;
    }

