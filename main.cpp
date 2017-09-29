/*
 *
 *
 *  Created on: Sep 28, 2017
 *      Author: dh4gan
 *
 *  tests single LKSystem object
 *
 */

#include "LKVertex.h"
#include "Graph.h"
#include <iostream>
using namespace std;

int main()
    {

    // Set up single LK System

    double preyGrow = 0.6666;
    double preyDeath = 1.333;

    double predGrow = 1.0;
    double predDeath = 1.0;

    double mutate = 0.0;
    double outflow = 0.1;
    double velocity = 0.5;
    double t0 = 0.0;

    double initialPrey = 1.8;
    double initialPred = 1.8;

    double dt =0.0001;

    double tmax = 5.0;
    double t = 0.0;
    int ID = 1;

    Vector3D position1(0.0,0.0,0.0);

    LKVertex test(ID,position1,initialPrey, initialPred, preyGrow,
	    preyDeath,predGrow, predDeath, mutate,
	    outflow, velocity, t0);

    // do test integration on single system

    test.initialiseLKSystem(t0, dt, initialPrey, initialPred);

    while(t < tmax)
	{

	test.updateLKSystem(t);
	cout << t << endl;
	test.writeToFile(t);
	t = t+dt;
	}


    vector<Vertex*> vertices;

    Vector3D position2(10.0,0.0,0.0);

    vertices.push_back(new LKVertex(ID,initialPrey, initialPred, preyGrow,
    	    preyDeath,predGrow, predDeath, mutate,
    	    outflow, velocity, t0));

    ID++;
    initialPrey = 0.0;
    initialPred = 0.0;

    vertices.push_back(new LKVertex(ID,initialPrey, initialPred, preyGrow,
    	    preyDeath,predGrow, predDeath, mutate,
    	    outflow, velocity, t0));


    // Create new graph

    vector<Edge*> edges;

    double range = 40.0;
    Graph graph(vertices,edges);
    graph.createNeighbourNetwork(range);


    // Attempt coupled calculation

    t = 0;


    graph.initialiseLKSystems(t,dt);

    // TODO - how to set t0 when empty site begins receiving inflow
    while(t<tmax)
	{

	graph.updateLKSystems(t);

	t = t+dt;
	}


    }
