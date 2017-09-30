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
    double outflow = 0.01;
    double velocity = 0.5;
    double t0 = 0.0;

    double initialPrey = 1.8;
    double initialPred = 1.8;

    double dt =0.0001;

    double tmax = 10.0;
    double t = 0.0;
    int ID = 1;

    Vector3D position1(0.0,0.0,0.0);

    LKVertex test(ID,position1,initialPrey, initialPred, preyGrow,
	    preyDeath,predGrow, predDeath, mutate,
	    outflow, velocity, t0);

    // do test integration on single system

    printf("Testing single system \n");
    test.initialiseLKSystem(t0, dt, initialPrey, initialPred);

    while(t < tmax)
	{

	test.updateLKSystem(t);
	cout << t << endl;
	test.writeToFile(t);
	t = t+dt;

	}

    printf("Single system test complete \n");


    printf("Attempting two-system test \n");
    vector<Vertex*> vertices;

    Vector3D position2(1.0,0.0,0.0);

    vertices.push_back(new LKVertex(ID,position1,initialPrey, initialPred, preyGrow,
    	    preyDeath,predGrow, predDeath, mutate,
    	    outflow, velocity, t0));

    ID++;
    initialPrey = 0.0;
    initialPred = 0.0;

    vertices.push_back(new LKVertex(ID,position2,initialPrey, initialPred, preyGrow,
    	    preyDeath,predGrow, predDeath, mutate,
    	    outflow, velocity, t0));


    // Create new graph

    vector<Edge*> edges;

    double range = 40.0;
    Graph graph(vertices,edges);
    graph.createNeighbourNetwork(range);

    string graphFile = "testgraph.txt";
    graph.writeToFile(graphFile);

    // Attempt coupled calculation

    t = 0;


    printf("Initialising all LK Systems \n");
    graph.initialiseLKSystems(t,dt);


    while(t<tmax)
	{

	printf("Time: %f \n",t);
	graph.updateLKSystems(t);

	t = t+dt;
	}


    }
