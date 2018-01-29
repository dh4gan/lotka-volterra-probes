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
#include "Constants.h"
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

    double dt =0.001;

    double tmax = 100.0;
    double t = 0.0;

    int nVertices = 100;

    int iseed = -37;
    double rmax = 10.0;
    double range = 10.0;

    Graph fullgraph;


    // Generate a cluster run
    fullgraph.generateCluster(iseed, nVertices,rmax);
    fullgraph.createNeighbourNetwork(range);

    string fullGraphFile = "fulltest.graph";
    fullgraph.writeToFile(fullGraphFile);

    //Use only minimum spanning Forest for computation
    Graph graph = fullgraph.minimumSpanningForest();

    // Generate parameters for LK systems

    graph.generateConstantLKParameters(initialPrey, initialPred, preyGrow,
		preyDeath, predGrow, predDeath, mutate, outflow, velocity, t0);



    // Write graph to file
    string graphFile = "test.graph";
    graph.writeToFile(graphFile);

    // Begin Calculation

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






