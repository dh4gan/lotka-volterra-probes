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
#include "parFile.h"

using namespace std;

int main()
    {

    double t = 0.0;
    double t0 = t;
    double dt = 1.0e-3;

    // Read in parameter file

    parFile input;

    input.readParams();

    // Generate graph

    Graph fullgraph;

    // Generate either a cluster run

    if(input.icChoice==1)
      {
	fullgraph.generateCluster(input.iseed, input.nVertices,input.rmax);
      }
    else
      {
	fullgraph.generateGHZ(input.iseed, input.nVertices, input.rmin, input.rmax, input.scaleLength);
      }

    fullgraph.createNeighbourNetwork(input.range);

    string fullGraphFile = "fulltest.graph";
    fullgraph.writeToFile(fullGraphFile);

    //Use only minimum spanning Forest for computation
    Graph graph = fullgraph.minimumSpanningForest();

    // Generate parameters for LK systems


    if(input.parChoice=="uniform")
      {

	graph.generateUniformLKParameters(input.initialPrey, input.initialPred,
	                          	input.preyGrow1, input.preyGrow2,
	                          	input.preyDeath1, input.preyDeath2,
	                          	input.predGrow1, input.predGrow2,
	                          	input.predDeath1, input.predDeath2,
	                          	input.mutationRate1, input.mutationRate2,
	                          	input.outflowRate1, input.outflowRate2,
	                          	input.velocity1, input.velocity2);

      }
    else if(input.parChoice=="gaussian")
      {

// TODO - write gaussian sampling of LK parameters
      }
    else
      {
	graph.generateConstantLKParameters(input.initialPrey, input.initialPred, input.preyGrow1,
			input.preyDeath1, input.predGrow1, input.predDeath1, input.mutationRate1, input.outflowRate1, input.velocity, t0);

      }




    // Write graph to file
    string graphFile = "test.graph";
    graph.writeToFile(graphFile);

    // Begin Calculation

    t = 0;

    printf("Initialising all LK Systems \n");
    graph.initialiseLKSystems(t,dt);

    while(t<input.tmax)
	{

	printf("Time: %f \n",t);
	graph.updateLKSystems(t);

	t = t+dt;
	}


    }






