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
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
    {

    double t = 0.0;
    double t0 = t;

    // Display header

    printf("  \n");
    printf("*********************************************** \n");
    printf("    Lotka-Volterra Probes interacting on a Stellar Network \n");
    printf("    Date Created : 1st January 2018 \n");
    printf("    Current Version: 22nd February 2018 \n");
    printf("*********************************************** \n");
    printf("  \n");

    // Check if parameter file in args

    string fileString;

    if (argc == 2)
    	{
    	fileString = string(argv[1]);

    	}
        else
    	{
         printf("Reading in the Users body data \n ");
         cout << "What is the input file? " << endl;

         getline(cin, fileString);

    	}

    // Read in parameter file

    parFile input(fileString);
    input.readParams();

    input.writeParamsToScreen();
    input.writeParamsToFile();

    // Generate graph

    Graph fullgraph;

    // Generate either a cluster run

    if(input.icChoice=="cluster")
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
	graph.generateGaussianLKParameters(input.initialPrey, input.initialPred,
		                          	input.preyGrow1, input.preyGrow2,
		                          	input.preyDeath1, input.preyDeath2,
		                          	input.predGrow1, input.predGrow2,
		                          	input.predDeath1, input.predDeath2,
		                          	input.mutationRate1, input.mutationRate2,
		                          	input.outflowRate1, input.outflowRate2,
		                          	input.velocity1, input.velocity2);
      }
    else
      {
	graph.generateConstantLKParameters(input.initialPrey, input.initialPred, input.preyGrow1,
			input.preyDeath1, input.predGrow1, input.predDeath1, input.mutationRate1, input.outflowRate1, input.velocity1, t0);

      }


    // Write stellar network graph to file
    string graphFile = "test.graph";
    graph.writeToFile(graphFile);

    // Begin Lotka-Volterra calculation

    t = 0;

    printf("Initialising all LK Systems \n");
    graph.initialiseLKSystems(t,input.dt);

    while(t<input.tmax)
	{

	printf("Time: %f \n",t);
	graph.updateLKSystems(t);

	t = t+input.dt;
	}


    }






