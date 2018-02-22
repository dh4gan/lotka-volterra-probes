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
    double tnext = 0.0;
    bool writeSnapshot;

    // Display header

    printf("  \n");
    printf("*********************************************** \n");
    printf("    Lotka-Volterra Probes interacting on a Stellar Network \n");
    printf("    Date Created : 1st January 2018 \n");
    printf("    Current Version: 22nd February 2018 \n");
    printf("*********************************************** \n");
    printf("  \n");

    // Either read parameter filename from argv[] or command line

    string fileString;

    if (argc == 2)
    	{
    	fileString = string(argv[1]);

    	}
        else
    	{
         cout << "What is the input file? " << endl;
         getline(cin, fileString);

    	}

    // Read in parameter file

    parFile input(fileString);
    input.readParams();


    // Write parameter details to screen and store in file
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

    string fullGraphFile = "stars.graph";
    fullgraph.writeToFile(fullGraphFile);

    // Reduce graph to minimum spanning Forest for computation
    Graph graph = fullgraph.minimumSpanningForest();

    // Generate parameters for all Lotka-Volterra systems

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
    string graphFile = "stars.MSF.graph";
    graph.writeToFile(graphFile);

    // Begin Lotka-Volterra calculation

    t = 0;

    printf("Initialising all LK Systems \n");
    graph.initialiseLKSystems(t,input.dt);

    while(t<input.tmax)
	{

	tnext = tnext +input.dt;

	if(tnext<input.tsnap)
	  {
	    writeSnapshot = true;
	  }
	printf("Time: %f \n",t);
	graph.updateLKSystems(t,writeSnapshot);

	t = t+input.dt;
	}


    }






