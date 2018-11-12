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


    // Display header

    printf("  \n");
    printf("%s",screenBar.c_str());
    printf("\tLotka-Volterra Probes interacting on a Stellar Network \n");
    printf("\t\tVersion: %s\n", VERSION);
    printf("\t\tCompiled: %s\n", __DATE__);
    printf("\t\tgit commit: %s \n", GIT_HASH);
    printf("%s",screenBar.c_str());
    printf("  \n");

            
    // Record start time
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
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
	                          	input.preyCarry1, input.preyCarry2,
	                          	input.predCarry1, input.predCarry2,
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
		                          	input.preyCarry1, input.preyCarry2,
		                          	input.predCarry1, input.predCarry2,
		                          	input.mutationRate1, input.mutationRate2,
		                          	input.outflowRate1, input.outflowRate2,
		                          	input.velocity1, input.velocity2);
      }
    else
      {
	graph.generateConstantLKParameters(input.initialPrey, input.initialPred,
					   input.preyGrow1, input.preyDeath1,
					   input.predGrow1, input.predDeath1,
					   input.preyCarry1, input.predCarry1,
					   input.mutationRate1, input.outflowRate1,
					   input.velocity1, t0);

      }


    // Write stellar network graph to file
    string graphFile = "stars.MSF.graph";
    graph.writeToFile(graphFile);

    // Begin Lotka-Volterra calculation

    t = 0;
    tnext = 0.0;

    printf("Initialising all LK Systems \n");
    graph.initialiseLKSystems(t,input.dt);


    while(t<input.tmax)
	{



	graph.updateLKSystems(t);

	if(tnext>input.tsnap)
	  {
	    printf("Time: %f \n",t);
	    graph.writeLKSnapshots(t);
	    tnext = 0.0;
	  }

	t = t+input.dt;
	tnext = tnext +input.dt;

	}



    // Simulation has now ended
    // Write elapsed runtime to the screen
    
    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        
    printf("%s",screenBar.c_str());
    printf("Run complete \n");
    printf("Wall Clock Runtime: %f s \n", time_span.count());
    printf("%s",screenBar.c_str());
        

    
    }






