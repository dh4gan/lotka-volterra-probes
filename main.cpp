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


    double dt =0.001;

    double tmax = 50.0;
    double t = 0.0;
    int ID = 1;

    Vector3D position1(0.0,0.0,0.0);

//    LKVertex test(ID,position1,initialPrey, initialPred, preyGrow,
//	    preyDeath,predGrow, predDeath, mutate,
//	    outflow, velocity, t0);
//
//    // do test integration on single system
//
//    printf("Testing single system \n");
//    test.initialiseLKSystem(t0, dt, initialPrey, initialPred);
//
//    while(t < tmax)
//	{
//
//	test.updateLKSystem(t);
//	cout << t << endl;
//	test.writeToFile(t);
//	t = t+dt;
//
//	}
//
//    printf("Single system test complete \n");


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


    Vector3D position3(0.0,3.0,0.0);

    ID++;
    vertices.push_back(new LKVertex(ID,position3,initialPrey, initialPred, preyGrow,
      	    preyDeath,predGrow, predDeath, mutate,
      	    outflow, velocity, t0));

    Vector3D position4(3.0,0.0,0.0);

    ID++;
    vertices.push_back(new LKVertex(ID,position4,initialPrey, initialPred, preyGrow,
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


double randomReal()
    {
    // Returns a random real between zero and 1
    double number = (double) rand() / (double) RAND_MAX;
    return number;
    }


Graph generateGraphConstantParameters(){

    // TODO - write function that generates a graph from list of vertex Positions (assuming constant LV parameters)
    Graph graph;
    return graph;
}

Graph generateGraphRandomParameters(){

    // TODO - write function that generates a graph from list of vertex Positions (assuming randomised LV parameters)
    Graph graph;
    return graph;
}


// TODO - write function that generates a set of positions according to e.g. a cluster or GHZ formation
vector<Vector3D> generateGHZ(){

    vector<Vector3D> positions;

    return positions;
}

vector<Vector3D> generateCluster(){

    vector<Vector3D> positions;

    return positions;
}


vector<Body*> placeStarsGHZExponentialDensity(double &totalmass, int &nStars, double &innerRadius,
	double &outerRadius, double &scale, double &G, int &iseed)
    {

    /*
     * Written by dh4gan, 2/12/13
     * Place stars such that they form an annular GHZ with an exponential surface density profile
     * Periastra and apastra must fall within the GHZ annulus
     */

    int istar, success;
    double eccmax, eccmax2,incmax;
    double mass, radius, inc, semimaj, ecc, meanAnom, argPer, longAscend;
    double periastron, apastron, percent, timecounter;
    string bodyType;
    vector<Body*> stars;

    // All stars have the same mass

    mass = totalmass / float(nStars);
    radius = 1.0;
    srand(iseed);

    // Assert maximum eccentricity and inclination
    eccmax = 0.7;
    incmax = 0.5;


    // Scale surface density distribution for the range of radii involved

    double sigma_0 = (exp(-innerRadius/scale) - exp(-outerRadius/scale));

    // Randomly assign star positions

    istar = 0;

    percent = 0.0;
    timecounter = 0.0;

    for (istar = 0; istar < nStars; istar++)
	{

	percent = percent + 1 / float(nStars);
	if (percent > 0.1)
	    {
	    timecounter = timecounter + 10.0;
	    printf("%.0f %% complete \n", timecounter);
	    percent = 0.0;
	    }

	// Assign semimajor axis such that the surface density profile is correct

	semimaj = exp(-innerRadius/scale) -randomReal()*sigma_0;
	semimaj = -log(semimaj)*scale;


	// This defines maximum eccentricity the body can have to avoid going inside inner radius
	eccmax = outerRadius/semimaj - 1.0;
	eccmax2 = 1.0 - innerRadius/semimaj;

	if(eccmax2 < eccmax) eccmax = eccmax2;


	success = 0;

	while (success == 0)
	    {
	    ecc = randomReal() *eccmax;

	    // Check that orbit falls within the GHZ

	    periastron = semimaj * (1.0 - ecc);
	    apastron = semimaj * (1.0 + ecc);


	    success = 1;

	    if (periastron < innerRadius)
		{
		success = 0;
		}

	    if (apastron > outerRadius)
		{
		success = 0;
		}

	    }

	inc = -incmax + randomReal() * 2.0 * incmax;
	meanAnom = randomReal() * 2.0 * pi;
	argPer = randomReal() * 2.0 * pi;
	longAscend = randomReal() * 2.0 * pi;

	bodyType = "Star";

	stars.push_back(
		new Body(bodyType, bodyType, mass, radius, semimaj, ecc, inc,
			meanAnom, argPer, longAscend, G, totalmass));
	}

    for (istar = 0; istar < nStars; istar++)
	{
	stars[istar]->calcVectorFromOrbit(G, totalmass);
	}

    return stars;
    }

