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


double randomReal();
double uniformSample(double min,double max);

Graph generateGraphConstantParameters(vector<Vector3D> positions, double initialPrey, double initialPred, double preyGrow,
double preyDeath,double predGrow, double predDeath, double mutate, double outflow, double velocity, double t0);

Graph generateGraphUniformParameters(vector<Vector3D> positions, double initialPrey, double initialPred,
	double preyGrowMin, double preyGrowMax,
	double preyDeathMin, double preyDeathMax,
	double predGrowMin, double predGrowMax,
	double predDeathMin, double predDeathMax,
	double mutateMin, double mutateMax,
	double outflowMin, double outflowMax,
	double velocityMin, double velocityMax);

vector<Vector3D> generateCluster(int nVertices);

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

    int nVertices = 10;


    vector<Vector3D> positions = generateCluster(nVertices);



    // Create new graph

    Graph graph = generateGraphConstantParameters(positions, initialPrey, initialPred, preyGrow,
		preyDeath, predGrow, predDeath, mutate, outflow, velocity, t0);


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


double randomReal()
    {
    // Returns a random real between zero and 1
    double number = (double) rand() / (double) RAND_MAX;
    return number;
    }

double uniformSample(double min, double max)
    {
    return min + (max-min)*randomReal();
    }


Graph generateGraphUniformParameters(vector<Vector3D> positions, double initialPrey, double initialPred,
	double preyGrowMin, double preyGrowMax,
	double preyDeathMin, double preyDeathMax,
	double predGrowMin, double predGrowMax,
	double predDeathMin, double predDeathMax,
	double mutateMin, double mutateMax,
	double outflowMin, double outflowMax,
	double velocityMin, double velocityMax){

    int ID;
    int nVertices = int(positions.size());
    double t0=0.0;

    printf("Generating graph with randomised parameters: %i vertices \n",nVertices);

    vector<Vertex*> vertices;


    for (int iVertex=0; iVertex<nVertices; iVertex++)
	{

	ID = iVertex+1;

	double preyGrow = uniformSample(preyGrowMin, preyGrowMax);
	double preyDeath = uniformSample(preyDeathMin, preyDeathMax);

	double predGrow = uniformSample(predGrowMin, predGrowMax);
	double predDeath = uniformSample(predDeathMin, predDeathMax);

	double mutate = uniformSample(mutateMin, mutateMax);
	double outflow = uniformSample(outflowMin, outflowMax);
	double velocity = uniformSample(velocityMin, velocityMax);


        vertices.push_back(new LKVertex(ID,positions[iVertex],initialPrey, initialPred, preyGrow,
        	    preyDeath,predGrow, predDeath, mutate,
        	    outflow, velocity, t0));

        // Only first vertex is given non-zero starting values
        if(iVertex==0)
            {
            initialPrey = 0.0;
            initialPred = 0.0;
            }
	}


        vector<Edge*> edges;

        double range = 40.0;
        Graph graph(vertices,edges);
        graph.createNeighbourNetwork(range);

    return graph;

}


vector<Vector3D> generateGHZ(int nVertices, double &innerRadius, double &outerRadius, double &scale){

    vector<Vector3D> positions;

    int istar, success;
       double eccmax, eccmax2,incmax;
       double mass, radius, inc, semimaj, ecc, meanAnom, argPer, longAscend;
       double periastron, apastron, percent, timecounter;
       string bodyType;
       vector<Body*> stars;

       // All stars have the same mass
       double totalmass = float(nVertices);
       double G = 4.49e-12; // G in kpc^3 Msol^{-1} kpc^{-2}
       radius = 1.0;
       srand(iseed);

       // Assert maximum eccentricity and inclination
       eccmax = 0.7;
       incmax = 0.5;


       // Scale surface density distribution for the range of radii involved

       double sigma_0 = (exp(-innerRadius/scale) - exp(-outerRadius/scale));

       // Randomly assign vertex positions

       percent = 0.0;
       timecounter = 0.0;

       for (int i = 0; i < nVertices; i++)
   	{

   	percent = percent + 1 / float(nVertices);
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

   	// ensure that eccentricity guarantees orbit entirely within GHZ
   	while (success == 0)
   	    {
   	    ecc = randomReal() *eccmax;

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

       for (int i = 0; i < nVertices; i++)
   	{
   	stars[istar]->calcVectorFromOrbit(G, totalmass);
   	}

       return stars;


    return positions;
}

vector<Vector3D> generateCluster(int nVertices){


    vector<Vector3D> positions;

    double rmax = 10.0;
    double rmin = 0.0;


    double thetamin = pi/2;
    double thetamax = pi/2;

    double phimin = 0.0;
    double phimax = 2.0*pi;
    int success = 0;

    double percent = 0.0;
    double timecounter= 0.0;

    for (int i=0; i<nVertices; i++)
	{

	percent = percent + 1 / float(nVertices);
   	if (percent > 0.1)
   	    {
   	    timecounter = timecounter + 10.0;
   	    printf("%.0f %% complete \n", timecounter);
   	    percent = 0.0;
   	    }

	success = 0;
	double r = 0;
	while(success==0)
	    {
	    r = uniformSample(rmin,rmax);

	    // (Normalised) density profile of a Plummer sphere (1,0.17)
	    double rho = pow(1 + r*r/(rmax*rmax),-2.5);

	    double randtest = uniformSample(0.0,1);

	    if (rho<randtest)
		{
		success=1;
		}
	}


	double theta = uniformSample(thetamin,thetamax);
	double phi = uniformSample(phimin,phimax);

	double x = r*sin(theta)*cos(phi);
	double y = r*sin(theta)*sin(phi);
	double z = 0.0;

	cout << "x " << x << endl;
	cout << "y " << y << endl;

	Vector3D nextVector(x,y,z);
	nextVector.printVector();
	positions.push_back(nextVector);

	}

    return positions;
}


