/*
 * LKVertex.cpp
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "LKVertex.h"
#include <string>
#include <sstream>
#include <iomanip>

LKVertex::LKVertex(int ID) :
    Vertex(ID)
    {
    nPrey = 1;
    nPredator = 0;

    preyGrowth = 10.0;
    preyDeath = 1.0;

    predatorGrowth = 10.0;
    predatorDeath = 1.0;

    mutationRate = 0.0;
    outflowRate = 0.0;
    probeVelocity = 1.0;
    tzero = 0.0;
    timestep = 0.0;

    ratePrey = 0.0;
    ratePredator = 0.0;
    preyIn = 0.0;
    preyOut = 0.0;
    predatorOut = 0.0;
    predatorIn = 0.0;
    outputFile = NULL;


    }

LKVertex::LKVertex(int ID, Vector3D pos) :
    Vertex(ID,pos)
    {
    nPrey = 1;
    nPredator = 0;

    preyGrowth = 10.0;
    preyDeath = 1.0;

    predatorGrowth = 10.0;
    predatorDeath = 1.0;

    mutationRate = 0.0;
    outflowRate = 0.0;
    probeVelocity = 1.0;
    tzero = 0.0;
    timestep = 0.0;

    ratePrey = 0.0;
    ratePredator = 0.0;
    preyIn = 0.0;
    preyOut = 0.0;
    predatorOut = 0.0;
    predatorIn = 0.0;
    outputFile = NULL;

    }

LKVertex::LKVertex(int ID, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex(ID)
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;
    predatorGrowth = predGrow;
    predatorDeath = predDeath;
    mutationRate = mutate;
    outflowRate = outflow;
    probeVelocity = vel;
    tzero = t0;
    timestep = t0/100.0;

    ratePrey = 0.0;
    ratePredator = 0.0;
    preyIn = 0.0;
    preyOut = 0.0;
    predatorOut = 0.0;
    predatorIn = 0.0;
    outputFile = NULL;

    }

LKVertex::LKVertex(int ID, Vector3D pos, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex(ID, pos)
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;
    predatorGrowth = predGrow;
    predatorDeath = predDeath;
    mutationRate = mutate;
    outflowRate = outflow;
    probeVelocity = vel;
    tzero = t0;
    timestep = t0/100.0;

    ratePrey = 0.0;
    ratePredator = 0.0;
    preyIn = 0.0;
    preyOut = 0.0;
    predatorOut = 0.0;
    predatorIn = 0.0;
    outputFile = NULL;
    }

void LKVertex::initialiseSystem(double time, double dt, double initialPrey, double initialPredator)

    {

    /*
     * Written 28/9/17 by dh4gan
     * sets up the LK system for integration
     */

    tzero = time;
    timestep = dt;
    nPrey = initialPrey;
    nPredator = initialPredator;


    // Open file for writing
    ostringstream ss;

    ss << setw(5) << setfill('0') << ident;
    string fileNumber = ss.str();
    string logFileName = "LKSystem_"+fileNumber+".log";
    outputFile = fopen(logFileName.c_str(), "w");

    }

void LKVertex::writeToFile(double time)
    {
    /*
     * written 28/9/17 by dh4gan
     * Writes the state of the LK system to file
     *
     */

	string formatString = "%f  %.4E  %.4E  %.4E  %.4E \n";
	fprintf(outputFile,formatString.c_str(), time, nPrey, nPredator, ratePrey, ratePredator);
    }

void LKVertex::updateLKSystem(double t)
    {
    /*
     * Written 28/9/17 by dh4gan
     * Integrates the LK system by one timestep
     *
     */

    // Calculate the rate of change of prey = birth rate - death rate - prey leaving + prey arriving

    ratePrey = preyGrowth*nPrey - preyDeath*nPrey*nPredator - preyOut + preyIn;

    // same for predators (but also include a mutation term from prey into predators)
    ratePredator = predatorGrowth*nPredator*nPrey - predatorDeath*nPredator + mutationRate*nPrey - predatorOut + predatorIn;

    nPrey = nPrey + ratePrey*timestep;
    nPredator = nPredator + ratePredator*timestep;

    }

void LKVertex::computeOutwardFlux(double t)

    {

    /*
     * Written 29/9/17 by dh4gan
     * computes the flux of predators/prey to connected LKVertex objects
     * Also updates inward fluxes?
     */

    double distance, outwardPrey, outwardPredator;
    Vertex* v;
    // Get list of all connected vertices

    vector<Vertex*> connected = getConnectedVertices();

    int nConnected = connected.size();

    // loop over each connected vertex
    for (int i=0; i<nConnected; i++)
	{

	// Flux only active if distance/speed < t-tzero

	v = connected[i];
	distance = calcVertexSeparation(v);

	if(t-tzero < distance/probeVelocity)
	    {
	    outwardPrey = outflowRate*nPrey*probeVelocity/distance;
	    outwardPredator = outflowRate*nPredator*probeVelocity/distance;
	    preyOut = preyOut + outwardPrey;
	    predatorOut = predatorOut + outwardPredator;

	    v->addInwardPrey(outwardPrey);
	    v->addInwardPredator(outwardPredator);
	    }

	}

    }

