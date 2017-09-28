/*
 * LKVertex.cpp
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "LKVertex.h"
#include <string>

LKVertex::LKVertex() :
    Vertex()
    {
    nPrey = 1;
    nPredator = 0;

    preyGrowth = 10.0;
    preyDeath = 1.0;

    predatorGrowth = 10.0;
    predatorDeath = 1.0;

    mutationRate = 0.0;
    outflowRate = 0.0;
    velocity = 1.0;
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

LKVertex::LKVertex(Vector3D pos) :
    Vertex(pos)
    {
    nPrey = 1;
    nPredator = 0;

    preyGrowth = 10.0;
    preyDeath = 1.0;

    predatorGrowth = 10.0;
    predatorDeath = 1.0;

    mutationRate = 0.0;
    outflowRate = 0.0;
    velocity = 1.0;
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

LKVertex::LKVertex(double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex()
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;
    predatorGrowth = predGrow;
    predatorDeath = predDeath;
    mutationRate = mutate;
    outflowRate = outflow;
    velocity = vel;
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

LKVertex::LKVertex(Vector3D pos, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex(pos)
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;
    predatorGrowth = predGrow;
    predatorDeath = predDeath;
    mutationRate = mutate;
    outflowRate = outflow;
    velocity = vel;
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
    string logFileName = "LKSystem.log";
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

