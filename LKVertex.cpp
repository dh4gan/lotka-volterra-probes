/*
 * LKVertex.cpp
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "LKVertex.h"

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

    }

LKVertex::LKVertex(int initialPrey, int initialPredator, double preyGrow,
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

    }

LKVertex::LKVertex(Vector3D pos, int initialPrey, int initialPredator, double preyGrow,
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

    }

void LKVertex::initialiseSystem(double time, double dt, int initialPrey, int initialPredator)

    {
    tzero = time;
    timestep = dt;
    nPrey = initialPrey;
    nPredator = initialPredator;
    }





