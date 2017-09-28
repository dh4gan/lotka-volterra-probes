/*
 * LKVertex.h
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "Vertex.h"

#ifndef LKVERTEX_H_
#define LKVERTEX_H_


class LKVertex: public Vertex {

public:
    LKVertex();
    LKVertex(Vector3D pos);

    LKVertex(double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);

    LKVertex(Vector3D pos, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);

    // TODO - set virtual methods for LKVertex in Vertex

    void initialiseSystem(double time, double dt, double initialPrey, double initialPredator);
    void computeFluxes(double t);
    void updateLKSystem(double t);
    void writeToFile(double time);

protected:

    double nPrey, nPredator;
    double preyGrowth, preyDeath, mutationRate;
    double predatorGrowth, predatorDeath;
    double velocity, outflowRate;
    double tzero, timestep;

    double ratePrey, ratePredator, preyOut, preyIn, predatorOut, predatorIn;

    FILE* outputFile;


    };



#endif /* LKVERTEX_H_ */
