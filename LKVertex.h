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

    LKVertex(int initialPrey, int initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);

    LKVertex(Vector3D pos, int initialPrey, int initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);

    // TODO - set virtual methods for LKVertex in Vertex

    void initialiseSystem(double time, double dt, int initialPrey, int initialPredator);
    void computeFluxes(double t);
    void updateLKSystem(double t);

protected:

    int nPrey, nPredator;
    double preyGrowth, preyDeath, mutationRate;
    double predatorGrowth, predatorDeath;
    double velocity, outflowRate;
    double tzero, timestep;


    };



#endif /* LKVERTEX_H_ */
