/*
 * LKVertex.h
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "Vertex.h"
#include <stdio.h>
#ifndef LKVERTEX_H_
#define LKVERTEX_H_


class LKVertex: public Vertex {

public:
    LKVertex(int ID);
    LKVertex(int ID, Vector3D pos);

    LKVertex(int ID, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);

    LKVertex(int ID, Vector3D pos, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);


    // set and get methods

    void setNPrey(double n){nPrey=n;}
    double getNPrey(){return nPrey;}

    void setNPredator(double n){nPredator=n;}
    double getNPredator(){return nPredator;}

    void setPreyGrowth(double a){preyGrowth=a;}
    double getPreyGrowth(){return preyGrowth;}

    void setPreyDeath(double a){preyDeath=a;}
    double getPreyDeath(){return preyDeath;}

    void setPredatorGrowth(double a){predatorGrowth=a;}
    double getPredatorGrowth(){return predatorGrowth;}

    void setPredatorDeath(double a){predatorDeath=a;}
    double getPredatorDeath(){return predatorDeath;}

    void setMutationRate(double m){mutationRate=m;}
    double getMutationRate(){return mutationRate;}

    void setOutflowRate(double o){outflowRate=o;}
    double getOutflowRate(){return outflowRate;}

    void setProbeVelocity(double v){probeVelocity=v;}
    double getProbeVelocity(){return probeVelocity;}

    void setTZero(double t){tzero=t;}
    double getTZero(double t){return tzero;}

    // Increment methods for inward/outward fluxes of predators/prey
    void addOutwardPrey(double increment){preyOut = preyOut+increment;}
    void addInwardPrey(double increment){preyIn = preyIn+increment;}

    void addOutwardPredator(double increment){predatorOut = predatorOut+increment;}
    void addInwardPredator(double increment){predatorIn = predatorIn+increment;}

    void initialiseLKSystem(double time, double dt);
    void initialiseLKSystem(double time, double dt, double initialPrey, double initialPredator);
    void setLKParameters(double initialPrey, double initialPred,double preyGrow,
	    double preyDie, double predGrow, double predDeath,
	    double mutate, double outflow, double vel, double t0);
    void determineTZero(double time);
    void computeOutwardFlux(double t);
    void updateLKSystem(double t);
    void writeToFile(double time);

protected:

    double nPrey, nPredator;
    double preyGrowth, preyDeath, mutationRate;
    double predatorGrowth, predatorDeath;
    double probeVelocity, outflowRate;
    double tzero, timestep;

    double ratePrey, ratePredator, preyOut, preyIn, predatorOut, predatorIn;

    FILE* outputFile;


    };



#endif /* LKVERTEX_H_ */
