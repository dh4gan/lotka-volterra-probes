/*
 * LKVertex.cpp
 *
 *  Created on: 28 Sep 2017
 *      Author: dhf
 */

#include "LKVertex.h"
#include <string>
#include "Constants.h"
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

    preyCapacity = 1.0e30;
    predCapacity = 1.0e30;

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
    infoFile = NULL;


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

    preyCapacity = 1.0e30;
    predCapacity = 1.0e30;

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
    infoFile = NULL;

    }

LKVertex::LKVertex(int ID, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath,double preyCarry, double predCarry,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex(ID)
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;

    preyCapacity = preyCarry;
    predCapacity = predCarry;

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
    infoFile = NULL;

    }

LKVertex::LKVertex(int ID, Vector3D pos, double initialPrey, double initialPredator, double preyGrow,
	    double preyDie, double predGrow, double predDeath, double preyCarry, double predCarry,
	    double mutate, double outflow, double vel, double t0) :
		    Vertex(ID, pos)
    {
    nPrey = initialPrey;
    nPredator = initialPredator;

    preyGrowth = preyGrow;
    preyDeath = preyDie;

    predatorGrowth = predGrow;
    predatorDeath = predDeath;

    preyCapacity = preyCarry;
    predCapacity = predCarry;

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
    infoFile = NULL;
    }

void LKVertex::initialiseLKSystem(double time, double dt)

    {

    /*
     * Written 28/9/17 by dh4gan
     * sets up the LK system for integration
     */

    if(nPrey+nPredator > 0.0)
	{
	tzero = time;
	}
    else
	{
	tzero = -1.0e30;
	}
    timestep = dt;

    // Open file for writing
    ostringstream ss;

    ss << setw(5) << setfill('0') << ident;
    string fileNumber = ss.str();
    string logFileName = "LKSystem_"+fileNumber+".log";
    string infoFileName = "LKSystem_"+fileNumber+".info";

    printf("Initialising output %s \n", logFileName.c_str());

    outputFile = fopen(logFileName.c_str(), "w");
    infoFile = fopen(infoFileName.c_str(),"w");

    }

void LKVertex::initialiseLKSystem(double time, double dt, double initialPrey, double initialPredator)

    {
    initialiseLKSystem(time, dt);
    nPrey = initialPrey;
    nPredator = initialPredator;
    }

void LKVertex::determineTZero(double time)
    {

    //printf("ID %i, tzero %f, N+P %f \n", ident,tzero,nPrey+nPredator);
    if(tzero<0.0 and nPrey+nPredator>0.0)
	{
	tzero = time;
	}

    }


void LKVertex::setLKParameters(double initialPrey, double initialPred,double preyGrow,
	    double preyDie, double predGrow, double predDeath,double preyCarry,double predCarry,
	    double mutate, double outflow, double vel, double t0)

    {

    /*
     * Written 29/1/18 by dh4gan
     * Sets LK system to the following parameters
     *
     */

    nPrey = initialPrey;
    nPredator = initialPred;
    preyGrowth = preyGrow;

    predatorGrowth = predGrow;
    predatorDeath = predDeath;

    preyCapacity = preyCarry;
    predCapacity = predCarry;

    mutationRate = mutate;
    outflowRate = outflow;
    probeVelocity = vel;
    tzero = t0;

    }

void LKVertex::writeLKParametersToFile()

{
  /*
   * Written 2/8/17 by dh4gan
   * Writes the LK system parameters to an information file
   *
   */

    time_t rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    fprintf(infoFile,"Parameters for LKVertex %i \n", ident);
    fprintf(infoFile, "Time: %s \n", asctime(timeinfo));
    fprintf(infoFile,"Position: %f %f %f \n", position.elements[0],position.elements[1],position.elements[2]);
    fprintf(infoFile, "%s",screenBar.c_str());

    // Print position, and number of connections

    fprintf(infoFile,"Initial Prey Count: %.4E \n", nPrey);
    fprintf(infoFile,"Initial Predator Count: %.4E \n", nPredator);

    fprintf(infoFile,"Prey Growth Rate: %.4E \n", preyGrowth);
    fprintf(infoFile,"Prey Death Rate: %.4E \n", preyDeath);
    fprintf(infoFile,"Prey Carrying Capacity: %.4E \n", preyCapacity);

    fprintf(infoFile,"Predator Growth Rate: %.4E \n", predatorGrowth);
    fprintf(infoFile,"Predator Death Rate: %.4E \n", predatorDeath);
    fprintf(infoFile,"Predator Carrying Capacity: %.4E \n", predCapacity);

    fprintf(infoFile, "Mutation Rate: %.4E \n", mutationRate);
    fprintf(infoFile, "Outflow Rate: %.4E \n", outflowRate);
    fprintf(infoFile, "Probe Velocity: %.4E c\n", probeVelocity/c_kpc_Myr);



    fclose(infoFile);

}

void LKVertex::writePopulationsToFile(double time)
    {
    /*
     * written 28/9/17 by dh4gan
     * Writes the state of the LK system to file
     *
     */

    int nVariables = 7;
    string formatString = "";

	for (int i=0; i<nVariables; i++)
	    {

	    formatString = formatString + "%.4E  ";
	    }
	formatString = formatString + "\n";

	fprintf(outputFile,formatString.c_str(), time, nPrey, nPredator, ratePrey, ratePredator, preyOut, predatorOut);
    }

void LKVertex::updateLKSystem(double t)
    {
    /*
     * Written 28/9/17 by dh4gan
     * Integrates the LK system by one timestep
     * Assumes both predators and prey practise logistic growth
     */

    // Calculate the rate of change of prey = birth rate - death rate - prey leaving + prey arriving

    ratePrey = preyGrowth*nPrey*(1.0 - nPrey/preyCapacity) - preyDeath*nPrey*nPredator - preyOut + preyIn;

    // same for predators (but also include a mutation term from prey into predators)
    ratePredator = predatorGrowth*nPredator*nPrey*(1.0 - nPredator/predCapacity) - predatorDeath*nPredator + mutationRate*nPrey - predatorOut + predatorIn;

    nPrey = nPrey + ratePrey*timestep;
    nPredator = nPredator + ratePredator*timestep;

    if(nPrey < 0.0) nPrey = 0.0;
    if(nPredator < 0.0) nPredator = 0.0;
    determineTZero(t);

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

    if (tzero >= 0.0) // outward flux only happens if prey/predators present (t>tzero)
	{
	// loop over each connected vertex
	for (int i = 0; i < nConnected; i++)
	    {

	    // Flux only active if distance/speed < t-tzero

	    v = connected[i];
	    distance = calcVertexSeparation(v);
	    //printf("Distance %f, velocity %f, time-tzero %f \n", distance, probeVelocity, t-tzero);

	    if (t - tzero > distance / probeVelocity)
		{
		outwardPrey = outflowRate * nPrey * probeVelocity / distance;
		outwardPredator = outflowRate * nPredator * probeVelocity
			/ distance;
		preyOut = preyOut + outwardPrey;
		predatorOut = predatorOut + outwardPredator;

		v->addInwardPrey(outwardPrey);
		v->addInwardPredator(outwardPredator);
		}

	    }

	}
    }

