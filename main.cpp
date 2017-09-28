/*
 * lotka_volterra_probes
 *
 *  Created on: Sep 28, 2017
 *      Author: dh4gan
 *
 *  carries out coupled Lotka Volterra equations on stellar networks
 *
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Body.h"
#include "Constants.h"
#include "Graph.h"
#include <omp.h>

#include <fstream>
#include <sstream>
using namespace std;


// TODO - write functions to read input parameters from file

int readParameters(ifstream &inputStream, string &prefixString, int &nRuns,
	int &nStars,double &dt, double &tMax,int &iseed, string &GHZChoice, double &innerRadius, double &outerRadius);


vector<Body*> placeStarsRandom(double &totalmass, int &nStars, double &semimajMax, double &G,
	int &iseed);
vector<Body*> placeStarsGHZUniformDensity(double &totalmass, int &nStars, double &innerRadius,
	double &outerRadius, double &G, int &iseed);

vector<Body*> placeStarsGHZExponentialDensity(double &totalmass, int &nStars, double &innerRadius,
	double &outerRadius, double &scale, double &G, int &iseed);

int main()
    {

    int irun, istar, jstar, itime;
    int nStars, nRuns, nTime;
    int iseed, ntransit;
    double time, tMax, dt;
    double percent, timecounter;
    double innerRadius = 6.0;
    double outerRadius = 10.0;
    double scale = 3.5; // Scale radius of surface density distribution
    double G = 4.49e-12; // G in kpc^3 Msol^{-1} kpc^{-2}
    double totalMass = 9.76e10; // Set such that stars at 8 kpc have period of 230 Myr (like the Sun)

    vector<Body*> stars;

    string GHZChoice="exponential";
    string inputFileString, runNumString, timeNumString, prefixString;
    string logFileString, systemsFileString, graphFileString, treeFileString, runFileString, timeFileString, outputString;

    FILE *outputLog, *outputSystem;

    /*
     * 1. Read in input parameters
     *
     */

    inputFileString = "lotka_volterra_probes.params";

    printf("Reading in main parameter file %s \n", inputFileString.c_str());


    ifstream myfile(inputFileString.c_str());
    readParameters(myfile, prefixString, nRuns, nStars, dt, tMax,iseed,GHZChoice, innerRadius,outerRadius);

    printf("Star Number: %i \n", nStars);
    printf("GHZ defined to be %s \n", GHZChoice.c_str());
    printf("Inner Radius: %f kpc\nOuter Radius: %f kpc\nScale Length: %f kpc\n", innerRadius, outerRadius, scale);


    /*
     * 2. set up stellar network graph
     */



    /*
     * 3. Set up Lotka-Volterra systems (global fixed timestep)?
     */


    /*
     * 4. Begin time evolution
     */


	/*
	 * 4a. Compute LK system evolution
	 */

	/*
	 * 4b. Move stars according to orbits
	 */


	/*
	 * 4c. Repeat
	 */

    /*
     * 2.  Define planet properties
     *
     */

    vector<double>planetrad(nStars,0.0);
    vector<double>planeta(nStars,0.0);
    vector<Vector3D>starspin(nStars);
    vector<Vector3D>planetavector(nStars);
    vector<double> planeti(nStars,0.0);

    // Matrix element (i,j) is non-zero if i can observe j's transit - in general NOT SYMMETRIC
    float transit[nStars][nStars];

    // This matrix records all observed transits over all timesteps
    int cumulativeTransit[nStars][nStars];

    nTime = int(tMax / dt);

    double nRunZeros = int(log10(nRuns) + 1);
    double nTimeZeros = int(log10(nTime) + 1);

    /*
     * 3. Begin MCR process
     */

    // Begin loop over runs

printf("Beginning MCR \n");
    for (irun = 0; irun < nRuns; irun++)
	{

	iseed = iseed - 2.0 * irun;

	// Place Stars
	printf("-----\n");
	printf("Placing Stars for Run %i, seed %i \n", irun + 1, iseed);

	if (GHZChoice == "random")
	    {
	    stars = placeStarsRandom(totalMass, nStars, outerRadius, G, iseed);
	    }
	else if (GHZChoice == "uniform")
	    {
	    stars = placeStarsGHZUniformDensity(totalMass, nStars, innerRadius,
		    outerRadius, G, iseed);
	    }
	else
	    {
	    stars = placeStarsGHZExponentialDensity(totalMass, nStars,
		    innerRadius, outerRadius, scale, G, iseed);
	    }


	// Set up planet properties
	setPlanetPropertiesRandom(nStars, planeta, planeti, planetrad,planetavector, starspin);
	//setPlanetPropertiesFromData(nStars, planeta, planeti, planetrad,planetavector, starspin);

	// Set up output file names

	// Log file

	ostringstream convert;
	convert << irun + 1;

	runNumString = convert.str();
	while (runNumString.length() < nRunZeros)
	    {
	    runNumString = "0" + runNumString;
	    }

	logFileString = prefixString + "_log."+ runNumString;

	printf("Opening log file for run %i, %s \n ", irun, logFileString.c_str());

	outputLog = fopen(logFileString.c_str(), "w");

	// Write header data
	fprintf(outputLog, "%+.4E  %i \n", dt, nStars);

	// Write system properties to file

	systemsFileString = prefixString+"_systems."+runNumString;

	//strcpy(outputFile,systemsFileString.c_str());

	printf("Write system properties to file %s \n", systemsFileString.c_str());

	outputSystem = fopen(systemsFileString.c_str(), "w");

	for (istar=0; istar < nStars; istar++)
	    {
	    outputString = "";

	    // Start with star's position


	    // Print star's spin vector and planet properties
		fprintf(outputSystem, "%+.4E   %+.4E   %+.4E  ",
			starspin[istar].elements[0],
			starspin[istar].elements[1],
			starspin[istar].elements[2]);

		fprintf(outputSystem, "%+.4E   %+.4E \n  ",
					planeta[istar], planeti[istar]);
	    }
	    fclose(outputSystem);


	// Loop over time

	time = 0.0;
	percent = 0.0;
	timecounter = 0.0;

	// Set Cumulative Matrix to Zero

#pragma omp parallel default(none) \
	shared(transit) \
		private(istar,jstar)
		{
#pragma omp for schedule(runtime) ordered
		for (istar = 0; istar < nStars; istar++)
		    {
		    for(jstar=0;jstar<nStars; jstar++)
			{
			cumulativeTransit[istar][jstar]=0;
			}
		    }
		}


	for (itime=0; itime< nTime; itime++)

	    {
	    time = time + dt;
	    printf("Time: %f \n", time);
	    percent = percent + dt / float(tMax);
	    if (percent > 0.1)
		{
		timecounter = timecounter + 10.0;
		printf("%.0f %% complete \n", timecounter);
		percent = 0.0;
		}


	    // Reset graph objects for this timestep
	    Graph fullGraph;
	    Graph minSpanTree;

	    // Update mean anomalies of all bodies, and recalculate positions

#pragma omp parallel default(none) \
	shared(stars,G,totalMass,dt,nStars) \
		private(istar)
		{
#pragma omp for schedule(runtime) ordered
		for (istar = 0; istar < nStars; istar++)
		    {
		    stars[istar]->moveAlongOrbit(G, totalMass, dt);

		    // Add to the graph object
		    Vertex* v = new Vertex(stars[istar]->getPosition());
		    fullGraph.addVertex(v);

		    // Clear transit matrix for next run

		    for(jstar=0;jstar<nStars; jstar++)
			{
			transit[istar][jstar]=0.0;
			}

		    }
		}


	    ntransit = 0;

	    // Begin loop over all star pairs

	    for (istar = 0; istar < nStars; istar++)
		{

		// Calculate intersections

#pragma omp parallel default(none) \
	shared(stars,beams,dt,nStars,transit,ibeam) \
		private(jstar) \
		    {
#pragma omp for schedule(runtime) ordered
	    for(jstar=0;jstar<nStars; jstar++)
		{

		// Skip if istar = jstar
		if(istar==jstar)
		    {
		    continue;
		    }

		// Test for transit connection

		bool connected = checkTransitBridge(stars[istar],stars[jstar],starspin[jstar],planeta[jstar],planetrad[jstar]);

		// If connected, then add to list

		if(connected)
		    {
		    transit[istar][jstar] = stars[istar]->calcSeparation(stars[jstar]);
		    ntransit +=1;

		    cumulativeTransit[istar][jstar] = 1;

		    // Add Edge to the graph
		    Edge* e = new Edge(fullGraph.getVertex(istar),fullGraph.getVertex(jstar), stars[istar]->calcSeparation(stars[jstar]));

		    // Add Edge to each vertex's catalogue of edges

		    fullGraph.getVertex(istar)->addConnectedEdge(e, fullGraph.getVertex(jstar));
		    fullGraph.getVertex(jstar)->addConnectedEdge(e, fullGraph.getVertex(istar));

		    // Finally, add edge to graph
		    fullGraph.addEdge(e);

		    }


		}

	    } // End of loop over star pairs

	// Strip away isolated vertices
	//printf("Removing Isolated Vertices \n");
	int nIsolated = fullGraph.removeIsolatedVertices();
	//printf("%i Vertices remaining: %i edges \n", fullGraph.getNVertices(), fullGraph.getTotalEdges());
	fullGraph.findConnectedComponents();

	// Calculate the MSF

	//printf("Calculating Minimum Spanning Forest \n");
	minSpanTree = fullGraph.minimumSpanningForest();
	minSpanTree.findConnectedComponents();

	// Open files to hold timestep data
	    ostringstream convert;
	    convert << itime + 1;

	    timeNumString = convert.str();
	    while (timeNumString.length() < nTimeZeros)
		{
		timeNumString = "0" + timeNumString;
		}

	    graphFileString = prefixString + "_run_" + runNumString + "."
		    + timeNumString;

	    treeFileString = prefixString + "_MSF_run_" + runNumString + "."
	    		    + timeNumString;

	    /*printf("Writing graph and MSF data for timestep %i, %s %s \n ", itime,
		    graphFileString.c_str(), treeFileString.c_str());*/


	// Write the full graph to file
	fullGraph.writeToFile(graphFileString);

	// Write MSF to file
	minSpanTree.writeToFile(treeFileString);

	//printf("Graph Files written: updating log file \n");
	// Write data about this timestep to the run log file

	fprintf(outputLog,"%+.4E  %i %i %i %f \n", time, ntransit, nIsolated, fullGraph.getNConnectedComponents(), fullGraph.calculateTotalEdgeWeight());
	fflush(outputLog);
	    }

	// End of loop over time

	// Construct a cumulative graph using final positions of vertices

	Graph cumulativeGraph;


	printf("Calculating cumulative network for run %i \n", irun);
	for(istar=0; istar<nStars; istar++)
	    {
	    Vertex* v = new Vertex(stars[istar]->getPosition());
	    cumulativeGraph.addVertex(v);
	    }

	int nCumulativeTransit=0;
	   for (istar = 0; istar < nStars; istar++)
			{

			// Calculate intersections

	#pragma omp parallel default(none) \
		shared(stars,beams,dt,nStars,transit,ibeam) \
			private(jstar) \
			    {
	#pragma omp for schedule(runtime) ordered
		    for(jstar=0;jstar<nStars; jstar++)
			{

			// Skip if istar = jstar
			if(istar==jstar)
			    {
			    continue;
			    }

			// Test for transit connection

			bool connected = cumulativeTransit[istar][jstar]>0.0 or cumulativeTransit[jstar][istar]>0.0;

			// If connected, then add to list

			if(connected)
			    {

			    nCumulativeTransit +=1;

			    cumulativeTransit[istar][jstar] = 1;

			    // Add Edge to the graph
			    Edge* e = new Edge(cumulativeGraph.getVertex(istar),cumulativeGraph.getVertex(jstar), stars[istar]->calcSeparation(stars[jstar]));

			    // Add Edge to each vertex's catalogue of edges

			    cumulativeGraph.getVertex(istar)->addConnectedEdge(e, cumulativeGraph.getVertex(jstar));
			    cumulativeGraph.getVertex(jstar)->addConnectedEdge(e, cumulativeGraph.getVertex(istar));

			    // Finally, add edge to graph
			    cumulativeGraph.addEdge(e);

			    }


			}

		    } // End of loop over star pairs

	   // Write cumulative graph (and its MSF) to file

		// Strip away isolated vertices
		//printf("Removing Isolated Vertices \n");
		int nIsolated = cumulativeGraph.removeIsolatedVertices();
		//printf("%i isolated vertices in cumulative graph \n",nIsolated);
		//printf("%i Vertices remaining: %i edges \n", cumulativeGraph.getNVertices(), cumulativeGraph.getTotalEdges());
		cumulativeGraph.findConnectedComponents();

		// Calculate the MSF

		printf("Calculating Minimum Spanning Forest \n");
		Graph minSpanTree = cumulativeGraph.minimumSpanningForest();
		minSpanTree.findConnectedComponents();

		// Open files to hold timestep data
		    //ostringstream convert;
		    convert << itime + 1;

		    timeNumString = convert.str();
		    while (timeNumString.length() < nTimeZeros)
			{
			timeNumString = "0" + timeNumString;
			}

		    graphFileString = prefixString + "_cumulative_" + runNumString;

		    treeFileString = prefixString + "_MSF_cumulative_" + runNumString;

		    /*printf("Writing cumulative graph and MSF data for timestep %i, %s %s \n ", itime,
			    graphFileString.c_str(), treeFileString.c_str());*/


		// Write the full graph to file
		cumulativeGraph.writeToFile(graphFileString);

		// Write MSF to file
		minSpanTree.writeToFile(treeFileString);

	// Reset arrays for next MCR

	printf("Run %i complete \n", irun + 1);
	stars.clear();

	fclose(outputLog);
	// End of loop over runs
	}

    }
    // End of program

double randomReal()
    {
    // Returns a random real between zero and 1
    double number = (double) rand() / (double) RAND_MAX;
    return number;
    }

void setPlanetPropertiesRandom(int &nStars, vector<double> &planeta, vector<double> &planeti, vector<double> &planetrad, vector<Vector3D> &planetavector, vector<Vector3D>&starspin)
    {
    /*
     * Written by dh4gan
     * This method sets up uniform distributions of planet semimajor axis
     * and inclination
     *
     */

    double startheta, starphi;

    for (int i=0; i<nStars; i++)

	{

	// Randomly sample planet semimajor axis and inclination

	planeta[i] = randomReal()*100.0*AUtokpc;
	planeti[i] = randomReal()*pi;

	// Randomly orient stellar spin vector
	startheta = randomReal()*pi;
	starphi = randomReal()*2.0*pi;

	starspin[i].defineFromSpherical(1.0, startheta,starphi);

	}

    }

void setPlanetPropertiesFromData(int &nStars, vector<double> &planeta,
	vector<double>&planeti, vector<double>&planetrad,
	vector<Vector3D>&planetavector, vector<Vector3D>&starspin)
    {

    // Read exoplanet data (histogram format)

    int i, j;

    vector<float> exo_a_hist, exo_i_hist;
    vector<float> exo_a_edges, exo_i_edges;

    string aFileString = "ahistogram.dat";
    string iFileString = "ahistogram.dat";
    FILE *aFile = fopen(aFileString.c_str(), "r");
    FILE *iFile = fopen(iFileString.c_str(), "r");

    printf("Reading Exoplanet Data from Files: \n %s\n %s\n",aFileString.c_str(), iFileString.c_str());
    readExoplanetData(aFile, iFile, exo_a_hist, exo_a_edges, exo_i_hist,
	    exo_i_edges);

    printf("File Read complete \n");

    double max_a = exo_a_edges[exo_a_edges.size()-1];
    double max_i = exo_i_edges[exo_i_edges.size()-1];

    printf("Maximum a: %+.4E\n", max_a);
    printf("Maximum i: %+.4E\n",max_i);

    double randtest, histvalue, startheta, starphi;
    // Use accept-reject to sample these distributions

    // Do Semimajor Axis first
    for (i = 0; i < nStars; i++)
	{

	randtest = 1.0e30;
	histvalue = 0.0;

	while (randtest > histvalue)
	    {
	    planeta[i] = randomReal() * max_a;

	    randtest = randomReal();

	    // Find planeta in exo_a_edges
	    j = 0;
	    while (planeta[i] > exo_a_edges[j])
		{
		j++;
		}

	    // Compare histogram value with randtest: if randtest greater than histogram value, repeat
	    histvalue = exo_a_hist[j];

	    printf("%+.4E %+.4E %i\n", planeta[i], exo_a_edges[j], j);

	    }

	planeta[i] = planeta[i]*AUtokpc;

	printf("%i: planeta: %+.4E\n",i, planeta[i]);

	}

    // Now Do Inclination
    for (i = 0; i < nStars; i++)
	{

//	randtest = 0.0;
//	histvalue = 1.0e30;
//
//	while (randtest > histvalue)
//	    {
//	    planeti[i] = randomReal() * max_i;
//	    randtest = randomReal();
//
//	    // Find planeta in exo_a_edges
//	    j = 0;
//	    while (planeti[i] < exo_i_edges[j])
//		{
//		j++;
//		}
//
//	    // Compare histogram value with randtest: if randtest greater than histogram value, repeat
//	    histvalue = exo_i_hist[j];

//	    }

    planeti[i] = randomReal()*pi;

	// Randomly orient stellar spin vector
	startheta = randomReal() * pi;
	starphi = randomReal() * 2.0 * pi;

	starspin[i].defineFromSpherical(1.0, startheta, starphi);

	}
    }

void readExoplanetData(FILE* aFile, FILE* iFile, vector<float> &exo_a_hist,
	vector<float> &exo_a_edges, vector<float>&exo_i_hist,
	vector<float> &exo_i_edges)
    {

    /*
     * Written 29/10/14 by dh4gan
     * This method reads in histogram data for exoplanet semimajor axis
     * and orbital inclination
     * This method requires the FILE pointers to already be initialised
     * (i.e. the files must already be open)
     *
     */

    int nSemiMajorAxis, nInclination;
    float temp1=0.0f;
    float temp2 = 0.0f;

    // Read number of entries in each histogram File

    fscanf(aFile, "%i\n", &nSemiMajorAxis);
    fscanf(iFile, "%i\n", &nInclination);


    printf("Semi-Major Axis File has %i entries \n",nSemiMajorAxis);
    printf("Inclination file has %i entries \n", nInclination);

    // Read entirety of semimajor axis file

    for (int i = 0; i < nSemiMajorAxis; i++)
	{
	fscanf(aFile, "%F %F", &temp1, &temp2);

	exo_a_edges.push_back(temp1);
	exo_a_hist.push_back(temp2);
	}

    fclose(aFile);

    printf("Semi Major axis data read \n");
    for (int i = 0; i < nInclination; i++)
	{
	fscanf(iFile, "%F %F", &temp1, &temp2);
	exo_i_edges.push_back(temp1);
	exo_i_hist.push_back(temp2);
	}
    printf("Inclination data read \n");

    fclose(iFile);
    }




bool checkTransitBridge(Body* ibody, Body* jbody, Vector3D starspin, double planeta, double planetrad)
    {

    // Calculate r_ij = rj - ri (line of sight vector)

    Vector3D rij = jbody->getPosition().subtractVector(ibody->getPosition());

    // This forms unit vector for plane to project jbody system onto (assumes ibody's planet a negligible)
    rij = rij.unitVector();

    // calculate angle between line of sight and stellar spin vector

    double angle = -rij.dotProduct(starspin.unitVector());
    angle = acos(angle);

    // Subtract pi/2 to get angle between planet orbit and line of sight (for zero inclination orbits)
    angle = fabs(angle - pi/2.0);

    // Calculate projected distance from this angle
    double projection = planeta*angle;

    double proj_crit = planetrad + jbody->getRadius();
    proj_crit = proj_crit*rsoltokpc;

    //printf ("Projections: %+.4E   %+.4E   %+.4E \n", angle, projection, proj_crit);
    if (projection < proj_crit)
	{
	return true;
	}
    else
	{
	return false;
	}
    }

vector<Body*> placeStarsRandom(double &totalmass, int &nStars,
	double &semimajMax, double &G, int &iseed)
    {

    int istar;
    double mass, radius, inc, semimaj, ecc, meanAnom, argPer, longAscend;
    string bodyType;
    vector<Body*> stars;

    // All stars have the same mass

    mass = totalmass / float(nStars);
    radius = 1.0;
    srand(iseed);

    // Randomly assign star positions
    for (istar = 0; istar < nStars; istar++)
	{
	semimaj = randomReal() * semimajMax;
	ecc = randomReal();
	inc = randomReal() * 2.0 * pi;
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

vector<Body*> placeStarsGHZUniformDensity(double &totalmass, int &nStars, double &innerRadius,
	double &outerRadius, double &G, int &iseed)
    {

    /*
     * Written by dh4gan, 2/12/13
     * Place stars such that they form an annular GHZ
     * Periastra and apastra must fall within the GHZ annulus
     */

    int istar, success;
    double eccmax, incmax, dsemimaj;
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

    // This constrains maximum and minimum semimajor axes

    dsemimaj = outerRadius - innerRadius;

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

	semimaj = innerRadius + randomReal() * dsemimaj;

	success = 0;

	while (success == 0)
	    {
	    ecc = randomReal() * eccmax;

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


int readParameters(ifstream &inputStream, string &prefixString, int &nRuns,
	int &nStars, double &dt, double &tMax,int &iseed, string &GHZChoice, double &innerRadius, double &outerRadius)
    {

    string par, line;

    // check that the file exists as cpp does not check
    //and it just returns 0 if it doesn't exist
    if (inputStream == 0)
	{
	cout << "No Input file found" << endl;
	return -1;
	}

    // Then loop through each line using getline and then
    //assign to vectors
    while (getline(inputStream, line))
	{
	istringstream iss(line);
	iss >> par;

	cout << "Reading " << par << endl;

	if (par == "OutputPrefix")
	    {

	    iss >> prefixString;
	    cout << prefixString << endl;
	    }

	if (par == "nStars")
	    {
	    iss >> nStars;

	    }

	if (par == "nRuns")
	    {
	    iss >> nRuns;

	    }

	if (par == "GHZChoice")
	    {
	    iss >> GHZChoice;
	    }

	if (par == "rInner")
	    {
	    iss >> innerRadius;
	    }

	if (par == "rOuter")
	    {
	    iss >> outerRadius;
	    }

	if (par == "Timestep")
	    {
	    iss >> dt;

	    }

	if (par == "MaximumTime")
	    {
	    iss >> tMax;
	    }

	if (par == "RandomSeed")
	    {
	    iss >> iseed;

	    iseed = -abs(iseed);
	    }



	}
    inputStream.close();
    return 0;
    }
