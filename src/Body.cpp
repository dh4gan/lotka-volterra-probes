/*
 * Body.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "Body.h"
#include <iostream> // Included for debug lines only
#include <math.h>

Body::Body()
    {
    name = "Body";
    type = "Body";
    mass = 1.0;
    radius = 1.0;
    collisionBounce = true;

    Vector3D zero;
    position = zero;
    velocity = zero;
    acceleration = zero;
    jerk = zero;
    snap = zero;
    crackle = zero;

    semiMajorAxis = 0.0;

    eccentricityVector = zero; // Eccentricity Vector
    eccentricity = 0.0; // Eccentricity = Absolute Magnitude of Eccentricity Vector

    orbitalAngularMomentum = zero;
    magOrbitalAngularMomentum = 0.0;

    inclination = 0.0;
    trueAnomaly = 0.0;
    meanAnomaly = 0.0;
    eccentricAnomaly = 0.0;
    longitudePeriapsis = 0.0;
    longitudeAscendingNode = 0.0;
    argumentPeriapsis = 0.0;

    timestep =0.0;

    }

Body::Body(string &namestring, string &typestring, double &m, double &rad, Vector3D &pos,
	Vector3D &vel)
    {
    Vector3D zero;
    name = namestring;
    mass = m;
    type = typestring;
    radius = rad;
    collisionBounce = true;

    position = pos;
    velocity = vel;

    acceleration = zero;
    jerk = zero;
    snap = zero;
    crackle = zero;

    semiMajorAxis = 0.0;

    eccentricityVector = zero;
    eccentricity = 0.0;

    orbitalAngularMomentum = zero;
    magOrbitalAngularMomentum = 0.0;

    inclination = 0.0;
    trueAnomaly = 0.0;
    meanAnomaly = 0.0;
    eccentricAnomaly = 0.0;

    longitudePeriapsis = 0.0;
    argumentPeriapsis = 0.0;
    longitudeAscendingNode = 0.0;

    timestep = 0.0;

    }


Body::Body(string &namestring, string &typestring, double &m, double &rad, double semimaj, double ecc, double inc, double longascend,
				double argper, double meananom, double G, double totalMass)
    {

    Vector3D zero;

    name = namestring;
    type = typestring;
    mass = m;
    radius = rad;
    collisionBounce = true;

    semiMajorAxis = semimaj;

    eccentricityVector = zero; // Eccentricity Vector
    eccentricity = ecc; // Eccentricity = Absolute Magnitude of Eccentricity Vector

    orbitalAngularMomentum = zero;
    magOrbitalAngularMomentum = 0.0;

    inclination = inc;

    position = zero;
    velocity = zero;
    meanAnomaly = meananom;
    eccentricAnomaly = 0.0;

    argumentPeriapsis = argper;

    longitudeAscendingNode = longascend;
    longitudePeriapsis = argper + longitudeAscendingNode;

    calcTrueAnomaly();
    calcVectorFromOrbit(G, totalMass);
  }

Body::~Body()
    {
    }

/* Calculation Methods */


double Body::calcSeparation(Body* other)
    {
    /*
     * Written 21/8/14 by dh4gan
     * Calculates the separation between the Body and another Body
     *
     */

    return getPosition().relativeVector(other->getPosition()).magVector();

    }

Vector3D Body::calcSeparationVector(Body* other)
    {
    /*
     * Written 21/8/14 by dh4gan
     * Calculates the separation vector between the Body and another Body
     * Origin is the Body object this method acts on (not the Body pointer argument)
     *
     */

    return getPosition().relativeVector(other->getPosition());

    }

void Body::calcOrbitalAngularMomentum()
    {

    /* Author: dh4gan
     * Calculates orbital angular momentum vector and its magnitude, given positions and velocities
     Note that this method assumes that the rotation axis is based at the origin
     */

    orbitalAngularMomentum = position.crossProduct(velocity);
    magOrbitalAngularMomentum = orbitalAngularMomentum.magVector();

    }

void Body::calcEccentricity(double G, double totmass)
    {
    /* Author: dh4gan
     * Calculates orbital eccentricity vector
     * Assumes orbitalAngularMomentum is up to date
     */

    int i;
    double magvel, magpos;
    double gravparam, vdotr;
    Vector3D zerovector;


    gravparam = G * totmass;
    magpos = position.magVector();
    magvel = velocity.magVector();
    vdotr = velocity.dotProduct(position);

    if (magpos == 0.0)
	{
	eccentricityVector = zerovector;
	}
    else
	{
	for (i = 0; i < 3; i++)
	    {
	    eccentricityVector.elements[i] = (magvel * magvel
		    * position.elements[i] - vdotr * velocity.elements[i])
		    / (gravparam) - position.elements[i] / magpos;
	    }

	}

    eccentricity = eccentricityVector.magVector();

    }

void Body::calcOrbitFromVector(double G, double totmass)
    {

    /* Author:dh4gan
     * Calculates all the orbital elements of the Body, using its position and velocity
     * Calls calcOrbitalAngularMomentum, calcEccentricity in the process
     */

    double pi = 3.141592658285;
    Vector3D nplane; // Vector pointing toward Ascending Node
    double edotR, ndotR, edotn, ndotV, rdotV;
    double magpos, nscalar;

    calcOrbitalAngularMomentum();
    calcEccentricity(G, totmass);

    // Calculate semi-major axis from |h| and |e|

    semiMajorAxis = magOrbitalAngularMomentum * magOrbitalAngularMomentum / (G
	    * (totmass) * (1 - eccentricity * eccentricity));

    // Calculate Orbital Inclination

    if (magOrbitalAngularMomentum > 0.0)
	{
	inclination = acos(orbitalAngularMomentum.elements[2]
		/ magOrbitalAngularMomentum);
	}
    else
	{
	inclination = 0.0;
	}

    // Calculate Longitude of the Ascending Node

    if (inclination == 0.0)
	{

	longitudeAscendingNode = 0.0;

	nplane.elements[0] = magOrbitalAngularMomentum;
	nplane.elements[1] = 0.0;
	nplane.elements[2] = 0.0;

	}
    else
	{

	nplane.elements[0] = -orbitalAngularMomentum.elements[2];
	nplane.elements[1] = orbitalAngularMomentum.elements[1];
	nplane.elements[2] = 0.0;

	nscalar = nplane.magVector();
	longitudeAscendingNode = acos(nplane.elements[0] / nscalar);

	if (nplane.elements[1] < 0.0)
	    {
	    longitudeAscendingNode = 2.0 * pi - longitudeAscendingNode;
	    }
	}

    // Calculate true anomaly

    magpos = position.magVector();

    // If orbit circular, no inclination, then use the position vector itself

    if (eccentricity == 0.0 and inclination == 0.0)
	{
	trueAnomaly = acos(position.elements[0] / magpos);
	if (velocity.elements[0] < 0.0)
	    {
	    trueAnomaly = 2.0 * pi - trueAnomaly;
	    }

	}

    // If orbit circular and inclination non-zero, then use the orbital plane vector
    else if (eccentricity == 0.0)
	{
	ndotR = nplane.dotProduct(position);
	ndotR = ndotR / (magpos * nscalar);

	ndotV = nplane.dotProduct(velocity);

	trueAnomaly = acos(ndotR);

	if (ndotV > 0.0)
	    {
	    trueAnomaly = 2.0 * pi - trueAnomaly;
	    }
	}
    // For non-circular orbits use the eccentricity vector
    else
	{
	edotR = eccentricityVector.dotProduct(position);
	edotR = edotR / (magpos * eccentricity);

	rdotV = velocity.dotProduct(position);

	trueAnomaly = acos(edotR);

	if (rdotV < 0.0)
	    {
	    trueAnomaly = 2.0 * pi - trueAnomaly;
	    }

	}

    // Finally, calculate the longitude of periapsis - first calculate the argument of periapsis

    edotn = eccentricityVector.dotProduct(nplane);
    edotn = edotn / (nscalar * eccentricity);

    argumentPeriapsis = acos(edotn);
    if (eccentricityVector.elements[2] < 0.0)
	{
	argumentPeriapsis = 2.0 * pi - argumentPeriapsis;
	}

    longitudePeriapsis = argumentPeriapsis + longitudeAscendingNode;

    }

double Body::calcPeriod(double G, double totalMass)
	    {
    double pi = 3.1415926585;
    double period;

    period = sqrt(4.0 * pi * pi * semiMajorAxis * semiMajorAxis * semiMajorAxis / (G
    	    * totalMass));
    return period;

	    }

void Body::calcEccentricAnomaly()
{

	int ncalc;
	double tolerance;
	double Eold, Enext, fE, fdashE;

// If orbit is circular, meanAnomaly and eccentricAnomaly are equal

    if(eccentricity==0.0)
	{
    eccentricAnomaly = meanAnomaly;
	}
    else
	{
    // calculate eccentric anomaly (Newton-Raphson iteration)

    Eold = eccentricAnomaly;

    fE = Eold - eccentricity * sin(Eold) - meanAnomaly;
    fdashE = 1 - eccentricity * cos(Eold);
    tolerance = 1.0e30;
    ncalc = 0;
    while (fabs(tolerance) > 1.0e-3 and ncalc < 100)
	{

	if (fdashE != 0.0)
	    {
	    Enext = Eold - fE / fdashE;
	    }
	else
	    {
	    Enext = Eold * 1.05;
	    }

	fE = Enext - eccentricity * sin(Enext) - meanAnomaly;
	fdashE = 1 - eccentricity * cos(Enext);

	tolerance = (Enext - Eold);
	Eold = Enext;
	ncalc++;
	}

    eccentricAnomaly = Enext;
	}

}

void Body::calcTrueAnomaly(double G, double totalMass, double time)
    {
    /* Author: dh4gan 7/8/13
     * Calculates the true anomaly of a Body, given the time and its orbital parameters
     * Assumes t=0 corresponds to true anomaly = 0
     *
     *
     */

    double period;

    double pi = 3.1415926585;


    // Calculate period

    period = calcPeriod(G,totalMass);

    // Calculate mean anomaly
    meanAnomaly = fmod(2.0*pi*time / period, 2.0*pi);

    // Calculate Eccentric Anomaly
    calcEccentricAnomaly();

    // Finally, calculate true anomaly

    trueAnomaly = 2.0*atan2(sqrt(1.0+eccentricity)*sin(eccentricAnomaly/2.0),sqrt(1.0-eccentricity)*cos(eccentricAnomaly/2.0));

    }

void Body::calcTrueAnomaly()
    {
    /* Author: dh4gan 7/8/13
     * Calculates the true anomaly of a Body, given the time and its orbital parameters
     * This overloaded method assumes mean anomaly already known
     *
     */

    calcEccentricAnomaly();

    // Finally, calculate true anomaly

    trueAnomaly = 2.0*atan2(sqrt(1.0+eccentricity)*sin(eccentricAnomaly/2.0),sqrt(1.0-eccentricity)*cos(eccentricAnomaly/2.0));

    }


void Body::calcVectorFromOrbit(double G, double totmass)
    {

    /* Author:dh4gan 1/8/13
     *
     * Calculates a Body's position and velocity, given its orbital elements (relative to CoM)
     * Uses separation and true anomaly to calculate position in frame coplanar with orbit
     * Then uses a series of rotations to give correct inclination and longitudes of periapsis and ascending nodes
     *
     */

    double magpos, magvel;
    double semiLatusRectum, gravparam;

    /* 1. calculate distance from CoM using semimajor axis, eccentricity and true anomaly*/

    magpos = semiMajorAxis * (1.0 - eccentricity * eccentricity) / (1.0
	    + eccentricity * cos(trueAnomaly));

    /* 2. Calculate position vector in orbital plane */

    position.elements[0] = magpos * cos(trueAnomaly);
    position.elements[1] = magpos * sin(trueAnomaly);
    position.elements[2] = 0.0;

    /* 3. Calculate velocity vector in orbital plane */
    semiLatusRectum = fabs(semiMajorAxis * (1.0 - eccentricity * eccentricity));
    gravparam = G * totmass;

    if (semiLatusRectum != 0.0)
	{
	magvel = sqrt(gravparam / semiLatusRectum);
	}
    else
	{
	magvel = 0.0;
	}

    velocity.elements[0] = -magvel * sin(trueAnomaly);
    velocity.elements[1] = magvel * (cos(trueAnomaly) + eccentricity);
    velocity.elements[2] = 0.0;

    /* 4. Begin rotations:
     * Firstly, Rotation around z axis by -argumentPeriapsis */

    if(argumentPeriapsis!=0.0)
	{
    position.rotateZ(-1 * argumentPeriapsis);
    velocity.rotateZ(-1 * argumentPeriapsis);
	}

    /* Secondly, Rotate around x by -inclination */

    if(inclination !=0.0)
	{
    position.rotateX(-1 * inclination);
    velocity.rotateX(-1 * inclination);
	}

    /* Lastly, Rotate around z by longitudeAscendingNode */

    if(longitudeAscendingNode !=0.0)
	{
    position.rotateZ(-1 * longitudeAscendingNode);
    velocity.rotateZ(-1 * longitudeAscendingNode);
	}

    }

void Body::moveAlongOrbit(double G, double totalMass, double dt)
{
	/* Author: dh4gan 8/11/13
	 * Given timestep dt, recalculates mean anomaly and updates position of the Body
	 */

	double pi = 3.1415926585;

	double period = calcPeriod(G,totalMass);
	meanAnomaly = meanAnomaly + (2.0*pi*dt/period);
	meanAnomaly = fmod(meanAnomaly, 2.0*pi);

	calcTrueAnomaly();
	calcVectorFromOrbit(G,totalMass);

}

void Body::calcTimestep(double greekEta)
    {
    /* Author: Joao Ferreira
     * Calculate the preferred timestep for the Body given its acceleration, jerk, snap and crackle
     *
     *
     * UPDATED : DRH 14/3/2013
     */

    double tolerance = 1e-20;

    double normJ = jerk.magVector();
    double normA = acceleration.magVector();
    double normC = crackle.magVector();
    double normS = snap.magVector();

    // If numerator zero, give a warning

    if (normA * normS + normJ * normJ < tolerance)
	{
	cout << "warning in calcTimestep: numerator zero for Body " << name
		<< endl;
	setTimestep(0.0);
	}

    // If denominator zero, give a warning - set timestep very large

    else if (normC * normJ + normS * normS < tolerance)
	{
	cout << "warning in calcTimestep: denominator zero for Body " << name
		<< endl;
	setTimestep(1.0e30);
	}

    // Otherwise calculate timestep as normal
    else
	{
	setTimestep(sqrt(greekEta * (normA * normS + normJ * normJ) / (normC
		* normJ + normS * normS)));
	}

    }

void Body::calcAccelJerk(double G, vector<Body*> bodyarray,
	double softening_length)
    {
    /* Author: David the legend Harvey
     * Calculate the gravitational acceleration on a Body given an array of **pointers** to bodies
     * Method :
     * 		Loop over each body in the system and calculate the vector r
     * 		from the body in question, then work out the force vector
     * 		and therefore the acceleration vector
     *
     * 		a_i = -G * M_i
     * 		    ------------- * r^hat_i
     * 		         |r_i|**2
     *
     * 		j_i = -G * M_i
     * 		    ------------- * v^hat_i - 3 alpha_i a_i
     * 		         |r_i|**2
     *
     *		alpha_i = (r_i.v_i)/|r_i|**2
     *
     *
     *		r^hat is the unit vector in the direction of the vector connecting
     *		M1 is the mass of the body
     *		M2 is the mass of the other body
     *		|r| is the distance between the two bodies
     *
     *
     * Argument :
     * 		bodyarray : a vector containing all the Body objects
     * 					in the system that the Body in question belongs to
     *
     */
    int b;
    int number_bodies = bodyarray.size();

    double alpha, factor, mass_b;
    double rmag, r2, r3, r2_1;

    Vector3D position_j;
    Vector3D velocity_j;
    Vector3D rel_velocity;
    Vector3D rel_position;

    Vector3D accelterm;
    Vector3D jerkterm;
    Vector3D jerkterm1;
    Vector3D jerkterm2;

    for (b = 0; b < number_bodies; b++)
	{
	//Get the position and velocity of the body in question

	position_j = bodyarray[b]->getPosition();
	velocity_j = bodyarray[b]->getVelocity();

	// Get relative position and velocity

	rel_position = position.relativeVector(position_j);
	rel_velocity = velocity.relativeVector(velocity_j);

	rmag = rel_position.magVector();

	// Since the body in question is inside the body array
	// I need to make sure the body I am looping through
	// isn't itself so skip the loop if distance ==0

	if (rmag < 1.0e-2 * softening_length)
	    {
	    continue;
	    }

	// Add gravitational softening
	r2 = rmag * rmag + softening_length * softening_length;
	r2_1 = 1.0 / r2;
	rmag = sqrt(r2);
	r3 = rmag * rmag * rmag;

	// Define this factor, as it's useful
	mass_b = bodyarray[b]->getMass();
	factor = -G * mass_b / r3;

	// Calculate acceleration term
	accelterm = rel_position.scaleVector(factor);

	// acceleration = acceleration - accelterm

	acceleration = accelterm.relativeVector(acceleration);

	// now jerk - calculate alpha term
	alpha = rel_velocity.dotProduct(rel_position);

	alpha = alpha * r2_1;
	jerkterm1 = rel_velocity.scaleVector(factor);
	jerkterm2 = accelterm.scaleVector(-3 * alpha);

	// jerkterm = jerkterm2 - jerkterm1
	jerkterm = jerkterm1.addVector(jerkterm2);

	// jerk = jerk - jerkterm

	jerk = jerkterm.relativeVector(jerk);

	} // End of loop
    // End of method

    }

void Body::calcSnapCrackle(double G, vector<Body*> bodyarray,
	double softening_length)
    {
    /* Author: dh4gan 11/3/13
     * Calculate the snap and crackle on a Body given an array of **pointers** to bodies
     * WARNING: This method does not work unless calcAccelJerk is called prior to this
     * This method will calculate a lot of the same terms as calcAccelJerk, but
     * this is the only way these calculations can be done on a Body by Body basis
     */

    int b;
    int number_bodies = bodyarray.size();

    double alpha, beta, gamma, factor;
    double rmag, r2, r3, r2_1;
    double vmag, v2;

    Vector3D zerovector;
    Vector3D position_j;
    Vector3D velocity_j;
    Vector3D acceleration_j;
    Vector3D jerk_j;

    Vector3D rel_velocity;
    Vector3D rel_position;
    Vector3D rel_acceleration;
    Vector3D rel_jerk;

    Vector3D accelterm;

    Vector3D jerkterm;
    Vector3D jerkterm1;
    Vector3D jerkterm2;

    Vector3D snapterm;
    Vector3D snapterm1;
    Vector3D snapterm2;
    Vector3D snapterm3;

    Vector3D crackleterm;
    Vector3D crackleterm1;
    Vector3D crackleterm2;
    Vector3D crackleterm3;
    Vector3D crackleterm4;

    for (b = 0; b < number_bodies; b++)
	{
	//Get the position and velocity of the body in question

	position_j = bodyarray[b]->getPosition();
	velocity_j = bodyarray[b]->getVelocity();
	acceleration_j = bodyarray[b]->getAcceleration();
	jerk_j = bodyarray[b]->getJerk();

	//bodyarray[b]->setSnap(zerovector);
	//bodyarray[b]->setCrackle(zerovector);

	// Get relative position, velocity, acceleration and jerk

	rel_position = position.relativeVector(position_j);
	rel_velocity = velocity.relativeVector(velocity_j);
	rel_acceleration = acceleration.relativeVector(acceleration_j);
	rel_jerk = jerk.relativeVector(jerk_j);

	rmag = rel_position.magVector();
	vmag = rel_velocity.magVector();

	v2 = vmag * vmag;

	//Since the body in question is inside the body array
	//I need to make sure the body I am looping through
	//isn't itself so skip the loop if distance ==0

	if (rmag < 1.0e-2 * softening_length)
	    {
	    continue;
	    }

	// Add gravitational softening
	r2 = rmag * rmag + softening_length * softening_length;
	r2_1 = 1.0 / r2;
	rmag = sqrt(r2);
	r3 = rmag * rmag * rmag;

	// Define this factor, as it's useful
	factor = G * bodyarray[b]->getMass() / r3;

	// Calculate acceleration term
	accelterm = rel_position.scaleVector(factor);

	// now jerk - calculate alpha term
	alpha = rel_velocity.dotProduct(rel_position) * r2_1;

	//alpha = alpha*r2_1;
	jerkterm1 = rel_velocity.scaleVector(factor);
	jerkterm2 = accelterm.scaleVector(3 * alpha);

	// jerkterm = jerkterm2 - jerkterm1
	jerkterm = jerkterm2.relativeVector(jerkterm1);

	// calculate snap terms
	beta = (v2 + rel_position.dotProduct(rel_acceleration)) * r2_1 + alpha
		* alpha;

	snapterm1 = rel_acceleration.scaleVector(factor);
	snapterm2 = jerkterm.scaleVector(-6 * alpha);
	snapterm3 = accelterm.scaleVector(-3 * beta);

	snapterm = snapterm1.addVector(snapterm2, snapterm3);

	snap = snapterm.relativeVector(snap);

	// Finally crackle terms

	gamma = (3.0 * rel_velocity.dotProduct(rel_acceleration)
		+ rel_position.dotProduct(rel_jerk)) * r2_1;
	gamma = gamma + alpha * (3.0 * beta - 4.0 * alpha * alpha);

	crackleterm1 = rel_jerk.scaleVector(factor);
	crackleterm2 = snapterm.scaleVector(-9.0 * alpha);
	crackleterm3 = jerkterm.scaleVector(-9.0 * beta);
	crackleterm4 = accelterm.scaleVector(-3.0 * gamma);

	crackleterm = crackleterm1.addVector(crackleterm2, crackleterm3,
		crackleterm4);

	crackle = crackleterm.relativeVector(crackle);

	} // End of loop
    // End of method
    }

