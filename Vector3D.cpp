/*
 * Vector3D.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: dhf
 */


#include "Vector3D.h"
#include <iostream>

using namespace std;

Vector3D::Vector3D()
    {
    elements[0] = elements[1] = elements[2] = 0.0;
    }

Vector3D::Vector3D(float x, float y, float z)
	    {
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	    }


Vector3D::~Vector3D()
    {
    }

float Vector3D::dotProduct(Vector3D other)
// Written by dh4gan 11/3/13
// Simple method to take the dotProduct of two vectors
    {
    float dotProduct = 0.0;
    for (int i =0; i< 3; i++)
	{
	dotProduct += elements[i] *other.elements[i];
	}
    return dotProduct;
    }

Vector3D Vector3D::scaleVector(float scale)
// Written by dh4gan 11/3/13
// Simple method to scale a vector a by float b
    {
    Vector3D scaleVector;
    for (int i=0; i<3; i++)
	{
	scaleVector.elements[i] = elements[i] * scale;
	}

    return scaleVector;
    }

float Vector3D::magVector()
// Written by dh4gan 11/3/13
// Simple method to take the magnitude of a vector
    {

    float mag=0.0;

    for (int i=0; i<3; i++ )
	{
	mag += elements[i] * elements[i];
	}
    mag = sqrt(mag);
    return mag;
    }

Vector3D Vector3D::unitVector()
// Written by dh4gan 12/6/13
// Returns the unit vector
    {
    float mag = magVector();
    Vector3D unit;
    for (int i=0; i<3; i++)
	{
	unit.elements[i] = elements[i]/mag;
	}
    return unit;
    }

Vector3D Vector3D::addVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i];
	}

    return add;
    }

Vector3D Vector3D::addVector(Vector3D b, Vector3D c)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i] + c.elements[i];
	}

    return add;
    }

Vector3D Vector3D::addVector(Vector3D b, Vector3D c, Vector3D d)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i] + c.elements[i] + d.elements[i];
	}

    return add;
    }


Vector3D Vector3D::relativeVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to find the relative vector b-a
    {

    Vector3D relative;

    for (int i=0; i<3; i++)
	{
	relative.elements[i] = b.elements[i] - elements[i];
	}

    return relative;
    }

Vector3D Vector3D::subtractVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to find the relative vector b-a
    {

    Vector3D relative;

    for (int i=0; i<3; i++)
	{
	relative.elements[i] = elements[i] - b.elements[i];
	}

    return relative;
    }

Vector3D Vector3D::crossProduct(Vector3D other)
    {
    // Written by dh4gan 16/3/13
    // Calculates the cross product between a and b

    Vector3D cross;

    cross.elements[0] = elements[1] * other.elements[2] - elements[2] * other.elements[1];
    cross.elements[1] = elements[2] * other.elements[0] - elements[0] * other.elements[2];
    cross.elements[2] = elements[0] * other.elements[1] - elements[1] * other.elements[0];

    return cross;

    }

void Vector3D::defineFromSpherical(float r, float theta, float phi)
    {
    // Written by dh4gan 13/8/14
    // Defines vector elements using spherical coordinates

    elements[0] = r*sin(theta)*cos(phi);
    elements[1] = r*sin(theta)*sin(phi);
    elements[2] = r*cos(theta);

    }


void Vector3D::rotateX(float angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the x axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0];
    elements[1] = oldvec.elements[1]*cos(angle) - oldvec.elements[2]*sin(angle);
    elements[2] = oldvec.elements[1]*sin(angle) + oldvec.elements[2]*cos(angle);

    }

void Vector3D::rotateY(float angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the y axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0]*cos(angle) + oldvec.elements[2]*sin(angle);
    elements[1] = oldvec.elements[1];
    elements[2] = -oldvec.elements[0]*sin(angle) + oldvec.elements[2]*cos(angle);

    }

void Vector3D::rotateZ(float angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the z axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0]*cos(angle) - oldvec.elements[1]*sin(angle);
    elements[1] = oldvec.elements[0]*sin(angle) + oldvec.elements[1]*cos(angle);
    elements[2] = oldvec.elements[2];

    }

void Vector3D::printVector()
    {

    cout << elements[0] << "  " << elements[1] << "   " << elements[2] << endl;

    }

