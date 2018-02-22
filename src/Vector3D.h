/*
 * Vector3D.h
 *
 *  Created on: Mar 14, 2013
 *      Author: dhf
 *
 *      A simple 3D Vector Class (holds useful vector operations)
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include "math.h"
#include <vector>


class Vector3D{
    public:
	Vector3D();
	Vector3D(float x,float y,float z);
	virtual ~Vector3D();

	float dotProduct(Vector3D other);
	Vector3D scaleVector(float scale);
	float magVector();
	Vector3D unitVector();
	Vector3D addVector(Vector3D b);
	Vector3D addVector(Vector3D b, Vector3D c);
	Vector3D addVector(Vector3D b, Vector3D c, Vector3D d);
	Vector3D relativeVector(Vector3D b);
	Vector3D subtractVector(Vector3D b);
	Vector3D crossProduct(Vector3D other);

	void defineFromSpherical(float r, float theta, float phi);

	void rotateX(float angle);
	void rotateY(float angle);
	void rotateZ(float angle);

	void printVector();

	float elements [3];

};




#endif /* VECTOR3D_H_ */
