#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <m3Vector.h>
#include <vector>

class Particle
{
public:
	m3Vector 	pos;			// Position of the particle
	m3Vector 	vel;			// Velocity of the particle
	m3Vector	predicted_vel;	// Predicted velocity: calculated with all forces but viscoelastic and pressure ones
	m3Vector	inter_vel;		// Intermediate velocity of the particle
	m3Vector	corrected_vel;	// Corrected velocity using SM
	m3Vector 	acc;			// Acceleration of the particle
	float		mass;

	m3Vector 	mOriginalPos;	// Original positions of the mesh points
	m3Vector 	mGoalPos;		// Goal positions
	bool		mFixed;			// Whether de particle is fixed in place

	float 		dens;			// density
	float 		pres;			// pressure
};

class Cell
{
public:
	std::vector<Particle*> contained_particles;
};

#endif
