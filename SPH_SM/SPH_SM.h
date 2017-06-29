/*
This class implements the SPH [5] and SM [2] algorithms with the needed modifications to correct velocity as presented in [1]. 

@author Octavio Navarro
@version 1.0
*/

#ifndef __SPH_SM_H__
#define __SPH_SM_H__

#include <m3Vector.h>
#include <m3Bounds.h>
#include <m3Real.h>
#include <m3Matrix.h>
#include <m9Matrix.h>

#include <Particle.h>

#include <vector>
#include <map>

#define PI 3.141592f
#define INF 1E-12f

class SPH_SM
{
	private:
	
		m3Real 		kernel;					// kernel or h in kernel function
		int 		Max_Number_Paticles;	// initial array for particles
		int 		Number_Particles;		// paticle number

		m3Vector 	Grid_Size;				// Size of a size of each grid voxel
		m3Vector 	World_Size;				// screen size
		m3Real 		Cell_Size;				// Size of the divisions in the grid; used to determine the cell position for the has grid; kernel or h
		int 		Number_Cells;			// Number of cells in the hash grid

		m3Vector 	Gravity;
		m3Real 		K;						// ideal pressure formulation k; Stiffness of the fluid
											// The lower the value, the stiffer the fluid
		m3Real 		Stand_Density;			// ideal pressure formulation p0
		m3Real 		Time_Delta;			
		m3Real 		Wall_Hit;				// To manage collisions with the environment.
		m3Real 		mu;						// Viscosity.

		// Smoothing kernel constants for SPH
		m3Real 		Poly6_constant, Spiky_constant, Visco_Constant;

		// SM Parameters
		m3Bounds bounds;					// Controls the bounds of the simulation.
		m3Real alpha;						// alpha[0...1] Simulates stiffness.
		m3Real beta;						// Not entirely sure about this parameter, but the larger, the less rigid the object.
		bool quadraticMatch;				// Linear transformations can only represent shear and stretch. To extend the range of motion by twist and bending modes, we move from linear to quadratic transformations.
		bool volumeConservation;			// Allows the object to conserve its volume.
		bool allowFlip;

		// Max velocity allowed for a particle.
		m3Vector 	max_vel = m3Vector(INF, INF, INF);	

		Particle*	Particles;
		Cell*		Cells;

	public:
		SPH_SM();
		~SPH_SM();

		void Init_Fluid(std::vector<m3Vector> positions);	// initialize fluid
		void Init_Particle(m3Vector pos, m3Vector vel);		// initialize particle system
		m3Vector Calculate_Cell_Position(m3Vector pos);		// get cell position
		int Calculate_Cell_Hash(m3Vector pos);				// get cell hash number
		void add_viscosity(m3Real);

		// Kernel function
		m3Real Poly6(m3Real r2);		// for density
		m3Real Spiky(m3Real r);			// for pressure
		m3Real Visco(m3Real);			// for viscosity

		/// Hashed the particles into a grid
		void Find_neighbors();

		/// SPH SM methods
		/// Calculates the predicted velocity, and the corrected velocity using SM, in order to
		/// obtain the intermediate velocity that is input to SPH. Taken from 2014 - A velocity correcting method
		/// for volume preserving viscoelastic fluids
		void calculate_velocity_correction();
		void projectPositions();

		/// Applies external forces for F-adv, including gravity
		void apply_external_forces(m3Vector* forcesArray, int* indexArray, int size);

		void calculate_corrected_velocity();
		void calculate_intermediate_velocity();

		void Compute_Density_SingPressure();
		void Compute_Force();
		void Update_Pos_Vel();

		void compute_SPH_SM();
		void Animation();

		inline int Get_Particle_Number() { return Number_Particles; }
		inline m3Vector Get_World_Size() { return World_Size; }
		inline Particle* Get_Paticles()	 { return Particles; }
		inline Cell* Get_Cells()		 { return Cells; }
};


#endif