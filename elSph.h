#pragma once

#include "elObject.h"

inline void el_init(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb, Elvari& elvari) {

	elFluid.rest_volume_all().setConstant(restVolume);

	elFluid.rest_compressionRate_all().setOnes();
	elBound.rest_compressionRate_all().setOnes();

	elFluid.rest_density_all().setConstant(1000);
	elBound.rest_density_all().setConstant(1000);

	elFluid.vel_all().colwise() = Array3_f(0, -3, 0);
}

inline void el_updateBoundWeight(Elbound& elBound, Elneighb& elNeighb) {

	// [Step 1] clean memnory (left in previous time step) before new summation
	elBound.weight_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		// [Step 2] get uid and neighb info
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			// [Step 3] neighbour is fluid or bound?
			// here we only need boundary neighbours
			if (!isFluid(neighbs.uid[n])) {
				elBound.weight(i) += neighbs.W[n];
			}
		}
		elBound.weight(i) *= restVolume;
		elBound.weight(i) = 0.7 / elBound.weight(i); // Eqn.12 of Pressure Boundaries 
		
		elBound.rest_volume(i) = elBound.weight(i) * restVolume;
	}
}

inline void el_neighbourSearch(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	elNeighb.positioning(elFluid, elBound);
	elNeighb.updateRelations(elFluid, elBound);
}

inline void el_prepareSphAttribute(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

	/*** LOOP 1 ***/
	// Object: Calculate mass() from {rest_density(), rest_volume()}
	elFluid.mass_all() = elFluid.rest_volume_all().array() * elFluid.rest_density_all().array();
	elBound.mass_all() = elBound.rest_volume_all().array() * elBound.rest_density_all().array();

	/*** LOOP 2 ***/
	// Object: Calculate Psi(), sph_volume(), alpha_term1(), alpha_term2()
	elFluid.sph_Psi_all().setZero();
	elBound.sph_Psi_all().setZero();
	elFluid.alpha_term1_all().setZero();
	elFluid.alpha_term2_all().setZero();
	elBound.alpha_term2_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.sph_Psi(i) += elFluid.X(neighbs.uid[n]) * neighbs.W[n];
				elFluid.alpha_term1(i) += neighbs.grad_W_vec[n] * elFluid.X(neighbs.uid[n]);
				elFluid.alpha_term2(i) += pow(elFluid.X(neighbs.uid[n]), 2) * pow(neighbs.grad_W[n], 2) / elFluid.mass(neighbs.uid[n]);
			}
			else {
				elFluid.sph_Psi(i) += elBound.weight(bid(neighbs.uid[n])) * elFluid.X(i) * neighbs.W[n];
				elFluid.alpha_term1(i) += neighbs.grad_W_vec[n] * elBound.weight(bid(neighbs.uid[n])) * elFluid.X(i);
			}
		}
		elFluid.sph_volume(i) = elFluid.X(i) / elFluid.sph_Psi(i);
		elFluid.alpha(i) = (elFluid.alpha_term1(i).transpose().matrix() *
			elFluid.alpha_term1(i).matrix())(0, 0) / elFluid.mass(i) + elFluid.alpha_term2(i);
		if (elFluid.alpha(i) < 10e-6) { elFluid.alpha(i) = 10e-6; }
	}
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.sph_Psi(i) += elFluid.X(neighbs.uid[n]) * neighbs.W[n];
				elBound.alpha_term2(i) += pow(elFluid.X(neighbs.uid[n]), 2) * pow(neighbs.grad_W[n], 2) / elFluid.mass(neighbs.uid[n]); // Number Density
			}
			else {
				elBound.sph_Psi(i) += elBound.X(bid(neighbs.uid[n])) * neighbs.W[n];
			}
		}
		elBound.sph_volume(i) = elBound.X(i) / elBound.sph_Psi(i);
		elBound.alpha(i) = elBound.alpha_term2(i);
		if (elBound.alpha(i) < 10e-6) { elBound.alpha(i) = 10e-6; }
	}
}

inline void el_advection(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb, Elvari& elvari) {

	// Object: Compute vel_adv()
	elFluid.adv_vel_all() = elFluid.vel_all();
	elFluid.adv_acce_all().colwise() = elvari.gravity;
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			Array3_f Aij, xij;
			if (isFluid(neighbs.uid[n])) {
				Aij = elFluid.vel(i) - elFluid.vel(neighbs.uid[n]);
				xij = elFluid.pos(i) - elFluid.pos(neighbs.uid[n]);
				elFluid.adv_acce(i) += viscosity_nu * laplacian_W(Aij, xij, neighbs.dis[n],
					neighbs.grad_W_vec[n], elFluid.sph_volume(neighbs.uid[n]));
			}
			else {
				Aij = elFluid.vel(i) - elBound.vel(bid(neighbs.uid[n]));
				xij = elFluid.pos(i) - elBound.pos(bid(neighbs.uid[n]));
				elFluid.adv_acce(i) += viscosity_nu * laplacian_W(Aij, xij, neighbs.dis[n],
					neighbs.grad_W_vec[n], elBound.sph_volume(bid(neighbs.uid[n])));
			}
		}
		elFluid.adv_vel(i) += elFluid.adv_acce(i) * timeStep;
	}
}

inline void el_loop_adv_Psi(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	elFluid.adv_Psi_all().setZero();
	elBound.adv_Psi_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.adv_Psi(i) += ((elFluid.adv_vel(i) - elFluid.adv_vel(neighbs.uid[n])).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elFluid.X(neighbs.uid[n]);
			}
			else {
				elFluid.adv_Psi(i) += ((elFluid.adv_vel(i) - elBound.vel(bid(neighbs.uid[n]))).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elBound.weight(bid(neighbs.uid[n])) * elFluid.X(i);
			}
		}
	}
	elFluid.adv_Psi_all() *= timeStep;
	elFluid.adv_Psi_all() += elFluid.sph_Psi_all();
	elFluid.adv_Psi_all() -= elFluid.rest_Psi_all();
	elFluid.adv_Psi_all() = (elFluid.adv_Psi_all().array() < 0).select(0, elFluid.adv_Psi_all());
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.adv_Psi(i) += ((elBound.vel(i) - elFluid.adv_vel(neighbs.uid[n])).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elFluid.X(neighbs.uid[n]); // Number Density
			}
			else {
				elBound.adv_Psi(i) += ((elBound.vel(i) - elBound.vel(bid(neighbs.uid[n]))).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elBound.X(bid(neighbs.uid[n]));
			}
		}
	}
	elBound.adv_Psi_all() *= timeStep;
	elBound.adv_Psi_all() += elBound.sph_Psi_all();
	elBound.adv_Psi_all() -= elBound.rest_Psi_all();
	elBound.adv_Psi_all() = (elBound.adv_Psi_all().array() < 0).select(0, elBound.adv_Psi_all());
}

inline void el_loop_adv_Psi_changeRate(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	el_loop_adv_Psi(elFluid, elBound, elNeighb);
	elFluid.adv_Psi_changeRate_all() = elFluid.adv_Psi_all() * timeStep_1;
	elBound.adv_Psi_changeRate_all() = elBound.adv_Psi_all() * timeStep_1;
}

inline void el_incompressibleSolver(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

	val_f compressibleState = 1;
	val_i currentIteration = 0;

	/*** LOOP 1 ***/
	// Object: Compute adv_Psi()
	el_loop_adv_Psi(elFluid, elBound, elNeighb);

	while (compressibleState > incompressibleThreshold || currentIteration < minIncompressibleSolverIter)
	{
		if (currentIteration > maxIncompressibleSolverIter) { break; }
		currentIteration++;

		/*** LOOP 2 ***/
		// Object: Compute pressure force
		elFluid.pressure_force_all().setZero();
		elBound.pressure_force_all().setZero();
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			val_i& uid = elFluid.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi(i) * elFluid.X(neighbs.uid[n]) / elFluid.alpha(i))
							+ (elFluid.adv_Psi(neighbs.uid[n]) * elFluid.X(i) / elFluid.alpha(neighbs.uid[n])));
				}
				else {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi(i) * elBound.X(bid(neighbs.uid[n])) / elFluid.alpha(i))
							+ (elBound.adv_Psi(bid(neighbs.uid[n])) * elFluid.X(i) / elBound.alpha(bid(neighbs.uid[n]))));
				}
			}
		}
		elFluid.pressure_force_all() *= -timeStep2_1;
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			val_i& uid = elBound.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {}else {}
			}
		}
		elBound.pressure_force_all() *= -timeStep2_1;

		elFluid.pressure_force_all() *= timeStep;
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			elFluid.adv_vel(i) += elFluid.pressure_force(i).array() / elFluid.mass(i);
		}

		/*** LOOP 3 ***/
		// Object: Compute adv_Psi()
		el_loop_adv_Psi(elFluid, elBound, elNeighb);

		compressibleState = (elFluid.adv_Psi_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart;
		cout << "Iter incom: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
	elFluid.pos_all() += elFluid.vel_all() * timeStep;
}

inline void el_divergenceFreeSolver(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	val_f compressibleState = 1;
	val_i currentIteration = 0;

	/*** LOOP 1 ***/
	// Object: Compute adv_Psi_changeRate()
	el_loop_adv_Psi_changeRate(elFluid, elBound, elNeighb);

	while (compressibleState > divergenceFreeThreshold || currentIteration < minDivergenceFreeSolverIter)
	{
		if (currentIteration > maxDivergenceFreeSolverIter) { break; }
		currentIteration++;

		/*** LOOP 2 ***/
		// Object: Compute pressure force
		elFluid.pressure_force_all().setZero();
		elBound.pressure_force_all().setZero();
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			val_i& uid = elFluid.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi_changeRate(i) * elFluid.X(neighbs.uid[n]) / elFluid.alpha(i))
							+ (elFluid.adv_Psi_changeRate(neighbs.uid[n]) * elFluid.X(i) / elFluid.alpha(neighbs.uid[n])));
				}
				else {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi_changeRate(i) * elBound.X(bid(neighbs.uid[n])) / elFluid.alpha(i))
							+ (elBound.adv_Psi_changeRate(bid(neighbs.uid[n])) * elFluid.X(i) / elBound.alpha(bid(neighbs.uid[n]))));
				}
			}
		}
		elFluid.pressure_force_all() *= -timeStep_1;
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			val_i& uid = elBound.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {}
				else {}
			}
		}
		elBound.pressure_force_all() *= -timeStep_1;
		
		elFluid.pressure_force_all() *= timeStep;
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			elFluid.adv_vel(i) += elFluid.pressure_force(i).array() / elFluid.mass(i);
		}

		/*** LOOP 3 ***/
		// Object: Compute adv_Psi_changeRate()
		el_loop_adv_Psi_changeRate(elFluid, elBound, elNeighb);

		compressibleState = (elFluid.adv_Psi_changeRate_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart * timeStep;
		cout << "Iter div-f: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
}

inline void el_update_adv_Psi_pj(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

	/*** LOOP 1 ***/
	// Object: Compute devi_pi()
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		elFluid.devi_pi(i) = elFluid.pressure(i) * elFluid.devi_pi_tmp(i);
	}

	/*** LOOP 2 ***/
	// Object: Compute devi_pj()
	elFluid.devi_pj_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.devi_pj(i) += elFluid.X(neighbs.uid[n]) * elFluid.pressure(neighbs.uid[n])
					/ pow(elFluid.sph_Psi(neighbs.uid[n]), 2) * neighbs.grad_W_vec[n];
			}
			else {
				elFluid.devi_pj(i) += elBound.X(bid(neighbs.uid[n])) * elBound.pressure(bid(neighbs.uid[n]))
					/ pow(elBound.sph_Psi(bid(neighbs.uid[n])), 2) * neighbs.grad_W_vec[n];
			}
			elFluid.devi_pj(i) *= elFluid.X(i) / elFluid.mass(i);
		}
	}
	elFluid.devi_pj_all() *= -timeStep2;

#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.adv_Psi_pj(i) += elFluid.X(neighbs.uid[n]) * ((elFluid.devi_pj(i) - elFluid.devi_pi(neighbs.uid[n])
					- elFluid.devi_pj(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
			else {
				elFluid.adv_Psi_pj(i) += elFluid.X(i) * elBound.weight(bid(neighbs.uid[n]))
					* (elFluid.devi_pj(i).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
		}
		elFluid.adv_Psi_pj(i) += timeStep2 * elFluid.pressure(i) * elFluid.X(i) / pow(elFluid.sph_Psi(i), 2) * elFluid.alpha_term2(i);
	}
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.adv_Psi_pj(i) += elFluid.X(neighbs.uid[n]) * ((-elFluid.devi_pi(neighbs.uid[n])
					- elFluid.devi_pj(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
			else {

			}
		}
		elBound.adv_Psi_pj(i) += timeStep2 * elBound.pressure(i) * elBound.X(i) / pow(elBound.sph_Psi(i), 2) * elBound.alpha_term2(i);
	}
}
	
inline void el_incompressibleSolver_II(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	/*** LOOP 1 ***/
	// Object: Compute adv_Psi()
	el_loop_adv_Psi(elFluid, elBound, elNeighb);

	val_f compressibleState = 1;
	val_i currentIteration = 0;

	elFluid.pressure_all().setZero();
	elBound.pressure_all().setZero();
	elFluid.adv_Psi_pj_all().setZero();
	elBound.adv_Psi_pj_all().setZero();

	/*** LOOP 2 ***/
	// Object: devi_pi_tmp()
	elFluid.devi_pi_tmp_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.devi_pi_tmp(i) += elFluid.X(neighbs.uid[n]) * neighbs.grad_W_vec[n];
			}
			else {
				elFluid.devi_pi_tmp(i) += elBound.X(bid(neighbs.uid[n])) * neighbs.grad_W_vec[n];
			}
		}
	}
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		elFluid.devi_pi_tmp(i) *= elFluid.X(i) / elFluid.mass(i) / pow(elFluid.sph_Psi(i), 2);
	}
	elFluid.devi_pi_tmp_all() *= -timeStep2;

	elFluid.vel_all() = elFluid.adv_vel_all();
	/*** IMPLICIT ITERATION ***/
	while (compressibleState > incompressibleThreshold_ii || currentIteration < minIncompressibleSolverIter)
	{
		if (currentIteration > maxIncompressibleSolverIter) { break; }
		currentIteration++;

		elFluid.adv_vel_all() = elFluid.vel_all();

		/*** LOOP 3 ***/
		// Object: Compute pressure()
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			elFluid.pressure(i) = (1 - relaxingFactor) * elFluid.pressure(i)
				+ relaxingFactor * pow(elFluid.sph_Psi(i), 2) * timeStep2_1 / elFluid.X(i)
				/ elFluid.alpha(i) * (elFluid.adv_Psi(i) + elFluid.adv_Psi_pj(i));
		}
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			elBound.pressure(i) = (1 - relaxingFactor) * elBound.pressure(i)
				+ relaxingFactor * pow(elBound.sph_Psi(i), 2) * timeStep2_1 / elBound.X(i)
				/ elBound.alpha(i) * (elBound.adv_Psi(i) + elBound.adv_Psi_pj(i));
		}
		elFluid.pressure_all() = (elFluid.pressure_all().array() < 0).select(0, elFluid.pressure_all());
		elBound.pressure_all() = (elBound.pressure_all().array() < 0).select(0, elBound.pressure_all());
		
		/*** LOOP 4 ***/
		// Object: Compute adv_Psi_pj()
		elFluid.adv_Psi_pj_all().setZero();
		elBound.adv_Psi_pj_all().setZero();
		el_update_adv_Psi_pj(elFluid, elBound, elNeighb);

		/*** LOOP 5 ***/
		// Object: Compute pressure_force()
		elFluid.pressure_force_all().setZero();
		elBound.pressure_force_all().setZero();
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			val_i& uid = elFluid.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {
					elFluid.pressure_force(i) += -elFluid.X(i) * elFluid.X(neighbs.uid[n]) *
						(elFluid.pressure(i) / pow(elFluid.sph_Psi(i), 2) + elFluid.pressure(neighbs.uid[n]) / pow(elFluid.sph_Psi(neighbs.uid[n]), 2))
						* neighbs.grad_W_vec[n];
				}
				else {
					elFluid.pressure_force(i) += -elFluid.X(i) * elBound.X(bid(neighbs.uid[n])) *
						(elFluid.pressure(i) / pow(elFluid.sph_Psi(i), 2) + elBound.pressure(bid(neighbs.uid[n])) / pow(elBound.sph_Psi(bid(neighbs.uid[n])), 2))
						* neighbs.grad_W_vec[n];
				}
			}
		}
		/*** LOOP 6 ***/
		// Object: Compute vel_adv()
		elFluid.pressure_force_all() *= timeStep;
#pragma omp parallel for
		for (val_i i = 0; i < elFluid.numPart; i++) {
			elFluid.adv_vel(i) += elFluid.pressure_force(i).array() / elFluid.mass(i);
		}
		el_loop_adv_Psi(elFluid, elBound, elNeighb);
		compressibleState = (elFluid.adv_Psi_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart;
		cout << "Iter incom: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
	elFluid.pos_all() += elFluid.vel_all() * timeStep;
}