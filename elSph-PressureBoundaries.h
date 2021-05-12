#pragma once

#include "elObject.h"

/* content intro */
// Common
	// el_PB_prepareSphAttribute(elFluid, elBound, elNeighb);
// IISPH
	// el_PB_incompressibleSolver_II(elFluid, elBound, elNeighb);
// DFSPH
	// el_PB_incompressibleSolver(elFluid, elBound, elNeighb);
	// el_PB_divergenceFreeSolver(elFluid, elBound, elNeighb);


inline void el_PB_prepareSphAttribute(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

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
				elFluid.sph_Psi(i) += elFluid.rest_volume(neighbs.uid[n]) * neighbs.W[n];
			}
			else {
				elFluid.sph_Psi(i) += elBound.rest_volume(bid(neighbs.uid[n])) * neighbs.W[n];
			}
		}
		elFluid.sph_volume(i) = elFluid.rest_volume(i) / elFluid.sph_Psi(i);
	}
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.sph_Psi(i) += elFluid.rest_volume(neighbs.uid[n]) * neighbs.W[n];
			}
			else {
				elBound.sph_Psi(i) += elBound.rest_volume(bid(neighbs.uid[n])) * neighbs.W[n];
			}
		}
		elBound.sph_volume(i) = elBound.rest_volume(i) / elBound.sph_Psi(i);
	}

#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.alpha_term1(i) += neighbs.grad_W_vec[n] * elFluid.sph_volume(neighbs.uid[n]);
				elFluid.alpha_term2(i) += pow(elFluid.sph_volume(neighbs.uid[n]), 2) * pow(neighbs.grad_W[n], 2) / elFluid.mass(neighbs.uid[n]);
			}
			else {
				elFluid.alpha_term1(i) += neighbs.grad_W_vec[n] * elBound.sph_volume(bid(neighbs.uid[n]));
			}
		}
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
				elBound.alpha_term2(i) += pow(elFluid.sph_volume(neighbs.uid[n]), 2) * pow(neighbs.grad_W[n], 2) / elFluid.mass(neighbs.uid[n]); // Number Density
			}
			else {
			}
		}
		elBound.alpha(i) = elBound.alpha_term2(i);
		if (elBound.alpha(i) < 10e-6) { elBound.alpha(i) = 10e-6; }
	}
}


inline void el_PB_loop_adv_Psi(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	elFluid.adv_Psi_all().setZero();
	elBound.adv_Psi_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.adv_Psi(i) += ((elFluid.adv_vel(i) - elFluid.adv_vel(neighbs.uid[n])).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elFluid.sph_volume(neighbs.uid[n]);
			}
			else {
				elFluid.adv_Psi(i) += ((elFluid.adv_vel(i) - elBound.vel(bid(neighbs.uid[n]))).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elBound.sph_volume(bid(neighbs.uid[n]));
			}
		}
	}
	elFluid.adv_Psi_all() *= timeStep;
	elFluid.adv_Psi_all() += elFluid.sph_Psi_all();
	elFluid.adv_Psi_all() -= elFluid.rest_compressionRate_all();
	elFluid.adv_Psi_all() = (elFluid.adv_Psi_all().array() < 0).select(0, elFluid.adv_Psi_all());
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.adv_Psi(i) += ((elBound.vel(i) - elFluid.adv_vel(neighbs.uid[n])).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elFluid.sph_volume(neighbs.uid[n]); // Number Density
			}
			else {
				elBound.adv_Psi(i) += ((elBound.vel(i) - elBound.vel(bid(neighbs.uid[n]))).transpose().matrix() *
					neighbs.grad_W_vec[n].matrix())(0, 0) * elBound.sph_volume(bid(neighbs.uid[n]));
			}
		}
	}
	elBound.adv_Psi_all() *= timeStep;
	elBound.adv_Psi_all() += elBound.sph_Psi_all();
	elBound.adv_Psi_all() -= (elBound.rest_compressionRate_all() * elBound.weight_all());
	elBound.adv_Psi_all() = (elBound.adv_Psi_all().array() < 0).select(0, elBound.adv_Psi_all());
}

inline void el_PB_loop_adv_Psi_changeRate(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	el_PB_loop_adv_Psi(elFluid, elBound, elNeighb);
	elFluid.adv_Psi_changeRate_all() = elFluid.adv_Psi_all() * timeStep_1;
	elBound.adv_Psi_changeRate_all() = elBound.adv_Psi_all() * timeStep_1;
}

inline void el_PB_update_adv_Psi_pj(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

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
				elFluid.devi_pj(i) += elFluid.sph_volume(neighbs.uid[n]) * elFluid.pressure(neighbs.uid[n]) * neighbs.grad_W_vec[n];
			}
			else {
				elFluid.devi_pj(i) += elBound.sph_volume(bid(neighbs.uid[n])) * elBound.pressure(bid(neighbs.uid[n])) * neighbs.grad_W_vec[n];
			}
			elFluid.devi_pj(i) *= elFluid.sph_volume(i) / elFluid.mass(i);
		}
	}
	elFluid.devi_pj_all() *= -timeStep2;

#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		val_i& uid = elFluid.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elFluid.adv_Psi_pj(i) += elFluid.sph_volume(neighbs.uid[n]) * ((elFluid.devi_pj(i) - elFluid.devi_pi(neighbs.uid[n])
					- elFluid.devi_pj(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
			else {
				elFluid.adv_Psi_pj(i) += elBound.sph_volume(bid(neighbs.uid[n]))
					* (elFluid.devi_pj(i).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
		}
		elFluid.adv_Psi_pj(i) += timeStep2 * elFluid.pressure(i) * elFluid.sph_volume(i) * elFluid.alpha_term2(i);
	}
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		val_i& uid = elBound.uid(i);
		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
		for (val_i n = 0; n < neighbs.neighbNum; n++) {
			if (isFluid(neighbs.uid[n])) {
				elBound.adv_Psi_pj(i) += elFluid.sph_volume(neighbs.uid[n]) * ((-elFluid.devi_pi(neighbs.uid[n])
					- elFluid.devi_pj(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0);
			}
			else {

			}
		}
		elBound.adv_Psi_pj(i) += timeStep2 * elBound.pressure(i) * elBound.sph_volume(i) * elBound.alpha_term2(i);
	}
}

inline void el_PB_incompressibleSolver_II(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	/*** LOOP 1 ***/
	// Object: Compute adv_Psi()
	el_PB_loop_adv_Psi(elFluid, elBound, elNeighb);

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
				elFluid.devi_pi_tmp(i) += elFluid.sph_volume(neighbs.uid[n]) * neighbs.grad_W_vec[n];
			}
			else {
				elFluid.devi_pi_tmp(i) += elBound.sph_volume(bid(neighbs.uid[n])) * neighbs.grad_W_vec[n];
			}
		}
	}
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		elFluid.devi_pi_tmp(i) *= elFluid.sph_volume(i) / elFluid.mass(i);
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
				+ relaxingFactor * timeStep2_1 / elFluid.sph_volume(i)
				/ elFluid.alpha(i) * (elFluid.adv_Psi(i) + elFluid.adv_Psi_pj(i));
		}
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			elBound.pressure(i) = (1 - relaxingFactor) * elBound.pressure(i)
				+ relaxingFactor * timeStep2_1 / elBound.sph_volume(i)
				/ elBound.alpha(i) * (elBound.adv_Psi(i) + elBound.adv_Psi_pj(i));
		}
		elFluid.pressure_all() = (elFluid.pressure_all().array() < 0).select(0, elFluid.pressure_all());
		elBound.pressure_all() = (elBound.pressure_all().array() < 0).select(0, elBound.pressure_all());

		/*** LOOP 4 ***/
		// Object: Compute adv_Psi_pj()
		elFluid.adv_Psi_pj_all().setZero();
		elBound.adv_Psi_pj_all().setZero();
		el_PB_update_adv_Psi_pj(elFluid, elBound, elNeighb);

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
					elFluid.pressure_force(i) += -elFluid.sph_volume(i) * elFluid.sph_volume(neighbs.uid[n]) *
						(elFluid.pressure(i) + elFluid.pressure(neighbs.uid[n]))
						* neighbs.grad_W_vec[n];
				}
				else {
					elFluid.pressure_force(i) += -elFluid.sph_volume(i) * elBound.sph_volume(bid(neighbs.uid[n])) *
						(elFluid.pressure(i) + elBound.pressure(bid(neighbs.uid[n])))
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
		el_PB_loop_adv_Psi(elFluid, elBound, elNeighb);
		compressibleState = (elFluid.adv_Psi_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart;
		cout << "Iter incom: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
	elFluid.pos_all() += elFluid.vel_all() * timeStep;
}

inline void el_PB_incompressibleSolver(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {

	val_f compressibleState = 1;
	val_i currentIteration = 0;

	/*** LOOP 1 ***/
	// Object: Compute adv_Psi()
	el_PB_loop_adv_Psi(elFluid, elBound, elNeighb);

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
						* ((elFluid.adv_Psi(i) * elFluid.sph_volume(neighbs.uid[n]) / elFluid.alpha(i))
							+ (elFluid.adv_Psi(neighbs.uid[n]) * elFluid.sph_volume(i) / elFluid.alpha(neighbs.uid[n])));
				}
				else {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi(i) * elBound.sph_volume(bid(neighbs.uid[n])) / elFluid.alpha(i))
							+ (elBound.adv_Psi(bid(neighbs.uid[n])) * elFluid.sph_volume(i) / elBound.alpha(bid(neighbs.uid[n]))));
				}
			}
		}
		elFluid.pressure_force_all() *= -timeStep2_1;
#pragma omp parallel for
		for (val_i i = 0; i < elBound.numPart; i++) {
			val_i& uid = elBound.uid(i);
			NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
			for (val_i n = 0; n < neighbs.neighbNum; n++) {
				if (isFluid(neighbs.uid[n])) {}
				else {}
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
		el_PB_loop_adv_Psi(elFluid, elBound, elNeighb);

		compressibleState = (elFluid.adv_Psi_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart;
		cout << "Iter incom: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
	elFluid.pos_all() += elFluid.vel_all() * timeStep;
}

inline void el_PB_divergenceFreeSolver(Elfluid& elFluid, Elbound& elBound, Elneighb& elNeighb) {
	val_f compressibleState = 1;
	val_i currentIteration = 0;

	/*** LOOP 1 ***/
	// Object: Compute adv_Psi_changeRate()
	el_PB_loop_adv_Psi_changeRate(elFluid, elBound, elNeighb);

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
						* ((elFluid.adv_Psi_changeRate(i) * elFluid.sph_volume(neighbs.uid[n]) / elFluid.alpha(i))
							+ (elFluid.adv_Psi_changeRate(neighbs.uid[n]) * elFluid.sph_volume(i) / elFluid.alpha(neighbs.uid[n])));
				}
				else {
					elFluid.pressure_force(i) += neighbs.grad_W_vec[n]
						* ((elFluid.adv_Psi_changeRate(i) * elBound.sph_volume(bid(neighbs.uid[n])) / elFluid.alpha(i))
							+ (elBound.adv_Psi_changeRate(bid(neighbs.uid[n])) * elFluid.sph_volume(i) / elBound.alpha(bid(neighbs.uid[n]))));
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
		el_PB_loop_adv_Psi_changeRate(elFluid, elBound, elNeighb);

		compressibleState = (elFluid.adv_Psi_changeRate_all() / elFluid.rest_Psi_all()).sum() / elFluid.numPart * timeStep;
		cout << "Iter div-f: " << currentIteration << endl;
		cout << "compressibleState: " << compressibleState << endl;
	}
	elFluid.vel_all() = elFluid.adv_vel_all();
}