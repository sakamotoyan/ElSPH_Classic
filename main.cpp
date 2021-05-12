/*** TO RUN THIS PARGRAMME, USE "g++-9 -fopenmp main.cpp -o main" TO COMPILE ***/

#include "elSph.h"

void test()
{
	Elfluid elFluid(MaxPartNum_fluid, fluidGridConfig);
	Elbound elBound(MaxPartNum_bound, boundGridConfig);
	Elneighb elNeighb;
	Elvari elvari;

	FloatBlock refPosBlock_fluid, refPosBlock_bound;
	val_i a = 5;
	val_i b = 5;
	Array3_f sizeCuboid_fluid(diamPart * a, diamPart * a, diamPart * a);
	Array3_f sizeCuboid_bound(diamPart * b, diamPart * b, diamPart * b);
	Array3_f posCenter_fluid(0, 0, 0);
	Array3_f posCenter_bound(0, -diamPart * a, 0);
	refPosBlock_fluid.asPosBlock(Array3_f(diamPart), posCenter_fluid, sizeCuboid_fluid);
	refPosBlock_bound.asPosBlock(Array3_f(diamPart), posCenter_bound, sizeCuboid_bound);

	// [Step 2] import data
	// [Step 2.1] get the number of particles of fluid and boundary 
	// [note] upper limits of numPart is in globaValue.h
	elFluid.numPart = refPosBlock_fluid.getNumMembers();
	elBound.numPart = refPosBlock_bound.getNumMembers();
	// [Step 2.2] feed data into elFluid and elBound respectively
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		elFluid.pos(i) = refPosBlock_fluid.dataMat.col(i);
	}
#pragma omp parallel for
	for (val_i i = 0; i < elBound.numPart; i++) {
		elBound.pos(i) = refPosBlock_bound.dataMat.col(i);
	}

	// [Step 3] neighbour search
	// [Step 3.1] feed position data into background grid
	elNeighb.positioning(elFluid, elBound);
	// [Step 3.2] build up the reference struture for all particles
	elNeighb.updateRelations(elFluid, elBound);

	// [step 4] SPH step
	el_init(elFluid, elBound, elNeighb, elvari);
	el_updateBoundWeight(elBound, elNeighb);
	el_prepareSphAttribute(elFluid, elBound, elNeighb);


	/*cout << "before: " << endl << elFluid.adv_vel_all() << endl;*/
	el_advection(elFluid, elBound, elNeighb, elvari);
	//cout << "before: " << endl << elFluid.adv_vel_all() << endl;
	//el_incompressibleSolver(elFluid, elBound, elNeighb);
	//el_divergenceFreeSolver(elFluid, elBound, elNeighb);
	el_incompressibleSolver_II(elFluid, elBound, elNeighb);
	PAUSE;
	
	//cout << "adv_Psi_all: " << endl << elBound.adv_Psi_all() << endl; 

	
	//elFluid.vel_all();
	

	//cout << "elFluid.alpha_all:" << endl << elFluid.alpha_all() << endl;
	//cout << "elBound.alpha_all:" << endl << elBound.alpha_all() << endl;
	// elNeighb.loopNonZeros();
	// cout << "grid total:" << elNeighb.gridTotal << endl;

	// NeighbRelation& re = elNeighb.neighbRelations[elBound.uid(0)];
	// cout << "total neighbs " << re.neighbNum << endl;
	// for (int i = 0;i < re.neighbNum;i++) {
	// 	cout << "neighb " << i << ", uid is: " << re.uid[i] << ", dis is " << re.dis[i] << endl;
	// }


	// #pragma omp parallel for
	// 	for (int i = 0;i < elFluid.numPart;i++) {
	// 		val_i& uid = elFluid.uid(i);
	// 		NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
	// 		for (int n = 0;n < neighbs.neighbNum;n++) {
	// 			if (neighbs.dis[n] == 0)break;
	// 			Array3_f dis = elFluid.pos(i) - elFluid.pos(neighbs.uid[n]);
	// 			elFluid.gradTest(i) += (dis / neighbs.dis[n]) * neighbs.grad_W[n];
	// 		}
	// 	}


		// cout << "all grad: " << endl << elFluid.gradTest_all().colwise().sum() << endl;

		// PAUSE;
		// show things
		// elFluid.show(elFluid.numPart);
		// elBound.show(elBound.numPart);
		// cout << "all n nums: " << endl << elFluid.neighbNum_all() << endl;
		// elFluid.cleanSpecificIntAttributes(elFluid.neighbNum_i);
		// cout << "all n nums: " << endl << elFluid.neighbNum_all() << endl;
		// cout << "all n nums: " << endl << elBound.neighbNum_all() << endl;
		// elNeighb.loopNonZeros();
		// cout << "backGroundGrid row/col: " << elNeighb.backGroundGrid.rows() << "/" << elNeighb.backGroundGrid.cols() << endl;

	// cout << "test: " << elFluid.pos(0).transpose() << " / " << elFluid.pos(1).transpose()
	// 		<< " norm: " << (elFluid.pos(1) - elFluid.pos(0)).matrix().norm() << endl;
}

int main()
{
	CALL_TIME(test(), "programme");
	/*Array<val_f,Dynamic,Dynamic> a;
	a.resize(3, 1);
	a << 0, 1, 2;
	cout << "test " << a/2 << endl;*/
	// for (int i = 0;i < 0;i++) {
	// 	cout << "shouldn't happen!" << endl;
	// }
	return 0;
}

