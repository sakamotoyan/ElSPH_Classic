#pragma once

#include "elObject.h"

class ElPhase : public DataGrid
{
private:
    val_i viscosity_nu_1f[3] = { 1,0,1 };
    val_i restDensity_1f[3] = { 1,1,2 };
    val_i rgba_4f[3] = { 1,3,4 };
public:
    val_i phaseNum;
    inline ElPhase(val_i num = phaseNumber, const val_i* config = phaseGridConfig) : DataGrid(num, config) {
        phaseNum = phaseNumber;
    }

    inline val_f& viscosity_nu(val_i i)
    {
        val_i(&ref)[3] = viscosity_nu_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef viscosity_nu_all()
    {
        val_i(&ref)[3] = viscosity_nu_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& restDensity(val_i i)
    {
        val_i(&ref)[3] = restDensity_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef restDensity_all()
    {
        val_i(&ref)[3] = restDensity_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef rgba(val_i i)
    {
        val_i(&ref)[3] = rgba_4f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline val_f& colorA(val_i i)
    {
        return rgba(i)(3, 0);
    }
    inline FloatBlockRef rgba_all()
    {
        val_i(&ref)[3] = rgba_4f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
};

inline void el_mixColor(Elfluid& elFluid, ElPhase& elPhase) {
	elFluid.color_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
			val_f singelA = elPhase.colorA(ph) * elFluid.volumeFraction(i, ph);
			if (singelA < 10e-6)continue;
			elFluid.colorA(i) += singelA;
			elFluid.color(i).col(0).head(3) = (elFluid.color(i).col(0).head(3) * elFluid.colorA(i) * (val_f(1.0) - singelA)
				+ (elPhase.rgba(ph).col(0).head(3) * singelA)) / elFluid.colorA(i);
		}
	}
}

inline void driftToVolume(Elfluid& elFluid, Elneighb& elNeighb, ElPhase& elPhase) {
    bool hasNegative = true;
    elFluid.transact_all().setOnes();
    while (hasNegative) {
        hasNegative = false;
        elFluid.volumeFractionCache_all() = elFluid.volumeFraction_all();
#pragma omp parallel for
        for (val_i i = 0; i < elFluid.numPart; i++) {
            if (elFluid.transact(i) < 0)continue;
            for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
                val_i& uid = elFluid.uid(i);
                NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
                val_f tmpTran = 0;
                val_f tmpDiff = 0;
                for (val_i n = 0; n < neighbs.neighbNum; n++) {
                    if (isFluid(neighbs.uid[n])) {
                        if (elFluid.transact(neighbs.uid[n]) > 0) {
                            tmpTran += -elFluid.rest_volume(neighbs.uid[n]) * ((elFluid.volumeFraction(i, ph) * elFluid.drift_vel(i, ph)
                                + elFluid.volumeFraction(neighbs.uid[n], ph) * elFluid.drift_vel(neighbs.uid[n], ph)).transpose().matrix()
                                * neighbs.grad_W_vec[n].matrix())(0, 0);
                            tmpDiff += elFluid.rest_volume(neighbs.uid[n]) * (elFluid.volumeFraction(i, ph) - elFluid.volumeFraction(neighbs.uid[n], ph))
                                * ((elFluid.pos(i) - elFluid.pos(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0)
                                / (pow(neighbs.dis[n], 2) + 10e-6);
                            tmpTran *= transactionCoefficient;
                            tmpDiff *= diffusionCoefficinet;
                            tmpTran = timeStep * (tmpTran + tmpDiff);
                            elFluid.volumeFractionCache(i, ph) += tmpTran;
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for (val_i i = 0; i < elFluid.numPart; i++) {
            for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
                if (elFluid.volumeFractionCache(i, ph) < 0) {
                    hasNegative = true;
                    elFluid.transact(i) = -1;
                }
            }
        }
    }

    elFluid.volumeFraction_all() = elFluid.volumeFractionCache_all();
}

inline void el_multiphaseSolver(Elfluid& elFluid, Elneighb& elNeighb, ElPhase& elPhase, Elvari& elvari) {

    // Compute zeta & adv_acce
    elFluid.mt_zeta_all().setZero();
#pragma omp parallel for
	for (val_i i = 0; i < elFluid.numPart; i++) {
		for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
            elFluid.mt_zeta(i) += elFluid.volumeFraction(i, ph) * (elPhase.restDensity(ph) - elFluid.rest_density(i)) / elPhase.restDensity(ph);
		}
        if (elFluid.mt_zeta(i) == 1) { elFluid.mt_zeta(i) = 1 - 10e-6; }
        elFluid.intermediate_adv_acce(i) = (elFluid.adv_acce(i) - (elFluid.mt_zeta(i) * elvari.gravity)) / (1 - elFluid.mt_zeta(i));
	}

    // Compute intermediate adv_acce
#pragma omp parallel for
    for (val_i i = 0; i < elFluid.numPart; i++) {
        for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
            elFluid.drift_vel(i, ph) = timeStep
                * ((elFluid.intermediate_adv_acce(i) * elFluid.rest_density(i) / elPhase.restDensity(ph))
                    + (elvari.gravity * (elPhase.restDensity(ph) - elFluid.rest_density(i)) / elPhase.restDensity(ph))
                    - elFluid.adv_acce(i));
        }  
    }

    // Compute volumeFraction change
    driftToVolume(elFluid, elNeighb, elPhase);
}

inline void el_multiphaseSolver_SCA21(Elfluid& elFluid, Elneighb& elNeighb, ElPhase& elPhase, Elvari& elvari) {
#pragma omp parallel for
    for (val_i i = 0; i < elFluid.numPart; i++) {
        for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
            elFluid.drift_vel(i, ph) = timeStep
                * ((elPhase.restDensity(ph) - elFluid.rest_density(i)) / elFluid.rest_density(i))
                * (elvari.gravity - elFluid.adv_acce(i));
        }
    }

    bool hasNegative = true;
    elFluid.transact_all().setOnes();
    while (hasNegative) {
        hasNegative = false;
        elFluid.volumeFractionCache_all() = elFluid.volumeFraction_all();
#pragma omp parallel for
        for (val_i i = 0; i < elFluid.numPart; i++) {
            if (elFluid.transact(i) < 0)continue;
            for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
                val_i& uid = elFluid.uid(i);
                NeighbRelation& neighbs = elNeighb.neighbRelations[uid];
                val_f tmpTran = 0;
                val_f tmpDiff = 0;
                for (val_i n = 0; n < neighbs.neighbNum; n++) {
                    if (isFluid(neighbs.uid[n])) {
                        if (elFluid.transact(neighbs.uid[n]) > 0) {
                            tmpTran += -elFluid.sph_volume(neighbs.uid[n]) * ((elFluid.volumeFraction(i, ph) * elFluid.drift_vel(i, ph)
                                + elFluid.volumeFraction(neighbs.uid[n], ph) * elFluid.drift_vel(neighbs.uid[n], ph)).transpose().matrix()
                                * neighbs.grad_W_vec[n].matrix())(0, 0);
                            tmpDiff += elFluid.sph_volume(neighbs.uid[n]) * (elFluid.volumeFraction(i, ph) - elFluid.volumeFraction(neighbs.uid[n], ph))
                                * ((elFluid.pos(i) - elFluid.pos(neighbs.uid[n])).transpose().matrix() * neighbs.grad_W_vec[n].matrix())(0, 0)
                                / (pow(neighbs.dis[n], 2) + 10e-6);
                            tmpTran *= transactionCoefficient;
                            tmpDiff *= diffusionCoefficinet;
                            tmpTran = timeStep * (tmpTran + tmpDiff);
                            elFluid.volumeFractionCache(i, ph) += tmpTran;
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for (val_i i = 0; i < elFluid.numPart; i++) {
            for (val_i ph = 0; ph < elPhase.phaseNum; ph++) {
                if (elFluid.volumeFractionCache(i, ph) < 0) {
                    hasNegative = true;
                    elFluid.transact(i) = -1;
                }
            }
        }
    }

    elFluid.volumeFraction_all() = elFluid.volumeFractionCache_all();
}