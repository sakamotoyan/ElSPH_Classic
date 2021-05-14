#pragma once

#include "elDataStructureTemplate.h"
/****************** OBJECT GRID ******************/
class ObjectGrid : public DataGrid
{
private:
    /*** Addressing Hint***/
    val_i universalId_1i[3] = {0, 0, 1}; // [object type] e.g. for all fluid oType is FLUID (defined in util.h)
    val_i neighbNum_1i[3] = {0, 1, 1};   // [object Id] e.g. the first fluid oId is 1 and also to the first boundary object

    val_i position_3f[3] = {0, 0, 3}; // postion is in posIn[0] value block, start from posIn[1] row, takes up posIn[2]*sizeof(val) mem
    val_i velocity_3f[3] = {1, 0, 3}; // statement of velocity, same as postion
public:
    inline ObjectGrid(val_i num = 0, const val_i *config = objectGridConfig) : DataGrid(num, config) {}

    /*** Addressing ***/
    inline FloatBlockRef pos(val_i i)
    {
        val_i(&ref)[3] = position_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline Array3_f pos_const(val_i i)
    {
        val_i(&ref)[3] = position_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef pos_all()
    {
        val_i(&ref)[3] = position_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef vel(val_i i)
    {
        val_i(&ref)[3] = velocity_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef vel_all()
    {
        val_i(&ref)[3] = velocity_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_i &uid(val_i i)
    {
        val_i(&ref)[3] = universalId_1i; // change this por different attributes
        return intBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }

    inline val_i &neighbNum(val_i i)
    {
        val_i(&ref)[3] = neighbNum_1i; // change this por different attributes
        return intBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline IntegerBlockRef neighbNum_all()
    {
        val_i(&ref)[3] = neighbNum_1i; // change this por different attributes
        return intBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f rest_compressionRate(val_i i)
    {
        return 1;
    }
};

/****************** ELFLUID ******************/
class Elfluid : public ObjectGrid
{
private:
    FloatBlock multiphaseBlock;

    val_i mass_1f[3] = {1, 3, 1}; // X
    val_i sph_volume_1f[3] = {1, 4, 1};
    val_i rest_density_1f[3] = {1, 5, 1}; // Psi0
    val_i sph_Psi_1f[3] = {1, 6, 1};      // Psi
    val_i adv_vel_3f[3] = {1, 7, 3};
    val_i rest_volume_1f[3] = {1, 10, 1}; // X
    val_i alpha_term1_3f[3] = {1, 11, 3};
    val_i alpha_term2_1f[3] = {1, 14, 1};
    val_i alpha_1f[3] = {1, 15, 1};
    val_i adv_acce_3f[3] = { 1,16,3 };
    val_i adv_Psi_1f[3] = { 1, 19, 1 };
    val_i pressure_1f[3] = { 1,20,1 };
    val_i pressure_force_3f[3] = { 1,21,3 };
    val_i rest_compressionRate_1f[3] = { 1,25,1 };
    val_i adv_Psi_changeRate_1f[3] = { 1,26,1 };
    val_i adv_Psi_pj_1f[3] = { 1,27,1 };
    val_i devi_pi_3f[3]={ 1,28,3 };
    val_i devi_pj_3f[3]={ 1,31,3 };
    val_i devi_pi_tmp_3f[3] = { 1,34,3 };
    val_i volumeFraction_5f[3] = { 2,0,5 };
    val_i drift_vel_15f[3] = { 2,5,15 };
    val_i color_4f[3] = { 2,20,4 };
    val_i mt_zeta_1f[3] = { 2,24,1 };
    val_i intermediate_adv_acce_3f[3] = { 2,25,3 };
    val_i volumeFractionCache_5f[3] = { 2,28,5 };
    val_i transaction_1f[3] = { 2,33,1 };
    val_i mt_gamma_1f[3] = { 2,34,1 };
    // rest_compressionRate_1f Psi 0

public:
    inline Elfluid(val_i num = MaxPartNum_fluid, const val_i* config = fluidGridConfig) : ObjectGrid(num, config) {
        val_i phaseBlockElementNumber = 37;
        floatBlocks.push_back(&multiphaseBlock);
        multiphaseBlock.resize(phaseBlockElementNumber, num);
        multiphaseBlock.dataMat.setZero();
    }

    inline val_f &rest_density(val_i i)
    {
        val_i(&ref)[3] = rest_density_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_density_all()
    {
        val_i(&ref)[3] = rest_density_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& rest_compressionRate(val_i i)
    {
        val_i(&ref)[3] = rest_compressionRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_compressionRate_all()
    {
        val_i(&ref)[3] = rest_compressionRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &rest_volume(val_i i)
    {
        val_i(&ref)[3] = rest_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_volume_all()
    {
        val_i(&ref)[3] = rest_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& adv_Psi(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_all()
    {
        val_i(&ref)[3] = adv_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &sph_volume(val_i i)
    {
        val_i(&ref)[3] = sph_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef sph_volume_all()
    {
        val_i(&ref)[3] = sph_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &mass(val_i i)
    {
        val_i(&ref)[3] = mass_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef mass_all()
    {
        val_i(&ref)[3] = mass_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &sph_Psi(val_i i)
    {
        val_i(&ref)[3] = sph_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef sph_Psi_all()
    {
        val_i(&ref)[3] = sph_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef adv_vel(val_i i)
    {
        val_i(&ref)[3] = adv_vel_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef adv_vel_all()
    {
        val_i(&ref)[3] = adv_vel_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef alpha_term1(val_i i)
    {
        val_i(&ref)[3] = alpha_term1_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef alpha_term1_all()
    {
        val_i(&ref)[3] = alpha_term1_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &alpha_term2(val_i i)
    {
        val_i(&ref)[3] = alpha_term2_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef alpha_term2_all()
    {
        val_i(&ref)[3] = alpha_term2_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f &alpha(val_i i)
    {
        val_i(&ref)[3] = alpha_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef alpha_all()
    {
        val_i(&ref)[3] = alpha_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef adv_acce(val_i i)
    {
        val_i(&ref)[3] = adv_acce_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef adv_acce_all()
    {
        val_i(&ref)[3] = adv_acce_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& pressure(val_i i)
    {
        val_i(&ref)[3] = pressure_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef pressure_all()
    {
        val_i(&ref)[3] = pressure_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef pressure_force(val_i i)
    {
        val_i(&ref)[3] = pressure_force_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef pressure_force_all()
    {
        val_i(&ref)[3] = pressure_force_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& adv_Psi_changeRate(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_changeRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_changeRate_all()
    {
        val_i(&ref)[3] = adv_Psi_changeRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& adv_Psi_pj(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_pj_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_pj_all()
    {
        val_i(&ref)[3] = adv_Psi_pj_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef devi_pi(val_i i)
    {
        val_i(&ref)[3] = devi_pi_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef devi_pi_all()
    {
        val_i(&ref)[3] = devi_pi_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef devi_pj(val_i i)
    {
        val_i(&ref)[3] = devi_pj_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef devi_pj_all()
    {
        val_i(&ref)[3] = devi_pj_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef devi_pi_tmp(val_i i)
    {
        val_i(&ref)[3] = devi_pi_tmp_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef devi_pi_tmp_all()
    {
        val_i(&ref)[3] = devi_pi_tmp_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& volumeFraction(val_i i, val_i ph)
    {
        val_i(&ref)[3] = volumeFraction_5f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1] + ph, i, 1, 1)(0, 0);
    }
    inline FloatBlockRef volumeFractions(val_i i)
    {
        val_i(&ref)[3] = volumeFraction_5f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef volumeFraction_all()
    {
        val_i(&ref)[3] = volumeFraction_5f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef drift_vel(val_i i, val_i ph)
    {
        val_i(&ref)[3] = drift_vel_15f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1] + (ph * 3), i, 3, 1);
    }
    inline FloatBlockRef drift_vels(val_i i)
    {
        val_i(&ref)[3] = drift_vel_15f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef drift_vel_all()
    {
        val_i(&ref)[3] = drift_vel_15f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline FloatBlockRef color(val_i i)
    {
        val_i(&ref)[3] = color_4f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline val_f& colorA(val_i i) {return color(i)(3, 0);}
    
    inline FloatBlockRef color_all()
    {
        val_i(&ref)[3] = color_4f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& mt_zeta(val_i i)
    {
        val_i(&ref)[3] = mt_zeta_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef mt_zeta_all()
    {
        val_i(&ref)[3] = mt_zeta_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef intermediate_adv_acce(val_i i)
    {
        val_i(&ref)[3] = intermediate_adv_acce_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef intermediate_adv_acce_all()
    {
        val_i(&ref)[3] = intermediate_adv_acce_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& volumeFractionCache(val_i i, val_i ph)
    {
        val_i(&ref)[3] = volumeFractionCache_5f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1] + ph, i, 1, 1)(0, 0);
    }
    inline FloatBlockRef volumeFractionCache_all()
    {
        val_i(&ref)[3] = volumeFractionCache_5f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& transact(val_i i)
    {
        val_i(&ref)[3] = transaction_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef transact_all()
    {
        val_i(&ref)[3] = transaction_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& mt_gamma(val_i i)
    {
        val_i(&ref)[3] = mt_gamma_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef mt_gamma_all()
    {
        val_i(&ref)[3] = mt_gamma_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
};

/****************** ELBOUND ******************/
class Elbound : public ObjectGrid
{

private:
    val_i contributionWeight_1f[3] = {1, 3, 1}; // [3] of floatBlock[1]
    val_i rest_volume_1f[3] = {1, 4, 1};        // [4] of floatBlock[1]
    val_i rest_density_1f[3] = {1, 5, 1};
    val_i mass_1f[3] = {1, 6, 1};
    val_i sph_Psi_1f[3] = {1, 7, 1};
    val_i alpha_term2_1f[3] = {1, 8, 1};
    val_i alpha_1f[3] = {1, 9, 1};
    val_i sph_volume_1f[3] = { 1, 10, 1 };
    val_i adv_Psi_1f[3] = { 1, 11, 1 };
    val_i pressure_1f[3] = { 1,12,1 };
    val_i pressure_force_3f[3] = { 1,13,3 };
    val_i rest_compressionRate_1f[3] = { 1,16,1 };
    val_i adv_Psi_changeRate_1f[3] = { 1,17,1 };
    val_i adv_Psi_pj_1f[3] = { 1,18,1 };

public:
	inline Elbound(val_i num = MaxPartNum_bound, const val_i* config = boundGridConfig) : ObjectGrid(num, config) {}

    inline val_f &weight(val_i i)
    {
        val_i(&ref)[3] = contributionWeight_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef weight_all()
    {
        val_i(&ref)[3] = contributionWeight_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &rest_volume(val_i i)
    {
        val_i(&ref)[3] = rest_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_volume_all()
    {
        val_i(&ref)[3] = rest_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &rest_density(val_i i)
    {
        val_i(&ref)[3] = rest_density_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_density_all()
    {
        val_i(&ref)[3] = rest_density_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& rest_compressionRate(val_i i)
    {
        val_i(&ref)[3] = rest_compressionRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef rest_compressionRate_all()
    {
        val_i(&ref)[3] = rest_compressionRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &mass(val_i i)
    {
        val_i(&ref)[3] = mass_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef mass_all()
    {
        val_i(&ref)[3] = mass_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &sph_Psi(val_i i)
    {
        val_i(&ref)[3] = sph_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef sph_Psi_all()
    {
        val_i(&ref)[3] = sph_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& alpha_term2(val_i i)
    {
        val_i(&ref)[3] = alpha_term2_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef alpha_term2_all()
    {
        val_i(&ref)[3] = alpha_term2_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f &alpha(val_i i)
    {
        val_i(&ref)[3] = alpha_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef alpha_all()
    {
        val_i(&ref)[3] = alpha_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& sph_volume(val_i i)
    {
        val_i(&ref)[3] = sph_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef sph_volume_all()
    {
        val_i(&ref)[3] = sph_volume_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }

    inline val_f& adv_Psi(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_all()
    {
        val_i(&ref)[3] = adv_Psi_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& pressure(val_i i)
    {
        val_i(&ref)[3] = pressure_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef pressure_all()
    {
        val_i(&ref)[3] = pressure_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline FloatBlockRef pressure_force(val_i i)
    {
        val_i(&ref)[3] = pressure_force_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline FloatBlockRef pressure_force_all()
    {
        val_i(&ref)[3] = pressure_force_3f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& adv_Psi_changeRate(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_changeRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_changeRate_all()
    {
        val_i(&ref)[3] = adv_Psi_changeRate_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
    inline val_f& adv_Psi_pj(val_i i)
    {
        val_i(&ref)[3] = adv_Psi_pj_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline FloatBlockRef adv_Psi_pj_all()
    {
        val_i(&ref)[3] = adv_Psi_pj_1f;
        return floatBlocks[ref[0]]->dataMat.block(ref[1], 0, ref[2], numPart);
    }
};

/****************** NeighbAux ******************/
//Aux dense data grid for Neighb module
class NeighbAux : public DataGrid
{
public:
    inline NeighbAux(val_i num = 0, const val_i *config = neighbAuxConfig) : DataGrid(num, config) {}

    val_i indexLocation_i[3] = {0, 0, 3};
    val_i codedIndex_i[3] = {0, 3, 1};
    val_i particleId_i[3] = {0, 4, 1};

    inline IntegerBlockRef bgIndex_vec(val_i i)
    {
        val_i(&ref)[3] = indexLocation_i;
        return intBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1);
    }
    inline val_i &bgIndex_int(val_i i)
    {
        val_i(&ref)[3] = codedIndex_i;
        return intBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
    inline val_i &pid(val_i i)
    {
        val_i(&ref)[3] = particleId_i;
        return intBlocks[ref[0]]->dataMat.block(ref[1], i, ref[2], 1)(0, 0);
    }
};

/****************** NEIGHBOUR RELATION ******************/
class NeighbRelation
{
public:
    val_i neighbNum = 0;
    vector<val_i> uid;
    vector<val_f> dis;
    vector<val_f> W;
    vector<val_f> grad_W;
    vector<Array<val_f,3,1>> grad_W_vec;

    inline void clear()
    {
        this->neighbNum = 0;
        this->uid.clear();
        this->dis.clear();
        this->W.clear();
        this->grad_W.clear();
        this->grad_W_vec.clear();
    }

    inline void push(val_i uid, val_f dis, val_f W, val_f grad_W, Array<val_f, 3, 1>& dis_vec)
    {
        this->neighbNum++;
        this->uid.push_back(uid);
        this->dis.push_back(dis);
        this->W.push_back(W);
        this->grad_W.push_back(grad_W); 
        if (dis == 0) {
            this->grad_W_vec.push_back(dis_vec);
        }
        else {
            dis_vec *= (grad_W / dis);
            this->grad_W_vec.push_back(dis_vec);
        }
    }
};

/****************** NEIGHBOUR MODULE ******************/
class Elneighb
{
    Elvari elvari;
private:
    inline Array3_i locateBgIndex(FloatBlockRef pos)
    {
        return ((pos.array() - elvari.spaceMin) / neighbSearchGridSize).floor().cast<val_i>();
    }
    inline const val_i codeBgIndex(IntegerBlockRef bgIndex)
    {
        return bgIndex(0, 0) + (gridNum(0) * bgIndex(1, 0)) + (gridNum(0) * gridNum(1) * bgIndex(2, 0));
    }
    inline const val_i codeBgIndex(Array3_i bgIndex)
    {
        return bgIndex(0, 0) + (gridNum(0) * bgIndex(1, 0)) + (gridNum(0) * gridNum(1) * bgIndex(2, 0));
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    Array<val_i, Dynamic, Dynamic> backGroundGrid;
    vector<NeighbRelation> neighbRelations;

    Array3_i gridNum;
    val_i gridTotal;
    NeighbAux aux;
    omp_lock_t *lock;
    int ii;
    val_i searchRange = 2; // +- 1block

    inline Elneighb()
    {
        aux.reset(MaxPartNum);
        neighbRelations.resize(MaxPartNum);
        gridNum = ((elvari.spaceMax - elvari.spaceMin) / neighbSearchGridSize).ceil().cast<val_i>();
        gridTotal = iterProduct(gridNum);
        backGroundGrid.resize(gridMembers, gridTotal);
        backGroundGrid.setZero();
        lock = new omp_lock_t[gridTotal];
        for (ii = 0; ii < gridTotal; ii++) {
            omp_init_lock(&(lock[ii]));
        }
        
    }

    inline void loopNonZeros()
    {
        val_i sum = 0;
        for (val_i col = 0; col < backGroundGrid.cols(); col++)
        {
            val_i &counter = backGroundGrid(nonZeros, col);
            sum += counter;
            if (counter > 0)
                cout << "has " << counter << " non zeros" << endl;
            for (val_i row = 0; row < counter; row++)
            {
                cout << "\t value: " << backGroundGrid(row, col) << "\t row: " << row << "\t\t col: " << col << endl;
            }
        }
        cout << "summary: " << sum << " non zeros in all" << endl;
    }

    inline void push(val_i col, val_i value)
    {
        omp_set_lock(&(lock[col]));
        //cout << "value being pushed to: " << col << endl;
        val_i &counter = backGroundGrid(nonZeros, col);
        backGroundGrid(counter, col) = value;
        counter++;
        omp_unset_lock(&(lock[col]));
    }

    inline void positioning(Elfluid &elFluid, Elbound &elBound)
    {
        // 1. [elFluid] [elBound] destribution of uid
        // 2. [backGroundGrid] logging thoes uid
        // 3. [aux] logging cel_number(vec)(coded) for all elFluid and elBound
        backGroundGrid.row(nonZeros).setZero();
#pragma omp parallel for // for each col in elFluid
        for (val_i i = 0; i < elFluid.numPart; i++)
        {
            elFluid.uid(i) = i;
            val_i &uid = elFluid.uid(i);
            aux.bgIndex_vec(uid) = locateBgIndex(elFluid.pos(i));
            aux.bgIndex_int(uid) = codeBgIndex(aux.bgIndex_vec(uid));
            //#pragma omp critical
            this->push(aux.bgIndex_int(uid), uid);
            // cout << "listen push from fluid" << aux.bgIndex_vec(i).transpose() << "coded as :" << aux.bgIndex_int(i) << endl;
        }

#pragma omp parallel for // for each col in elBound
        for (val_i i = 0; i < elBound.numPart; i++)
        {
            elBound.uid(i) = i + MaxPartNum_fluid;
            val_i &uid = elBound.uid(i);
            aux.bgIndex_vec(uid) = locateBgIndex(elBound.pos(i));
            aux.bgIndex_int(uid) = codeBgIndex(aux.bgIndex_vec(uid));
            //#pragma omp critical
            this->push(aux.bgIndex_int(uid), uid);
            // cout << "listen push from bound" << aux.bgIndex_vec(elBound.uid(i)).transpose() << "coded as :" << aux.bgIndex_int(elBound.uid(i)) << endl;
        }
    }

    inline void updateRelations(Elfluid &elFluid, Elbound &elBound)
    {
        val_i upperLimit = backGroundGrid.cols();
#pragma omp parallel for // for each col in elFluid
        for (val_i i = 0; i < elFluid.numPart; i++)
        {
            val_i &uid = elFluid.uid(i);
            neighbRelations[uid].clear();
            for (val_i x = scanRange[0]; x < scanRange[1]; x++)
            {
                for (val_i y = scanRange[0]; y < scanRange[1]; y++)
                {
                    for (val_i z = scanRange[0]; z < scanRange[1]; z++)
                    {
                        val_i codedBgIndex = codeBgIndex(aux.bgIndex_vec(uid) + Array3_i(x, y, z));
                        if (codedBgIndex > -1 && codedBgIndex < upperLimit)
                        {
                            for (val_i nz = 0; nz < backGroundGrid(nonZeros, codedBgIndex); nz++)
                            {
                                Array<val_f, 3, 1> dis_vec;
                                val_f dis;
                                val_i &neighb_uid = backGroundGrid(nz, codedBgIndex);
                                if (isFluid(neighb_uid))
                                {
                                    dis_vec = elFluid.pos(i) - elFluid.pos(neighb_uid);
                                    dis = dis_vec.matrix().norm();
                                    if (dis < smoothRadius)
                                    {
                                        neighbRelations[uid].push(neighb_uid, dis, W(dis), grad_W(dis), dis_vec);
                                    }
                                }
                                else
                                {
                                    dis_vec = elFluid.pos(i) - elBound.pos(bid(neighb_uid));
                                    dis = dis_vec.matrix().norm();
                                    if (dis < smoothRadius)
                                    {
                                        neighbRelations[uid].push(neighb_uid, dis, W(dis), grad_W(dis), dis_vec);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            elFluid.neighbNum(i) = neighbRelations[uid].neighbNum;
        }
#pragma omp parallel for // for each col in elBound
        for (val_i i = 0; i < elBound.numPart; i++)
        {
            val_i &uid = elBound.uid(i);
            neighbRelations[uid].clear();
            for (val_i x = scanRange[0]; x < scanRange[1]; x++)
            {
                for (val_i y = scanRange[0]; y < scanRange[1]; y++)
                {
                    for (val_i z = scanRange[0]; z < scanRange[1]; z++)
                    {
                        val_i codedBgIndex = codeBgIndex(aux.bgIndex_vec(uid) + Array3_i(x, y, z));
                        if (codedBgIndex > -1 && codedBgIndex < upperLimit)
                        {
                            for (val_i nz = 0; nz < backGroundGrid(nonZeros, codedBgIndex); nz++)
                            {
                                Array<val_f, 3, 1> dis_vec;
                                val_f dis;
                                val_i &neighb_uid = backGroundGrid(nz, codedBgIndex);
                                if (isFluid(neighb_uid))
                                {
                                    dis_vec = elBound.pos(i) - elFluid.pos(neighb_uid);
                                    dis = dis_vec.matrix().norm();
                                    if (dis < smoothRadius)
                                    {
                                        neighbRelations[uid].push(neighb_uid, dis, W(dis), grad_W(dis), dis_vec);
                                    }
                                }
                                else
                                {
                                    dis_vec = elBound.pos(i) - elBound.pos(bid(neighb_uid));
                                    val_f dis = dis_vec.matrix().norm();
                                    if (dis < smoothRadius)
                                    {
                                        neighbRelations[uid].push(neighb_uid, dis, W(dis), grad_W(dis), dis_vec);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            elBound.neighbNum(i) = neighbRelations[uid].neighbNum;
        }
    }
};