#!/usr/bin/env python3

"""
# pytochondrion
This module contains the functions required to run mitochondrial simulations
"""

import numpy as np # type: ignore
from params import *

def mitochondrial_simulator(time: float, initial_values: list[float], *args, **kwargs) -> list:
    
    """
    A simulator of mitochondrial energetics
    Returns a list of vectors, each containing the estimated values for each function with the range of t
    Parameters
    ----------
    Y : list(float)
    A vector of initial values
    t : list(float)
    A range of times
    """
    
    # initial values
    nadh_0, uqh_0, c2plus_0, o2_0, h_i_0, atp_ti_0, pi_ti_0, atp_te_0, adp_te_0, pi_te_0, pcr_0 = initial_values
    
    # additional kinetic equations & calculations
    adp_ti = a_iSUM - atp_ti_0
    atp_fe = atp_te_0 / (1 + mg_fe / k_DTe)
    nad = n_t - nadh_0
    v_DH = k_DH * (1 / ((1 + k_mN) / ((nad/nadh_0)**p_D)) ) # (*Substrate dehydrogenation*)
    uq = ubq_t - uqh_0
    c3plus = c_t - c2plus_0
    e_mN = e_mN0 + ((Z/2) * np.log10(nad / nadh_0))
    e_mU = e_mU0 + ((Z/2) * np.log10(uq / uqh_0))
    e_mc = e_mc0 + Z * np.log10(c3plus / c2plus_0)
    pH_i = -np.log10(h_i_0 / 1000000)
    deltapH = Z*(pH_i - pH_e)
    deltap = (1 / (1 - u)) * deltapH
    e_ma = e_mc + deltap*(2 + (2 * u)) / 2
    a_3over2 = 10**((e_ma - e_ma0) / Z) # (*a3plus/a2plus ratio*)
    a2plus = a_t / (1 + a_3over2)
    a3plus = a_t - a2plus
    delta_psi = -(deltap - deltapH)
    deltaG_C1 = e_mU - e_mN - deltap * 4/2
    deltaG_C3 = e_mc - e_mU - deltap * (4 - 2*u) / 2
    deltaG_P = deltaG_P0 / F + Z * np.log10(1000000*atp_ti_0 / (adp_ti * pi_ti_0))
    deltaG_SN = n_A * deltap - deltaG_P
    v_C1 = k_C1 * deltaG_C1 # (*Complex I*)
    v_C3 = k_C3 * deltaG_C3 #
    v_C4 = k_C4 * a2plus * c2plus_0 * 1/(1 + k_mO/o2_0)
    gamma = 10**(deltaG_SN/Z)
    v_SN = k_SN * (gamma - 1)/(gamma + 1) # (*ATP Synthase*)
    v_LK = k_LK1 * (np.exp(k_LK2 * deltap) - 1) # (*Proton leak*)
    atp_me = atp_te_0 - atp_fe
    adp_fe = adp_te_0 /(1 + mg_fe / k_DDe)
    adp_me = adp_te_0 - adp_fe
    atp_fi = atp_ti_0 / (1 + mg_fi / k_DTi)
    atp_mi = atp_ti_0 - atp_fi
    amp_e = a_eSUM - atp_te_0 - adp_te_0
    cr = c_SUM - pcr_0
    psi_i = 0.65 * delta_psi
    psi_e = -0.35 * delta_psi
    pi_je = pi_te_0 / (1 + 10**(pH_e - pK_a))
    pi_ji = pi_ti_0/(1 + 10**(pH_i - pK_a))
    adp_fi = adp_ti / (1 + mg_fi / k_DDi)
    adp_mi = adp_ti - adp_fi
    v_AK = (k_fAK * adp_fe * adp_me) - (k_bAK * atp_me * amp_e) # (*adenylate kinase*)
    v_CK = (k_fCK * adp_te_0 * pcr_0 * h_e) - (k_bCK * atp_te_0 * cr) # (*creatine kinase*)
    v_EX = k_EX * (adp_fe / (adp_fe + atp_fe * 10**(-psi_e / Z)) - adp_fi / (adp_fi + atp_fi * 10**(-psi_i / Z))) * (1/(1 + k_mADP/adp_fe)) # (*ANT*)
    v_PI = k_PI * (pi_je * h_e - pi_ji * h_i_0) # (*phosphate transporter*)
    v_UT = k_UT * 1/(1 + k_mA / atp_te_0) # (*ATP usage*)
    atp_tot = (atp_ti_0 / r_cm) + atp_te_0 # (*calculates total cell [ATP]*)
    pi_tot = (pi_ti_0 / r_cm) + pi_te_0 # (*calculates total cell [Pi]*)
    adp_tot = (adp_ti / r_cm) + adp_te_0 # (*calculates total cell [ADP]*)
    c_0i = (10**(-pH_i) - 10**((-pH_i) - dpH)) / dpH # (*natural buffer capacity in matrix*)
    c_0e = (10**(-pH_e) - 10**((-pH_e) - dpH)) / dpH # (*natural buffer capacity in cytosol*)
    r_buffi = c_buffi / c_0i # (*buffering capacity coefficient for matrix*)
    
    # differential equations
    f0 = (v_DH - v_C1) * (r_cm / b_N) # nadh_of_t
    f1 = (v_C1 - v_C3) * r_cm # uqh_of_t
    f2 = (v_C3 - 2 * v_C4) * 2 * r_cm # c2plus_of_t
    f3 = 0 # -v_C4 # (*or 0*) # o2_of_t
    f4 = -((2 * (2 + 2 * u) * v_C4) + ((4 - 2*u) * v_C3) + (4 * v_C1) - (n_A * v_SN) - (u * v_EX) - ((1 - u) * v_PI) - (v_LK)) * r_cm / r_buffi # h_i_of_t 
    f5 = (v_SN - v_EX) * r_cm # atp_ti_of_t
    f6 = (v_PI - v_SN) * r_cm # pi_ti_of_t
    f7 = (v_EX - v_UT + v_AK + v_CK) # atp_te_of_t
    f8 = (v_UT - v_EX - (2 * v_AK) - v_CK) # adp_te_of_t
    f9 = (v_UT - v_PI) # pi_te_of_t
    f10 = -v_CK # pcr_of_t
    
    return [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]

def calculate_values(ode_solution: list[float]) -> dict[str, float]:
    """
    """
    # unpack values
    nadh_0, uqh_0, c2plus_0, o2_0, h_i_0, atp_ti_0, pi_ti_0, atp_te_0, adp_te_0, pi_te_0, pcr_0 = ode_solution
    
    # additional kinetic equations & calculations
    adp_ti = a_iSUM - atp_ti_0
    atp_fe = atp_te_0 / (1 + mg_fe / k_DTe)
    nad = n_t - nadh_0
    v_DH = k_DH * (1 / ((1 + k_mN) / ((nad/nadh_0)**p_D)) ) # (*Substrate dehydrogenation*)
    uq = ubq_t - uqh_0
    c3plus = c_t - c2plus_0
    e_mN = e_mN0 + ((Z/2) * np.log10(nad / nadh_0))
    e_mU = e_mU0 + ((Z/2) * np.log10(uq / uqh_0))
    e_mc = e_mc0 + Z * np.log10(c3plus / c2plus_0)
    pH_i = -np.log10(h_i_0 / 1000000)
    deltapH = Z*(pH_i - pH_e)
    deltap = (1 / (1 - u)) * deltapH
    e_ma = e_mc + deltap*(2 + (2 * u)) / 2
    a_3over2 = 10**((e_ma - e_ma0) / Z) # (*a3plus/a2plus ratio*)
    a2plus = a_t / (1 + a_3over2)
    a3plus = a_t - a2plus
    delta_psi = -(deltap - deltapH)
    deltaG_C1 = e_mU - e_mN - deltap * 4/2
    deltaG_C3 = e_mc - e_mU - deltap * (4 - 2*u) / 2
    deltaG_P = deltaG_P0 / F + Z * np.log10(1000000*atp_ti_0 / (adp_ti * pi_ti_0))
    deltaG_SN = n_A * deltap - deltaG_P
    v_C1 = k_C1 * deltaG_C1 # (*Complex I*)
    v_C3 = k_C3 * deltaG_C3 #
    v_C4 = k_C4 * a2plus * c2plus_0 * 1/(1 + k_mO/o2_0)
    gamma = 10**(deltaG_SN/Z)
    v_SN = k_SN * (gamma - 1)/(gamma + 1) # (*ATP Synthase*)
    v_LK = k_LK1 * (np.exp(k_LK2 * deltap) - 1) # (*Proton leak*)
    atp_me = atp_te_0 - atp_fe
    adp_fe = adp_te_0 /(1 + mg_fe / k_DDe)
    adp_me = adp_te_0 - adp_fe
    atp_fi = atp_ti_0 / (1 + mg_fi / k_DTi)
    atp_mi = atp_ti_0 - atp_fi
    amp_e = a_eSUM - atp_te_0 - adp_te_0
    cr = c_SUM - pcr_0
    psi_i = 0.65 * delta_psi
    psi_e = -0.35 * delta_psi
    pi_je = pi_te_0 / (1 + 10**(pH_e - pK_a))
    pi_ji = pi_ti_0/(1 + 10**(pH_i - pK_a))
    adp_fi = adp_ti / (1 + mg_fi / k_DDi)
    adp_mi = adp_ti - adp_fi
    v_AK = (k_fAK * adp_fe * adp_me) - (k_bAK * atp_me * amp_e) # (*adenylate kinase*)
    v_CK = (k_fCK * adp_te_0 * pcr_0 * h_e) - (k_bCK * atp_te_0 * cr) # (*creatine kinase*)
    v_EX = k_EX * (adp_fe / (adp_fe + atp_fe * 10**(-psi_e / Z)) - adp_fi / (adp_fi + atp_fi * 10**(-psi_i / Z))) * (1/(1 + k_mADP/adp_fe)) # (*ANT*)
    v_PI = k_PI * (pi_je * h_e - pi_ji * h_i_0) # (*phosphate transporter*)
    v_UT = k_UT * 1/(1 + k_mA / atp_te_0) # (*ATP usage*)
    atp_tot = (atp_ti_0 / r_cm) + atp_te_0 # (*calculates total cell [ATP]*)
    pi_tot = (pi_ti_0 / r_cm) + pi_te_0 # (*calculates total cell [Pi]*)
    adp_tot = (adp_ti / r_cm) + adp_te_0 # (*calculates total cell [ADP]*)
    c_0i = (10**(-pH_i) - 10**((-pH_i) - dpH)) / dpH # (*natural buffer capacity in matrix*)
    c_0e = (10**(-pH_e) - 10**((-pH_e) - dpH)) / dpH # (*natural buffer capacity in cytosol*)
    r_buffi = c_buffi / c_0i # (*buffering capacity coefficient for matrix*)
    return {"adp_ti" : adp_ti, 
            "atp_fe" : atp_fe}