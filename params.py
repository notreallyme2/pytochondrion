#!/usr/bin/env python3

"""
# params
This module contains all the paramaters and constants used by the model
"""

# ****CONSTANTS****
T = 298 #(*Cell temp. in Kelvin*)
R = 0.0083 #(*kJ/M/K*)
F = 0.0965 #(*kJ/M/mV*)
S = 2.303*R*T
Z = 2.303*R*T/F
u = 0.861 #(* u=\[CapitalDelta]\[Psi]/\[CapitalDelta]p *)

# ****PARAMETERS****
r_cm = 3.35 # (*cell volume / mt volume ratio*)
b_N = 5 # (*buffering capacity coefficient for NAD*)
c_t = 270 # (*total concentration of cytochrome c in uM*)
ubq_t = 1350 # (*total concentration of ubiquinone in uM*)
n_t = 2970 # (*total concentration of NAD in uM*)

p_D = 0.8 # (*Dehydrogenation sensitivity to NAD ratio*)

k_DH = 96293 # (*Substrate dehydrogenation parameters in uM/min*)
k_mN = 100
k_C1 = 819.61 # (*in uM/mV/min*)
k_C3 = 467.90 # (*in uM/mV/min*)
k_C4 = 12.348 # (*in uM/min*)
k_mO = 120 # (*uM_apparent k_mO]=0.8 uM*)
k_SN = 117706 # (*in uM/min*)
k_EX = 187185 # (*in uM/min*)
k_mADP = 3.5 # (*uM*)
k_PI = 238.11 # (*in /uM/min*)
k_UT = 12244 # (*ATP usage rate in uM/min# 12244(low work)-61220(high work)*)
k_mA = 150 # (*uM*)
k_LK1 = 8.5758 # (*normally 8.5758 in uM/min*)
k_LK2 = 0.038
k_fAK = 862.10 # (*adenylate kinase forward rate constant in uM/min*)
k_bAK = 22.747 # (*adenylate kinase reverse rate constant in uM/min*)
k_fCK = 1.9258 # (*creatine kinase forward rate constant in uM/min# times \ 0.656 for CHF (Weiss et al. (2005)*)
k_bCK = 0.00087538 # (*creatine kinase reverse rate constant in uM/min*)

a_eSUM = 6700.2 # (*uM_total cytosolic adenine nucleotide conc.*)
a_iSUM = 16260 # (*uM_total mitochondrial adenine nucleotide conc.*)
c_SUM = 25000 # (*total creatine concentration_Normal: 25000 (orig. model)_\ CHF: 16100 (Weiss et al. (2005)*)
p_SUM = 45582#(*total phosphate pool_not actually used!!*)
mg_fe = 4000 # (*uM_free cytosolic magnesium conc.*)
mg_fi = 380 # (*uM_free mt. magnesium conc.*)
k_DTe = 24 # (*uM_Mg diss. const. for cyt. ATP*)
k_DDe = 347 # (*uM_Mg diss. const. for cyt. ADP*)
k_DTi = 17 # (*uM_Mg diss. const. for mt. ATP*)
k_DDi = 282 # (*uM_Mg diss. const. for mt. ADP*)
pH_e = 7.0
h_e = 1000000 * 10**(-pH_e)
dpH = 0.001
c_buffi = 0.022 # (* moles of protons/pH unit_buffer capacity of matrix*)
pK_a = 6.8 # (*phosphate proton-buffering midpoint*)
n_A = 2.5 # (*H+/ATP stoichiometry of the synthase*)
deltaG_P0 = 31.9 # (*kJ/mol*)
e_mN0 = -320 # (*standard redox potential of NAD couple_in mV*)
e_mU0 = 85 # (*standard redox potential of UQ couple_in mV*)
e_mc0 = 250 # (*standard redox potential of cyt. c couple_in mV*)
a_t = 135 # (*total conc. of cyt. a_in uM*)
e_ma0 = 540 # (*standard redox potential of cyt. a couple_in mV*)