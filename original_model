(****CONSTANTS****)
T = 298; (*Cell temp. in Kelvin*)
R = 0.0083;(*kJ/M/K*)
F = 0.0965;(*kJ/M/mV*)
S = 2.303*R*T;
Z = 2.303*R*T/F;
u = 0.861;(* u=\[CapitalDelta]\[Psi]/\[CapitalDelta]p *)

(****PARAMETERS****)
Subscript[r, cm] = 3.35; (*cell volume / mt volume ratio*)
Subscript[b, N] = 5; (*buffering capacity coefficient for NAD*)
Subscript[c, t] = 270; (*total concentration of cytochrome c in uM*)
Subscript[ubq, t] = 1350; (*total concentration of ubiquinone in uM*)
Subscript[n, t] = 2970; (*total concentration of NAD in uM*)
Subscript[k, DH] = 
  96293;(*Substrate dehydrogenation parameters in uM/min*)
Subscript[k, mN] = 100;
Subscript[p, D] = 0.8;(*Dehydrogenation sensitivity to NAD ratio*)
Subscript[k, C1] = 819.61;(*in uM/mV/min*)
Subscript[k, C3] = 467.90;(*in uM/mV/min*)
Subscript[k, C4] = 12.348;(*in uM/min*)
Subscript[k, mO] = 120; (*uM, apparent Subscript[k, mO]=0.8 uM*)
Subscript[k, SN] = 117706;(*in uM/min*)
Subscript[k, EX] = 187185;(*in uM/min*)
Subscript[k, mADP] = 3.5; (*uM*)
Subscript[k, PI] = 238.11;(*in /uM/min*)
Subscript[k, UT] := 
  12244;(*ATP usage rate in uM/min; 12244(low work)-61220(high work)*)
Subscript[k, mA] = 150;(*uM*)
Subscript[k, LK1] = 8.5758;(*normally 8.5758 in uM/min*)
Subscript[k, LK2] = 0.038;
Subscript[k, fAK] = 
  862.10; (*adenylate kinase forward rate constant in uM/min*)
Subscript[k, bAK] = 
  22.747;(*adenylate kinase reverse rate constant in uM/min*)
Subscript[k, fCK] = 
  1.9258;(*creatine kinase forward rate constant in uM/min; times \
0.656 for CHF (Weiss et al. (2005)*)
Subscript[k, bCK] = 
  0.00087538;(*creatine kinase reverse rate constant in uM/min*)
Subscript[a, eSUM] = 
  6700.2; (*uM, total cytosolic adenine nucleotide conc.*)
Subscript[a, iSUM] = 
  16260;(*uM, total mitochondrial adenine nucleotide conc.*)
Subscript[c, SUM] = 
  25000; (*total creatine concentration, Normal: 25000 (orig. model), \
CHF: 16100 (Weiss et al. (2005)*)
Subscript[p, SUM] = 45582;(*total phosphate pool, not actually used!!*)
Subscript[mg, fe] = 4000;(*uM, free cytosolic magnesium conc.*)
Subscript[mg, fi] = 380;(*uM, free mt. magnesium conc.*)
Subscript[k, DTe] = 24;(*uM, Mg diss. const. for cyt. ATP*)
Subscript[k, DDe] = 347;(*uM, Mg diss. const. for cyt. ADP*)
Subscript[k, DTi] = 17;(*uM, Mg diss. const. for mt. ATP*)
Subscript[k, DDi] = 282;(*uM, Mg diss. const. for mt. ADP*)
Subscript[pH, e] = 7.0;
Subscript[h, e] = 1000000*10^-Subscript[pH, e];
dpH = 0.001;
Subscript[c, buffi] = 
  0.022; (* moles of protons/pH unit, buffer capacity of matrix*)
Subscript[pK, a] = 6.8;(*phosphate proton-buffering midpoint*)
Subscript[n, A] = 2.5;(*H+/ATP stoichiometry of the synthase*)
Subscript[deltaG, P0] = 31.9;(*kJ/mol*)
Subscript[e, 
  mN0] = -320;(*standard redox potential of NAD couple, in mV*)
Subscript[e, mU0] = 85;(*standard redox potential of UQ couple, in mV*)
Subscript[e, mc0] = 
  250;(*standard redox potential of cyt. c couple, in mV*)
Subscript[a, t] = 135;(*total conc. of cyt. a, in uM*)
Subscript[e, ma0] = 
  540;(*standard redox potential of cyt. a couple, in mV*)

  (****KINETIC EQUATIONS AND CALCULATIONS****)
Subscript[v, DH] := 
  Subscript[k, DH]
    1/(1 + Subscript[k, mN]/(nad/nadh[t]))^Subscript[p, 
    D]; (*Substrate dehydrogenation*)
Subscript[v, C1] := 
  Subscript[k, C1]*Subscript[deltaG, C1]; (*Complex I*)
Subscript[v, C3] := Subscript[k, C3]*Subscript[deltaG, C3];
Subscript[v, C4] := 
  Subscript[k, C4]*a2plus*c2plus[t]*1/(1 + Subscript[k, mO]/o2[t]);
Subscript[v, SN] := 
  Subscript[k, SN] (\[Gamma] - 1)/(\[Gamma] + 1);(*ATP Synthase*)
\[Gamma] := 10^(Subscript[deltaG, SN]/Z);
Subscript[v, LK] := 
  Subscript[k, LK1]*(Exp[Subscript[k, LK2]*deltap] - 1);(*Proton leak*)
Subscript[v, 
  AK] := (Subscript[k, fAK]*Subscript[adp, fe]*Subscript[adp, 
     me]) - (Subscript[k, bAK]*Subscript[atp, me]*Subscript[amp, 
     e]);(*adenylate kinase*)
Subscript[v, 
  CK] := (Subscript[k, fCK]*Subscript[adp, te][t]*pcr[t]*Subscript[h, 
     e]) - (Subscript[k, bCK]*Subscript[atp, te][t]*
     cr);(*creatine kinase*)
Subscript[v, EX] := 
  Subscript[k, 
   EX]*(Subscript[adp, fe]/(
     Subscript[adp, fe] + 
      Subscript[atp, fe]*10^(-Subscript[\[Psi], e]/Z)) - Subscript[
     adp, fi]/(
     Subscript[adp, fi] + 
      Subscript[atp, fi]*10^(-Subscript[\[Psi], i]/Z))) (1/(
    1 + Subscript[k, mADP]/Subscript[adp, fe]));(*ANT*)
Subscript[v, PI] := 
  Subscript[k, 
   PI] (Subscript[pi, je]*Subscript[h, e] - 
     Subscript[pi, ji]*Subscript[h, i][t]);(*phosphate transporter*)
Subscript[v, UT] := 
  Subscript[k, UT]*1/(
   1 + Subscript[k, mA]/Subscript[atp, te][t]);(*ATP usage*)
Subscript[deltaG, C1] := 
  Subscript[e, mU] - Subscript[e, mN] - deltap*4/2;
Subscript[deltaG, C3] := 
  Subscript[e, mc] - Subscript[e, mU] - deltap*(4 - 2 u)/2;
Subscript[deltaG, SN] := Subscript[n, A]*deltap - Subscript[deltaG, P];
Subscript[deltaG, P] := 
  Subscript[deltaG, P0]/F + 
   Z*Log[10, 
     1000000*Subscript[atp, ti][
        t]/(Subscript[adp, ti]*Subscript[pi, ti][t])];
Subscript[adp, ti] := Subscript[a, iSUM] - Subscript[atp, ti][t];
Subscript[atp, fe] := 
  Subscript[atp, te][
    t]/(1 + Subscript[mg, fe]/Subscript[k, 
      DTe]);(*the following equations describe free and \
magnesium-bound adenine nucleotides*)
Subscript[atp, me] := Subscript[atp, te][t] - Subscript[atp, fe];
Subscript[adp, fe] := 
  Subscript[adp, te][t]/(1 + Subscript[mg, fe]/Subscript[k, DDe]);
Subscript[adp, me] := Subscript[adp, te][t] - Subscript[adp, fe];
Subscript[atp, fi] := 
  Subscript[atp, ti][t]/(1 + Subscript[mg, fi]/Subscript[k, DTi]);
Subscript[atp, mi] := Subscript[atp, ti][t] - Subscript[atp, fi];
Subscript[atp, 
  tot] := (Subscript[atp, ti][t]/Subscript[r, cm]) +Subscript[atp, 
    te][t];(*calculates total cell [ATP]*)
Subscript[pi, 
  tot] := (Subscript[pi, ti][t]/Subscript[r, cm]) + 
   Subscript[pi, te][t];(*calculates total cell [Pi]*)
Subscript[adp, 
  tot] := (Subscript[adp, ti]/Subscript[r, cm]) + 
   Subscript[adp, te][t];(*calculates total cell [ADP]*)
Subscript[adp, fi] := 
  Subscript[adp, ti]/(1 + Subscript[mg, fi]/Subscript[k, DDi]);
Subscript[adp, mi] := Subscript[adp, ti] - Subscript[adp, fi];
Subscript[amp, e] := 
  Subscript[a, eSUM] - Subscript[atp, te][t] - Subscript[adp, te][t];
cr := Subscript[c, SUM] - pcr[t];
Subscript[e, mN] := Subscript[e, mN0] + ((Z/2)*Log[10, (nad/nadh[t])]);
Subscript[e, mU] := Subscript[e, mU0] + ((Z/2)*Log[10, (uq/uqh[t])]);
Subscript[e, mc] := 
  Subscript[e, mc0] + Z*Log[10, (c3plus/c2plus[t])];
Subscript[e, ma] := Subscript[e, mc] + deltap*(2 + (2*u))/2;
Subscript[a, 3 over2] := 
  10^((Subscript[e, ma] - Subscript[e, ma0])/Z);(*a3plus/a2plus ratio*)
a2plus := Subscript[a, t]/(1 + Subscript[a, 3 over2]);
a3plus := Subscript[a, t] - a2plus;
Subscript[c, 
  0 i] := (10^-Subscript[pH, i] - 10^(-Subscript[pH, i] - dpH))/
   dpH;(*natural buffer capacity in matrix*)
Subscript[c, 
  0 e] := (10^-Subscript[pH, e] - 10^(-Subscript[pH, e] - dpH))/
   dpH;(*natural buffer capacity in cytosol*)
Subscript[r, buffi] := 
  Subscript[c, buffi]/Subscript[c, 
   0 i];(*buffering capacity coefficient for matrix*)
Subscript[pH, i] := -Log[10, (Subscript[h, i][t]/1000000)];
deltapH := Z*(Subscript[pH, i] - Subscript[pH, e]);
deltap := (1/(1 - u)) deltapH;
delta\[Psi] := -(deltap - deltapH);
Subscript[\[Psi], i] := 0.65*delta\[Psi];
Subscript[\[Psi], e] := -0.35*delta\[Psi];
Subscript[pi, je] := 
  Subscript[pi, te][
    t]/(1 + 10^(Subscript[pH, e] - Subscript[pK, a]));
Subscript[pi, ji] := 
  Subscript[pi, ti][
    t]/(1 + 10^(Subscript[pH, i] - Subscript[pK, a]));
nad := Subscript[n, t] - nadh[t];
uq := Subscript[ubq, t] - uqh[t];
c3plus := Subscript[c, t] - c2plus[t]

(****DIFFERENTIAL EQUATIONS****)
odeSolution =
 NDSolve[
  {nadh'[
     t] == (Subscript[v, DH] - Subscript[v, C1]) Subscript[r, cm]/
      Subscript[b, N],
   uqh'[t] == (Subscript[v, C1] - Subscript[v, C3]) Subscript[r, cm],
   c2plus'[t] == (Subscript[v, C3] - 2*Subscript[v, C4])*2*Subscript[
     r, cm],
   o2'[t] == 0,(*or o2'[t]\[Equal]-Subscript[v, C4]*)
   Subscript[h, i]'[
     t] == -((2*(2 + 2*u)*Subscript[v, C4]) + ((4 - 2*u)*Subscript[v, 
          C3]) + (4*Subscript[v, C1]) - (Subscript[n, A]*Subscript[v, 
          SN]) - (u*Subscript[v, EX]) - ((1 - u)*Subscript[v, 
          PI]) - (Subscript[v, LK]))*
     Subscript[r, cm]/Subscript[r, buffi],
   Subscript[atp, ti]'[t] == (Subscript[v, SN] - Subscript[v, EX])*
     Subscript[r, cm],
   Subscript[pi, ti]'[t] == (Subscript[v, PI] - Subscript[v, SN])*
     Subscript[r, cm],
   Subscript[atp, te]'[
     t] == (Subscript[v, EX] - Subscript[v, UT] + Subscript[v, AK] + 
      Subscript[v, CK]),
   Subscript[adp, te]'[
     t] == (Subscript[v, UT] - Subscript[v, 
      EX] - (2*Subscript[v, AK]) - Subscript[v, CK]),
   Subscript[pi, te]'[t] == (Subscript[v, UT] - Subscript[v, PI]),
   pcr'[t] == -Subscript[v, CK],
   nadh[0] == 822,(*starting values from BK*)
   uqh[0] == 1143,
   c2plus[0] == 60.89,
   o2[0] == 240,(*from BK (1991)*)
   Subscript[h, i][0] == 0.037,
   Subscript[atp, ti][0] == 6987,
   Subscript[pi, ti][0] == 7020,
   Subscript[atp, te][0] == 6668,
   Subscript[adp, te][0] == 31.7,
   Subscript[pi, te][0] == 2599,
   pcr[0] == 12219},
  {nadh[t], uqh[t], c2plus[t], o2[t], Subscript[h, i][t], 
   Subscript[atp, ti][t], Subscript[pi, ti][t], Subscript[atp, te][t],
    Subscript[adp, te][t], Subscript[pi, te][t], pcr[t]},
  {t, 0, 10}]