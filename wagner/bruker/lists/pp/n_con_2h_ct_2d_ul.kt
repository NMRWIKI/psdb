# 1 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt"
;n_con_2h_ct_2d_ul.kt
;for uniform labeled sample
;cpd decoupling version
;avance-version (04/09/08)
;CaN nitrogen detection
;3D sequence with
;   15N detected correlation for double resonance
;
;      F2(CO, t1) -> F1(N,t2)
;
;on/off resonance 13C pulses using shaped pulses
;phase sensitive (t1)
;with 2H decoupling during relaxation delay
;(use parameterset )
;
;coded K Takeuchi 121808, tested OK


prosol relations=<triple_c>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"


"d11=30m"
"d12=20u"

;"d21=16.7m" ;1/4JCON total delay for CO evolution can be determined by paropt
;"d22=16.7m" ;1/4JCON total delay for N evolution can be determined by paropt


;definition for constant time for t1
"d31=d21/2-p12/2" ;delay 1 for CO constant time
"d30=d21/2-p12/2" ;delay 2 for CO constant time


;definition for nitrogen transverse delays
"d27=d22/2-p12/2" ;delay for N reforcus
"d26=d27-d12"

	"in31=in0/4"
	"in30=in0/4"


;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)


"spoff23=0"
"spoff24=0"
"spoff25=0"
"spoff26=bf3*((cnst22-cnst21)/1000000)"
"spoff15=bf3*((cnst22-cnst21)/2000000)"


# 1 "mc_line 60 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt expanding definition part of mc command before ze"
define delay MCWRK
define delay MCREST
define loopcounter ST1CNT
"ST1CNT = td1 / (2)"
"MCWRK = 0.333333*d11"
"MCREST = d11 - d11"
# 60 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt"
1 ze
# 1 "mc_line 60 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt expanding definition of mc command after ze"
# 61 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt"


;f1:N, f2:H, f3:C, f4:D

  d11 pl16:f3
  d11 setnmr4|28
  d11 pl17:f4
# 1 "mc_line 68 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt expanding start label for mc command"
2 MCWRK  do:f2 do:f3 do:f4
LBLSTS1, MCWRK  * 2
LBLF1, MCREST
# 69 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt"
3 d11 setnmr4^24
  9m setnmr3^0
  d1 pl1:f1 pl12:f2 pl0:f3
  50u setnmr3|0 setnmr0|34|32|33
  d12 setnmr4|24
  5u
  d12
  (p11:sp23 ph3):f3 ; CO90 Q5	;CO indirect_____________
   d31
  (p12:sp26 ph1):f3 ;Ca180 Q3
   d31
  (center (p2 ph1):f1 (p12:sp24 ph1):f3) ; N180 CO180 Q3
   d30
  (p12:sp26 ph1):f3 ;Ca180 Q3
  d30
  (p11:sp25 ph2):f3 ; CO90tr Q5tr	;________________________

  4u
  p16:gp1
  d16
  20u cpd2:f2 cpd4:f4

  (p1 ph1):f1			;Ny i->i+/-1 transfer____
  d27
  (p12:sp26 ph1):f3  ;Ca180
  d27
  (center (p2 ph1):f1 (p12:sp24 ph1):f3) ;N180 CO180
  d27
  (p12:sp26 ph1):f3 ;Ca180
  d26
  d12 pl16:f3

  4u setnmr0^34^32^33
  go=2 ph31 cpd3:f3

# 1 "mc_line 104 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt expanding mc command in line"
  MCWRK  do:f2 do:f3 do:f4 wr #0 if #0 zd ip3
  lo to LBLSTS1 times 2
  MCWRK dd31  MCWRK  id30
  lo to LBLF1 times ST1CNT
# 106 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_con_2h_ct_2d_ul.kt"



  d11 do:f2 do:f3 do:f4
  d11 setnmr4^24
  d11 setnmr3^0
  d11 setnmr4^28


exit


ph1=0
ph2=1
ph3=0 2
ph5=0
ph31=1 3 ;ph31 follow ph3


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;pl16: f3 channel - power level for CPD/BB decoupling
;sp3: f2 channel - shaped pulse 180 degree (adiabatic)
;sp23: f1 channel - shaped pulse  90 degree  (on resonance)
;sp24: f1 channel - shaped pulse 180 degree  (on resonance)
;sp25: f1 channel - shaped pulse  90 degree  (on resonance)
;                   for time reversed pulse
;sp26: f1 channel - shaped pulse 180 degree  (C=O off resonance)
;sp27: f1 channel - shaped pulse 180 degree  (Ca off resonance)
;sp28: f1 channel - shaped pulse 180 degree  (Ca on resonance)
;p11: f1 channel -  90 degree shaped pulse
;p12: f1 channel - 180 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;d0 : incremented delay (F1 in 2D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d23: 1/(2J(NCa))                                      [25 msec]
;d22: 1/4JCON total delay for N evolution can be determined by paropt
;d21: 1/4JCON total delay for CO evolution can be determined by paropt
;o1p: Calpha chemical shift (cnst22)
;in0: 1/(1 * SW(Ca)) =  DW(Ca)
;nd0: 1
;NS: 8 * n
;DS: >= 32
;td1: number of experiments in F1
;FnMODE: States-TPPI (or TPPI) in F1
;cpd2: decoupling according to sequence defined by cpdprg2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)



;$Id: c_can_mq.2,v 1.1.4.1 2004/11/23 15:08:14 ber Exp $
