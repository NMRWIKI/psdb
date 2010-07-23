# 1 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt"
;n_can_2h_ct_2d_ul.kt
;for uniform labeled sample
;cpd decoupling version
;avance-version (04/09/08)
;CaN nitrogen detection
;3D sequence with
;   15N detected correlation for double resonance
;
;      F2(Ca, t1) -> F1(N,t2)
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

"d21=11.3m" ;1/8JaN  delay for Ca evolution can be determined by paropt
"d22=14.3m" ;1/(2*JAB)
;"d23=11.3m" ;1/8JaN  delay for N evolution can be determined by paropt

;definition for constant time for t1
"d31=d22/2-p12/2" ;delay 1 for Ca constant time
"d30=d22-d21-p1" ;delay 2 for Ca constant time
"d29=d21-d22/2-p12/2-p1" ;delay 3 for Ca constant time
"d28=d22/2-p12/2" ;delay 4 for Ca constant time


;definition for nitrogen transverse delays
"d27=d23/2-p12/2" ;delay for N reforcus
"d26=d27-d12-4u"

	"in31=in0/4"
	"in30=0"
	"in29=in0/4"
	"in28=in0/4"


;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)


"spoff23=0"
"spoff24=0"
"spoff25=0"
"spoff26=bf3*((cnst21-cnst22)/1000000)"
"spoff15=bf3*((cnst21-cnst22)/2000000)"


# 1 "mc_line 64 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt expanding definition part of mc command before ze"
define delay MCWRK
define delay MCREST
define loopcounter ST1CNT
"ST1CNT = td1 / (2)"
"MCWRK = 0.200000*d11"
"MCREST = d11 - d11"
# 64 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt"
1 ze
# 1 "mc_line 64 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt expanding definition of mc command after ze"
# 65 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt"


;f1:N, f2:H, f3:C, f4:D

  d11 pl16:f3
  d11 setnmr4|28
  d11 pl17:f4
# 1 "mc_line 72 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt expanding start label for mc command"
2 MCWRK  do:f2 do:f3 do:f4
LBLSTS1, MCWRK  * 4
LBLF1, MCREST
# 73 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt"
3 d11 setnmr4^24
  9m setnmr3^0
  d1 pl1:f1 pl12:f2 pl0:f3
  50u setnmr3|0 setnmr0|34|32|33
  d12 setnmr4|24
  5u
  d12 fq1:f4  ;D-alpha
  d12 cpd4:f4
  5u
  (p11:sp23 ph3:r):f3 ; a90 Q5	;Ca indirect_____________
   d31
  (p12:sp26 ph1):f3 ;CO180 Q3
   d31
  (p12:sp24 ph1):f3 ; Ca180 Q3
   d30
  (p2 ph1):f1 ; N180
  d29
  (p12:sp26 ph1):f3 ;CO180 Q3
  d28
  (p11:sp25 ph2):f3 ; a90tr Q5tr	;________________________

  4u do:f4
  5u
  p16:gp1
  d16 fq1:f4  ;D-N
  20u cpd2:f2 cpd4:f4
  5u
  
  (p1 ph1):f1			;Ny i->i+/-1 transfer____
  d27
  (p12:sp26 ph1):f3
  d27
  (center (p2 ph1):f1 (p12:sp24 ph1):f3)
  d27
  (p12:sp26 ph1):f3
  d26
  d12 pl16:f3

  4u setnmr0^34^32^33
  go=2 ph31 cpd3:f3

# 1 "mc_line 114 file /opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt expanding mc command in line"
  MCWRK  do:f2 do:f3 do:f4 wr #0 if #0 zd ip3
  lo to LBLSTS1 times 2
  MCWRK dd31  MCWRK  id30  MCWRK  id29  MCWRK  id28
  lo to LBLF1 times ST1CNT
# 116 "/opt/xwinnmr/exp/stan/nmr/lists/pp/n_can_2h_ct_2d_ul.kt"



  d11 do:f2 do:f3 do:f4
  d11 setnmr4^24
  d11 setnmr3^0
  d11 setnmr4^28


exit


ph1=0
ph2=1
ph3=0 2
ph5=0
ph31=0 2 ;ph31 follow ph3


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
;d21: 11.3m 1/8JaN  delay for Ca evolution can be determined by paropt
;d22: 14.3m 1/(2*JAB)
;d23: 11.3m 1/8JaN  delay for N evolution can be determined by paropt
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
