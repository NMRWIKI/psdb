#ifndef LINT
#endif

/*

Restore, but use V13 to avoid clash qwith pahse cycling
REVERT then add d3 and d4
GAM shift J maximum half a chunk later by unbalancing halves of t1
GAM swap order of 180s
GAM use gzlvl1/gt1 for CTP pulses
GAM 21viii07 add a little diffusion weighting
GAM version 21viii07 : use 2D-style acquisition with more secure timing
GAM version 24 iii 06 :  start data collection one dwell time early by reducing gstab6

Pulse sequence for collecting a homodecoupled (pureshift) spectrum
Refs:
Sterk 1996 ?????



Instructions:
Set up an array of d2 values.

np_fid - number of points per fid
at_fid - acquisition time used per fid
dw     - dwell time (1/sw*2)

d2 should be arrayed in at_fid

at_fid = dw * np_fid

example: using 24 datapoints per fid at spectral width (sw) of 3000 Hz
at_fid=1/sw * 48= 16 ms

The d2 array should be set up as:
16, 32, 48 etc to the desired fid lengt in the resulting assembled spectrum


pulse seq:
soft270--A--hard180--B--soft180--C--acquire

the delays A,B,C must be set up so that A+C=B to ensure correct phase.



MN 26iv06 now aquiring before refocusing to half the amount of d2 values necessary
          of the large frequence dependent phase errors.
MN 20Nov06  Adding phase cycling
MN 20Nov06  setting the delays so that A+C=B is automatically true
MN 28Jan07  updated abort to abort_message
MN 15Feb07  updated pbox_pulse to use &

Juan Aguilar 05_2008  using hard90 instead of soft90
Juan Aguilar 05_2008  soft180 pulse generation supports new Varian soft pulse generation recommendations



  */

#include <standard.h>
#include <Pbox_psg.h>


pulsesequence()
{
double 	
	gzlvl1=getval("gzlvl1"),		/* diffusion en/decoding */
	gzlvl2=getval("gzlvl2"),		/* homdecouple level*/
	gzlvl_ss=getval("gzlvl_ss"),
	gt1=getval("gt1"),				/*gradient encoding time*/
	corr=getval("corr"),				/*gradient encoding time*/
  	dosyfrq=getval("sfrq"),
	dosytimecubed,
	Ddelta,
    	gstab1=getval("gstab1"),
    	gstab2=getval("gstab2"),
    	gstab5=getval("gstab5"),
    	gstab6=getval("gstab6"),
    	gdel=getval("gdel"),
    	gtss=getval("gtss"),
    	del = getval("del"),
    	d3 = getval("d3"),
    	d4 = getval("d4"),
	selpw=getval("selpw"),		/* pulse length for the soft 180  */
	selpwr=getval("selpwr"),	/* power level  for the soft 180  */
    	nchunk = getval("nchunk")
	;

char delflag[MAXSTR],
     sspul[MAXSTR],
     selshape[MAXSTR],		/* pulse file for the soft 180 */
     lockcomp[MAXSTR];
	

getstr("delflag",delflag);
getstr("sspul",sspul);
getstr("lockcomp",lockcomp);
getstr("selshape",selshape);


/*Check conditions*/
  if ( (gzlvl2*selpw/gdel)> 32767) /*should be generalized for other systems*/
   {
      printf("Either gdel is too short or gzlvl2 too strong\n");
      abort_message("Either gdel is too short or gzlvl2 too strong\n");
   }


	/*phase cycling*/

        mod2(ct,v5);
        dbl(v5,v5);
        hlv(ct,v6);
        hlv(v6,v6);
        hlv(v6,v7);
        hlv(v7,v8);
        mod2(v6,v6);
        add(v5,v6,v1);   /* 270-pulse 0202 1313*/
        mod4(v1,v1);
        hlv(ct,v5);
        mod2(v5,v5);
        mod2(v7,v7);
        dbl(v7,v7);
        add(v5,v7,v3);   /* second 180  0011 0011 2233 2233*/
        mod4(v8,v2);      /* first 180   (0)16 (1)16 (2)16 (3)16*/
        dbl(v2,v9);
        add(v1,v9,v10);
        dbl(v3,v9);
        add(v10,v9,v11);
        mod4(v11,oph);    /* oph = v1 + 2*v2 + 2*v3 */

	del=2.0*(gstab2+gt1)+2.0*(gt1+pw)+rof1+rof2+gstab6+gstab1+2.0*gstab5+selpw;
/* Need to check  */

    Ddelta=gt1;/*the diffusion-encoding pulse width is gt1*/
	dosytimecubed=Ddelta*Ddelta*(del - (Ddelta/3.0)); 
    putCmd("makedosyparams(%e,%e)\n",dosytimecubed,dosyfrq);



   /* equilibrium period */
   status(A);


	if (sspul[0] == 'y')
      {  
		 zgradpulse(gzlvl_ss,gtss);	
		 obspower(tpwr);
         rgpulse(pw, zero, 0.0, 0.0);
         rgpulse(pw, one , 0.0, 0.0);
		 zgradpulse(gzlvl_ss,gtss);
      }  
	delay(d1);
	

   status(B);

	obspower(tpwr);					/* POWER CHANGE (HARD) */
	delay(rof1);					/* POWER CHANGE (HARD) */
	
    	rgpulse(pw, v1, rof1, rof2);			/* HARD 90 */
   
	
	
	if (d2_index>0)					/* INCREMENTED DELAY */
	{						/* INCREMENTED DELAY */
	initval(d2_index*nchunk*0.25+0.1,v12);		/* INCREMENTED DELAY */
	starthardloop(v12);				/* INCREMENTED DELAY */
		delay(1.0/sw);				/* INCREMENTED DELAY */
	endhardloop(); 					/* INCREMENTED DELAY */
	}						 
	
	zgradpulse(gzlvl1*0.5,gt1);			/* diffusion encoding and CTP selection */
	obspower(selpwr);				/* POWER CHANGE (SOFT) */
	delay(gstab2);					/* */
	
	rgradient('z',gzlvl2);					/* SLICE selection*/
	shaped_pulse(selshape,selpw, v3, rof1, rof2);	   	/* SLICE selection - Soft 180 */
	rgradient('z',0.0);					/* SLICE selection*/
	
	delay(gt1/2.0+gstab2);

	zgradpulse(-gzlvl1*0.5,gt1);				/* second CTP pulse*/
	obspower(tpwr);						/* Power level change (HARD*/
	delay(gt1/2.0+gstab5+gstab6);
   	rgpulse(pw*2,v2,rof1,rof2);				/* Hard 180*/
	delay(gstab5);
	zgradpulse(-gzlvl1,gt1);				/* third CTP pulse*/
	delay(gstab6-corr-1.0/sw);
delay(d4);							/* */
	if (d2_index>1)
	{
	initval(d2_index*nchunk*0.25-nchunk*0.25+0.1,v13);
	starthardloop(v13);
		delay(1.0/sw);
	endhardloop(); 
	}
	
/* Start acquisition before complete refocusing */

   /* --- observe period --- */
   status(C);
}
