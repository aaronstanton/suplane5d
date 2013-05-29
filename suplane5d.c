/* Copyright (c) Signal Analysis and Imaging Group (SAIG), University of Alberta, 2013.*/
/* All rights reserved.                       */
/* suplane5d  :  $Date:May     2013- Last version May 2013  */

#include "su.h"
#include "cwp.h"
#include "segy.h"
#include "header.h"
#include <time.h>
#include "plane5d.h"

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   							                                  ",
  " SUPLANE5D  model linear or curved events in 4 spatial dimensions  ",
  "                                                                   ",
  "           User provides:                                          ",
  "                                                                   ",
  "           Other parameters:                                       ",
  "           verbose=0; (=1 to show messages)                        ",
  "           dt=0.004; (sample rate)                                 ",
  "           tmax=1; (end time)                                      ",
  "           fmin = 0; (min frequency)                               ",
  "           fmax = 0.5/dt; (max frequency)                          ",
  "           f0 = 30; (dominant frequency of ricker wavelet)         ",
  "           nevent=1; (number of events)                            ",
  "                               ",
  "           To acheive shot/receiver geometry code:                    ",
  "           sxmin=50;           ",
  "           sxmax=50;           ",
  "           symin=50;           ",
  "           symax=50;           ",
  "           gxmin=50;           ",
  "           gxmax=50;           ",
  "           gymin=1;            ",
  "           gymax=100;          ",
  "           dsx=1;              ",
  "           dsy=1;              ",
  "           dgx=1;              ",
  "           dgy=1;              ",
  "                               ", 
  "           To acheive midpoint/offset geometry code:                  ", 
  "           mxmin=0;            ",
  "           mxmax=0;            ",
  "           mymin=0;            ",
  "           mymax=0;            ",
  "           hxmin=0;            ",
  "           hxmax=0;            ",
  "           hymin=0;            ",
  "           hymax=0;            ",
  "           dmx=1;              ",
  "           dmy=1;              ",
  "           dhx=1;              ",
  "           dhy=1;              ",
  "                               ",
  "                               ",
  "           For these parameters separate multiple values using a comma: ",
  "           amp = 1; (amplitude for each event)                          ",
  "           t0 = 0.2; (zero offset time for each event )                 ",
  "           v_mx = 10000000; (velocity for each event in mx direction)   ",
  "           v_my = 10000000; (velocity for each event in my direction)   ",
  "           v_hx = 5000; (velocity for each event in hx direction)       ",
  "           v_hy = 5000; (velocity for each event in hy direction)       ",
  "           curve_mx = 2; (power multiplied against the distance to achieve curvature in mx direction) ",
  "           curve_my = 2; (power multiplied against the distance to achieve curvature in my direction) ",
  "           curve_hx = 2; (power multiplied against the distance to achieve curvature in hx direction) ",
  "           curve_hy = 2; (power multiplied against the distance to achieve curvature in hy direction)",
  "                                                                   ",
  "           It is also possible to perturb the shot/receiver or midpoint/offset coordinates ",
  "           using the following parameters:                         ",
  "           when working with shot/receiver geometry code:                  ", 
  "           sx_std_dev = 0; (standard deviation of the sx corrdinate from the regular position) ",
  "           sy_std_dev = 0; (standard deviation of the sy corrdinate from the regular position) ",
  "           gx_std_dev = 0; (standard deviation of the gx corrdinate from the regular position) ",
  "           gy_std_dev = 0; (standard deviation of the gy corrdinate from the regular position) ",
  "           when working with midpoint/offset geometry code:                  ", 
  "           mx_std_dev = 0; (standard deviation of the mx corrdinate from the regular position) ",
  "           my_std_dev = 0; (standard deviation of the my corrdinate from the regular position) ",
  "           hx_std_dev = 0; (standard deviation of the hx corrdinate from the regular position) ",
  "           hy_std_dev = 0; (standard deviation of the hy corrdinate from the regular position) ",
  "                                                                   ",
  "           Example coding:                                         ",
  "           # makes a single 2D shot gather with irregular spacing       ",
  "           suplane5d verbose=1 gy_std_dev=5 > d.su                      ",
  "           suxwigb wclip=-1 bclip=1 key=gy < d.su &                     ",
 NULL};
/* Credits:
 * Aaron Stanton
 * Trace header fields accessed:
 * Last changes: May : 2013 
 */
/**************** end self doc ***********************************/

segy tr;

int main(int argc, char **argv)
{
  int verbose;
  time_t start,finish;
  double elapsed_time;
  int it,ih;
  int nt, ntr, nh; 
  int ievent,nevent;    /* number of events */
  float *amp=NULL;	    /* array of amplitudes for each event */
  float *t0=NULL;	    /* array of time at zero offset for each event */
  float *v_mx=NULL;	    /* array of dips in the mx direction for each event */
  float *v_my=NULL;	    /* array of dips in the my direction for each event */
  float *v_hx=NULL;	    /* array of dips in the hx direction for each event */
  float *v_hy=NULL;	    /* array of dips in the hy direction for each event */
  float *curve_mx=NULL;	/* array of curvatures in the mx direction for each event */
  float *curve_my=NULL;	/* array of curvatures in the my direction for each event */
  float *curve_hx=NULL;	/* array of curvatures in the hx direction for each event */
  float *curve_hy=NULL;	/* array of curvatures in the hy direction for each event */
  float dt,tmax;
  float f0,fmin,fmax;
  float **d;  
  float *sx,*sy,*gx,*gy,*mx,*my,*hx,*hy,*h,*az;
  float sx_std_dev,sy_std_dev,gx_std_dev,gy_std_dev,mx_std_dev,my_std_dev,hx_std_dev,hy_std_dev;
  float *sx_dev,*sy_dev,*gx_dev,*gy_dev,*mx_dev,*my_dev,*hx_dev,*hy_dev; 
  float sxmin,sxmax,dsx;
  float symin,symax,dsy;
  float gxmin,gxmax,dgx;
  float gymin,gymax,dgy;
  float mxmin,mxmax,dmx;
  float mymin,mymax,dmy;
  float hxmin,hxmax,dhx;
  float hymin,hymax,dhy;
  int isx,nsx;
  int isy,nsy;
  int igx,ngx;
  int igy,ngy;
  int imx,nmx;
  int imy,nmy;
  int ihx,nhx;
  int ihy,nhy;
  int mode;
  
  fprintf(stderr,"*******SUPLANE5D*********\n");
  /* Initialize */
  initargs(argc, argv);
  requestdoc(0); /* stdin not used */

  start=time(0);    
  /* Get parameters */
  if (!getparint("verbose", &verbose)) verbose=0;
  if (!getparint("ntr",&ntr)){ 
    ntr = 1000000;
    if (verbose) fprintf(stderr,"warning: ntr paramater not set; using ntr=%d\n",ntr);
  }
  if (!getparfloat("dt", &dt)) dt=0.004;
  if (!getparfloat("tmax", &tmax)) tmax=1;
  nt = (int) (tmax/dt) + 1;
  if (!getparfloat("fmin",&fmin)) fmin = 0;
  if (!getparfloat("fmax",&fmax)) fmax = 0.5/dt;
  if (!getparfloat("sx_std_dev",&sx_std_dev)) sx_std_dev = 0;
  if (!getparfloat("sy_std_dev",&sy_std_dev)) sy_std_dev = 0;
  if (!getparfloat("gx_std_dev",&gx_std_dev)) gx_std_dev = 0;
  if (!getparfloat("gy_std_dev",&gy_std_dev)) gy_std_dev = 0;
  if (!getparfloat("mx_std_dev",&mx_std_dev)) mx_std_dev = 0;
  if (!getparfloat("my_std_dev",&my_std_dev)) my_std_dev = 0;
  if (!getparfloat("hx_std_dev",&hx_std_dev)) hx_std_dev = 0;
  if (!getparfloat("hy_std_dev",&hy_std_dev)) hy_std_dev = 0;
  fmax = MIN(fmax,0.5/dt);
  if (!getparfloat("f0",&f0)) f0 = 30; /* 30 Hz dominant freq for ricker wavelet */
  if (!getparint("nevent",&nevent)){ 
    nevent=1; 
    amp = ealloc1float(nevent); amp[0] = 1;
    t0 = ealloc1float(nevent); t0[0] = 0.2;
    v_mx = ealloc1float(nevent); v_mx[0] = 10000000;
    v_my = ealloc1float(nevent); v_my[0] = 10000000;
    v_hx = ealloc1float(nevent); v_hx[0] = 5000;
    v_hy = ealloc1float(nevent); v_hy[0] = 5000;
    curve_mx = ealloc1float(nevent); curve_mx[0] = 2;
    curve_my = ealloc1float(nevent); curve_my[0] = 2;
    curve_hx = ealloc1float(nevent); curve_hx[0] = 2;
    curve_hy = ealloc1float(nevent); curve_hy[0] = 2;
  }
  else{
    if (!(nevent == countparval("amp"))) err("must give amp= vector");
    amp = ealloc1float(nevent);	getparfloat("amp", amp);
    
    if (!(nevent == countparval("t0"))) err("must give t0= vector");
    t0 = ealloc1float(nevent);	getparfloat("t0", t0);

    if (!(nevent == countparval("v_mx"))) err("must give v_mx= vector");
    if (!(nevent == countparval("v_my"))) err("must give v_my= vector");
    if (!(nevent == countparval("v_hx"))) err("must give v_hx= vector");
    if (!(nevent == countparval("v_hy"))) err("must give v_hy= vector");
    v_mx = ealloc1float(nevent);	getparfloat("v_mx", v_mx);
    v_my = ealloc1float(nevent);	getparfloat("v_my", v_my);
    v_hx = ealloc1float(nevent);	getparfloat("v_hx", v_hx);
    v_hy = ealloc1float(nevent);	getparfloat("v_hy", v_hy);
    
    if (!(nevent == countparval("curve_mx"))) err("must give curve_mx= vector");
    if (!(nevent == countparval("curve_my"))) err("must give curve_my= vector");
    if (!(nevent == countparval("curve_hx"))) err("must give curve_hx= vector");
    if (!(nevent == countparval("curve_hy"))) err("must give curve_hy= vector");
    curve_mx = ealloc1float(nevent);	getparfloat("curve_mx", curve_mx);
    curve_my = ealloc1float(nevent);	getparfloat("curve_mx", curve_mx);
    curve_hx = ealloc1float(nevent);	getparfloat("curve_hx", curve_hx);
    curve_hy = ealloc1float(nevent);	getparfloat("curve_hy", curve_hy);
  }
  if (verbose){
    for (ievent=0;ievent<nevent;ievent++){ 
      fprintf(stderr,"amp[%d]      =%f\n",ievent,amp[ievent]);
      fprintf(stderr,"t0[%d]       =%f\n",ievent,t0[ievent]);
      fprintf(stderr,"v_mx[%d]   =%f\n",ievent,v_mx[ievent]);
      fprintf(stderr,"v_my[%d]   =%f\n",ievent,v_my[ievent]);
      fprintf(stderr,"v_hx[%d]   =%f\n",ievent,v_hx[ievent]);
      fprintf(stderr,"v_hy[%d]   =%f\n",ievent,v_hy[ievent]);
      fprintf(stderr,"curve_mx[%d] =%f\n",ievent,curve_mx[ievent]);
      fprintf(stderr,"curve_my[%d] =%f\n",ievent,curve_my[ievent]);
      fprintf(stderr,"curve_hx[%d] =%f\n",ievent,curve_hx[ievent]);
      fprintf(stderr,"curve_hy[%d] =%f\n",ievent,curve_hy[ievent]);
    } 
  }
  
  if (!getparfloat("sxmin", &sxmin)) sxmin=50;
  if (!getparfloat("sxmax", &sxmax)) sxmax=50;
  if (!getparfloat("symin", &symin)) symin=50;
  if (!getparfloat("symax", &symax)) symax=50;
  if (!getparfloat("gxmin", &gxmin)) gxmin=50;
  if (!getparfloat("gxmax", &gxmax)) gxmax=50;
  if (!getparfloat("gymin", &gymin)) gymin=1;
  if (!getparfloat("gymax", &gymax)) gymax=100;
  if (!getparfloat("dsx", &dsx)) dsx=1;
  if (!getparfloat("dsy", &dsy)) dsy=1;
  if (!getparfloat("dgx", &dgx)) dgx=1;
  if (!getparfloat("dgy", &dgy)) dgy=1;
  /*  */
  if (!getparfloat("mxmax", &mxmax)) mxmax=0;
  if (!getparfloat("mymin", &mymin)) mymin=0;
  if (!getparfloat("mymax", &mymax)) mymax=0;
  if (!getparfloat("hxmin", &hxmin)) hxmin=0;
  if (!getparfloat("hxmax", &hxmax)) hxmax=0;
  if (!getparfloat("hymin", &hymin)) hymin=0;
  if (!getparfloat("hymax", &hymax)) hymax=0;
  if (!getparfloat("dmx", &dmx)) dmx=1;
  if (!getparfloat("dmy", &dmy)) dmy=1;
  if (!getparfloat("dhx", &dhx)) dhx=1;
  if (!getparfloat("dhy", &dhy)) dhy=1;
  
  mode = 1;
  nsx = (int) truncf((sxmax - sxmin)/dsx) + 1;
  nsy = (int) truncf((symax - symin)/dsy) + 1;
  ngx = (int) truncf((gxmax - gxmin)/dgx) + 1;
  ngy = (int) truncf((gymax - gymin)/dgy) + 1;
  nmx = (int) truncf((mxmax - mxmin)/dmx) + 1;
  nmy = (int) truncf((mymax - mymin)/dmy) + 1;
  nhx = (int) truncf((hxmax - hxmin)/dhx) + 1;
  nhy = (int) truncf((hymax - hymin)/dhy) + 1;
  if (!getparfloat("mxmin",&mxmin)){
  	mode = 1;
    mxmin = 0;
    if (verbose) fprintf(stderr,"using surface coordinates (sx,sy,gx,gy)\n");
  }
  else{
  	mode = 2;
    if (verbose) fprintf(stderr,"using subsurface coordinates (mx,my,hx,hy)\n");
  }

  if (mode==1) nh = nsx*nsy*ngx*ngy;
  else         nh = nmx*nmy*nhx*nhy;
  /* Allocate memory for data */
  /* fprintf(stderr,"nt=%d\n",nt);
  fprintf(stderr,"nh=%d\n",nh); */

  d  = ealloc2float(nt,nh);
  sx = ealloc1float(nh);
  sy = ealloc1float(nh);
  gx = ealloc1float(nh);
  gy = ealloc1float(nh);
  mx = ealloc1float(nh);
  my = ealloc1float(nh);
  hx = ealloc1float(nh);
  hy = ealloc1float(nh);
  h  = ealloc1float(nh);
    
  sx_dev = ealloc1float(nsx);
  sy_dev = ealloc1float(nsy);
  gx_dev = ealloc1float(ngx);
  gy_dev = ealloc1float(ngy);
  mx_dev = ealloc1float(nmx);
  my_dev = ealloc1float(nmy);
  hx_dev = ealloc1float(nhx);
  hy_dev = ealloc1float(nhy);
  
  if (mode==1){
  	for (isx=0;isx<nsx;isx++) sx_dev[isx] = sx_std_dev*frannor();
  	for (isy=0;isy<nsy;isy++) sy_dev[isy] = sy_std_dev*frannor();
  	for (igx=0;igx<ngx;igx++) gx_dev[igx] = gx_std_dev*frannor();
  	for (igy=0;igy<ngy;igy++) gy_dev[igy] = gy_std_dev*frannor();
    ih = 0;   
    for (isx=0;isx<nsx;isx++){ 
      for (isy=0;isy<nsy;isy++){ 
        for (igx=0;igx<ngx;igx++){ 
          for (igy=0;igy<ngy;igy++){ 
            sx[ih] = isx*dsx + sxmin + sx_dev[isx];
            sy[ih] = isy*dsy + symin + sy_dev[isy]; 
            gx[ih] = igx*dgx + gxmin + gx_dev[igx]; 
            gy[ih] = igy*dgy + gymin + gy_dev[igy]; 
            ih++;
          }
        }
      }
    } 
  } 
  else{ 
  	for (imx=0;imx<nmx;imx++) mx_dev[imx] = mx_std_dev*frannor();
  	for (imy=0;imy<nmy;imy++) my_dev[imy] = my_std_dev*frannor();
  	for (ihx=0;ihx<nhx;ihx++) hx_dev[ihx] = hx_std_dev*frannor();
  	for (ihy=0;ihy<nhy;ihy++) hy_dev[ihy] = hy_std_dev*frannor();
    ih = 0;   
    for (imx=0;imx<nmx;imx++){ 
      for (imy=0;imy<nmy;imy++){ 
        for (ihx=0;ihx<nhx;ihx++){ 
          for (ihy=0;ihy<nhy;ihy++){ 
            mx[ih] = imx*dmx + mxmin + mx_dev[imx];
            my[ih] = imy*dmy + mymin + my_dev[imy]; 
            hx[ih] = ihx*dhx + hxmin + hx_dev[ihx]; 
            hy[ih] = ihy*dhy + hymin + hy_dev[ihy]; 
            ih++;
          }
        }
      }
    } 
  } 
  if (mode==1){ 
    for (ih=0;ih<nh;ih++){
      mx[ih] = (gx[ih] + sx[ih])/2;
      my[ih] = (gy[ih] + sy[ih])/2;
      hx[ih] = gx[ih] - sx[ih];
      hy[ih] = gy[ih] - sy[ih];      
      h[ih]  = sqrt(hx[ih]*hx[ih] + hy[ih]*hy[ih]);      
    } 
  } 
  else{ 
    for (ih=0;ih<nh;ih++){
      sx[ih] = mx[ih] - 0.5*hx[ih];
      sy[ih] = my[ih] - 0.5*hy[ih];
      gx[ih] = sx[ih] + hx[ih];
      gy[ih] = sy[ih] + hy[ih];
      h[ih]  = sqrt(hx[ih]*hx[ih] + hy[ih]*hy[ih]);      
    } 
  } 
  
  if (verbose) fprintf(stderr,"creating %d traces \n", nh);
  for (ih=0;ih<nh;ih++){ 
    for (it=0;it<nt;it++){   
      d[ih][it] = 0;
    }
  }
  plane5d(d,
          nt,nh,
          mx,my,hx,hy,
          nevent,amp,t0,
          v_mx,v_my,v_hx,v_hy,
          curve_mx,curve_my,curve_hx,curve_hy,
          dt,fmin,fmax,f0);
  
 /* ***********************************************************************
 outputting data:
 *********************************************************************** */
  for (ih=0;ih<nh;ih++){ 
    memcpy((void *) tr.data,(const void *) d[ih],nt*sizeof(float));
    tr.ntr = nh;
    tr.sx  = (int) sx[ih];
    tr.sy  = (int) sy[ih];
    tr.gx  = (int) gx[ih];
    tr.gy  = (int) gy[ih];
    tr.gelev  = (int) mx[ih];
    tr.selev  = (int) my[ih];
    tr.gdel   = (int) hx[ih];
    tr.sdel   = (int) hy[ih];
    tr.offset = (int) h[ih];
    tr.ns = nt;
    tr.dt = NINT(dt*1000000.);
    tr.tracl = tr.tracr = ih + 1;
    fputtr(stdout,&tr);
  }


 /* ***********************************************************************
 end outputting data
 *********************************************************************** */
 free2float(d);

 /******** End of output **********/
 finish=time(0);
 elapsed_time=difftime(finish,start);
 fprintf(stderr,"Total time required: %6.2f \n", elapsed_time);
 
 return EXIT_SUCCESS;
}




