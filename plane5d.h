/* Copyright (c) Signal Analysis and Imaging Group (SAIG), University of Alberta, 2013.*/
/* All rights reserved.                       */
/* plane5d  :  $Date: May 2013- Last version May    2013  */

#include "su.h"
#include "cwp.h"
#include "fftw3.h"

void plane5d(float **d,
             int nt,int nh,
             float *mx,float *my,float *hx,float *hy,
             int nevent,float *amp,float *t0,
             float *v_mx,float *v_my,float *v_hx,float *v_hy,
             float *curve_mx,float *curve_my,float *curve_hx,float *curve_hy,
             float dt,float fmin,float fmax,float f0);
void ricker_wavelet(float *w, float f,float dt);

void plane5d(float **d,
             int nt,int nh,
             float *mx,float *my,float *hx,float *hy,
             int nevent,float *amp,float *t0,
             float *v_mx,float *v_my,float *v_hx,float *v_hy,
             float *curve_mx,float *curve_my,float *curve_hx,float *curve_hy,
             float dt,float fmin,float fmax,float f0)
{
  complex shift;
  complex czero;
  int padfactor;
  int ntfft;
  int nw;
  int ih,it,iw;
  int ifmin,ifmax;
  complex Complex_I;
  complex **D;
  float omega;
  float* w;
  complex* W;
  float delay;
  float shift_mx, shift_my, shift_hx, shift_hy;
  float* in1;
  complex* out1;
  fftwf_plan p1;
  complex* in2;
  float* out2;
  fftwf_plan p2;
  int N;
  int ievent;
  
  czero.r=0;czero.i=0;  
  padfactor = 2;
  ntfft = npfar(padfactor*nt);
  nw = (int) ntfft/2 + 1; 
  /* fprintf(stderr,"nw=%d\n",nw);
  fprintf(stderr,"nh=%d\n",nh); */
  D = ealloc2complex(ntfft,nh);
  Complex_I.r = 0;
  Complex_I.i = 1;

  w = ealloc1float(ntfft);
  for (it=0;it<ntfft;it++) w[it] = 0;
  ricker_wavelet(w,f0,dt);
  
  W = ealloc1complex(nw);
  
  
  /**********************************************************************************************/
  N = ntfft;
  in1 = ealloc1float(N);
  out1 = ealloc1complex(nw);
  p1 = fftwf_plan_dft_r2c_1d(N, in1, (fftwf_complex*)out1, FFTW_ESTIMATE);
  for(it=0;it<ntfft;it++){
    in1[it] = w[it];
  }
  fftwf_execute(p1); 
  for(iw=0;iw<nw;iw++){
    W[iw] = out1[iw]; 
  }
  /**********************************************************************************************/
  delay = dt*(trunc((2.2/f0/dt))+1)/2;
  
  for (ih=0;ih<nh;ih++){ 
    for (iw=0;iw<nw;iw++){ 
      D[ih][iw] =  czero;
    }
  } 


  if(fmin>0){ 
    ifmin = (int) truncf(fmin*dt*ntfft);
  }
  else{
    ifmin = 0;
  }
  if(fmax*dt*ntfft<nw){ 
    ifmax = (int) truncf(fmax*dt*ntfft);
  }
  else{
    ifmax = 0;
  }

  for (ievent=0;ievent<nevent;ievent++){ 
    for (ih=0;ih<nh;ih++){ 
      for (iw=ifmin;iw<ifmax;iw++){
        omega = (float) 2*PI*iw/ntfft/dt;
        shift_mx = powf(fabsf(mx[ih]),curve_mx[ievent])/(v_mx[ievent]*fabsf(v_mx[ievent]));
        shift_my = powf(fabsf(my[ih]),curve_my[ievent])/(v_my[ievent]*fabsf(v_my[ievent]));
        shift_hx = powf(fabsf(hx[ih]),curve_hx[ievent])/(v_hx[ievent]*fabsf(v_hx[ievent]));
        shift_hy = powf(fabsf(hy[ih]),curve_hy[ievent])/(v_hy[ievent]*fabsf(v_hy[ievent]));
        /*
        fprintf(stderr,"shift_mx=%f shift_my=%f shift_hx=%f shift_hy=%f\n",shift_mx,shift_my,shift_hx,shift_hy);
		*/
        shift.r =   cos(omega*(t0[ievent] + shift_mx + shift_my + shift_hx + shift_hy - delay));
        shift.i =  -sin(omega*(t0[ievent] + shift_mx + shift_my + shift_hx + shift_hy - delay));
        D[ih][iw] = cadd(D[ih][iw],crmul(cmul(W[iw],shift),amp[ievent]));
      }
    } 
  }
  /**********************************************************************************************/
  N = ntfft;
  in2 = ealloc1complex(N);
  out2 = ealloc1float(N);
  p2 = fftwf_plan_dft_c2r_1d(N, (fftwf_complex*)in2, out2, FFTW_ESTIMATE);
  for (ih=0;ih<nh;ih++){
    for(iw=0;iw<nw;iw++){
      in2[iw] = D[ih][iw];
    }
    fftwf_execute(p2); 
    for(it=0;it<nt;it++){
      d[ih][it] = out2[it]; 
    }
  }
  /**********************************************************************************************/
  for (ih=0;ih<nh;ih++) for (it=0; it<nt; it++) d[ih][it]=d[ih][it]/ntfft;

  free2complex(D);
  fftwf_destroy_plan(p1);
  free1float(in1);free1complex(out1);
  free1complex(in2);free1float(out2);
  
  return; 
}
void ricker_wavelet(float *w, float f,float dt)
{
  int iw, nw, nc;
  float alpha, beta;  
  nw = (int) 2*trunc((float) (2.2/f/dt)/2) + 1;
  nc = (int) trunc((float) nw/2);
 
  for (iw=0;iw<nw-2;iw++){
    alpha = (nc-iw+1)*f*dt*PI;
  	beta = alpha*alpha;
    w[iw] = (1-beta*2)*exp(-beta);
  }

  return;
}


