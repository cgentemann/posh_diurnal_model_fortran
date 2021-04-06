      subroutine get_diurnal_sst_models(rgeo,wnd,precip,diurnal_sst,diurnal_sst2,diurnal_sst3)
      implicit none

c	c.gentemann gentemann@remss.com 6/2009

c	returns diurnal_sst = webster96
c			diurnal_sst2=kk02 skin
c			diurnal_sst3=kk02 bulk
      real(4), parameter :: pi=3.141592654
	real(4), parameter :: solar_constant=1360   !watts/m**2
	real(4) wnd,precip,peak_solar,rgeo
      real(4) diurnal_sst
      real(4) a0,b0,c0,d0,e0,f0	 
	real diurnal_sst2,diurnal_sst3

c     calculate solar insolation
 
	peak_solar=rgeo
c	clayson
	if(wnd>2)then
	a0=0.00265
	b0=0.028
	c0=-0.838
	d0=-0.00105
	e0=0.158
	f0=0.262
      else
	a0=0.002
	b0=0.041
	c0=0.212
	d0=-0.000185
	e0=-0.329
	f0=0.328
	endif

	if(wnd<.01)wnd=.01
	diurnal_sst=f0+a0*peak_solar+b0*precip+c0*log(wnd)+d0*peak_solar*log(wnd)+e0*wnd
	if(rgeo<=0)diurnal_sst=0

c skin algorithm kk02
cu>2.5ms-1	3.2708x10-6	-0.079982	-1.3329x10-6 	0.073287
cu<=2.5ms-1	5.6814x10-6 	0.40052	-3.9637x10-6 	-0.367
	if(wnd>2.5)then
	a0=3.2708e-6
	b0=-0.079982
	c0=-1.3329e-6
	d0=0.073287
c	a0=0.000003049;b0=-0.028258;c0=-0.0000011987;d0=-.025893
      else
	a0=5.6814e-6
	b0=0.40052
	c0=-3.9637e-6
	d0=-0.367
c	a0=0.0000050109;b0=0.22063;c0=-0.0000033394;d0=-0.20216
	endif
	diurnal_sst2=a0*peak_solar**2+b0*log(wnd)+c0*peak_solar**2*log(wnd)+d0
	if(diurnal_sst2<0)diurnal_sst2=0;
	if(rgeo<=0)diurnal_sst2=0

c bulk algorithm kko2
c	u>2.5ms-1	2.3989x10-6	0.057289	-9.2463x10-7 	-0.14236
c	u<=2.5ms-1	1.9361x10-6 	0.014576	-4.1966x10-7 	-.10322

	if(wnd>2.5)then
	a0=2.3989e-6
	b0=0.057289
	c0=-9.2463e-7
	d0=-0.14236
      else
	a0=1.9361e-6
	b0=0.014576
	c0=-4.1966e-7
	d0=-.10322
	endif
	diurnal_sst3=a0*peak_solar**2+b0*log(wnd)+c0*peak_solar**2*log(wnd)+d0
	if(diurnal_sst3<0)diurnal_sst3=0;
	if(rgeo<=0)diurnal_sst3=0

	RETURN
	END
