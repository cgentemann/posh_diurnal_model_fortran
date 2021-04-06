	program fluxes
c	c.gentemann 6/2009 gentemann@remss.com
c	gentemann, minnett, and ward, (in press), "Profiles of Ocean Surface Heating (POSH):a new model of upper ocean diurnal warming", JGR


	real hr_dawn,hr_dylen
      real*8 hUm,hTm,hUs,hTs
	real*8 ws_h,ta_h,qq_h,rf,tau,ustar,tstar,qstar,cd,ch,ce,rr,rt
	real*8 rq,zl,zo,zot,zoq,dter,t0,wg,tau_r,ef_webb,qh,qe
      real*8 ws,sst,atb,qq,pp,zi,rain
      real*8 rl,rs,ts_depth
      real*8 dt_wrm
      real*8 glat,glon
      real*8 intime,sol_time      
      real*8 time_old,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old
      real*8 tk_pwp,fxp
	real*4 precip
      integer mn,yy,dd,jamset,jump
      COMMON /old/time_old,qcol_ac,tau_ac,tau_old,rf_old,hf_old,
     &            ef_old,jamset,jump,fxp,tk_pwp   
      real wnd_in,air_in,q_in,pp_in,
	1	rs_in,rl_in,rain_in,xlat,xlon,sst_in,skin_in,rey_in                          
	real avwnd,avwndc,pkinso
	real pi,sun_theta,dt_warm_profile(100)
	real*4,dimension(24):: wnd_in1,air_in1,q_in1,pp_in1,rs_in1,rl_in1,rain_in1,sst_in1,skin_in1,rey_in1
	pi=3.14159
      qcol_ac=0.
      tau_ac=0.
      time_old=0.
      jamset=0.
      tau_old=0.
      hf_old=0.
      ef_old=0.
      rf_old=0.
      jamset=0
      jump=0
      hum=10.                 !height of wind measurement
      htm=10.                 !height of air temp. and RH
      hUs=10.                 !10m standard levels
      hTs=10.
      ts_depth=3.0            !6m data
      pp=1008.
      zi=600.
      fxp=0.5       
      tk_pwp=19.0  
      jwarm=1
      jcool=1
	jump=0;jamset=0;fxp=0.5;tk_pwp=19.0;tau_ac=0.0;qcol_ac=0.0;dt_wrm=0.0  
	idyjlsv=0

	open(unit=34,file='d:\diurnal\posh\posh_test_v2\input_data.txt',form='formatted')
	open(unit=35,file='d:\diurnal\posh\posh_test_v2\output_data.txt',form='formatted')
	open(unit=36,file='d:\diurnal\posh\posh_test_v2\output_data_profile.txt',form='formatted')

      qcol_ac=0.
      tau_ac=0.
      time_old=0.
      jamset=0.
      tau_old=0.
      hf_old=0.
      ef_old=0.
      rf_old=0.
      jamset=0
      jump=0
      hum=10.                 !height of wind measurement
      htm=10.                 !height of air temp. and RH
      hUs=10.                 !10m standard levels
      hTs=10.
      ts_depth=3.0            !6m data
      pp=1008.
      zi=600.
      fxp=0.5       
      tk_pwp=19.0  
      jwarm=1
      jcool=1
	jump=0;jamset=0;fxp=0.5;tk_pwp=19.0;tau_ac=0.0;qcol_ac=0.0;dt_wrm=0.0  
	idyjlsv=0

	avwnd=0;avwndc=0;

	do k=1,24
	read(34,'(6I5,20F12.2)') yy,mn,dd,ihh,imm,iss,wnd_in1(k),air_in1(k),q_in1(k),pp_in1(k),
 	1	rs_in1(k),rl_in1(k),rain_in1(k),xlat,xlon,sst_in1(k),skin_in1(k),rey_in1(k)
	avwnd=avwnd+wnd_in1(k)
	avwndc=avwndc+1
	enddo
	close(34)
	avwnd=avwnd/avwndc !some of the diurnal models input daily average wind, so that is calculated here

c	add some tests to make sure input sign conventions are correct
	if (maxval(q_in1)<0) then
		print*, 'possible sign error in q_in1, humidity (should be positive)'
		stop
	endif
	if (maxval(pp_in1)<0) then
		print*, 'possible sign error in pp_in1, pressure (should be positive)'
		stop
	endif
	if (maxval(rs_in1)<0) then
		print*, 'possible sign error in rs_in1, shortwave (in) (should be positive)'
		stop
	endif
	if (maxval(rl_in1)<0) then
		print*, 'possible sign error in rl_in1, longwave (in) (should be positive)'
		stop
	endif


	do ik=1,137 !the data is hourly, here I interpolate the data 	to ~10 minutes
	k1=floor(ik/6.)+1
	k2=k1+1
	f1=1.-(ik-(k1-1)*6.)/6.
	f2=1.-f1
	ihh=k1
	imm=(ik-6.*floor(ik/6.))*10
	wnd_in=f1*wnd_in1(k1)+f2*wnd_in1(k2)
	air_in=f1*air_in1(k1)+f2*air_in1(k2)
	q_in=f1*q_in1(k1)+f2*q_in1(k2)
	pp_in=f1*pp_in1(k1)+f2*pp_in1(k2)
	rs_in=f1*rs_in1(k1)+f2*rs_in1(k2)
	rl_in=f1*rl_in1(k1)+f2*rl_in1(k2)
	rain_in=0
 	sst_in=f1*sst_in1(k1)+f2*sst_in1(k2)
	skin_in=f1*skin_in1(k1)+f2*skin_in1(k2)
	rey_in=f1*rey_in1(k1)+f2*rey_in1(k2)

c	print*, yy,ihh,imm,iss,xhour,isecdy
c	print*, xlat,xlon

	if(xlon<0)xlon=xlon+360	  !make sure 0-360 rather than -180 to 180
	if(xlon<0)xlon=xlon+360
	lyr=yy
c	if(lyr<1980)lyr=yy+2000
	ws_test=wnd_in 
	pp=pp_in
	if(rl_in<0)rl_in=380      !longwave is missing sometimes, set to 380 if missing.  380 was common value for explorer data.

	call dy2jul(lyr,idyjl,mn,dd)  !this subroutine goes back and forth between month/day and day of year
c	print*, ikk,idyjl
      index=index+1              !count data records (hours) 

	u=ws_test
	ws=u;
	sst=sst_in
	atb=air_in
	q=q_in/100.	  !humidity input is 0-100, % of relative humidity
	qq=q;
	hsb=10.
	hlb=50.
	tub=0.03
	rs=rs_in
	rl=rl_in
	rain=rain_in
	ts=sst;                          
	glon=xlon;glat=xlat;

	xhour=ihh+imm/60.+iss/3600.
	isecdy=nint(xhour*3600.)

c	CLG here geometric cal of TOA insolation, dawn, and length of day.  dusk =hr_dawn+hr_dylen
	call fd_insolation_inst(lyr,idyjl,isecdy,xlat,xlon, hr_dawn, hr_dylen, rgeo)
      call fd_solar_elevation(lyr,idyjl,isecdy,xlat,xlon, sun_theta) !sun elevation in deg

	if(idyjlsv==0)idyjlsv=idyjl
	rlocal=xhour+(xlon/360)*24.
	if(rlocal>24)rlocal=rlocal-24

ccc
c	R.Weihs & CLG update, do not reset at 6AM, reset at 1 hr before hr_dawn when possible
c	as the length of day increases, changes, hr_dawn = 0 and hr_dylen=24
c	so when the length of day < 23 hours, reset 1 hr before dawn
c	when the lenght of day >=23 hours, reset at dawn
c	R.Weihs comment: careful here if you are not running a timeseries, but instead NWP maps.
c	dawn could occur on previous day
	if(idyjlsv==0)idyjlsv=idyjl  
	if(idyjl/=idyjlsv)then
		if(hr_dylen<23 .and. rlocal>hr_dawn-1)then
			jump=0;jamset=0;fxp=0.5;tk_pwp=19.0;tau_ac=0.0;qcol_ac=0.0;dt_wrm=0.0  
			idyjlsv=idyjl
		endif
		if(hr_dylen>=23 .and. rlocal>hr_dawn)then
			jump=0;jamset=0;fxp=0.5;tk_pwp=19.0;tau_ac=0.0;qcol_ac=0.0;dt_wrm=0.0  
			idyjlsv=idyjl
		endif
	endif
c end of code to reset ~dawn

	if(Jwarm.gt.0) then
		if(index.eq.1) then
		  Jwarm=2
		else
		  Jwarm=1
		endif
	endif
c
      intime=(float(ihh)+float(imm)/60.+iss/3600.) !eg 18.4833

	xhour=ihh+imm/60.+iss/3600.
	isecdy=nint(xhour*3600.)
	wnd=ws;
	pkinso=maxval(rs_in1)		 !peak insolation value used in one diurnal model
	precip=0;

c	print*, lyr,idyjl,isecdy
c	CALL SUNLOC1(lyr,idyjl,isecdy, SUNLAT,SUNLON,SUNDIS)
	CALL GET_DIURNAL_SST_MW(lyr,idyjl,XHOUR,xlat,XLON,wnd, DSST)   
      call get_diurnal_sst_models(pkinso,avwnd,precip,dsst2,dsst3,dsst4)
	call get_diurnal_sst_seviri2(lyr,idyjl,xhour,xlat,xlon,wnd, dsst5)

c and convert to local solar time in seconds

      sol_time=mod(glon/15+intime+24,24.)*3600   !eg 17580
c
c call bulk flux routine
c                    
	rwnd_diss=.00005
	rsol_diss=.00008
	ss=0

c	write(*,'(20F7.1)') sol_time,rs,rl,dt_wrm,dter,tk_pwp
      call bulk_flux(sol_time,glat,hUm,hTm,hUs,hTs,ws,sst,atb,qq,ws_h,Ta_h,qq_h,
     &   rs,rl,rain,pp,zi,Jcool,Jwarm,QH,QE,RF,TAU,Ustar,Tstar,Qstar,
     &   CD,CH,CE,RR,RT,RQ,ZL,ZO,zot,zoq,dt_wrm,dter,
     &   T0,ts_depth,wg,TAU_r,EF_webb,sun_theta,rwnd_diss,rsol_diss,
     &   dt_warm_profile)

	write(35,'(6I5,30F12.2)') yy,mn,dd,ihh,imm,iss,wnd_in,air_in,q_in,pp_in,
 	1	rs_in,rl_in,rain_in,xlat,xlon,sst_in,skin_in,rey_in,dt_wrm,dsst,dsst2,dsst3,dsst4,dsst5
      write(36,203)yy,mn,dd,ihh,imm,isecdy,tk_pwp,dt_wrm,dt_warm_profile

203   format(6i5,102f16.3)

	
	enddo

	close(35)
	close(36)
	print*,'done'

	end


 	include 'd:\diurnal\posh\posh_test_v2\subroutines\openbig.for'
	include 'd:\diurnal\posh\posh_test_v2\subroutines\bulf12.f'  !changed profile ot all exp & .0001 for wind diss & morning exp profile
	include 'd:\diurnal\posh\posh_test_v2\subroutines\asl.f'  
	include 'd:\diurnal\posh\posh_test_v2\subroutines\zeta.f'
	include 'd:\diurnal\posh\posh_test_v2\subroutines\h_adjust.f'
	include 'd:\diurnal\posh\posh_test_v2\subroutines\gravity.f'
 	INCLUDE 'd:\diurnal\posh\posh_test_v2\subroutines\SUNLOC1.F'            
 	INCLUDE	'd:\diurnal\posh\posh_test_v2\subroutines\GET_DIURNAL_SST6.F'	
 	INCLUDE	'd:\diurnal\posh\posh_test_v2\subroutines\GET_DIURNAL_clayson2.F'	
	include 'd:\diurnal\posh\posh_test_v2\subroutines\jul2dy.f'
	include 'd:\diurnal\posh\posh_test_v2\subroutines\zip.f'
      INCLUDE 'd:\diurnal\posh\posh_test_v2\subroutines\fd_insolation_inst.f'
      INCLUDE 'd:\diurnal\posh\posh_test_v2\subroutines\fd_solar_elevation.f'
      INCLUDE 'd:\diurnal\posh\posh_test_v2\subroutines\fd_insolation.f'

