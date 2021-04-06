c*****************************************************************
c	c.gentemann gentemann@remss.com 6/2009

c	cg change, no reset at midnight
c	to change so that there was no reset at midnight i commented out
c	lines in this program and then added lines in the main program.
c	it wasn't possible to put the reset easily within this code, so it 
c	is external and ***must be added to main program***
c	the second fix was to fix dtime so that when the day was
c	crossed the dtime calcualtion still made sense 


c	qcol_ac*.96 if wind >2, increased absorption
c	
c	equations numbers next to some code references paper:
c	gentemann, minnett, and ward, (in press), "Profiles of Ocean Surface Heating (POSH):a new model of upper ocean diurnal warming", JGR


      subroutine bulk_flux(sol_time,glat,hUm,hTm,hUs,hTs,                  
     & ws,sst,atb,qq,ws_h,Ta_h,qq_h,rs,rl,rainx,pp,zix,Jcoolx,Jwarmx,
     & HF,EF,RF,TAU,Ustar,Tstar,Qstar,
     & CD,CH,CE,RRx,RTx,RQx,ZLx,ZOx,zotx,zoqx,
     & dt_wrmx,dterx,T0,ts_depthx,wgx,TAU_r,Hl_webb,sun_theta,rwnd_diss,rsol_diss,dt_warm_profile)
c
	real rwt,rws
      real*8 hUm,hTm,hUs,hTs,ws_h,Ta_h,qq_h,ts_depthx
      real*8 ws,sst,atb,qq,pp,zix,rainx,glat
      real*8 HF,EF,TAU,Ustar,Qstar,Tstar,TAU_r,Hl_webb
      real*8 rl,rs,RF,T0,sol_time,wgx,ee
      real*8 CD,CE,CH,RRx,RTx,RQx,Zlx,ZOx
      real*8 zotx,zoqx,dt_wrmx,dterx,qjoule,qr_out,dtime      
      real*8 ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,ctd1,ctd2,rich
      real*8 dqs_dt,dwat,dtmp
c
      real*8 U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,QA
      real*8 USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      real*8 al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,visa
      real*8 visw,von,fdg,DU_Wg,Wg,newtime      
      real*8 time_old,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old
      real*8 fxp,tk_pwp,dsea,tk_pwp2
	real sun_theta,dt_warm_profile(100)
      integer ID,jump,jamset,ihumid
      COMMON /old/time_old,qcol_ac,tau_ac,tau_old,rf_old,hf_old,
     &            ef_old,jamset,jump,fxp,tk_pwp    
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      COMMON/const/al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,
     &             visa,visw,von,fdg
      COMMON/wgust/DU_Wg,Wg    


C
C *
c
c initialize variables that appear in COMMON
c
      Jcool=Jcoolx
      Jwarm=Jwarmx
      ZU=hUm       !height of wind measurement
      ZT=hTm       !height of temperature measurement
      ZQ=hTm       !height of water vapor measurement
      ZUs=hUs      !standard height of wind measurement
      ZTs=hTs      !standard height of temperature measurement
      ZQs=hTs      !standard height of water vapor measurement
      U=ws         !wind speed m/s
	rws=ws
      TS=sst       !surface temp. Celsius
      T=atb        !air temp. Celsius
      P=pp         !pressure mb
      zi=zix       !atmospheric boundary layer depth
      toK=273.16
c	CLG add diurnal warming to sst for outgoing radiation
c 	R.Weil correction, use dt_wrmx rather than dt_wrm because it is previous time step values
      Rnl= 0.97*(5.67e-8*(TS+toK+dt_wrmx)**4-rl)    !Net longwave (up = +)		  CLG equation 10
      Rns=0.945*rs                          !Net shortwave (into water)		  CLG equation 11
      rain=rainx                            !rainfall
      ts_depth=ts_depthx                    !depth of sst measurement
	dt_wrm=0

c artifically increase winds used to cal profile in morning s

	r1=0;r2=0;r3=0;r4=0;r5=0;r6=0;						   !clg equation 17
	if(u.le.1.5)r2=1;
	if(u.gt.1.5 .and. u.le.3) then
	r2=1.-(u-1.5)/1.5;r3=1.-r2
	endif
	if(u.gt.3 .and. u.le.4.5) then
	r3=1.-(u-3)/1.5;r4=1.-r3
	endif
	if(u.gt.4.5 .and. u.le.6) then
	r4=1.-(u-4.5)/1.5;r5=1.-r4
	endif
	if(u.gt.6 .and. u.le.7.5) then
	r5=1.-(u-6.)/1.5;r6=1.-r5
	endif
	if(u.gt.7 ) then
	r6=1.
	endif

      call gravity(glat,grav)     
c
c Warm Layer
c
	dt_warm_profile=0;
      if(Jwarm.ne.0) then    
        newtime=sol_time                    !run time in secs 
        if(Jwarm.eq.2) then                    
          jump=1
          goto 16                                 !set time_old and pass thru' ASL
c       elseif(newtime.gt.21600.and.jump.eq.1) then
c        goto 16                                 !6 am too late to start
        elseif(u/=ws)then         !reset all var. at midnight
        else
		jump=0
          rich=.65                                    !critical Rich. No.
          ctd1=sqrt(2*rich*cpw/(al*grav*rhow))        !u*^2 integrated so
          ctd2=sqrt(2*al*grav/(rich*rhow))/(cpw**1.5) !has /rhow in both
          dtime=newtime-time_old                      !delta time
		if(dtime<0)dtime=dtime+86400
          qr_out=rnl+hf_old+ef_old+rf                 !flux out from previous pass
          q_pwp=fxp*rns-qr_out                        !effective net warming           CLG Equation 9
c          print*, fxp*rns,qr_out
c          print*, rnl,hr_old,ef_old,rf
c          if(q_pwp.lt.50.and.jamset.eq.0) print*, 'intg', q_pwp,jamset    !integration threshold
          if(q_pwp.lt.50.and.jamset.eq.0) go to 16    !integration threshold
          jamset=1
c	dissipate wind
c	function of wind stress for as square of wind gets higher, rate of dissipation increases
c	increase in turbulence rise to more viscous dissipations, look for papers on google scholar for this
c	wind cubed!because rate of deeping mixed layer wind cubed
c	rates .0009 m2/s3 clayson
		rtem_tau=(1.-rwnd_diss*dtime) !			  CLG equation 16
		if(rtem_tau<0)rtem_tau=0;
c	print*, dtime

          tau_ac=rtem_tau*tau_ac+tau_old*dtime  !tau from previous pass  maybe even have coeff of dissipation dep on wndspeed	CLG Equation 7
          if(qcol_ac+q_pwp*dtime.gt.0) then
            do 10 index=1,5                           !iterate for warm layer thicknes
c			tk_pwp2=tk_pwp*cosd(sun_theta)
			tk_pwp2=tk_pwp*cosd(sun_theta) !corrected by r.weihs 1.2011
              fxp=1.-(0.2370*dexp(-tk_pwp2/34.840)+ 0.3600*dexp(-tk_pwp2/2.2660)		  !nine-band model		CLG equation 14
	1			  + 0.1790*dexp(-tk_pwp2/0.0315)+ 0.0870*dexp(-tk_pwp2/0.0055)
	2			  + 0.0800*dexp(-tk_pwp2/8.32E-4)+ 0.0250*dexp(-tk_pwp2/1.26E-4)
	3			  + 0.0250*dexp(-tk_pwp2/3.13E-4)+ 0.0070*dexp(-tk_pwp2/7.82E-4)
	4			  + 0.0004*dexp(-tk_pwp2/1.44E-5))										 
			fxp=fxp*1.2																	!CLG equation 15
			if(fxp>1)fxp=1			!don't allow *1.2 to create absorption >1
              qjoule=(fxp*rns-qr_out)*dtime
              if((qcol_ac+qjoule.gt.0.0))
     &         tk_pwp=min(19.,ctd1*tau_ac/sqrt(qcol_ac+qjoule))			!CLG Equation 6
   10       continue
          else
            fxp=.76
            tk_pwp=19
            qjoule=(fxp*rns-qr_out)*dtime  
          endif

		rtem_qcol=(1.-rsol_diss*dtime)					  !CLG equation 16
		if(rtem_qcol<0)rtem_qcol=0;
          qcol_ac=rtem_qcol*qcol_ac+qjoule               !integrate heat input	CLG equation 5

          if(qcol_ac.gt.0) then
c			if(tau_ac<.0001)print*, 'tau',tau_ac
			if(tau_ac<.0001)tau_ac=.01
			dt_wrm=ctd2*(qcol_ac)**1.5/tau_ac  !pwp model warming   CLG Equation (4) 
c CLG R.W add test to ensure the dt_wrm does not explode because of ASL subroutine problems
			if(dt_wrm>10)dt_wrm=10
          else
            dt_wrm=0.
		  dt_warm_profile=0;
          endif      
		endif
		if(dt_wrm>0) then
			rsm1=0;rsm2=0;
			do i=1,100
				rz=i/100.
				rsm1=rsm1+dt_wrm*fun1(rz)
				rsm2=rsm2+r2*fun2(rz)+r3*fun3(rz)+r4*fun4(rz)+r5*fun5(rz)+r6*fun6(rz)
			enddo
			rsm1=rsm1/100.;rsm2=rsm2/100.;
			rwt=rsm1/rsm2
c return 100 observations of dt_warm within the warm layer (depth scaled by tk_wp)
			do i=1,100
				rz=i/100.
				dt_warm_profile(i)=rwt*(r2*fun2(rz)+r3*fun3(rz)+r4*fun4(rz)+r5*fun5(rz)+r6*fun6(rz))
			enddo
			dt_wrm=dt_warm_profile(1)
		endif
		if(tk_pwp.lt.ts_depth) then            !sensor deeper than pwp layer
		dsea=dt_wrm                          !all warming must be added to ts
		else                                   !warming deeper than sensor
		rz=ts_depth/tk_pwp
		iz=nint(rz*100)
		dt_wrm_sensor=dt_warm_profile(iz)
		dsea=dt_wrm-dt_wrm_sensor
		endif
        ts=ts+dsea                             !add warming above sensor for new ts
   16   time_old=newtime
        endif

c end of warm layer
c
   15 call humidity(T,P,QA)         !Teten's formula returns sat. air in mb
      if(qq.lt.2.) then             !checks whether humidity in g/Kg or RH      
         R=qq
         ee=QA*R                    !convert from RH using vapour pressure      
         Q=.62197*(ee/(P-0.378*ee)) !Spec. humidity kg/kg
      else
         Q=qq/1000.                 !g/kg to kg/kg
      endif
      QA=.62197*(QA/(P-0.378*QA))   !convert from mb to spec. humidity  kg/kg
      call humidity(TS,P,QS)        !sea QS returned in mb      
      QS=QS*0.98                    !reduced for salinity Kraus 1972 p. 46
      QS=.62197*(QS/(P-0.378*QS))   !convert from mb to spec. humidity  kg/kg
c
c calculate atmospheric surface layer      
c
      call ASL(Jcool,IER)
c	print*, 'ier >=0',ier
      if(IER.ge.0) then
c
c compute surface stress (TAU), sensible heat flux (HF),  
c latent heat flux (EF) & other parameters
c
      S=sqrt(u*u + wg*wg)           !velocity incl. gustiness param.
      TAU=rhoa*USR*usr*u/S          !kinematic units					   CLG equation 8
      HF=-cpa*rhoa*USR*TSR											   !CLG equation 12
      EF=-xlv*rhoa*USR*QSR											   !CLG equation 13
c
c compute heat flux due to rainfall
c
      dwat=2.11e-5*((T+toK)/toK)**1.94                    !water vapour diffusivity
      dtmp=(1.+3.309e-3*T-1.44e-6*T*T)*0.02411/(rhoa*cpa) !heat diffusivity
      dqs_dt=QA*xlv/(rgas*(T+toK)**2)                     !Clausius-Clapeyron
      alfac= 1/(1+0.622*(dqs_dt*xlv*dwat)/(cpa*dtmp))     !wet bulb factor
      W=-1.61*USR*QSR-(1+1.61*Q)*USR*TSR/(T+toK)
             
      rF= rain*alfac*cpw*((TS-T)+(QS-Q)*xlv/cpa)/3600.
c
c Webb correction to latent heat flux 
c
      Hl_webb=rhoa*xlv*W*Q
c
c compute momentum flux due to rainfall
c      
      TAU_r=0.85*rain/3600*u  
c
      dterx=dter               !cool skin parameters
      tktx=tkt          
      T0=ts-dter
      tau_old=tau 
c      ef_old=ef-Hl_webb
	ef_old=ef
      hf_old=hf
      dt_wrmx=dt_wrm           !warm layer parameter
c      
c compute transfer coefficients
c
         CD=(USR/S)**2
         CH=USR*TSR/(S*(T-TS+.0098*zt+dter)) !revise 2e to 2f to include '+dter'
         CE=USR*QSR/(S*(Q-QS+dqer))                                      
         Ustar=USR
         Tstar=TSR
         Qstar=QSR
         RRx=RR
         RTx=RT
         RQx=RQ
         ZLx=ZL
         ZOx=ZO
         zotx=zot
         zoqx=zoq
         wgx=wg
         ihumid=0
         if(qq .lt. 2) ihumid=1
c	r.weihs added check, no need to adjust heights if already at 10m heights
		IF (ZU .NE. ZUs .OR. ZT .NE. ZTs .OR. ZQ .NE. ZQs) then
			  call h_adjust(ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,ihumid)
			ws_h=U_hs
			Ta_h=T_hs
			qq_h=Q_hs 
		ELSE
			ws_h = ws
			Ta_h = atb
			qq_h = qq
		ENDIF
      else                           !input parameters out of range
         EF=-999.
         HF=-999.
         TAU=-999.
         TAUr=-999.
         EF_webb=-999.
         Ustar=-999.
         Tstar=-999.
         Qstar=-999.
         RRx=-999.
         RTx=-999.
         RQx=-999.
         ZLx=-999.
         ZOx=-999.
         ws_h=-999.
         Ta_h=-999.
         qq_h=-999.   
         wg=-999.
      endif
      return      !return to main program
      end
c
c ------------------------------------------------------------------

       
	real function fun1(z)
	real z
	fun1=(1.-z)
	end
	real function fun2(z)
	real z
c	fun2=exp(-15*z**2)
	fun2=exp(-9.5*z**2)
	end
	real function fun3(z)
	real z
	fun3=exp(-9.5*z**3)
	end
	real function fun4(z)
	real z
	fun4=exp(-9.5*z**5)
	end
	real function fun5(z)
	real z
	fun5=exp(-9.5*z**7)
	end
	real function fun6(z)
	real z
	fun6=exp(-9.5*z**9)
	end

