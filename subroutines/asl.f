      subroutine ASL(Jcool,IER)
c
c TO EVALUATE SURFACE FLUXES, SURFACE ROUGHNESS AND STABILITY OF
c THE ATMOSPHERIC SURFACE LAYER FROM BULK PARAMETERS BASED ON
c LIU ET AL. (79) JAS 36 1722-1735 
c    
      real*8 U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,D,S
      real*8 USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      real*8 al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,visa
      real*8 visw,von,fdg,DU_Wg,Wg
      real*8 PTZ,PQZ,PUZ,ZTl,ZLN,ZQl,wetc,bigc,cpv,dq,dqq,dt,dtt
      integer ID
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      COMMON/const/al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,
     &             visa,visw,von,fdg
      COMMON/wgust/DU_Wg,Wg     

c Factors
      Beta=1.2     ! evaluated from Fairall's low windspeed turbulence data
      Von=0.4      ! von Karman's "constant"
c      fdg=1.00     ! Fairall's LKB rr to von karman adjustment
      fdg=1.00     !based on results from Flux workshop August 1995
      toK=273.16   ! Celsius to Kelvin
c      grav=9.72    ! gravity equatorial value (ref. IGPP-SIO)
 
c Air constants and coefficients
      Rgas=287.1                  !J/kg/K     gas const. dry air
      xlv=(2.501-0.00237*TS)*1e+6  !J/kg  latent heat of vaporization at TS
      Cpa=1004.67                 !J/kg/K specific heat of dry air (Businger 1982)
      Cpv=Cpa*(1+0.84*Q)          !Moist air - currently not used (Businger 1982)
      rhoa=P*100./(Rgas*(T+toK)*(1.+.61*Q)) !kg/m3  Moist air density ( " )
      visa=1.326e-5*(1+6.542e-3*T+8.301e-6*T*T-4.84e-9*T*T*T)   !m2/s
          !Kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
 
c Cool skin constants
      al=2.1e-5*(ts+3.2)**0.79     !water thermal expansion coefft.
      be=0.026                     !salinity expansion coefft.
      cpw=4000.                    !J/kg/K specific heat water
      rhow=1022.                   !kg/m3  density water
      visw=1.e-6                   !m2/s kinematic viscosity water
      tcw=0.6                      !W/m/K   Thermal conductivity water
      bigc=16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
      wetc=.622*xlv*QS/(rgas*(TS+toK)**2) !correction for dq;slope of sat. vap.
 
c Initialise everything
      IER=0
      ZL=0.                         
      Dter=0.                        !cool skin Dt
      Dqer=0.                        !cool skin Dq

c Initial guesses
c     next line originally required for surface current, input now requires relative wind
c      US=0.                        !surface current = 0.
c
      Wg=0.5                        !Gustiness factor initial guess
      ZO=.0005                      !roughness initial guess
      tkt=.001                      !guess sublayer thickness
      DU=U                          !assumes U is measured rel. to current
      DU_Wg=(DU**2.+Wg**2.)**.5     !include gustiness in wind spd. difference
                                    !equivalent to S in definition of fluxes
      DT=T-TS+.0098*zt              !potential temperature diff        
      DQ=Q-QS
      USR=.04*DU_Wg
      TSR=.04*DT 
      QSR=.04*DQ
      if(DU_Wg.ne.0.)then
         TA=T+toK
c         print*, grav,ZU,dt,ta,dq,du_wg,grav*zu*(DT+0.61*TA*DQ),(TA*DU_Wg**2)
         RI=grav*ZU*(DT+0.61*TA*DQ)/(TA*DU_Wg**2)
      else
c		print*, 'du_wg/=0',du_Wg
         IER=-1
         RI=-999.
      endif
c     if(RI.gt.0.25) print*, 'RI > .25'
      if(RI.gt.0.25) IER=-1
      do 200 index=1,20        !iterate 20 times 
         call ZETA(T,Q,USR,TSR,QSR,ZU,ZLN)
         ZL=ZLN
         PUZ=PSI(1,ZL)
         ZTL=ZL*ZT/ZU
         ZQL=ZL*ZQ/ZU
         PTZ=PSI(2,ZTL)
         PQZ=PSI(2,ZQL)
         ZO=0.011*USR*USR/grav + 0.11*visa/USR    !after Smith 1988 
         USR=DU_Wg*von/(dlog(ZU/ZO)-PUZ)          !Gustiness effect incl.
         RR=ZO*USR/VISA 
         call LKB(RR,RT,1)      
         if(RT.eq.-999.) then
           IER = -2                         !error - return
          return
         endif  
   21    call LKB(RR,RQ,2)
         if(RQ.eq.-999.) then
           IER = -2                         !error - return
           return
         endif
   22    zot=rt*visa/usr
         zoq=rq*visa/usr
         S = (dlog(ZT/zot)-PTZ)/(von*fdg) !coeff fdg=1.04 included following
         D = (dlog(ZQ/zoq)-PQZ)/(von*fdg) !Fairall observations during COARE. 
                                          !NOTE coeff changed to 1.
                                          !This has been adjusted to 0.94
                                          !following discussion during FLux
                                          !workshop August 1995.
         dtt=(dt+dter)
         dqq=(dq+dqer)               ! or we could calculate new sat. hum.
         tsr=dtt/S                              !! modify
         qsr=dqq/D                              !! fluxes
         TVSR=TSR*(1.+0.61*Q)+(0.61*TA*QSR)
         Bf=-grav/TA*USR*TVSR
         if(Bf.gt.0) then
           Wg=Beta*(Bf*zi)**0.333  
         else
           Wg=0.
         endif
         DU_Wg=(DU**2.+Wg**2.)**.5        !include gustiness in wind spd.
         if(Jcool.ne.0) then              !Cool skin
           hsb=-rhoa*cpa*usr*tsr
           hlb=-rhoa*xlv*usr*qsr
c             write(*,'(A,5F16.8)')  'hlb',-rhoa,xlv,usr,qsr
           qout=rnl+hsb+hlb
           dels=rns*(.137+11*tkt-6.6e-5/tkt*(1-dexp(-tkt/8.0e-4))) ! Eq.16 
                                                    ! Shortwave, correction
                                                    ! from version 2e to 2f
           qcol=qout-dels
           if(qcol.gt.0.) then
             alq=Al*qcol+be*hlb*cpw/xlv                ! Eq. 7 Buoy flux water
c             write(*,'(A,5F16.8)')  ' alq',be,hlb,cpw,xlv
             xlamx=6/(1+abs(bigc*alq/usr**4)**.75)**.333  ! Eq 13 Saunders coeff.
             tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)      ! Eq.11 Sublayer thickness
             dter=qcol*tkt/tcw                         ! Eq.12 Cool skin
           else
             dter=0.
           endif    
           dqer=wetc*dter
         endif             !end cool skin
  200 continue             !end iterations                 
      return               !to main subroutine, bulk_flux
      end
c
c------------------------------------------------------------------
      subroutine humidity(T,P,Qsat)                                 
c
c Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
c     
      real*8 T,P,Qsat
c     
      Qsat = (1.0007+3.46e-6*P)*6.1121*dexp(17.502*T/(240.97+T))
      return
      end
c
c------------------------------------------------------------------
      subroutine LKB(RR,RT,IFLAG)
c
c TO DETERMINE THE LOWER BOUNDARY VALUE RT OF THE 
c LOGARITHMIC PROFILES OF TEMPERATURE (IFLAG=1) 
c OR HUMIDITY (IFLAG=2) IN THE ATMOSPHERE FROM ROUGHNESS 
c REYNOLD NUMBER RR BETWEEN 0 AND 1000.  OUT OF RANGE
c RR INDICATED BY RT=-999. BASED ON LIU ET AL.(1979)
c JAS 36 1722-1723 
c-------------------------------------------------------------------
c Scalar RR relation from Moana Wave data.
c
c      real*8 A(9,2),B(9,2),RAN(9),RR,RT
c      integer iflag
c      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
c     &       0.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
c      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
c     &       0.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
c      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./        
c-------------------------------------------------------------------
c
c Scalar RR relation from Liu et al. 
c
      real*8 a(8,2),b(8,2),ran(8),RR,RT
      integer iflag                                            
c
      data a/0.177,1.376,1.026,1.625,4.661,34.904,1667.19,5.88e5,               
     & 0.292,1.808,1.393,1.956,4.994,30.709,1448.68,2.98e5/                     
      data b/0.,0.929,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,                
     & 0.,0.826,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/                      
      data ran/0.11,0.825,3.0,10.0,30.0,100.,300.,1000./                        
c    
      I=1
      if((RR.le.0.).or.(RR.ge.1000.)) goto 90
   10 continue
      if(RR.le.RAN(I)) goto 20
      I=I+1
      goto 10
   20 RT=A(I,IFLAG)*RR**B(I,IFLAG)
      goto 99
   90 RT=-999. 
   99 return
      end
c
c-----------------------------------------------------------------
      function PSI(ID,ZL)
C
C TO EVALUATE THE STABILITY FUNCTION PSI FOR WIND SPEED (IFLAG=1)
C OR FOR TEMPERATURE AND HUMIDITY PROFILES FROM STABILITY 
C PARAMETER ZL. SEE LIU ET AL (1979).
c Modified to include convective form following Fairall (Unpublished)
C     
      real*8 Zl,chik,psik,f,chic,psic
      integer ID
c      
      if(ZL)10,20,30                                                 
   10 F=1./(1+zl*zl) 
      CHIK=(1.-16.*ZL)**0.25
      if(ID.eq.1) goto 11
      PSIK=2.*dlog((1.+CHIK*CHIK)/2.)
      goto 12
   11 PSIK=2.*dlog((1.+CHIK)/2.)+dlog((1.+CHIK*CHIK)/2.)
     & -2.*ATAN(CHIK)+2.*ATAN(1.)
   12 CHIC=(1.-12.87*ZL)**.333    !for very unstable conditions
      PSIC=1.5*dlog((CHIC*CHIC+CHIC+1.)/3.)
     &      -(3.**.5)*ATAN((2*CHIC+1.)/(3.**.5))
     &      +4.*ATAN(1.)/(3.**0.5) 
c                     
c match Kansas and free-conv. forms with weighting F
c
      PSI= F*PSIK+(1-F)*PSIC                                        
      goto 99
   20 PSI=0.
      goto 99
   30 continue
      PSI=-4.7*ZL
   99 return
      end
c
