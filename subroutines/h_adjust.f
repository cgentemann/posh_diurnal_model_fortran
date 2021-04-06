c-------------------------------------------------------------------------
      subroutine H_ADJUST(ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,IHUMID)
c
c This subroutine adjusts the U,T,Q variables to the specified
c standard height (ZUs,ZTs,ZQs) using the loglayer profiles.
c The DELTA correction (adjustment) is relative to the surface
c measurement.             Cronin 4/13/94
c     
      real*8 ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,ZUsL,ZTsL,ZQsL
      real*8 U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,PUZs,PTZs,PQZs
      real*8 U_wg_hs,Rho_hs,Rho_avg,QA,Rho,P_hs,ee
      real*8 USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      real*8 al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,visa
      real*8 visw,von,fdg,DU_Wg,Wg,S,D
      integer ID,ihumid 
c      
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      COMMON/const/al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,
     &             visa,visw,von,fdg
      COMMON/wgust/DU_Wg,Wg

      call ZETA(T,Q,USR,TSR,QSR,ZUs,ZUsL)
      call ZETA(T,Q,USR,TSR,QSR,ZTs,ZTsL)
      call ZETA(T,Q,USR,TSR,QSR,ZQs,ZQsL)
      PUZs= PSI(1,ZUsL)
      PTZs= PSI(2,ZTsL)
      PQZs= PSI(2,ZQsL)
 
      S = (dlog(ZTs*USR/VISA/RT)-PTZs)/(von*fdg)
      D = (dlog(ZQs*USR/VISA/RQ)-PQZs)/(von*fdg)
      T_hs =TSR*S +TS - dter -.0098*ZTs
      Q_hs =(QSR*D + QS - dqer)*1000
      U_wg_hs = USR*(dlog(ZUs/ZO) - PUZs)/0.4
      if(U_wg_hs.ge.Wg) then
         U_hs = SQRT(U_wg_hs**2 - Wg**2)
      else
         U_hs = U_wg_hs
      endif
c
      if(IHUMID.eq.1) then    ! then need to convert sp hum into rh
         Q_hs = Q_hs/1000     ! sh kg/kg
         RHO=1./(287.*(T+273.16)*(1.+.61*Q))*P*100.
         P_hs = P - (RHO*grav*(ZTs - ZT))/100 !Approx hydrost.Pressure mb
         RHO_hs=1./(287.*(T_hs+273.16)*(1.+.61*Q_hs))*P_hs*100
         RHO_avg = (RHO + RHO_hs)/2
         P_hs = P -(RHO_avg*grav*(ZTs - ZT))/100 !hydrostatic Pressure
         call humidity(T_hs,P_hs,QA)         !Teten's formula for Pvap,sat
         ee=Q_hs*P_hs/(.62197 + .378*Q_hs)   !to get vapor pressure
         Q_hs = ee/QA                        !to get relative humidity
      endif
      return
      end
