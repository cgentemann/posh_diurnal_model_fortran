c------------------------------------------------------------
      subroutine ZETA(T,Q,USR,TSR,QSR,Z,ZL)
C
C TO EVALUATE OBUKHOVS STABILITY PARAMETER Z/L FROM AVERAGE
C TEMP T IN DEG C, AVERAGE HUMIDITY Q IN GM/GM, HEIGHT IN M,
C AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS
C SEE LIU ET AL. (1979)
C     
      real*8 T,Q,OB,TVSR,TV,TA,sgn
      real*8 USR,TSR,QSR,Z,ZL
      real*8 al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,visa
      real*8 visw,von,fdg
      COMMON/const/al,beta,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,
     &             visa,visw,von,fdg
c
      TA=T+toK
      TV=TA*(1.+0.61*Q)
      TVSR=TSR*(1.+0.61*Q)+0.61*TA*QSR    
      sgn=sign(1.,tvsr)               !added this to avoid program
      if(abs(tvsr) .lt. 1.e-3) then   !failure when TVSR is very small
         tvsr=sgn*tvsr
      endif
      OB=TV*USR*USR/(grav*VON*TVSR) 
      ZL=Z/OB                     
      if(ZL .gt. 1000) ZL=1000.
      goto 99
   10 ZL=0. 
   99 return
      end
