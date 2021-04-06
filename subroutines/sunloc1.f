C     12/12/2003 VERSION CHANGED ON AUG 26 2005.  IMPLICIT NONE WAS ADDED.
C     Same as SUNLOC August 01, 1998, 01:04:31 PM, except this also returns R.
	SUBROUTINE SUNLOC1(LYEAR,IDAYJL,ISECDY, SUNLAT,SUNLON,R)
	IMPLICIT NONE

	REAL(8) XDAY,TX
      REAL(4) Z(3),X(3),TOSUN(3),PROJECT(3)
	real(8) D01,D11,S11,C11,S21,C21
	real(8) D02,D12,S12,C12,S22,C22
	real(8) slope,period,ecc,xinc,rad,xfac,alpha,t0,t,ang,cosinc
	real(8) amag,xphi,eot,hang,sinang,cosang,sinhang,coshang,corrlat,corrlon
	real(4) phi
	integer(4) istart,lyear,idayjl,isecdy
	real(4) sunlat,sunlon,r

	DATA D01,D11,S11,C11,S21,C21/ -1.15538E-3, -6.86302E-8, -8.48082E-1, -1.29184E-1, -3.18681E-2, -5.59290E-4/
	DATA D02,D12,S12,C12,S22,C22/ -1.57115E-3, -3.82973E-5, -9.91281E-4, -2.75336E-2,  1.83422E-2,  1.09567E-3/
	DATA SLOPE/1.81818E-5/ !=0.2/11000.0
	DATA PERIOD,ECC,XINC/365.2422,0.016708617,23.43929111/
      DATA RAD/0.01745329252/
	DATA ISTART/1/

	IF(ISTART.EQ.1) THEN
      XFAC=360./PERIOD                 
      ALPHA=XFAC*(355.5674 - 0.3210) !WINTER SOLSTICES DEC 21 13:37, 2000; ADDITION OFFSET OF .321 TO COMPENSTATE FOR ECC
	T0=XFAC*2.2083 !PERIHELION JAN 3, 5 HOURS, 2000
	Z(1)=SIND(XINC)*COSD(ALPHA)
	Z(2)=SIND(XINC)*SIND(ALPHA)
	Z(3)=COSD(XINC)
	X(1)=COSD(XINC)*COSD(ALPHA+180)
	X(2)=COSD(XINC)*SIND(ALPHA+180)
	X(3)=SIND(XINC)
	ENDIF

	IF(LYEAR.LT.1950 .OR. LYEAR.GT.2099) STOP 'ERROR IN LYEAR'
	XDAY=365.D0*(LYEAR-1950) + (IDAYJL-1) + ISECDY/86400.D0   + INT((LYEAR-1949)/4) - 18262.D0   !TDAY=0 FOR JAN 1 0Z, 2000 
      TX=XFAC*XDAY
      T=MOD(TX,360.D0)
      IF(T.LT.0) T=T+360. 
                                                 
      ANG=T+2.*ECC*(SIND(T-T0)+SIND(T0))/RAD                                           
      R=1.-ECC*COSD(T-T0)

      TOSUN(1)=-COSD(ANG)                                                          
      TOSUN(2)=-SIND(ANG)                                                  
      TOSUN(3)=0
      COSINC=DOT_PRODUCT(TOSUN,Z) 
      SUNLAT=90-ACOSD(COSINC) 

      PROJECT=TOSUN - COSINC*Z 
	AMAG=SQRT(PROJECT(1)**2 + PROJECT(2)**2 + PROJECT(3)**2)
	PROJECT=PROJECT/AMAG
      CALL FINDANG(X,PROJECT,Z, PHI)                             
      XPHI=PHI-9.6 -0.326
	EOT=T-XPHI
	IF(EOT.LT.-180) EOT=EOT+360
	IF(EOT.GT. 180) EOT=EOT-360

	SUNLON=180-360*ISECDY/86400. - EOT
	IF(SUNLON.LT.0) SUNLON=SUNLON+360.

	ANG=XFAC*XDAY
	HANG=0.5*ANG
	SINANG=SIND(ANG)
	COSANG=COSD(ANG)
	SINHANG=SIND(HANG)
	COSHANG=COSD(HANG)
	CORRLAT = D01 + D11*XDAY + SLOPE*XDAY*(S11*SINANG + C11*COSANG + S21*SINHANG + C21*COSHANG)
	CORRLON = D02 + D12*XDAY +             S12*SINANG + C12*COSANG + S22*SINHANG + C22*COSHANG

	SUNLAT=SUNLAT + CORRLAT
	SUNLON=SUNLON + CORRLON
	IF(SUNLON.LT.  0) SUNLON=SUNLON+360
	IF(SUNLON.GE.360) SUNLON=SUNLON-360

      RETURN
      END 


                                                           
      SUBROUTINE FINDANG(X,Y,Z, ANG) 
	implicit none       
       
	REAL(4) X(3),Y(3),Z(3),D(3)
	real(4) sinang,cosang,ang

      D(1)=X(2)*Y(3)-X(3)*Y(2)                                                    
      D(2)=X(3)*Y(1)-X(1)*Y(3)                                                    
      D(3)=X(1)*Y(2)-X(2)*Y(1)                                                    
      SINANG=SQRT(D(1)**2+D(2)**2+D(3)**2) 
	IF(DOT_PRODUCT(D,Z).LT.0) SINANG=-SINANG 
      COSANG=DOT_PRODUCT(X,Y) 
	ANG=ATAN2D(SINANG,COSANG)
	IF(ANG.LT.0) ANG=ANG+360.
      RETURN                                                            
      END                                                           


