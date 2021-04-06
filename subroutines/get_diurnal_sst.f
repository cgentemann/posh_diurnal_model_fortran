      subroutine get_diurnal_sst(xhour,xlon,wind, Q,diurnal_sst)
      implicit none

      real(4), parameter :: pi=3.141592654
	real(4), parameter :: solar_constant=1360   !watts/m**2

c	real(4) sunlat,sundis,xhour,xlon,xlat,wind
	real(4) xhour,wind,xlon
      real(4) diurnal_sst, Q, rlmt
      real(4) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,w	 


      data a0,a1,a2,a3,a4,a5/4.118E-3, -4.132E-3, 0.8746E-3, -0.2460E-3, 0.2762E-3,-0.0609E-3/
      data    b1,b2,b3,b4,b5/          -5.093E-3, 2.5830E-3, -0.5143E-3,-0.3355E-3, 0.2269E-3/
      data w/0.261799388/  !2*pi/24.

      rlmt=xhour + (xlon/360.)*24.

      if (Q>154) then
      diurnal_sst=(1.23459589668*(Q-154.)-1.56088414175378e-3*(Q-154.)**2)*exp(-.44*wind)*
     &				(a0 + a1*cos(   rlmt*w) + b1*sin(   rlmt*w) + a2*cos(2.*rlmt*w) + b2*sin(2.*rlmt*w)+
     &				      a3*cos(3.*rlmt*w) + b3*sin(3.*rlmt*w) + a4*cos(4.*rlmt*w) + b4*sin(4.*rlmt*w)+
     &                      a5*cos(5.*rlmt*w) + b5*sin(5.*rlmt*w))
	else
	diurnal_sst=0;
	endif
	if(diurnal_sst<0)diurnal_sst=0;

	RETURN
	END


	subroutine get_diurnal_sst_seviri2(xhour,xlon,wind, rsol,dsst)


c	c.gentemann 4.2008
c	inputs:
c	lyr=year
c	idyjl=ordinal day
c	xhour = LOCAL TIME IN HOURS
c	xlat=latitude in deg
c	wind = wind speed in m/s

c	output:
c	dsst = diurnal amplitude in deg

c	integer lyr,idyjl
	real xhour,wind,dsst
	real w,rsol
      real(4), parameter :: pi=3.141592654
	real*8 xpar(11)
	real diurnal_cos

	xpar(1)=0.90779440686153
	xpar(2)=-5.654018322117584E-001 
	xpar(3)=-9.109286437084861E-001
	xpar(4)=4.938677367379611E-002
	xpar(5)=2.921478004135862E-001  
	xpar(6)=-8.767453446309354E-004 
	xpar(7)= -6.803299659334103E-002  
	xpar(8)=-4.456366487917374E-003   
	xpar(9)=2.214046984740299E-003
	xpar(10)=2.166914823201959E-002  
	xpar(11)=-2.116565831620519E-003

c	a0=4.118;a1=-4.132;b1=-5.093;a2=.8746;b2=2.583;
c	a3=-.246;b3=-.5143;a4=.2762;b4=-.3355;a5=-.0609;b5=.2269;
	w=2.*pi/24.;

c	call fd_insolation(lyr,idyjl,43200,xlat,rsol)
      rlmt=xhour + (xlon/360.)*24.

	dsst=0;
	if (rsol>97) then
	diurnal_cos=xpar(1)+
	.			xpar(2)*cos(rlmt*w)+xpar(3)*sin(rlmt*w)+
     .			xpar(4)*cos(2*rlmt*w)+xpar(5)*sin(2*rlmt*w)+
     .			xpar(6)*cos(3*rlmt*w)+xpar(7)*sin(3*rlmt*w)+
     .			xpar(8)*cos(4*rlmt*w)+xpar(9)*sin(4*rlmt*w)+
     .			xpar(10)*cos(5*rlmt*w)+xpar(11)*sin(5*rlmt*w)


	dsst=((0.00224946*(rsol-97.))*exp(-0.27*wind) +
	+	 (2.01422e-6*(rsol-97.)**2)*exp(-0.27*wind))*diurnal_cos
     	endif


	end

