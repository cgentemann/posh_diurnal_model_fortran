	subroutine get_diurnal_sst_mw(lyr,idyjl,xhour,xlat,xlon,wind, dsst)

c	c.gentemann gentemann@remss.com 
c	c.gentemann 10/2004
c	inputs:
c	lyr=year
c	idyjl=ordinal day
c	xhour = LOCAL TIME IN HOURS
c	xlat=latitude in deg
c	wind = wind speed in m/s

c	output:
c	dsst = diurnal amplitude in deg

	integer lyr,idyjl
	real xhour,xlat,wind,dsst
	real a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,w,rsol
      real(4), parameter :: pi=3.141592654

	a0=4.118;a1=-4.132;b1=-5.093;a2=.8746;b2=2.583;
	a3=-.246;b3=-.5143;a4=.2762;b4=-.3355;a5=-.0609;b5=.2269;
	w=2*pi/24.;

	rhr=xhour+xlon/360.*24.
	if(rhr>24.) rhr=rhr-24
	if(rhr==0) rhr=24

	call fd_insolation(lyr,idyjl,43200,xlat,rsol)
	dsst=0;
	if (rsol>134) then
	dsst=(0.00102838361636*(rsol-134.)-9.423218330084637e-7*(rsol-134.)**2)*exp(-.44*wind)*
     +				(a0+a1*cos(rhr*w)+b1*sin(rhr*w)+a2*cos(2.*rhr*w)+b2*sin(2.*rhr*w)+
     +				 a3*cos(3.*rhr*w)+b3*sin(3.*rhr*w)+a4*cos(4.*rhr*w)+b4*sin(4.*rhr*w)+
     +			  	 a5*cos(5.*rhr*w)+b5*sin(5.*rhr*w))
	endif
	if(dsst<0)dsst=0;
      end	subroutine



	subroutine get_diurnal_sst_ir(lyr,idyjl,xhour,xlat,xlon,wind, dsst)

c	c.gentemann 10.2004
c	inputs:
c	lyr=year
c	idyjl=ordinal day
c	xhour = LOCAL TIME IN HOURS
c	xlat=latitude in deg
c	wind = wind speed in m/s

c	output:
c	dsst = diurnal amplitude in deg

	integer lyr,idyjl
	real xhour,xlat,wind,dsst
	real a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,w,rsol
      real(4), parameter :: pi=3.141592654

      a0 =       1.169/1000.*3.163003;
      a1 =     -0.6638/1000.*3.163003;
      b1 =        -1.2/1000.*3.163003;
      a2 =     0.03301/1000.*3.163003;
      b2 =      0.3939/1000.*3.163003;
      a3 =     0.02399/1000.*3.163003;
      b3 =    -0.04651/1000.*3.163003;
      a4 =     0.00446/1000.*3.163003;
      b4 =   -0.005024/1000.*3.163003;
      a5 =    0.008777/1000.*3.163003;
      b5 =   -0.008232/1000.*3.163003;
      
	w=2*pi/24.;
	rhr=xhour+xlon/360.*24.
	if(rhr>24.) rhr=rhr-24
	if(rhr==0) rhr=24

	call fd_insolation(lyr,idyjl,43200,xlat,rsol)

	dsst=0;
	if (rsol>88) then
	dsst=(0.9767161*(rsol-88.)-.0004813153*(rsol-88.)**2)*exp(-.28*wind)*
     +				(a0+a1*cos(rhr*w)+b1*sin(rhr*w)+a2*cos(2.*rhr*w)+b2*sin(2.*rhr*w)+
     +				 a3*cos(3.*rhr*w)+b3*sin(3.*rhr*w)+a4*cos(4.*rhr*w)+b4*sin(4.*rhr*w)+
     +				 a5*cos(5.*rhr*w)+b5*sin(5.*rhr*w))
	endif
	if(dsst<0)dsst=0;

      end	subroutine


	subroutine get_diurnal_sst_seviri2(lyr,idyjl,xhour,xlat,xlon,wind, dsst)

c	c.gentemann 4.2008
c	inputs:
c	lyr=year
c	idyjl=ordinal day
c	xhour = LOCAL TIME IN HOURS
c	xlat=latitude in deg
c	wind = wind speed in m/s

c	output:
c	dsst = diurnal amplitude in deg

	integer lyr,idyjl
	real xhour,xlat,wind,dsst
c	real a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5
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
	rhr=xhour+xlon/360.*24.
	if(rhr>24.) rhr=rhr-24
	if(rhr==0) rhr=24

	call fd_insolation(lyr,idyjl,43200,xlat,rsol)

	dsst=0;
	if (rsol>97) then

	diurnal_cos=xpar(1)+
     +			xpar(2)*cos(rhr*w)+xpar(3)*sin(rhr*w)+
     +			xpar(4)*cos(2*rhr*w)+xpar(5)*sin(2*rhr*w)+
     +			xpar(6)*cos(3*rhr*w)+xpar(7)*sin(3*rhr*w)+
     +			xpar(8)*cos(4*rhr*w)+xpar(9)*sin(4*rhr*w)+
     +			xpar(10)*cos(5*rhr*w)+xpar(11)*sin(5*rhr*w)


	dsst=((0.00224946*(rsol-97.))*exp(-0.27*wind)+(2.01422e-6*(rsol-97.)**2)*exp(-0.27*wind))*diurnal_cos
     	endif

	if(dsst<0)dsst=0;

	end


