      subroutine fd_insolation_inst(lyear,idayjl,isecdy,xlat,xlon, hr_dawn, hr_dylen, rsol)
      implicit none

c Now, during the time period of 1978 through 1998, 
c the mean value of daily averages for the solar constant from six different satellites 
c yielded a solar constant of 1366.1 Wm2. This same source offers a minimum - maximum range of 
c the readings for 1363 - 1368 Wm2. Adjustments yielded "a solar constant" calculated to 1366.22 Wm2. 
c [Royal Meteorological Institute of Belgium: Department of Aerology; 
c http://remotesensing.oma.be/RadiometryPapers/article2.html.]. 
c	c.gentemann gentemann@remss.com 6/2009

	real(4), parameter :: solar_constant=1366.2   !watts/m**2
      real(4), parameter :: pi=3.141592654
c     inputs and outputs
	integer(4) lyear,idayjl,isecdy  ! year, julian day, seconds in day
	real(4)    xlat                 ! latitude of observation
	real(4)    rsol                 ! average insolation watts/m**2   
c	real(4) sunlat,sunlon,sundis,rr,H,sinh,Q,cosZ,hh,xlon,hr_dylen,hr_dawn,local_time,hr_dusk
	real(4) sunlat,sunlon,sundis,rr,H,cosZ,hh,xlon,hr_dylen,hr_dawn,local_time,hr_dusk

	call sunloc1(lyear,idayjl,isecdy, sunlat,sunlon,sundis)	  !sundis is radius/(mean radius)

      rr=-tand(xlat)*tand(sunlat);
      if(rr> 1) rr= 1.
      if(rr<-1) rr=-1.
      H = acos(rr);
	hr_dylen=H*24./pi
	hr_dawn=12.-hr_dylen/2.
	hr_dusk=12.+hr_dylen/2.
	local_time=(isecdy/86400.*24.+xlon/360.*24.)
	if(local_time>24)local_time=local_time-24.
	hh=local_time
	if((hh>hr_dawn .and. hh<hr_dusk) .or. hr_dylen==24) then !inbetween dawn and dusk
	cosZ=sind(xlat)*sind(sunlat) + cosd(xlat)*cosd(sunlat)*cos((hh+12)/24.*2*pi)
	rsol=solar_constant*cosZ/sundis**2	
	else
	rsol=0
	endif
	end

      subroutine fd_insolation_peak(lyear,idayjl,xlat,xlon, rsol)
      implicit none

c Now, during the time period of 1978 through 1998, 
c the mean value of daily averages for the solar constant from six different satellites 
c yielded a solar constant of 1366.1 Wm2. This same source offers a minimum - maximum range of 
c the readings for 1363 - 1368 Wm2. Adjustments yielded "a solar constant" calculated to 1366.22 Wm2. 
c [Royal Meteorological Institute of Belgium: Department of Aerology; 
c http://remotesensing.oma.be/RadiometryPapers/article2.html.]. 

	real(4), parameter :: solar_constant=1366.2   !watts/m**2
      real(4), parameter :: pi=3.141592654
c     inputs and outputs
	integer(4) lyear,idayjl  ! year, julian day
	real(4)    xlat                 ! latitude of observation
	real(4)    rsol                 ! average insolation watts/m**2   
	real(4) sunlat,sunlon,sundis,rr,H,cosZ,hh,xlon,hr_dylen,hr_dawn,local_time,hr_dusk
c	real(4) sunlat,sunlon,sundis,rr,H,sinh,Q,cosZ,hh,xlon,hr_dylen,hr_dawn,local_time,hr_dusk

	call sunloc1(lyear,idayjl,43200, sunlat,sunlon,sundis)	  !sundis is radius/(mean radius)
      rr=-tand(xlat)*tand(sunlat);
      if(rr> 1) rr= 1.
      if(rr<-1) rr=-1.
      H = acos(rr);
	hr_dylen=H*24./pi
	hr_dawn=12.-hr_dylen/2.
	hr_dusk=12.+hr_dylen/2.
	local_time=(43200/86400.*24.+xlon/360.*24.)
	if(local_time>24)local_time=local_time-24.
	hh=local_time
	if((hh>hr_dawn .and. hh<hr_dusk) .or. hr_dylen==24) then !inbetween dawn and dusk
	cosZ=sind(xlat)*sind(sunlat) + cosd(xlat)*cosd(sunlat)*cos((hh+12)/24.*2*pi)
	rsol=solar_constant*cosZ/sundis**2	
	else
	rsol=0
	endif
      return
      end

c      INCLUDE 'O:\sun_moon\SUNLOC1.F'

