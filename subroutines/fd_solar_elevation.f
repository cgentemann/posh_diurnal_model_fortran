      subroutine fd_solar_elevation(lyear,idayjl,isecdy,xlat,xlon, sun_theta)
      implicit none

	real(4), parameter :: solar_constant=1366.2   !watts/m**2
      real(4), parameter :: pi=3.141592654
c     inputs and outputs
	integer(4) lyear,idayjl,isecdy  ! year, julian day, seconds in day
	real(4)    xlat                 ! latitude of observation
	real(4)    sun_theta                 ! sun elevation  
	real(4) sunlat,sunlon,sundis,hh,xlon,local_time,cosZ

	call sunloc1(lyear,idayjl,isecdy, sunlat,sunlon,sundis)	  !sundis is radius/(mean radius)
	local_time=(isecdy/86400.*24.+xlon/360.*24.)
	if(local_time>24)local_time=local_time-24.
	hh=local_time
	cosZ=sind(xlat)*sind(sunlat) + cosd(xlat)*cosd(sunlat)*cos((hh+12)/24.*2*pi)
	sun_theta=asin(cosZ)*180./pi
	if(sun_theta<-90)sun_theta=-90
	if(sun_theta>90)sun_theta=90
      return
      end
