      subroutine fd_insolation(lyear,idayjl,isecdy,xlat, rsol)
      implicit none

c Now, during the time period of 1978 through 1998, 
c the mean value of daily averages for the solar constant from six different satellites 
c yielded a solar constant of 1366.1 Wm2. This same source offers a minimum - maximum range of 
c the readings for 1363 - 1368 Wm2. Adjustments yielded "a solar constant" calculated to 1366.22 Wm2. 
c [Royal Meteorological Institute of Belgium: Department of Aerology; 
c http://remotesensing.oma.be/RadiometryPapers/article2.html.]. 
c	c.gentemann gentemann@remss.com
c	6/2009

	real(4), parameter :: solar_constant=1366.2   !watts/m**2
      real(4), parameter :: pi=3.141592654

c     inputs and outputs
	integer(4) lyear,idayjl,isecdy  ! year, julian day, seconds in day
	real(4)    xlat                 ! latitude of observation
	real(4)    rsol                 ! average insolation watts/m**2   
	real(4) sunlat,sunlon,sundis,rr,H,sinh,Q

	call sunloc1(lyear,idayjl,isecdy, sunlat,sunlon,sundis)	  !sundis is radius/(mean radius)
        
      rr=-tand(xlat)*tand(sunlat);
      if(rr> 1) rr= 1.
      if(rr<-1) rr=-1.

      H = acos(rr);	 
      sinH=sqrt(1 - rr**2) !sqrt(1-cos(h)**2) =sin(h)
	Q=sind(xlat)*sind(sunlat)*H + cosd(xlat)*cosd(sunlat)*sinH
      rsol = (solar_constant/pi)*Q/sundis**2

      return
      end
