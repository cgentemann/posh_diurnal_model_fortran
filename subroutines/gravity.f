      Subroutine gravity(lat,g)
c       calculates g as a funciton of latitude using the 1980 IUGG formula
c         
c       Bulletin Geodesique, Vol 62, No 3, 1988 (Geodesist's Handbook)
c       p 356, 1980 Gravity Formula (IUGG, H. Moritz)
c       units are in m/sec^2 and have a relative precision of 1 part
c       in 10^10 (0.1 microGal)
c       code by M. Zumberge.
c
c       check values are:
c
c        g = 9.780326772 at latitude  0.0
c        g = 9.806199203 at latitude 45.0
c        g = 9.832186368 at latitude 90.0
c
      real*8 gamma, c1, c2, c3, c4, phi, lat, g
 
      gamma = 9.7803267715
      c1 = 0.0052790414
      c2 = 0.0000232718
      c3 = 0.0000001262
      c4 = 0.0000000007
      phi = lat * 3.14159265358979 / 180.0
      g = gamma * (1.0 
     $ + c1 * ((sin(phi))**2)
     $ + c2 * ((sin(phi))**4)
     $ + c3 * ((sin(phi))**6)
     $ + c4 * ((sin(phi))**8))
c
      return
      end
