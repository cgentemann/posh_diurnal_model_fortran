
	subroutine jul2dy(lyr,idyjl,imon,idy)

	INTEGER totday(0:12),mnth(365),totdayly(0:12),mnthly(366)
	data totday/0,31,59,90,120,151,181,212,243,273,304,334,365/ 
	data totdayly/0,31,60,91,121,152,182,213,244,274,305,335,366/ 
	mnth=(/(1,i=1,31),(2,i=1,28),(3,i=1,31),(4,i=1,30),(5,i=1,31),(6,i=1,30),
	1	(7,i=1,31),(8,i=1,31),(9,i=1,30),(10,i=1,31),(11,i=1,30),(12,i=1,31)/)
	mnthly=(/(1,i=1,31),(2,i=1,29),(3,i=1,31),(4,i=1,30),(5,i=1,31),(6,i=1,30),
	1	(7,i=1,31),(8,i=1,31),(9,i=1,30),(10,i=1,31),(11,i=1,30),(12,i=1,31)/)
	
	if(nint(lyr/4.)*4==lyr) then
		imon=mnthly(idyjl)
		idy=idyjl-totdayly(imon-1)
	else
		imon=mnth(idyjl)
		idy=idyjl-totday(imon-1)
	endif

	return
	end


	subroutine dy2jul(lyr,idyjl,imon,idy)

	INTEGER totday(0:12),totdayly(0:12)
	data totday/0,31,59,90,120,151,181,212,243,273,304,334,365/ 
	data totdayly/0,31,60,91,121,152,182,213,244,274,305,335,366/ 

	if(nint(lyr/4.)*4==lyr) then
		idyjl=totdayly(imon-1)+idy
	else
		idyjl=totday(imon-1)+idy
	endif
	
	return
	end



	real function dyinyr(lyr)
	dyinyr=365
	if(nint(lyr/4.)*4==lyr) dyinyr=366
	return
	end

	subroutine day_increment(lyr,idyjl,increment,iyr,idy)

	idy=idyjl+increment;
	iyr=lyr;
	if(idy<1) then
		iyr=iyr-1;idy=idy+dyinyr(iyr)
	endif
	if(idy>dyinyr(iyr)) then
		idy=idy-dyinyr(iyr);iyr=iyr+1;
	endif

	end subroutine

	subroutine get_name_mon(imon,amon)
	integer imon
	character*3 amon

	if(imon==1) amon='JAN'
	if(imon==2) amon='FEB'
	if(imon==3) amon='MAR'
	if(imon==4) amon='APR'
	if(imon==5) amon='MAY'
	if(imon==6) amon='JUN'
	if(imon==7) amon='JUL'
	if(imon==8) amon='AUG'
	if(imon==9) amon='SEP'
	if(imon==10) amon='OCT'
	if(imon==11) amon='NOV'
	if(imon==12) amon='DEC'

	end subroutine
