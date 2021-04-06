	subroutine unzipfile(foutZ)
	USE MSFLIB
	character*(*) foutZ
	character*250 fout,unzip
	fout=foutZ
	ilen=LEN_TRIM(foutZ)
	fout(ilen-2:ilen)='   '
	write(unzip,1016) trim(foutZ), trim(fout)
	IOK=SYSTEMQQ(unzip)
1016	format('C:\GZIP.EXE -d -f -c ',A, ' > ', A)
	ilen=LEN_TRIM(foutZ)
	foutZ(ilen-2:ilen)='   '
	return
	end

	subroutine zipfile(foutZ)
	USE MSFLIB
      CHARACTER*(*) foutZ
	character*250 zip
	write(zip,1016) trim(foutZ)
	IOK=SYSTEMQQ(zip)
	ilen=LEN_TRIM(foutZ)
	foutZ(ilen+1:ilen+3)='.gz'
1016	format('C:\GZIP.EXE -f ',A)
	return
	end

	subroutine del_unzipfile(fname)
	USE MSFLIB
	character*(*) fname
	character*250 foutZ
	write(foutZ,1401)trim(fname)
1401	format('del ',A)
	IOK=SYSTEMQQ(foutZ)
	return
	end

	subroutine delfile(fname)
	USE MSFLIB
	character*(*) fname
	character*250 foutZ
	write(foutZ,1401)trim(fname)
1401	format('del ',A)
	IOK=SYSTEMQQ(foutZ)
	return
	end
	
	subroutine copyfile(fout,fdir)
	USE MSFLIB
	character*(*) fout,fdir
	character*400 fstr
	write(fstr,1016) trim(fout),trim(fdir)
c	print*, trim(fstr)
	IOK=SYSTEMQQ(fstr)
1016	format('copy /Y ',A,' ',A)
	return
	end

	subroutine unbzipfile(foutZ)
	USE MSFLIB
	character*(*) foutZ
	character*250 fout,unzip
	fout=foutZ
	ilen=LEN_TRIM(foutZ)
	fout(ilen-3:ilen)='   '
	write(unzip,1016) trim(foutZ), trim(fout)
	IOK=SYSTEMQQ(unzip)
1016	format('C:\bzip.exe -d -f -c ',A, ' > ', A)
	ilen=LEN_TRIM(foutZ)
	foutZ(ilen-3:ilen)='   '
	return
	end

	subroutine bzipfile(foutZ)
	USE MSFLIB
      CHARACTER*(*) foutZ
	character*250 zip
	write(zip,1016) trim(foutZ)
	IOK=SYSTEMQQ(zip)
1016	format('C:\bzip.exe -f ',A)
	return
	end

	

	subroutine unzipfileZ(foutZ)
	USE MSFLIB
	character*(*) foutZ
	character*250 fout,unzip
	fout=foutZ
	ilen=LEN_TRIM(foutZ)
	fout(ilen-1:ilen)='   '
	write(unzip,1016) trim(foutZ), trim(fout)
	IOK=SYSTEMQQ(unzip)
1016	format('C:\GZIP.EXE -d -f -c ',A, ' > ', A)
	ilen=LEN_TRIM(foutZ)
	foutZ(ilen-1:ilen)='   '
	return
	end

	subroutine movefile(fout,fdir)
	USE MSFLIB
	character*(*) fout,fdir
	character*420 fstr
	write(fstr,1016) trim(fout),trim(fdir)
	IOK=SYSTEMQQ(fstr)
1016	format('move /Y ',A,' ',A)
	return
	end
