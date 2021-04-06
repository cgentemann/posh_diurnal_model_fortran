C     The 12/15/2000 version was changed on 2/22/2002.
C     The old version is saved as OPENBIG_V01
C     The changes are:
C     1.  No more Lahey code
C     2.  SHARE='DENYNONE' is added, which is the default for MS Fortran but not for new Compaq
C     3. 	IOSTAT=IERROR is added to give more info when error occurs
C     4.  Routine will retry if user requests.

      SUBROUTINE OPENBIG(IUNIT,FILENAME,S)

      CHARACTER*(*) FILENAME,S

    5 CONTINUE

      OPEN(IUNIT,FILE=FILENAME,STATUS=S,BLOCKSIZE=32000,ACCESS='SEQUENTIAL',FORM='BINARY',
     &     SHARE='DENYNONE',IOSTAT=IERROR,ERR=10)
      RETURN
C
   10 CONTINUE
      WRITE(*,1001) FILENAME,IUNIT,S,IERROR
 1001 FORMAT('ERROR IN OPENBIG ',A50,1X,I3,1X,A10,I7)
      WRITE(*,*) 'ENTER 0 TO TRY AGAIN, ENTER 1 TO STOP'
	READ(*,*) IGO
	IF(IGO.EQ.1) STOP 'PROGRAM STOPPED BY USER'
	GOTO 5
      END