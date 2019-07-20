#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
c-----------------------------------------------------------------------
      subroutine ppiclf_gop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      real*8 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ie
!
      if (op.eq.'+  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,ppiclf_comm,ie)
      elseif (op.EQ.'M  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_max,ppiclf_comm,ie)
      elseif (op.EQ.'m  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_min,ppiclf_comm,ie)
      elseif (op.EQ.'*  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,ppiclf_comm,ie)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_igop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      integer*4 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ierr
!
      if     (op.eq.'+  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
      elseif (op.EQ.'M  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
      elseif (op.EQ.'m  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
      elseif (op.EQ.'*  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglsum(a,n)
! 
      implicit none
! 
! Input:
! 
      integer*4 a(1)
      integer*4 n
! 
! Internal:
! 
      integer*4 tsum
      integer*4 tmp(1),work(1)
      integer*4 i
!
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call ppiclf_igop(tmp,work,'+  ',1)
      ppiclf_iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      real*8 function ppiclf_glsum (x,n)
! 
      implicit none
! 
! Vars:
! 
      real*8 x(1),tsum
      real*8 tmp(1),work(1)
      integer*4 i,n
!
      TSUM = 0.0d0
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL ppiclf_GOP(TMP,WORK,'+  ',1)
      ppiclf_GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      REAL*8 A(1),tmax
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL ppiclf_GOP(TMP,WORK,'M  ',1)
      ppiclf_GLMAX=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      integer*4 a(1),tmax
      integer*4 tmp(1),work(1)
      integer*4 i,n
!
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call ppiclf_igop(tmp,work,'M  ',1)
      ppiclf_iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      REAL*8 A(1),tmin
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL ppiclf_GOP(TMP,WORK,'m  ',1)
      ppiclf_GLMIN = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      integer*4 a(1),tmin
      integer*4 tmp(1),work(1)
      integer*4 i, n
!
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call ppiclf_igop(tmp,work,'m  ',1)
      ppiclf_iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmin(vec,n)
! 
      implicit none
! 
! Vars:
! 
      REAL*8 VEC(1),tmin
      integer*4 i, n
!
      TMIN = 99.0E20
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      ppiclf_VLMIN = TMIN
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmax(vec,n)
! 
      implicit none
! 
! Vars:
! 
      REAL*8 VEC(1),tmax
      integer*4 i, n
!
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      ppiclf_VLMAX = TMAX
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_copy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      real*8 a(1),b(1)
      integer*4 i,n
!

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_icopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      INTEGER*4 A(1), B(1)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_chcopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      CHARACTER*1 A(1), B(1)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_exittr(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
! 
! Vars:
! 
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, ierr, k, idata
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ERROR ',132a1)

c     call mpi_finalize (ierr)
      call mpi_abort(ppiclf_comm, 1, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsri(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!

      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsi(stringi,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*10 s10
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s10,10) idata
   10 format(1x,i9)
      call ppiclf_chcopy(stringo(ilen),s10,10)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+9)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsr(stringi,rdata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*15 s15
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s15,15) rdata
   15 format(1x,1p1e14.6)
      call ppiclf_chcopy(stringo(ilen),s15,15)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+14)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_prints(stringi)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen-1)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE PPICLF_BLANK(A,N)
! 
      implicit none
! 
! Vars:
!
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
      integer*4 i,n
!
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER*4 FUNCTION PPICLF_INDX1(S1,S2)
! 
      implicit none
! 
! Vars:
!
      CHARACTER*1 S1(132),S2
      integer*4 n1, i
!
      N1=132
      PPICLF_INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         IF (S1(I).EQ.S2) THEN
            PPICLF_INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      character*132 FUNCTION PPICLF_CHSTR(S1)
! 
      implicit none
! 
! Vars:
!
      CHARACTER*132 S1
!
      PPICLF_CHSTR = S1
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character fnamei*(*)
      logical ifro
      CHARACTER*1 BLNK
      DATA BLNK/' '/
      character*132 fname
      character*1   fname1(132)
      equivalence  (fname1,fname)
      integer*4 imode, ierr, mpi_fh
!
      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY 
      endif

      call MPI_file_open(ppiclf_comm,fnamei,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_byte_read_mpi(buf,icount,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 iout, icount, mpi_fh, ierr
!

      iout = icount ! icount is in 4-byte words
      call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
     &                       MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 icount, iorank, mpi_fh, ierr, iout
!

      iout = icount ! icount is in 4-byte words
      if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                        MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_close_mpi(mpi_fh,ierr)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*4 mpi_fh, ierr
!

      call MPI_file_close(mpi_fh,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_set_view(ioff_in,mpi_fh)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*8 ioff_in
      integer*4 mpi_fh, ierr
!
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                       'native',MPI_INFO_NULL,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_bcast(buf,len)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)
      integer*4 len, ierr
!

      call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

      return
      end
C--------------------------------------------------------------------------
