c-----------------------------------------------------------------------
      subroutine ppiclf_gop( x, w, op, n)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      include 'mpif.h'

      real x(n), w(n)
      character*3 op

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
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      include 'mpif.h'

      integer x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function ppiclf_iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
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
      function ppiclf_glsum (x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL ppiclf_GOP(TMP,WORK,'+  ',1)
      ppiclf_GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function ppiclf_glmax(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
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
      function ppiclf_iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
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
      function ppiclf_glmin(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
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
      function ppiclf_iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
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
      real function ppiclf_vlmin(vec,n)
      REAL VEC(1)
      TMIN = 99.0E20
C
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      VLMIN = TMIN
      return
      END
c-----------------------------------------------------------------------
      real function ppiclf_vlmax(vec,n)
      REAL VEC(1)
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      VLMAX = TMAX
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      include 'mpif.h'

      character fnamei*(*)
      logical ifro
      CHARACTER*1 BLNK
      DATA BLNK/' '/
 
      character*132 fname
      character*1   fname1(132)
      equivalence  (fname1,fname)

c     write(6,*) fnamei

      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY 
      endif

      call MPI_file_open(ppiclf_comm,fnamei,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)

      return
      end
C--------------------------------------------------------------------------
c     subroutine ppiclf_byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)
c     include 'mpif.h'

c     real*4 buf(1)          ! buffer

c     iout = icount ! icount is in 4-byte words
c     call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
c    &                       MPI_STATUS_IGNORE,ierr)
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      include 'mpif.h'

      real*4 buf(1)          ! buffer

      iout = icount ! icount is in 4-byte words
      if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                        MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_close_mpi(mpi_fh,ierr)
      include 'mpif.h'

      call MPI_file_close(mpi_fh,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_set_view(ioff_in,mpi_fh)
      include 'mpif.h'
      integer*8 ioff_in
    
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                       'native',MPI_INFO_NULL,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_bcast(buf,len)
#include "ppiclf_user.h"
#include "ppiclf.h"
#include "PPICLF"
      include 'mpif.h'
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

      return
      end
C--------------------------------------------------------------------------
