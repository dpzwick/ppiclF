c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      real x(n), w(n)
      character*3 op

      if (op.eq.'+  ') then
      call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,lpm_comm,ie)
      elseif (op.EQ.'M  ') then
      call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_max,lpm_comm,ie)
      elseif (op.EQ.'m  ') then
      call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_min,lpm_comm,ie)
      elseif (op.EQ.'*  ') then
      call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,lpm_comm
     >                                                              ,ie)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine igop( x, w, op, n)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      integer x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,lpm_comm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,lpm_comm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,lpm_comm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,lpm_comm,ierr)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call igop(tmp,work,'+  ',1)
      iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      function glsum (x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function glmax(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLMAX=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call igop(tmp,work,'M  ',1)
      iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function glmin(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL GOP(TMP,WORK,'m  ',1)
      GLMIN = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call igop(tmp,work,'m  ',1)
      iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
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

      call MPI_file_open(lpm_comm,fnamei,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)

      return
      end
C--------------------------------------------------------------------------
c     subroutine byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)
c     include 'mpif.h'

c     real*4 buf(1)          ! buffer

c     iout = icount ! icount is in 4-byte words
c     call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
c    &                       MPI_STATUS_IGNORE,ierr)
c--------------------------------------------------------------------------
      subroutine byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      real*4 buf(1)          ! buffer

      iout = icount ! icount is in 4-byte words
      if(iorank.ge.0 .and. lpm_nid.ne.iorank) iout = 0
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                        MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine byte_close_mpi(mpi_fh,ierr)
      include 'mpif.h'

      call MPI_file_close(mpi_fh,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine byte_set_view(ioff_in,mpi_fh)
      include 'mpif.h'
      integer*8 ioff_in
    
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                       'native',MPI_INFO_NULL,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine bcast(buf,len)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,0,lpm_comm,ierr)

      return
      end
C--------------------------------------------------------------------------
