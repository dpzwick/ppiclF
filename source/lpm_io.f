!-----------------------------------------------------------------------
      subroutine lpm_io_vtu_write_grd(filein1,iobig)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      real*4  rout_pos(3      *LPM_LPART) 
     >       ,rout_sln(LPM_LRS*LPM_LPART)
     >       ,rout_lrp(LPM_LRP*LPM_LPART)
     >       ,rout_lip(3      *LPM_LPART)

      character*5 sprop1
      character*9 rprop1

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile
      character*13 vtufile1
      character*50 dumstr
      character*6  prostr

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         
      integer vtu,vtu1,pth, nvtx_total, ncll_total
      integer*4 iint
      integer*8 idisp_pos,idisp_cll,idisp_lrp,idisp_lip
      integer*8 stride_lenv, stride_lenc

      integer icount_pos(LPM_BX1, LPM_BY1, LPM_BZ1)

      real*4 rpoint(3)

      call lpm_comm_ghost_create
      call lpm_comm_ghost_send
      call lpm_solve_project_bins

      icalld1 = icalld1+1

      nnp   = lpm_np
      nxx   = LPM_NPART

      if (lpm_rparam(12) .gt. 2) then
         ncll_total = lpm_ndxgp*(lpm_bx-1)
     >               *lpm_ndygp*(lpm_by-1)
     >               *lpm_ndzgp*(lpm_bz-1)
         nvtx_total = (lpm_ndxgp*(lpm_bx-1)+1)
     >               *(lpm_ndygp*(lpm_by-1)+1)
     >               *(lpm_ndzgp*(lpm_bz-1)+1)
      else

        ! 2d here ....
      endif

      lpm_ndxgpp1 = lpm_ndxgp + 1
      lpm_ndygpp1 = lpm_ndygp + 1
      lpm_ndzgpp1 = lpm_ndzgp + 1
      lpm_ndxygpp1 = lpm_ndxgpp1*lpm_ndygpp1


      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'grd'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize  = 4

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------

      ! get which bin this processor holds
      ibin = modulo(lpm_nid,lpm_ndxgp)
      jbin = modulo(lpm_nid/lpm_ndxgp,lpm_ndygp)
      kbin = lpm_nid/(lpm_ndxgp*lpm_ndygp)

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

      if (lpm_nid .eq. 0) then

      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (iobig .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (iobig .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') lpm_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') lpm_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call lpm_io_vtu_data(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      do ie=1,LPM_LRP_PRO
         write(prostr,'(A4,I2.2)') "PRO-",ie
         call lpm_io_vtu_data(vtu,prostr,1,iint)
         iint = iint + 1*isize*nvtx_total + isize
      enddo
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      ndumx = lpm_ndxgp*(lpm_bx-1) + 1
      ndumy = lpm_ndygp*(lpm_by-1) + 1
      ndumz = lpm_ndzgp*(lpm_bz-1) + 1
      do ie=0,lpm_np-1
         i = modulo(ie,lpm_ndxgp)
         j = modulo(ie/lpm_ndxgp,lpm_ndygp)
         k = ie/(lpm_ndxgp*lpm_ndygp)


c        il = i-1
c        ir = i+1
c        jl = j-1
c        jr = j+1
c        kl = k-1
c        kr = k+1
c        nmult = lpm_ndxgp*lpm_ndygp
c        nid1  = il + lpm_ndxgp*jl + nmult*kl
c        nid2  = i  + lpm_ndxgp*jl + nmult*kl
c        nid3  = ir + lpm_ndxgp*jl + nmult*kl
c        nid4  = il + lpm_ndxgp*j  + nmult*kl
c        nid5  = i  + lpm_ndxgp*j  + nmult*kl
c        nid6  = ir + lpm_ndxgp*j  + nmult*kl
c        nid7  = il + lpm_ndxgp*jr + nmult*kl
c        nid8  = i  + lpm_ndxgp*jr + nmult*kl
c        nid9  = ir + lpm_ndxgp*jr + nmult*kl
c        nid10 = il + lpm_ndxgp*jl + nmult*k 
c        nid11 = i  + lpm_ndxgp*jl + nmult*k 
c        nid12 = ir + lpm_ndxgp*jl + nmult*k 
c        nid13 = il + lpm_ndxgp*j  + nmult*k 
c        nid14 = i  + lpm_ndxgp*j  + nmult*k 
c        nid15 = ir + lpm_ndxgp*j  + nmult*k 
c        nid16 = il + lpm_ndxgp*jr + nmult*k 
c        nid17 = i  + lpm_ndxgp*jr + nmult*k 
c        nid18 = ir + lpm_ndxgp*jr + nmult*k 
c        nid19 = il + lpm_ndxgp*jl + nmult*kr
c        nid20 = i  + lpm_ndxgp*jl + nmult*kr
c        nid21 = ir + lpm_ndxgp*jl + nmult*kr
c        nid22 = il + lpm_ndxgp*j  + nmult*kr
c        nid23 = i  + lpm_ndxgp*j  + nmult*kr
c        nid24 = ir + lpm_ndxgp*j  + nmult*kr
c        nid25 = il + lpm_ndxgp*jr + nmult*kr
c        nid26 = i  + lpm_ndxgp*jr + nmult*kr
c        nid27 = ir + lpm_ndxgp*jr + nmult*kr

      do kk=1,lpm_bz-1
      do jj=1,lpm_by-1
      do ii=1,lpm_bx-1

         itmp = i*(lpm_bx-1) + (ii-1)
         jtmp = j*(lpm_by-1) + (jj-1)
         ktmp = k*(lpm_bz-1) + (kk-1)

         kl = ktmp
         kr = ktmp+1
         jl = jtmp
         jr = jtmp+1
         il = itmp
         ir = itmp+1

c        ndum = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         npa = il + ndumx*jl + ndumx*ndumy*kl
         npb = ir + ndumx*jl + ndumx*ndumy*kl
         npc = il + ndumx*jr + ndumx*ndumy*kl
         npd = ir + ndumx*jr + ndumx*ndumy*kl
         npe = il + ndumx*jl + ndumx*ndumy*kr
         npf = ir + ndumx*jl + ndumx*ndumy*kr
         npg = il + ndumx*jr + ndumx*ndumy*kr
         nph = ir + ndumx*jr + ndumx*ndumy*kr
         write(vtu,'(I0,A)',advance='no')  npa, ' '
         write(vtu,'(I0,A)',advance='no')  npb, ' '
         write(vtu,'(I0,A)',advance='no')  npc, ' '
         write(vtu,'(I0,A)',advance='no')  npd, ' '
         write(vtu,'(I0,A)',advance='no')  npe, ' '
         write(vtu,'(I0,A)',advance='no')  npf, ' '
         write(vtu,'(I0,A)',advance='no')  npg, ' '
         write(vtu,'(I0)'  ,advance='yes') nph
      enddo
      enddo
      enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') 8*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') 11
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(lpm_comm,ierr)

      ! write points first
      call byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,lpm_bz
      do j=1,lpm_by
      do i=1,lpm_bx
         icount_pos(i,j,k) = 0
      enddo
      enddo
      enddo
      if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         do k=1,lpm_bz
         do j=1,lpm_by
         do i=1,lpm_bx
            if (i .ne. lpm_bx .and.
     >          j .ne. lpm_by .and.
     >          k .ne. lpm_bz) then
                  icount_pos(i,j,k) = 3
            endif

            if (i .eq. lpm_bx) then
            if (j .ne. lpm_by) then
            if (k .ne. lpm_bz) then
            if (ibin .eq. lpm_ndxgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (j .eq. lpm_by) then
            if (i .ne. lpm_bx) then
            if (k .ne. lpm_bz) then
            if (jbin .eq. lpm_ndygp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (k .eq. lpm_bz) then
            if (i .ne. lpm_bx) then
            if (j .ne. lpm_by) then
            if (kbin .eq. lpm_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (i .eq. lpm_bx) then
            if (j .eq. lpm_by) then
            if (k .ne. lpm_bz) then
            if (ibin .eq. lpm_ndxgp-1) then
            if (jbin .eq. lpm_ndygp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. lpm_bx) then
            if (k .eq. lpm_bz) then
            if (j .ne. lpm_by) then
            if (ibin .eq. lpm_ndxgp-1) then
            if (kbin .eq. lpm_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (j .eq. lpm_by) then
            if (k .eq. lpm_bz) then
            if (i .ne. lpm_bx) then
            if (jbin .eq. lpm_ndygp-1) then
            if (kbin .eq. lpm_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. lpm_bx) then
            if (j .eq. lpm_by) then
            if (k .eq. lpm_bz) then
            if (ibin .eq. lpm_ndxgp-1) then
            if (jbin .eq. lpm_ndygp-1) then
            if (kbin .eq. lpm_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif
            endif

         enddo
         enddo
         enddo
      endif

      do k=1,lpm_bz
      do j=1,lpm_by
      do i=1,lpm_bx
         stride_lenv = lpm_grid_i(i,j,k)
         idisp_pos   = ivtu_size + isize*(3   *stride_lenv + 1)
         icount_dum  = icount_pos(i,j,k)
         rpoint(1)   = lpm_grid_x(i,j,k)
         rpoint(2)   = lpm_grid_y(i,j,k)
         rpoint(3)   = lpm_grid_z(i,j,k)
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)

      enddo
      enddo
      enddo

      call byte_close_mpi(pth,ierr)


      do ie=1,LPM_LRP_PRO

      if_pos = 1*isize*nvtx_total

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(lpm_comm,ierr)

      call byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,lpm_bz
      do j=1,lpm_by
      do i=1,lpm_bx
            stride_lenv = lpm_grid_i(i,j,k)
            idisp_pos   = ivtu_size + isize*(3*nvtx_total ! position fld
     >                    + (ie-1)*nvtx_total ! prev fields
     >                    + 1*stride_lenv    ! this fld
     >                    + 1 + ie)          ! ints
            icount_dum  = icount_pos(i,j,k)/3 ! either zero or 1
            rpoint(1)   = lpm_grid_fld(i,j,k,ie)
            call byte_set_view(idisp_pos,pth)
            call byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)
      enddo
      enddo
      enddo

      call byte_close_mpi(pth,ierr)

      enddo
            

      ! still need to add 2d


      call mpi_barrier(lpm_comm,ierr)

      if (lpm_nid .eq. 0) then
      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_io_vtu_write_bins(filein1,iobig)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      real*4  rout_pos(3      *LPM_LPART) 
     >       ,rout_sln(LPM_LRS*LPM_LPART)
     >       ,rout_lrp(LPM_LRP*LPM_LPART)
     >       ,rout_lip(3      *LPM_LPART)

      character*5 sprop1
      character*9 rprop1

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile
      character*13 vtufile1
      character*50 dumstr

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         
      integer vtu,vtu1,pth, nvtx_total, ncll_total
      integer*4 iint
      integer*8 idisp_pos,idisp_cll,idisp_lrp,idisp_lip
      integer*8 stride_lenv(8), stride_lenc

      real*4 rpoint(3)

      icalld1 = icalld1+1

      nnp   = lpm_np
      nxx   = LPM_NPART

      nvtx_total = (lpm_ndxgp+1)*(lpm_ndygp+1)
      if (lpm_rparam(12) .gt. 2) nvtx_total = nvtx_total*(lpm_ndzgp+1)
      ncll_total = lpm_ndxgp*lpm_ndygp
      if (lpm_rparam(12) .gt. 2) ncll_total = ncll_total*lpm_ndzgp


      lpm_ndxgpp1 = lpm_ndxgp + 1
      lpm_ndygpp1 = lpm_ndygp + 1
      lpm_ndzgpp1 = lpm_ndzgp + 1
      lpm_ndxygpp1 = lpm_ndxgpp1*lpm_ndygpp1


      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'bin'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize  = 4

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------

      ! get which bin this processor holds
      ibin = modulo(lpm_nid,lpm_ndxgp)
      jbin = modulo(lpm_nid/lpm_ndxgp,lpm_ndygp)
      kbin = lpm_nid/(lpm_ndxgp*lpm_ndygp)

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin+1

      nbinpa = il + lpm_ndxgpp1*jl + lpm_ndxygpp1*kl
      nbinpb = ir + lpm_ndxgpp1*jl + lpm_ndxygpp1*kl
      nbinpc = il + lpm_ndxgpp1*jr + lpm_ndxygpp1*kl
      nbinpd = ir + lpm_ndxgpp1*jr + lpm_ndxygpp1*kl
      nbinpe = il + lpm_ndxgpp1*jl + lpm_ndxygpp1*kr
      nbinpf = ir + lpm_ndxgpp1*jl + lpm_ndxygpp1*kr
      nbinpg = il + lpm_ndxgpp1*jr + lpm_ndxygpp1*kr
      nbinph = ir + lpm_ndxgpp1*jr + lpm_ndxygpp1*kr

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         stride_lenv(5) = nbinpe
         stride_lenv(6) = nbinpf
         stride_lenv(7) = nbinpg
         stride_lenv(8) = nbinph
         
         stride_lenc    = lpm_nid
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

      if (lpm_nid .eq. 0) then

      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (iobig .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (iobig .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') lpm_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') lpm_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call lpm_io_vtu_data(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      do ii=0,ncll_total-1
         i = modulo(ii,lpm_ndxgp)
         j = modulo(ii/lpm_ndxgp,lpm_ndygp)
         k = ii/(lpm_ndxgp*lpm_ndygp)
          
c     do K=0,lpm_ndzgp-1
         kl = K
         kr = K+1
c     do J=0,lpm_ndygp-1
         jl = J
         jr = J+1
c     do I=0,lpm_ndxgp-1
         il = I
         ir = I+1
         npa = il + lpm_ndxgpp1*jl + lpm_ndxygpp1*kl
         npb = ir + lpm_ndxgpp1*jl + lpm_ndxygpp1*kl
         npc = il + lpm_ndxgpp1*jr + lpm_ndxygpp1*kl
         npd = ir + lpm_ndxgpp1*jr + lpm_ndxygpp1*kl
         npe = il + lpm_ndxgpp1*jl + lpm_ndxygpp1*kr
         npf = ir + lpm_ndxgpp1*jl + lpm_ndxygpp1*kr
         npg = il + lpm_ndxgpp1*jr + lpm_ndxygpp1*kr
         nph = ir + lpm_ndxgpp1*jr + lpm_ndxygpp1*kr
         write(vtu,'(I0,A)',advance='no')  npa, ' '
         write(vtu,'(I0,A)',advance='no')  npb, ' '
         write(vtu,'(I0,A)',advance='no')  npc, ' '
         write(vtu,'(I0,A)',advance='no')  npd, ' '
         write(vtu,'(I0,A)',advance='no')  npe, ' '
         write(vtu,'(I0,A)',advance='no')  npf, ' '
         write(vtu,'(I0,A)',advance='no')  npg, ' '
         write(vtu,'(I0)'  ,advance='yes') nph
c     enddo
c     enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') 8*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') 11
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total


      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(lpm_comm,ierr)

      ! write points first
      call byte_open_mpi(vtufile,pth,.false.,ierr)

      ! point A
      icount_pos = 0
      if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3   *stride_lenv(1) + 1)
      rpoint(1)  = lpm_binx(1,1)
      rpoint(2)  = lpm_biny(1,1)
      rpoint(3)  = lpm_binz(1,1)
      call byte_set_view(idisp_pos,pth)
      call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (lpm_rparam(12) .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (ibin .eq. lpm_ndxgp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(2,1)
            rpoint(2)  = lpm_biny(1,1)
            rpoint(3)  = lpm_binz(1,1)
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (jbin .eq. lpm_ndygp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(1,1)
            rpoint(2)  = lpm_biny(2,1)
            rpoint(3)  = lpm_binz(1,1)
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (kbin .eq. lpm_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(1,1)
            rpoint(2)  = lpm_biny(1,1)
            rpoint(3)  = lpm_binz(2,1)
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (ibin .eq. lpm_ndxgp-1) then
         if (jbin .eq. lpm_ndygp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(2,1)
            rpoint(2)  = lpm_biny(2,1)
            rpoint(3)  = lpm_binz(1,1)
         endif
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (ibin .eq. lpm_ndxgp-1) then
         if (kbin .eq. lpm_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(2,1)
            rpoint(2)  = lpm_biny(1,1)
            rpoint(3)  = lpm_binz(2,1)
         endif
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (jbin .eq. lpm_ndygp-1) then
         if (kbin .eq. lpm_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(1,1)
            rpoint(2)  = lpm_biny(2,1)
            rpoint(3)  = lpm_binz(2,1)
         endif
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         if (ibin .eq. lpm_ndxgp-1) then
         if (jbin .eq. lpm_ndygp-1) then
         if (kbin .eq. lpm_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = lpm_binx(2,1)
            rpoint(2)  = lpm_biny(2,1)
            rpoint(3)  = lpm_binz(2,1)
         endif
         endif
         endif
         endif
         call byte_set_view(idisp_pos,pth)
         call byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

      endif

      call byte_close_mpi(pth,ierr)

      call mpi_barrier(lpm_comm,ierr)

      if_cll = 1*isize*ncll_total

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_cll
        close(vtu)
      endif

      call mpi_barrier(lpm_comm,ierr)

      ! write points first
      call byte_open_mpi(vtufile,pth,.false.,ierr)

      !
      ! cell values
      idisp_cll = ivtu_size + isize*(3  *(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (lpm_nid .le. lpm_ndxgp*lpm_ndygp*lpm_ndzgp-1) then
         icount_cll = 1
      endif
      rpoint(1)  = real(nxx)
      call byte_set_view(idisp_cll,pth)
      call byte_write_mpi(rpoint,icount_cll,iorank,pth,ierr)

      call byte_close_mpi(pth,ierr)

      if (lpm_nid .eq. 0) then
      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_io_vtu_write(filein1,iobig)
#include "lpm_user.h"
#include "lpm.h"
#include "LPM"
      include 'mpif.h'

      real*4  rout_pos(3      *LPM_LPART) 
     >       ,rout_sln(LPM_LRS*LPM_LPART)
     >       ,rout_lrp(LPM_LRP*LPM_LPART)
     >       ,rout_lip(3      *LPM_LPART)

      character*5 sprop1
      character*9 rprop1

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile
      character*13 vtufile1
      character*50 dumstr

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         
      integer vtu,vtu1,pth,prevs(2,lpm_np)
      integer*4 iint
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip
      integer*8 stride_len

      icalld1 = icalld1+1

      nnp   = lpm_np
      nxx   = LPM_NPART

      npt_total = iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 3

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'par'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      iadd = 0
      if_pos = 3      *isize*npt_total
      if_sln = LPM_LRS*isize*npt_total
      if_lrp = LPM_LRP*isize*npt_total
      if_lip = 3      *isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = lpm_y(jx,i)
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = lpm_y(jy,i)
c        if (if3d) then
            ic_pos = ic_pos + 1
            rout_pos(ic_pos) = lpm_y(jz,i)
c        else
c           ic_pos = ic_pos + 1
c           rout_pos(ic_pos) = 0.0
c        endif

         do j=1,LPM_LRS
            ic_sln = ic_sln + 1
            rout_sln(ic_sln) = lpm_y(j,i)
         enddo

         do j=1,LPM_LRP
            ic_lrp = ic_lrp + 1
            rout_lrp(ic_lrp) = lpm_rprop(j,i)
         enddo

         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = lpm_iprop(5,i)
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = lpm_iprop(6,i)
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = lpm_iprop(7,i)

      enddo

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call fgslib_crystal_ituple_transfer(i_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call fgslib_crystal_ituple_sort(i_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (lpm_nid .ne. 0) then
      do i=1,lpm_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (lpm_nid .eq. 0) then

      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (iobig .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (iobig .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') lpm_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') lpm_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call lpm_io_vtu_data(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'

      call lpm_io_vtu_data(vtu,'lpm-y',LPM_LRS,iint)
      iint = iint + LPM_LRS*isize*npt_total + isize

      call lpm_io_vtu_data(vtu,'lpm-rprop',LPM_LRP,iint)
      iint = iint + LPM_LRP*isize*npt_total + isize

      call lpm_io_vtu_data(vtu,'lpm-iprop',3,iint)
      iint = iint + 3*isize*npt_total + isize

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3   *stride_len + 1)
      idisp_sln = ivtu_size + isize*(3   *npt_total + LPM_LRS*stride_len
     >                      + 2)
      idisp_lrp = ivtu_size + isize*(3   *npt_total  + LPM_LRS*npt_total
     >                      + LPM_LRP*stride_len + 3)
      idisp_lip = ivtu_size + isize*(3   *npt_total  + LPM_LRS*npt_total
     >                      + LPM_LRP*npt_total + 3*stride_len + 4 )

      ! how much to write
      icount_pos = 3      *nxx
      icount_sln = LPM_LRS*nxx
      icount_lrp = LPM_LRP*nxx
      icount_lip = 3      *nxx

      iorank = -1

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(lpm_comm,ierr)

      ! write
      call byte_open_mpi(vtufile,pth,.false.,ierr)
      call byte_set_view(idisp_pos,pth)
      call byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      call mpi_barrier(lpm_comm,ierr)

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_sln
        close(vtu)
      endif

      call mpi_barrier(lpm_comm,ierr)

      ! write
      call byte_open_mpi(vtufile,pth,.false.,ierr)
      call byte_set_view(idisp_sln,pth)
      call byte_write_mpi(rout_sln,icount_sln,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_lrp
        close(vtu)
      endif

      call mpi_barrier(lpm_comm,ierr)

      ! write
      call byte_open_mpi(vtufile,pth,.false.,ierr)
      call byte_set_view(idisp_lrp,pth)
      call byte_write_mpi(rout_lrp,icount_lrp,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      ! integer write
      if (lpm_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_lip
        close(vtu)
      endif

      ! write
      call byte_open_mpi(vtufile,pth,.false.,ierr)
      call byte_set_view(idisp_lip,pth)
      call byte_write_mpi(rout_lip,icount_lip,iorank,pth,ierr)
      call byte_close_mpi(pth,ierr)

      if (lpm_nid .eq. 0) then
      vtu=867+lpm_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_io_vtu_data(vtu,dataname,ncomp,idist)

      integer vtu,ncomp
      integer*4 idist
      character (len = *) dataname
      character*50 dumstr

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="'
      write(vtu,'(A)',advance='no') dataname
      write(vtu,'(A)',advance='no') '" NumberOfComponents="'
      write(vtu,'(I0)',advance='no') ncomp
      write(vtu,'(A)',advance='no') '" format="append" '
      write(vtu,'(A)',advance='no') 'offset="'
      write(vtu,'(I0)',advance='no') idist
      write(vtu,'(A)',advance='yes') '"/>'

      return
      end
!-----------------------------------------------------------------------
