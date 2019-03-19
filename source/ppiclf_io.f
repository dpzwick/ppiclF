!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteSubBinVTU(filein1)
#include "PPICLF"
      include 'mpif.h'

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile
      character*6  prostr

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      integer vtu,pth, nvtx_total, ncll_total
      integer*4 iint
      integer*8 idisp_pos
      integer*8 stride_lenv

      integer icount_pos(PPICLF_BX1, PPICLF_BY1, PPICLF_BZ1)

      real*4 rpoint(3)

      call ppiclf_printsi(' *Begin WriteSubBinVTU$',ppiclf_cycle)

      call ppiclf_prints('  -Begin CreateSubBin$')
         call ppiclf_comm_CreateSubBin
      call ppiclf_prints('   End CreateSubBin$')

      call ppiclf_prints('  -Begin CreateGhost$')
         call ppiclf_comm_CreateGhost
      call ppiclf_prints('   End CreateGhost$')

      call ppiclf_prints('  -Begin MoveGhost$')
         call ppiclf_comm_MoveGhost
      call ppiclf_prints('   End MoveGhost$')

      call ppiclf_prints('  -Begin ProjectParticleSubBin$')
         call ppiclf_solve_ProjectParticleSubBin
      call ppiclf_prints('   End ProjectParticleSubBin$')

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      if (ppiclf_ndim .gt. 2) then
         ncll_total = ppiclf_ndxgp*(ppiclf_bx-1)
     >               *ppiclf_ndygp*(ppiclf_by-1)
     >               *ppiclf_ndzgp*(ppiclf_bz-1)
         nvtx_total = (ppiclf_ndxgp*(ppiclf_bx-1)+1)
     >               *(ppiclf_ndygp*(ppiclf_by-1)+1)
     >               *(ppiclf_ndzgp*(ppiclf_bz-1)+1)
      else

        ! 2d here ....
      endif

      ppiclf_ndxgpp1 = ppiclf_ndxgp + 1
      ppiclf_ndygpp1 = ppiclf_ndygp + 1
      ppiclf_ndzgpp1 = ppiclf_ndzgp + 1
      ppiclf_ndxygpp1 = ppiclf_ndxgpp1*ppiclf_ndygpp1


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
      ibin = modulo(ppiclf_nid,ppiclf_ndxgp)
      jbin = modulo(ppiclf_nid/ppiclf_ndxgp,ppiclf_ndygp)
      kbin = ppiclf_nid/(ppiclf_ndxgp*ppiclf_ndygp)

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      do ie=1,PPICLF_LRP_PRO
         write(prostr,'(A4,I2.2)') "PRO-",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
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
      ndumx = ppiclf_ndxgp*(ppiclf_bx-1) + 1
      ndumy = ppiclf_ndygp*(ppiclf_by-1) + 1
      ndumz = ppiclf_ndzgp*(ppiclf_bz-1) + 1
      do ie=0,ppiclf_np-1
         i = modulo(ie,ppiclf_ndxgp)
         j = modulo(ie/ppiclf_ndxgp,ppiclf_ndygp)
         k = ie/(ppiclf_ndxgp*ppiclf_ndygp)


c        il = i-1
c        ir = i+1
c        jl = j-1
c        jr = j+1
c        kl = k-1
c        kr = k+1
c        nmult = ppiclf_ndxgp*ppiclf_ndygp
c        nid1  = il + ppiclf_ndxgp*jl + nmult*kl
c        nid2  = i  + ppiclf_ndxgp*jl + nmult*kl
c        nid3  = ir + ppiclf_ndxgp*jl + nmult*kl
c        nid4  = il + ppiclf_ndxgp*j  + nmult*kl
c        nid5  = i  + ppiclf_ndxgp*j  + nmult*kl
c        nid6  = ir + ppiclf_ndxgp*j  + nmult*kl
c        nid7  = il + ppiclf_ndxgp*jr + nmult*kl
c        nid8  = i  + ppiclf_ndxgp*jr + nmult*kl
c        nid9  = ir + ppiclf_ndxgp*jr + nmult*kl
c        nid10 = il + ppiclf_ndxgp*jl + nmult*k 
c        nid11 = i  + ppiclf_ndxgp*jl + nmult*k 
c        nid12 = ir + ppiclf_ndxgp*jl + nmult*k 
c        nid13 = il + ppiclf_ndxgp*j  + nmult*k 
c        nid14 = i  + ppiclf_ndxgp*j  + nmult*k 
c        nid15 = ir + ppiclf_ndxgp*j  + nmult*k 
c        nid16 = il + ppiclf_ndxgp*jr + nmult*k 
c        nid17 = i  + ppiclf_ndxgp*jr + nmult*k 
c        nid18 = ir + ppiclf_ndxgp*jr + nmult*k 
c        nid19 = il + ppiclf_ndxgp*jl + nmult*kr
c        nid20 = i  + ppiclf_ndxgp*jl + nmult*kr
c        nid21 = ir + ppiclf_ndxgp*jl + nmult*kr
c        nid22 = il + ppiclf_ndxgp*j  + nmult*kr
c        nid23 = i  + ppiclf_ndxgp*j  + nmult*kr
c        nid24 = ir + ppiclf_ndxgp*j  + nmult*kr
c        nid25 = il + ppiclf_ndxgp*jr + nmult*kr
c        nid26 = i  + ppiclf_ndxgp*jr + nmult*kr
c        nid27 = ir + ppiclf_ndxgp*jr + nmult*kr

      do kk=1,ppiclf_bz-1
      do jj=1,ppiclf_by-1
      do ii=1,ppiclf_bx-1

         itmp = i*(ppiclf_bx-1) + (ii-1)
         jtmp = j*(ppiclf_by-1) + (jj-1)
         ktmp = k*(ppiclf_bz-1) + (kk-1)

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

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         icount_pos(i,j,k) = 0
      enddo
      enddo
      enddo
      if (ppiclf_nid .le. ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by .and.
     >          k .ne. ppiclf_bz) then
                  icount_pos(i,j,k) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_ndxgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (k .ne. ppiclf_bz) then
            if (jbin .eq. ppiclf_ndygp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (kbin .eq. ppiclf_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_ndxgp-1) then
            if (jbin .eq. ppiclf_ndygp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (k .eq. ppiclf_bz) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_ndxgp-1) then
            if (kbin .eq. ppiclf_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_ndygp-1) then
            if (kbin .eq. ppiclf_ndzgp-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (ibin .eq. ppiclf_ndxgp-1) then
            if (jbin .eq. ppiclf_ndygp-1) then
            if (kbin .eq. ppiclf_ndzgp-1) then
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

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         stride_lenv = ppiclf_grid_i(i,j,k)
         idisp_pos   = ivtu_size + isize*(3   *stride_lenv + 1)
         icount_dum  = icount_pos(i,j,k)
         rpoint(1)   = ppiclf_grid_x(i,j,k)
         rpoint(2)   = ppiclf_grid_y(i,j,k)
         rpoint(3)   = ppiclf_grid_z(i,j,k)
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)

      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)


      do ie=1,PPICLF_LRP_PRO

      if_pos = 1*isize*nvtx_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
           stride_lenv = ppiclf_grid_i(i,j,k)
           idisp_pos   = ivtu_size + isize*(3*nvtx_total ! position fld
     >                   + (ie-1)*nvtx_total ! prev fields
     >                   + 1*stride_lenv    ! this fld
     >                   + 1 + ie)          ! ints
           icount_dum  = icount_pos(i,j,k)/3 ! either zero or 1
           rpoint(1)   = ppiclf_grid_fld(i,j,k,ie)
           call ppiclf_byte_set_view(idisp_pos,pth)
           call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)
      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      enddo
            

      ! still need to add 2d


      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteSubBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteBinVTU(filein1)
#include "PPICLF"
      include 'mpif.h'

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      integer vtu,pth, nvtx_total, ncll_total
      integer*4 iint
      integer*8 idisp_pos,idisp_cll
      integer*8 stride_lenv(8), stride_lenc

      real*4 rpoint(3)

      call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      nvtx_total = (ppiclf_ndxgp+1)*(ppiclf_ndygp+1)
      if (ppiclf_ndim .gt. 2) 
     >    nvtx_total = nvtx_total*(ppiclf_ndzgp+1)
      ncll_total = ppiclf_ndxgp*ppiclf_ndygp
      if (ppiclf_ndim .gt. 2) ncll_total = ncll_total*ppiclf_ndzgp


      ppiclf_ndxgpp1 = ppiclf_ndxgp + 1
      ppiclf_ndygpp1 = ppiclf_ndygp + 1
      ppiclf_ndzgpp1 = ppiclf_ndzgp + 1
      ppiclf_ndxygpp1 = ppiclf_ndxgpp1*ppiclf_ndygpp1


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
      ibin = modulo(ppiclf_nid,ppiclf_ndxgp)
      jbin = modulo(ppiclf_nid/ppiclf_ndxgp,ppiclf_ndygp)
      kbin = ppiclf_nid/(ppiclf_ndxgp*ppiclf_ndygp)

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin+1

      nbinpa = il + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kl
      nbinpb = ir + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kl
      nbinpc = il + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kl
      nbinpd = ir + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kl
      nbinpe = il + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kr
      nbinpf = ir + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kr
      nbinpg = il + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kr
      nbinph = ir + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kr

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (ppiclf_nid .le. ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         stride_lenv(5) = nbinpe
         stride_lenv(6) = nbinpf
         stride_lenv(7) = nbinpg
         stride_lenv(8) = nbinph
         
         stride_lenc    = ppiclf_nid
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
      iint = iint + 1   *isize*ncll_total + isize
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
         i = modulo(ii,ppiclf_ndxgp)
         j = modulo(ii/ppiclf_ndxgp,ppiclf_ndygp)
         k = ii/(ppiclf_ndxgp*ppiclf_ndygp)
          
c     do K=0,ppiclf_ndzgp-1
         kl = K
         kr = K+1
c     do J=0,ppiclf_ndygp-1
         jl = J
         jr = J+1
c     do I=0,ppiclf_ndxgp-1
         il = I
         ir = I+1
         npa = il + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kl
         npb = ir + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kl
         npc = il + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kl
         npd = ir + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kl
         npe = il + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kr
         npf = ir + ppiclf_ndxgpp1*jl + ppiclf_ndxygpp1*kr
         npg = il + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kr
         nph = ir + ppiclf_ndxgpp1*jr + ppiclf_ndxygpp1*kr
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

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total


      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      ! point A
      icount_pos = 0
      if (ppiclf_nid .le. 
     >    ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3   *stride_lenv(1) + 1)
      rpoint(1)  = ppiclf_binx(1,1)
      rpoint(2)  = ppiclf_biny(1,1)
      rpoint(3)  = ppiclf_binz(1,1)
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (ppiclf_ndim .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (ibin .eq. ppiclf_ndxgp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(2,1)
            rpoint(2)  = ppiclf_biny(1,1)
            rpoint(3)  = ppiclf_binz(1,1)
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (jbin .eq. ppiclf_ndygp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(1,1)
            rpoint(2)  = ppiclf_biny(2,1)
            rpoint(3)  = ppiclf_binz(1,1)
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (kbin .eq. ppiclf_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(1,1)
            rpoint(2)  = ppiclf_biny(1,1)
            rpoint(3)  = ppiclf_binz(2,1)
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (ibin .eq. ppiclf_ndxgp-1) then
         if (jbin .eq. ppiclf_ndygp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(2,1)
            rpoint(2)  = ppiclf_biny(2,1)
            rpoint(3)  = ppiclf_binz(1,1)
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (ibin .eq. ppiclf_ndxgp-1) then
         if (kbin .eq. ppiclf_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(2,1)
            rpoint(2)  = ppiclf_biny(1,1)
            rpoint(3)  = ppiclf_binz(2,1)
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (jbin .eq. ppiclf_ndygp-1) then
         if (kbin .eq. ppiclf_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(1,1)
            rpoint(2)  = ppiclf_biny(2,1)
            rpoint(3)  = ppiclf_binz(2,1)
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         if (ibin .eq. ppiclf_ndxgp-1) then
         if (jbin .eq. ppiclf_ndygp-1) then
         if (kbin .eq. ppiclf_ndzgp-1) then
            icount_pos = 3
            rpoint(1)  = ppiclf_binx(2,1)
            rpoint(2)  = ppiclf_biny(2,1)
            rpoint(3)  = ppiclf_binz(2,1)
         endif
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

      endif

      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      if_cll = 1*isize*ncll_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_cll
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      !
      ! cell values
      idisp_cll = ivtu_size + isize*(3  *(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (ppiclf_nid .le. ppiclf_ndxgp*ppiclf_ndygp*ppiclf_ndzgp-1) then
         icount_cll = 1
      endif
      rpoint(1)  = real(nxx)
      call ppiclf_byte_set_view(idisp_cll,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_cll,iorank,pth,ierr)

      call ppiclf_byte_close_mpi(pth,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleVTU(filein1)
#include "PPICLF"
      include 'mpif.h'

      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)

      character (len = *)  filein1
      character*3 filein
      character*12 vtufile

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      integer vtu,pth,prevs(2,ppiclf_np)
      integer*4 iint
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip
      integer*8 stride_len

      call ppiclf_printsi(' *Begin WriteParticleVTU$',ppiclf_cycle)

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      npt_total = ppiclf_iglsum(nxx,1)

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
      if_sln = PPICLF_LRS*isize*npt_total
      if_lrp = PPICLF_LRP*isize*npt_total
      if_lip = 3      *isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = ppiclf_y(jx,i)
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = ppiclf_y(jy,i)
c        if (if3d) then
            ic_pos = ic_pos + 1
            rout_pos(ic_pos) = ppiclf_y(jz,i)
c        else
c           ic_pos = ic_pos + 1
c           rout_pos(ic_pos) = 0.0
c        endif

         do j=1,PPICLF_LRS
            ic_sln = ic_sln + 1
            rout_sln(ic_sln) = ppiclf_y(j,i)
         enddo

         do j=1,PPICLF_LRP
            ic_lrp = ic_lrp + 1
            rout_lrp(ic_lrp) = ppiclf_rprop(j,i)
         enddo

         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = ppiclf_iprop(5,i)
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = ppiclf_iprop(6,i)
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = ppiclf_iprop(7,i)

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
      call fgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call fgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3   ,iint)
      iint = iint + 3   *isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'

      call ppiclf_io_WriteDataArrayVTU(vtu,'ppiclf-y',PPICLF_LRS,iint)
      iint = iint + PPICLF_LRS*isize*npt_total + isize

      call ppiclf_io_WriteDataArrayVTU(vtu,'ppiclf-rprop',PPICLF_LRP
     >                                ,iint)
      iint = iint + PPICLF_LRP*isize*npt_total + isize

      call ppiclf_io_WriteDataArrayVTU(vtu,'ppiclf-iprop',3,iint)
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

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
      idisp_sln = ivtu_size + isize*(3*npt_total + PPICLF_LRS*stride_len
     >                      + 2)
      idisp_lrp = ivtu_size + isize*(3*npt_total  + PPICLF_LRS*npt_total
     >                      + PPICLF_LRP*stride_len + 3)
      idisp_lip = ivtu_size + isize*(3*npt_total  + PPICLF_LRS*npt_total
     >                      + PPICLF_LRP*npt_total + 3*stride_len + 4 )

      ! how much to write
      icount_pos = 3      *nxx
      icount_sln = PPICLF_LRS*nxx
      icount_lrp = PPICLF_LRP*nxx
      icount_lip = 3      *nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_sln
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_sln,pth)
      call ppiclf_byte_write_mpi(rout_sln,icount_sln,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_lrp
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_lrp,pth)
      call ppiclf_byte_write_mpi(rout_lrp,icount_lrp,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_lip
        close(vtu)
      endif

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_lip,pth)
      call ppiclf_byte_write_mpi(rout_lip,icount_lip,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteDataArrayVTU(vtu,dataname,ncomp,idist)

      integer vtu,ncomp
      integer*4 idist
      character (len = *) dataname

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
      subroutine ppiclf_io_OutputDiagAll
#include "PPICLF"

      call ppiclf_io_OutputDiagGen
      call ppiclf_io_OutputDiagGhost
      if (ppiclf_lsubbin) call ppiclf_io_OutputDiagSubBin
      if (ppiclf_overlap) call ppiclf_io_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGen
#include "PPICLF"

      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = PPICLF_NDXGP*PPICLF_NDYGP*PPICLF_NDZGP

      call ppiclf_printsi('  -Cycle                  :$',ppiclf_cycle)
      call ppiclf_printsi('  -Output Freq.           :$',ppiclf_iostep)
      call ppiclf_printsr('  -Time                   :$',ppiclf_time)
      call ppiclf_printsr('  -dt                     :$',ppiclf_dt)
      call ppiclf_printsi('  -Global particles       :$',npart_tot)
      call ppiclf_printsi('  -Local particles (Max)  :$',npart_max)
      call ppiclf_printsi('  -Local particles (Min)  :$',npart_min)
      call ppiclf_printsi('  -Local particles (Ideal):$',npart_ide)
      call ppiclf_printsi('  -Total ranks            :$',ppiclf_np)
      call ppiclf_printsi('  -Problem dimensions     :$',ppiclf_ndim)
      call ppiclf_printsi('  -Integration method     :$',ppiclf_imethod)
      call ppiclf_printsi('  -Number of bins total   :$',nbin_total)
      call ppiclf_printsi('  -Number of bins (x)     :$',PPICLF_NDXGP)
      call ppiclf_printsi('  -Number of bins (y)     :$',PPICLF_NDYGP)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$',PPICLF_NDZGP)
      call ppiclf_printsr('  -Bin xl coordinate      :$',ppiclf_binb(1))
      call ppiclf_printsr('  -Bin xr coordinate      :$',ppiclf_binb(2))
      call ppiclf_printsr('  -Bin yl coordinate      :$',ppiclf_binb(3))
      call ppiclf_printsr('  -Bin yr coordinate      :$',ppiclf_binb(4))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zl coordinate      :$',ppiclf_binb(5))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zr coordinate      :$',ppiclf_binb(6))

      call ppiclf_prints('  End General Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGrid
#include "PPICLF"

      call ppiclf_prints(' *Begin Grid Info$')

         nel_max_orig   = ppiclf_iglmax(ppiclf_nee,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nee,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nee,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_neltb,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_neltb,1)
         nel_total_map = ppiclf_iglsum(ppiclf_neltb,1)

      call ppiclf_printsi('  -Orig. Global cells     :$',nel_total_orig)
      call ppiclf_printsi('  -Orig. Local cells (Max):$',nel_max_orig)
      call ppiclf_printsi('  -Orig. Local cells (Min):$',nel_min_orig)
      call ppiclf_printsi('  -Map Global cells       :$',nel_total_map)
      call ppiclf_printsi('  -Map Local cells (Max)  :$',nel_max_map)
      call ppiclf_printsi('  -Map Local cells (Min)  :$',nel_min_map)

      call ppiclf_prints('  End Grid Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGhost
#include "PPICLF"

      call ppiclf_prints(' *Begin Ghost Info$')

         npart_max = ppiclf_iglmax(ppiclf_npart_gp,1)
         npart_min = ppiclf_iglmin(ppiclf_npart_gp,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart_gp,1)

      call ppiclf_printsi('  -Global ghosts          :$',npart_tot)
      call ppiclf_printsi('  -Local ghosts (Max)     :$',npart_max)
      call ppiclf_printsi('  -Local ghosts (Min)     :$',npart_min)

      call ppiclf_prints('  End Ghost Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagSubBin
#include "PPICLF"

      call ppiclf_prints(' *Begin SubBin Info$')

      nbin_total = ppiclf_bx*ppiclf_by*ppiclf_bz

      call ppiclf_printsi('  -Number of local bin    :$',nbin_total)
      call ppiclf_printsi('  -Number of local bin (x):$',PPICLF_BX)
      call ppiclf_printsi('  -Number of local bin (y):$',PPICLF_BY)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of local bin (z):$',PPICLF_BZ)
      call ppiclf_printsr('  -Bin width (x)          :$',ppiclf_rdx)
      call ppiclf_printsr('  -Bin width (y)          :$',ppiclf_rdy)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin width (z)          :$',ppiclf_rdz)
      call ppiclf_printsr('  -Filter width           :$',ppiclf_filter)
      call ppiclf_printsr('  -Filter cut-off         :$'
     >                                                 ,ppiclf_d2chk(2))
      call ppiclf_printsi('  -SubBins per filter res.:$',ppiclf_ngrids)

      call ppiclf_prints('  End SubBin Info$')

      return
      end
!-----------------------------------------------------------------------
