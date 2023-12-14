!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin)
     > bind(C, name="ppiclc_io_ReadParticleVTU")
#else
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      character*1 dum_read
      character*132 filein2
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 vtu, isize, jx, jy, jz, ivtu_size, ifound, i, npt_total,
     >          npart_min, npmax, npart, ndiff, iorank, icount_pos,
     >          icount_sln, icount_lrp, icount_lip, j, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, pth, ierr, indx1
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      EXTERNAL PPICLF_CHSTR
      ! Sam - this modifies the interface for ppiclC. Will now need to
      ! include NULL for optional arguments
      integer*4, optional :: istartoutin
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!
      call ppiclf_prints(' *Begin ReadParticleVTU$')
      
      call ppiclf_solve_InitZero
      PPICLF_RESTART = .true.

      if (present(istartoutin)) then
        istartout = istartoutin
      else
        istartout = 0
      end if


      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) u (3)
      filein2 = ppiclf_chstr(filein1(1:indx1))

      isize = 4
      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=trim(filein2)
     >    ,access='stream',form="unformatted")

      ivtu_size = -1
      ifound = 0
      do i=1,1000000
      read(vtu) dum_read
      if (dum_read == '_') ifound = ifound + 1
      if (ifound .eq. 2) then
         ivtu_size = i
         exit
      endif
      enddo
      read(vtu) npt_total
      close(vtu)
      npt_total = npt_total/isize/3
      endif

      call ppiclf_bcast(npt_total,isize)
      call ppiclf_bcast(ivtu_size,isize)


      npart_min = npt_total/ppiclf_np+1
      if (npart_min*ppiclf_np .gt. npt_total) npart_min = npart_min-1

      if (npt_total .gt. PPICLF_LPART*ppiclf_np) 
     >   call ppiclf_exittr('Increase LPART to at least$',0.0d0
     >    ,npart_min)


      npmax = min(npt_total/PPICLF_LPART+1,ppiclf_np)
      stride_len = 0
      if (ppiclf_nid .le. npmax-1 .and. ppiclf_nid. ne. 0) 
     >stride_len = ppiclf_nid*PPICLF_LPART

      npart = PPICLF_LPART
      if (ppiclf_nid .gt. npmax-1) npart = 0

      ndiff = npt_total - (npmax-1)*PPICLF_LPART
      if (ppiclf_nid .eq. npmax-1) npart = ndiff

      iorank = -1

      call ppiclf_byte_open_mpi(trim(filein2),pth,.true.,ierr)

      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
      icount_pos = npart*3   
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_read_mpi(rout_pos,icount_pos,pth,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)
         j = 1 + (i-1)*npart
         icount_sln = npart

         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_read_mpi(rout_sln(j),icount_sln
     >                            ,pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         j = 1 + (i-1)*npart
         icount_lrp = npart

         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_read_mpi(rout_lrp(j),icount_lrp
     >                            ,pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)

         j = 1 + (i-1)*npart
         icount_lip = npart

         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_read_mpi(rout_lip(j),icount_lip
     >                            ,pth,ierr)
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      ic_pos = 0
      ic_sln = 0
      ic_lrp = 0
      ic_lip = 0
      do i=1,npart
         ic_pos = ic_pos + 1
         ppiclf_y(jx,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         ppiclf_y(jy,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
         ppiclf_y(jz,i) = rout_pos(ic_pos)
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,npart
         ic_sln = ic_sln + 1
         ppiclf_y(j,i) = rout_sln(ic_sln)
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,npart
         ic_lrp = ic_lrp + 1
         ppiclf_rprop(j,i) = rout_lrp(ic_lrp)
      enddo
      enddo
      do j=5,7
      do i=1,npart
         ic_lip = ic_lip + 1
         ppiclf_iprop(j,i) = int(rout_lip(ic_lip))
      enddo
      enddo

      ppiclf_npart = npart

      call ppiclf_printsi('  End ReadParticleVTU$',npt_total)

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadWallVTK(filein1)
     > bind(C, name="ppiclc_io_ReadWallVTK")
#else
      subroutine ppiclf_io_ReadWallVTK(filein1)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*8 points(3,4*PPICLF_LWALL)
      integer*4 fid, nmax, i, j, i1, i2, i3, isize, irsize, ierr,
     >          npoints, nwalls, indx1
      character*1000 text
      character*132 filein2
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      external PPICLF_CHSTR
!
      !! THROW ERRORS HERE IN FUTURE
      call ppiclf_prints(' *Begin ReadWallVTK$')
      
      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) k (3)
      filein2 = ppiclf_chstr(filein1(1:indx1))

      if (ppiclf_nid .eq. 0) then

      fid = 432

      open (unit=fid,file=trim(filein2),action="read")
      
      nmax = 10000
      do i=1,nmax
         read (fid,*,iostat=ierr) text,npoints

         do j=1,npoints
            read(fid,*) points(1,j),points(2,j),points(3,j)
         enddo

         read (fid,*,iostat=ierr) text,nwalls

         do j=1,nwalls
            if (ppiclf_ndim .eq. 2) then
               read(fid,*) i1,i2
    
               i1 = i1 + 1
               i2 = i2 + 1

               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1)/),
     >                 (/points(1,i2),points(2,i2)/),
     >                 (/points(1,i1),points(2,i1)/))  ! dummy 2d
    
            elseif (ppiclf_ndim .eq. 3) then
               read(fid,*) i1,i2,i3

               i1 = i1 + 1
               i2 = i2 + 1
               i3 = i3 + 1
    
               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1),points(3,i1)/),
     >                 (/points(1,i2),points(2,i2),points(3,i2)/),
     >                 (/points(1,i3),points(2,i3),points(3,i3)/))

            endif
         enddo
            
         exit
      enddo

      close(fid)

      endif

      isize  = 4
      irsize = 8
      call ppiclf_bcast(ppiclf_nwall, isize)
      call ppiclf_bcast(ppiclf_wall_c,9*PPICLF_LWALL*irsize)
      call ppiclf_bcast(ppiclf_wall_n,4*PPICLF_LWALL*irsize)

      call ppiclf_printsi('  End ReadWallVTK$',nwalls)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteSubBinVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos
      integer*8 stride_lenv
      integer*4 icount_pos(PPICLF_BX1, PPICLF_BY1, PPICLF_BZ1)
      real*4 rpoint(3)
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz,
     >          isize, ibin, jbin, kbin, ndumx, ndumy, i, j, k,
     >          ie, ndum, kmax, ii, jj, kk, itmp, jtmp, ktmp, il, ir,
     >          jl, jr, kl, kr, npa, npb, npc, npd, npe, npf, npg, nph,
     >          ioff_dum, itype, ivtu_size, if_pos, ierr, icount_dum,
     >          iorank
!

      call ppiclf_printsi(' *Begin WriteSubBinVTU$',ppiclf_cycle)

      call ppiclf_prints(' *Begin InitSolve$')
         call ppiclf_solve_InitSolve
      call ppiclf_prints('  End InitSolve$')

      call ppiclf_prints(' *Begin CreateSubBin$')
         call ppiclf_comm_CreateSubBin
      call ppiclf_prints('  End CreateSubBin$')

      call ppiclf_prints(' *Begin ProjectParticleSubBin$')
         call ppiclf_solve_ProjectParticleSubBin
      call ppiclf_prints('  End ProjectParticleSubBin$')

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      ncll_total = ppiclf_n_bins(1)*(ppiclf_bx-1)
     >            *ppiclf_n_bins(2)*(ppiclf_by-1)
      if (ppiclf_ndim.eq.3) ncll_total = ncll_total
     >            *ppiclf_n_bins(3)*(ppiclf_bz-1)
      nvtx_total = (ppiclf_n_bins(1)*(ppiclf_bx-1)+1)
     >            *(ppiclf_n_bins(2)*(ppiclf_by-1)+1)
      if (ppiclf_ndim.eq.3) nvtx_total = nvtx_total
     >            *(ppiclf_n_bins(3)*(ppiclf_bz-1)+1)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'grd'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

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
      write(vtu,'(A)',advance='no') 'type="ProjectParticleSubBinInt32" '
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
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
c     call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
c     call ppiclf_io_WriteDataArrayVTU(vtu,"NID",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
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
      ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
      ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
      ndum = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
      do ie=0,ndum-1
         i = modulo(ie,ppiclf_n_bins(1))
         j = modulo(ie/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = ie/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      kmax = 1
      if (ppiclf_ndim .eq. 3) kmax = ppiclf_bz-1

      do kk=1,kmax
      do jj=1,ppiclf_by-1
      do ii=1,ppiclf_bx-1

         itmp = i*(ppiclf_bx-1) + (ii-1)
         jtmp = j*(ppiclf_by-1) + (jj-1)
         ktmp = 0
         if (ppiclf_ndim .eq. 3)
     >   ktmp = k*(ppiclf_bz-1) + (kk-1)

         kl = ktmp
         kr = ktmp+1
         jl = jtmp
         jr = jtmp+1
         il = itmp
         ir = itmp+1

c        ndum = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         if (ppiclf_ndim .eq. 3) then
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
         else
            npa = il + ndumx*jl + ndumx*ndumy*kl
            npb = ir + ndumx*jl + ndumx*ndumy*kl
            npc = il + ndumx*jr + ndumx*ndumy*kl
            npd = ir + ndumx*jr + ndumx*ndumy*kl

            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
      enddo
      enddo
      enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      ! write offsetts here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
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
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ppiclf_ndim .eq. 3) then
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
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (k .ne. ppiclf_bz) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (k .eq. ppiclf_bz) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
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
         ! 2D
         else
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by) then
                  icount_pos(i,j,1) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif
            endif
         enddo
         enddo

         endif

      endif

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         stride_lenv = ppiclf_grid_i(i,j,k)
         idisp_pos   = ivtu_size + isize*(3*stride_lenv + 1)
         icount_dum  = icount_pos(i,j,k)
         rpoint(1)   = ppiclf_grid_x(i,j,k)
         rpoint(2)   = ppiclf_grid_y(i,j,k)
         rpoint(3)   = 0.0d0
         if (ppiclf_ndim .eq. 3)
     >   rpoint(3)   = ppiclf_grid_z(i,j,k)
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)

      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)


      ! projected fields
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
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos,idisp_cll,stride_lenv(8),stride_lenc
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz, ibin,
     >          jbin, kbin, il, ir, jl, jr, kl, kr, nbinpa, nbinpb,
     >          nbinpc, nbinpd, nbinpe, nbinpf, nbinpg, nbinph, i, j, k,
     >          ii, npa, npb, npc, npd, npe, npf, npg, nph, 
     >          ioff_dum, itype, iorank, if_cll, if_pos, icount_pos, 
     >          icount_cll, ierr, isize, ivtu_size
      real*4 rpoint(3)
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) then
          icalld1 = istartout
          if (istartout .eq. 0) icalld1 = 1
      end if

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      nvtx_total = (ppiclf_n_bins(1)+1)*(ppiclf_n_bins(2)+1)
      if (ppiclf_ndim .gt. 2) 
     >    nvtx_total = nvtx_total*(ppiclf_n_bins(3)+1)
      ncll_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)
      if (ppiclf_ndim .gt. 2) ncll_total = ncll_total*ppiclf_n_bins(3)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'bin'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin
      if (ppiclf_ndim .eq. 3) then
         kl = kbin
         kr = kbin+1
      endif

      nbinpa = il + ndxgpp1*jl + ndxygpp1*kl
      nbinpb = ir + ndxgpp1*jl + ndxygpp1*kl
      nbinpc = il + ndxgpp1*jr + ndxygpp1*kl
      nbinpd = ir + ndxgpp1*jr + ndxygpp1*kl
      if (ppiclf_ndim .eq. 3) then
         nbinpe = il + ndxgpp1*jl + ndxygpp1*kr
         nbinpf = ir + ndxgpp1*jl + ndxygpp1*kr
         nbinpg = il + ndxgpp1*jr + ndxygpp1*kr
         nbinph = ir + ndxgpp1*jr + ndxygpp1*kr
      endif

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         if (ppiclf_ndim .eq. 3) then
            stride_lenv(5) = nbinpe
            stride_lenv(6) = nbinpf
            stride_lenv(7) = nbinpg
            stride_lenv(8) = nbinph
         endif
         
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
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
         i = modulo(ii,ppiclf_n_bins(1))
         j = modulo(ii/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = 0
         if (ppiclf_ndim .eq. 3)
     >   k = ii/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
          
c     do K=0,ppiclf_n_bins(3)-1
         kl = K
         kr = K
         if (ppiclf_ndim .eq. 3) then
            kl = K
            kr = K+1
         endif
c     do J=0,ppiclf_n_bins(2)-1
         jl = J
         jr = J+1
c     do I=0,ppiclf_n_bins(1)-1
         il = I
         ir = I+1

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            npe = il + ndxgpp1*jl + ndxygpp1*kr
            npf = ir + ndxgpp1*jl + ndxygpp1*kr
            npg = il + ndxgpp1*jr + ndxygpp1*kr
            nph = ir + ndxgpp1*jr + ndxygpp1*kr
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
c     enddo
c     enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
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
     >    ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3*stride_lenv(1) + 1)
      rpoint(1)  = sngl(ppiclf_binx(1,1))
      rpoint(2)  = sngl(ppiclf_biny(1,1))
      rpoint(3)  = 0.0
      if (ppiclf_ndim .eq. 3)
     >rpoint(3)  = sngl(ppiclf_binz(1,1))
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (ppiclf_ndim .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = 0.0
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

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
      idisp_cll = ivtu_size + isize*(3*(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
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
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      character*3 filein
      character*12 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
     >          ierr, ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteParticleVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

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
      if_pos = 3*isize*npt_total
      if_sln = 1*isize*npt_total
      if_lrp = 1*isize*npt_total
      if_lip = 1*isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jy,i))
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
            rout_pos(ic_pos) = sngl(ppiclf_y(jz,i))
         else
            rout_pos(ic_pos) = 0.0
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,nxx
         ic_sln = ic_sln + 1
         rout_sln(ic_sln) = sngl(ppiclf_y(j,i))
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,nxx
         ic_lrp = ic_lrp + 1
         rout_lrp(ic_lrp) = sngl(ppiclf_rprop(j,i))
      enddo
      enddo
      do j=5,7
      do i=1,nxx
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(ppiclf_iprop(j,i))
      enddo
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
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
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
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do ie=1,PPICLF_LRS
         write(prostr,'(A1,I2.2)') "y",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,PPICLF_LRP
         write(prostr,'(A4,I2.2)') "rprop",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,3
         write(prostr,'(A3,I2.2)') "tag",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

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

      ! how much to write
      icount_pos = 3*nxx
      icount_sln = 1*nxx
      icount_lrp = 1*nxx
      icount_lip = 1*nxx

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

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_sln
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lrp
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)
         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lip
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

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
!
      implicit none
!
! Input:
!
      integer*4 vtu,ncomp
      integer*4 idist
      character (len = *) dataname
!
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
!
      implicit none
!
      include "PPICLF"
!
      call ppiclf_prints('*********** PPICLF OUTPUT *****************$')
      call ppiclf_io_OutputDiagGen
      call ppiclf_io_OutputDiagGhost
      if (ppiclf_lsubbin) call ppiclf_io_OutputDiagSubBin
      if (ppiclf_overlap) call ppiclf_io_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGen
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot, npart_ide, nbin_total
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

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
      call ppiclf_printsi('  -Number of bins (x)     :$',
     >                    ppiclf_n_bins(1))
      call ppiclf_printsi('  -Number of bins (y)     :$'
     >                    ,ppiclf_n_bins(2))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$'
     >                    ,ppiclf_n_bins(3))
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
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nel_max_orig, nel_min_orig, nel_total_orig, nel_max_map,
     >          nel_min_map, nel_total_map
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
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
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
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
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin_total
!

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
