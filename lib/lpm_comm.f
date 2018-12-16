!-----------------------------------------------------------------------
      subroutine lpm_comm_setup
#include "LPM"

c     call lpm_comm_interp_setup(i_fp_hndl,0.0,idum,lpm_nelt)
      call fgslib_crystal_setup(i_cr_hndl,lpm_comm,lpm_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_comm_interp_setup(ih,tolin,nmsh,nelm)
#include "LPM"

      ! need to add xm1, etc in LPM for this case

c     common /intp_h/ ih_intp(2,1)
c     common /intp/   tol

c     data ihcounter /0/
c     save ihcounter

c     real xmi, ymi, zmi
c     common /CBXMI/ xmi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT),
c    $               ymi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT),
c    $               zmi(LPM_LX1*LPM_LY1*LPM_LZ1*LPM_LELT)
c
c     real w(2*LPM_LX1*LPM_LY1*LPM_LZ1)
c     logical if3d

c     tol = max(5e-13,tolin)
c     npt_max = 128
c     bb_t    = 0.01

c     ! fix below
c     if3d = .true.

c     if (nmsh.ge.1 .and. nmsh.lt.LPM_LX1-1) then
c        nxi = nmsh+1
c        nyi = nxi
c        nzi = 1
c        if (if3d) nzi = nxi
c        do ie = 1,nelm
c          call map_m_to_n(xmi((ie-1)*nxi*nyi*nzi + 1),nxi,xm1(1,1,1,ie)
c    $                     ,lx1,if3d,w,size(w))
c          call map_m_to_n(ymi((ie-1)*nxi*nyi*nzi + 1),nyi,ym1(1,1,1,ie)
c    $                     ,ly1,if3d,w,size(w))
c          if (if3d) 
c    $     call map_m_to_n(zmi((ie-1)*nxi*nyi*nzi + 1),nzi,zm1(1,1,1,ie)
c    $                     ,lz1,if3d,w,size(w))
c        enddo

c        n = nelm*nxi*nyi*nzi 
c        call fgslib_findpts_setup(ih_intp1,nekcomm,npp,ldim,
c    &                             xm1,ym1,zm1,nx1,ny1,nz1,
c    &                             nelm,nx1,ny1,nz1,bb_t,nelm+2,nelm+2,
c    &                             npt_max,tol)

c        call fgslib_findpts_setup(ih_intp2,nekcomm,npp,ldim,
c    $                             xmi,ymi,zmi,nxi,nyi,nzi,
c    $                             nelm,2*nxi,2*nyi,2*nzi,bb_t,n,n,
c    $                             npt_max,tol)
c     else
c        n = nelm*nx1*ny1*nz1
c        call fgslib_findpts_setup(ih_intp1,nekcomm,npp,ldim,
c    &                             xm1,ym1,zm1,nx1,ny1,nz1,
c    &                             nelm,2*nx1,2*ny1,2*nz1,bb_t,n,n,
c    &                             npt_max,tol)
c        ih_intp2 = ih_intp1
c     endif

c     ihcounter = ihcounter + 1 
c     ih = ihcounter
c     if (ih .gt. INTP_HMAX)
c    $   call exitti('Maximum number of handles exceeded!$',INTP_HMAX)
c     ih_intp(1,ih) = ih_intp1
c     ih_intp(2,ih) = ih_intp2

      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_comm_bin_setup
#include "LPM"

      integer  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer  ifac(3), icount(3)
      real     d2new(3)
      logical  partl

      real                      lpm_xerange(2,3,lpm_lbmax)
      common /lpm_elementrange/ lpm_xerange

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (lpm_rparam(12) .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 3

      iperiodicx = int(lpm_rparam(8))
      iperiodicy = int(lpm_rparam(9))
      iperiodicz = int(lpm_rparam(10))
      ndim       = int(lpm_rparam(12))

      ! compute binb
      xmin = 1E8
      ymin = 1E8
      zmin = 1E8
      xmax = 0.
      ymax = 0.
      zmax = 0.
      do i=1,lpm_npart
         rduml = lpm_y(ix,i) - lpm_d2chk(2)
         rdumr = lpm_y(ix,i) + lpm_d2chk(2)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = lpm_y(iy,i) - lpm_d2chk(2)
         rdumr = lpm_y(iy,i) + lpm_d2chk(2)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         rduml = lpm_y(iz,i) - lpm_d2chk(2)
         rdumr = lpm_y(iz,i) + lpm_d2chk(2)
         if (rduml .lt. zmin) zmin = rduml
         if (rdumr .gt. zmax) zmax = rdumr
      enddo

      lpm_binb(1) = glmin(xmin,1)
      lpm_binb(2) = glmax(xmax,1)
      lpm_binb(3) = glmin(ymin,1)
      lpm_binb(4) = glmax(ymax,1)
      lpm_binb(5) = 0.0
      lpm_binb(6) = 0.0
      if(lpm_rparam(12) .gt. 2) lpm_binb(5) = glmin(zmin,1)
      if(lpm_rparam(12) .gt. 2) lpm_binb(6) = glmax(zmax,1)


      lpm_binb(1) = max(lpm_binb(1),lpm_xdrange(1,1))
      lpm_binb(2) = min(lpm_binb(2),lpm_xdrange(2,1))
      lpm_binb(3) = max(lpm_binb(3),lpm_xdrange(1,2))
      lpm_binb(4) = min(lpm_binb(4),lpm_xdrange(2,2))
      if(lpm_rparam(12) .gt. 2) lpm_binb(5) = 
     >                           max(lpm_binb(5),lpm_xdrange(1,3))
      if(lpm_rparam(12) .gt. 2) lpm_binb(6) = 
     >                           min(lpm_binb(6),lpm_xdrange(2,3))


      if (iperiodicx .eq. 0) then
         lpm_binb(1) = lpm_xdrange(1,1)
         lpm_binb(2) = lpm_xdrange(2,1)
      endif
      if (iperiodicy .eq. 0) then
         lpm_binb(3) = lpm_xdrange(1,2)
         lpm_binb(4) = lpm_xdrange(2,2)
      endif
      if (iperiodicz .eq. 0 .and. lpm_rparam(12) .gt. 2) then
         lpm_binb(5) = lpm_xdrange(1,3)
         lpm_binb(6) = lpm_xdrange(2,3)
      endif

      ifac(1) = 1
      ifac(2) = 1
      ifac(3) = 1
      icount(1) = 0
      icount(2) = 0
      icount(3) = 0
      d2new(1) = lpm_d2chk(2)
      d2new(2) = lpm_d2chk(2)
      d2new(3) = lpm_d2chk(2)

      lpm_ndxgp = floor( (lpm_binb(2) - lpm_binb(1))/d2new(1))
      lpm_ndygp = floor( (lpm_binb(4) - lpm_binb(3))/d2new(2))
      lpm_ndzgp = 1
      if (lpm_rparam(12) .gt. 2) lpm_ndzgp = 
     >                      floor( (lpm_binb(6) - lpm_binb(5))/d2new(3))


      if (lpm_ndxgp*lpm_ndygp*lpm_ndzgp .gt. lpm_np .or. 
     >    int(lpm_rparam(4)) .eq. 1) then
         nmax = 1000
         d2chk_save = lpm_d2chk(2)
         
         do i=1,nmax
         do j=0,ndim-1
            ifac(j+1) = 1 + i
            d2new(j+1) = (lpm_binb(2+2*j) - lpm_binb(1+2*j))/ifac(j+1)
            nbb = ifac(1)*ifac(2)*ifac(3)

            if( nbb .gt. lpm_np ) then
            if( int(lpm_rparam(4)) .eq. 1 .or.
     >          int(lpm_rparam(4)) .eq. 0 .and.d2new(j+1).lt.d2chk_save)
     >          then
               icount(j+1) = icount(j+1) + 1
               ifac(j+1) = ifac(j+1) - icount(j+1)
               d2new(j+1) = (lpm_binb(2+2*j) -lpm_binb(1+2*j))/ifac(j+1)
            endif
            endif
         enddo
            if (icount(1) .gt. 0) then
            if (icount(2) .gt. 0) then
            if (icount(3) .gt. 0) then
               exit
            endif
            endif
            endif
         enddo
      endif

! -------------------------------------------------------
c SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! how many spacings in each direction
      lpm_ndxgp = floor( (lpm_binb(2) - lpm_binb(1))/d2new(1))
      lpm_ndygp = floor( (lpm_binb(4) - lpm_binb(3))/d2new(2))
      lpm_ndzgp = 1
      if (lpm_rparam(12) .gt. 2) lpm_ndzgp = 
     >                      floor( (lpm_binb(6) - lpm_binb(5))/d2new(3))

      ! grid spacing for that many spacings
      lpm_rdxgp = (lpm_binb(2) - lpm_binb(1))/real(lpm_ndxgp)
      lpm_rdygp = (lpm_binb(4) - lpm_binb(3))/real(lpm_ndygp)
      lpm_rdzgp = 1.
      if (lpm_rparam(12) .gt. 2) lpm_rdzgp = (lpm_binb(6) - lpm_binb(5))
     >                           /real(lpm_ndzgp)

      ninc = 2
      rxlbin = lpm_binb(1)
      rxrbin = lpm_binb(2)
      rylbin = lpm_binb(3)
      ryrbin = lpm_binb(4)
      rzlbin = lpm_binb(5)
      rzrbin = lpm_binb(6)
      if (iperiodicx .ne. 0) then
         rxlbin = rxlbin - ninc/2*lpm_rdxgp
         rxrbin = rxrbin + ninc/2*lpm_rdxgp
         rxlbin = max(rxlbin,lpm_xdrange(1,1))
         rxrbin = min(rxrbin,lpm_xdrange(2,1))
      endif
      if (iperiodicy .ne. 0) then
         rylbin = rylbin - ninc/2*lpm_rdygp
         ryrbin = ryrbin + ninc/2*lpm_rdygp
         rylbin = max(rylbin,lpm_xdrange(1,2))
         ryrbin = min(ryrbin,lpm_xdrange(2,2))
      endif
      if (iperiodicz .ne. 0) then
      if (lpm_rparam(12) .gt. 2) then
         rzlbin = rzlbin - ninc/2*lpm_rdzgp
         rzrbin = rzrbin + ninc/2*lpm_rdzgp
         rzlbin = max(rzlbin,lpm_xdrange(1,3))
         rzrbin = min(rzrbin,lpm_xdrange(2,3))
      endif
      endif

      nbin_now = lpm_ndxgp*lpm_ndygp*lpm_ndzgp

      if (int(lpm_rparam(4)) .eq. 1) return ! only for projection

! see which bins are in which elements
c     lpm_neltb = 0
c     do ie=1,lpm_nelt
c     do k=1,lpm_lz1
c     do j=1,lpm_ly1
c     do i=1,lpm_lx1
c        rxval = xm1(i,j,k,ie) 
c        ryval = ym1(i,j,k,ie) 
c        rzval = 0.
c        if(if3d) rzval = zm1(i,j,k,ie)

c        if (rxval .gt. lpm_binb(2)) goto 1233
c        if (rxval .lt. lpm_binb(1)) goto 1233
c        if (ryval .gt. lpm_binb(4)) goto 1233
c        if (ryval .lt. lpm_binb(3)) goto 1233
c        if (if3d .and. rzval .gt. lpm_binb(6)) goto 1233
c        if (if3d .and. rzval .lt. lpm_binb(5)) goto 1233

c        ii    = floor((rxval-lpm_binb(1))/lpm_rdxgp) 
c        jj    = floor((ryval-lpm_binb(3))/lpm_rdygp) 
c        kk    = floor((rzval-lpm_binb(5))/lpm_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
c        if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
c        if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
c        ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk
c        nrank = modulo(ndum,np)

c        lpm_neltb = lpm_neltb + 1
c        if(lpm_neltb .gt. lpm_lbmax) then
c          write(6,*) 'increase lbmax',nid,lpm_neltb,lpm_lbmax
c          call exitt
c        endif

c        lpm_er_map(1,lpm_neltb) = ie
c        lpm_er_map(2,lpm_neltb) = nid
c        lpm_er_map(3,lpm_neltb) = ndum
c        lpm_er_map(4,lpm_neltb) = nrank
c        lpm_er_map(5,lpm_neltb) = nrank
c        lpm_er_map(6,lpm_neltb) = nrank

c        if (lpm_neltb .gt. 1) then
c        do il=1,lpm_neltb-1
c           if (lpm_er_map(1,il) .eq. ie) then
c           if (lpm_er_map(4,il) .eq. nrank) then
c              lpm_neltb = lpm_neltb - 1
c              goto 1233
c           endif
c           endif
c        enddo
c        endif
c1233 continue
c     enddo
c     enddo
c     enddo
c     enddo

c     nxyz = lx1*ly1*lz1
c     do ie=1,lpm_neltb
c        iee = lpm_er_map(1,ie)
c        call copy(lpm_xm1b(1,1,1,1,ie), xm1(1,1,1,iee),nxyz)
c        call copy(lpm_xm1b(1,1,1,2,ie), ym1(1,1,1,iee),nxyz)
c        call copy(lpm_xm1b(1,1,1,3,ie), zm1(1,1,1,iee),nxyz)
c     enddo

c     lpm_neltbb = lpm_neltb

c     do ie=1,lpm_neltbb
c        call icopy(lpm_er_maps(1,ie),lpm_er_map(1,ie),LPM_LRMAX)
c     enddo


c     nl   = 0
c     nii  = LPM_LRMAX
c     njj  = 6
c     nrr  = nxyz*3
c     nkey = 3
c     call fgslib_crystal_tuple_transfer(i_cr_hndl,lpm_neltb,lpm_lbmax
c    >                  , lpm_er_map,nii,partl,nl,lpm_xm1b,nrr,njj)
c     call fgslib_crystal_tuple_sort    (i_cr_hndl,lpm_neltb
c    $              , lpm_er_map,nii,partl,nl,lpm_xm1b,nrr,nkey,1)


c     do ie=1,lpm_neltb
c     do k=1,nz1
c     do j=1,ny1
c     do i=1,nx1
c        rxval = lpm_xm1b(i,j,k,1,ie)
c        ryval = lpm_xm1b(i,j,k,2,ie)
c        rzval = 0.
c        if(if3d) rzval = lpm_xm1b(i,j,k,3,ie)
c        
c        ii    = floor((rxval-lpm_binb(1))/lpm_rdxgp) 
c        jj    = floor((ryval-lpm_binb(3))/lpm_rdygp) 
c        kk    = floor((rzval-lpm_binb(5))/lpm_rdzgp) 
c        if (.not. if3d) kk = 0
c        if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
c        if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
c        if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
c        ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk

c        lpm_modgp(i,j,k,ie,1) = ii
c        lpm_modgp(i,j,k,ie,2) = jj
c        lpm_modgp(i,j,k,ie,3) = kk
c        lpm_modgp(i,j,k,ie,4) = ndum
c  
c     enddo
c     enddo
c     enddo
c     enddo

c     do ie=1,lpm_neltb
c        lpm_xerange(1,1,ie) = vlmin(lpm_xm1b(1,1,1,1,ie),nxyz)
c        lpm_xerange(2,1,ie) = vlmax(lpm_xm1b(1,1,1,1,ie),nxyz)
c        lpm_xerange(1,2,ie) = vlmin(lpm_xm1b(1,1,1,2,ie),nxyz)
c        lpm_xerange(2,2,ie) = vlmax(lpm_xm1b(1,1,1,2,ie),nxyz)
c        lpm_xerange(1,3,ie) = vlmin(lpm_xm1b(1,1,1,3,ie),nxyz)
c        lpm_xerange(2,3,ie) = vlmax(lpm_xm1b(1,1,1,3,ie),nxyz)

c        ilow  = floor((lpm_xerange(1,1,ie) - lpm_binb(1))/lpm_rdxgp)
c        ihigh = floor((lpm_xerange(2,1,ie) - lpm_binb(1))/lpm_rdxgp)
c        jlow  = floor((lpm_xerange(1,2,ie) - lpm_binb(3))/lpm_rdygp)
c        jhigh = floor((lpm_xerange(2,2,ie) - lpm_binb(3))/lpm_rdygp)
c        klow  = floor((lpm_xerange(1,3,ie) - lpm_binb(5))/lpm_rdzgp)
c        khigh = floor((lpm_xerange(2,3,ie) - lpm_binb(5))/lpm_rdzgp)
c        if (.not. if3d) then
c           klow = 0
c           khigh = 0
c        endif

c        lpm_el_map(1,ie) = ilow  + lpm_ndxgp*jlow  
c    >                            + lpm_ndxgp*lpm_ndygp*klow
c        lpm_el_map(2,ie) = ihigh + lpm_ndxgp*jhigh 
c    >                            + lpm_ndxgp*lpm_ndygp*khigh
c        lpm_el_map(3,ie) = ilow
c        lpm_el_map(4,ie) = ihigh
c        lpm_el_map(5,ie) = jlow
c        lpm_el_map(6,ie) = jhigh
c        lpm_el_map(7,ie) = klow
c        lpm_el_map(8,ie) = khigh
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine lpm_comm_findpts
#include "LPM"

c     common /intp_h/ ih_intp(2,1)

c     ih_intp2 = ih_intp(2,i_fp_hndl)

      ix = 1
      iy = 2
      iz = 3

c     call fgslib_findpts(ih_intp2           !   call fgslib_findpts( ihndl,
c    $        , lpm_iprop (1 ,1),LPM_LIP        !   $             rcode,1,
c    $        , lpm_iprop (3 ,1),LPM_LIP        !   &             proc,1,
c    $        , lpm_iprop (2 ,1),LPM_LIP        !   &             elid,1,
c    $        , lpm_rprop2(1 ,1),LPM_LRP2       !   &             rst,ndim,
c    $        , lpm_rprop2(4 ,1),LPM_LRP2       !   &             dist,1,
c    $        , lpm_y     (ix,1),LPM_LRS        !   &             pts(    1),1,
c    $        , lpm_y     (iy,1),LPM_LRS        !   &             pts(  n+1),1,
c    $        , lpm_y     (iz,1),LPM_LRS ,LPM_NPART) !   &             pts(2*n+1),1,n)

      do i=1,lpm_npart
         lpm_iprop(4,i) = lpm_iprop(3,i)
      enddo

      ! instead get which bin it is in
      do i=1,lpm_npart
         ! check if particles are greater or less than binb bounds....
         ! add below....

         ii    = floor((lpm_y(ix,i)-lpm_binb(1))/lpm_rdxgp) 
         jj    = floor((lpm_y(iy,i)-lpm_binb(3))/lpm_rdygp) 
         kk    = floor((lpm_y(iz,i)-lpm_binb(5))/lpm_rdzgp) 
         if (lpm_rparam(12) .lt. 3) kk = 0
         if (ii .eq. lpm_ndxgp) ii = lpm_ndxgp - 1
         if (jj .eq. lpm_ndygp) jj = lpm_ndygp - 1
         if (kk .eq. lpm_ndzgp) kk = lpm_ndzgp - 1
         ndum  = ii + lpm_ndxgp*jj + lpm_ndxgp*lpm_ndygp*kk
         nrank = modulo(ndum, lpm_np)

         lpm_iprop(8,i)  = ii
         lpm_iprop(9,i)  = jj
         lpm_iprop(10,i) = kk
         lpm_iprop(11,i) = ndum

         lpm_iprop(4,i)  = nrank ! where particle is actually moved
      enddo


      return
      end
!-----------------------------------------------------------------------
      subroutine lpm_comm_crystal
#include "LPM"

      logical partl    
      integer lpm_ipmap1(1,LPM_LPART)
     >       ,lpm_ipmap2(1,LPM_LPART)
     >       ,lpm_ipmap3(1,LPM_LPART)
     >       ,lpm_ipmap4(1,LPM_LPART)

      parameter(lrf = LPM_LRS*4 + LPM_LRP + LPM_LRP2)
      real rwork(lrf,LPM_LPART)

      do i=1,lpm_npart
         ic = 1
         call copy(rwork(ic,i),lpm_y(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_y1((i-1)*LPM_LRS+1),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_ydot(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_ydotc(1,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(rwork(ic,i),lpm_rprop(1,i),LPM_LRP)
         ic = ic + LPM_LRP
         call copy(rwork(ic,i),lpm_rprop2(1,i),LPM_LRP2)
      enddo

      j0 = 4
      call fgslib_crystal_tuple_transfer(i_cr_hndl,lpm_npart ,LPM_LPART
     $           ,lpm_iprop ,LPM_LIP,partl,0,rwork,lrf ,j0)

      do i=1,lpm_npart
         ic = 1
         call copy(lpm_y(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_y1((i-1)*LPM_LRS+1),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_ydot(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_ydotc(1,i),rwork(ic,i),LPM_LRS)
         ic = ic + LPM_LRS
         call copy(lpm_rprop(1,i),rwork(ic,i),LPM_LRP)
         ic = ic + LPM_LRP
         call copy(lpm_rprop2(1,i),rwork(ic,i),LPM_LRP2)
      enddo
        
      return
      end
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
