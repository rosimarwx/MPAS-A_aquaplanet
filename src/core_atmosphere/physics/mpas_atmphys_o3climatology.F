! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_o3climatology
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_atmphys_date_time
 use mpas_atmphys_constants
 use mpas_atmphys_utilities

!wrf physics:
 use module_ra_cam_support, only: r8, getfactors

 implicit none
 private
 public:: init_o3climatology,          &
          update_o3climatology,        &
          o3climatology_from_MPAS,     &
          init_o3_ape !added by RRB to use hemispheric symmetric ozone

 integer,parameter:: latsiz = 64
 integer,parameter:: lonsiz = 1

!mpas_atmphys_o3climatology contains the subroutines needed to initialize,interpolate,and update
!the climatological monthly-mean ozone volume mixing ratios o3clim to the MPAS grid. Input data
!files are the same as the ones used in the CAM long- and short-wave radiation codes.
!when namelist parameter config_o3climatology is set to true, o3clim is used in the RRTMG long-
!wave and short-wave radiation codes,and replaces the annual-mean ozone sounding used by default.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-07-03.
!
! subroutines in mpas_atmphys_o3climatology:
! ------------------------------------------
! init_o3climatology     : read the CAM ozone data files.
! update_o3climatology   : interpolates the ozone volume mixing ratio to the current Julian day
!                          as done for the greeness fraction in the MPAS time manager.
! o3climatology_from_MPAS: interpolates the ozone volume mixing ratio to the current Julian day
!                          as in the CAM radiation codes.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the subroutine vinterp_ozn to its own module module_ra_rrtmg_vinterp.F in physics_wrf.
!   Laura D. Fowler (laura@ucar.edu) / 2017-01-27.


 contains


!=================================================================================================================
 subroutine init_o3climatology(mesh,atm_input)
!=================================================================================================================

!This subroutine assumes a uniform distribution of ozone concentration. It should be replaced
!with monthly climatology varying ozone distribution.

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: atm_input

!local pointers:
 integer, pointer:: nCells,num_months,levsiz
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: pin
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

!local variables:
 integer,parameter:: pin_unit = 27
 integer,parameter:: lat_unit = 28
 integer,parameter:: oz_unit  = 29
 integer,parameter:: open_ok  = 0

 integer:: i,i1,i2,istat,k,j,m
 integer:: iCell
 
 real(kind=RKIND):: lat,lon,dlat,dlatCell
 real(kind=RKIND),dimension(latsiz):: lat_ozone
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: ozmixin

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine physics_init_o3:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nMonths',num_months)
 call mpas_pool_get_dimension(mesh,'nOznLevels',levsiz)

 call mpas_pool_get_array(atm_input,'pin',pin)
 call mpas_pool_get_array(atm_input,'ozmixm',ozmixm)
 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

!-- read in ozone pressure data:
 open(pin_unit,file='OZONE_PLEV.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_PLEV.TBL')
 do k = 1,levsiz
    read(pin_unit,*) pin(k)
 enddo
 close(pin_unit)

!-- read in ozone lat data:
 open(lat_unit, file='OZONE_LAT.TBL',action='READ',status='OLD',iostat=istat) 
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                             'failure opening OZONE_LAT.TBL')
 do j = 1, latsiz
    read(lat_unit,*) lat_ozone(j)
!   call mpas_log_write('$i $r', intArgs=(/j/), realArgs=(/lat_ozone(j)/))
 enddo
 close(lat_unit)

!-- read in ozone data:
 open(oz_unit,file='OZONE_DAT.TBL',action='READ',status='OLD',iostat=istat)
 if(istat /= open_ok) &
    call physics_error_fatal('subroutine oznini: ' // &
                                'failure opening OZONE_DAT.TBL')

 allocate(ozmixin(lonsiz,levsiz,latsiz,num_months))
 do m=1,num_months
 do j=1,latsiz ! latsiz=64
 do k=1,levsiz ! levsiz=59
 do i=1,lonsiz ! lonsiz=1
    read(oz_unit,*) ozmixin(i,k,j,m)
 enddo
 enddo
 enddo
 enddo
 close(oz_unit)

!INTERPOLATION OF INPUT OZONE DATA TO MPAS GRID:
!call mpas_log_write('max latCell=$r', realArgs=(/maxval(latCell)/degrad/))
!call mpas_log_write('min latCell=$r', realArgs=(/minval(latCell)/degrad/))
!call mpas_log_write('max lonCell=$r', realArgs=(/maxval(lonCell)/degrad/))
!call mpas_log_write('min lonCell=$r', realArgs=(/minval(lonCell)/degrad/))
!call mpas_log_write('')
!call mpas_log_write('max lat_ozone=$r', realArgs=(/maxval(lat_ozone)/))
!call mpas_log_write('min lat_ozone=$r', realArgs=(/minval(lat_ozone)/))
 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    lon = lonCell(iCell)/degrad
    if(lat .gt. lat_ozone(latsiz)) then
     i1 = latsiz
     i2 = latsiz
    elseif(lat .lt. lat_ozone(1)) then
       i1 = 1
       i2 = 1
   else
       do i = 1, latsiz
          if(lat.ge.lat_ozone(i) .and. lat.lt.lat_ozone(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do m = 1,num_months
    do k = 1,levsiz
    do j = 1,lonsiz
       dlat     = lat_ozone(i2)-lat_ozone(i1)
       dlatCell = lat-lat_ozone(i1)
       if(dlat == 0.) then
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m)
       else
          ozmixm(m,k,iCell) = ozmixin(j,k,i1,m) &
                     + (ozmixin(j,k,i2,m)-ozmixin(j,k,i1,m))*dlatCell/dlat
       endif
    enddo 
    enddo       
    enddo
!   do k = 1, levsiz
!      call mpas_log_write('$i $i $i $r $r $r $r $r $r', intArgs=(/iCell,i1,i2/), &
!                          realArgs=(/lat_ozone(i1),lat,lat_ozone(i2),ozmixin(1,k,i1,1),ozmixm(1,k,iCell),ozmixin(1,k,i2,1)/))
!   enddo
 enddo
 deallocate(ozmixin)

! call mpas_log_write('--- end subroutine physics_init_o3.')

 end subroutine init_o3climatology

!=================================================================================================================
 subroutine update_o3climatology(current_date,mesh,atm_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in) :: mesh
 character(len=*),intent(in):: current_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer, pointer:: nCellsSolve,nOznLevels
 real(kind=RKIND),dimension(:,:),pointer  :: o3clim
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

!local variables:
 integer:: iCell,iLev

!-----------------------------------------------------------------------------------------------------------------
! call mpas_log_write('')
! call mpas_log_write('--- enter subroutine physics_update_o3:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nOznLevels',nOznLevels)

 call mpas_pool_get_array(atm_input,'ozmixm',ozmixm)
 call mpas_pool_get_array(diag_physics,'o3clim',o3clim)

 do iLev = 1,nOznLevels 
    call monthly_interp_to_date(nCellsSolve,current_date,ozmixm(:,iLev,:),o3clim(iLev,:))
 enddo
! call mpas_log_write('--- end subroutine physics_update_o3:')

 end subroutine update_o3climatology

!=================================================================================================================
 subroutine o3climatology_from_MPAS(julian,mesh,atm_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in):: julian
 type(mpas_pool_type),intent(in):: atm_input 

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,pointer:: nOznLevels,nVertLevels,nCellsSolve,nMonths
 real(kind=RKIND),dimension(:,:),pointer  :: o3clim
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmixm

!local variables:
 logical:: finddate
 logical:: ozncyc

 integer:: iCell,k
 integer:: ijul,m,nm,np,np1
 integer, dimension(12) ::  date_oz
 data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/

 real(kind=r8):: intjulian_r8
 real(kind=r8):: cdayozp_r8, cdayozm_r8
 real(kind=r8):: fact1_r8, fact2_r8
 real(kind=RKIND):: fact1,fact2

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine o3climatology_from_MPAS:')

 call mpas_pool_get_dimension(mesh,'nOznLevels',nOznLevels)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)
 
 call mpas_pool_get_array(atm_input,'ozmixm',ozmixm)
 call mpas_pool_get_array(diag_physics,'o3clim',o3clim)

 ozncyc = .true.

!julian starts from 0.0 at 0Z on 1 Jan.
 intjulian_r8 = real(julian + 1.0_RKIND, r8)    ! offset by one day

!jan 1st 00z is julian=1.0 here
 ijul=int(intjulian_r8)
!note that following will drift. need to use actual month/day info to compute julian.
 intjulian_r8 = intjulian_r8 - real(ijul, r8)
 ijul = mod(ijul,365)
 if(ijul .eq. 0) ijul=365
 intjulian_r8 = intjulian_r8 + real(ijul, r8)
 np1=1

 finddate=.false.
 do m = 1, nMonths
    if(date_oz(m).gt.intjulian_r8 .and. .not.finddate) then
       np1 = m
       finddate = .true.
    endif
 enddo
 cdayozp_r8=date_oz(np1)

 if(np1 .gt. 1) then
    cdayozm_r8 = date_oz(np1-1)
    np = np1
    nm = np-1
 else
    cdayozm_r8 = date_oz(12)
    np = np1
    nm = nMonths
 endif
 call getfactors(ozncyc,np1,cdayozm_r8,cdayozp_r8,intjulian_r8,fact1_r8,fact2_r8) 
 fact1 = real(fact1_r8, RKIND)
 fact2 = real(fact2_r8, RKIND)
 
!call mpas_log_write('--- end subroutine getfactors:')
!call mpas_log_write('fact1 =$r', realArgs=(/fact1/))
!call mpas_log_write('fact2 =$r', realArgs=(/fact2/))

!Time interpolation.
 do k = 1, nOznLevels
 do iCell = 1, nCellsSolve
    o3clim(k,iCell) = fact1 * ozmixm(nm,k,iCell) + fact2 * ozmixm(np,k,iCell)
 end do
!call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/k/), realArgs=(/o3clim(k,1),ozmixm(nm,k,1),ozmixm(np,k,1), &
!             o3clim(k,nCellsSolve),ozmixm(nm,k,nCellsSolve),ozmixm(np,k,nCellsSolve)/))
 end do

!call mpas_log_write('--- end subroutine o3climatology_from_MPAS')

 end subroutine o3climatology_from_MPAS

!!!! The following subroutine was added by RRB to add zonally symmetric ozone
!=================================================================================================================
 subroutine init_o3_ape(mesh,atm_input,diag_physics)
!=================================================================================================================

!This subroutine reads in a zonally symmetric ozone distribution
!for aquaplanet simulations

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer, pointer:: nCells,levsiz
 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: pin
 real(kind=RKIND),dimension(:,:),pointer:: o3clim

!local variables:
 integer,parameter:: open_ok  = 0

 integer:: i,i1,i2,istat,k,j,m
 integer:: iCell
 
 real(kind=RKIND):: lat,lon,dlat,dlatCell
 real(kind=RKIND),dimension(latsiz):: lat_ozone
 real(kind=RKIND),dimension(:,:),allocatable:: ozmixin

!Variables to read in ozone
 character finame*80,label*120
 integer :: lato3d  != 64  ! number of data latitudinal grid
 integer :: layo3d  != 59  ! number of data vertical layers
 integer :: lvlo3d  != 60  ! number of data vertical layer interfaces
! real, allocatable :: glat(:)
 real, allocatable :: gprs(:), galt(:)
 real, allocatable :: gpri(:), gali(:)
 real, allocatable :: gdat(:,:)
 real :: o3vubc, prsubc, altubc

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine physics_init_o3:'

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nOznLevels',levsiz)

 call mpas_pool_get_array(atm_input,'pin',pin)
 call mpas_pool_get_array(diag_physics,'o3clim',o3clim)
 call mpas_pool_get_array(mesh,'latCell',latCell)

 !!!!Read in ozone data. File name: apeo3.dat, source: http://climate.ncas.ac.uk/ape/ape_ozone.html
      
 open (unit=45,file='apeo3.dat',status='old')
 read (45,'(i3)') lato3d

 read (45,'(i3)') layo3d
 lvlo3d = 1 + layo3d

 !File has data we don't need...
 read (45,'(3(1pe12.5))') o3vubc,prsubc,altubc

 allocate(galt(layo3d))
 allocate(gpri(lvlo3d))
 allocate(gali(lvlo3d))
 allocate(ozmixin(latsiz,levsiz))


 write(0,*) 'Reading ozone data...'
 read (45,'(a)') label   ! grid latitudes (deg) --------->
 read (45,'(10(1pe12.5))') (lat_ozone(j),j=1,lato3d)
 read (45,'(a)') label   ! layer pressure (mb) ---------->
 read (45,'(10(1pe12.5))') (pin(k),k=1,layo3d)

 read (45,'(a)') label   ! layer altitude (km) ---------->
 read (45,'(10(1pe12.5))') (galt(k),k=1,layo3d)

 read (45,'(a)') label   ! interface pressure (mb) ------>
 read (45,'(10(1pe12.5))') (gpri(k),k=1,lvlo3d)

 read (45,'(a)') label   ! interface altitude (km) ------>
 read (45,'(10(1pe12.5))') (gali(k),k=1,lvlo3d)

 read (45,'(a)') label   ! o3 vmr (ppmv) lat-alt
 read (45,'(10(1pe12.5))') ((ozmixin(j,k),j=1,latsiz),k=1,levsiz)
 close(45)

 !Convert from ppmv to mol/mol; 1 ppmv = 1e-6 mol/mol
 ozmixin = ozmixin*1.e-6

! write(0,*) gdat


!INTERPOLATION OF INPUT OZONE DATA TO MPAS GRID:
!write(0,*) 'max latCell=', maxval(latCell)/degrad
!write(0,*) 'min latCell=', minval(latCell)/degrad
!write(0,*) 'max lonCell=', maxval(lonCell)/degrad
!write(0,*) 'min lonCell=', minval(lonCell)/degrad
!write(0,*)
!write(0,*) 'max lat_ozone=',maxval(lat_ozone)
!write(0,*) 'min lat_ozone=',minval(lat_ozone)
 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    if(lat .gt. lat_ozone(latsiz)) then
     i1 = latsiz
     i2 = latsiz
    elseif(lat .lt. lat_ozone(1)) then
       i1 = 1
       i2 = 1
   else
       do i = 1, latsiz
          if(lat.ge.lat_ozone(i) .and. lat.lt.lat_ozone(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do k = 1,levsiz
       dlat     = lat_ozone(i2)-lat_ozone(i1)
       dlatCell = lat-lat_ozone(i1)
       if(dlat == 0.) then
          o3clim(k,iCell) = ozmixin(i1,k)
       else
          o3clim(k,iCell) = ozmixin(i1,k) &
                     + (ozmixin(i2,k)-ozmixin(i1,k))*dlatCell/dlat
       endif
    enddo 
!   do k = 1, levsiz
!      write(0,102) iCell,i1,i2,lat_ozone(i1),lat,lat_ozone(i2),ozmixin(1,k,i1,1), &
!                   ozmixm(1,k,iCell),ozmixin(1,k,i2,1)
!   enddo
 enddo
 deallocate(ozmixin)

! write(0,*) '--- end subroutine init_o3_ape.'

!formats:
 101 format(i3,12(1x,e15.8))
 102 format(i6,i6,i6,6(1x,e15.8))

 end subroutine init_o3_ape

!=================================================================================================================
 end module mpas_atmphys_o3climatology
!=================================================================================================================
