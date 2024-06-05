module output_routes
use constants
use basin_grid
use ocean_variables
use ocean_bc
use iodata_routes
use rw_ctl_routes
use ocalg_routes

implicit none

include 'reclen.fi'

contains
!-------------------------------------------------------------------------
subroutine local_output(path2data,  &
                        nrec,       &
                        year,       &
                       month,       & 
                         day,       &
                        hour,       &
                      minute,       &
                      tstep,        &
                      calendar  )

include 'locout.fi'

integer,intent(in):: nrec, year, month, day, hour, minute, calendar
character*(*),intent(in)::  path2data
real(4),intent(in)::  tstep

character fname*256
integer m,n,k,ierr
real(8) z8(nz), zw8(nz+1), z0(1), z1(1), zlev8(nlev_loc)

z8 = dble(z )*1000.0d0
zw8= dble(zw)*1000.0d0
z0 = 0.0d0
z1 = 1.0d0
zlev8=dble(zlev_loc)

if (rank == 0) then
  write(*,*) 'Writing local output (instantaneous on time), record number ', nrec
endif

!$omp parallel do private(m,n,k)
do n=bnd_y1, bnd_y2
 do m=bnd_x1, bnd_x2 
  
  aux_array2d_01(m,n) = 0.0
  
  do k=1, nz 
   aux_array3d_tgr1(m,n,k)=0.0
  enddo

  do k=1, nlev_loc
   aux_array3d_zgr_loc(m,n,k)=0.0
  enddo

 enddo
enddo
!$omp end parallel do

if(nrec==1) then
  !writing HHQ
  ierr=0
  call pwdstd(path2data,'LOCAL/hhq.dat',nrec,hhq_rest,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'HHQ is written'
  
  call fulfname(fname,path2data,'LOCAL/hhq.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                      'HHQ, m',    &     !title of dataset
                          'hhq'   )      !variable name
  endif

  !writing geographical longitude on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=geo_lon_t(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/geolon_t.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Geolon is written'

  call fulfname(fname,path2data,'LOCAL/geolon_t.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Geographical longitude on t-grid, degree',    &     !title of dataset
                          'glon'   )      !variable name
  endif

  !writing geographical latitude on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=geo_lat_t(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/geolat_t.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Geolat is written'

  call fulfname(fname,path2data,'LOCAL/geolat_t.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Geographical latitude on t-grid, degree',    &     !title of dataset
                          'glat'   )      !variable name
  endif

!writing cosine of angles between model and geographical grids on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=rotvec_coeff(m,n,1)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/cosrot.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Cosrot is written'

  call fulfname(fname,path2data,'LOCAL/cosrot.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Cosine of angle between model and geographical grid lines',    &     !title of dataset
                          'cosrot'   )      !variable name
  endif

!writing sine of angles between model and geographical grids on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=rotvec_coeff(m,n,2)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/sinrot.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Sinrot is written'

  call fulfname(fname,path2data,'LOCAL/sinrot.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Sine of angle between model and geographical grid lines',    &     !title of dataset
                          'sinrot'   )      !variable name
  endif
     
!writing DX grid step

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=dx(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/dx.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'DX is written'

  call fulfname(fname,path2data,'LOCAL/dx.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
              'X grid step, m',    &     !title of dataset
                          'dx'   )      !variable name
  endif

!writing DY grid step

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=dy(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/dy.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'DY is written'

  call fulfname(fname,path2data,'LOCAL/dy.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
              'Y grid step, m',    &     !title of dataset
                          'dy'   )      !variable name
  endif

endif

if(ssh_output_loc>0) then

!writing SSH
!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=ssh(m,n)+mistot(m,n)/RefDen*float(variable_volume_budget)
   enddo
  enddo
!$omp end parallel do
  
  ierr=0
  
  call pwdstd(path2data,'LOCAL/ssh.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSH is written'

  call fulfname(fname,path2data,'LOCAL/ssh.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                      'SSH, m',    &     !title of dataset
                          'ssh'   )      !variable name
  endif

!writing UB
!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     aux_array2d_01(m,n)=(ubrtr(m,n)*dyh(m,n)+ubrtr(m-1,n)*dyh(m-1,n))/(2.0*dy(m,n)*hhq(m,n))
    endif
   enddo
  enddo
!$omp end parallel do
  
  ierr=0
  
  call pwdstd(path2data,'LOCAL/ub.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'UB is written'

  call fulfname(fname,path2data,'LOCAL/ub.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                     'UB, m/s',    &     !title of dataset
                            'u'   )      !variable name
  endif

!writing VB
!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     aux_array2d_01(m,n)=(vbrtr(m,n)*dxh(m,n)+vbrtr(m,n-1)*dxh(m,n-1))/(2.0*dx(m,n)*hhq(m,n))
    endif
   enddo
  enddo
!$omp end parallel do
  
  ierr=0
  
  call pwdstd(path2data,'LOCAL/vb.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'VB is written'

  call fulfname(fname,path2data,'LOCAL/vb.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                     'VB, m/s',    &     !title of dataset
                            'v'   )      !variable name
  endif

endif

if(mld_output_loc>0) then
!mixed layer depths writing

  ierr=0
  call pwdstd(path2data,'LOCAL/mld_dens.dat',nrec,mld_dens,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'MLD-DENS is written'

  call fulfname(fname,path2data,'LOCAL/mld_dens.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
 'Mixed layer depth based on density difference, m',    &     !title of dataset
                          'mld'   )      !variable name
  endif

endif

if(sls_output_loc>0) then
!mixed layer depths writing

  ierr=0
  call pwdstd(path2data,'LOCAL/sls.dat',nrec,sls,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SLS is written'

  call fulfname(fname,path2data,'LOCAL/sls.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
         'Steric sea level, m',    &     !title of dataset
                          'sls'   )      !variable name
  endif

endif

if(uv_output_loc>0) then 
!-----------------------------------------------------------------------------------------------------
  if(grid_shift_loc==0) then !writing on the model grid
  
  ierr=0
  !writing zonal velocity
    call pwdstd(path2data,'LOCAL/uu.dat',nrec,uu,llu,nx,ny,nz,m1loc-1,m2loc,n1loc,n2loc,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'U at u-grid s-levels is written'

    call fulfname(fname,path2data,'LOCAL/uu.dat',ierr)
    
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                        nx_loc+1,    &     !x-dimension
                          ny_loc,    &     !y-dimension
                              nz,    &     !z-dimension
                            nrec,    &     !t-dimension
                        xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xu(m1loc-1:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                              1,     &     !z-grid type (0 - linear, 1 - levels)
                              z8,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
            'zonal velocity, m/s',   &     !title of dataset
                              'u'   )      !variable name
    endif
  
  ierr=0
  !writing meridional velocity
    call pwdstd(path2data,'LOCAL/vv.dat',nrec,vv,llv,nx,ny,nz,m1loc,m2loc,n1loc-1,n2loc,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'V at v-grid s-levels is written'

    call fulfname(fname,path2data,'LOCAL/vv.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                        ny_loc+1,   &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),    &     !first x-value (if linear) or x-array (if levels)
                          dxst,     &     !x-step (if linear)
                      ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
              yv(n1loc-1:n2loc),    &     !first y-value (if linear) or x-array (if levels)
                          dyst,     &     !y-step (if linear)
                          1,        &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
      'meridional velocity, m/s',   &     !title of dataset
                            'v'    )      !variable name
    endif

  else !writing on T-grid
  
  !writing zonal velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
          do k=1,nz
            aux_array3d_tgr1(m,n,k)= (uu(m  ,n,k)*dyh(m  ,n)*hhu(m  ,n)   &
                                     +uu(m-1,n,k)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) 
          enddo
        endif
      enddo
    enddo
    !$omp end parallel do
  
  if(sig_output_loc>0) then    !writing at s-levels
   
    ierr=0

    call pwdstd(path2data,'LOCAL/uu.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'U at t-grid s-levels is written'

    call fulfname(fname,path2data,'LOCAL/uu.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                          ny_loc,   &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                       xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                 xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
'zonal velocity at s-levels, m/s',  &     !title of dataset
                            'u'    )      !variable name
    endif
  
  endif !writing at s-levels

  if(z_output_loc>0) then !writing at z-levels
   
    call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
             aux_array3d_zgr_loc,   &       !output 3d field (on z-levels)
                     hhq_rest,      &       !2d field of bottom topography (in metres),
                           lu,      &       !temperature mask
                            z,      &       !array of s-levels
                     zlev_loc,      &       !array of z-levels (in metres)
                bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
                bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                           nz,      &       !number of s-levels
                     nlev_loc,      &       !number of z-levels
                            0,      &       !parameter of task
                       undef) 
    ierr=0

    call pwdstd(path2data,'LOCAL/uz.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'U at t-grid z-levels is written'

    call fulfname(fname,path2data,'LOCAL/uz.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                          ny_loc,   &     !y-dimension
                        nlev_loc,   &     !z-dimension
                            nrec,   &     !t-dimension
                       xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                 xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                           zlev8,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
'zonal velocity at z-levels, m/s',  &     !title of dataset
                            'u'    )      !variable name
           
    endif

  endif !writing at z-levels

  !writing meridional velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
          do k=1,nz
           aux_array3d_tgr1(m,n,k)= (vv(m,n  ,k)*dxh(m,n  )*hhv(m,n  )    &
                                    +vv(m,n-1,k)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n)
          enddo
        endif
      enddo
    enddo
    !$omp end parallel do

  if(sig_output_loc>0) then    !writing at s-levels
    
    ierr=0

    call pwdstd(path2data,'LOCAL/vv.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'V at t-grid s-levels is written'

    call fulfname(fname,path2data,'LOCAL/vv.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                          nx_loc,    &     !x-dimension
                          ny_loc,    &     !y-dimension
                              nz,    &     !z-dimension
                            nrec,    &     !t-dimension
                        xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                                1,   &     !z-grid type (0 - linear, 1 - levels)
                              z8,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
'meridional velocity at s-levels, m/s',  &     !title of dataset
                              'v'   )      !variable name
    endif
  
  endif !writing at s-levels

  if(z_output_loc>0) then !writing at z-levels
   
    call s2z(aux_array3d_tgr1   ,   &       !input  3d field (on s-levels)
             aux_array3d_zgr_loc,   &       !output 3d field (on z-levels)
                        hhq_rest,   &       !2d field of bottom topography (in metres),
                              lu,   &       !temperature mask
                               z,   &       !array of s-levels
                        zlev_loc,   &       !array of z-levels (in metres)
                bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
                bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                           nz,      &       !number of s-levels
                     nlev_loc,      &       !number of z-levels
                            0,      &       !parameter of task
                          undef) 
    ierr=0

    call pwdstd(path2data,'LOCAL/vz.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'V at t-grid z-levels is written'

    call fulfname(fname,path2data,'LOCAL/vz.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                          ny_loc,   &     !y-dimension
                        nlev_loc,   &     !z-dimension
                            nrec,   &     !t-dimension
                       xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                 xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                           zlev8,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
'meridional velocity at z-levels, m/s',  &     !title of dataset
                              'v'   )      !variable name
           
    endif

  endif !writing at z-levels

  endif !gridtype

endif !write or not UV

!--------------------------------------------------------------------------------
if(ts_output_loc>0) then
  !writing temperature
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/tt.dat',nrec,tt,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'T at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/tt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
'Temperature at s-levels, °C',    &     !title of dataset
                          'tt'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(       tt,      &       !input  3d field (on s-levels)
  aux_array3d_zgr_loc,      &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
             zlev_loc,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
             nlev_loc,      &       !number of z-levels
                    0,      &       !parameter of task
               undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/tz.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'T at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/tz.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
 'Temperature at z-levels, °C',    &     !title of dataset
                           'tt'   )      !variable name
          
   endif
 
  endif !writing at z-levels

  !writing salinity

  aux_array3d_tgr1=ss+salref

 if(sig_output_loc>0) then    !writing at s-levels
 
  ierr=0
  call pwdstd(path2data,'LOCAL/ss.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'S at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/ss.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
              'Salinity, PSU',    &     !title of dataset
                          'ss'   )      !variable name
  endif

 endif !writing at s-levels
 
 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
        aux_array3d_zgr_loc,       &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
             zlev_loc,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
             nlev_loc,      &       !number of z-levels
                    0,      &       !parameter of task
               undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/sz.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'S at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/sz.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
   'Salinity at z-levels, PSU',    &     !title of dataset
                           'ss'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!--------------------------------------------------------------------------------
if(den_sgt_output_loc>0) then
  !writing in situ density
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/den_sgt.dat',nrec,den_sgt,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-t at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/den_sgt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    'Density sigma-t, kg/m^3',    &     !title of dataset
                          'den'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(  den_sgt,      &       !input  3d field (on s-levels)
  aux_array3d_zgr_loc,      &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
             zlev_loc,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
             nlev_loc,      &       !number of z-levels
                    0,      &       !parameter of task
               undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/den_sgt_z.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-t at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/den_sgt_z.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
     'Density sigma-t, kg/m^3',    &     !title of dataset
                           'den'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!--------------------------------------------------------------------------------
if(den_sg0_output_loc>0) then
  !writing potential density (p=0)
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/den_sg0.dat',nrec,den_sg0,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-0 at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/den_sg0.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    'Density sigma-0, kg/m^3',    &     !title of dataset
                          'den'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(  den_sg0,      &       !input  3d field (on s-levels)
  aux_array3d_zgr_loc,      &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
             zlev_loc,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
             nlev_loc,      &       !number of z-levels
                    0,      &       !parameter of task
               undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/den_sg0_z.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-0 at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/den_sg0_z.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
     'Density sigma-0, kg/m^3',    &     !title of dataset
                           'den'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!--------------------------------------------------------------------------------
if(den_sgm_output_loc>0) then
  !writing potential density (p=p(zmean))
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/den_sgm.dat',nrec,den_sgm,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-zmean/2 at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/den_sgm.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
'Density sigma-zmean/2, kg/m^3',    &     !title of dataset
                          'den'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(  den_sgm,      &       !input  3d field (on s-levels)
  aux_array3d_zgr_loc,      &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
             zlev_loc,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
             nlev_loc,      &       !number of z-levels
                    0,      &       !parameter of task
               undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/den_sgm_z.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'Density potential sigma-zmean/2 at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/den_sgm_z.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
 'Density sigma-zmean/2, kg/m^3',    &     !title of dataset
                           'den'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!--------------------------------------------------------------------------------
if(pt_output_loc>0) then
  !writing passive tracer
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/pt.dat',nrec,pass_tracer,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Passive tracer at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/pt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
 'Passive tracer at s-levels',    &     !title of dataset
                          'pt'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(pass_tracer,    &       !input  3d field (on s-levels)
    aux_array3d_zgr_loc,    &       !output 3d field (on z-levels)
               hhq_rest,    &       !2d field of bottom topography (in metres),
                     lu,    &       !temperature mask
                      z,    &       !array of s-levels
               zlev_loc,    &       !array of z-levels (in metres)
          bnd_x1,bnd_x2,    &       !dimension on x     !dimension on x
          bnd_y1,bnd_y2,    &       !dimension on y     !dimension on y
                     nz,    &       !number of s-levels
               nlev_loc,    &       !number of z-levels
                      0,    &       !parameter of task
                 undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/ptz.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'Passive tracer at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/ptz.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
  'Passive tracer at z-levels',    &     !title of dataset
                           'pt'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!--------------------------------------------------------------------------------
if(age_output_loc>0) then
  !writing water ideal age
 
 if(sig_output_loc>0) then    !writing at s-levels
    
  ierr=0
  call pwdstd(path2data,'LOCAL/age.dat',nrec,age,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Age at s-levels is written'

  call fulfname(fname,path2data,'LOCAL/age.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
 'Ideal age at s-levels, days',   &     !title of dataset
                         'age'   )      !variable name
  endif
 
 endif !writing at s-levels

 if(z_output_loc>0) then !writing at z-levels
  
   call s2z(  age,    &       !input  3d field (on s-levels)
aux_array3d_zgr_loc,    &       !output 3d field (on z-levels)
         hhq_rest,    &       !2d field of bottom topography (in metres),
               lu,    &       !temperature mask
                z,    &       !array of s-levels
         zlev_loc,    &       !array of z-levels (in metres)
    bnd_x1,bnd_x2,    &       !dimension on x     !dimension on x
    bnd_y1,bnd_y2,    &       !dimension on y     !dimension on y
               nz,    &       !number of s-levels
         nlev_loc,    &       !number of z-levels
                0,    &       !parameter of task
           undef) 
   ierr=0

   call pwdstd(path2data,'LOCAL/agez.dat',nrec,aux_array3d_zgr_loc,lu,nx,ny,nlev_loc,m1loc,m2loc,n1loc,n2loc,1,nlev_loc,ierr)
   if(rank==0.and.deb_out_scr>0) write(*,*) 'Ideal age at z-levels is written'

   call fulfname(fname,path2data,'LOCAL/agez.dat',ierr)
   if (rank == 0) then
     call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                       nlev_loc,   &     !z-dimension
                           nrec,   &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                             1,    &     !z-grid type (0 - linear, 1 - levels)
                          zlev8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,    &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                         month,    &     !month  of the first field
                           day,    &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                         tstep,    &     !time step (in seconds
       'Ideal age at z-levels',    &     !title of dataset
                          'age'   )      !variable name
          
   endif
 
  endif !writing at z-levels

endif

!-----------------------------------------------------------------------------------------------------
if(wind_output_loc>0) then
 !writing zonal wind stress

 if(grid_shift_loc==0) then !writing on the model grid

  ierr=0  
  call pwdstd(path2data,'LOCAL/taux.dat',nrec,taux_oc,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Taux at u-grid is written'

  call fulfname(fname,path2data,'LOCAL/taux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                      nx_loc+1,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
              xu(m1loc-1:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds
         'zonal wind stress, Pa',  &     !title of dataset
                            'tx'   )     !variable name
  endif 
  
  !writing meridional wind stress
  
  ierr=0

  call pwdstd(path2data,'LOCAL/tauy.dat',nrec,tauy_oc,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Tauy at v-grid is written'

  call fulfname(fname,path2data,'LOCAL/tauy.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                      ny_loc+1,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                       xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),    &     !first x-value (if linear) or x-array (if levels)
                         dxst,     &     !x-step (if linear)
                     ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
             yv(n1loc-1:n2loc),    &     !first y-value (if linear) or x-array (if levels)
                         dyst,     &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds
    'Meridional wind stress, Pa',  &     !title of dataset
                            'ty'   )     !variable name
  endif
 
  !writing zonal wind speed 
  ierr=0
  call pwdstd(path2data,'LOCAL/uwnd.dat',nrec,windx,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Uwnd at u-grid is written'

  call fulfname(fname,path2data,'LOCAL/uwnd.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                      nx_loc+1,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            1,    &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
             xu(m1loc-1:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        1,        &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds
        'zonal wind speed, m/s',  &     !title of dataset
                          'u'    )      !variable name
  endif

  !writing meridional wind speed   
  ierr=0
  call pwdstd(path2data,'LOCAL/vwnd.dat',nrec,windy,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Vwnd at v-grid is written'

  call fulfname(fname,path2data,'LOCAL/vwnd.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                      ny_loc+1,   &     !y-dimension
                            1,    &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
              xt(m1loc:m2loc),    &     !first x-value (if linear) or x-array (if levels)
                        dxst,     &     !x-step (if linear)
                    ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
            yv(n1loc-1:n2loc),    &     !first y-value (if linear) or x-array (if levels)
                        dyst,     &     !y-step (if linear)
                        1,        &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds
  'meridional wind speed, m/s',   &     !title of dataset
                          'v'    )      !variable name
  endif

 else !writing on T-grid
 
 !writing zonal wind stress
  !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= ( taux_oc(m  ,n)*dyh(m  ,n)    &
                                  +taux_oc(m-1,n)*dyh(m-1,n) )/2.0/dy(m,n)
        endif
      enddo
    enddo
  !$omp end parallel do
  
  ierr=0
  
  call pwdstd(path2data,'LOCAL/taux.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Taux at t-grid is written'

  call fulfname(fname,path2data,'LOCAL/taux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds
         'zonal wind stress, Pa',  &     !title of dataset
                            'tx'   )     !variable name
  endif
 
 !writing meridional wind stress
  !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= ( tauy_oc(m,n  )*dxh(m,n  )    &
                                  +tauy_oc(m,n-1)*dxh(m,n-1) )/2.0/dx(m,n)
        endif
      enddo
    enddo
  !$omp end parallel do
 
  ierr=0

  call pwdstd(path2data,'LOCAL/tauy.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Tauy at t-grid is written'
  
  call fulfname(fname,path2data,'LOCAL/tauy.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_loc,    &     !x-dimension
                        ny_loc,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds
    'Meridional wind stress, Pa',  &     !title of dataset
                            'ty'   )     !variable name
   endif

!writing zonal wind speed
  ierr=0
  call pwdstd(path2data,'LOCAL/uwnd.dat',nrec,uwnd,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Uwnd at t-grid is written'

  call fulfname(fname,path2data,'LOCAL/uwnd.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            1,    &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        1,        &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds
        'zonal wind speed, m/s',  &     !title of dataset
                          'u'    )      !variable name
  endif

!writing meridional wind speed  
  ierr=0
  call pwdstd(path2data,'LOCAL/vwnd.dat',nrec,vwnd,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Vwnd at t-grid is written'

  call fulfname(fname,path2data,'LOCAL/vwnd.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            1,    &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        1,        &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds
  'meridional wind speed, m/s',   &     !title of dataset
                          'v'    )      !variable name
  endif
       
 endif

endif

!--------------------------------------------------------------------------------------------------
!writing surface preset data

if (ss_output_loc>0) then
  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if (lu(m,n)>0.5) then
        ssd(m,n,1) =tatm(m,n)
        ssd(m,n,2) =qatm(m,n)
        ssd(m,n,3) = lwr(m,n)
        ssd(m,n,4) = swr(m,n)
        ssd(m,n,5) =slpr(m,n)
        ssd(m,n,6) =rain(m,n)
        ssd(m,n,7) =snow(m,n)
        ssd(m,n,8) =runoff(m,n)
        ssd(m,n,9) =runoff_solid(m,n)
        ssd(m,n,10)=sst_obs(m,n)
        ssd(m,n,11)=sss_obs(m,n)   
      endif
    enddo
  enddo
  !$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/ssdata.dat',nrec,ssd,lu,nx,ny,11,m1loc,m2loc,n1loc,n2loc,1,11,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSdata is written'
  
  call fulfname(fname,path2data,'LOCAL/ssdata.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            11,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        0,        &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    '1-TA, 2-QA, 3-LWdw, 4-SWdw, 5-SLP, 6-rain, 7-snow, 8-runoff, 9-runoff solid, 10-sst obs, 11-sss obs',   &     !title of dataset
                        'data'    )      !variable name
  endif

  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if (lu(m,n)>0.5) then
        ssf(m,n,1) = sensheat(m,n)
        ssf(m,n,2) =  latheat(m,n)
        ssf(m,n,3) =   lw_bal(m,n)
        ssf(m,n,4) = sw_bal_atm(m,n)
        ssf(m,n,5) =  sw_bal_oc(m,n)
        ssf(m,n,6) = hf_tot_atm(m,n)
        ssf(m,n,7) =  hf_tot_oc(m,n)
        ssf(m,n,8) = wf_tot_atm(m,n)
        ssf(m,n,9) =  wf_tot_oc(m,n)
        ssf(m,n,10)= sf_tot_oc(m,n)   
      endif
    enddo
  enddo
  !$omp end parallel do

  ierr=0

  call pwdstd(path2data,'LOCAL/ssflux.dat',nrec,ssf,lu,nx,ny,10,m1loc,m2loc,n1loc,n2loc,1,10,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSflux is written'
  
  call fulfname(fname,path2data,'LOCAL/ssflux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            10,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        0,        &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    '1-SensHeat, 2-LatHeat, 3-LongWBal, 4-ShortWBalAtm, 5-SWBalOc, 6-HeatBalAtm, 7-HeatBalOc, 8-WaterBalAtm, 9-WaterBalOc, 10-SaltBalOc',   &     !title of dataset
                        'data'    )      !variable name
  endif

  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if (lu(m,n)>0.5) then
        vfl(m,n,1) = vol_flux(m,n)*RefDen
        vfl(m,n,2) = (tflux_adv(m,n)+tflux_dif(m,n))*RefDen*HeatCapWater
        vfl(m,n,3) =   swflux(m,n)*RefDen*HeatCapWater
        vfl(m,n,4) = hf_sugar(m,n)
        vfl(m,n,5) = (sflux_adv(m,n)+sflux_dif(m,n))*RefDen
      endif
    enddo
  enddo
  !$omp end parallel do

  ierr=0

  call pwdstd(path2data,'LOCAL/vflux.dat',nrec,vfl,lu,nx,ny,5,m1loc,m2loc,n1loc,n2loc,1,5,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Vflux is written'

  call fulfname(fname,path2data,'LOCAL/vflux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                             5,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear))
                        0,        &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    '1-VolFlux, 2-TemFlux, 3-SwFlux, 4-SugarFlux, 5-SalFlux',   &     !title of dataset
                        'data'    )      !variable name
  endif

endif

!--------------------------------------------------------------------------------
if(amuv_output_loc>0) then
  !writing lateral diffusivity
  ierr=0
  call pwdstd(path2data,'LOCAL/amts.dat',nrec,amts,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'AMTS is written'
  
  call fulfname(fname,path2data,'LOCAL/amts.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
  'Lateral diffusivity, m^2/s',   &     !title of dataset
                          'mu'   )      !variable name
  endif

  !writing lateral viscosity
  ierr=0
  call pwdstd(path2data,'LOCAL/amuv.dat',nrec,amuv,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'AMUV is written'
  
  call fulfname(fname,path2data,'LOCAL/amuv.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    'Lateral viscosity, m^2/s',   &     !title of dataset
                          'mu'   )      !variable name
  endif

  !writing lateral biharmonic viscosity

  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if (lu(m,n)>0.5) then
        do k=1,nz
         aux_array3d_tgr1(m,n,k)= amuv4(m,n,k)**2
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

  ierr=0
  call pwdstd(path2data,'LOCAL/amuv4.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'AMUV4 is written'
  
  call fulfname(fname,path2data,'LOCAL/amuv4.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    'Lateral biharmonic viscosity, m^4/s',   &     !title of dataset
                          'mu'   )      !variable name
  endif

endif

!--------------------------------------------------------------------------------
if(anzu_output_loc>0) then

  !writing vertical diffusivity
  ierr=0
  call pwdstd(path2data,'LOCAL/anzt.dat',nrec,anzt,lu,nx,ny,nz+1,m1loc,m2loc,n1loc,n2loc,1,nz+1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'ANZT is written'

  call fulfname(fname,path2data,'LOCAL/anzt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                          nz+1,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                          zw8,    &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
  'vertical diffusivity, m^2/s',  &     !title of dataset
                          'nu'   )      !variable name
  endif

  !writing vertical viscosity
  ierr=0
  call pwdstd(path2data,'LOCAL/anzu.dat',nrec,anzu,lu,nx,ny,nz+1,m1loc,m2loc,n1loc,n2loc,1,nz+1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'ANZU is written'

  call fulfname(fname,path2data,'LOCAL/anzu.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                          nz+1,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                          zw8,    &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
    'vertical viscosity, m^2/s',  &     !title of dataset
                          'nu'   )      !variable name
  endif
endif

!--------------------------------------------------------------------------------
if(ww_output_loc>0) then
  !writing vertical velocity
  ierr=0
  call pwdstd(path2data,'LOCAL/ww.dat',nrec,ww,lu,nx,ny,nz+1,m1loc,m2loc,n1loc,n2loc,1,nz+1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'WW is written'

  call fulfname(fname,path2data,'LOCAL/ww.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                          nz+1,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                          zw8,    &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
      'vertical velocity, m/s',   &     !title of dataset
                          'w'    )      !variable name
  endif
endif

!--------------------------------------------------------------------------------
if(ahice_output_loc>0) then
  !writing ice compactness
  
  ierr=0
  call pwdstd(path2data,'LOCAL/aice.dat',nrec,aistot,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Aice is written'

  call fulfname(fname,path2data,'LOCAL/aice.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                     xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
       'Ice compactness, 0-1',    &     !title of dataset
                          'ai'   )      !variable name
  endif

  !writing ice volume
  
  ierr=0
  call pwdstd(path2data,'LOCAL/hice.dat',nrec,hitot,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Hice is written'

  call fulfname(fname,path2data,'LOCAL/hice.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),  &     !first x-value (if linear) or x-array (if levels)
                          dxst,   &     !x-step (if linear)
                      ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc:n2loc),  &     !first y-value (if linear) or x-array (if levels)
                          dyst,   &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
   'Ice volume per square, m',    &     !title of dataset
                          'hi'   )      !variable name
  endif

  !writing snow volume
  
  ierr=0
  call pwdstd(path2data,'LOCAL/hsnow.dat',nrec,hstot,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Hsnow is written'

  call fulfname(fname,path2data,'LOCAL/hsnow.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                        nx_loc,   &     !x-dimension
                        ny_loc,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
               xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,   &     !x-step (if linear)
                      ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
               yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,   &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
  'Snow volume per square, m',    &     !title of dataset
                          'hs'   )      !variable name
  endif

endif

if(uvice_output_loc>0) then 
!-----------------------------------------------------------------------------------------------------
  if(grid_shift_loc==0) then !writing on the model grid
  
  ierr=0
  !writing zonal velocity
    call pwdstd(path2data,'LOCAL/uice.dat',nrec,uice,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Uice at u-grid is written'

    call fulfname(fname,path2data,'LOCAL/uice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                        nx_loc+1,    &     !x-dimension
                          ny_loc,    &     !y-dimension
                               1,    &     !z-dimension
                            nrec,    &     !t-dimension
                        xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                xu(m1loc-1:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                              1,     &     !z-grid type (0 - linear, 1 - levels)
                              z0,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
        'zonal ice velocity, m/s',   &     !title of dataset
                              'u'   )      !variable name
    endif
  
  ierr=0
  !writing meridional velocity
    call pwdstd(path2data,'LOCAL/vice.dat',nrec,vice,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Vice at v-grid is written'

    call fulfname(fname,path2data,'LOCAL/vice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                        ny_loc+1,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc:m2loc),    &     !first x-value (if linear) or x-array (if levels)
                          dxst,     &     !x-step (if linear)
                      ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
              yv(n1loc-1:n2loc),    &     !first y-value (if linear) or x-array (if levels)
                          dyst,     &     !y-step (if linear)
                          1,        &     !z-grid type (0 - linear, 1 - levels)
                             z0,    &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
  'meridional ice velocity, m/s',   &     !title of dataset
                            'v'    )      !variable name
    endif

  else !writing on T-grid
  
  !writing zonal velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (uice(m  ,n)*dyh(m  ,n)   &
                                 +uice(m-1,n)*dyh(m-1,n) )/2.0/dy(m,n) 
        endif
      enddo
    enddo
    !$omp end parallel do
    
    ierr=0

    call pwdstd(path2data,'LOCAL/uice.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Uice at t-grid is written'

    call fulfname(fname,path2data,'LOCAL/uice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                          nx_loc,   &     !x-dimension
                          ny_loc,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                       xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                 xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                              z0,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
        'ice zonal velocity, m/s',  &     !title of dataset
                            'u'    )      !variable name
    endif

  !writing meridional velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (vice(m,n  )*dxh(m,n  )    &
                                 +vice(m,n-1)*dxh(m,n-1) )/2.0/dx(m,n)
        endif
      enddo
    enddo
    !$omp end parallel do
  
    ierr=0

    call pwdstd(path2data,'LOCAL/vice.dat',nrec,aux_array2d_01,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Vice at t-grid is written'

    call fulfname(fname,path2data,'LOCAL/vice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                          nx_loc,    &     !x-dimension
                          ny_loc,    &     !y-dimension
                               1,    &     !z-dimension
                            nrec,    &     !t-dimension
                        xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc:m2loc),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc:n2loc),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                                1,   &     !z-grid type (0 - linear, 1 - levels)
                              z0,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
    'ice meridional velocity, m/s',  &     !title of dataset
                              'v'   )      !variable name
    endif

  endif

endif

if (rank == 0) then
 write(*,*) 'Local output is finished'
endif

endsubroutine local_output

!===================================================================================================================
subroutine cpwrite(path2data,      &
                       year,       &
                      month,       & 
                        day,       &
                       hour,       &
                     minute,       &
                     tstep,        &
                     calendar )

integer,intent(in):: year, month, day, hour, minute, calendar
character*(*),intent(in):: path2data
real(4),intent(in):: tstep

real(8) z8(nz), zw8(nz+1), z0(1), z1(1)
character fname*256
integer m, n, ierr

z8 = dble(z )*1000.0d0
zw8= dble(zw)*1000.0d0
z0 = 0.0d0
z1 = 1.0d0

 if(rank==0) write(*,*) 'Writing checkpoints for restart'

!$omp parallel do private(m,n)
do n=bnd_y1, bnd_y2
 do m=bnd_x1, bnd_x2 
  aux_array2d_01(m,n) = 0.0
 enddo
enddo
!$omp end parallel do

!----------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cptt.dat',1,tt,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cptt is written'

call fulfname(fname,path2data,'cptt.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
  'Potential temperature, C',   &     !title of dataset
                         'tt'   )     !variable name
endif
!-------------------------------------------------------------------------------------

ierr=0
call pwdstd(path2data,'cpss.dat',1,ss,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpss is written'

call fulfname(fname,path2data,'cpss.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
      'Salinity-salref, PSU',   &     !title of dataset
                         'ss'   )     !variable name
endif

ierr=0
call pwdstd(path2data,'cpden0.dat',1,den_ini,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpden0 is written'

call fulfname(fname,path2data,'cpden0.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
     'Density sigma-hmean/2',   &     !title of dataset
                        'den'   )     !variable name
endif

!----------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cpage.dat',1,age,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpage is written'

call fulfname(fname,path2data,'cpage.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
          'Ideal age, years',   &     !title of dataset
                        'age'   )     !variable name
endif

!----------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cppt.dat',1,pass_tracer,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cppt is written'

call fulfname(fname,path2data,'cppt.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
            'Passive tracer',   &     !title of dataset
                        'pt'   )     !variable name
endif

!-----------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cpuu.dat',1,uu,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpuu is written'

call fulfname(fname,path2data,'cpuu.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob+1, &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
       'Zonal velocity, m/s',   &     !title of dataset
                         'u'   )     !variable name
endif
!----------------------------------------------------------------------------------

ierr=0
call pwdstd(path2data,'cpvv.dat',1,vv,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpvv is written'

call fulfname(fname,path2data,'cpvv.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob+1, &     !y-dimension
                          nz,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
  'Meridional velocity, m/s',   &     !title of dataset
                         'v'   )     !variable name
endif
!-----------------------------------------------------------------------------------------

ierr=0
call pwdstd(path2data,'cpq2.dat',1,q2,lu,nx,ny,nz+1,mmm,mm,nnn,nn,1,nz+1,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpq2 is written'

call fulfname(fname,path2data,'cpq2.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                        nz+1,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                         zw8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
  'turbulent kunetic energy, (m/s)^2',   &     !title of dataset
                         'q2'   )     !variable name
endif
!-----------------------------------------------------------------------------------------

ierr=0
call pwdstd(path2data,'cpq2l.dat',1,q2l,lu,nx,ny,nz+1,mmm,mm,nnn,nn,1,nz+1,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpq2l is written'

call fulfname(fname,path2data,'cpq2l.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                        nz+1,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                         zw8,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
  'turbulent kunetic energy by length scale, (m/s)^2*m',   &     !title of dataset
                        'q2l'   )     !variable name
endif
!-----------------------------------------------------------------------------------------

ierr=0
call pwdstd8(path2data,'cpssh8.dat',1, ssh, lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpssh is written'

ierr=0
aux_array2d_01=sngl(ssh)
call pwdstd(path2data,'cpssh.dat',1,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
call fulfname(fname,path2data,'cpssh.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
                    'SSH, m',   &     !title of dataset
                        'ssh'   )     !variable name
endif

!--------------------------------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cpaice.dat',1,aice,lu,nx,ny,mgrad,mmm,mm,nnn,nn,1,mgrad,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpaice is written'

call fulfname(fname,path2data,'cpaice.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                       mgrad,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           0,   &     !z-grid type (0 - linear, 1 - levels)
                          z1,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
      'ice compactness, 0-1',   &     !title of dataset
                         'ai'   )     !variable name
endif

!-------------------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cphice.dat',1,hice,lu,nx,ny,mgrad,mmm,mm,nnn,nn,1,mgrad,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cphice is written'

call fulfname(fname,path2data,'cphice.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                       mgrad,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           0,   &     !z-grid type (0 - linear, 1 - levels)
                          z1,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
     'ice volume per area, m',   &     !title of dataset
                         'hi'   )     !variable name
endif

!------------------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cphsnow.dat',1,hsnow,lu,nx,ny,mgrad,mmm,mm,nnn,nn,1,mgrad,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cphsnow is written'

call fulfname(fname,path2data,'cphsnow.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                       mgrad,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           0,   &     !z-grid type (0 - linear, 1 - levels)
                          z1,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'snow volume per area, m',   &     !title of dataset
                         'hs'   )     !variable name
endif

!----------------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cpsice.dat',1,sice,lu,nx,ny,mgrad,mmm,mm,nnn,nn,1,mgrad,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpsice is written'

call fulfname(fname,path2data,'cpsice.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                       mgrad,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           0,   &     !z-grid type (0 - linear, 1 - levels)
                          z1,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
  'ice salt per area, PSU*m',   &     !title of dataset
                         'si'   )     !variable name
endif

!--------------------------------------------------------------------------------------
ierr=0
call pwdstd(path2data,'cpuice.dat',1,uice,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpuice is written'

call fulfname(fname,path2data,'cpuice.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob+1, &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                  yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'Zonal ice velocity, m/s',   &     !title of dataset
                         'u'   )     !variable name
endif
!----------------------------------------------------------------------------------

ierr=0
call pwdstd(path2data,'cpvice.dat',1,vice,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
if(rank==0.and.deb_out_scr>0) write(*,*) 'cpvice is written'

call fulfname(fname,path2data,'cpvice.dat',ierr)
if (rank == 0) then
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob+1, &     !y-dimension
                           1,   &     !z-dimension
                           1,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                           1,   &     !z-grid type (0 - linear, 1 - levels)
                          z0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds 
'Meridional ice velocity, m/s',   &     !title of dataset
                         'v'   )     !variable name
endif

if(rank==0) write(*,*) 'Checkpoints are written'

endsubroutine cpwrite

!============================================================================
subroutine data_calc

integer m,n,k

  meancalc = meancalc+1

  !$omp parallel do
  do n=ny_start-1,ny_end+1
    do m=nx_start-1,nx_end+1
      
      do k=1,nz
        tt_calc(m,n,k)=  tt_calc(m,n,k) + tt(m,n,k)*hhq(m,n)
        ss_calc(m,n,k)=  ss_calc(m,n,k) + ss(m,n,k)*hhq(m,n)
        uu_calc(m,n,k)=  uu_calc(m,n,k) + uu(m,n,k)*hhu(m,n)
        vv_calc(m,n,k)=  vv_calc(m,n,k) + vv(m,n,k)*hhv(m,n)
        pt_calc(m,n,k)=  pt_calc(m,n,k) + pass_tracer(m,n,k)*hhq(m,n)
       age_calc(m,n,k)= age_calc(m,n,k) + age(m,n,k)*hhq(m,n)
       den_sgt_calc(m,n,k) = den_sgt_calc(m,n,k) + den_sgt(m,n,k)*hhq(m,n)
       den_sg0_calc(m,n,k) = den_sg0_calc(m,n,k) + den_sg0(m,n,k)*hhq(m,n)
       den_sgm_calc(m,n,k) = den_sgm_calc(m,n,k) + den_sgm(m,n,k)*hhq(m,n)
      enddo
       sls_calc(m,n)= sls_calc(m,n) + sls(m,n)
       mld_dens_calc(m,n)= mld_dens_calc(m,n) + mld_dens(m,n)
       hhq_calc(m,n)=  hhq_calc(m,n)+ hhq(m,n)
       hhu_calc(m,n)=  hhu_calc(m,n)+ hhu(m,n)
       hhv_calc(m,n)=  hhv_calc(m,n)+ hhv(m,n)
       ssh_calc(m,n)=  ssh_calc(m,n)+ ssh(m,n) +mistot(m,n)/RefDen*float(variable_volume_budget)
      ssh_hhq_calc(m,n)=  ssh_hhq_calc(m,n)+ ssh(m,n)            
     uwnd_calc(m,n)= uwnd_calc(m,n)+uwnd(m,n)
     vwnd_calc(m,n)= vwnd_calc(m,n)+vwnd(m,n)
      txo_calc(m,n)=  txo_calc(m,n)+taux_oc(m,n)
      tyo_calc(m,n)=  tyo_calc(m,n)+tauy_oc(m,n)

      ssdata_calc(m,n,1)  = ssdata_calc(m,n,1)  +   tatm(m,n)
      ssdata_calc(m,n,2)  = ssdata_calc(m,n,2)  +   qatm(m,n)
      ssdata_calc(m,n,3)  = ssdata_calc(m,n,3)  +    lwr(m,n)
      ssdata_calc(m,n,4)  = ssdata_calc(m,n,4)  +    swr(m,n)
      ssdata_calc(m,n,5)  = ssdata_calc(m,n,5)  +   slpr(m,n)
      ssdata_calc(m,n,6)  = ssdata_calc(m,n,6)  +   rain(m,n)
      ssdata_calc(m,n,7)  = ssdata_calc(m,n,7)  +   snow(m,n)
      ssdata_calc(m,n,8)  = ssdata_calc(m,n,8)  + runoff(m,n)
      ssdata_calc(m,n,9)  = ssdata_calc(m,n,9)  + runoff_solid(m,n)
      ssdata_calc(m,n,10) = ssdata_calc(m,n,10)  +sst_obs(m,n)
      ssdata_calc(m,n,11) = ssdata_calc(m,n,11) +sss_obs(m,n)      

      ssflux_calc(m,n,1)  = ssflux_calc(m,n,1)  +  sensheat(m,n)
      ssflux_calc(m,n,2)  = ssflux_calc(m,n,2)  +   latheat(m,n)
      ssflux_calc(m,n,3)  = ssflux_calc(m,n,3)  +    lw_bal(m,n)
      ssflux_calc(m,n,4)  = ssflux_calc(m,n,4)  +sw_bal_atm(m,n)
      ssflux_calc(m,n,5)  = ssflux_calc(m,n,5)  + sw_bal_oc(m,n)
      ssflux_calc(m,n,6)  = ssflux_calc(m,n,6)  +hf_tot_atm(m,n)
      ssflux_calc(m,n,7)  = ssflux_calc(m,n,7)  + hf_tot_oc(m,n)
      ssflux_calc(m,n,8)  = ssflux_calc(m,n,8)  +wf_tot_atm(m,n)
      ssflux_calc(m,n,9)  = ssflux_calc(m,n,9)  + wf_tot_oc(m,n)
      ssflux_calc(m,n,10) = ssflux_calc(m,n,10) + sf_tot_oc(m,n) 

      vflux_calc(m,n,1)  = vflux_calc(m,n,1)  + vol_flux(m,n)*RefDen
      vflux_calc(m,n,2)  = vflux_calc(m,n,2)  + (tflux_adv(m,n)+tflux_dif(m,n))*RefDen*HeatCapWater
      vflux_calc(m,n,3)  = vflux_calc(m,n,3)  +   swflux(m,n)*RefDen*HeatCapWater
      vflux_calc(m,n,4)  = vflux_calc(m,n,4)  + hf_sugar(m,n)
      vflux_calc(m,n,5)  = vflux_calc(m,n,5)  + (sflux_adv(m,n)+sflux_dif(m,n))*RefDen

      aice_calc(m,n) =aice_calc(m,n) +aistot(m,n)
      hice_calc(m,n) =hice_calc(m,n) + hitot(m,n)
      hsnow_calc(m,n)=hsnow_calc(m,n)+ hstot(m,n)
      uice_calc(m,n) =uice_calc(m,n) +  uice(m,n)
      vice_calc(m,n) =vice_calc(m,n) +  vice(m,n)
    enddo
  enddo
  !$omp end parallel do

endsubroutine data_calc

!==============================================================================================
subroutine global_output(path2data,  &
                         nrec,       &
                         year,       &
                        month,       & 
                          day,       &
                         hour,       &
                       minute,       &
                       tstep,        &
                       calendar  )

include 'globout.fi'

integer, intent(in):: nrec, year, month, day, hour, minute, calendar
character fname*256
character*(*), intent(in):: path2data

real(4) tstep
integer m,n,k, ierr
real(8) z8(nz), zw8(nz+1), z0(1), z1(1), zlev8(nlev_glob)

z8 = dble(z )*1000.0d0
zw8= dble(zw)*1000.0d0
z0 = 0.0d0
z1 = 1.0d0
zlev8=dble(zlev_glob)

!$omp parallel do private(m,n,k)
do n=bnd_y1, bnd_y2
 do m=bnd_x1, bnd_x2 
  
  aux_array2d_01(m,n) = 0.0
  aux_array2d_02(m,n) = 0.0
    
  do k=1, nz 
   aux_array3d_tgr1(m,n,k)=0.0
  enddo

  do k=1, nlev_glob
   aux_array3d_zgr_glob(m,n,k)=0.0
  enddo

 enddo
enddo
!$omp end parallel do

if (rank == 0) then
  write(*,*) 'Writing global output (averaged for the period), record number ', nrec
endif

if(nrec==1) then
  !writing HHQ
  ierr=0
  call pwdstd(path2data,'GLOBAL/hhq.dat',nrec,hhq_rest,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'HHQ is written'

  call fulfname(fname,path2data,'GLOBAL/hhq.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
                      'HHQ, m',   &     !title of dataset
                         'hhq'   )      !variable name
  endif

  !writing geographical longitude on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=geo_lon_t(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/geolon_t.dat',nrec,aux_array2d_01,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Geolon is written'

  call fulfname(fname,path2data,'GLOBAL/geolon_t.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                       nx_glob,    &     !x-dimension
                       ny_glob,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Geographical longitude on t-grid, degree',    &     !title of dataset
                          'glon'   )      !variable name
  endif

  !writing geographical latitude on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=geo_lat_t(m,n)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/geolat_t.dat',nrec,aux_array2d_01,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Geolat is written'

  call fulfname(fname,path2data,'GLOBAL/geolat_t.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                       nx_glob,    &     !x-dimension
                       ny_glob,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Geographical latitude on t-grid, degree',    &     !title of dataset
                          'glat'   )      !variable name
  endif

!writing cosine of angles between model and geographical grids on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=rotvec_coeff(m,n,1)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/cosrot.dat',nrec,aux_array2d_01,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Cosrot is written'

  call fulfname(fname,path2data,'GLOBAL/cosrot.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                       nx_glob,    &     !x-dimension
                       ny_glob,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Cosine of angle between model and geographical grid lines',    &     !title of dataset
                          'cosrot'   )      !variable name
  endif

!writing sine of angles between model and geographical grids on t-grid

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
   do m=nx_start, nx_end
    aux_array2d_01(m,n)=rotvec_coeff(m,n,2)
   enddo
  enddo
!$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/sinrot.dat',nrec,aux_array2d_01,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Sinrot is written'

  call fulfname(fname,path2data,'GLOBAL/sinrot.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                       nx_glob,    &     !x-dimension
                       ny_glob,    &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                      xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
    'Sine of angle between model and geographical grid lines',    &     !title of dataset
                          'sinrot'   )      !variable name
  endif

endif

if(ssh_output_glob>0) then
  !writing SSH
  ierr=0

  !$omp parallel do private(m,n,k)
  do n=bnd_y1, bnd_y2
   do m=bnd_x1, bnd_x2 
    aux_array2d_01(m,n)=sngl(ssh_calc(m,n)/dfloat(meancalc))    
   enddo
  enddo
  !$omp end parallel do


  call pwdstd(path2data,'GLOBAL/ssh.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSH is written'

  call fulfname(fname,path2data,'GLOBAL/ssh.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
                      'SSH, m',   &     !title of dataset
                         'ssh'   )      !variable name
  endif

endif

!------------------------------------------------------------------------------------------
if(mld_output_glob>0) then
!mixed layer depths writing

  !$omp parallel do private(m,n,k)
  do n=bnd_y1, bnd_y2
   do m=bnd_x1, bnd_x2   
    aux_array2d_01(m,n) = mld_dens_calc(m,n)/float(meancalc)
   enddo
  enddo
  !$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/mld_dens.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'MLD-DENS is written'

  call fulfname(fname,path2data,'GLOBAL/mld_dens.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_glob,   &     !x-dimension
                        ny_glob,   &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                       xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
 'Mixed layer depth based on density difference, m',    &     !title of dataset
                          'mld'   )      !variable name
  endif

endif

if(sls_output_glob>0) then
!mixed layer depths writing

  ierr=0
  !$omp parallel do private(m,n,k)
  do n=bnd_y1, bnd_y2
   do m=bnd_x1, bnd_x2   
    aux_array2d_01(m,n) = sls_calc(m,n)/float(meancalc)
   enddo
  enddo
  !$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/sls.dat',nrec,sls,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SLS is written'

  call fulfname(fname,path2data,'GLOBAL/sls.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,     &     !file name
                        undef,     &     !value for undefined points
                        nx_glob,   &     !x-dimension
                        ny_glob,   &     !y-dimension
                              1,   &     !z-dimension
                          nrec,    &     !t-dimension
                       xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                      calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,    &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                          hour,    &     !hour   of the first field
                        minute,    &     !minute of the first field
                          tstep,   &     !time step (in seconds)
         'Steric sea level, m',    &     !title of dataset
                          'sls'   )      !variable name
  endif

endif
  
!--------------------------------------------------------------------------------
if(ts_output_glob>0) then
  !writing temperature

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=tt_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/tt.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'T at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/tt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
'Temperature at s-levels, °C',   &     !title of dataset
                         'tt'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/tz.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'T at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/tz.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
'Temperature at z-levels, °C',   &     !title of dataset
                         'tt'   )      !variable name
  endif
 
 endif

 !writing salinity

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
          aux_array3d_tgr1(m,n,k)=ss_calc(m,n,k)/hhq_calc(m,n)+salref
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/ss.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'S at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/ss.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                            nz,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         1,       &     !z-grid type (0 - linear, 1 - levels)
                            z8,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
   'Salinity at s-levels, PSU',   &     !title of dataset
                          'ss'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
                    hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/sz.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'S at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/sz.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
  'Salinity at z-levels, PSU',   &     !title of dataset
                         'ss'   )      !variable name
  endif
 
 endif

endif

!--------------------------------------------------------------------------------
if(den_sgt_output_glob>0) then
  !writing density in situ

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=den_sgt_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sgt.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-t at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sgt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-t, kg/m^3',    &     !title of dataset
                         'den'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sgt_z.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-t at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sgt_z.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-t, kg/m^3',    &     !title of dataset
                        'den'   )      !variable name
  endif
 
 endif

endif

!--------------------------------------------------------------------------------
if(den_sg0_output_glob>0) then
  !writing density SLP potential-1000

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=den_sg0_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sg0.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-0 at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sg0.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-0, kg/m^3',    &     !title of dataset
                         'den'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sg0_z.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-0 at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sg0_z.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-0, kg/m^3',    &     !title of dataset
                        'den'   )      !variable name
  endif
 
 endif

endif

!--------------------------------------------------------------------------------
if(den_sgm_output_glob>0) then
  !writing density potential for z=hmean/2

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=den_sgm_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sgm.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-hmean/2 at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sgm.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-hmean/2, kg/m^3',    &     !title of dataset
                         'den'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
             hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/den_sgm_z.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Density sigma-hmean/2 at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/den_sgm_z.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
   'Density sigma-hmean/2, kg/m^3',    &     !title of dataset
                        'den'   )      !variable name
  endif
 
 endif

endif

!--------------------------------------------------------------------------------
if(pt_output_glob>0) then
  !writing passive tracer

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=pt_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/pt.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Passive tracer at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/pt.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
 'Passive tracer at s-levels',   &     !title of dataset
                         'pt'   )      !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
                    hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/ptz.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Passive tracer at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/ptz.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
 'Passive tracer at z-levels',   &     !title of dataset
                         'pt'   )      !variable name
  endif
 
 endif

endif

!--------------------------------------------------------------------------------
if(age_output_glob>0) then
  !writing ideal age

  !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
        do k=1,nz
           aux_array3d_tgr1(m,n,k)=age_calc(m,n,k)/hhq_calc(m,n)
        enddo
      endif
    enddo
  enddo
  !$omp end parallel do

 if(sig_output_glob>0) then
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/age.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Ideal age at s-levels is written'

  call fulfname(fname,path2data,'GLOBAL/age.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                           nz,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                           z8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
      'Ideal age at s-levels',   &     !title of dataset
                         'age'   )     !variable name
  endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
                    hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0
  call pwdstd(path2data,'GLOBAL/agez.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Ideal age at z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/agez.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
      'Ideal age at z-levels',   &     !title of dataset
                         'age'   )      !variable name
  endif
 
 endif

endif

if(uv_output_glob>0) then 
!-----------------------------------------------------------------------------------------------------

  if (grid_shift_glob==0) then !writing on the model grid
  
  !writing zonal velocity
  ierr=0
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start-1, nx_end
        if (llu(m,n)>0.5) then
          do k=1,nz
            aux_array3d_tgr1(m,n,k)=uu_calc(m,n,k)/hhu_calc(m,n)
          enddo
        endif
      enddo
    enddo
    !$omp end parallel do
  
    call pwdstd(path2data,'GLOBAL/uu.dat',nrec,aux_array3d_tgr1,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'U at u-grid s-levels is written'

    call fulfname(fname,path2data,'GLOBAL/uu.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                        nx_glob+1,  &     !x-dimension
                        ny_glob,    &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
            'zonal velocity, m/s',  &     !title of dataset
                             'u'   )      !variable name
    endif

!writing meridional velocity
  ierr=0
     !$omp parallel do
     do n=ny_start-1, ny_end
       do m=nx_start, nx_end
         if(llv(m,n)>0.5) then
           do k=1,nz
             aux_array3d_tgr1(m,n,k)=vv_calc(m,n,k)/hhv_calc(m,n)
           enddo
         endif
        enddo
    enddo
    !$omp end parallel do

    call pwdstd(path2data,'GLOBAL/vv.dat',nrec,aux_array3d_tgr1,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'V at v-grid s-levels is written'

    call fulfname(fname,path2data,'GLOBAL/vv.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                       ny_glob+1,   &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
       'meridional velocity, m/s',  &     !title of dataset
                             'v'   )      !variable name
    endif

  else !writing on T-grid

!writing zonal velocity 

    !$omp parallel do
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if(lu(m,n)>0.5) then
          do k=1,nz
            aux_array3d_tgr1(m,n,k)= ( uu_calc(m  ,n,k)*dyh(m  ,n)  &
                                      +uu_calc(m-1,n,k)*dyh(m-1,n) )/2.0/hhq_calc(m,n)/dy(m,n)
          enddo
        endif
      enddo
    enddo
 !$omp end parallel do

 if(sig_output_glob>0) then

    ierr=0 
    call pwdstd(path2data,'GLOBAL/uu.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'U at t-grid s-levels is written'

    call fulfname(fname,path2data,'GLOBAL/uu.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                         ny_glob,   &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
'zonal velocity at s-levels, m/s',  &     !title of dataset
                             'u'   )      !variable name
    endif
 
 endif

 if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
                    hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0

  call pwdstd(path2data,'GLOBAL/uz.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'U at t-grid z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/uz.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
'zonal velocity at z-levels, m/s',  &     !title of dataset
                          'u'   )      !variable name
  endif
 
 endif

!writing meridional velocity

    !$omp parallel do
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if(lu(m,n)>0.5) then
          do k=1,nz
            aux_array3d_tgr1(m,n,k)=  ( vv_calc(m,n  ,k)*dxh(m,n  )  &
                                       +vv_calc(m,n-1,k)*dxh(m,n-1) )/2.0/hhq_calc(m,n)/dx(m,n)
          enddo
        endif
      enddo
    enddo
    !$omp end parallel do

 if(sig_output_glob>0) then

  ierr=0  
    call pwdstd(path2data,'GLOBAL/vv.dat',nrec,aux_array3d_tgr1,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'V at t-grid s-levels is written'

    call fulfname(fname,path2data,'GLOBAL/vv.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                         ny_glob,   &     !y-dimension
                              nz,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z8,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
 'meridional velocity at s-levels, m/s',  &     !title of dataset
                             'v'   )      !variable name
    endif
   
  endif
  
  if(z_output_glob>0) then
   
   call s2z(aux_array3d_tgr1,      &       !input  3d field (on s-levels)
            aux_array3d_zgr_glob,  &       !output 3d field (on z-levels)
                    hhq_rest,      &       !2d field of bottom topography (in metres),
                   lu,      &       !temperature mask
                    z,      &       !array of s-levels
            zlev_glob,      &       !array of z-levels (in metres)
        bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
        bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
                   nz,      &       !number of s-levels
            nlev_glob,      &       !number of z-levels
                    0,      &       !parameter of task
               undef)   
  
  ierr=0

  call pwdstd(path2data,'GLOBAL/vz.dat',nrec,aux_array3d_zgr_glob,lu,nx,ny,nlev_glob,mmm,mm,nnn,nn,1,nlev_glob,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'V at t-grid z-levels is written'

  call fulfname(fname,path2data,'GLOBAL/vz.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,   &     !file name
                       undef,    &     !value for undefined points
                      nx_glob,   &     !x-dimension
                      ny_glob,   &     !y-dimension
                    nlev_glob,   &     !z-dimension
                         nrec,   &     !t-dimension
                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,    &     !x-step (if linear)
                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,    &     !y-step (if linear)
                        1,       &     !z-grid type (0 - linear, 1 - levels)
                        zlev8,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,   &     !z-step (if linear)
                     calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                         year,   &     !year   of the first field
                        month,   &     !month  of the first field
                          day,   &     !day    of the first field
                         hour,   &     !hour   of the first field
                       minute,   &     !minute of the first field
                        tstep,   &     !time step (in seconds)
 'meridional velocity at z-levels, m/s',  &     !title of dataset
                          'v'   )      !variable name
  endif
 
 endif ! z-output

 endif !t-grid

endif !uv

!-----------------------------------------------------------------------------------------------
!write wind speed
if(wind_output_glob>0) then

  if(grid_shift_glob==0) then
  !writing zonal velocity
  ierr=0
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start-1, nx_end
        if (llu(m,n)>0.5) then
            aux_array2d_01(m,n)= (uwnd_calc(m  ,n)*dy(m  ,n)   &
                                 +uwnd_calc(m+1,n)*dy(m+1,n) )/2.0/dyh(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do
  
  call pwdstd(path2data,'GLOBAL/uwnd.dat',nrec,aux_array2d_01,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Uwnd at u-grid is written'

  call fulfname(fname,path2data,'GLOBAL/uwnd.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                     nx_glob+1,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                            1,    &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                        1,        &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds
        'zonal wind speed, m/s',  &     !title of dataset
                          'u'    )      !variable name
  endif

  !writing meridional velocity
  ierr=0
    !$omp parallel do 
    do n=ny_start-1, ny_end
      do m=nx_start, nx_end
        if (llv(m,n)>0.5) then
            aux_array2d_01(m,n)= (vwnd_calc(m,n  )*dx(m,n  )   &
                                 +vwnd_calc(m,n+1)*dx(m,n+1) )/2.0/dxh(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do
  
    call pwdstd(path2data,'GLOBAL/vwnd.dat',nrec,aux_array2d_01,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Vwnd at v-grid is written'

    call fulfname(fname,path2data,'GLOBAL/vwnd.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                       ny_glob+1,   &     !y-dimension
                              1,    &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                          1,        &     !z-grid type (0 - linear, 1 - levels)
                              z0,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
     'meridional wind speed, m/s',  &     !title of dataset
                            'v'    )      !variable name
    endif

!writing zonal wind stress
  ierr=0
  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start-1, nx_end
      if (llu(m,n)>0.5) then
          aux_array2d_01(m,n)= txo_calc(m,n)/float(meancalc)
      endif
    enddo
  enddo
  !$omp end parallel do

  call pwdstd(path2data,'GLOBAL/taux.dat',nrec,aux_array2d_01,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Taux at u-grid is written'

  call fulfname(fname,path2data,'GLOBAL/taux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                     nx_glob+1,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds
        'zonal wind stress, Pa',  &     !title of dataset
                           'tx'   )      !variable name
  endif

!------------------------------------------------------------------------------------------------------
!writing meridional wind stress
  ierr=0
  !$omp parallel do 
  do n=ny_start-1, ny_end
    do m=nx_start, nx_end
      if (llv(m,n)>0.5) then
          aux_array2d_01(m,n)= tyo_calc(m,n)/float(meancalc)
      endif
    enddo
  enddo
  !$omp end parallel do

  call pwdstd(path2data,'GLOBAL/tauy.dat',nrec,aux_array2d_01,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Tauy at v-grid is written'

  call fulfname(fname,path2data,'GLOBAL/tauy.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                     ny_glob+1,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         1,       &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds
   'meridional wind stress, Pa',  &     !title of dataset
                           'ty'   )      !variable name
  endif

  else
  
   !writing zonal velocity
   ierr=0
   !$omp parallel do 
   do n=ny_start, ny_end
     do m=nx_start, nx_end
       if (lu(m,n)>0.5) then
           aux_array2d_01(m,n)= uwnd_calc(m,n)/float(meancalc)
       endif
     enddo
   enddo
   !$omp end parallel do
   
     call pwdstd(path2data,'GLOBAL/uwnd.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr) 
     if(rank==0.and.deb_out_scr>0) write(*,*) 'Uwnd at t-grid is written'
     
     call fulfname(fname,path2data,'GLOBAL/uwnd.dat',ierr)
     if (rank == 0) then
       call ctl_file_write(fname,    &     !file name
                           undef,    &     !value for undefined points
                          nx_glob,   &     !x-dimension
                          ny_glob,   &     !y-dimension
                                1,   &     !z-dimension
                             nrec,   &     !t-dimension
                         xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                       xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                       yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                            1,       &     !z-grid type (0 - linear, 1 - levels)
                               z0,   &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                         calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                             year,   &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                             hour,   &     !hour   of the first field
                           minute,   &     !minute of the first field
                            tstep,   &     !time step (in seconds
           'zonal wind speed, m/s',  &     !title of dataset
                              'u'   )      !variable name
     endif

  !writing meridional velocity
   ierr=0
   !$omp parallel do 
   do n=ny_start, ny_end
     do m=nx_start, nx_end
       if (lu(m,n)>0.5) then
           aux_array2d_01(m,n)= vwnd_calc(m,n)/float(meancalc)
       endif
     enddo
   enddo
   !$omp end parallel do
   
     call pwdstd(path2data,'GLOBAL/vwnd.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
     if(rank==0.and.deb_out_scr>0) write(*,*) 'Vwnd at t-grid is written'

     call fulfname(fname,path2data,'GLOBAL/vwnd.dat',ierr)
     if (rank == 0) then
       call ctl_file_write(fname,    &     !file name
                           undef,    &     !value for undefined points
                          nx_glob,   &     !x-dimension
                          ny_glob,   &     !y-dimension
                                1,   &     !z-dimension
                             nrec,   &     !t-dimension
                         xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                       xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                       yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                            1,       &     !z-grid type (0 - linear, 1 - levels)
                               z0,   &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                         calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                             year,   &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                             hour,   &     !hour   of the first field
                           minute,   &     !minute of the first field
                            tstep,   &     !time step (in seconds
      'meridional wind speed, m/s',  &     !title of dataset
                              'v'   )      !variable name
     endif

  !writing zonal wind stress
    ierr=0
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (txo_calc(m  ,n)*dyh(m  ,n)   &
                                 +txo_calc(m-1,n)*dyh(m-1,n) )/2.0/dy(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do

    call pwdstd(path2data,'GLOBAL/taux.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Taux at t-grid is written'

    call fulfname(fname,path2data,'GLOBAL/taux.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                         ny_glob,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z0,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
          'zonal wind stress, Pa',  &     !title of dataset
                             'tx'   )      !variable name
    endif
  
  !------------------------------------------------------------------------------------------------------
  !writing meridional wind stress
    ierr=0
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (tyo_calc(m,n  )*dxh(m,n  )   &
                                 +tyo_calc(m,n-1)*dxh(m,n-1) )/2.0/dx(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do

    call pwdstd(path2data,'GLOBAL/tauy.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Tauy at t-grid is written'

    call fulfname(fname,path2data,'GLOBAL/tauy.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                         ny_glob,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                           1,       &     !z-grid type (0 - linear, 1 - levels)
                              z0,   &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                           month,   &     !month  of the first field
                             day,   &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                           tstep,   &     !time step (in seconds
     'meridional wind stress, Pa',  &     !title of dataset
                             'ty'   )      !variable name
    endif

  endif
endif

!--------------------------------------------------------------------------------------------------
!writing surface data
if(ss_output_glob>0) then
  !$omp parallel do 
  do n=ny_start, ny_end
    do m=nx_start, nx_end
      if(lu(m,n)>0.5) then
          ssdata_calc(m,n,:)=ssdata_calc(m,n,:)/float(meancalc)      
          ssflux_calc(m,n,:)=ssflux_calc(m,n,:)/float(meancalc)      
           vflux_calc(m,n,:)= vflux_calc(m,n,:)/float(meancalc)      
      endif
    enddo
  enddo
  !$omp end parallel do

  ierr=0
  call pwdstd(path2data,'GLOBAL/ssdata.dat',nrec,ssdata_calc,lu,nx,ny,11,mmm,mm,nnn,nn,1,11,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSdata is written'

  call fulfname(fname,path2data,'GLOBAL/ssdata.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                            11,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         0,       &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
    '1-TA, 2-QA, 3-LWdw, 4-SWdw, 5-SLP, 6-rain, 7-snow, 8-runoff, 9-runoff solid, 10-sst obs, 11-sss obs',   &     !title of dataset
                         'data'   )      !variable name
  endif

  ierr=0
  call pwdstd(path2data,'GLOBAL/ssflux.dat',nrec,ssflux_calc,lu,nx,ny,10,mmm,mm,nnn,nn,1,10,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'SSflux is written'

  call fulfname(fname,path2data,'GLOBAL/ssflux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                            10,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         0,       &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
    '1-SensHeat, 2-LatHeat, 3-LongWBal, 4-ShortWBalAtm, 5-SWBalOc, 6-HeatBalAtm, 7-HeatBalOc, 8-WaterBalAtm, 9-WaterBalOc, 10-SaltBalOc',   &     !title of dataset
                         'data'   )      !variable name
  endif

  ierr=0
  call pwdstd(path2data,'GLOBAL/vflux.dat',nrec,vflux_calc,lu,nx,ny,5,mmm,mm,nnn,nn,1,5,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Vflux is written'

  call fulfname(fname,path2data,'GLOBAL/vflux.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             5,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                         0,       &     !z-grid type (0 - linear, 1 - levels)
                            z1,   &     !first z-value (if linear) or x-array (if levels)
                         1.0d0,   &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                         month,   &     !month  of the first field
                           day,   &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                         tstep,   &     !time step (in seconds)
    '1-VolFlux, 2-TemFlux, 3-SwFlux, 4-SugarFlux, 5-SalFlux',   &     !title of dataset
                         'data'   )      !variable name
  endif

endif

!--------------------------------------------------------------------------------
if(ahice_output_glob>0) then
  !writing ice compactness
  ierr=0
  
  !$omp parallel do private(m,n)
    do n=ny_start,ny_end
     do m=nx_start, nx_end
      aux_array2d_01(m,n)=aice_calc(m,n)/float(meancalc)
     enddo
    enddo
  !$omp end parallel do

  call pwdstd(path2data,'GLOBAL/aice.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Aice is written'

  call fulfname(fname,path2data,'GLOBAL/aice.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
       'Ice compactness, 0-1',    &     !title of dataset
                          'ai'   )      !variable name
  endif

  !writing ice volume
  
  ierr=0
  
  !$omp parallel do private(m,n)
    do n=ny_start,ny_end
     do m=nx_start, nx_end
      aux_array2d_01(m,n)=hice_calc(m,n)/float(meancalc)
     enddo
    enddo
  !$omp end parallel do

  call pwdstd(path2data,'GLOBAL/hice.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Hice is written'

  call fulfname(fname,path2data,'GLOBAL/hice.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
   'Ice volume per square, m',    &     !title of dataset
                          'hi'   )      !variable name
  endif

  !writing snow volume
  
  ierr=0
  
  !$omp parallel do private(m,n)
    do n=ny_start,ny_end
     do m=nx_start, nx_end
      aux_array2d_01(m,n)=hsnow_calc(m,n)/float(meancalc)
     enddo
    enddo
  !$omp end parallel do

  call pwdstd(path2data,'GLOBAL/hsnow.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  if(rank==0.and.deb_out_scr>0) write(*,*) 'Hsnow is written'

  call fulfname(fname,path2data,'GLOBAL/hsnow.dat',ierr)
  if (rank == 0) then
    call ctl_file_write(fname,    &     !file name
                        undef,    &     !value for undefined points
                       nx_glob,   &     !x-dimension
                       ny_glob,   &     !y-dimension
                             1,   &     !z-dimension
                          nrec,   &     !t-dimension
                      xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                         dxst,    &     !x-step (if linear)
                     ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                         dyst,    &     !y-step (if linear)
                            1,    &     !z-grid type (0 - linear, 1 - levels)
                            z0,   &     !first z-value (if linear) or x-array (if levels)
                        1.0d0,    &     !z-step (if linear)
                      calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                          year,   &     !year   of the first field
                        month,    &     !month  of the first field
                          day,    &     !day    of the first field
                          hour,   &     !hour   of the first field
                        minute,   &     !minute of the first field
                        tstep,    &     !time step (in seconds)
  'Snow volume per square, m',    &     !title of dataset
                          'hs'   )      !variable name
  endif
endif

if(uvice_output_glob>0) then 
!-----------------------------------------------------------------------------------------------------
  if(grid_shift_glob==0) then !writing on the model grid
  
  uice_calc=uice_calc/float(meancalc)
  ierr=0
  !writing zonal velocity
    call pwdstd(path2data,'GLOBAL/uice.dat',nrec,uice_calc,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Uice at u-grid is written'

    call fulfname(fname,path2data,'GLOBAL/uice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                       nx_glob+1,    &     !x-dimension
                         ny_glob,    &     !y-dimension
                               1,    &     !z-dimension
                            nrec,    &     !t-dimension
                         xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xu(mmm-1:mm),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                       yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                              1,     &     !z-grid type (0 - linear, 1 - levels)
                              z0,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
        'zonal ice velocity, m/s',   &     !title of dataset
                              'u'   )      !variable name
    endif

  vice_calc=vice_calc/float(meancalc)  
  ierr=0
  !writing meridional velocity
    call pwdstd(path2data,'GLOBAL/vice.dat',nrec,vice_calc,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Vice at v-grid is written'

    call fulfname(fname,path2data,'GLOBAL/vice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                       ny_glob+1,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yv(nnn-1:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                          1,        &     !z-grid type (0 - linear, 1 - levels)
                             z0,    &     !first z-value (if linear) or x-array (if levels)
                           1.0d0,   &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
  'meridional ice velocity, m/s',   &     !title of dataset
                            'v'    )      !variable name
    endif

  else !writing on T-grid
  
  !writing zonal velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (uice_calc(m  ,n)*dyh(m  ,n)   &
                                 +uice_calc(m-1,n)*dyh(m-1,n) )/2.0/dy(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do
    
    ierr=0

    call pwdstd(path2data,'GLOBAL/uice.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Uice at t-grid is written'

    call fulfname(fname,path2data,'GLOBAL/uice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                         nx_glob,   &     !x-dimension
                         ny_glob,   &     !y-dimension
                               1,   &     !z-dimension
                            nrec,   &     !t-dimension
                        xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                      xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                           dxst,    &     !x-step (if linear)
                       ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                      yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                           dyst,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                              z0,   &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                        calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,   &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                            hour,   &     !hour   of the first field
                          minute,   &     !minute of the first field
                          tstep,    &     !time step (in seconds
        'ice zonal velocity, m/s',  &     !title of dataset
                            'u'    )      !variable name
    endif

  !writing meridional velocity
    !$omp parallel do 
    do n=ny_start, ny_end
      do m=nx_start, nx_end
        if (lu(m,n)>0.5) then
            aux_array2d_01(m,n)= (vice_calc(m,n  )*dxh(m,n  )    &
                                 +vice_calc(m,n-1)*dxh(m,n-1) )/2.0/dx(m,n)/float(meancalc)
        endif
      enddo
    enddo
    !$omp end parallel do
  
    ierr=0

    call pwdstd(path2data,'GLOBAL/vice.dat',nrec,aux_array2d_01,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'Vice at t-grid is written'

    call fulfname(fname,path2data,'GLOBAL/vice.dat',ierr)
    if (rank == 0) then
      call ctl_file_write(fname,     &     !file name
                          undef,     &     !value for undefined points
                         nx_glob,    &     !x-dimension
                         ny_glob,    &     !y-dimension
                               1,    &     !z-dimension
                            nrec,    &     !t-dimension
                         xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                       xt(mmm:mm),   &     !first x-value (if linear) or x-array (if levels)
                            dxst,    &     !x-step (if linear)
                        ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                       yt(nnn:nn),   &     !first y-value (if linear) or x-array (if levels)
                            dyst,    &     !y-step (if linear)
                                1,   &     !z-grid type (0 - linear, 1 - levels)
                              z0,    &     !first z-value (if linear) or x-array (if levels)
                            1.0d0,   &     !z-step (if linear)
                        calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                            year,    &     !year   of the first field
                            month,   &     !month  of the first field
                              day,   &     !day    of the first field
                            hour,    &     !hour   of the first field
                          minute,    &     !minute of the first field
                            tstep,   &     !time step (in seconds
    'ice meridional velocity, m/s',  &     !title of dataset
                              'v'   )      !variable name
    endif

  endif

endif



!$omp parallel do 
do n=bnd_y1, bnd_y2
  do m=bnd_x1, bnd_x2
      tt_calc(m,n,:)=0.0
      ss_calc(m,n,:)=0.0
      pt_calc(m,n,:)=0.0
     age_calc(m,n,:)=0.0
      uu_calc(m,n,:)=0.0
      vv_calc(m,n,:)=0.0
 den_sgt_calc(m,n,:)=0.0
 den_sg0_calc(m,n,:)=0.0
 den_sgm_calc(m,n,:)=0.0
      uwnd_calc(m,n)=0.0
      vwnd_calc(m,n)=0.0
       txo_calc(m,n)=0.0
       tyo_calc(m,n)=0.0
  ssdata_calc(m,n,:)=0.0
  ssflux_calc(m,n,:)=0.0
   vflux_calc(m,n,:)=0.0
       hhq_calc(m,n)=0.0
       hhu_calc(m,n)=0.0
       hhv_calc(m,n)=0.0
       ssh_calc(m,n)=0.0
   ssh_hhq_calc(m,n)=0.0
       sls_calc(m,n)=0.0
  mld_dens_calc(m,n)=0.0
      aice_calc(m,n)=0.0
      hice_calc(m,n)=0.0
     hsnow_calc(m,n)=0.0
      uice_calc(m,n)=0.0
      vice_calc(m,n)=0.0
  enddo
enddo
!$omp end parallel do

meancalc=0

if (rank == 0) then
 write(*,*) 'Global output is finished'
endif

endsubroutine global_output
!===================================================================
subroutine integral(path2data,  &
                        nrec,   &
                        year,   &
                       month,   & 
                         day,   &
                        hour,   &
                      minute,   &
                      tstep,    &
                      calendar  )

integer nrec, year, month, day, hour, minute, calendar, ierr
character fname*256
character*(*) path2data
real(4) integr_data(22)
real(4) tstep
integer m,n,k,grad

real(8) buf_real8
real(8) calc_tt,calc_tw,calc_t2d, calc_ut, calc_u2d, calc_s2d
real(8) t_ave,t_surf,t_mid,t_bot,s_ave,s_surf,s_mid,s_bot
real(8) ebcl,ebtr,ew,e1lyr,mssh,mssh2, mes, ub
real(4) north_hem, south_hem

real(8)  calc_aice_north,  calc_aice_south,       &
         calc_hice_north,  calc_hice_south,       &
        calc_hsnow_north, calc_hsnow_south,       &
         calc_sice_north,  calc_sice_south

real(8) x0(1),y0(1),z0(1)

x0=1.0d0
y0=0.0d0
z0=0.0d0

t_ave=0.0; t_surf=0.0; t_mid=0.0; t_bot=0.0
s_ave=0.0; s_surf=0.0; s_mid=0.0; s_bot=0.0
ebcl=0.0;  ebtr=0.0; ew=0.0; e1lyr=0.0; mssh=0.0; mssh2=0.0

calc_tt=0.0; calc_tw=0.0; calc_t2d=0.0
calc_ut=0.0; calc_u2d=0.0; calc_s2d=0.0

 calc_aice_north = 0.0;  calc_aice_south = 0.0
 calc_hice_north = 0.0;  calc_hice_south = 0.0
calc_hsnow_north = 0.0; calc_hsnow_south = 0.0
 calc_sice_north = 0.0;  calc_sice_south = 0.0

!$omp parallel do private(m,n,k,grad,mes,ub,north_hem, south_hem)   &
!$omp reduction(+: t_ave,t_surf,t_mid,t_bot,          &
!$omp              s_ave,s_surf,s_mid,s_bot,          &
!$omp            ebcl,ebtr,ew,e1lyr,mssh,mssh2,       &
!$omp            calc_tt, calc_tw, calc_t2d,          &
!$omp            calc_ut, calc_u2d, calc_s2d,         &
!$omp            calc_aice_north,  calc_aice_south,   &
!$omp            calc_hice_north,  calc_hice_south,   &
!$omp           calc_hsnow_north, calc_hsnow_south,   &
!$omp            calc_sice_north,  calc_sice_south)
do n=ny_start, ny_end
 do m=nx_start, nx_end
   if(lu(m,n)>0.5) then
    
    mes=sqt(m,n)
    calc_s2d=calc_s2d+mes
    mssh =mssh  + mes*(ssh(m,n) + mistot(m,n)/RefDen*float(variable_volume_budget))
    mssh2=mssh2 + mes*(ssh(m,n) + mistot(m,n)/RefDen*float(variable_volume_budget))**2

    north_hem = ( 1.0 + sign( 1.0,sngl(geo_lat_t(m,n)) ) )/2.0
    south_hem = 1.0-north_hem
    do grad=1,mgrad
      calc_aice_north = calc_aice_north + aice(m,n,grad)*mes*north_hem
      calc_hice_north = calc_hice_north + hice(m,n,grad)*mes*north_hem
      calc_sice_north = calc_sice_north + sice(m,n,grad)*mes*north_hem
     calc_hsnow_north =calc_hsnow_north +hsnow(m,n,grad)*mes*north_hem
      calc_aice_south = calc_aice_south + aice(m,n,grad)*mes*south_hem
      calc_hice_south = calc_hice_south + hice(m,n,grad)*mes*south_hem
      calc_sice_south = calc_sice_south + sice(m,n,grad)*mes*south_hem
     calc_hsnow_south =calc_hsnow_south +hsnow(m,n,grad)*mes*south_hem
    enddo

    mes=sqt(m,n)*hhq(m,n)
    calc_t2d = calc_t2d + mes
    t_surf=t_surf + mes*tt(m,n,1)
    t_mid =t_mid  + mes*tt(m,n,nz/2)
    t_bot =t_bot  + mes*tt(m,n,nz)
    s_surf=s_surf + mes*ss(m,n,1)
    s_mid =s_mid  + mes*ss(m,n,nz/2)
    s_bot =s_bot  + mes*ss(m,n,nz)

    do k=1,nz
     mes=sqt(m,n)*hhq(m,n)*dz(k)
     calc_tt = calc_tt + mes
     t_ave=t_ave + mes*tt(m,n,k)
     s_ave=s_ave + mes*ss(m,n,k)

     mes=sqt(m,n)*hhq(m,n)*hzt(k)
     calc_tw = calc_tw + mes
     ew = ew + mes*ww(m,n,k)**2
    enddo

   endif

   if(lcu(m,n)>0.5) then

    ub=ubrtr(m,n)/hhu(m,n)    
    mes=squ(m,n)*hhu(m,n)
    calc_u2d=calc_u2d+mes
    e1lyr=e1lyr + mes*uu(m,n,1)**2
    ebtr = ebtr + mes*ub**2
    
    do k=1,nz
     mes=squ(m,n)*hhu(m,n)*dz(k)
     calc_ut = calc_ut + mes
     ebcl = ebcl + mes*uu(m,n,k)**2
    enddo

   endif

   if(lcv(m,n)>0.5) then

    ub=vbrtr(m,n)/hhv(m,n)    
    mes=sqv(m,n)*hhv(m,n)
    calc_u2d = calc_u2d + mes
    ebtr = ebtr + mes*ub**2
    e1lyr=e1lyr + mes*vv(m,n,1)**2

    do k=1,nz
     mes=sqv(m,n)*hhv(m,n)*dz(k)
     calc_ut = calc_ut + mes
     ebcl = ebcl + mes*vv(m,n,k)**2
    enddo

   endif

 enddo
enddo
!$omp end parallel do

    buf_real8 = t_ave
    call mpi_allreduce(buf_real8, t_ave, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = t_surf
    call mpi_allreduce(buf_real8, t_surf, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = t_mid
    call mpi_allreduce(buf_real8, t_mid, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = t_bot
    call mpi_allreduce(buf_real8, t_bot, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = s_ave
    call mpi_allreduce(buf_real8, s_ave, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = s_surf
    call mpi_allreduce(buf_real8, s_surf, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = s_mid
    call mpi_allreduce(buf_real8, s_mid, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = s_bot
    call mpi_allreduce(buf_real8, s_bot, 1, mpi_real8, mpi_sum, cart_comm, ierr)

    buf_real8 = ebcl
    call mpi_allreduce(buf_real8, ebcl, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = ebtr
    call mpi_allreduce(buf_real8, ebtr, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = ew
    call mpi_allreduce(buf_real8, ew, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = e1lyr
    call mpi_allreduce(buf_real8, e1lyr, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = mssh
    call mpi_allreduce(buf_real8, mssh, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = mssh2
    call mpi_allreduce(buf_real8, mssh2, 1, mpi_real8, mpi_sum, cart_comm, ierr)

    buf_real8 = calc_tt
    call mpi_allreduce(buf_real8, calc_tt, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = calc_tw
    call mpi_allreduce(buf_real8, calc_tw, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = calc_t2d
    call mpi_allreduce(buf_real8, calc_t2d, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = calc_ut
    call mpi_allreduce(buf_real8, calc_ut, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = calc_u2d
    call mpi_allreduce(buf_real8, calc_u2d, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = calc_s2d
    call mpi_allreduce(buf_real8, calc_s2d, 1, mpi_real8, mpi_sum, cart_comm, ierr)

    buf_real8 =  calc_aice_north
    call mpi_allreduce(buf_real8, calc_aice_north, 1, mpi_real8, mpi_sum, cart_comm, ierr)    
    buf_real8 =  calc_hice_north
    call mpi_allreduce(buf_real8, calc_hice_north, 1, mpi_real8, mpi_sum, cart_comm, ierr) 
    buf_real8 =  calc_sice_north
    call mpi_allreduce(buf_real8, calc_sice_north, 1, mpi_real8, mpi_sum, cart_comm, ierr)    
    buf_real8 = calc_hsnow_north
    call mpi_allreduce(buf_real8,calc_hsnow_north, 1, mpi_real8, mpi_sum, cart_comm, ierr)     
    buf_real8 =  calc_aice_south
    call mpi_allreduce(buf_real8, calc_aice_south, 1, mpi_real8, mpi_sum, cart_comm, ierr)     
    buf_real8 =  calc_hice_south
    call mpi_allreduce(buf_real8, calc_hice_south, 1, mpi_real8, mpi_sum, cart_comm, ierr)     
    buf_real8 =  calc_sice_south
    call mpi_allreduce(buf_real8, calc_sice_south, 1, mpi_real8, mpi_sum, cart_comm, ierr)     
    buf_real8 = calc_hsnow_south
    call mpi_allreduce(buf_real8,calc_hsnow_south, 1, mpi_real8, mpi_sum, cart_comm, ierr) 

 if(rank==0) then
  integr_data(1)=t_ave/calc_tt
  integr_data(2)=t_surf/calc_t2d
  integr_data(3)=t_mid/calc_t2d
  integr_data(4)=t_bot/calc_t2d
  integr_data(5)=s_ave/calc_tt   +salref
  integr_data(6)=s_surf/calc_t2d +salref
  integr_data(7)=s_mid/calc_t2d  +salref
  integr_data(8)=s_bot/calc_t2d  +salref
  integr_data(9)= sqrt(ebcl/calc_ut) - sqrt(ebtr/calc_u2d)
  integr_data(10)= sqrt(ew/calc_tw)
  integr_data(11)= sqrt(ebtr/calc_u2d)
  integr_data(12)= sqrt(e1lyr/calc_u2d)
  integr_data(13)= mssh/calc_s2d
  integr_data(14)= sqrt(mssh2/calc_s2d)
  integr_data(15)= calc_aice_north/1e6
  integr_data(16)= calc_hice_north/1e9
  integr_data(17)=calc_hsnow_north/1e9
  integr_data(18)= calc_sice_north/1e9
  integr_data(19)= calc_aice_south/1e6
  integr_data(20)= calc_hice_south/1e9
  integr_data(21)=calc_hsnow_south/1e9
  integr_data(22)= calc_sice_south/1e9

  write(*,*) 'Writing integral output, record number ', nrec
   
   call fulfname(fname,path2data,'LOCAL/integrals.dat',ierr)
   open (20,file=fname, access='direct', form='unformatted',     &
                          recl=22*lrecl,err=120)
    write(20,rec=nrec,err=122) integr_data
   close(20) 

      call ctl_file_write(fname,    &     !file name
                          undef,    &     !value for undefined points
                             22,    &     !x-dimension
                              1,    &     !y-dimension
                              1,    &     !z-dimension
                           nrec,    &     !t-dimension
                              0,    &     !x-grid type (0 - linear, 1 - levels)
                             x0,    &     !first x-value (if linear) or x-array (if levels)
                          1.0d0,    &     !x-step (if linear)
                              1,    &     !y-grid type (0 - linear, 1 - levels)
                             y0,    &     !first y-value (if linear) or x-array (if levels)
                          1.0d0,    &     !y-step (if linear)
                              1,    &     !z-grid type (0 - linear, 1 - levels)
                             z0,    &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,    &     !z-step (if linear)
                       calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,    &     !year   of the first field
                          month,    &     !month  of the first field
                            day,    &     !day    of the first field
                           hour,    &     !hour   of the first field
                         minute,    &     !minute of the first field
                          tstep,    &     !time step (in seconds
'01-t_ave, 02-t_surf, 03-t_mid, 04-t_bot, 05-s_ave, 06-s_surf, 07-s_mid, 08-s_bot, 09-Ebcl, 10-Ew, &
& 11-Ebtr, 12-Ebcl1, 13-mssh, 14-ssh-rms, 15-NH ice area, 16-NH ice volume, 17-NH snow volume, 18-NH salt ice, &
& 19-SH ice area, 20-SH ice volume, 21-SH snow volume, 22-SH salt ice ',  &     !title of dataset
                            'data'   )     !variable name

! printout of averaged temperatute and salinity
   write(*,110) integr_data(1:8)
! printout velocity
   write(*,111) integr_data(9:14)
! printout ice
   write(*,112) integr_data(15:22)

 endif

return 

110   format(3x,'tsum;',5x,'t_1;',6x,'t_nz/2;',3x,'t_nz;',5x,           &    
                'ssum;',5x,'s_1;',6x,'s_nz/2;',3x,'s_nz;'/8f10.6)
111  format(2x,'vel: bcl ',5x,'ww',9x, 'btr',8x,'ek1',8x,'ssh',6x,'ssh2'/  &  
                f11.7,      f12.8,    f11.7,   f11.8   ,f11.7 , f11.7)
112  format(3x,'AR_N;',5x,'ViN;',6x,'V_sn_N;',3x,'Si_N;',5x,           &    
               'AR_S;',5x,'ViS;',6x,'V_sn_S;',3x,'Si_S;'/8es10.2)

120 write(*,'(1x,a)')'   Error in opening file integrals.dat'
    call mpi_abort(cart_comm, 1, ierr)
    stop
122 write(*,'(1x,a)')'   Error in writing in file integrals.dat'
    call mpi_abort(cart_comm, 1, ierr)
    stop   
endsubroutine integral

endmodule output_routes