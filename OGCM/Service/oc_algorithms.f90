module ocalg_routes
use constants
implicit none

contains

!=====================================================
subroutine cyclize_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines:
!do n=1,ny do k=1,nz
!  ff(mmm-1,n,k) = ff(mm ,n,k)
!  ff(mm +1,n,k) = ff(mmm,n,k)
!enddo enddo
    integer, intent(in):: mmm, mm, nx, ny, nz
    integer n, k
    real(4), intent(inout):: ff(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),     &
                          (bnd_y2 - bnd_y1 + 1)*nz,         &
                          mpi_real4, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),   &
                          (bnd_y2 - bnd_y1 + 1)*nz,         &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,             &
                          mpi_real4, dist_rank, 1,              &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,             &
                          mpi_real4, dist_rank, 1,              &
                          cart_comm, stat, ierr)
    endif

endsubroutine cyclize_x

!======================================================
subroutine cyclize8_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
!do n=1,ny do k=1,nz
!  ff(mmm-1,n,k) = ff(mm ,n,k)
!  ff(mm +1,n,k) = ff(mmm,n,k)
!enddo enddo
    integer, intent(in):: mmm, mm, nx, ny, nz
    integer n, k
    real(8), intent(inout) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
    endif

endsubroutine cyclize8_x

!=====================================================
subroutine cyclize_y(ff,nx,ny,nz,nnn,nn)
    use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
!do m=1,nx do k=1,nz
!  ff(m,nnn-1,k) = ff(m,nn ,k)
!  ff(m, nn+1,k) = ff(m,nnn,k)
!enddo enddo
    integer, intent(in):: nnn, nn, nx, ny, nz
    integer m, k
    real(4), intent(inout):: ff(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(2) .eq. 0) then
  !-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_coord(1) 
        p_dist(2) = p_size(2) - 1
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(bnd_x1:bnd_x2, nnn, 1:nz),     &
                          (bnd_x2 - bnd_x1 + 1)*nz,         &
                          mpi_real4, dist_rank, 1,          &
                          ff(bnd_x1:bnd_x2, nnn-1, 1:nz),   &
                          (bnd_x2 - bnd_x1 + 1)*nz,         &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
    endif

    if (p_coord(2) .eq. (p_size(2) - 1)) then
  !-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = 0
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(bnd_x1:bnd_x2, nn, 1:nz),          &
                          (bnd_x2 - bnd_x1 + 1)*nz,             &
                          mpi_real4, dist_rank, 1,              &
                          ff(bnd_x1:bnd_x2, nn+1, 1:nz),        &
                          (bnd_x2 - bnd_x1 + 1)*nz,             &
                          mpi_real4, dist_rank, 1,              &
                          cart_comm, stat, ierr)
    endif

endsubroutine cyclize_y

!=====================================================
subroutine cyclize8_y(ff,nx,ny,nz,nnn,nn)
    use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
!do m=1,nx do k=1,nz
!  ff(m,nnn-1,k) = ff(m,nn ,k)
!  ff(m, nn+1,k) = ff(m,nnn,k)
!enddo enddo
    integer, intent(in):: nnn, nn, nx, ny, nz
    integer m, k
    real(8), intent(inout):: ff(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(2) .eq. 0) then
  !-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_coord(1) 
        p_dist(2) = p_size(2) - 1
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(bnd_x1:bnd_x2, nnn, 1:nz),     &
                          (bnd_x2 - bnd_x1 + 1)*nz,         &
                          mpi_real8, dist_rank, 1,          &
                          ff(bnd_x1:bnd_x2, nnn-1, 1:nz),   &
                          (bnd_x2 - bnd_x1 + 1)*nz,         &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
    endif

    if (p_coord(2) .eq. (p_size(2) - 1)) then
  !-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = 0
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(bnd_x1:bnd_x2, nn, 1:nz),          &
                          (bnd_x2 - bnd_x1 + 1)*nz,             &
                          mpi_real8, dist_rank, 1,              &
                          ff(bnd_x1:bnd_x2, nn+1, 1:nz),        &
                          (bnd_x2 - bnd_x1 + 1)*nz,             &
                          mpi_real8, dist_rank, 1,              &
                          cart_comm, stat, ierr)
    endif

endsubroutine cyclize8_y

!======================================================
subroutine cyclize_int4_x(ff,nx,ny,nz,mmm,mm)
  use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
!do n=1,ny do k=1,nz
!  ff(mmm-1,n,k) = ff(mm ,n,k)
!  ff(mm +1,n,k) = ff(mmm,n,k)
!enddo enddo
  integer, intent(in):: mmm, mm, nx, ny, nz
  integer n, k
  integer(4), intent(inout):: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

  integer, dimension(2) :: p_dist
  integer :: dist_rank
  integer :: ierr
  integer stat(MPI_status_size)

  if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
      p_dist(1) = p_size(1) - 1
      p_dist(2) = p_coord(2)
      call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

      call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                        (bnd_y2 - bnd_y1 + 1)*nz,              &
                        mpi_integer4, dist_rank, 1,          &
                        ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                        (bnd_y2 - bnd_y1 + 1)*nz,              &
                        mpi_integer4, dist_rank, 1,          &
                        cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
  endif

  if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
      p_dist(1) = 0
      p_dist(2) = p_coord(2)
      call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

      call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                        (bnd_y2 - bnd_y1 + 1)*nz,              &
                        mpi_integer4, dist_rank, 1,          &
                        ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                        (bnd_y2 - bnd_y1 + 1)*nz,              &
                        mpi_integer4, dist_rank, 1,          &
                        cart_comm, stat, ierr)
  endif

endsubroutine cyclize_int4_x

!=====================================================
subroutine cyclize_int4_y(ff,nx,ny,nz,nnn,nn)
  use mpi_parallel_tools
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
!do m=1,nx do k=1,nz
!  ff(m,nnn-1,k) = ff(m,nn ,k)
!  ff(m, nn+1,k) = ff(m,nnn,k)
!enddo enddo
  integer, intent(in):: nnn, nn, nx, ny, nz
  integer m, k
  integer(4), intent(inout):: ff(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

  integer, dimension(2) :: p_dist
  integer :: dist_rank
  integer :: ierr
  integer stat(MPI_status_size)

  if (p_coord(2) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
      p_dist(1) = p_coord(1) 
      p_dist(2) = p_size(2) - 1
      call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

      call mpi_sendrecv(ff(bnd_x1:bnd_x2, nnn, 1:nz),     &
                        (bnd_x2 - bnd_x1 + 1)*nz,         &
                        mpi_integer4, dist_rank, 1,          &
                        ff(bnd_x1:bnd_x2, nnn-1, 1:nz),   &
                        (bnd_x2 - bnd_x1 + 1)*nz,         &
                        mpi_integer4, dist_rank, 1,          &
                        cart_comm, stat, ierr)
  endif

  if (p_coord(2) .eq. (p_size(2) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
      p_dist(1) = p_coord(1)
      p_dist(2) = 0
      call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

      call mpi_sendrecv(ff(bnd_x1:bnd_x2, nn, 1:nz),          &
                        (bnd_x2 - bnd_x1 + 1)*nz,             &
                        mpi_integer4, dist_rank, 1,              &
                        ff(bnd_x1:bnd_x2, nn+1, 1:nz),        &
                        (bnd_x2 - bnd_x1 + 1)*nz,             &
                        mpi_integer4, dist_rank, 1,              &
                        cart_comm, stat, ierr)
  endif

endsubroutine cyclize_int4_y
!=================================================================
  ! calculate a distance in 3-d between points on unit sphere
  function distance(lon1, lat1, lon2, lat2)
    use math_tools
    
    real(4) lon1,lat1,lon2,lat2 
    real(4) distance
    real(4) x1,y1,z1, x2,y2,z2

    z1 = sind(lat1)            
    x1 = cosd(lat1)*cosd(lon1)
    y1 = cosd(lat1)*sind(lon1)
   
    z2 = sind(lat2)
    x2 = cosd(lat2)*cosd(lon2)
    y2 = cosd(lat2)*sind(lon2)

    distance = sqrt((z1-z2)**2+(x1-x2)**2+(y1-y2)**2)
  
  end function distance

!======================================================================
! module contains the sweep procedure needed for different cases
subroutine factor(nd,a,b,c,eta,rksi,ii,jj)
!----------------------------------------------------------------------
! up-down : a(n)*f(n-1) + b(n)*f(n) + c(n)*f(n+1) = eta(n) ; n=ii,jj
! a(ii) and c(jj) are not used, so they may be set = 0.
! optimized by dianski n.a.
      integer, intent(in):: ii, jj, nd
      real(4), intent(in):: a(nd),b(nd),c(nd), eta(nd)
      real(4), intent(inout):: rksi(nd)
      real(4), allocatable:: x(:), y(:)
      integer j1, jj1
      real(4) scr

      allocate(x(nd), y(nd))

      x(ii) = - c(ii) / b(ii)
      y(ii) = eta(ii) / b(ii)
      jj1 = ii+1
      do  j1=jj1,jj
       scr=1.0/(a(j1) * x(j1-1) + b(j1))
       x(j1) = -c(j1) * scr
       y(j1) = (eta(j1) - a(j1)*y(j1-1)) * scr
      end do
      rksi(jj)=y(jj)
      jj1 = jj - 1
      do j1=jj1,ii,-1
       rksi(j1) = x(j1) * rksi(j1+1) + y(j1)
      end do
      
      deallocate(y, x)

endsubroutine factor
!======================================================================
! module contains the sweep procedure needed for different cases
subroutine factor8(nd,a,b,c,eta,rksi,ii,jj)
!----------------------------------------------------------------------
! up-down : a(n)*f(n-1) + b(n)*f(n) + c(n)*f(n+1) = eta(n) ; n=ii,jj
! a(ii) and c(jj) are not used, so they may be set = 0.
! optimized by dianski n.a.
      integer, intent(in):: ii, jj, nd
      real(8), intent(in):: a(nd),b(nd),c(nd), eta(nd)
      real(8), intent(inout):: rksi(nd)
      real(8), allocatable:: x(:), y(:)
      integer j1, jj1
      real(8) scr

      allocate(x(nd), y(nd))

      x(ii) = - c(ii) / b(ii)
      y(ii) = eta(ii) / b(ii)
      jj1 = ii+1
      do  j1=jj1,jj
       scr=1.0/(a(j1) * x(j1-1) + b(j1))
       x(j1) = -c(j1) * scr
       y(j1) = (eta(j1) - a(j1)*y(j1-1)) * scr
      end do
      rksi(jj)=y(jj)
      jj1 = jj - 1
      do j1=jj1,ii,-1
       rksi(j1) = x(j1) * rksi(j1+1) + y(j1)
      end do
      
      deallocate(y, x)

endsubroutine factor8

!================================================================
subroutine hbal_snow(hbal,  &
                     wnd,   &
                     slp,   &
                    tair,   &
                    qair,   &
                     lwd,   &
                     swd,   &
                      hi,   &
                   spart,   &
               swpen_top,   &
               swpen_bot,   &
                     kph,   &
                    tsrf,   &
                    tbot,   &
                      sh,   &
                      ev,   &
                      lh,   &
                    cdai,   &
                      lw,   &
                     swb,   &
                     swi,   &
                     ow_alb)

include 'atmforcing.fi'

real(4), intent(in):: tsrf, tbot, wnd, slp, tair, qair, lwd, swd, hi, spart, swpen_top, swpen_bot, kph, ow_alb
real(4), intent(inout):: hbal, sh, ev, lh, cdai, lw, swb, swi
real ai, as

  call air_ice_turbulent_fluxes(wnd,  &   ! wind modulo, m/s
                                slp,  &   ! sea level pressure, Pa
                               tair,  &   ! air temperature,  �C
                               tsrf,  &   ! sea surface temp, �C
                               qair,  &   ! air specific humidity, kg/kg
                          u_height,   &   ! height of wind datasets, m
                          t_height,   &   ! height of tair datasets, m
                          q_height,   &   ! height of qair datasets, m
                                 sh,  &   ! sensible heat flux, W/m^2
                                 ev,  &   ! evaporation rate, kg/m^2/s
                                 lh,  &   ! latent heat flux, W/m^2
                             cdai   )     ! drag coefficient [kg/m^2/s]

  lw = Emissice*( lwd - StephBoltz*(tsrf+273.15)**4 )
  
  ai=albedo_ice(tbot, hi, ow_alb)
  as=albedo_snow(tsrf)
  
  swb=swd * (spart*(1.0-as) + (1.0-spart)*(1.0-ai)*(1.0-swpen_top))
  swi=swd * (1.0-spart)*(1.0-ai)*swpen_bot

  hbal = sh + lh + lw + swb + kph*(tbot-tsrf)

endsubroutine hbal_snow

!====================================================================
subroutine  hbal_ice(hbal,  &
                     wnd,   &
                     slp,   &
                    tair,   &
                    qair,   &
                     lwd,   &
                     swd,   &
                      hi,   &
               swpen_top,   &
               swpen_bot,   &
                     kph,   &
                    tsrf,   &
                    tbot,   &
                      sh,   &
                      ev,   &
                      lh,   &
                    cdai,   &
                      lw,   &
                     swb,   &
                     swi,   &
                     ow_alb)

include 'atmforcing.fi'

real(4), intent(in):: tsrf, tbot, wnd, slp, tair, qair, lwd, swd, hi, swpen_top, swpen_bot, kph, ow_alb
real(4), intent(inout):: hbal, sh, ev, lh, cdai, lw, swb, swi

real(4) ai

  call air_ice_turbulent_fluxes(wnd,  &   ! wind modulo, m/s
                                slp,  &   ! sea level pressure, Pa
                               tair,  &   ! air temperature,  �C
                               tsrf,  &   ! sea surface temp, �C
                               qair,  &   ! air specific humidity, kg/kg
                          u_height,   &   ! height of wind datasets, m
                          t_height,   &   ! height of tair datasets, m
                          q_height,   &   ! height of qair datasets, m
                                 sh,  &   ! sensible heat flux, W/m^2
                                 ev,  &   ! evaporation rate, kg/m^2/s
                                 lh,  &   ! latent heat flux, W/m^2
                             cdai   )     ! drag coefficient [kg/m^2/s]

  lw = Emissice*(lwd-StephBoltz*(tsrf+273.15)**4)

  ai=albedo_ice(tsrf, hi, ow_alb)
  
  swb=swd * (1.0-ai)*(1.0-swpen_top)
  swi=swd * (1.0-ai)*swpen_bot

  hbal = sh + lh + lw + swb + kph*(tbot-tsrf)

endsubroutine hbal_ice

!==========ice albedo as a function of temperature and height===========
function albedo_ice(tice, hi, ow_alb)

 real(4) albedo_ice, tice, hi, ow_alb
 real(4) tem

  tem=max(min(tice,0.0),-1.0)

  albedo_ice = ow_alb + (alb_ice_dry-alb_ice_decr*(tem+1.0)-ow_alb) * min(hi,hice_max_alb)/hice_max_alb

endfunction albedo_ice

!==========snow albedo as a function of temperature===========
function albedo_snow(tsnow)

 real(4) albedo_snow, tsnow
 real(4) tem
 
  tem=max(min(tsnow,0.0),-1.0)
 
  albedo_snow = alb_snow_dry - alb_snow_decr * (tem+1.0)

endfunction albedo_snow

!========program for interpolation from z-horizons to s-levels========
subroutine z2s(f1,      &     !input  3d field (on z-levels)
               f2,      &     !output 3d field (on s-levels)
              hhq,      &     !2d field of bottom topography (in metres)
              msk,      &     !temperature mask
           zsigma,      &     !array of s-levels
             zlvl,      &     !array of z-levels (in metres)
    bnd_x1,bnd_x2,      &     !dimension on x
    bnd_y1,bnd_y2,      &     !dimension on y
               nz,      &     !number of s-levels
             nlvl,      &     !number of z-levels
             lev1,      &     !parameter of task
             over,      &     !undefined value
             ierr)            !error index

! if lev1=1 then first level of f1 is definited as fist level of f1
! if lev1=0 then first level of f1 is definited correspondly to its level position

integer,intent(in):: bnd_x1,bnd_x2, bnd_y1,bnd_y2, nz, nlvl

real(4),intent(in):: f1(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlvl),       &
                    hhq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),            &
                    msk(bnd_x1:bnd_x2,bnd_y1:bnd_y2),            &
                    zsigma(nz), zlvl(nlvl)

real(4),intent(inout):: f2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

real(4) fz1,fz2,z1,z2,zs(150),fz(150),ds,deep
real(4),intent(in):: over
real(4) over5
integer i,j,k,kdeep,kups,nocegrid,lev1,ierr, kr

! over5 - 5% vicinity to over
  over5=abs(over)-abs(over)/20.0
  ierr=0
  write(*,*)' interpolation from z-levels to s-levels'
      
  if(nz>150) then
     write(*,*)' error in routine z2s:'
     write(*,*)' number of s-levels is greater than 150'
     ierr=1
     return
  endif
  
  if(nlvl>150) then
     write(*,*)' error in routine z2s:'
     write(*,*)' number of z-levels is greater than 150'
     ierr=1
     return
  end if

!  interpolation
  nocegrid=0
  
  do j=bnd_y1,bnd_y2
   do i=bnd_x1,bnd_x2
    if(msk(i,j)>0.5) then
      nocegrid=nocegrid+1
      deep=hhq(i,j)
      
      do k=1,nz
         zs(k)=zsigma(k)*deep
      end do
     
      fz(1)=f1(i,j,1)
      
      if(abs(fz(1))>over5) then
         write(*,*)' error in routine z2s:'
         write(*,'(a,2i4,a)') ' in point ',i,j,' input value of upper level is undefined!'
         ierr=1
         return
      end if

!--------- making profile without bottom-------------
      ierr=0
      
      do k=2,nlvl
       fz(k)=f1(i,j,k)
       if(abs(fz(k))>over5) then
         fz(k)=fz(k-1)
         ierr=ierr+1
       end if
      end do

      if(ierr==nlvl-1) then
        write(*,*)' warning in routine z2s:'
        write(*,1000) i,j,deep,zlvl(1)
      end if
      
      ierr=0

!  searching number of upper sigma level
      kups=1
 
       do while(zs(kups)<zlvl(1).and.kups<nz)
        kups=kups+1
       end do
 
      kups=kups-1
 
      if(kups>=1) then

!   for upper sigma levels of f2
         do k=1,kups
           f2(i,j,k)=f1(i,j,1)
         enddo
    
      else
         if (lev1/=0) then
          f2(i,j,1)=f1(i,j,1)
          kups=1
         end if
      end if


! searching the deepest z level
      kdeep=1
      
      do while(zlvl(kdeep)<deep.and.kdeep<nlvl)
       kdeep=kdeep+1
      end do
      
      kr=1
      fz1=fz(1)
      fz2=fz1
      z1=zlvl(1)-zlvl(1)/20.0
      z2=zlvl(1)

      do k=kups+1,nz
       
       ds=zs(k)
       
       do while((ds<=z1.or.ds>z2).and.kr<kdeep)
          kr=kr+1
          fz1=fz2
          z1=z2
          fz2=fz(kr)
          z2=zlvl(kr)
       end do
       
       if(ds>=zlvl(kdeep)) then
         f2(i,j,k)=fz(kdeep)
       else
         f2(i,j,k)=(fz1*(z2-ds)+fz2*(ds-z1))/(z2-z1)
       end if

      end do

      end if
    enddo
 enddo

 write(*,*) ' for control: number of ocean grid in interpolation is ',nocegrid
 return

1000  format(' in point i=',i5,',',' j=',i5,' deep=',f9.2,'m'/    &
            ' there is only one level for interpolation'/         &
            ' at z-level of',f9.2,'m')

endsubroutine z2s
!======================================================================
subroutine s2z(f1,      &       !input  3d field (on s-levels)
               f2,      &       !output 3d field (on z-levels)
              hhq,      &       !2d field of bottom topography (in metres),
              msk,      &       !temperature mask
           zsigma,      &       !array of s-levels
             zlvl,      &       !array of z-levels (in metres)
    bnd_x1,bnd_x2,      &       !dimension on x     !dimension on x
    bnd_y1,bnd_y2,      &       !dimension on y     !dimension on y
               nz,      &       !number of s-levels
             nlvl,      &       !number of z-levels
             lev1,      &       !parameter of task
             over)              !undefined value
            
! Ь diansky n.a. 18.03.98 15:43
! program for interpolation from sigma levels on common levels
! f1 - first  3-d field (input data on sigma levels)
! f2 - second 3-d field (output data on common levels)
! hhh - 2-d field of bottom topography (IN CENTIMETERS!!!)
! msk - mask of ocean grids
! zsigma - values of sigma levels
! zlvl - values of common levels (IN METERS!!!)
! nx,ny - x,y dimension
! nz - number of sigma levels
! nlvl - number of common levels
! lev1 - parameter of task:
! if lev1=1 then first level of f2 is definited as fist level of f1
! if lev1=0 then first level of f2 is definited correspondly to its
!           level position
! over - undefinite value
   
      integer,intent(in)::  bnd_x1,bnd_x2,bnd_y1,bnd_y2, nz, nlvl, lev1
      integer k, i, j, kr, kdeep, kshal

      real(4),intent(in):: f1(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),        &
                          hhq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &
                          msk(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &
                          zsigma(nz),zlvl(nlvl)

      real(4),intent(inout):: f2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlvl)
      real(4) fz1,fz2,zs1,zs2,zs(nz), dlvl, zdeep, deep
      real(4),intent(in):: over

! interpolation
      do j=bnd_y1,bnd_y2
         do i=bnd_x1,bnd_x2
           if(msk(i,j)>0.5) then
           deep=hhq(i,j)
           do k=1,nz
              zs(k)=zsigma(k)*deep
           end do
!  finding levels nearest to shalow and deep sigma levels
           zdeep=1.0e+10
           kshal=0
           do k=1,nlvl
                 dlvl=zlvl(k)
              if(dlvl<=zs(1)) then
                 kshal=k
              end if
              if(abs(zs(nz)-dlvl)<=zdeep) then
                 kdeep=k
                 zdeep=abs(zs(nz)-dlvl)
              end if
           end do

           if(kshal>=1) then
!   for upper layers of f2
              do k=1,kshal
                 f2(i,j,k)=f1(i,j,1)
              enddo
           else
              if (lev1/=0) then
                 f2(i,j,1)=f1(i,j,1)
                 kshal=1
              end if
           end if

           kr=1
           fz1=f1(i,j,1)
           fz2=fz1
           zs1=zs(1)-zs(1)/20.0
           zs2=zs(1)

           do k=kshal+1,kdeep
              dlvl=zlvl(k)
              do while((dlvl<=zs1.or.dlvl>zs2).and.kr<nz)
                 kr=kr+1
                 fz1=fz2
                 zs1=zs2
                 fz2=f1(i,j,kr)
                 zs2=zs(kr)
              end do

              if(dlvl>=zs(nz)) then
                f2(i,j,k)=f1(i,j,nz)
              else
                f2(i,j,k)=(fz1*(zs2-dlvl)+fz2*(dlvl-zs1))/(zs2-zs1)
              end if
           end do

           if(kdeep<nlvl) then
              do k=kdeep+1,nlvl
              f2(i,j,k)=over
              enddo
           end if

           end if
         enddo
      enddo

endsubroutine s2z

!=======10 m neutral stability bulk-coefficients
subroutine bulk_n10_ocean(u10,      &  !10m wind
                         zeta,      &  !stability function
                           cd,      &  !drag coefficient
                           cd_rt,   &  !drag coeff. sq.root
                           ch,      &  !heat coefficient
                           ce)         !vapour coefficient
real(4),intent(in):: u10, zeta
real(4),intent(inout):: cd, cd_rt, ch, ce
real(4) stab, wnd

 wnd= max(min(u10,33.0),1.0)

 cd = ( 2.7/wnd + 0.142 + wnd/13.09 - 3.14807e-10*wnd**6 )/1000.0
 cd_rt = sqrt(cd) 

 ce    = cd_rt * 34.6/1000.0                        ! L-Y eqn. 6b again

 stab  = 0.5 + sign(0.5,zeta)
 ch    = cd_rt*(18.0*stab+32.7*(1.0-stab))/1000.0   ! L-Y eqn. 6c again  

endsubroutine bulk_n10_ocean

!=============================================================================
subroutine air_sea_turbulent_fluxes(wnd,        &   ! wind modulo, m/s
                                    slp,        &   ! sea level pressure, Pa
                                    tair,       &   ! air temperature,  °C
                                    tsea,       &   ! sea surface temp, °C
                                    ssea,       &   ! sea surface salinity, psu
                                    qair,       &   ! air specific humidity, kg/kg
                                    u_hgt,      &   ! height of wind datasets, m
                                    t_hgt,      &   ! height of tair datasets, m
                                    q_hgt,      &   ! height of qair datasets, m
                                    sens_heat,  &   ! sensible heat flux, W/m^2
                                    evap_rate,  &   ! evaporation rate, kg/m^2/s
                                    lat_heat,   &   ! latent heat flux, W/m^2
                                    cdao   )        ! drag coefficient [kg/m^2/s]      
  
 real(4),intent(in):: wnd, slp, tair, tsea, ssea, qair, u_hgt, t_hgt, q_hgt  !input data
 real(4),intent(inout):: sens_heat, evap_rate, lat_heat, cdao                   !output data
 
 integer, parameter:: niter=2
 integer iter

 real(4) rho_a, tv, qsat, qv, psat, cpa, lhv
 real(4) u, u10, t_zu, q_zu, tstar, qstar, xx, stab
 real(4) cd, ch, ce, ustar, bstar   
 real(4) cd_n10, ce_n10, ch_n10, cd_n10_rt   ! neutral 10m drag coefficients                         
 real(4) zeta_u, zeta_t, zeta_q, x2, x,  &
             psi_m_u, psi_m_t, psi_m_q,  &
             psi_h_u, psi_h_t, psi_h_q       ! stability parameters
 
     t_zu = tair
     q_zu = qair

     u = max(wnd, 1.0)         ! 0.5 m/s floor on wind
     u10 = u

 call bulk_n10_ocean(u10,      &  !10m wind
               t_zu-tsea,      &  !stability function
                  cd_n10,      &  !drag coefficient
                  cd_n10_rt,   &  !drag coeff. sq.root
                  ch_n10,      &  !heat coefficient
                  ce_n10)   

! first guess for exchange coeff's at z
      cd    = cd_n10                                          
      ch    = ch_n10
      ce    = ce_n10

      tv    = virt_temp(t_zu,q_zu)
      qv    = virt_hum(q_zu)
      psat  = svp_as(tsea,ssea)
      qsat  = spec_hum(slp,psat)

      do iter=1,niter

       call atm_turb_scales(cd, ch, ce, u,       &
                            t_zu, tsea, tv,      &
                            q_zu, qsat, qv,      & 
                            ustar, tstar,        & 
                            qstar, bstar  )    

! stability for U-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         u_hgt,     &  !heihgt of the variable
                        zeta_u,     &  !stability function
                       psi_m_u,     &  !momentum stability 
                       psi_h_u)        !heat stability
 
 ! stability for T-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         t_hgt,     &  !heihgt of the variable
                        zeta_t,     &  !stability function
                       psi_m_t,     &  !momentum stability 
                       psi_h_t)        !heat stability
! stability for Q-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         q_hgt,     &  !heihgt of the variable
                        zeta_q,     &  !stability function
                       psi_m_q,     &  !momentum stability 
                       psi_h_q)        !heat stability
      
        u10 = u/(1.0+cd_n10_rt*(log(u_hgt/10.0)-psi_m_u)/vonkarman) ! L-Y eqn. 9
        t_zu= tair  - tstar/vonkarman*(log(t_hgt/u_hgt) +psi_h_u - psi_h_t )
        q_zu= qair  - qstar/vonkarman*(log(q_hgt/u_hgt) +psi_h_u - psi_h_q )      

        call bulk_n10_ocean(u10,      &  !10m wind
                         zeta_t,      &  !stability function
                         cd_n10,      &  !drag coefficient
                         cd_n10_rt,   &  !drag coeff. sq.root
                         ch_n10,      &  !heat coefficient
                         ce_n10)  

        xx = (log(u_hgt/10.0)-psi_m_u)/vonkarman 
        cd = cd_n10/(1.0+cd_n10_rt*xx)**2                        ! L-Y 10a

        xx = (log(u_hgt/10.0)-psi_h_u)/vonkarman 
        ch = ch_n10/(1.0+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10b (corrected code aug2007)
        ce = ce_n10/(1.0+ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10c (corrected code aug2007)

        tv    = virt_temp(t_zu,q_zu)
        qv    = virt_hum(q_zu)
      
      end do

      rho_a = air_den(tv,slp)
        cpa = air_spec_heat(q_zu)
        lhv = lat_heat_vapor(tsea)

      sens_heat = cpa*rho_a*ch*u*(t_zu-tsea)
      evap_rate =     rho_a*ce*u*(q_zu-qsat)
      lat_heat  = lhv*evap_rate
      cdao = rho_a*cd*u

endsubroutine air_sea_turbulent_fluxes

!=============================================================================
subroutine air_ice_turbulent_fluxes(wnd,        &   ! wind modulo, m/s
                                    slp,        &   ! sea level pressure, Pa
                                    tair,       &   ! air temperature,  °C
                                    tsrf,       &   ! sea surface temp, °C
                                    qair,       &   ! air specific humidity, kg/kg
                                    u_hgt,      &   ! height of wind datasets, m
                                    t_hgt,      &   ! height of tair datasets, m
                                    q_hgt,      &   ! height of qair datasets, m
                                    sens_heat,  &   ! sensible heat flux, W/m^2
                                    subl_rate,  &   ! evaporation rate, kg/m^2/s
                                    lat_heat,   &   ! latent heat flux, W/m^2
                                    cdai   )        ! drag coefficient [kg/m^2/s]    
  
  real(4),intent(in):: wnd, slp, tair, tsrf, qair, u_hgt, t_hgt, q_hgt       !input data
  real(4),intent(inout):: sens_heat, subl_rate, lat_heat, cdai                  !output data
  
  integer, parameter:: niter=5
  integer iter

  real(4), parameter:: cd_n10_i    = 1.63e-3,           &
                       cd_n10_i_rt = sqrt(cd_n10_i),    &
                       ce_n10_i    = 1.63e-3,           &
                       ch_n10_i    = 1.63e-3    

  real(4) rho_a, tv, qsat, qv, psat, cpa, lhs
  real(4) u, t_zu, q_zu, tstar, qstar, xx, stab
  real(4) cd, ch, ce, ustar, bstar                             
  real(4) zeta_u, zeta_t, zeta_q, x2, x,  &
              psi_m_u, psi_m_t, psi_m_q,  &
              psi_h_u, psi_h_t, psi_h_q       ! stability parameters
  
      t_zu = tair
      q_zu = qair

      u = max(wnd, 1.0)         ! 0.5 m/s floor on wind

! first guess for exchange coeff's at z
      cd    = cd_n10_i                                          
      ch    = ch_n10_i
      ce    = ce_n10_i

      tv    = virt_temp(t_zu,q_zu)
      qv    = virt_hum(q_zu)
      psat  = svp_ai(tsrf)
      qsat  = spec_hum(slp,psat)

      do iter=1,niter

       call atm_turb_scales(cd, ch, ce, u,       &
                            t_zu, tsrf, tv,      &
                            q_zu, qsat, qv,      & 
                            ustar, tstar,        & 
                            qstar, bstar  )    

! stability for U-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         u_hgt,     &  !heihgt of the variable
                        zeta_u,     &  !stability function
                       psi_m_u,     &  !momentum stability 
                       psi_h_u)        !heat stability
 
 ! stability for T-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         t_hgt,     &  !heihgt of the variable
                        zeta_t,     &  !stability function
                       psi_m_t,     &  !momentum stability 
                       psi_h_t)        !heat stability
! stability for Q-height
      call atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                         q_hgt,     &  !heihgt of the variable
                        zeta_q,     &  !stability function
                       psi_m_q,     &  !momentum stability 
                       psi_h_q)        !heat stability
      
        t_zu= tair  - tstar/vonkarman*(log(t_hgt/u_hgt) +psi_h_u - psi_h_t )
        q_zu= qair  - qstar/vonkarman*(log(q_hgt/u_hgt) +psi_h_u - psi_h_q )      

        xx = (log(u_hgt/10.0)-psi_m_u)/vonkarman 
        cd = cd_n10_i/(1.0+cd_n10_i_rt*xx)**2                        ! L-Y 10a

        xx = (log(u_hgt/10.0)-psi_h_u)/vonkarman 
        ch = ch_n10_i/(1.0+ch_n10_i*xx/cd_n10_i_rt)*sqrt(cd/cd_n10_i) ! 10b (corrected code aug2007)
        ce = ce_n10_i/(1.0+ce_n10_i*xx/cd_n10_i_rt)*sqrt(cd/cd_n10_i) ! 10c (corrected code aug2007)

        tv    = virt_temp(t_zu,q_zu)
        qv    = virt_hum(q_zu)
      
      end do

      rho_a = air_den(tv,slp)
        cpa = air_spec_heat(q_zu)
        lhs = lat_heat_subl(tsrf)

      sens_heat = cpa*rho_a*ch*u*(t_zu-tsrf)
      subl_rate =     rho_a*ce*u*(q_zu-qsat)
      lat_heat  = lhs*subl_rate
      cdai = rho_a*cd*u

endsubroutine air_ice_turbulent_fluxes

!=============================================================
subroutine atm_turb_scales(cd,      &  !drag   coefficient
                           ch,      &  !heat   coefficient
                           ce,      &  !vapour coefficient
                            u,      &  !wind speed
                           ta,      &  !air temperature
                           ts,      &  !surface temperature
                           tv,      &  !virtual temperature
                           qa,      &  !air humidity
                           qs,      &  !surface humidity
                           qv,      &  !virtual humidity
                           ustar,   &  !turbulent u-scale
                           tstar,   &  !turbulent t-scale
                           qstar,   &  !turbulent q-scale
                           bstar  )    !turbulent bouyancy scale

real(4),intent(in):: cd, ch, ce, u, ta, ts, tv, qa, qs, qv
real(4),intent(inout):: ustar, tstar, qstar, bstar
real(4) cd_rt
  
  cd_rt = sqrt(cd)
  ustar = cd_rt*u                                ! L-Y eqn. 7a
  tstar = (ch/cd_rt)*(ta-ts)                     ! L-Y eqn. 7b
  qstar = (ce/cd_rt)*(qa-qs)                     ! L-Y eqn. 7c
  bstar    = FreeFallAcc*(tstar/tv+qstar/qv)   

endsubroutine atm_turb_scales

!=====================================================================
subroutine atm_stability(bstar,     &  !bouyancy turb. scale
                         ustar,     &  !wind turbylent scale
                           hgt,     &  !heihgt of the variable
                          zeta,     &  !stability function
                         psi_m,     &  !momentum stability 
                         psi_h)        !heat stability

real(4),intent(in):: bstar, ustar, hgt
real(4),intent(inout)::zeta, psi_m, psi_h

real(4) x2, x

! stability for the given height
        zeta   = vonkarman*bstar*hgt/(ustar*ustar)      ! L-Y eqn. 8a
        zeta   = max(min(zeta,10.0),-10.0)              ! undocumented NCAR
    
        if (zeta > 0.0) then
          psi_m = -5.0*zeta                             ! L-Y eqn. 8c
          psi_h = -5.0*zeta                             ! L-Y eqn. 8c
        else
          x2 = sqrt(1.0-16.0*zeta)                     ! L-Y eqn. 8b
          x = sqrt(x2);
          psi_m = log( (1.0+2.0*x+x2)*(1.0+x2)/8.0 ) - 2.0*( atan(x)-atan(1.0) ) ! L-Y eqn. 8d
          psi_h = 2.0*log((1.0+x2)/2.0)                 ! L-Y eqn. 8e
        end if

endsubroutine atm_stability

!========Specific humidity from air and vapour pressures===================================
function spec_hum(p_air, p_vap)
 
 real(4) p_air,    &   !air pressure(Pa)
         p_vap,    &   !vapour pressure(Pa)
         spec_hum      !specific humidity (undim)

  spec_hum = mw2ma*p_vap/(p_air-(1.0-mw2ma)*p_vap)

endfunction spec_hum

!=======Saturated vapour pressure at air-sea interface==================
function svp_as(sst,sss)
 
 real(4) sst,          &   !Surface temperature (°C)
         sss,          &   !Sea surface salinity (°C)
         svp_as            !Saturated vapour pressure at air-sea (Pa)
 
 real(4) pow, fs, t, s
 
 t = max(min(sst,40.0),-40.0)
 s = max(min(sss,35.0),  0.0)/35.0
 pow = (0.7859 + 0.03477*t)/(1.0 + 0.00412*t)   ! power index
 fs=  1.0 - 0.02*s
 
 svp_as = fs * 10.0**(pow + 2.0)

endfunction svp_as

!=======Saturated vapour pressure at air-ice interface==================
function svp_ai(tsurf)

 real(4) tsurf,          &   !Surface temperature (°C)  
         svp_ai              !Saturated vapour pressure at air-ice (Pa)

 real(4) pow, t
 
 t = max(min(tsurf,40.0),-40.0)
 pow = (0.7859 + 0.03477*t)/(1.0 + 0.00412*t) + 0.00422*min(t,0.0)   ! power index
 svp_ai = 10.0**(pow + 2.0)

endfunction svp_ai

!=====Saturated vapour factor for air========================
function svp_air_factor(tem,p_air)
 
 real(4) tem,          &   !Air temperature (°C)
         p_air,        &   !Air pressure(Pa)
         svp_air_factor
 real(4) t

 t = max(min(tem,40.0),-40.0)
 svp_air_factor = 1.0 + 1.0e-8*p_air*(4.5 + 0.0006*t*t)      

endfunction svp_air_factor

!=======Air specific heat
function air_spec_heat(q_air)
real(4) q_air,           &   !Air specific humidity(kg/kg)
        air_spec_heat        !Air specific heat (J/kg/K)
 
 air_spec_heat = 3.5 * rgas_air * (1.0 - q_air + 8.0*q_air/7.0/mw2ma)

endfunction air_spec_heat

!========Latent heat of vaporisation
function lat_heat_vapor(tsurf)
real(4) tsurf,    &      !Surface temperature (°C)
        lat_heat_vapor   !Latent heat of vaporisation (J/kg)
real(4) t
 
 t = max(min(tsurf,40.0),-2.0)

 lat_heat_vapor= 2.5008e+06 - 2.3e+03*t

endfunction lat_heat_vapor

!========Latent heat of sublimation
function lat_heat_subl(tsurf)
real(4) tsurf,    &     !Surface temperature (°C)
        lat_heat_subl   !Latent heat of sublimation (J/kg)
real(4) t
 
 t = max(min(tsurf,0.0),-35.0)

 lat_heat_subl= 2.839e+06 - 3.6*(t+35.0)**2

endfunction lat_heat_subl

!=======virtual temperature
function virt_temp(t_air,q_air)
real(4) t_air,    &     !Air temperature (°C)
        q_air,    &     !Air specific humidity (kg/kg)
     virt_temp          !Air virtual temperature (K)

 virt_temp = (t_air + 273.15) * (1.0 - q_air + q_air/mw2ma)

endfunction virt_temp

!=======virtual humidity
function virt_hum(q_air)
real(4) q_air,    &     !Air specific humidity (kg/kg)
      virt_hum          !Air virtual humidity (kg/kg)

 virt_hum = q_air + mw2ma/(1.0 - mw2ma)

endfunction virt_hum

!=======air density
function air_den(virt_temp,p_air)
real(4) virt_temp,   &     !Air virtual temperature (K)
           p_air,    &     !Air pressure (Pa)
         air_den           !Air density (kg/m^3)

  air_den = p_air/rgas_air/virt_temp

endfunction air_den

!=======sea water freesing point as a function of salinity===========
function freeze_temp(sal)
 real(4) sal, freeze_temp, s
    s=max(sal,0.0)
    freeze_temp= (-0.0575 + 1.710523e-3 * sqrt(s) - 2.154996e-4 * s) *s
endfunction freeze_temp

!==========================================================================
function density_brydon_narrow_subtracted(theta,sal,pres)
! Sea water density deviation from density 1000.0 [kg/m**3]
! as a function of theta[°C](potential),sal[PSU], pres[Pa]
! By Brydon et al.:"A new approximation of the equation of state for
!                   seawater, suitable for numerical ocean models."
! In: J.Geophis.Res.,v.104,No.C1, p.1537-1540, 1999.
! Pressure profile is subtracted
!
! Pressure in equations is in MPa (1.0e+06 Pa)
!
! Case 1: More accurate, but in a less wide variable range:
! -2<theta<30 °C; 30<s<38 PSU; 0<p<50 MPa

real(4), parameter:: a1=-1.36471e-01, b1= 5.06423e-01, c1=-5.52640e-04,       &
                     a2= 4.68181e-02, b2=-3.57109e-03, c2= 4.88584e-06,       &
                     a3= 8.07004e-01, b3=-8.76148e-04, c3= 9.96027e-07,       &
                     a4=-7.45353e-03, b4= 5.25243e-05, c4=-7.25139e-08,       &
                     a5=-2.94418e-03, b5= 1.57976e-05, c5=-3.98736e-09,       &
                     a6= 3.43570e-05, b6=-3.46686e-07, c6= 4.00631e-10,       &
                     a7= 3.48658e-05, b7=-1.68764e-07, c7= 8.26368e-11
      
      real(4) theta, sal, pres, density_brydon_narrow_subtracted
      real(4) t, s, p
      real(4) f1,f2,f3,f4,f5,f6,f7
             
      t=min(max(theta,-2.0),50.0)
      s=min(max(sal,0.0),50.0)
      p=max(pres*1.0e-06, 0.0)  !Transition of Pa to MPa

!     f1 = a1 + p*(b1 + c1*p)
      f1= 0.0 
      f2 = a2 + p*(b2 + c2*p)
      f3 = a3 + p*(b3 + c3*p)
      f4 = a4 + p*(b4 + c4*p)
      f5 = a5 + p*(b5 + c5*p)
      f6 = a6 + p*(b6 + c6*p)
      f7 = a7 + p*(b7 + c7*p)

      density_brydon_narrow_subtracted =        &
                f1 + t*(f2 + t*(f4 + f6*t))    &
                   + s*(f3 + t*(f5 + f7*t))
!A check value is density_brydon_narrow( 1°C, 35 PSU, 50*1e+6 Pa) =  50.402782 kg/m^3
endfunction density_brydon_narrow_subtracted

!==========================================================================
function density_brydon_wide_subtracted(theta,sal,pres)
! Sea water density deviation from density 1000.0 [kg/m**3]
! as a function of theta[°C](potential),sal[PSU], pres[Pa]
! By Brydon et al.:"A new approximation of the equation of state for
!                   seawater, suitable for numerical ocean models."
! In: J.Geophis.Res.,v.104,No.C1, p.1537-1540, 1999.
! Pressure profile is subtracted
!
! Pressure in equations is in MPa (1.0e+06 Pa)

! Case 2: Less accurate, but in a more wide variable range:
! -2<theta<40 °C; 0<s<42 PSU; 0<p<100 MPa                    
real(4), parameter:: a1=-9.20601e-02, b1= 5.07043e-01, c1=-5.43283e-04,       &
                     a2= 5.10768e-02, b2=-3.69119e-03, c2= 6.54837e-06,       &
                     a3= 8.05999e-01, b3=-9.34012e-04, c3= 1.38777e-06,       &
                     a4=-7.40849e-03, b4= 5.33243e-05, c4=-1.01563e-07,       &
                     a5=-3.01036e-03, b5= 1.75145e-05, c5=-2.34892e-08,       &
                     a6= 3.32267e-05, b6=-3.25887e-07, c6= 4.98612e-10,       &
                     a7= 3.21931e-05, b7=-1.65849e-07, c7= 2.17612e-10
      
      real(4) theta, sal, pres, refden, density_brydon_wide_subtracted
      real(4) t, s, p
      real(4) f1,f2,f3,f4,f5,f6,f7
             
      t=min(max(theta,-2.0),50.0)
      s=min(max(sal,0.0),50.0)
      p=max(pres*1.0e-06, 0.0)  !Transition of Pa to MPa

!     f1 = a1 + p*(b1 + c1*p)
      f1= 0.0 
      f2 = a2 + p*(b2 + c2*p)
      f3 = a3 + p*(b3 + c3*p)
      f4 = a4 + p*(b4 + c4*p)
      f5 = a5 + p*(b5 + c5*p)
      f6 = a6 + p*(b6 + c6*p)
      f7 = a7 + p*(b7 + c7*p)

      density_brydon_wide_subtracted =          &
                f1 + t*(f2 + t*(f4 + f6*t))    &
                   + s*(f3 + t*(f5 + f7*t))

endfunction density_brydon_wide_subtracted

!==========================================================================
function density_brydon_narrow_full(theta,sal,pres)
implicit none
! Sea water density deviation from density 1000.0 [kg/m**3]
! as a function of theta[°C](potential),sal[PSU], pres[Pa]
! By Brydon et al.:"A new approximation of the equation of state for
!                   seawater, suitable for numerical ocean models."
! In: J.Geophis.Res.,v.104,No.C1, p.1537-1540, 1999.
! Pressure profile is subtracted
!
! Pressure in equations is in MPa (1.0e+06 Pa)
!
! Case 1: More accurate, but in a less wide variable range:
! -2<theta<30 °C; 30<s<38 PSU; 0<p<50 MPa

real(4), parameter:: a1=-1.36471e-01, b1= 5.06423e-01, c1=-5.52640e-04,       &
                     a2= 4.68181e-02, b2=-3.57109e-03, c2= 4.88584e-06,       &
                     a3= 8.07004e-01, b3=-8.76148e-04, c3= 9.96027e-07,       &
                     a4=-7.45353e-03, b4= 5.25243e-05, c4=-7.25139e-08,       &
                     a5=-2.94418e-03, b5= 1.57976e-05, c5=-3.98736e-09,       &
                     a6= 3.43570e-05, b6=-3.46686e-07, c6= 4.00631e-10,       &
                     a7= 3.48658e-05, b7=-1.68764e-07, c7= 8.26368e-11
      
      real(4) theta, sal, pres, density_brydon_narrow_full
      real(4) t, s, p
      real(4) f1,f2,f3,f4,f5,f6,f7
       
      t=min(max(theta,-2.0),50.0)
      s=min(max(sal,0.0),50.0)
      p=max(pres*1.0e-06, 0.0)  !Transition of Pa to MPa

      f1 = a1 + p*(b1 + c1*p)
      f2 = a2 + p*(b2 + c2*p)
      f3 = a3 + p*(b3 + c3*p)
      f4 = a4 + p*(b4 + c4*p)
      f5 = a5 + p*(b5 + c5*p)
      f6 = a6 + p*(b6 + c6*p)
      f7 = a7 + p*(b7 + c7*p)

      density_brydon_narrow_full =        &
          f1 + t*(f2 + t*(f4 + f6*t))    &
             + s*(f3 + t*(f5 + f7*t))
!A check value is density_brydon_narrow( 1°C, 35 PSU, 50*1e+6 Pa) =  50.402782 kg/m^3

endfunction density_brydon_narrow_full

!==========================================================================
function density_brydon_wide_full(theta,sal,pres)
implicit none
! Sea water density deviation from density 1000.0 [kg/m**3]
! as a function of theta[°C](potential),sal[PSU], pres[Pa]
! By Brydon et al.:"A new approximation of the equation of state for
!                   seawater, suitable for numerical ocean models."
! In: J.Geophis.Res.,v.104,No.C1, p.1537-1540, 1999.
! Pressure profile is subtracted
!
! Pressure in equations is in MPa (1.0e+06 Pa)

! Case 2: Less accurate, but in a more wide variable range:
! -2<theta<40 °C; 0<s<42 PSU; 0<p<100 MPa                    
real(4), parameter:: a1=-9.20601e-02, b1= 5.07043e-01, c1=-5.43283e-04,       &
                     a2= 5.10768e-02, b2=-3.69119e-03, c2= 6.54837e-06,       &
                     a3= 8.05999e-01, b3=-9.34012e-04, c3= 1.38777e-06,       &
                     a4=-7.40849e-03, b4= 5.33243e-05, c4=-1.01563e-07,       &
                     a5=-3.01036e-03, b5= 1.75145e-05, c5=-2.34892e-08,       &
                     a6= 3.32267e-05, b6=-3.25887e-07, c6= 4.98612e-10,       &
                     a7= 3.21931e-05, b7=-1.65849e-07, c7= 2.17612e-10
      
      real(4) theta, sal, pres, refden, density_brydon_wide_full
      real(4) t, s, p
      real(4) f1,f2,f3,f4,f5,f6,f7

      t=min(max(theta,-2.0),50.0)
      s=min(max(sal,0.0),50.0)
      p=max(pres*1.0e-06, 0.0)  !Transition of Pa to MPa

      f1 = a1 + p*(b1 + c1*p)
 !    f1= 0.0 
      f2 = a2 + p*(b2 + c2*p)
      f3 = a3 + p*(b3 + c3*p)
      f4 = a4 + p*(b4 + c4*p)
      f5 = a5 + p*(b5 + c5*p)
      f6 = a6 + p*(b6 + c6*p)
      f7 = a7 + p*(b7 + c7*p)

      density_brydon_wide_full =          &
          f1 + t*(f2 + t*(f4 + f6*t))    &
             + s*(f3 + t*(f5 + f7*t))

endfunction density_brydon_wide_full

!====Potential temperature as a function of T,S,P ======================
function potential_temperature(tem,sal,pres)
implicit none

! potential temperature as a function of tem[°C](in situ), sal[PSU], pres[Pa]
!Salinity in equation is in deviation from 35 PSU
!Pressure in equation is in bar (1.0e+05 Pa)
real(4), parameter:: a001=-3.6504e-04, a101=-8.3198e-05, a201= 5.4065e-07, a301=-4.0274e-09,    &
                     a011=-1.7439e-05, a111= 2.9778e-07,                                        &
                     a002=-8.9309e-07, a102= 3.1628e-08, a202=-2.1987e-10, a012= 4.1057e-09,    &
                     a003= 1.6056e-10, a103=-5.0484e-12

real(4) tem, sal, pres, potential_temperature
real(4) t,s,p

t=min(max(tem,-2.0),50.0)
s=min(max(sal, 0.0),50.0) - 35.0
p=max(pres*1.0e-05,0.0)   !Transition p from Pa to bar

potential_temperature = t + p*( a001 + t*(a101 + t*(a201 + t*a301)) + s*(a011 + t*a111)   &
                          + p*( a002 + t*(a102 + t*a202) + a012*s   + p*(a003 + t*a103) ) )
!A check value is potential_temperature(10°, 25 PSU, 1e+08 Pa) = 8.4678516 °C
endfunction potential_temperature

!==========density for input in situ temperature (Unesco 1981)=======================================+++
function density_unesco_slp(tem,sal)
! tem  - temperature [°C] (in situ)
! sal  - salinity [PSU]
! density_unesco_slp - deviation of density by sea level pressure from 1000.0 [kg/m**3]
! in equations pressure is in bars
implicit none

! polynom coefficients for the fresh water state equation 
real(4), parameter:: fw0=    -0.157406,  fw1= 6.793952e-02, fw2=-9.095290e-03,     &
                     fw3= 1.001685e-04,  fw4=-1.120083e-06, fw5= 6.536332e-09

! polynom coefficients for the surface sea water state equation 
real(4), parameter:: sw01  = 0.824493,    sw11  =-4.0899e-03,   sw21= 7.6438e-05,   &
                     sw31  =-8.2467e-07,  sw41  = 5.3875e-09,                       &
                     sw03p2=-5.72466e-03, sw13p2= 1.0227e-04, sw23p2=-1.6546e-06,   &
                     sw02  = 4.8314e-04

!Coefficients for surface pressure density are same for inputs in temperatures in situ and potential
!The coefficients for Young module are different


real(4) tem, sal, density_unesco_slp
real(4) t, s

t=min(max(tem,-2.0),50.0)
s=min(max(sal, 0.0),50.0)

! density of sea water (p=0) [cm/s**3]
density_unesco_slp = fw0 + t*( fw1  + t*(fw2    + t*(fw3    + t*(fw4    + t*fw5 ) ) ) )   &
                         + s*( sw01 + t*(sw11   + t*(sw21   + t*(sw31   + t*sw41) ) )     &
                           +    sqrt(s)*(sw03p2 + t*(sw13p2 + t *sw23p2) )   +  s*sw02 )

!Values for checking the formula are 
!density_unesco_slp( 5°C, 0 PSU,     0 Pa) = -0.03325 kg/m^3, 
!density_unesco_slp( 5°C,35 PSU,     0 Pa) = 27.67547 kg/m^3,
 
 endfunction density_unesco_slp

!==========density for input in situ temperature (Unesco 1981)=======================================+++
function young_module_slp_tins(tem,sal)
! tem  - temperature [°C] (in situ)
! sal  - salinity [PSU]
! young_module_slp_tins - young module at surface (for in situ temperature input)

implicit none

!Here are ones for in situ temperature (Unesco 1981)

! polynom coefficients for the fresh water Young module
real(4), parameter:: pfw0= 19652.21,     pfw1= 148.4206,    pfw2=-2.327105,       &
                     pfw3= 1.360477e-02, pfw4=-5.155288e-05

! polynom coefficients for the surface sea water Young module
real(4), parameter:: psw01  = 54.6746,   psw11  =-0.603459,    psw21  = 1.09987e-02, psw31=-6.1670e-05,        &
                     psw03p2= 7.944e-02, psw13p2= 1.64830e-02, psw23p2=-5.30090e-04

real(4) tem, sal, young_module_slp_tins
real(4) t, s

t=min(max(tem,-2.0),50.0)
s=min(max(sal, 0.0),50.0)

young_module_slp_tins = pfw0 + t*( pfw1 + t*( pfw2   + t*( pfw3 + t*pfw4 ) ) )      &
                             + s*(psw01 + t*(psw11   + t*(psw21 + t*psw31) )        &
                      +sqrt(s)*(psw03p2 + t*(psw13p2 + t* psw23p2            )  )  )

 endfunction young_module_slp_tins

!==========density for input in potential temperature (Jackett McDougall 1995)===================================+++
function young_module_slp_tpot(theta,sal)
! theta  - temperature [°C] (potential)
! sal  - salinity [PSU]
! young_module_slp_tpot - young module at surface (for in potential temperature input)

implicit none

!Here are ones for in potential temperature (Jackett McDougall 1995)

! polynom coefficients for the fresh water Young module
real(4), parameter:: pfw0= 19659.33,     pfw1= 144.4304,    pfw2=-1.706103,       &
                     pfw3= 9.648704e-03, pfw4=-4.190253e-05

! polynom coefficients for the surface sea water Young module
real(4), parameter:: psw01  = 52.84855,  psw11 =-3.101089e-01, psw21  = 6.283263e-03, psw31=-5.084188e-05,        &
                     psw03p2= 3.88664e-01, psw13p2= 9.085835e-03, psw23p2=-4.619924e-04

real(4) theta, sal, young_module_slp_tpot
real(4) t, s

t=min(max(theta,-2.0),50.0)
s=min(max(sal,0.0),50.0)

young_module_slp_tpot = pfw0 + t*( pfw1 + t*( pfw2   + t*( pfw3 + t*pfw4 ) ) )      &
                             + s*(psw01 + t*(psw11   + t*(psw21 + t*psw31) )        &
                      +sqrt(s)*(psw03p2 + t*(psw13p2 + t* psw23p2            )  )  )

 endfunction young_module_slp_tpot

!==========density for input in situ temperature (Unesco 1981)=======================================+++
function density_unesco_tins(tem,sal,pres,den_surf,k_surf)
! tem  - temperature [°C] (in situ)
! sal  - salinity [PSU]
! pres - pressure[Pa]
! density_unesco_tins - density [kg/m**3]
! in equations pressure is in bars
implicit none

!Here are ones for in situ temperature (Unesco 1981)

!polynom coefficients for the full Young module 
real(4), parameter:: pr001= 3.239908,    pr101= 1.43713e-03, pr201= 1.16092e-04,   pr301=-5.77905e-07,       &
                     pr011= 2.28380e-03, pr111=-1.09810e-05, pr211=-1.60780e-06, pr03p21= 1.910750e-04,      &
                     pr002= 8.50935e-05, pr102=-6.12293e-06, pr202= 5.27870e-08,                             &
                     pr012=-9.93480e-07, pr112= 2.08160e-08, pr212= 9.16970e-10

real(4) tem, sal, pres,den_surf,k_surf,density_unesco_tins
real(4) t, s, p
real(4) k_pres

t=min(max(tem,-2.0),50.0)
s=min(max(sal, 0.0),50.0)
p=max(pres*1.0e-05,0.0)   !Transition p from Pa to bar

k_pres= k_surf + p*( pr001 + t*(pr101 + t*(pr201 + t*pr301) ) + s*(pr011 + t*(pr111 + t* pr211 ) + sqrt(s)*pr03p21 )     &
               + p*( pr002 + t*(pr102 + t* pr202            ) + s*(pr012 + t*(pr112 + t* pr212 )                   )   )   )

density_unesco_tins = sngl( (dble(den_surf) + 1000.0d0) /(1.0d0-dble(p)/dble(k_pres)) - 1000.0d0 )

!Values for checking the formula are 
!density_unesco_tins( 5°C, 0 PSU,     0 Pa) = -0.03325 kg/m^3, 
!density_unesco_tins( 5°C,35 PSU,     0 Pa) = 27.67547 kg/m^3,
!density_unesco_tins(25°C,35 PSU, 1e+08 Pa) = 62.5381  kg/m^3
 
 endfunction density_unesco_tins

 !==========density for input inpotential temperature (Jackett McDougall 1995)=======================================+++
function density_unesco_tpot(theta,sal,pres,den_surf,k_surf)
! tem  - temperature [°C] (potential)
! sal  - salinity [PSU]
! pres - pressure[Pa]
! density_unesco_tins - density [kg/m**3]
! in equations pressure is in bars
implicit none

!Here are ones for in potential temperature (Jackett McDougall 1995)

!polynom coefficients for the full Young module 
real(4), parameter:: pr001= 3.186519,     pr101= 2.212276e-02, pr201=-2.984642e-04, pr301  = 1.956415e-06,    &
                     pr011= 6.704388e-03, pr111=-1.847318e-04, pr211= 2.059331e-07, pr03p21= 1.480266e-04,    &
                     pr002= 2.102898e-04, pr102=-1.202016e-05, pr202= 1.394680e-07,                           &
                     pr012=-2.040237e-06, pr112= 6.128773e-08, pr212= 6.207323e-10

real(4) theta, sal, pres,den_surf,k_surf,density_unesco_tpot
real(4) t, s, p
real(4) k_pres

t=min(max(theta,-2.0),50.0)
s=min(max(sal,0.0),50.0)
p=max(pres*1.0e-05,0.0)   !Transition p from Pa to bar

k_pres= k_surf + p*( pr001 + t*(pr101 + t*(pr201 + t*pr301) ) + s*(pr011 + t*(pr111 + t* pr211 ) + sqrt(s)*pr03p21 )     &
               + p*( pr002 + t*(pr102 + t* pr202            ) + s*(pr012 + t*(pr112 + t* pr212 )                   )   )   )

density_unesco_tpot = sngl( (dble(den_surf) + 1000.0d0) /(1.0d0-dble(p)/dble(k_pres)) - 1000.0d0 )

!Values for checking the formula are 
!density_unesco_tins( 5°C, 0 PSU,     0 Pa) = -0.03325 kg/m^3, 
!density_unesco_tins( 5°C,35 PSU,     0 Pa) = 27.67547 kg/m^3,
!density_unesco_tins(25°C,35 PSU, 1e+08 Pa) = 62.5381  kg/m^3
 
 endfunction density_unesco_tpot

!=======================================================================
subroutine steric_level(den0,den,sls)
use basin_grid

real(4),intent(in):: den0(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &
                      den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
real(4),intent(inout):: sls(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

integer m,n,k

real(4) denfrac

!$omp parallel do private(m,n,k,denfrac)
  do n=ny_start,ny_end
    do m=nx_start,nx_end
     if(lu(m,n)>0.5) then      
       sls(m,n)=0.0
       do k=1,nz
        denfrac = (den0(m,n,k)+1000.0)/(den(m,n,k)+1000.0)
        sls(m,n)=sls(m,n)+log(denfrac)*dz(k)*hhq(m,n)
       enddo
     endif    
    enddo
  enddo
!$omp end parallel do 

endsubroutine steric_level

endmodule ocalg_routes