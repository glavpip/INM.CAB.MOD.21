module ocpar_routes
use constants
use basin_grid
use ocean_variables
use key_switches
use iodata_routes
use rwpar_routes
use gridcon_routes
use basinpar_routes
use ocalg_routes
!use sea_level_no_split
use parallel_sea_level_no_split
use depth_routes
use velocity_routes

implicit none

contains

subroutine ocean_model_parameters(tau)
use basin_grid
use ocean_variables
use key_switches
use iodata_routes
use rwpar_routes
use gridcon_routes
use basinpar_routes
use ocalg_routes

real(4), intent(in):: tau

character(256)    t_mask_file,       &  !name of file with temperature point sea-land mask
       bottom_topography_file,       &  !name of file with bottom topography
       help_string
integer  m, n, k, ierr
real(4) difrate

! define parameters of task
! description of parameters see in file with mame filepar

if (rank == 0) then
 open (90, file='phys_proc.par', status='old')

  read(90,*)  ksw_ts,  ksw_uv                       !TS (left) and UV (right) equation solving (0-no, 1-yes)                                                                                              
  read(90,*)  ksw_age, ksw_pt                       !Ideal age (left) and PassTracer (right) equation solving (0-no, 1-yes)                                                                                                         
  read(90,*)  ksw_lat, ksw_lat4                     !Lateral 2nd (left, 0-no, 1-sigma, 2-z, 3-rho (2&3 for TS)) and 4th (right, 0-no, 1-yes, UV only) order mix
  read(90,*)  ksw_vert                              !vertical mix (0-constant coeff, 1-PP, 2-MY)                                                                                            
  read(90,*)  ksw_dens                              !pressure gradient computation (0-no (constant density), 1-local potential density)
  read(90,*)  ksw_ith, ksw_itr, ksw_idn             !sea ice thermodynamics, transport and dynamics using (0-no, 1-yes)                                                                                                                      
  read(90,*)  ksw_ssbc                              !Type of surface boundary conditions (1-surface T&S and wind stress are prescribed; 2-T&S fluxes and wind stress are prescribed; 3-T&S fluxes and wind stress are simulated
  read(90,*)  ksw_wflux                             !normalize global mean salt balance (0-no, 1-normalize water/salt flux (depending on variable volume))                                                               
  read(90,*)  ksw_lbc_ts, ksw_lbc_uv, ksw_lbc_ssh   !open boundary conditions for T&S, U&V and SSH using (0-no, 1-yes)                                                                                                         

  read(90,*) sst_relax,    sss_relax            !Surface nudging rate for temperature and salinity [m/s]            
  read(90,*) ldiff_ts_ref, ldiff_ts_upw, ldiff_ts_smg         !lateral diffusivity ratio (undim) for temperature: background, upwind and smagorinsky                               
  read(90,*) lvisc_ref,       lvisc_upw, lvisc_smg            !lateral vicosity(2nd order) ratio (undim): background, upwind and smagorinsky
  read(90,*) lvisc_4_ref,   lvisc_4_upw, lvisc_4_smg          !lateral vicosity(4th order) ratio (undim): background, upwind and smagorinsky
  read(90,*) tsfrac_lat                         !fraction of salinity lateral diffusivity due to one for temperature
  read(90,*) vdiff_ts_min, vdiff_ts_max         !vertical background and top diffusivity for T [m**2/s]
  read(90,*) vvisc_min,    vvisc_max            !vertical background and top viscosity for UV [m**2/s]        
  read(90,*) tsfrac_vert                        !fraction of salinity lateral diffusion due to one for temperature
  read(90,*) rfric_min, rfric_max               !minimum and maximum rayleigh dissipation rate (1/s)
  
  help_string =' '
  read (90,'(a)') help_string   ! file with t-mask'
  call get_first_lexeme(help_string ,t_mask_file   )

  help_string =' '      
  read (90,'(a)') help_string  ! file with bottom topography'
  call get_first_lexeme(help_string ,bottom_topography_file  )

 close(90)
endif

call mpi_bcast(ksw_ts      , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_uv      , 1, mpi_integer, 0, cart_comm, ierr) 
call mpi_bcast(ksw_age     , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_pt      , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_lat     , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_lat4    , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_vert    , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_dens    , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_ith     , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_itr     , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_idn     , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_ssbc    , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_wflux   , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_lbc_ts  , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_lbc_uv  , 1, mpi_integer, 0, cart_comm, ierr)   
call mpi_bcast(ksw_lbc_ssh , 1, mpi_integer, 0, cart_comm, ierr)   

call mpi_bcast(sst_relax   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(sss_relax   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(ldiff_ts_ref, 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(ldiff_ts_upw, 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(ldiff_ts_smg, 1, mpi_real4,   0, cart_comm, ierr) 
call mpi_bcast(lvisc_ref   , 1, mpi_real4,   0, cart_comm, ierr)
call mpi_bcast(lvisc_upw   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(lvisc_smg   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(lvisc_4_ref , 1, mpi_real4,   0, cart_comm, ierr)
call mpi_bcast(lvisc_4_upw , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(lvisc_4_smg , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(tsfrac_lat  , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(vdiff_ts_min, 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(vdiff_ts_max, 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(vvisc_min   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(vvisc_max   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(tsfrac_vert , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(rfric_min   , 1, mpi_real4,   0, cart_comm, ierr)   
call mpi_bcast(rfric_max   , 1, mpi_real4,   0, cart_comm, ierr) 

call mpi_bcast(t_mask_file, 256, mpi_character, 0, cart_comm, ierr)            ! file with t-mask
call mpi_bcast(bottom_topography_file, 256, mpi_character, 0, cart_comm, ierr) ! file with bottom topography

if (rank == 0) then
  write(*,'(2i4,a)') ksw_ts,  ksw_uv,          ' - TS (left) and UV (right) equation'
  write(*,'(2i4,a)') ksw_age, ksw_pt,          ' - Ideal age (left) and PassTracer (right) equation'
  write(*,'(2i4,a)') ksw_lat, ksw_lat4,        ' - Lateral 2nd and 4th  order mix' 
  write(*,'(i9,a)') ksw_vert,                  ' - Vertical mix '
  write(*,'(i9,a)') ksw_dens,                  ' - Pressure gradient '
  write(*,'(3i3,a)') ksw_ith,ksw_itr,ksw_idn,  ' - Sea ice thdyn, trans and dyn '
  write(*,'(i9,a)') ksw_ssbc,                  ' - Type of surface boundary conditions'
  write(*,'(i9,a)') ksw_wflux,                 ' - Normalize global mean water/salt balance'
  write(*,'(3i3,a)') ksw_lbc_ts,ksw_lbc_uv,ksw_lbc_ssh,   ' - Open boundary for T&S, U&V and SSH '
  write(*,'(2es12.4,a)') sst_relax, sss_relax,        ' - Nudging for temperature and salinity[m/s]'
  write(*,'(3es12.4,a)') ldiff_ts_ref, ldiff_ts_upw, ldiff_ts_smg,  ' - Lateral diffusivity for temperature '
  write(*,'(3es12.4,a)') lvisc_ref,    lvisc_upw,    lvisc_smg,     ' - Lateral vicosity(2nd order)'
  write(*,'(3es12.4,a)') lvisc_4_ref,  lvisc_4_upw,  lvisc_4_smg,   ' - Lateral vicosity(4th order)'
  write(*,'(es12.4,a)') tsfrac_lat,    ' - fraction of salinity lateral diffusion'
  write(*,'(2es12.4,a)') vdiff_ts_min, vdiff_ts_max,  ' - vertical background and top diffusivity for T [m**2/s]'
  write(*,'(2es12.4,a)') vvisc_min, vvisc_max,        ' - vertical background and top viscosity for UV [m**2/s]'
  write(*,'(es12.4,a)') tsfrac_vert,   ' - fraction of salinity vertical diffusion'
  write(*,'(2es12.4,a)') rfric_min, rfric_max,        ' - minimum and maximum Rayleigh dissipation rate (1/s)'
  write(*,'(a,a)')  ' file with T-point sea-land mask: ', t_mask_file(1:len_trim (t_mask_file))
  write(*,'(a,a)')  '     file with bottom topography: ', bottom_topography_file(1:len_trim (bottom_topography_file))
endif

igrzts = min(IABS(ksw_ssbc),2) ! type of condition for T and S on sea surface
igrzpt = 2
igrzage= 1

dkft = sst_relax
dkfs = sss_relax

! area mask initialization
call gridcon(t_mask_file)
! setting vertical t-,w- grid levels
call vgrid
! define grid geographical coordinates, steps and coriolis parameters
call basinpar

call prdstd(' ', bottom_topography_file, 1, hhq_rest, lu1, nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)

call syncborder_real(hhq_rest, 1)

if(periodicity_x/=0) then
  call cyclize_x(hhq_rest,nx,ny,1,mmm,mm)
end if

if(periodicity_y/=0) then
  call cyclize_y(hhq_rest,nx,ny,1,nnn,nn)
end if

!-------------define diffusive and viskosity coefficients-------------
!$omp parallel do private(m,n,k, difrate)
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1
         difrate=sqt(m,n)/8.0
         do k=1,nz    
!         for temperature at T-grid
         amts_ref(m,n,k) = ldiff_ts_ref*difrate/tau
!         for velocity at T-grid
         amuv_ref(m,n,k) = lvisc_ref*difrate/tau
         amuv4_ref(m,n,k) = sqrt(lvisc_4_ref/tau)*difrate
          if(m==nx/2.and.n==ny/2.and.k==1) then
           write(*,*) 'amts_ref(nx/2,ny/2,1)=',amts_ref(m,n,k)
           write(*,*) 'amuv_ref(nx/2,ny/2,1)=',amuv_ref(m,n,k)
           write(*,*) 'amuv4_ref(nx/2,ny/2,1)=',amuv4_ref(m,n,k)**2
          endif                  
         end do
       enddo
      enddo
!$omp end parallel do

! solar shortwave heating penetrate divergence function.
call shortwave_divergence

!$omp parallel do private(m,n)
    do n=ny_start-1,ny_end+1
     do m=nx_start-1,nx_end+1
! set upper coefficients for temperature vertical diffusion
            anzt(m,n,1) = vdiff_ts_max*lu(m,n)
! set upper coefficients for momentum vertical viscosity
            anzu(m,n,1) = vvisc_max*lu(m,n)       
! set coefficients for temperature vertical diffusion
            anzt(m,n,2:nz+1) = vdiff_ts_min*lu(m,n)
!c set coefficients for velocity vertical viscosity
            anzu(m,n,2:nz+1) = vvisc_min*lu(m,n)               
     enddo
  enddo
!$omp end parallel do

endsubroutine ocean_model_parameters

!====================================================================
subroutine shortwave_divergence

integer m,n,k, ierr
!-----------------------------------------------------------------------
!     Solar Shortwave energy penetrates below the ocean surface. Clear
!     water assumes energy partitions between two exponentials as
!     follows:
!
!     58% of the energy decays with a 35 cm e-folding scale
!     42% of the energy decays with a 23 m  e-folding scale
!
!     Paulson and Simpson (1977 Irradiance measurements in the upper
!                               ocean JPO 7, 952-956)
!     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
!                  A General Circulation Model for Upper Ocean
!                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
!-----------------------------------------------------------------------
!
! => Shortwave penetration is a double exponential as follows:
!     parameter(rpart = 0.58,efold1 = 0.35,efold2 = 23.0)
! new approuch of using short wave radiation: upper part
! of about 60% added to heat flux and residual part
! of 40% of the energy decays with a 20 m e-folding scale
      
      real(8) buf_real8
      real(8) swarg,pen1,pen2,spacesum,pointsum, sumk

      spacesum=0.0
      pointsum=0.0

      do n = ny_start,ny_end
        do m = nx_start,nx_end

        if(lu(m,n)>0.5) then

          pointsum=pointsum+1.0 
          pen2=1.0
          sumk = 0.0

          do k=1,nz
           swarg = -zw(k+1)*hhq_rest(m,n)/efold
           pen1 = pen2
           pen2 = dexp(swarg)           
           divswrad(m,n,k) = pen1 - pen2
           sumk = sumk + divswrad(m,n,k)       
          enddo

          do k=1,nz      
           divswrad(m,n,k) = divswrad(m,n,k)/sumk
           spacesum = spacesum + divswrad(m,n,k)   
          enddo
                     
        end if
       enddo
      enddo

    buf_real8 = pointsum
    call mpi_allreduce(buf_real8, pointsum, 1, mpi_real8, mpi_sum, cart_comm, ierr)
    buf_real8 = spacesum
    call mpi_allreduce(buf_real8, spacesum, 1, mpi_real8, mpi_sum, cart_comm, ierr)

      if (rank == 0) then
        write(*,'(2x,a,f10.2)') 'sum of SW divergence coefficient =',spacesum
        write(*,'(2x,a,f10.2,a)')'for ',pointsum,' t-grid points'
      endif

endsubroutine shortwave_divergence

!=======================================================
subroutine inicond(start_type,path2ocp,tau)

    integer, intent(in):: start_type    
    real(4), intent(in):: tau
    character*(*), intent(in):: path2ocp

    character(256) fname
    integer m, n, k, ierr, lqp, grad
    real(4) hx2, hy2, pgrz
    real*8 :: time_count

    if(rank==0) write(*,*) 'Reading checkpoints for initial conditions'

    ! initial conditions for temperature, salinity

    !read potential temperature    
    call prdstd(path2ocp,'cptt.dat', 1,tt,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
    !read salinity       
    call prdstd(path2ocp,'cpss.dat', 1,ss,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
    if(rank==0.and.deb_out_scr>0) write(*,*) 'T and S are read'

    call syncborder_real(ss, nz)
    call syncborder_real(tt, nz)

    if(periodicity_x/=0) then
        call cyclize_x(tt,nx,ny,nz,mmm,mm)
        call cyclize_x(ss,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize_y(tt,nx,ny,nz,nnn,nn)
        call cyclize_y(ss,nx,ny,nz,nnn,nn)
    end if

    if (start_type>0) then

        if(iabs(ksw_age)>0) then
        
        !read ideal age          
          call prdstd(path2ocp,'cpage.dat', 1, age,lu ,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
          if(rank==0.and.deb_out_scr>0) write(*,*) 'Ideal age is read'
          call syncborder_real(age, nz)
          
          if(periodicity_x/=0) then
              call cyclize_x(age,nx,ny,nz,mmm,mm)
          endif
          
          if(periodicity_y/=0) then
              call cyclize_y(age,nx,ny,nz,nnn,nn)
          endif
          
        endif

        !read zonal velocity
        call prdstd(path2ocp,'cpuu.dat', 1 ,uu, llu, nx,ny,nz, mmm-1,mm,nnn  ,nn,1,nz,ierr)
        !read meridional velocity
        call prdstd(path2ocp,'cpvv.dat', 1 ,vv, llv, nx,ny,nz, mmm  ,mm,nnn-1,nn,1,nz,ierr)
        if(rank==0.and.deb_out_scr>0) write(*,*) 'U and V are read'

        call syncborder_real(uu, nz)
        call syncborder_real(vv, nz)

        if(periodicity_x/=0) then
            call cyclize_x(uu ,nx,ny,nz  ,mmm,mm)
            call cyclize_x(vv ,nx,ny,nz  ,mmm,mm)
        end if

        if(periodicity_y/=0) then
            call cyclize_y(uu ,nx,ny,nz  ,nnn,nn)
            call cyclize_y(vv ,nx,ny,nz  ,nnn,nn)
        end if

        !  read turbulent kinetic energy and length scale       
        call prdstd(path2ocp,'cpq2.dat' , 1, q2  , lu, nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)
        call prdstd(path2ocp,'cpq2l.dat', 1, q2l , lu, nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)
        if(rank==0.and.deb_out_scr>0) write(*,*) 'K and KL are read'
        call syncborder_real(q2, nz+1)
        call syncborder_real(q2l, nz+1)

        if(periodicity_x/=0) then
            call cyclize_x(q2  ,nx,ny,nz+1,mmm,mm)
            call cyclize_x(q2l ,nx,ny,nz+1,mmm,mm)
        end if

        if(periodicity_y/=0) then
            call cyclize_y(q2  ,nx,ny,nz+1,nnn,nn)
            call cyclize_y(q2l ,nx,ny,nz+1,nnn,nn)
        end if       

        ! read sea surface heights for present(1) and previous(2) time steps
        call prdstd8(path2ocp,'cpssh8.dat', 1,ssh ,lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
        if(rank==0.and.deb_out_scr>0) write(*,*) 'SSH is read'
        call syncborder_real8(ssh, 1)

        if(periodicity_x/=0) then
            call cyclize8_x(ssh,nx,ny,1,mmm,mm)
        end if

        if(periodicity_y/=0) then
            call cyclize8_y(ssh , nx,ny,1,nnn,nn)
        end if

        ! reading hice - ice volume per square      
        call prdstd(path2ocp,'cphice.dat',1,hice,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
        ! reading aice - ice compactness
        call prdstd(path2ocp,'cpaice.dat',1,aice,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
        ! reading hsnow - snow volume per square
        call prdstd(path2ocp,'cphsnow.dat',1,hsnow,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
        ! reading sice - salt per square
        call prdstd(path2ocp,'cpsice.dat',1,sice,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
        if(rank==0.and.deb_out_scr>0) write(*,*) 'Aice, Hice, Hsnow and Sice are read'

        ! read zonal ice velocity
        call prdstd(path2ocp,'cpuice.dat', 1,uice, llu,nx,ny,1, mmm-1,mm,nnn  ,nn,1,1,ierr)
        ! read meridional ice velocity
        call prdstd(path2ocp,'cpvice.dat', 1,vice, llv,nx,ny,1, mmm  ,mm,nnn-1,nn,1,1,ierr)
        if(rank==0.and.deb_out_scr>0) write(*,*) 'Uice and Vice are read'

        call syncborder_real( aice , mgrad)
        call syncborder_real(  hice, mgrad)
        call syncborder_real( hsnow, mgrad)        
        call syncborder_real(  sice, mgrad)       
        call syncborder_real(  uice, 1)
        call syncborder_real(  vice, 1)
        
        if(periodicity_x/=0) then
            call cyclize_x( aice ,nx,ny,mgrad,mmm,mm)
            call cyclize_x(  hice,nx,ny,mgrad,mmm,mm)
            call cyclize_x( hsnow,nx,ny,mgrad,mmm,mm)        
            call cyclize_x(  sice,nx,ny,mgrad,mmm,mm)
            call cyclize_x(  uice,nx,ny,1,mmm,mm)
            call cyclize_x(  vice,nx,ny,1,mmm,mm)
        end if
       
        if(periodicity_y/=0) then
            call cyclize_y( aice ,nx,ny,mgrad,nnn,nn)
            call cyclize_y(  hice,nx,ny,mgrad,nnn,nn)
            call cyclize_y( hsnow,nx,ny,mgrad,nnn,nn)        
            call cyclize_y(  sice,nx,ny,mgrad,nnn,nn)
            call cyclize_y(  uice,nx,ny,1,nnn,nn)
            call cyclize_y(  vice,nx,ny,1,nnn,nn)
        end if
       
        !$omp parallel do private(m,n,grad) 
        do n = ny_start-1,ny_end+1
          do m = nx_start-1,nx_end+1
            if (lu(m,n)>0.5) then
             aistot(m,n) = 0.0
              hitot(m,n) = 0.0        
              hstot(m,n) = 0.0 
             do grad=1,mgrad
               aistot(m,n) = aistot(m,n) + aice(m,n,grad)
                hitot(m,n) =  hitot(m,n) + hice(m,n,grad)
                hstot(m,n) =  hstot(m,n) + hsnow(m,n,grad)
             enddo
             mistot(m,n) = den_ice*hitot(m,n)+den_snow*hstot(m,n)
            endif
          enddo
        enddo
        !$omp end parallel do    
    
    endif

    if(iabs(ksw_pt)>0) then
      
      call prdstd(path2ocp,'cppt.dat', 1, pass_tracer,lu ,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
      if(rank==0.and.deb_out_scr>0) write(*,*) 'Passive tracer is read'
      call syncborder_real(pass_tracer, nz)
      
      if(periodicity_x/=0) then
          call cyclize_x(pass_tracer,nx,ny,nz,mmm,mm)
      endif
      
      if(periodicity_y/=0) then
          call cyclize_y(pass_tracer,nx,ny,nz,nnn,nn)
      endif
      
    endif

    !initialize depth for internal mode      
    call hh_init(hhq, hhu, hhv, hhh,   &
                 ssh, hhq_rest, hhu_rest, hhv_rest)
      
    !$omp parallel do private(m,n,k) 
    do n = ny_start-1, ny_end+1
      do m = nx_start-1, nx_end+1
        do k = 1, nz
          uh(m,n,k)=uu(m,n,k)*hhu(m,n)
          vh(m,n,k)=vv(m,n,k)*hhv(m,n)
         aux_array3d_tgr1(m,n,k)=uh(m,n,k)
         aux_array3d_tgr2(m,n,k)=vh(m,n,k)        
        enddo
        hhq_n(m,n) = hhq(m,n)
        hhu_n(m,n) = hhu(m,n)
        hhv_n(m,n) = hhv(m,n)
        hhh_n(m,n) = hhh(m,n)
      enddo
    enddo
    !$omp end parallel do

!-----------------density definition-----------------------------------
    if (iabs(ksw_dens)>0) then
        !$omp parallel do private(m,n,k,pgrz) 
        do n=ny_start-1,ny_end+1
          do m=nx_start-1,nx_end+1
            if (lu(m,n)>0.5) then
              pgrz=FreeFallAcc*RefDen*hhq(m,n) 
              do k=1,nz
                den_sg0(m,n,k) = density_unesco_slp(tt(m,n,k),ss(m,n,k)+salref)                
                yng_slp(m,n,k) = young_module_slp_tpot(tt(m,n,k),ss(m,n,k)+salref)                
                den_sgt(m,n,k) = density_unesco_tpot(tt(m,n,k),ss(m,n,k)+salref,pgrz*z(k),      &
                                                den_sg0(m,n,k),yng_slp(m,n,k))
                den_sgm(m,n,k) = density_unesco_tpot(tt(m,n,k),ss(m,n,k)+salref,hmp2*1.0e+04,   &
                                                den_sg0(m,n,k),yng_slp(m,n,k))
              enddo
            endif
          enddo
        enddo
        !$omp end parallel do
    endif

    if(start_type==0) then
     den_ini = den_sg0
    else
     !read potential temperature    
     call prdstd(path2ocp,'cpden0.dat', 1,den_ini,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
    endif

    call steric_level(den_ini,den_sg0,sls)

    if(iabs(ksw_vert)>0) then
!Computing Richardson number terms
     call vertical_frequencies(uu,vv,tt,ss,den_sg0,yng_slp,vbf2,shf2)
!Mixed layer depth based on slp-referred potential density
     call mld_density(den_sg0,mld_dens)
    endif
    
    !storing old velocities

      !Removing barotropic components from the 3D velocity
      call depth_ave(aux_array3d_tgr1,ubrtr,llu,1) 
      call depth_ave(aux_array3d_tgr2,vbrtr,llv,1)

!Preliminary vertical velocity for momentum transport
      call vertical_velocity(aux_array3d_tgr1,aux_array3d_tgr2,ww,vol_flux)

    !--------------Rayleigh friction initialization
    !$omp parallel do private(m,n,k, hx2, hy2) 
    do n = ny_start,ny_end
      do m = nx_start,nx_end

        hx2=  ((hhq_rest(m+1,n)-hhq_rest(m  ,n))/dxt(m  ,n))**2     &
             +((hhq_rest(m  ,n)-hhq_rest(m-1,n))/dxt(m-1,n))**2 
        hy2=  ((hhq_rest(m,n+1)-hhq_rest(m,n  ))/dyt(m,n  ))**2     &
             +((hhq_rest(m,n  )-hhq_rest(m,n-1))/dyt(m,n-1))**2  
        r_diss(m,n)=rfric_min+(rfric_max-rfric_min)*sqrt((hx2+hy2)/2.0)
        AlbOpWater(m,n) = 0.069 - 0.011*cosd(2.0*sngl(geo_lat_t(m,n)))          
      
      enddo
    enddo
    !$omp end parallel do
  
    call syncborder_real(r_diss, 1)
    call syncborder_real(AlbOpWater, 1)

    if(periodicity_x/=0) then
        call cyclize_x( r_diss, nx,ny,1,mmm,mm)
        call cyclize_x( AlbOpWater, nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize_y( r_diss, nx,ny,1,nnn,nn)
        call cyclize_y( AlbOpWater, nx,ny,1,nnn,nn)
    end if

    ! Sea level initialization
    !CALL FORM_MATRICE_SSH(r_diss, tau)
    !CALL FORM_PRECONDITION() 
    call init_parallel_sea_level_solver()
    if (rank .eq. 0) then
        print *, "Begin Matrix Form!"
    endif
    call start_timer(time_count)
    call form_matrice_ssh(tau, r_diss)
    call end_timer(time_count)
    if (rank .eq. 0) then
        print *, "Matrix Form time is ", time_count
    endif
    call init_ksp_solver()

    if(rank==0) write(*,*) 'Checkpoints are read, initialization is done'

endsubroutine inicond

endmodule ocpar_routes
