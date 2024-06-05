module ocstep_routes
use constants
use basin_grid
use ocean_variables
use ocean_bc
use key_switches
use ocalg_routes
use depth_routes
use mixing_routes
use flux_routes
use ts_routes
use velocity_routes
!use sea_level_no_split
use parallel_sea_level_no_split
implicit none

contains

subroutine ocean_model_step(tau, nstep_icedyn)

 real(4),    intent(in):: tau
 integer,    intent(in):: nstep_icedyn
 
 integer m, n, k, mark

 real(4) pgrz

! Computing 3D depth array  

!$omp parallel do private(m, n, k) 
   do n=ny_start-1,ny_end+1
    do m=nx_start-1,nx_end+1    
      z3d(m,n,1:nz) = hhq(m,n) * z(1:nz)*lu(m,n)
    enddo
   enddo
!$omp end parallel do  

! Computing depth-monotonic potential (relatively to the basin mean depth) density array for lateral diffusion  
  if(iabs(ksw_lat)>=3) then
!$omp parallel do private(m, n, k) 
   do n=ny_start-1,ny_end+1
    do m=nx_start-1,nx_end+1    
       den_isopyc(m,n,1) = den_sgm(m,n,1)*lu(m,n)
      do k=2,nz
       den_isopyc(m,n,k) = max(den_sgm(m,n,k),den_isopyc(m,n,k-1)+epsrho)*lu(m,n)
      enddo
    enddo
   enddo
!$omp end parallel do    
  endif

! boundary condition computation
 call bndcond(tau,nstep_icedyn)

! mixing parameters computation
 call mixing_processes(tau)
 
! Hydrodynamics module (velocity + SSH)
 if (iabs(ksw_uv)>0) then
  call hydrodynamics(tau)
 endif 

!Turbulent variables transport-diffusion   
 if(iabs(ksw_vert)>=2) then
  call turbulence_tran_diff(tau)                                         
 endif

 if(iabs(ksw_ts)>0) then
!Temperature and salinity transport-diffusion    
   call tem_sal_tran_diff(tau)         
 endif
    
 if(iabs(ksw_age)>0) then
!ideal age transport-diffusion 
   call ideal_age_tran_diff(tau)         
 endif

!Passive tracer transport-diffusion
 if(iabs(ksw_pt)>0) then
   call passtracer_tran_diff(tau)
 endif

! Shifting HHQ on time
 if (nonlinear_free_surface>0) then
 !$omp parallel do private(m,n,k)
   do n=ny_start-1, ny_end+1
    do m=nx_start-1, nx_end+1
     hhq(m,n)=hhq_n(m,n)
     hhu(m,n)=hhu_n(m,n)
     hhv(m,n)=hhv_n(m,n)
     hhh(m,n)=hhh_n(m,n)
    enddo
   enddo
 !$omp end parallel do
 endif

!-----------------density definition-----------------------------------
 if (ksw_dens>0) then
   !$omp parallel do private(m,n,k,pgrz) 
   do n=ny_start-1,ny_end+1
     do m=nx_start-1,nx_end+1
       if(lu(m,n)>0.5) then
        pgrz=FreeFallAcc*RefDen*hhq(m,n) 
        do k=1,nz
         den_sg0(m,n,k) = density_unesco_slp(tt(m,n,k),ss(m,n,k)+salref)                
         yng_slp(m,n,k) = young_module_slp_tpot(tt(m,n,k),ss(m,n,k)+salref)                
         den_sgt(m,n,k) = density_unesco_tpot(tt(m,n,k),ss(m,n,k)+salref,pgrz*z(k),      &
                                         den_sg0(m,n,k),yng_slp(m,n,k))
         den_sgm(m,n,k) = density_unesco_tpot(tt(m,n,k),ss(m,n,k)+salref,hmp2*1.0e+04,        &
                                         den_sg0(m,n,k),yng_slp(m,n,k))
        enddo
       endif
     enddo
   enddo
   !$omp end parallel do
 
 endif

     if(iabs(ksw_vert)>0) then
!Computing Richardson number terms
     call vertical_frequencies(uu,vv,tt,ss,den_sg0,yng_slp,vbf2,shf2)
!Mixed layer depth based on slp-referred potential density
     call mld_density(den_sg0,mld_dens)
    endif

    call steric_level(den_ini,den_sg0,sls)

 call check_errors()

endsubroutine ocean_model_step

!=======================================================================
subroutine check_errors()
  
  integer m,n,k,mark, ierr, grad

  mark=0

  !$omp parallel do private(m,n,k,grad)
      do n=bnd_y1,bnd_y2
       do m=bnd_x1,bnd_x2

      if(lu(m,n)>0.5) then
          if(ssh(m,n)<20.0.and.ssh(m,n)>-20.0) then
            continue
          else
            write(*,*) rank, 'in the point m=',m,'n=',n,'ssh=',ssh(m,n)
            mark=1
          endif

          do grad=1,mgrad
           
           if(aice(m,n,grad)<=1.0.and.aice(m,n,grad)>=0.0) then
             continue
           else
             write(*,*) rank, 'in the point m=',m,'n=',n,'grad=',grad,'aice=',aice(m,n,grad)
             mark=1
           endif

           if(hice(m,n,grad)<=100.0.and.hice(m,n,grad)>=0.0) then
             continue
           else
             write(*,*) rank, 'in the point m=',m,'n=',n,'grad=',grad,'hice=',hice(m,n,grad)
             mark=1
           endif

           if(hsnow(m,n,grad)<=100.0.and.hsnow(m,n,grad)>=0.0) then
             continue
           else
             write(*,*) rank, 'in the point m=',m,'n=',n,'grad=',grad,'hsnow=',hsnow(m,n,grad)
             mark=1
           endif

           if(sice(m,n,grad)<=10000.0.and.sice(m,n,grad)>=0.0) then
             continue
           else
             write(*,*) rank, 'in the point m=',m,'n=',n,'grad=',grad,'sice=',sice(m,n,grad)
             mark=1
           endif
          
          enddo

          do k=1,nz
            if(tt(m,n,k)<100.0.and.tt(m,n,k)>-100.0) then
              continue
            else
              write(*,*) rank, 'in the point m=',m,'n=',n,'k=',k,'tt=',tt(m,n,k)
              mark=1
            endif

            if(ss(m,n,k)<100.0.and.ss(m,n,k)>-100.0) then
              continue
            else
              write(*,*) rank, 'in the point m=',m,'n=',n,'k=',k,'ss=',ss(m,n,k)
              mark=1
            endif
          enddo         
      endif

      if(llu(m,n)>0.5) then
        do k=1,nz
          if(uu(m,n,k)<100.0.and.uu(m,n,k)>-100.0) then
            continue
          else
            write(*,*) rank, 'in the point m=',m,'n=',n,'k=',k,'uu=',uu(m,n,k)
            mark=1
          endif

        enddo 
        
          if(uice(m,n)<=100.0.and.uice(m,n)>=-100.0) then
            continue
          else
            write(*,*) rank, 'in the point m=',m,'n=',n,'uice=',uice(m,n)
            mark=1
          endif
                   
      endif

      if(llv(m,n)>0.5) then
        do k=1,nz
          if(vv(m,n,k)<100.0.and.vv(m,n,k)>-100.0) then
            continue
          else
            write(*,*) rank, 'in the point m=',m,'n=',n,'k=',k,'vv=',vv(m,n,k)
            mark=1
          endif

        enddo
          
          if(vice(m,n)<=100.0.and.vice(m,n)>=-100.0) then
            continue
          else
            write(*,*) rank, 'in the point m=',m,'n=',n,'vice=',vice(m,n)
            mark=1
          endif                 
      endif

    enddo
  enddo
  !$omp end parallel do

  if (mark>0) then
    call mpi_abort(cart_comm, 1, ierr)
    stop
  endif

endsubroutine check_errors

!=======================================================================
subroutine bndcond(tau,nstep_icedyn)
    
    real(4), intent(in):: tau
    integer, intent(in):: nstep_icedyn
    integer m,n,k
    integer ierr, lqp


!!! Correcting open boundary data
    if(ksw_lbc_ts>0) then
     do lqp=1, numb_of_lqp
        m=lqpx(lqp)        
        n=lqpy(lqp)
        if(m>=bnd_x1.and.m<=bnd_x2.and.n>=bnd_y1.and.n<=bnd_y2) then
         do k=1,nz
          tt(m,n,k)=tt(m,n,k)+tlbc_relax*tau*(tlqbw(lqp,k)-tt(m,n,k))
          ss(m,n,k)=ss(m,n,k)+slbc_relax*tau*(slqbw(lqp,k)-ss(m,n,k))
         enddo
        endif
     enddo
    endif

    if(ksw_lbc_ssh>0) then
     
     do lqp=1, numb_of_lqp
        m=lqpx(lqp)        
        n=lqpy(lqp)
        if(m>=bnd_x1.and.m<=bnd_x2.and.n>=bnd_y1.and.n<=bnd_y2) then
          ssh(m,n)=sshlqbw(lqp)
        endif
     enddo
     !Updating depth functions  
     if(nonlinear_free_surface>0) then
       call hh_update(hhq, hhu, hhv, hhh, ssh, hhq_rest)
     endif
    
    endif
        
    ! Computation of sea surface boundary conditions and sea ice
    if(ksw_ssbc>0) then
     call sea_surface_fluxes(tau,nstep_icedyn)
    endif

    ! Computing bottom stresses
     call sea_bottom_fluxes

endsubroutine bndcond

!=======================================================================
subroutine mixing_processes(tau)

    real(4), intent(in):: tau

    if(iabs(ksw_lat)>0.or.iabs(ksw_lat4)>0) then
     ! Computing hyfrid background+smagorinsky dissipation
      call lateral_mix_coeff(ldiff_ts_upw, lvisc_upw, lvisc_4_upw,      &
                             ldiff_ts_smg, lvisc_smg, lvisc_4_smg,      &
                     stress_t,stress_s,amts_ref,amuv_ref,amuv4_ref,     &
                     amts,amuv,amuv4,uu,vv)
    end if

    ! vertical diffusion coefficients
    if(iabs(ksw_vert)==1) then
      ! Philander-Pacanowski mixing
      call ppmix(vbf2,shf2,anzu,vvisc_max,vvisc_min,anzt,vdiff_ts_max,vdiff_ts_min)
    end if

    if(iabs(ksw_vert)==2) then
      ! Mellor-Yamada mixing
      call Mellor_Yamada_gendis(q2,q2l,rhs_q2,rhs_q2l,coef_q2,coef_q2l,vbf2,shf2,    &
                                anzt,anzu,vvisc_min,vdiff_ts_min,vvisc_max,vdiff_ts_max,bot_fric)
    endif

endsubroutine mixing_processes

!=======================================================================
subroutine hydrodynamics(tau)

    real(4), intent(in):: tau
    integer m, n, k, ierr

    ! computing non-advective part of relative vorticity
    !$omp parallel do private(m,n,k) 
     do n=ny_start, ny_end
       do m=nx_start, nx_end
         do k=1,nz
          aux_array3d_tgr5(m,n,k)=( (vv(m+1,n,k)*dyt(m+1,n)-vv(m,n,k)*dyt(m,n))     &
                                   -(uu(m,n+1,k)*dxt(m,n+1)-uu(m,n,k)*dxt(m,n))     &
                                   -( (vv(m+1,n,k)-vv(m,n,k))*dyb(m,n)              &
                                     -(uu(m,n+1,k)-uu(m,n,k))*dxb(m,n) ) )/sqh(m,n)*luu(m,n)
         enddo
       enddo
     enddo
    !$omp end parallel do

    call syncborder_real(aux_array3d_tgr5, nz)

    if(periodicity_x/=0) then
       call cyclize_x(aux_array3d_tgr5,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
       call cyclize_y(aux_array3d_tgr5,nx,ny,nz,nnn,nn)
    end if

! storing old velocity

!$omp parallel do private(m,n,k) 
 do n=bnd_y1,bnd_y2
  do m=bnd_x1,bnd_x2  
    aux_array3d_tgr3(m,n,:) = uh(m,n,:)
    aux_array3d_tgr4(m,n,:) = vh(m,n,:)
   enddo
 enddo
!$omp end parallel do

 !Removing barotropic components from the 3D velocity
 call depth_ave(aux_array3d_tgr3,ubrtr,llu,1) 
 call depth_ave(aux_array3d_tgr4,vbrtr,llv,1)

!Preliminary vertical velocity for momentum transport
  call vertical_velocity(aux_array3d_tgr3,aux_array3d_tgr4,ww,vol_flux)

 !Updating depth functions  
 if(nonlinear_free_surface>0) then
   call ssh_pred(ubrtr,vbrtr,ssh,ssh_n,vol_flux,tau)
   call hh_update(hhq_n, hhu_n, hhv_n, hhh_n, ssh_n, hhq_rest)
 endif

!Transport-diffusion module

    !Modified Cabaret scheme
    call u_transport_cab_mod(uu,tau,aux_array3d_tgr1,aux_array3d_tgr2,aux_array3d_tgr3,aux_array3d_tgr4,     &
                             aux_array3d_wgr1,aux_array3d_wgr2,uh,vh,ww,tx_adv)
    call v_transport_cab_mod(vv,tau,aux_array3d_tgr1,aux_array3d_tgr2,aux_array3d_tgr3,aux_array3d_tgr4,     &
                             aux_array3d_wgr1,aux_array3d_wgr2,uh,vh,ww,ty_adv)

!Momentum lateral curvilinear correction    
   call uv_curv_rot(uu,vv,aux_array3d_tgr1,aux_array3d_tgr2,aux_array3d_tgr5,tau)
   
   if(iabs(ksw_lat)>0) then
!Momentum lateral harmonic mix 
    call visc2_rhs(uu,vv,stress_t,stress_s,amuv,aux_array3d_tgr1,aux_array3d_tgr2)
    call visc2_nosplit(uu,vv,aux_array3d_tgr1,aux_array3d_tgr2,tau)
   endif

   if(iabs(ksw_lat4)>0) then
!Momentum lateral biharmonic mix    
    call uv_lat_diff_4_nosplit(uu,vv,aux_array3d_tgr1,aux_array3d_tgr2,aux_array3d_tgr3,aux_array3d_tgr4,tau,stress_t,stress_s,amuv4)
   endif

!Momentum vertical mix 
   call u_vdiff(uu,tau,anzu,tx_dif,tx_coef,bot_fric)
   call v_vdiff(vv,tau,anzu,ty_dif,ty_coef,bot_fric)

 !Computing pressure gradients
  if (iabs(ksw_dens)>0) then
   call pressure_gradients_ploc(tt,ss,uu,vv,ssh,z3d,den_sg0,yng_slp,tau)
  endif

!forming full stream functions

 !$omp parallel do private(m,n,k) 
    do n = ny_start-1, ny_end+1
      do m = nx_start-1, nx_end+1
        do k = 1, nz
          uh(m,n,k)=uu(m,n,k)*hhu_n(m,n)
          vh(m,n,k)=vv(m,n,k)*hhv_n(m,n)
        enddo
      enddo
    enddo
 !$omp end parallel do

 ! Computing RHS for shallow water equations as SLP gradient
 ! (atmospheric + sea ice (if water budget is variable))

 !$omp parallel do private(m,n,k) 
 do n=ny_start,ny_end
   do m=nx_start,nx_end
     if (lcu(m,n)>0.5) then
         RHSx2d(m,n)= - ( (mistot(m+1,n)-mistot(m,n))*FreeFallAcc*float(variable_volume_budget)   &
                          + slpr(m+1,n)-slpr(m,n)   )*hhu_n(m,n)/dxt(m,n)/RefDen
     endif

     if (lcv(m,n)>0.5) then
         RHSy2d(m,n)= - ( (mistot(m,n+1)-mistot(m,n))*FreeFallAcc*float(variable_volume_budget)   &
                          + slpr(m,n+1)-slpr(m,n)   )*hhv_n(m,n)/dyt(m,n)/RefDen
     endif
   end do
 end do
 !$omp end parallel do

 !Removing barotropic components from the 3D velocity
 call depth_ave(uh,ubrtr,llu,1) 
 call depth_ave(vh,vbrtr,llv,1)

 call baroclinic_adaptation(tau,uh,vh,aux_array3d_tgr1,aux_array3d_tgr2)

!Ultimate vertical velocity for tracer transport
 call vertical_velocity(uh,vh,ww,vol_flux)
       
 !CALL FORM_RHS(UBRTR,VBRTR,SSH, hhu_n, hhu_rest,    &
 !                              hhv_n, hhv_rest, rhsx2d,rhsy2d, wf_tot, tau)                 
 !CALL SOLVE_MATRICE_ITER(UBRTR,VBRTR,SSH,niter) 
 call form_rhs(ubrtr, vbrtr, ssh, hhu_n, hhu_rest, hhv_n, hhv_rest, rhsx2d, rhsy2d, vol_flux, tau)
 call solve_system(ubrtr, vbrtr, ssh)

 !write(*,'(a,i4)') 'Number of iterations in SWEQ is', niter

 call syncborder_real8(ubrtr, 1)
 call syncborder_real8(vbrtr, 1)
 call syncborder_real8(ssh, 1)

 if(periodicity_x/=0) then
   call cyclize8_x(ubrtr,nx,ny,1,mmm,mm)
   call cyclize8_x(vbrtr,nx,ny,1,mmm,mm)
   call cyclize8_x(ssh,nx,ny,1,mmm,mm)
 end if

 if(periodicity_y/=0) then
   call cyclize8_y(ubrtr,nx,ny,1,nnn,nn)
   call cyclize8_y(vbrtr,nx,ny,1,nnn,nn)
   call cyclize8_y(ssh,nx,ny,1,nnn,nn)
 end if
 !!!!! End of barotropic equations

 !Updating depth functions  
 if(nonlinear_free_surface>0) then
   call hh_update(hhq_n, hhu_n, hhv_n, hhh_n, ssh, hhq_rest)
 endif

!Adding barotropic components to the 3D velocity
 !$omp parallel do private(m,n,k)
 do n=ny_start-1, ny_end+1
   do m=nx_start-1, nx_end+1
     if(lcu(m,n)>0.5) then
       uh(m,n,1:nz)=uh(m,n,1:nz)+ubrtr(m,n)
       uu(m,n,1:nz)=uh(m,n,1:nz)/hhu_n(m,n)
     endif

     if(lcv(m,n)>0.5) then
       vh(m,n,1:nz)=vh(m,n,1:nz)+vbrtr(m,n)
       vv(m,n,1:nz)=vh(m,n,1:nz)/hhv_n(m,n)
     endif    
   enddo
 enddo
 !$omp end parallel do

endsubroutine hydrodynamics

!=======================================================================
subroutine turbulence_tran_diff(tau)

  real(4),    intent(in):: tau    

   call turb_transport_cab_mod(q2 ,tau,aux_array3d_wgr1,aux_array3d_wgr2,      &
                                       aux_array3d_wgr3,aux_array3d_wgr4,      &
                                       aux_array3d_tgr1,aux_array3d_tgr2,uh,vh,ww)
   call turb_transport_cab_mod(q2l,tau,aux_array3d_wgr1,aux_array3d_wgr2,      &
                                       aux_array3d_wgr3,aux_array3d_wgr4,      &
                                       aux_array3d_tgr1,aux_array3d_tgr2,uh,vh,ww)
   
   call turb_sdiff(q2 ,aux_array3d_wgr1,aux_array3d_wgr2,tau,amuv,tur_factor_nu)
   call turb_sdiff(q2l,aux_array3d_wgr1,aux_array3d_wgr2,tau,amuv,tur_factor_nu)

   call turb_vdiff(q2 ,tau,anzu,tur_factor_nu,rhs_q2 ,coef_q2 )
   call turb_vdiff(q2l,tau,anzu,tur_factor_nu,rhs_q2l,coef_q2l)

endsubroutine turbulence_tran_diff

!==============================================
subroutine tem_sal_tran_diff(tau)

  real(4),    intent(in):: tau    

  integer m,n,k, ierr
   
   !Temperature and salinity transport
    !Modified Cabaret scheme
    call tracer_transport_cab_mod(tt,tau,aux_array3d_tgr1,aux_array3d_tgr2,  &
                                         aux_array3d_tgr3,aux_array3d_tgr4,  &
                                         aux_array3d_wgr1,aux_array3d_wgr2,  &
                                         uh,vh,ww,tflux_adv)

    call tracer_transport_cab_mod(ss,tau,aux_array3d_tgr1,aux_array3d_tgr2,  &
                                         aux_array3d_tgr3,aux_array3d_tgr4,  &
                                         aux_array3d_wgr1,aux_array3d_wgr2,  &
                                         uh,vh,ww,sflux_adv)

   !Temperature and salinity lateral diffusion

     if(iabs(ksw_lat)==1) then      !sigma-diffusion 
      call tracer_sdiff(tt,aux_array3d_tgr1,aux_array3d_tgr2,tau,amts,       1.0)
      call tracer_sdiff(ss,aux_array3d_tgr1,aux_array3d_tgr2,tau,amts,tsfrac_lat)
     elseif(iabs(ksw_lat)==2) then  !z-diffusion 
      call tracer_zdiff(tt,tau,amts,       1.0,z3d,aux_array3d_tgr5,       &
                                  aux_array3d_tgr4,aux_array3d_tgr3,       &
                                  aux_array3d_tgr2,aux_array3d_tgr1)
      call tracer_zdiff(ss,tau,amts,tsfrac_lat,z3d,aux_array3d_tgr5,       &
                                  aux_array3d_tgr4,aux_array3d_tgr3,       &
                                  aux_array3d_tgr2,aux_array3d_tgr1)     
     elseif(iabs(ksw_lat)==3) then  !rho-diffusion 
      call tracer_isopyc(tt,tau,amts,       1.0,den_isopyc,z3d,      &
            aux_array3d_tgr5,aux_array3d_tgr4,aux_array3d_tgr3,      &
            aux_array3d_tgr2,aux_array3d_tgr1)
      call tracer_isopyc(ss,tau,amts,tsfrac_lat,den_isopyc,z3d,      &
            aux_array3d_tgr5,aux_array3d_tgr4,aux_array3d_tgr3,      &
            aux_array3d_tgr2,aux_array3d_tgr1)
     endif

! Computing penetrating radiation influence

!$omp parallel do private(m, n, k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end   
     if(lu(m,n)>0.5) then
      do k=1,nz
       tt(m,n,k) = tt(m,n,k) + tau*swflux(m,n)*divswrad(m,n,k)/dz(k)/hhq_n(m,n)
      enddo
     endif
    enddo
   enddo
!$omp end parallel do  

   !Temperature and salinity vertical diffusion

     call tracer_vdiff(tt,tau,anzt,        1.0,tflux_dif,igrzts)
     call tracer_vdiff(ss,tau,anzt,tsfrac_vert,sflux_dif,igrzts)

endsubroutine tem_sal_tran_diff

!==============================================
subroutine ideal_age_tran_diff(tau)

  real(4),    intent(in):: tau    

  integer m,n,k

   !Ideal age transport
    !Modified Cabaret scheme
    call tracer_transport_cab_mod(age,tau,aux_array3d_tgr1,aux_array3d_tgr2,  &
                                          aux_array3d_tgr3,aux_array3d_tgr4,  &
                                          aux_array3d_wgr1,aux_array3d_wgr2,  &
                                          uh,vh,ww,age_forc_adv)
   !Ideal age lateral diffusion

     if(iabs(ksw_lat)==1) then      !sigma-diffusion 
      call tracer_sdiff(age,aux_array3d_tgr1,aux_array3d_tgr2,tau,amts,1.0)
     elseif(iabs(ksw_lat)==2) then  !z-diffusion 
      call tracer_zdiff(age,tau,amts,1.0,z3d,aux_array3d_tgr5,       &
                            aux_array3d_tgr4,aux_array3d_tgr3,       &
                            aux_array3d_tgr2,aux_array3d_tgr1)  
     elseif(iabs(ksw_lat)==3) then  !rho-diffusion 
      call tracer_isopyc(age,tau,amts,1.0,den_isopyc,z3d,      &
      aux_array3d_tgr5,aux_array3d_tgr4,aux_array3d_tgr3,      &
      aux_array3d_tgr2,aux_array3d_tgr1)
     endif

! Computing age olding (in years)

!$omp parallel do private(m, n, k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end   
     if(lu(m,n)>0.5) then
      do k=1,nz
       age(m,n,k) = age(m,n,k) + tau/86400.0/365.2425
      enddo
     endif
    enddo
   enddo
!$omp end parallel do  

   !Age vertical diffusion

     call tracer_vdiff(age,tau,anzt,1.0,age_forc_dif,igrzage)

endsubroutine ideal_age_tran_diff

!=======================================================================
subroutine passtracer_tran_diff(tau)

  real(4),    intent(in):: tau    

  integer m,n,k

   !Passive tracer transport

    call tracer_transport_cab_mod(pass_tracer,tau,aux_array3d_tgr1,aux_array3d_tgr2,  &
                                                  aux_array3d_tgr3,aux_array3d_tgr4,  &
                                                  aux_array3d_wgr1,aux_array3d_wgr2,  &
                                                  uh,vh,ww,pt_forc_adv)

   !Passive tracer lateral diffusion

     if(iabs(ksw_lat)==1) then      !sigma-diffusion 
      call tracer_sdiff(pass_tracer,aux_array3d_tgr1,aux_array3d_tgr2,tau,amts,1.0)
     elseif(iabs(ksw_lat)==2) then  !z-diffusion 
      call tracer_zdiff(pass_tracer,tau,amts,1.0,z3d,aux_array3d_tgr5,       &
                                    aux_array3d_tgr4,aux_array3d_tgr3,       &
                                    aux_array3d_tgr2,aux_array3d_tgr1)  
     elseif(iabs(ksw_lat)==3) then  !rho-diffusion 
      call tracer_isopyc(pass_tracer,tau,amts,1.0,den_isopyc,z3d,      &
              aux_array3d_tgr5,aux_array3d_tgr4,aux_array3d_tgr3,      &
              aux_array3d_tgr2,aux_array3d_tgr1)
     endif

   !Passive tracer vertical diffusion

     call tracer_vdiff(pass_tracer,tau,anzt,1.0,pt_forc_dif,igrzpt)

endsubroutine passtracer_tran_diff

endmodule ocstep_routes
