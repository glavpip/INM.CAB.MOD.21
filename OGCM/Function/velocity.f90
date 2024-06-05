module velocity_routes
use constants
use basin_grid
use ocalg_routes
use mixing_routes
implicit none

contains

!====================================================================
subroutine u_transport_cab_mod(ff,tau,fx,fx0,fy,fy0,fz,fz0,uh,vh,ww,flx_top)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: ww
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: flx_top
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ):: fx,fx0,fy,fy0
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1):: fz,fz0

integer m,n,k, ind1,ind2
real(8) hmid, rhs, up, um, vp, vm, flxz(nz+1)
real(4) u, v, w
real(8) cvert(nz+1)

cvert(1) = 0.0
cvert(nz+1) = 0.0
cvert(2:nz) = 1.0

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx(m,n,:) = 0.0
      fy(m,n,:) = 0.0
      fz(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!Predictor step

!$omp parallel do private(m,n,k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lu(m,n)>0.5) then
       do k=1,nz
        fx(m,n,k) = (ff(m-1,n,k) + ff(m,n,k))/2.0
       enddo
      endif
            
      if (luu(m,n)>0.5) then
       do k=1,nz
        fy(m,n,k) = (ff(m,n,k) + ff(m,n+1,k))/2.0
       enddo
      endif

      if (lcu(m,n)>0.5) then
       do k=2,nz
        fz(m,n,k) = (ff(m,n,k-1) + ff(m,n,k))/2.0
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz)
    call syncborder_real(fy, nz)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz,nnn,nn)            
    end if

!$omp parallel do private(m,n,k,hmid,rhs,flxz,up,um,vp,vm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lcu(m,n)>0.5) then
       hmid=(hhu(m,n)+hhu_n(m,n))/2.0
       
       flxz(1)=flx_top(m,n)*squ(m,n)
       do k=2,nz
       flxz(k) = ( ww(m,n,k)*sqt(m,n) + ww(m+1,n,k)*sqt(m+1,n) )/2.0 * fz(m,n,k)
       enddo
       flxz(nz+1)=0.0
       
       do k=1,nz
        up=(uh(m,n,k)*dyh(m,n) + uh(m+1,n,k)*dyh(m+1,n))/2.0
        um=(uh(m,n,k)*dyh(m,n) + uh(m-1,n,k)*dyh(m-1,n))/2.0
        vp=(vh(m,n  ,k)*dxh(m,n  ) + vh(m+1,n  ,k)*dxh(m+1,n  ))/2.0*luu(m,n  )
        vm=(vh(m,n-1,k)*dxh(m,n-1) + vh(m+1,n-1,k)*dxh(m+1,n-1))/2.0*luu(m,n-1)

        rhs=up*fx(m+1,n,k) - um*fx(m,n,k) + vp*fy(m,n,k) - vm*fy(m,n-1,k)     &
             +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hhu(m,n)/hmid) - rhs*tau/2.0/squ(m,n)/hmid
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz,nnn,nn)
     end if

!Crossbeam step

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx0(m,n,:) = fx(m,n,:)
      fy0(m,n,:) = fy(m,n,:)
      fz0(m,n,:) = fz(m,n,:)
    enddo
   enddo
!$omp end parallel do  


!$omp parallel do private(m,n,k,ind1,ind2,u,v,w) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (lu(m,n)>0.5) then
       do k=1,nz
        u=(uh(m,n,k)*dyh(m,n) + uh(m-1,n,k)*dyh(m-1,n))/2.0        
        ind1= nint( - 0.5 - sign(0.5,u) )
        ind2= nint(  - sign(1.0,u) )
        fx(m,n,k) = 2.0*ff(m+ind1,n,k) - fx0(m+ind2,n,k)
       enddo
      endif

      if (luu(m,n)>0.5) then
       do k=1,nz
        v=(vh(m,n  ,k)*dxh(m,n  ) + vh(m+1,n  ,k)*dxh(m+1,n  ))/2.0
        ind1= nint( 0.5 - sign(0.5,v) )
        ind2= nint(  - sign(1.0,v) )
        fy(m,n,k) = (1.0 + luu(m,n+ind2))*ff(m,n+ind1,k) - luu(m,n+ind2)*fy0(m,n+ind2,k)
       enddo
      endif

      if (lcu(m,n)>0.5) then
       do k=2,nz
        w = ( ww(m,n,k)*sqt(m,n) + ww(m+1,n,k)*sqt(m+1,n) )/2.0
        ind1= nint( - 0.5 - sign(0.5,w) )
        ind2= nint(  - sign(1.0,w) )
        fz(m,n,k) = (1.0 + cvert(k+ind2))*ff(m,n,k+ind1) - cvert(k+ind2)*fz0(m,n,k+ind2)
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz)
    call syncborder_real(fy, nz)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz,nnn,nn)            
    end if

!Corrector step

!$omp parallel do private(m,n,k,hmid,rhs,flxz,up,um,vp,vm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lcu(m,n)>0.5) then
       hmid=(hhu(m,n)+hhu_n(m,n))/2.0
       
       flxz(1)=flx_top(m,n)*squ(m,n)
       do k=2,nz
       flxz(k) = ( ww(m,n,k)*sqt(m,n) + ww(m+1,n,k)*sqt(m+1,n) )/2.0 * fz(m,n,k)
       enddo
       flxz(nz+1)=0.0
       
       do k=1,nz
        up=(uh(m,n,k)*dyh(m,n) + uh(m+1,n,k)*dyh(m+1,n))/2.0
        um=(uh(m,n,k)*dyh(m,n) + uh(m-1,n,k)*dyh(m-1,n))/2.0
        vp=(vh(m,n  ,k)*dxh(m,n  ) + vh(m+1,n  ,k)*dxh(m+1,n  ))/2.0*luu(m,n  )
        vm=(vh(m,n-1,k)*dxh(m,n-1) + vh(m+1,n-1,k)*dxh(m+1,n-1))/2.0*luu(m,n-1)

        rhs=up*fx(m+1,n,k) - um*fx(m,n,k) + vp*fy(m,n,k) - vm*fy(m,n-1,k)     &
             +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hmid/hhu_n(m,n)) - rhs*tau/2.0/squ(m,n)/hhu_n(m,n)
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz,nnn,nn)
     end if

endsubroutine u_transport_cab_mod

!====================================================================
subroutine v_transport_cab_mod(ff,tau,fx,fx0,fy,fy0,fz,fz0,uh,vh,ww,flx_top)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: ww
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: flx_top
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ):: fx,fx0,fy,fy0
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1):: fz,fz0

integer m,n,k, ind1,ind2
real(8) hmid, rhs, up, um, vp, vm, flxz(nz+1)
real(4) u, v, w 
real(8) cvert(nz+1)

cvert(1) = 0.0
cvert(nz+1) = 0.0
cvert(2:nz) = 1.0

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx(m,n,:) = 0.0
      fy(m,n,:) = 0.0
      fz(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!Predictor step

!$omp parallel do private(m,n,k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (luu(m,n)>0.5) then
       do k=1,nz
        fx(m,n,k) = (ff(m,n,k) + ff(m+1,n,k))/2.0
       enddo
      endif
      
      if (lu(m,n)>0.5) then
       do k=1,nz
        fy(m,n,k) = (ff(m,n-1,k) + ff(m,n,k))/2.0
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=2,nz
        fz(m,n,k) = (ff(m,n,k-1) + ff(m,n,k))/2.0
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz)
    call syncborder_real(fy, nz)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz,nnn,nn)            
    end if

!$omp parallel do private(m,n,k,hmid,rhs,flxz,up,um,vp,vm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lcv(m,n)>0.5) then
       hmid=(hhv(m,n)+hhv_n(m,n))/2.0
       
       flxz(1)=flx_top(m,n)*sqv(m,n)
       do k=2,nz
       flxz(k) = ( ww(m,n,k)*sqt(m,n) + ww(m,n+1,k)*sqt(m,n+1) )/2.0 * fz(m,n,k)
       enddo
       flxz(nz+1)=0.0
       
       do k=1,nz
        up=(uh(m  ,n,k)*dyh(m  ,n) + uh(m  ,n+1,k)*dyh(m  ,n+1))/2.0*luu(m  ,n)
        um=(uh(m-1,n,k)*dyh(m-1,n) + uh(m-1,n+1,k)*dyh(m-1,n+1))/2.0*luu(m-1,n)
        vp=(vh(m,n,k)*dxh(m,n) + vh(m,n+1,k)*dxh(m,n+1))/2.0
        vm=(vh(m,n,k)*dxh(m,n) + vh(m,n-1,k)*dxh(m,n-1))/2.0

        rhs=up*fx(m,n,k) - um*fx(m-1,n,k) + vp*fy(m,n+1,k) - vm*fy(m,n,k)     &
             +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hhv(m,n)/hmid) - rhs*tau/2.0/sqv(m,n)/hmid
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz,nnn,nn)
     end if

!Crossbeam step

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx0(m,n,:) = fx(m,n,:)
      fy0(m,n,:) = fy(m,n,:)
      fz0(m,n,:) = fz(m,n,:)
    enddo
   enddo
!$omp end parallel do  


!$omp parallel do private(m,n,k,ind1,ind2,u,v,w) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (luu(m,n)>0.5) then
       do k=1,nz
        u =(uh(m  ,n,k)*dyh(m  ,n) + uh(m  ,n+1,k)*dyh(m  ,n+1))/2.0   
        ind1= nint( 0.5 - sign(0.5,u) )
        ind2= nint(  - sign(1.0,u) )
        fx(m,n,k) = (1.0 + luu(m+ind2,n))*ff(m+ind1,n,k) - luu(m+ind2,n)*fx0(m+ind2,n,k)
       enddo
      endif

      if (lu(m,n)>0.5) then
       do k=1,nz
        v=(vh(m,n,k)*dxh(m,n) + vh(m,n-1,k)*dxh(m,n-1))/2.0
        ind1= nint( - 0.5 - sign(0.5,v) )
        ind2= nint(  - sign(1.0,v) )
        fy(m,n,k) = 2.0*ff(m,n+ind1,k) - fy0(m,n+ind2,k)
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=2,nz
        w = ( ww(m,n,k)*sqt(m,n) + ww(m,n+1,k)*sqt(m,n+1) )/2.0
        ind1= nint( - 0.5 - sign(0.5,w) )
        ind2= nint(  - sign(1.0,w) )
        fz(m,n,k) = (1.0 + cvert(k+ind2))*ff(m,n,k+ind1) - cvert(k+ind2)*fz0(m,n,k+ind2)
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz)
    call syncborder_real(fy, nz)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz,nnn,nn)            
    end if

!Corrector step

!$omp parallel do private(m,n,k,hmid,rhs,flxz,up,um,vp,vm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lcv(m,n)>0.5) then
       hmid=(hhv(m,n)+hhv_n(m,n))/2.0
       
       flxz(1)=flx_top(m,n)*sqv(m,n)
       do k=2,nz
       flxz(k) = ( ww(m,n,k)*sqt(m,n) + ww(m,n+1,k)*sqt(m,n+1) )/2.0 * fz(m,n,k)
       enddo
       flxz(nz+1)=0.0
       
       do k=1,nz
        up=(uh(m  ,n,k)*dyh(m  ,n) + uh(m  ,n+1,k)*dyh(m  ,n+1))/2.0*luu(m  ,n)
        um=(uh(m-1,n,k)*dyh(m-1,n) + uh(m-1,n+1,k)*dyh(m-1,n+1))/2.0*luu(m-1,n)
        vp=(vh(m,n,k)*dxh(m,n) + vh(m,n+1,k)*dxh(m,n+1))/2.0
        vm=(vh(m,n,k)*dxh(m,n) + vh(m,n-1,k)*dxh(m,n-1))/2.0

        rhs=up*fx(m,n,k) - um*fx(m-1,n,k) + vp*fy(m,n+1,k) - vm*fy(m,n,k)     &
             +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hmid/hhv_n(m,n)) - rhs*tau/2.0/sqv(m,n)/hhv_n(m,n)
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz,nnn,nn)
     end if

endsubroutine v_transport_cab_mod

!============================================================
subroutine uv_curv_rot(uu,vv,u0,v0,vort,tau)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: vort
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: uu,vv
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: u0,v0

 integer m,n,k
 real(8) denom

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
     u0(m,n,:) = uu(m,n,:)
     v0(m,n,:) = vv(m,n,:)
    enddo
   enddo
!$omp end parallel do   

!$omp parallel do private(m,n,k, denom)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
         
         if(lcu(m,n)>0.5) then
          do k=1,nz
           denom= 1.0 + tau**2 * (vort(m,n,k)**2+vort(m,n-1,k)**2)/2.0
           uu(m,n,k)=( u0(m,n,k) + (vort(m,n  ,k)*(v0(m  ,n  ,k)+v0(m+1,n  ,k))       &
                                   +vort(m,n-1,k)*(v0(m  ,n-1,k)+v0(m+1,n-1,k)) )*tau/4.0  ) / denom
          enddo
         endif

         if(lcv(m,n)>0.5) then
          do k=1,nz
           denom= 1.0 + tau**2 * (vort(m,n,k)**2+vort(m-1,n,k)**2)/2.0
           vv(m,n,k)=( v0(m,n,k) - (vort(m  ,n,k)*(u0(m  ,n  ,k)+u0(m  ,n+1,k))       &
                                   +vort(m-1,n,k)*(u0(m-1,n  ,k)+u0(m-1,n+1,k)) )*tau/4.0 ) / denom
          enddo
         endif

       enddo
      enddo
!$omp end parallel do

    call syncborder_real(uu, nz)
    call syncborder_real(vv, nz)

    if(periodicity_x/=0) then
      call cyclize_x(uu,nx,ny,nz,mmm,mm)
      call cyclize_x(vv,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(uu,nx,ny,nz,nnn,nn)
      call cyclize_y(vv,nx,ny,nz,nnn,nn)
    end if

endsubroutine uv_curv_rot

!============================================================
subroutine visc2_rhs(uu,vv,str_t,str_s,mu,rhsu,rhsv)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in ):: uu,vv,mu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: rhsu,rhsv
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: str_t,str_s

 integer m,n,k
 real(8) muh_p, muh_m, rhs

 !$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
     rhsu(m,n,:) = 0.0
     rhsv(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do   

!computing deformation rates
  call stress_tension(uu,vv,str_t)
  call stress_shear  (uu,vv,str_s)

!$omp parallel do private(m,n,k,muh_p,muh_m,rhs)      
  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity      
      if(lcu(m,n)>0.5) then
        do k=1,nz
         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0
         muh_m=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n-1,k)+mu(m+1,n-1,k))/4.0   
         rhs = ( dy(m+1,n)**2*mu(m+1,n,k)*hhq_n(m+1,n)*str_t(m+1,n,k)                         &
                -dy(m  ,n)**2*mu(m  ,n,k)*hhq_n(m  ,n)*str_t(m  ,n,k) )/(dyh(m,n)*squ(m,n))   &
             + (dxb(m,n  )**2*muh_p*hhh_n(m,n  )*str_s(m,n  ,k)                               &
               -dxb(m,n-1)**2*muh_m*hhh_n(m,n-1)*str_s(m,n-1,k) )/(dxt(m,n)*squ(m,n))                     
         rhsu(m,n,k)=rhs/hhu_n(m,n)
        end do
      end if

!meridional velocity     
      if(lcv(m,n)>0.5) then
        do k=1,nz
         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0
         muh_m=(mu(m,n,k)+mu(m-1,n,k)+mu(m,n+1,k)+mu(m-1,n+1,k))/4.0        
         rhs =-( dx(m,n+1)**2*mu(m,n+1,k)*hhq_n(m,n+1)*str_t(m,n+1,k)                         &
                -dx(m,n  )**2*mu(m,n  ,k)*hhq_n(m,n  )*str_t(m,n  ,k) )/(dxh(m,n)*sqv(m,n))   &
             + (dyb(m  ,n)**2*muh_p*hhh_n(m  ,n)*str_s(m  ,n,k)                               &
               -dyb(m-1,n)**2*muh_m*hhh_n(m-1,n)*str_s(m-1,n,k) ) /(dyt(m,n)*sqv(m,n))  
         rhsv(m,n,k)=rhs/hhv_n(m,n)     
        end do
      end if

    end do
  end do
!$omp end parallel do

    call syncborder_real(rhsu, nz)
    call syncborder_real(rhsv, nz)

    if(periodicity_x/=0) then
       call cyclize_x(rhsu,nx,ny,nz,mmm,mm)
       call cyclize_x(rhsv,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
       call cyclize_y(rhsu,nx,ny,nz,nnn,nn)
       call cyclize_y(rhsv,nx,ny,nz,nnn,nn)
    end if

endsubroutine visc2_rhs

!============================================================
subroutine visc2_nosplit(uu,vv,rhsu,rhsv,tau)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in ):: rhsu,rhsv 
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: uu,vv
real(4), intent(in):: tau

 integer m,n
 
!$omp parallel do private(m,n)      
  do n=bnd_y1,bnd_y2
   do m=bnd_x1,bnd_x2 
     uu(m,n,1:nz)=uu(m,n,1:nz) + tau*rhsu(m,n,1:nz)
     vv(m,n,1:nz)=vv(m,n,1:nz) + tau*rhsv(m,n,1:nz)
    end do
  end do
!$omp end parallel do

endsubroutine visc2_nosplit

!====================================================
subroutine uv_lat_diff_4_nosplit(uu,vv,u0,v0,u1,v1,tau,str_t,str_s,mu4)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in ):: mu4
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: uu,vv
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ):: u0,v0,u1,v1,str_t,str_s
real(4), intent(in):: tau

 integer m,n,k, ierr
 real(8) rhs

!second order viscosity operation computation
 call visc2_rhs(uu,vv,str_t,str_s,mu4,u0,v0)

!fourth order viscosity operation computation
 call visc2_rhs(u0,v0,str_t,str_s,mu4,u1,v1)

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
     uu(m,n,:) = uu(m,n,:)-u1(m,n,:)*tau
     vv(m,n,:) = vv(m,n,:)-v1(m,n,:)*tau
    enddo
   enddo
!$omp end parallel do

endsubroutine uv_lat_diff_4_nosplit

!================================================================================
subroutine pressure_gradients_ploc(tt,ss,uu,vv,ssh,z3d,den_s,k_s,tau)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in ):: tt,ss,z3d,den_s,k_s
real(8), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2   ), intent(in ):: ssh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: uu,vv
real(4), intent(in):: tau

 real(8) p_hint, p_ssh, rhs, z_pp, z_pm, z_mp, z_mm, rho_pm, rho_mm,  rho_pp, rho_mp
 real(4) pgrz
 integer m,n,k

 !$omp parallel do private(m,n,k,p_hint,p_ssh,pgrz,rhs,z_pp,z_pm,z_mp,z_mm,rho_pm,rho_mm,rho_pp,rho_mp)      
 do n=ny_start,ny_end
   do m=nx_start,nx_end

     if(lcu(m,n)>0.5) then
          pgrz=FreeFallAcc*RefDen*hhu(m,n)
          rhs= - FreeFallAcc*tau/RefDen/2.0/dxt(m,n)

          z_pm= 0.0
          z_mm= 0.0
          z_pp= z3d(m+1,n,1)
          z_mp= z3d(m  ,n,1)

          rho_pp = den_s(m+1,n,1)
          rho_mp = den_s(m  ,n,1)

          rho_pm = rho_pp
          rho_mm = rho_mp

          p_hint= (rho_pp-rho_mm)*(z_mp-z_pm) - (z_pp-z_mm)*(rho_mp-rho_pm)
          p_ssh = (ssh(m+1,n)-ssh(m,n))*(rho_pp+rho_mp+(2000.0-2.0*RefDen))
          uu(m,n,1)=uu(m,n,1)+rhs*(p_hint+p_ssh)
       
         do k=2,nz
          z_pm= z_pp
          z_mm= z_mp
          z_pp= z3d(m+1,n,k)
          z_mp= z3d(m  ,n,k)

          rho_pm= density_unesco_tpot(tt(m+1,n,k-1),ss(m+1,n,k-1)+salref,pgrz*zw(k),den_s(m+1,n,k-1),k_s(m+1,n,k-1))
          rho_mm= density_unesco_tpot(tt(m  ,n,k-1),ss(m  ,n,k-1)+salref,pgrz*zw(k),den_s(m  ,n,k-1),k_s(m  ,n,k-1))
          rho_pp= density_unesco_tpot(tt(m+1,n,k  ),ss(m+1,n,k  )+salref,pgrz*zw(k),den_s(m+1,n,k  ),k_s(m+1,n,k  ))
          rho_mp= density_unesco_tpot(tt(m  ,n,k  ),ss(m  ,n,k  )+salref,pgrz*zw(k),den_s(m  ,n,k  ),k_s(m  ,n,k  ))

          p_hint= p_hint + (rho_pp-rho_mm)*(z_mp-z_pm) - (z_pp-z_mm)*(rho_mp-rho_pm)
          p_ssh = (ssh(m+1,n)-ssh(m,n))*(rho_pp+rho_mp+(2000.0-2.0*RefDen))
          uu(m,n,k)=uu(m,n,k)+rhs*(p_hint+p_ssh)
         enddo
     endif

     if(lcv(m,n)>0.5) then
          pgrz=FreeFallAcc*RefDen*hhv(m,n)
          rhs= - FreeFallAcc*tau/RefDen/2.0/dyt(m,n)

          z_pm= 0.0
          z_mm= 0.0
          z_pp= z3d(m,n+1,1)
          z_mp= z3d(m,n  ,1)

          rho_pp = den_s(m,n+1,1)
          rho_mp = den_s(m,n  ,1)
          
          rho_pm = rho_pp
          rho_mm = rho_mp
          
          p_hint= (rho_pp-rho_mm)*(z_mp-z_pm) - (z_pp-z_mm)*(rho_mp-rho_pm)
          p_ssh = (ssh(m,n+1)-ssh(m,n))*(rho_pp+rho_mp+(2000.0-2.0*RefDen))
          vv(m,n,1)=vv(m,n,1)+rhs*(p_hint+p_ssh)
       
         do k=2,nz
          z_pm= z_pp
          z_mm= z_mp
          z_pp= z3d(m,n+1,k)
          z_mp= z3d(m,n  ,k)

          rho_pm= density_unesco_tpot(tt(m,n+1,k-1),ss(m,n+1,k-1)+salref,pgrz*zw(k),den_s(m,n+1,k-1),k_s(m,n+1,k-1))
          rho_mm= density_unesco_tpot(tt(m,n  ,k-1),ss(m,n  ,k-1)+salref,pgrz*zw(k),den_s(m,n  ,k-1),k_s(m,n  ,k-1))
          rho_pp= density_unesco_tpot(tt(m,n+1,k  ),ss(m,n+1,k  )+salref,pgrz*zw(k),den_s(m,n+1,k  ),k_s(m,n+1,k  ))
          rho_mp= density_unesco_tpot(tt(m,n  ,k  ),ss(m,n  ,k  )+salref,pgrz*zw(k),den_s(m,n  ,k  ),k_s(m,n  ,k  ))

          p_hint= p_hint + (rho_pp-rho_mm)*(z_mp-z_pm) - (z_pp-z_mm)*(rho_mp-rho_pm)
          p_ssh = (ssh(m,n+1)-ssh(m,n))*(rho_pp+rho_mp+(2000.0-2.0*RefDen))
          vv(m,n,k)=vv(m,n,k)+rhs*(p_hint+p_ssh)
         enddo
     endif

   enddo
 enddo
 !$omp end parallel do

   call syncborder_real(uu, nz)
   call syncborder_real(vv, nz)
  
   if(periodicity_x/=0) then
     call cyclize_x(uu,nx,ny,nz,mmm,mm)
     call cyclize_x(vv,nx,ny,nz,mmm,mm)
   end if
  
   if(periodicity_y/=0) then
     call cyclize_y(uu,nx,ny,nz,nnn,nn)
     call cyclize_y(vv,nx,ny,nz,nnn,nn)
   end if

endsubroutine pressure_gradients_ploc

!=====================================================================
subroutine baroclinic_adaptation(tau,uh,vh,u0,v0)

real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: u0,v0

real(8) denom
integer m,n,k,ierr

!$omp parallel do private(m,n,k)
      do n=bnd_y1,bnd_y2
       do m=bnd_x1,bnd_x2
          u0(m,n,1:nz)=uh(m,n,1:nz)
          v0(m,n,1:nz)=vh(m,n,1:nz)
       enddo
      enddo
!$omp end parallel do

!$omp parallel do private(m,n,k, denom)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
         
         if(lcu(m,n)>0.5) then
          denom= 1.0 + tau**2 * (rlh_s(m,n)**2+rlh_s(m,n-1)**2)/2.0
          do k=1,nz
           uh(m,n,k)=( u0(m,n,k) + ( rlh_s(m,n  )*(v0(m  ,n  ,k)+v0(m+1,n  ,k))       &
                                    +rlh_s(m,n-1)*(v0(m  ,n-1,k)+v0(m+1,n-1,k)) )     &
                                   *tau/4.0 ) / denom
          enddo
         endif

         if(lcv(m,n)>0.5) then
          denom= 1.0 + tau**2 * (rlh_s(m,n)**2+rlh_s(m-1,n)**2)/2.0
          do k=1,nz
           vh(m,n,k)=( v0(m,n,k) - ( rlh_s(m  ,n)*(u0(m  ,n  ,k)+u0(m  ,n+1,k))       &
                                    +rlh_s(m-1,n)*(u0(m-1,n  ,k)+u0(m-1,n+1,k)) )     &
                                   *tau/4.0 ) / denom
          enddo
         endif

       enddo
      enddo
!$omp end parallel do

    call syncborder_real(uh, nz)
    call syncborder_real(vh, nz)

    if(periodicity_x/=0) then
      call cyclize_x(uh,nx,ny,nz,mmm,mm)
      call cyclize_x(vh,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(uh,nx,ny,nz,nnn,nn)
      call cyclize_y(vh,nx,ny,nz,nnn,nn)
    end if

endsubroutine baroclinic_adaptation

!======================================================================
subroutine vertical_velocity(uh,vh,ww,vflux)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: vflux
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: ww

!----------------------------------------------------------------------
! computing w from continuity equation on grid "c".
! input:
! u - zonal velocity on uc-grid
! v - meridional velocity on vc-grid
! output:
! w - vertical velocity on t-grid

integer m, n, k

!$omp parallel do private(m,n,k)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
         ww(m,n,nz+1)=0.0
         do k=nz,1,-1
          ww(m,n,k) = ww(m,n,k+1) +dz(k)*( ( uh(m,n,k)*dyh(m,n) - uh(m-1,n,k)*dyh(m-1,n)                  &
                                           + vh(m,n,k)*dxh(m,n) - vh(m,n-1,k)*dxh(m,n-1))/sqt(m,n) + vflux(m,n))
         enddo
        endif
       enddo
      enddo
!$omp end parallel do

    call syncborder_real(ww, nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(ww,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(ww,nx,ny,nz+1,nnn,nn)
    end if
      
endsubroutine vertical_velocity

!======================================================================
subroutine ssh_pred(ubh,vbh,ssh,ssh_n,vflux,tau)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2), intent(in):: vflux
real(8), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2), intent(in):: ubh,vbh,ssh
real(4), intent(in):: tau
real(8), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2), intent(inout):: ssh_n

integer m, n

!$omp parallel do private(m,n)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
          ssh_n(m,n) = ssh(m,n)+tau*(vflux(m,n)                             &
               -( ubh(m,n)*dyh(m,n) - ubh(m-1,n)*dyh(m-1,n)                 &
                + vbh(m,n)*dxh(m,n) - vbh(m,n-1)*dxh(m,n-1))/sqt(m,n) )         
        endif
       enddo
      enddo
!$omp end parallel do

    call syncborder_real8(ssh_n,1)

    if(periodicity_x/=0) then
      call cyclize8_x(ssh_n,nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize8_y(ssh_n,nx,ny,1,nnn,nn)
    end if
      
endsubroutine ssh_pred

!====================================================================
subroutine u_vdiff(ff,tau,nu,rhs,coef,botfr)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: nu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: rhs,coef,botfr
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz) 
real(8) bp, dp, dm

!$omp parallel do private(m,n,k,bp,dp,dm,a,b,c,eta,rksi)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lcu(m,n)>0.5) then
         bp = hhu_n(m,n)/tau

!Upper layer          
          k=1

          dm = coef(m,n)
          dp = (nu(m,n,k+1)+nu(m+1,n,k+1))/2.0/hhu_n(m,n)/hzt(k+1)

          c(k) = -dp
          a(k) = 0.0
          b(k) = bp*dz(k) + dp + dm
          eta(k) = bp*ff(m,n,k)*dz(k) + rhs(m,n)

! internal points.
         do k=2,nz-1

          dm = dp
          dp = (nu(m,n,k+1)+nu(m+1,n,k+1))/2.0/hhu_n(m,n)/hzt(k+1)

          c(k) = -dp
          a(k) = -dm
          b(k) =  bp*dz(k) + dp + dm
          eta(k) = bp*ff(m,n,k)*dz(k)
         enddo

!Bottom layer
          k=nz

          dm = dp
          dp = (botfr(m,n)+botfr(m+1,n))/2.0

          c(k) = 0.0
          a(k) = -dm
          b(k) =  bp*dz(k) + dm + dp
          eta(k) = bp*ff(m,n,k)*dz(k)

         call factor8(nz,a,b,c,eta,rksi,1,nz)
         do k=1,nz
          ff(m,n,k)=rksi(k)
         enddo
        endif
       enddo
      enddo
!$omp end parallel do

    call syncborder_real(ff, nz)

    if(periodicity_x/=0) then
      call cyclize_x(ff,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(ff,nx,ny,nz,nnn,nn)
    end if

endsubroutine u_vdiff

!====================================================================
subroutine v_vdiff(ff,tau,nu,rhs,coef,botfr)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: nu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: rhs,coef,botfr
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz) 
real(8) bp, dp, dm

!$omp parallel do private(m,n,k,bp,dp,dm,a,b,c,eta,rksi)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lcv(m,n)>0.5) then
         bp = hhv_n(m,n)/tau

!Upper layer          
          k=1

          dm = coef(m,n)
          dp = (nu(m,n,k+1)+nu(m,n+1,k+1))/2.0/hhv_n(m,n)/hzt(k+1)

          c(k) = -dp 
          a(k) = 0.0
          b(k) = bp*dz(k) + dp + dm
          eta(k) = bp*ff(m,n,k)*dz(k) + rhs(m,n) 

! internal points.
         do k=2,nz-1

          dm = dp
          dp = (nu(m,n,k+1)+nu(m,n+1,k+1))/2.0/hhv_n(m,n)/hzt(k+1)

          c(k) = -dp
          a(k) = -dm
          b(k) =  bp*dz(k) + dp + dm
          eta(k) = bp*ff(m,n,k)*dz(k)
         enddo

!Bottom layer
          k=nz

          dm = dp
          dp = (botfr(m,n)+botfr(m,n+1))/2.0

          c(k) = 0.0
          a(k) = -dm
          b(k) =  bp*dz(k) + dm + dp
          eta(k) = bp*ff(m,n,k)*dz(k)

         call factor8(nz,a,b,c,eta,rksi,1,nz)
         do k=1,nz
          ff(m,n,k)=rksi(k)
         enddo
        endif
       enddo
      enddo
!$omp end parallel do

    call syncborder_real(ff, nz)

    if(periodicity_x/=0) then
      call cyclize_x(ff,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(ff,nx,ny,nz,nnn,nn)
    end if

endsubroutine v_vdiff

endmodule velocity_routes
