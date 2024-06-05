module ts_routes
use constants
use basin_grid
use ocalg_routes
implicit none

contains

!====================================================================
subroutine tracer_transport_cab_mod(ff,tau,fx,fx0,fy,fy0,fz,fz0,uh,vh,ww,flx_top)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: ww
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: flx_top
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ):: fx,fx0,fy,fy0
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1):: fz,fz0

integer m,n,k,ind1,ind2 
real(8) hmid, rhs, flxz(nz+1), cvert(nz+1)

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
      
      if (lcu(m,n)>0.5) then
       do k=1,nz
        fx(m,n,k) = (ff(m,n,k) + ff(m+1,n,k))/2.0
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=1,nz
        fy(m,n,k) = (ff(m,n,k) + ff(m,n+1,k))/2.0
       enddo
      endif

      if (lu(m,n)>0.5) then
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

!$omp parallel do private(m,n,k,hmid,rhs,flxz) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lu(m,n)>0.5) then
       
        hmid=(hhq(m,n)+hhq_n(m,n))/2.0
        flxz(1)=flx_top(m,n)
       
       do k=2,nz
        flxz(k) = ww(m,n,k)*fz(m,n,k)
       enddo
        flxz(nz+1)=0.0

       do k=1,nz
        rhs=(uh(m,n,k)*dyh(m,n)*fx(m,n,k) - uh(m-1,n,k)*dyh(m-1,n)*fx(m-1,n,k)            &
            +vh(m,n,k)*dxh(m,n)*fy(m,n,k) - vh(m,n-1,k)*dxh(m,n-1)*fy(m,n-1,k))/sqt(m,n)  &
            +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hhq(m,n)/hmid) - rhs*tau/2.0/hmid
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


!$omp parallel do private(m,n,k,ind1,ind2) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (lcu(m,n)>0.5) then
       do k=1,nz
        ind1= nint( 0.5 - sign(0.5,uh(m,n,k)) )
        ind2= nint(  - sign(1.0,uh(m,n,k)) )
        fx(m,n,k) = (1.0 + lcu(m+ind2,n))*ff(m+ind1,n,k) - lcu(m+ind2,n)*fx0(m+ind2,n,k)
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=1,nz
        ind1= nint( 0.5 - sign(0.5,vh(m,n,k)) )
        ind2= nint(  - sign(1.0,vh(m,n,k)) )
        fy(m,n,k) = (1.0 + lcv(m,n+ind2))*ff(m,n+ind1,k) - lcv(m,n+ind2)*fy0(m,n+ind2,k)
       enddo
      endif

      if (lu(m,n)>0.5) then
       do k=2,nz
        ind1= nint( - 0.5 - sign(0.5,ww(m,n,k)) )
        ind2= nint(  - sign(1.0,ww(m,n,k)) )
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

!$omp parallel do private(m,n,k,hmid,rhs,flxz) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lu(m,n)>0.5) then
       hmid=(hhq(m,n)+hhq_n(m,n))/2.0
       
       flxz(1)=flx_top(m,n)
       do k=2,nz
       flxz(k) = ww(m,n,k)*fz(m,n,k)
       enddo
       flxz(nz+1)=0.0
       
       do k=1,nz
        rhs=(uh(m,n,k)*dyh(m,n)*fx(m,n,k) - uh(m-1,n,k)*dyh(m-1,n)*fx(m-1,n,k)            &
            +vh(m,n,k)*dxh(m,n)*fy(m,n,k) - vh(m,n-1,k)*dxh(m,n-1)*fy(m,n-1,k))/sqt(m,n)  &
            +(flxz(k+1) - flxz(k))/dz(k)
        ff(m,n,k) = ff(m,n,k)*(hmid/hhq_n(m,n)) - rhs*tau/2.0/hhq_n(m,n)
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

endsubroutine tracer_transport_cab_mod

!==========================================================================================================
subroutine tracer_sdiff(ff,fx,fy,tau,mu,factor_mu)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in):: mu
real(4), intent(in):: tau, factor_mu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: fx,fy

integer m,n,k
real(8) rhs

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
     fx(m,n,:) = 0.0
     fy(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!$omp parallel do private(m,n,k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lcu(m,n)>0.5) then
      do k=1,nz      
       fx(m,n,k) = - (mu(m,n,k)+mu(m+1,n,k))/2.0*factor_mu     &
                *hhu_n(m  ,n)*dyh(m  ,n)/dxt(m  ,n)*(ff(m+1,n,k)-ff(m  ,n,k))
      enddo
     endif
     
     if(lcv(m,n)>0.5) then
      do k=1,nz 
       fy(m,n,k) = - (mu(m,n,k)+mu(m,n+1,k))/2.0*factor_mu     &
                *hhv_n(m,n  )*dxh(m,n  )/dyt(m,n  )*(ff(m,n+1,k)-ff(m,n  ,k))
      enddo     
     endif

    enddo
   enddo
!$omp end parallel do        

    call syncborder_real(fx,nz)
    call syncborder_real(fy,nz)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz,nnn,nn)            
    end if

!$omp parallel do private(m,n,k,rhs) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lu(m,n)>0.5) then
      do k=1,nz
       rhs=fx(m,n,k)-fx(m-1,n,k)+fy(m,n,k)-fy(m,n-1,k)
       ff(m,n,k) = ff(m,n,k) - tau*rhs/hhq_n(m,n)/sqt(m,n)   
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

endsubroutine tracer_sdiff

!====================================================================
subroutine tracer_zdiff(ff,tau,mu,factor_mu,z3d,rhs_cc,rhs_mc,rhs_pc,rhs_cm,rhs_cp)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in):: mu,z3d
real(4), intent(in):: tau, factor_mu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: rhs_cc,rhs_mc,rhs_pc,rhs_cm,rhs_cp

integer m,n,k, m1, n1, ktop, kbot, kmean, kk
real(8) rhs, depth, rate, fcenter, fside, coef, flx, wght_top, wght_bot
   
!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      rhs_cc(m,n,:) = 0.0
      rhs_cm(m,n,:) = 0.0
      rhs_cp(m,n,:) = 0.0
      rhs_mc(m,n,:) = 0.0
      rhs_pc(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!$omp parallel do private(m, n, k, m1, n1, ktop, kbot, kmean, kk, depth, rate, fcenter, fside,        &
!$omp                     coef,flx, wght_top, wght_bot) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end   
     if(lu(m,n)>0.5) then
      do k=1,nz

       depth=z3d(m,n,k)
       fcenter=ff(m,n,k)
       coef=mu(m,n,k)/2.0*factor_mu

!point m+1       
       m1=m+1
       if(lu(m1,n)>0.5.and.depth>=z3d(m1,n,1).and.depth<=z3d(m1,n,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(depth>z3d(m1,n,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(depth-z3d(m1,n,kk))/(z3d(m1,n,kk+1)-z3d(m1,n,kk))
        wght_top=1.0-rate
        wght_bot=rate
        fside=ff(m1,n,kk)*wght_top+ff(m1,n,kk+1)*wght_bot
        flx=coef*(fside-fcenter)/dxt(m,n)*dy(m,n)*hhq_n(m,n)*dz(k)
        
        rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
        rhs_pc(m,n,kk  )=rhs_pc(m,n,kk  ) - wght_top*flx
        rhs_pc(m,n,kk+1)=rhs_pc(m,n,kk+1) - wght_bot*flx
       endif
!end point m+1  

!point m-1       
       m1=m-1
       if(lu(m1,n)>0.5.and.depth>=z3d(m1,n,1).and.depth<=z3d(m1,n,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(depth>z3d(m1,n,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(depth-z3d(m1,n,kk))/(z3d(m1,n,kk+1)-z3d(m1,n,kk))
        wght_top=1.0-rate
        wght_bot=rate
        fside=ff(m1,n,kk)*wght_top+ff(m1,n,kk+1)*wght_bot
        flx=coef*(fside-fcenter)/dxt(m1,n)*dy(m,n)*hhq_n(m,n)*dz(k)

        rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
        rhs_mc(m,n,kk  )=rhs_mc(m,n,kk  ) - wght_top*flx
        rhs_mc(m,n,kk+1)=rhs_mc(m,n,kk+1) - wght_bot*flx
       endif
!end point m-1  

!point n+1       
       n1=n+1
       if(lu(m,n1)>0.5.and.depth>=z3d(m,n1,1).and.depth<=z3d(m,n1,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(depth>z3d(m,n1,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(depth-z3d(m,n1,kk))/(z3d(m,n1,kk+1)-z3d(m,n1,kk))
        wght_top=1.0-rate
        wght_bot=rate
        fside=ff(m,n1,kk)*wght_top+ff(m,n1,kk+1)*wght_bot
        flx=coef*(fside-fcenter)/dyt(m,n)*dx(m,n)*hhq_n(m,n)*dz(k)

        rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
        rhs_cp(m,n,kk  )=rhs_cp(m,n,kk  ) - wght_top*flx
        rhs_cp(m,n,kk+1)=rhs_cp(m,n,kk+1) - wght_bot*flx
       endif
!end point n+1     

!point n-1       
       n1=n-1
       if(lu(m,n1)>0.5.and.depth>=z3d(m,n1,1).and.depth<=z3d(m,n1,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(depth>z3d(m,n1,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(depth-z3d(m,n1,kk))/(z3d(m,n1,kk+1)-z3d(m,n1,kk))
        wght_top=1.0-rate
        wght_bot=rate
        fside=ff(m,n1,kk)*wght_top+ff(m,n1,kk+1)*wght_bot
        flx=coef*(fside-fcenter)/dyt(m,n1)*dx(m,n)*hhq_n(m,n)*dz(k)

        rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
        rhs_cm(m,n,kk  )=rhs_cm(m,n,kk  ) - wght_top*flx
        rhs_cm(m,n,kk+1)=rhs_cm(m,n,kk+1) - wght_bot*flx
       endif
!end point n-1 
      enddo
     endif
    enddo
   enddo
!$omp end parallel do    

    call syncborder_real(rhs_pc, nz)
    call syncborder_real(rhs_mc, nz)
    call syncborder_real(rhs_cp, nz)
    call syncborder_real(rhs_cm, nz)
                
    if(periodicity_x/=0) then
      call cyclize_x(rhs_pc,nx,ny,nz,mmm,mm)
      call cyclize_x(rhs_mc,nx,ny,nz,mmm,mm)  
    end if
    
    if(periodicity_y/=0) then
      call cyclize_y(rhs_cp,nx,ny,nz,nnn,nn)
      call cyclize_y(rhs_cm,nx,ny,nz,nnn,nn)
    end if

!$omp parallel do private(m, n, k, rhs) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lu(m,n)>0.5) then
      do k=1,nz
       rhs= ( rhs_cc(m,n  ,k) + rhs_mc(m+1,n  ,k) + rhs_pc(m-1,n  ,k)    &
                              + rhs_cm(m  ,n+1,k) + rhs_cp(m  ,n-1,k) )/(sqt(m,n)*dz(k)) 
       ff(m,n,k)=ff(m,n,k) + tau*rhs/hhq_n(m,n)
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

endsubroutine tracer_zdiff

!====================================================================
subroutine tracer_isopyc(ff,tau,mu,factor_mu,den1,z3d,rhs_cc,rhs_mc,rhs_pc,rhs_cm,rhs_cp)
 
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in):: mu,z3d,den1
real(4), intent(in):: tau, factor_mu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz):: rhs_cc,rhs_mc,rhs_pc,rhs_cm,rhs_cp

integer m,n,k, m1, n1, ktop, kbot, kmean, kk

real(8) rhs, dnst, rate, fcenter, fside, coef, flx, wght_top, wght_bot, angle

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      rhs_cc(m,n,:) = 0.0
      rhs_mc(m,n,:) = 0.0
      rhs_pc(m,n,:) = 0.0
      rhs_cm(m,n,:) = 0.0
      rhs_cp(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do   

!$omp parallel do private(m, n, k, m1, n1, ktop, kbot, kmean, kk, dnst, rate, fcenter, fside,        &
!$omp                     coef,flx, wght_top, wght_bot,angle) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end   
     if(lu(m,n)>0.5) then
      do k=1,nz

       dnst=den1(m,n,k)
       fcenter=ff(m,n,k)
       coef=mu(m,n,k)/2.0*factor_mu

!point m+1       
       m1=m+1
       if(lu(m1,n)>0.5.and.dnst>=den1(m1,n,1).and.dnst<=den1(m1,n,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(dnst>den1(m1,n,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(dnst-den1(m1,n,kk))/max(den1(m1,n,kk+1)-den1(m1,n,kk),epsrho)
        wght_top=1.0-rate
        wght_bot=rate
        angle=(z3d(m1,n,kk)*wght_top+z3d(m1,n,kk+1)*wght_bot - z3d(m,n,k))/dxt(m,n)
        if(abs(angle)<=angle_max) then
         fside=ff(m1,n,kk)*wght_top+ff(m1,n,kk+1)*wght_bot         
         flx=coef*(fside-fcenter)/dxt(m,n)*dy(m,n)*hhq_n(m,n)*dz(k)

         rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
         rhs_pc(m,n,kk  )=rhs_pc(m,n,kk  ) - wght_top*flx
         rhs_pc(m,n,kk+1)=rhs_pc(m,n,kk+1) - wght_bot*flx
        endif
       endif
!end point m+1  

!point m-1       
       m1=m-1
       if(lu(m1,n)>0.5.and.dnst>=den1(m1,n,1).and.dnst<=den1(m1,n,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(dnst>den1(m1,n,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(dnst-den1(m1,n,kk))/max(den1(m1,n,kk+1)-den1(m1,n,kk),epsrho)
        wght_top=1.0-rate
        wght_bot=rate
        angle=(z3d(m1,n,kk)*wght_top+z3d(m1,n,kk+1)*wght_bot - z3d(m,n,k))/dxt(m1,n)
        if(abs(angle)<=angle_max) then        
         fside=ff(m1,n,kk)*wght_top+ff(m1,n,kk+1)*wght_bot
         flx=coef*(fside-fcenter)/dxt(m1,n)*dy(m,n)*hhq_n(m,n)*dz(k)

         rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
         rhs_mc(m,n,kk  )=rhs_mc(m,n,kk  ) - wght_top*flx
         rhs_mc(m,n,kk+1)=rhs_mc(m,n,kk+1) - wght_bot*flx
        endif
       endif
!end point m-1  

!point n+1       
       n1=n+1
       if(lu(m,n1)>0.5.and.dnst>=den1(m,n1,1).and.dnst<=den1(m,n1,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(dnst>den1(m,n1,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(dnst-den1(m,n1,kk))/max(den1(m,n1,kk+1)-den1(m,n1,kk),epsrho)
        wght_top=1.0-rate
        wght_bot=rate
        angle=(z3d(m,n1,kk)*wght_top+z3d(m,n1,kk+1)*wght_bot - z3d(m,n,k))/dyt(m,n)
        if(abs(angle)<=angle_max) then 
         fside=ff(m,n1,kk)*wght_top+ff(m,n1,kk+1)*wght_bot
         flx=coef*(fside-fcenter)/dyt(m,n)*dx(m,n)*hhq_n(m,n)*dz(k)

         rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
         rhs_cp(m,n,kk  )=rhs_cp(m,n,kk  ) - wght_top*flx
         rhs_cp(m,n,kk+1)=rhs_cp(m,n,kk+1) - wght_bot*flx
        endif
       endif
!end point n+1     

!point n-1       
       n1=n-1
       if(lu(m,n1)>0.5.and.dnst>=den1(m,n1,1).and.dnst<=den1(m,n1,nz)) then
        ktop=1
        kbot=nz

        do while(kbot>ktop+1)
         kmean=(ktop+kbot)/2
         if(dnst>den1(m,n1,kmean)) then
          ktop=kmean
         else
          kbot=kmean
         endif
        enddo

        kk=ktop
        rate=(dnst-den1(m,n1,kk))/max(den1(m,n1,kk+1)-den1(m,n1,kk),epsrho)
        wght_top=1.0-rate
        wght_bot=rate
        angle=(z3d(m,n1,kk)*wght_top+z3d(m,n1,kk+1)*wght_bot - z3d(m,n,k))/dyt(m,n1)
        if(abs(angle)<=angle_max) then          
         fside=ff(m,n1,kk)*wght_top+ff(m,n1,kk+1)*wght_bot
         flx=coef*(fside-fcenter)/dyt(m,n1)*dx(m,n)*hhq_n(m,n)*dz(k)

         rhs_cc(m,n,k)=rhs_cc(m,n,k) + flx
         rhs_cm(m,n,kk  )=rhs_cm(m,n,kk  ) - wght_top*flx
         rhs_cm(m,n,kk+1)=rhs_cm(m,n,kk+1) - wght_bot*flx
        endif
       endif
!end point n-1 
      enddo
     endif
    enddo
   enddo
!$omp end parallel do    

    call syncborder_real(rhs_pc, nz)
    call syncborder_real(rhs_mc, nz)
    call syncborder_real(rhs_cp, nz)
    call syncborder_real(rhs_cm, nz)
           
    if(periodicity_x/=0) then
      call cyclize_x(rhs_pc,nx,ny,nz,mmm,mm)
      call cyclize_x(rhs_mc,nx,ny,nz,mmm,mm)
    end if
            
    if(periodicity_y/=0) then
      call cyclize_y(rhs_cp,nx,ny,nz,nnn,nn)
      call cyclize_y(rhs_cm,nx,ny,nz,nnn,nn)
    end if

!$omp parallel do private(m,n,k,rhs) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lu(m,n)>0.5) then
      do k=1,nz
       rhs= ( rhs_cc(m,n,k) + rhs_mc(m+1,n,k) + rhs_pc(m-1,n,k)      &
                            + rhs_cm(m,n+1,k) + rhs_cp(m,n-1,k) )/(sqt(m,n)*dz(k))
       ff(m,n,k)=ff(m,n,k) + tau*rhs/hhq_n(m,n)
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

endsubroutine tracer_isopyc

!====================================================================
subroutine tracer_vdiff(ff,tau,nu,factor_nu,flx_top,ig_top)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: nu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: flx_top
integer, dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in)::  ig_top
real(4), intent(in):: tau, factor_nu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: ff

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz) 
real(8) bp, dp, dm

!$omp parallel do private(m,n,k,dp,dm,a,b,c,eta,rksi)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
         bp = hhq_n(m,n)/tau

!Upper layer          
          k=1

          dm = nu(m,n,k  )*factor_nu/hhq_n(m,n)/hzt(k  )
          dp = nu(m,n,k+1)*factor_nu/hhq_n(m,n)/hzt(k+1)

          c(k) = -dp
          a(k) = 0.0

          if(ig_top(m,n)==1) then
           b(k) =  bp*dz(k) + dp + dm
           eta(k) = bp*ff(m,n,k)*dz(k) + dm*flx_top(m,n)
          elseif(ig_top(m,n)==2) then
           b(k) =  bp*dz(k) + dp 
           eta(k) = bp*ff(m,n,k)*dz(k) + flx_top(m,n)
          endif

! internal points.
         do k=2,nz-1

          dm = dp
          dp = nu(m,n,k+1)*factor_nu/hhq_n(m,n)/hzt(k+1)

          c(k) = -dp
          a(k) = -dm
          b(k) =  bp*dz(k) + dp + dm
          eta(k) = bp*ff(m,n,k)*dz(k)
         enddo

!Bottom layer
          k=nz

          dm = dp
          dp = 0.0

          c(k) = 0.0
          a(k) = -dm
          b(k) =  bp*dz(k) + dm
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

endsubroutine tracer_vdiff

endmodule ts_routes