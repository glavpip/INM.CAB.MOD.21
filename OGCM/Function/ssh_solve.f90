!================================================
! MODULES FOR SEA LEVEL CALCULATIONS PREPERED BY Rusakov. 
! No splitting here. We form full matrice and solve it
!======================================================================
      module sea_level_no_split
      use constants
      use basin_grid
      implicit none

!======================================================================
!      DECLARATIONS
!======================================================================

        !
        !   Matrice in CSR format
        !
        DOUBLE PRECISION, ALLOCATABLE::A(:)
        INTEGER, ALLOCATABLE::IA(:), JA(:)
!        DOUBLE PRECISION A(200000)
!        INTEGER IA(200000), JA(200000)      
        !
        !   Matrice of Precondition in CSR format
        !
        DOUBLE PRECISION, ALLOCATABLE::AU(:)
        INTEGER, ALLOCATABLE::IAU(:), JAU(:)
        INTEGER, ALLOCATABLE::PERMUTATION(:) !vector for reordering for precond
        
        INCLUDE 'ssh.fi'  !input parameters of solver and number for precondition 

      !
        !   RHS
      !
         DOUBLE PRECISION, ALLOCATABLE::guess_init(:) , RHS(:) 
!        DOUBLE PRECISION RHS(200000)

        INTEGER nnz,  &  ! number of nonzeroes
             maxnnz,  &  ! size of A and JA
             nvar        ! number of variables

        !
        !     Arrays of indexes. connect two dimensional problem with one dimensional 
        !      vector excluding the ground.   
        !
        INTEGER  UINDX(NX,NY), VINDX(NX,NY), SLHINDX(NX,NY)
      INTEGER  IS_FORMED, MATRICE_PERMUTED,  PRECONDITIONER_FORMED
      
      INTEGER, PARAMETER::DEBUG = 0, & !set debug to 1 to see some useful information
                          DO_REORDERING_FOR_SPARSITY = 1 

        INTEGER, SAVE :: factors(8), &! array of pointer for Super LU
                        isfactored  = 0

!======================================================================
!     END OF DECLARATIONS
!======================================================================

      contains      
!======================================================================
!     SUBROUTINES 
!======================================================================
      
         
!----------------------------------------------------------------------
!     Form matrice of sea level task
!
!    dU/dt + rb*U - l*V = mg d sl/dx
!
!    dU/dt + rb*V + l*U = ng d sl/dx
!
!             SL          d            d   n                      
! freesurface*-----  - [  --(HU) +    --(  - * HV ) ] = 0          
!              dt         dx           dy  m                      
!
!    Implicit time scheme:
!
!           U-Uo                             d SLH               
!           ----  + Rb * U - l * V = m*g*H * -----               
!           tau                              d x                 
!                                                                
!                                                                
!           V-Vo                             d SLH               
!           ----  + Rb * V + l * U = n*g*H * -----     (1)       
!           tau                              d y                 
!                                                                
!             SLH - SLHo       d          d   n                      
! freesurface*----------  - [  --(U) +    --( - * V ) ] = 0          
!              m*tau           dx         dy  m                      
!                                                                
!********************************************************************
!      remark:   freesurface =0 - rigid lid condition
!                freesurface =1 - long gravity waves are avaliable
!----------------------------------------------------------------------
!
!     Form matrice in compress sparse row format (CSR). 
!    Use common ordering. the first is u, the second is v , the third is slh. 
!    index m -changed first. 
!
!  side effect 
!     uncyclize operation on LCU and LCV
!
!----------------------------------------------------------------------
      SUBROUTINE FORM_MATRICE_SSH(RBOTTOM, tau)
!      IMPLICIT NONE
      INTEGER M,N,K, k1, k2, k3
      real   tau
      real*4 RBOTTOM(NX,NY)
      real*8 frictau
      INTEGER ileft, iright, jbot, jtop


!----------------------------------------------------------------------
!   
!      Estimate matrix size, Form numbering of variables.
!      
!----------------------------------------------------------------------
      PRECONDITIONER_FORMED = 0
      IS_FORMED = 0
      MATRICE_PERMUTED = 0
      
!
!     construct indexes
!
      UINDX   = 0
    VINDX   = 0
    SLHINDX = 0

      k  = 0 ! index of the variable
    k1 = 0
    k2 = 0
    k3 = 0
!    
!     U variables
! 
      do n = NNN, NN
          do m = MMM, MM
              if ( LCU(m,n).gt.0.5 ) then
                  k  = k + 1
                k1 = k1 + 1
                  UINDX(m,n) = k
            endif
        enddo
    enddo
          
      
!    
!     V variables
!
      do n = NNN, NN
          do m = MMM, MM
              if ( LCV(m,n).gt.0.5 ) then
                  k  = k + 1
                k2 = k2 + 1
                  VINDX(m,n) = k
            endif
        enddo
    enddo


!    
!     SLH variables
!
      do n = NNN, NN
          do m = MMM, MM
              if ( LU(m,n).gt.0.5 ) then
                  k  = k + 1
                k3 = k3 + 1
                  SLHINDX(m,n) = k
            endif
        enddo
    enddo

!    
!     estimate the number of nonzeroes in a matrice
!
      maxnnz = k1 * 7 + k2 * 7 + k3 * 5
!    
!     number of variables
!
    nvar = k
      if (debug.gt.0) then
        write (*,*) ' in form matrix number of variables = ', nvar
        write (*,*) ' maximum number of nonzeroes is '
          write (*,*) ' 7*nu_points+7*nv_points+5*nT_points', maxnnz
    endif
!    
!     allocate memory for matrice rhs and initial guess 
!
      ALLOCATE (A(maxnnz), JA(maxnnz), IA(nvar+1), RHS(nvar+1))
      ALLOCATE (guess_init(nvar+1))

!----------------------------------------------------------------------
!   
!      Form matrice
!      
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!   
!      First equation
!      
!----------------------------------------------------------------------
!   something like: see reports
!
!           U-Uo                             d SLH               
!           ----  + Rb * U - l * V = m*g*H * -----               
!           tau                              d x                 
!
!   coriolis on LUH grid
!----------------------------------------------------------------------

      nnz = 0
      do n = NNN, NN
        do m = MMM, MM
        
          ileft  = m-1
          iright = m+1

          if( periodicity_x>0 ) then
              if (m.eq.MMM) then 
                  ileft  = MM
              endif
              if (m.eq.MM) then   !cyclic right boundary
                  iright = MMM
              end if
          end if !if cyclic condition

          jbot = n-1
          jtop = n+1

          if( periodicity_y>0 ) then
              if (n.eq.nnn) then 
                  jbot = nn
              endif
              if (n.eq.nn) then   !cyclic right boundary
                  jtop = nnn
              end if
          end if !if cyclic condition

          if ( LCU(m,n).gt.0.5 ) then
            !diagonal
          nnz = nnz + 1
            frictau = 1.0/tau + (rbottom(m,n)+rbottom(m+1,n))/2.0 
            A(nnz)  = frictau
          JA(nnz) = UINDX(m,n)
            IA(UINDX(m,n)) = nnz
              
            !v elements
          if ( LCV(m,n).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = -rlh_s(m,n)/4.0
            JA(nnz) = VINDX(m,n)
          endif
            if ( LCV(iright,n).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = -rlh_s(m,n)/4.0 
            JA(nnz) = VINDX(iright,n)
          endif
          if ( LCV(m,jbot).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = -rlh_s(m,jbot)/4.0
            JA(nnz) = VINDX(m,jbot)
          endif
          if (LCV(iright,jbot).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = -rlh_s(m,jbot)/4.0
            JA(nnz) = VINDX(iright,jbot)
          endif

          !slh elements
          if ( LU(m,n).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  =-FreeFallAcc*hhu_rest(m,n)/DXT(M,N)
            JA(nnz) = SLHINDX(m,n)
          endif
          if ( LU(iright,n).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = FreeFallAcc*hhu_rest(m,n)/DXT(M,N)
            JA(nnz) = SLHINDX(iright,n)
          endif
          endif


      enddo
    enddo

!----------------------------------------------------------------------
!   
!      Second equation
!      
!----------------------------------------------------------------------
!   something like: see reports
!
!           V-Vo                             d SLH               
!           ----  + Rb * V + l * U = n*g*H * -----     (1)       
!           tau                              d y                 
!
!   coriolis on LUH grid
!----------------------------------------------------------------------
      do n = NNN, NN
        do m = MMM, MM

          ileft  = m-1
          iright = m+1

          if( periodicity_x>0 ) then
              if (m.eq.MMM) then 
                  ileft  = MM
              endif
              if (m.eq.MM) then   !cyclic right boundary
                  iright = MMM
              end if
          end if !if cyclic condition

          jbot = n-1
          jtop = n+1

          if( periodicity_y>0 ) then
              if (n.eq.nnn) then 
                  jbot = nn
              endif
              if (n.eq.nn) then   !cyclic right boundary
                  jtop = nnn
              end if
          end if !if cyclic condition

          if ( LCV(m,n).gt.0.5 ) then
            !diagonal
          nnz = nnz + 1
            frictau = 1.0/tau + (rbottom(m,n)+rbottom(m,n+1))/2.0
            A(nnz)  = frictau
          JA(nnz) = VINDX(m,n)

            IA(VINDX(m,n)) = nnz
            
            !u elements
          if ( LCU(m,n).gt.0.5 ) then
            nnz = nnz + 1
              A(nnz)  =rlh_s(m,n)/4.0
            JA(nnz) = UINDX(m,n)
          endif
            if ( LCU(ileft,n).gt.0.5 ) then
              nnz = nnz + 1
              A(nnz)  =rlh_s(ileft,n)/4.0
              JA(nnz) = UINDX(ileft,n)
          endif
          if ( LCU(m,jtop).gt.0.5 ) then
            nnz = nnz + 1
              A(nnz)  =rlh_s(m,n)/4.0
            JA(nnz) = UINDX(m,jtop)
          endif
            if ( LCU(ileft,jtop).gt.0.5 ) then
              nnz = nnz + 1
              A(nnz)  =rlh_s(ileft,n)/4.0
            JA(nnz) = UINDX(ileft,jtop)
            end if              

          !slh elements
          if ( LU(m,n).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  =-FreeFallAcc*hhv_rest(m,n)/DYT(M,N)
            JA(nnz) = SLHINDX(m,n)
          endif
          if ( LU(m,jtop).gt.0.5 ) then
            nnz = nnz + 1
            A(nnz)  = FreeFallAcc*hhv_rest(m,n)/DYT(M,N)
            JA(nnz) = SLHINDX(m,jtop)
          endif
          endif

      enddo
    enddo

!----------------------------------------------------------------------
!   
!      Third equation
!      
!----------------------------------------------------------------------
!   something like: see reports
!
!                                                                
!             SLH - SLHo                d           d    H                    
! freesurface*----------  - m(ij)n ([  --(HU/n) +   --(  - * V ) ] = 0          
!                  tau                 dx           dy   m                      
! 
!   coriolis on LUH grid
!----------------------------------------------------------------------

      do n = NNN, NN
        do m = MMM, MM

          ileft  = m-1
          iright = m+1

          if( periodicity_x>0 ) then
              if (m.eq.MMM) then 
                  ileft  = MM
              endif
              if (m.eq.MM) then   !cyclic right boundary
                  iright = MMM
              end if
          end if !if cyclic condition

          jbot = n-1
          jtop = n+1

          if( periodicity_y>0 ) then
              if (n.eq.nnn) then 
                  jbot = nn
              endif
              if (n.eq.nn) then   !cyclic right boundary
                  jtop = nnn
              end if
          end if !if cyclic condition

          if ( LU(m,n).gt.0.5 ) then
            
            !diagonal
          nnz = nnz + 1
            A(nnz)  = 1.0/TAU !the same have to be in RHS
          JA(nnz) = SLHINDX(m,n)
            IA(SLHINDX(m,n)) = nnz

            !u elements NAN
          if ( LCU(m,n).gt.0.5 ) then
              nnz = nnz + 1
            A(nnz)  = DYH(M,N)/sqt(M,N)
            JA(nnz) = UINDX(m,n)
          endif
          if ( LCU(ileft,n).gt.0.5 ) then
              nnz = nnz + 1
            A(nnz)  =-DYH(ileft,N)/sqt(M,N)
            JA(nnz) = UINDX(ileft,n)
          endif

            !v elements
          if ( LCV(m,n).gt.0.5 ) then
              nnz = nnz + 1
            A(nnz)  = DXH(M,N)/sqt(M,N)
            JA(nnz) = VINDX(m,n)
          endif
          if ( LCV(m,jbot).gt.0.5 ) then
              nnz = nnz + 1
            A(nnz)  =-DXH(M,jbot)/sqt(M,N)
            JA(nnz) = VINDX(m,jbot)
          endif

          endif
      enddo
      
      
    enddo
      IA(nvar+1) = nnz+1

      IS_FORMED = 1

      if (debug.gt.0) then
           write(*,*) 'Actual number on non-zeroes ',nnz
      end if 

!c
!c    test. sort columns in a matrice
!c
!      ALLOCATE (jau(nzmax))
!      call csort (nvar,a,ja,ia,jau,.true.) 
!      DEALLOCATE (jau)

!c
!c       Test dump
!c
!      OPEN (77,FILE='matrice.dat')
!     call dump (1,nvar,.true.,a,ja,ia,77)
      

!      write(77,*) 'A: nnz = ', nnz
!      do n = 1, nnz
!        write(77,*) 'n = ', n, ' A(n)', A(n)
!      enddo
!      write(77,*) '***************************************************'
!      write(77,*) 'JA: nnz = ', nnz
!      do n = 1, nnz
!        write(77,*) 'n = ', n, ' JA(n)', JA(n)
!      enddo
!      write(77,*) '***************************************************'
!      write(77,*) 'IA: nvar = ', nvar
!      do n = 1, nvar+1
!        write(77,*) 'n = ', n, ' IA(n)', IA(n)

!c
!c       dump in ijv (matlab format)
!c
!      write(77,*) 'IA: nvar = ', nvar
!      do n = 1, nvar+1
!          do k = IA(n), IA(n+1)-1
!              write(77,*) n, JA(k), A(k)
!          enddo       
!      enddo
!      CLOSE(77)
!      stop
      if (debug.gt.0) then
          do k = 1, nnz
              if (JA(k).le.0) then 
                   write(*,*) 'error in form matrice ', 'zero index',k, JA(k)
              end if 
              if (A(k).ne.A(k)) then 
                   write(*,*) 'error in form matrice ', 'NAN',k, JA(k)
              end if
      

          enddo       
      end if 
      write (*,*) ' periodicity_x = ', periodicity_x
      write (*,*) ' periodicity_y = ', periodicity_y
      END SUBROUTINE


!----------------------------------------------------------------------
!     SOLVE matrice of sea level task with iterative method based on krylov spaces
!     or direct method.
!     first try to use direct method 
!     if error then use preconditioned iterative method
!----------------------------------------------------------------------
!      SUBROUTINE SOLVE_MATRICE(U,V,SLH)
!      IMPLICIT NONE
!      
!      call SOLVE_MATRICE_SUPERLU(U,V,SLH, ierr)
!
!      if (ierr.ne.0) then
!
!      endif
!
!      END SUBROUTINE
!--------------------------------------------------------------
!
!     Form ILUT preconditioner for iterative method
!      Trying to find the stronger preconditioner for avalaible memory (set by nzmax)
!     if DO_REORDERING_FOR_SPARSITY then use amd befor ilut
!
!--------------------------------------------------------------
      SUBROUTINE FORM_PRECONDITION()
      integer ierr,lfil,nwk
    integer, allocatable :: iw(:) ! (2*nvar+1) temporary array for ilut
      real*8   tol
    real*8, allocatable  :: wk(:) ! (nvar+2)   temporary array for ilut
      REAL TIME_IN_SEC,TIME_IN_SEC0,TIME_IN_SEC1
    integer niter
    integer, parameter::max_iter=20
!     
!     data for preconditioner
!           
      if(DO_REORDERING_FOR_SPARSITY.gt.0) then
          call do_amd_reorder_of_matrice()    
      end if

      ierr = -1

      if (precond.gt.0) then

        ALLOCATE (au(nzmax), jau(nzmax),iau(nvar+1) )
      ALLOCATE (iw(2*nvar+1), wk(nvar+2))     

!
!         initial preconditioner parameters
!
        lfil = nvar
        tol = 1E-16
        niter = 0
        do while (ierr.ne.0.and.niter.lt.max_iter)        
          niter = niter + 1
          nwk = nzmax
          if (debug.gt.0) then
              write (*,*) ' *******************************************'
            write (*,*) ' Building of the precondition ' 
              write (*,*) ' lfil = ', lfil,' tol = ', tol
          endif         
          call cpu_time(TIME_IN_SEC0)
          call ilut (nvar,a,ja,ia,lfil,tol,au,jau,iau,nwk,wk,iw,ierr)

          call cpu_time(TIME_IN_SEC1) 
          TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
          write (*,*) ' Build ierr=', ierr
          write (*,*) 'End of Building of the precondition '
          write (*,*) ' Time on precond  ', TIME_IN_SEC
        write (*,*) 'Actual Number of non-zeroes if preconditioner  ', iau(nvar-1)
          write (*,*) ' ******************************************* '
      
!
!          if unsuccess then try  to use weaker preconditionier
!
          lfil = lfil / 2
          tol = tol * 7
          if(lfil.lt.50) lfil = 50
        if(tol.gt.1) tol = 0.5

        enddo !(while (ierr.lt.0)

      if (ierr.ne.0) then 
        DEALLOCATE (iau, jau, au)
        ALLOCATE (au(nnz), jau(nnz),iau(nvar+1) )

!         call ilu0 (nvar,a,ja,ia, au,jau,iau, iw,ierr)
           lfil = 0
          tol  = 0.0
         nwk = nnz
           call ilut (nvar,a,ja,ia,lfil,tol,au,jau,iau,nwk,wk,iw,ierr)


         write (*,*) 'End of Building of the precondition ILU0 '
         write (*,*) ' Build ierr=', ierr
         write (*,*) ' ******************************************* '
      end if
      
      DEALLOCATE(wk, iw)
      endif !if precond
        
    if (ierr.eq.0) then 
        PRECONDITIONER_FORMED = 1
    end if
    
      
      END SUBROUTINE

!----------------------------------------------------------------------
!     Form RHS
!     f1 - on LCU grid corresponds to equation for U
!     f2 - on LCV grid corresponds to equation for V
!     U0,V0,SLH0 are changed on exit ( 1/tau*[U0,V0,SLH0] are added)
!     RHS = [f1,f2,0] + 1/tau*[U0,V0, freesurface * SLH0]
!----------------------------------------------------------------------
      SUBROUTINE FORM_RHS(u0,v0,ssh0, hu, hu_r, hv, hv_r, f1, f2, vflux, tau)

      REAL*8   f1(NX,NY), f2(NX,NY)           
      REAL*8   u0(NX,NY), v0(NX,NY), ssh0(NX,NY),RHS_U(NX,NY),RHS_V(NX,NY),RHS_SSH(NX,NY)
      REAL(4)  tau, vflux(nx,ny), hu(nx,ny), hu_r(nx,ny), hv(nx,ny), hv_r(nx,ny)
      REAL*8   tau1
    INTEGER M,N
      tau1 = dble(1.0/tau)
      
!      RHS: 

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=NNN,NN
        DO M=MMM,MM
!   1-ST EQUATION
         IF(LCU(M,N).GT.0.5) THEN          
          U0(m,n)=U0(M,N)+tau*( f1(M,N)            &
          - dfloat(nonlinear_free_surface)*FreeFallAcc          &
          *(hu(m,n)-hu_r(m,n))*(ssh0(m+1,n)-ssh0(m,n))/dxt(m,n) )
        RHS_U(M,N)=U0(M,N)*tau1

         ELSE
          RHS_U(M,N)=0.0
         END IF

!   2-ND EQUATION
         IF(LCV(M,N).GT.0.5) THEN
          V0(m,n)=V0(m,n)+tau*(f2(M,N)             &
          - dfloat(nonlinear_free_surface)*FreeFallAcc          &
          *(hv(m,n)-hv_r(m,n))*(ssh0(m,n+1)-ssh0(m,n))/dyt(m,n) )

        RHS_V(M,N)=V0(M,N)*tau1
         ELSE
          RHS_V(M,N)=0.0
         END IF

      END DO
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(M,N)
      DO N=NNN,NN
        DO M=MMM,MM

!   3-RD EQUATION
         IF(LU(M,N).GT.0.5) THEN
        SSH0(M,N)=SSH0(M,N)+tau*vflux(m,n)*dfloat(variable_volume_budget)
          RHS_SSH(M,N)=SSH0(M,N) * tau1
         ELSE
          RHS_SSH(M,N)=0.0
         END IF

      END DO
    END DO
!$OMP END PARALLEL DO

!      put data from previous step to initial guess (in matrice order) 
!
      call reorderToMatriceOrder(U0, V0, SSH0, guess_init)
!     
!      put data to rhs convert all to matrice order
!     
      call reorderToMatriceOrder(RHS_U,RHS_V,RHS_SSH, RHS)

      END SUBROUTINE

!----------------------------------------------------------------------
!     SOLVE matrice of sea level task with iterative method based on krylov spaces
!----------------------------------------------------------------------
      SUBROUTINE SOLVE_MATRICE_ITER(U,V,SLH,NSLOR)
!      IMPLICIT NONE
      REAL*8   U(NX,NY),V(NX,NY),SLH(NX,NY) ! put  solution here
      integer  ipar(16),NSLOR
      real*8, allocatable :: vv(:), sol(:)
      integer maxits,iout
      real*8  eps,fpar(16)
      external gmres , bcg, tfqmr, fgmres, dqgmres, bcgstab
      integer, parameter:: nkrylov = 20
      integer  nwkspacesz 
      REAL TIME_IN_SEC,TIME_IN_SEC0,TIME_IN_SEC1

!
!       FORM_PRECONDITION has to be called before first FORM_RHS
!      if (PRECONDITIONER_FORMED.eq.0) THEN
!          CALL FORM_PRECONDITION() ! construct precondition , doesn't need if direct solver is used
!      end if 



      nwkspacesz = (nvar+3)*(nkrylov+2) + (nkrylov+1)*nkrylov/2 !for gmres
      nwkspacesz = nvar * 7  !for bcg
    nwkspacesz = nvar * 8  !for bcgstab
      nwkspacesz = nvar * 11 !for tqmr  

      ALLOCATE (vv(nwkspacesz), sol(nvar))


!     
!     data for GMRES
!     
      eps  = 1.0D-08
      maxits = 1000
      iout = 6          

!     ipar(2) -- status of the preconditioning:
!       0 == no preconditioning
!       1 == left preconditioning only
!       2 == right preconditioning only
      ipar(2) = precond * PRECONDITIONER_FORMED
      ipar(3) = 1
      ipar(4) = nwkspacesz
      ipar(5) = nkrylov
      ipar(6) = maxits

      fpar(1) = eps      !rtol
      fpar(2) = 2.22D-16 !atol
      
      if (debug.gt.0) then
        write (*,*) ' ******************************************* '
        write (*,*) ' Start matrix solution by iterative method'
        write (*,*) ' maximim iteration = ', maxits
        write (*,*) ' eps = ', eps
    endif


!c--------------------------------------------------------------
!c     call GMRES
!c--------------------------------------------------------------
!      call cpu_time(TIME_IN_SEC0)
!      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,
!     &                                       au,jau,iau,gmres)
!      call cpu_time(TIME_IN_SEC1) 
!      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0

!      write(*,*) ' Time tooked by GMRES ', TIME_IN_SEC
!      write(*,*) 'initial residual norm ', fpar(3)
!      write(*,*) 'current residual norm ', fpar(6)
!      write(*,*) 'target residual norm ', fpar(4)
!
!c--------------------------------------------------------------
!c     call bcg
!c--------------------------------------------------------------
!      call cpu_time(TIME_IN_SEC0)
!      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,
!     &                                       au,jau,iau,bcg)
!      call cpu_time(TIME_IN_SEC1) 
!      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!
!      write(*,*) ' Time tooked by bcg ', TIME_IN_SEC
!--------------------------------------------------------------
!c--------------------------------------------------------------
!c     call bcgstab
!c--------------------------------------------------------------
!      call cpu_time(TIME_IN_SEC0)
!      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,
!     &                                       au,jau,iau,bcgstab)
!      call cpu_time(TIME_IN_SEC1) 
!      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!
!      write(*,*) ' Time tooked by bcgstab ', TIME_IN_SEC
!--------------------------------------------------------------
!     call tfqmr
!--------------------------------------------------------------
      call cpu_time(TIME_IN_SEC0)
      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,au,jau,iau,tfqmr)
      call cpu_time(TIME_IN_SEC1) 
      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
      if (debug.gt.0) then
          write(*,*) ' Time tooked by tfqmr ', TIME_IN_SEC
      end if

!c--------------------------------------------------------------
!c     call fgmres
!c--------------------------------------------------------------
!      call cpu_time(TIME_IN_SEC0)
!      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,
!     &                                       au,jau,iau,fgmres)
!      call cpu_time(TIME_IN_SEC1) 
!      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!
!      write(*,*) ' Time tooked by fqmres ', TIME_IN_SEC
!
!c--------------------------------------------------------------
!c     call dqgmres
!--------------------------------------------------------------
!      call cpu_time(TIME_IN_SEC0)
!      call runrc(nvar,rhs,sol,ipar,fpar,vv,guess_init,A,JA,IA,
!     &                                       au,jau,iau,dqgmres)
!      call cpu_time(TIME_IN_SEC1) 
!      TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!
!      write(*,*) ' Time tooked by dqgmres ', TIME_IN_SEC
!
!--------------------------------------------------------------
!     put solution into original order
!--------------------------------------------------------------
      call reorderToOriginalOrder(U,V,SLH,sol)
      NSLOR=ipar(7)
      DEALLOCATE (sol,vv)

      END SUBROUTINE

!--------------------------------------------------------------
!--------------------------------------------------------------
      subroutine runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia, au,jau,iau,solver)
 !     implicit none
      integer n,ipar(16),ia(n+1),ja(*),iau(*),jau(*)
      real*8 fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*)
      external solver
!-----------------------------------------------------------------------
!   . It starts the iterative linear system solvers
!     with a initial guess supplied by the user.
!
!     The structure {au, jau, iau} is assumed to have the output from
!     the ILU* routines in ilut.f.
!
!-----------------------------------------------------------------------
!     local variables
!
      integer i, iou, its
      real*8 res, dnrm2
!     real dtime, dt(2), time
!     external dtime
      external dnrm2
      save its,res
!
!     ipar(2) can be 0, 1, 2, please don't use 3
!
      if (ipar(2).gt.2) then
         print *, 'I can not do both left and right preconditioning.'
         return
      endif
!
!     normal execution
!
      its = 0
      res = 0.0D0

!      sol = 1.0
!      call amux(n,sol,rhs,a,ja,ia)

!     set initial condition for iterative solver
      sol = guess
      if (debug.eq.1) then
          write (iou, *) '# the rhs norm is', dnrm2(n,rhs,1)
          write (iou, *) '# the sol norm is', dnrm2(n,sol,1)  
      endif



      iou = 6
      ipar(1) = 0
 10   call solver(n,rhs,sol,ipar,fpar,wk)
!
!     output the residuals
!
      if (ipar(7).ne.its) then     
         if (MOD(its,5).eq.0.and.debug.eq.1) THEN
              write (iou, *) its, real(res)
         endif
         its = ipar(7)         
      endif
      res = fpar(5)
!
      if (ipar(1).eq.1) then
         call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.2) then
         call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
         call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,iau)
         goto 10
      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,iau)
         goto 10
      else if (ipar(1).le.0) then
         if (ipar(1).eq.0) then
            if(debug.gt.0) then
              print *, 'Iterative sovler has satisfied convergence test'
            endif
         else if (ipar(1).eq.-1) then
            print *, 'Iterative solver has iterated too many times.'
         else if (ipar(1).eq.-2) then
            print *, 'Iterative solver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4), ' elements.'
            STOP 'gmres'
         else if (ipar(1).eq.-3) then
            print *, 'Iterative sovler is facing a break-down.'
         else
            print *, 'Iterative solver terminated. code =', ipar(1)
         endif
      endif
!     time = dtime(dt)
      if(debug.eq.1) then
          write (iou, *) ipar(7), real(fpar(6))
          write (iou, *) '# retrun code =', ipar(1),  '    convergence rate =', fpar(7)
      endif
      if (ipar(1).lt.0) then 
          write (iou, *) '# retrun code =', ipar(1)
          write (iou, *)  ' ipar(12) = ', ipar(12)
      end if

!
!     check the error
!
      if (debug.eq.1) then
          call amux(n,sol,wk,a,ja,ia)
          do i = 1, n
              wk(n+i) = rhs(i)
              wk(i) = wk(i) - rhs(i)
          enddo
          write (iou, *) '# the actual residual norm is', dnrm2(n,wk,1)
          write (iou, *) '# the rhs norm is', dnrm2(n,rhs,1)
      end if
!
      if (iou.ne.6) close(iou)
      return
      end subroutine


!--------------------------------------------------------------
!
!--------------------------------------------------------------
      subroutine reorderToMatriceOrder(U, V, SLH, OUT)
      REAL*8   U(NX,NY),V(NX,NY),SLH(NX,NY) 
      REAL*8   OUT(nvar)
      integer  m,n

!$OMP PARALLEL DO PRIVATE(M,N)
      do n = NNN, NN
          do m = MMM, MM
              
              if ( LCU(m,n).gt.0.5 ) then
                  OUT(UINDX(m,n)) = U(m,n)
            endif

              if ( LCV(m,n).gt.0.5 ) then
                  OUT(VINDX(m,n)) = V(m,n)
            endif

              if ( LU(m,n).gt.0.5 ) then
                  OUT(SLHINDX(m,n)) = SLH(m,n)
            endif

        enddo
    enddo
!$OMP END PARALLEL DO
      
      IF (MATRICE_PERMUTED.GT.0) THEN 
          call dvperm (NVAR, OUT, PERMUTATION)          
      ENDIF


      end subroutine
!--------------------------------------------------------------
!
!--------------------------------------------------------------
      subroutine reorderToOriginalOrder(U, V, SLH, IN)
      REAL*8    U(NX,NY),V(NX,NY),SLH(NX,NY) 
      REAL*8   IN(nvar)
      integer  m,n
      INTEGER,ALLOCATABLE :: inv_PERMUT(:) 

      IF (MATRICE_PERMUTED.GT.0) THEN 
        ALLOCATE  (inv_PERMUT(nvar))
          call invlpw(nvar, PERMUTATION, INV_PERMUT)
          call dvperm (NVAR, IN, inv_PERMUT)          
        DEALLOCATE  (inv_PERMUT)
!          call invert_dvperm (NVAR, IN, inv_PERMUT)          
      ENDIF

!$OMP PARALLEL DO PRIVATE(M,N)      
      do n = NNN, NN
         do m = MMM, MM
              
              if ( LCU(m,n).gt.0.5 ) then
                  U(m,n) = IN(UINDX(m,n))
            endif

              if ( LCV(m,n).gt.0.5 ) then
                  V(m,n) = IN(VINDX(m,n))
            endif

              if ( LU(m,n).gt.0.5 ) then
                  SLH(m,n) = IN(SLHINDX(m,n))
            endif

        enddo
    enddo
!$OMP END PARALLEL DO

      end subroutine

!--------------------------------------------------------------
!
!--------------------------------------------------------------
      subroutine test_reorder()
!      IMPLICIT NONE
      REAL*8   ff(nvar)
      REAL*8   U0(NX,NY), V0(NX,NY), S0(NX,NY) 
      integer m,n, err
       
      U0 = 1.0
      V0 = 2.0
       S0 = 3.0
      call reorderToMatriceOrder(U0, V0, S0, ff)

      U0 = -1.0
      V0 = -2.0
       S0 = -3.0
      call reorderToOriginalOrder(U0, V0, S0, ff)

      err = 0
      do m=1,nx
        do n=1,ny 
          if(U0(m,n).ne.1.0.and.LCU(m,n).gt.0.5) then 
              write(*,*) 'error: invalid u0 at m,n  ',m,n
          endif
          if(V0(m,n).ne.2.0.and.LCV(m,n).gt.0.5) then 
              write(*,*) 'error: invalid v0 at m,n  ',m,n
          endif
          if(S0(m,n).ne.3.0.and. LU(m,n).gt.0.5) then 
              write(*,*) 'error: invalid s0 at m,n  ',m,n
          endif
        enddo
      enddo
      if (err.gt.0) STOP 'test_reorder'

      write (*,*) 'test reorder simple ok'

      do m=1,nx
        do n=1,ny 
          U0(m,n) = m+5* n
          V0(m,n) = m+5* n*n
          S0(m,n) = m+5* n*n*n
        enddo
      enddo
      call reorderToMatriceOrder(U0, V0, S0, ff)

      U0 = 0
      V0 = 0
      S0 = 0

      call reorderToOriginalOrder(U0, V0, S0, ff)
      do m=1,nx
        do n=1,ny 
          if(U0(m,n).ne.(m+5* n).and.LCU(m,n).gt.0.5) then 
              write(*,*) 'error: invalid u0 at m,n  ',m,n
          endif
          if(V0(m,n).ne.(m+5* n*n).and.LCV(m,n).gt.0.5) then 
              write(*,*) 'error: invalid v0 at m,n  ',m,n
          endif
          if(S0(m,n).ne.(m+5* n*n*n).and. LU(m,n).gt.0.5) then 
              write(*,*) 'error: invalid s0 at m,n  ',m,n
          endif
        enddo
      enddo
      write (*,*) 'test reorder more complicated: ok'
      end subroutine
      
!      subroutine EnergyProduct(x1, x2)
!      REAL*8   x1(nvar), x2(nvar)
!      integer  m,n
!      DOUBLE PRECISION dot, dd
!   
!      dot = 0.0
!      do n = NNN, NN
!         do m = MMM, MM
!             if ( LCU(m,n).gt.0.5 ) then
!                  dot = x1(m,n)*x2(m,n)*(HHQ(m,n) + HHQ(m,n-1))/2
!            endif
!        enddo
!    enddo
!
!      do n = NNN, NN
!        do m = MMM, MM
!              if ( LCV(m,n).gt.0.5 ) then
!                  dot = x1(m,n)*x2(m,n)*(HHQ(m-1,n) + HHQ(m,n))/2
!            endif
!        enddo
!    enddo
!
!      do n = NNN, NN
!        do m = MMM, MM
!              if ( LU(m,n).gt.0.5 ) then
!                  dot = x1(m,n)*x2(m,n)
!            endif
!        enddo
!    enddo
!      end subroutine



!
!
!
      subroutine test_diagonal()
  !    IMPLICIT NONE
      real*8 diag(nvar)
      integer idiag(nvar)
      integer m, len 
      
      diag = 0.0 
      call getdia (nvar, nvar, 0,a,ja,ia,len,diag,idiag,0)
      
!      if (len.ne.nvar) THEN 
!          WRITE(*,*) 'there are zeroes on the main diagonal'
!          STOP 'test_is_diagonal'
!      endif

      do m=1,nvar
        if( diag(m) .ne. A( IA(m) )) THEN 
!          WRITE(*,*) 'incorrect diagonal m=', m
!          STOP 'test_is_diagonal'
        endif
        if( diag(m) .le. 0 ) THEN 
          WRITE(*,*) 'value <= zero on a diagonal m=', m
!          STOP 'test_is_diagonal'
        endif

      enddo

      WRITE(*,*) 'test diagonal: ok'
      end subroutine


!
!         check that all non-zeroes are place in proper place
!
      subroutine test_matrice()
!      IMPLICIT NONE
      DOUBLE PRECISION elm  
      double precision  getelm
      integer m, n
      integer i,j, iadd



      
      do m=1,nx
        do n=1,ny 
          !u row
          if(LCU(m,n).gt.0.5) then 

            i = UINDX(m,n)
            elm = getelm (i,i,a,ja,ia,iadd,.false.) 
            if(LCV(m,n).gt.0.5) then 
              j = VINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCV(m+1,n).gt.0.5) then 
              j = VINDX(m+1,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCV(m,n-1).gt.0.5) then 
              j = VINDX(m,n-1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCV(m+1,n-1).gt.0.5) then 
              j = VINDX(m+1,n-1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            
            !pattern to slh
            if(LU(m,n).gt.0.5) then 
              j = SLHINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LU(m+1,n).gt.0.5) then 
              j = SLHINDX(m+1,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            
          endif 

          !v row
          if(LCV(m,n).gt.0.5) then 

            i = VINDX(m,n)
            elm = getelm (i,i,a,ja,ia,iadd,.false.) 
            if(LCU(m,n).gt.0.5) then 
              j = UINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCU(m,n+1).gt.0.5) then 
              j = UINDX(m,n+1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCU(m-1,n).gt.0.5) then 
              j = UINDX(m-1,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCU(m-1,n+1).gt.0.5) then 
              j = UINDX(m-1,n+1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif

            !pattern to slh
            if(LU(m,n).gt.0.5) then 
              j = SLHINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LU(m,n+1).gt.0.5) then 
              j = SLHINDX(m,n+1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.le.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif


          endif
          !slh row
          if(LU(m,n).gt.0.5) then 

            i = SLHINDX(m,n)
            elm = getelm (i,i,a,ja,ia,iadd,.false.) 
            if(LCU(m,n).gt.0.5) then 
              j = UINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCU(m-1,n).gt.0.5) then 
              j = UINDX(m-1,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCV(m,n).gt.0.5) then 
              j = VINDX(m,n)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif
            if(LCV(m,n-1).gt.0.5) then 
              j = VINDX(m,n-1)
              elm = getelm (i,j,a,ja,ia,iadd,.false.) 
              if(iadd.eq.0)  then
                write(*,*) 'error: invalid matrice at m,n ',m,n
              endif
            endif

          endif
        enddo
      enddo


      WRITE(*,*) 'test matrice pattern: ok'
      end subroutine
!
!
!

      subroutine test_is_symmetric()
!      IMPLICIT NONE
      integer iadd1, iadd2
      integer m, n, issym, isantisymm
      double precision elm1, elm2
      logical sorted 
      double precision  getelm
      sorted = .false. 

      isantisymm = 1
      issym      = 1

      do m=1,nvar
        do n=1,m-1

           elm1 = getelm (m,n,a,ja,ia,iadd1,.false.) 
           elm2 = getelm (n,m,a,ja,ia,iadd2,.false.) 
           if (iadd1.ne.0.and.iadd2.eq.0) THEN
             WRITE(*,*) 'non-symmetric pattern m=', n, ' n =',n
             STOP 'test_is_symmettic'
           endif
           if (iadd2.ne.0.and.iadd1.eq.0) THEN
             WRITE(*,*) 'non-symmetric pattern m=', n, ' n =',n
             STOP 'test_is_symmettic'
           endif
           if (elm1.ne.elm2)  issym = 0
           if (elm1.ne.-elm2) isantisymm = 0
        enddo
      enddo
      
      WRITE(*,*) 'issym ', issym
      WRITE(*,*) 'issym ', isantisymm
      WRITE(*,*) 'test symmetric: ok'
      end subroutine

!c-----------------------------------------------------------------------------
!c
!c    super lu has to be called with symmetric mode flag and with amd reordering!!!!
!c
!c-----------------------------------------------------------------------------
!c-----------------------------------------------------------------------------
!      SUBROUTINE SOLVE_MATRICE_SUPERLU(U,V,SLH, ierr)
!      IMPLICIT NONE
!      REAL*8   U(NX,NY),V(NX,NY),SLH(NX,NY) ! put  solution here
!      INTEGER iopt, nrhs, ldb, info
!      DOUBLE PRECISION SOL(nvar)
!      DOUBLE PRECISION, ALLOCATABLE::A2(:)
!      INTEGER, ALLOCATABLE::IA2(:), JA2(:)
!      REAL TIME_IN_SEC,TIME_IN_SEC0,TIME_IN_SEC1
!      INTEGER ierr
!
!      ierr = -1
!
!      nrhs = 1 !one rhs vector at a time
!      ldb  = nvar
!      if (debug.gt.0) then
!        write (*,*) ' ******************************************* '
!        write (*,*) ' Start matrix solution by direct solver'
!    endif
!c
!c First, factorize the matrix. The factors are stored in factor() handle.
!c
!      if(isfactored.eq.0) THEN
!        write (*,*) ' ******************************************* '
!      write (*,*) ' Start matrix factorization by direct solver'
!        write (*,*) ' ******************************************* '
!
!c        convert matrice in csc format needed by super lu
! 
!        ALLOCATE (A2(maxnnz), JA2(maxnnz), IA2(nvar+1))
!        call csrcsc (nvar,1,1,a,ja,ia,a2,ja2,ia2)
!        A(1:maxnnz)  = A2
!        IA(1:nvar)   = IA2(1:nvar)
!        JA(1:maxnnz) = JA2
!        DEALLOCATE (IA2, JA2, A2)
!       
!        iopt = 1
!        call cpu_time(TIME_IN_SEC0)
!        call c_fortran_dgssv( iopt, nvar, nnz, nrhs, a, ja, ia,  
!     $                        rhs, ldb, factors, info )
!
!        isfactored = 1
!        call cpu_time(TIME_IN_SEC1) 
!        TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!        if (info .eq. 0) then
!c
!c         everything is ok
!c
!           if (debug.gt.0) then
!              write (*,*) 'LU decompositon takes ' ,   TIME_IN_SEC
!              write (*,*) 'Factorization succeeded'
!           endif
!           ierr = 0
!
!        else
!
!c
!c         error in super lu - free memory and go back
!c
!           write(*,*) 'INFO from triangular solve = ', info
!
!           iopt = 3 !free superlu memory
!           call c_fortran_dgssv( iopt, nvar, nnz, nrhs, a, ia, ja,  !attention transposed matrice
!     $                        rhs, ldb, factors, info )
!c
!c        convert matrice in back 
!c
!          ALLOCATE (A2(maxnnz), JA2(maxnnz), IA2(nvar+1))
!          call csrcsc (nvar,1,1,a,ja,ia,a2,ja2,ia2)
!          A(1:maxnnz)  = A2
!          IA(1:nvar)   = IA2
!          JA(1:maxnnz) = JA2
!          DEALLOCATE (IA2, JA2, A2)
!
!          STOP 'error in superLu'
!          return
!
!        endif
!
!      endif
!
!c
!c    Second, solve the system using the existing factors.
!      call cpu_time(TIME_IN_SEC0)
!
!      iopt = 2
!      call c_fortran_dgssv( iopt, nvar, nnz, nrhs, a, ia, ja,  !attention transposed matrice
!     $                        rhs, ldb, factors, info )
!      if (info .eq. 0) then
!       if (debug.gt.0) then
!           write (*,*) 'Solve succeeded'
!       end if
!      else
!         write(*,*) 'INFO from triangular solve = ', info
!      endif
!
!      if (debug.gt.0) then
!       call cpu_time(TIME_IN_SEC1) 
!         TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0
!         write (*,*) 'Triangle system solution takes ' ,   TIME_IN_SEC
!    end if
!      
!      SOL = RHS(1:nvar)
!
!      if (info .eq. 0) then 
!       if (debug.gt.0) then
!           write (*,*) 'Factorization succeeded'
!       end if
!      else
!         write(*,*) 'INFO from factorization = ', info
!      endif
!
!c--------------------------------------------------------------
!c     put solution into original order
!c--------------------------------------------------------------
!      call reorderToOriginalOrder(U,V,SLH,SOL)
!
!      end subroutine

!=================================================================
!
!     Free all allocatated arrays 
!
!----------------------------------------------------------------------     
      Subroutine Free_Matrice()
        INTEGER iopt, nrhs, info, ldb

        DEALLOCATE (guess_init, RHS, IA, JA, A)

        IS_FORMED   = 0

      if(PRECONDITIONER_FORMED.gt.0)  THEN
        DEALLOCATE (au, jau,iau )
      ENDIF
      PRECONDITIONER_FORMED = 0

        ldb  = nvar
        nrhs =1

!       free storage in super lu
      IF (ISFACTORED.GT.0) THEN
          iopt = 3
!          call c_fortran_dgssv( iopt, nvar, nnz, nrhs, a, ja, ia,  
!     $                        rhs, ldb, factors, info )
      end if
        ISFACTORED = 0

      end subroutine

!-------------------------------------------------------------- 
!
! form symbolic matrice for amd reordering . 
! see http://www.cise.ufl.edu/research/sparse/amd
! of Davis, Patrick R. Amestoy, and Iain S. Duff.  
! Note that the diagonal entries are not present,
!
      subroutine  do_amd_reorder_of_matrice()
 !     IMPLICIT NONE
      INTEGER J,K, IWLEN, PFREE, NCMPA, n_iw
      INTEGER, ALLOCATABLE :: LAST (:), PE (:), LENGTH (:), IW (:), DEGREE (:), NV (:), NEXT (:), HEAD (:), W (:)

      REAL*8, ALLOCATABLE :: AO(:)
      INTEGER, ALLOCATABLE ::IAO(:), JAO(:)

      REAL TIME_IN_SEC,TIME_IN_SEC0,TIME_IN_SEC1
      if (debug.eq.1) THEN 
          WRITE(*,*) 'begining of amd reordering '
      end if 
!
! see documemtation to amd for details
!      
      IWLEN = 1.2*NNZ+nvar !recommended value for array is 1.2 * nz+nvar
      ALLOCATE  (PERMUTATION(NVAR))
      ALLOCATE  ( LAST (NVAR)    )
      ALLOCATE  ( PE (NVAR)      )
      ALLOCATE  ( LENGTH (NVAR)  )
      ALLOCATE  ( IW (IWLEN)     )
      ALLOCATE  ( DEGREE (NVAR)  )
      ALLOCATE  ( NV (NVAR)      )
      ALLOCATE  ( NEXT (NVAR)    )
      ALLOCATE  ( HEAD (NVAR)    )
      ALLOCATE  ( W (NVAR)       )
     

      if (IS_FORMED.eq.0) THEN
          STOP 'do_amd_reorder_of_matrice: matrice is not formed yet'
      endif 

!
!       load the matrix into the AMD workspace. exclude diagonal
!
        N_IW = 1
        DO J = 1,NVAR

            PE (J) = N_IW

            LENGTH(J) = IA (J+1) - IA (J)
            LENGTH(J) = LENGTH(J) - 1

            DO  K = IA(J), IA(J+1)-1       
               IF( JA(K).NE.J ) THEN
                  IW (N_IW) = JA (K)
                  N_IW = N_IW + 1
               ENDIF
            END DO 
        END DO


        PFREE = NNZ + 1
!
!       order the matrix (destroys the copy of A in IW, PE, and LEN)
!
        call cpu_time(TIME_IN_SEC0)
        CALL AMD (NVAR, PE, IW, LENGTH, IWLEN, PFREE, NV, NEXT, LAST, HEAD, PERMUTATION, DEGREE, NCMPA, W)
        call cpu_time(TIME_IN_SEC1) 
        TIME_IN_SEC=TIME_IN_SEC1-TIME_IN_SEC0

        write (*,*) ' Time on amd reordering ', TIME_IN_SEC


        if (debug.eq.1) THEN 
          WRITE(*,*) 'do_amd: number of compression ',NCMPA
          WRITE(*,*) 'do_amd: the length of IW required for no ',  ' compressions to be needed ', PFREE
          WRITE(*,*) 'do_amd: IWLEN ', IWLEN
        endif

      DEALLOCATE  ( W        )
      DEALLOCATE  ( HEAD     )
      DEALLOCATE  ( NEXT     )
      DEALLOCATE  ( NV       )
      DEALLOCATE  ( DEGREE   )
      DEALLOCATE  ( IW       )
      DEALLOCATE  ( LENGTH   )
      DEALLOCATE  ( PE       )
      DEALLOCATE  ( LAST     )


 
     
!
!      reorder matrice. symmetric reordering
!
!       use sparskit subroutine

      ALLOCATE (AO(maxnnz), JAO(maxnnz), IAO(nvar+1))

      call dperm (NVAR,  A, JA, IA, AO, JAO, IAO, PERMUTATION,  PERMUTATION, 1)
      A  = AO 
      IA = IAO
      JA = JAO

      MATRICE_PERMUTED = 1.0
      if (debug.eq.1) THEN 
          WRITE(*,*) 'do_amd: matrice permuted succesfully'
      endif 

      DEALLOCATE (IAO, JAO, AO )

      end subroutine


!----------------------------------------------------------------------- 
!
!      stolen from sparskit and changed
!       Attention: is not tested
!----------------------------------------------------------------------- 
      subroutine invert_dvperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place inverse permutation of a real vector x 
! according to the permutation array perm(*), i.e., on return, 
! the vector x satisfies,
!
!    x(j) :== x(perm(j)), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n     = length of vector x.
! perm     = integer array of length n containing the permutation  array.
! x    = input vector
!
! on return:
!---------- 
! x    = vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables 
      real*8 tmp, tmp1
      integer ii, init, k, next, j
!
      init      = 1
      ii        = perm(init)
      tmp    = x(ii)    
      perm(init)= -perm(init)
      k         = 0
!     
! loop
! 
 6    k = k+1
!
! save the chased element --
! 
      tmp1      = x(init) 
      x(init)     = tmp
      next      = perm(ii) 
      if (next .lt. 0 ) goto 65
!     
! test for end 
!
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
!
! end loop 
!
      goto 6
!
! reinitilaize cycle --
!
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      ii    = perm(init)
      tmp    = x(ii)
      perm(init)=-perm(init)
      goto 6
!     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
!     
      return
!-------------------end-of-dvperm--------------------------------------- 
!-----------------------------------------------------------------------
      end subroutine
!======================================================================
      endmodule sea_level_no_split

!  Additional function
!-----------------------------------------------------------------------
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
!-----end-of-distdot
!-----------------------------------------------------------------------
      end function
