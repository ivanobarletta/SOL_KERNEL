
   SUBROUTINE standard_cg( kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_pcg  ***
      !!                        standard CG
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solpcg
      !!
      INTEGER  ::   ji, jj, jn,k   ! dummy loop indices
      REAL(wp) ::   zgcad        ! temporary scalars
      REAL(wp), DIMENSION(2) ::   zsum
      REAL(wp), POINTER, DIMENSION(:,:) ::   zgcr
      !!----------------------------------------------------------------------
      !
      !CALL wrk_alloc( jpi, jpj, zgcr )
      !
      ! Initialization of the algorithm with standard PCG
      ! -------------------------------------------------
      !zgcr = 0._wp
      gcr  = 0._wp

      CALL mpp_lnk_2d( gcx, c_solver_pt, 1. )   ! lateral boundary condition

      ! p0 = r_0 = b - A*x_0
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
               &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
               &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
               &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
               &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
            gcr  (ji,jj) = zgcad
            gcdes(ji,jj) = zgcad
         END DO
      END DO

      DO k = 1, nn_nmax                               ! Iterative loop

         ! rnorme = <r_k-1,r_k-1>
         rnorme = glob_sum_2d(  gcr(:,:) * gcr(:,:)  )
         IF ( mpprank == 0 ) PRINT *, 'pcg iter',k,'rnorme',rnorme

         CALL mpp_lnk_2d( gcdes, c_solver_pt, 1. )   ! lateral boundary condition

         ! d_k-1 = A*p_k-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
                  &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
                  &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
            END DO
         END DO 

         ! alph = <r_k-1|r_k-1>/<p_k-1|d_k-1>
         radd = glob_sum_2d(  gcdes(:,:)  * gccd(:,:)  )
         alph = rnorme /radd

         ! x_k = x_k-1 + alph_k-1 * p_k-1
         ! r_k = r_k-1 - alph_k-1 * d_k-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
               gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
            END DO
         END DO

         rr = rnorme
         ! <r_k|r_k>
         rnorme = glob_sum_2d( gcr(:,:) * gcr(:,:) )
         ! beta_k = <r_k|r_k> / <r_k-1|r_k-1>
         beta = rnorme / rr

         ! test of convergence
         IF( rnorme < epsr .OR. k == nn_nmax ) THEN
            IF ( mpprank == 0 ) PRINT *, k,'pcg iterations' 
            res = SQRT( rnorme )
            niter = k
            ncut = 999
         ENDIF

         ! p_k = r_k + beta_k * p_k-1 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gcdes(ji,jj) = gcr (ji,jj) + beta * gcdes(ji,jj) 
           END DO
         END DO

         ! indicator of non-convergence or explosion
         IF( k == nn_nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999

      END DO

999   CONTINUE
          
      CALL mpp_lnk_2d( gcx, c_solver_pt, 1. )      ! Output in gcx with lateral b.c. applied
      ! 
      !CALL wrk_dealloc( jpi, jpj, zgcr )

   END SUBROUTINE standard_cg

   SUBROUTINE standard_cg2( kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_pcg  ***
      !!                        standard CG
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solpcg
      !!
      INTEGER  ::   ji, jj, jn,k   ! dummy loop indices
      REAL(wp) ::   zgcad        ! temporary scalars
      REAL(wp), DIMENSION(2) ::   zsum
      REAL(wp), POINTER, DIMENSION(:,:) ::   zgcr
      !!----------------------------------------------------------------------
      !
      !CALL wrk_alloc( jpi, jpj, zgcr )
      !
      ! Initialization of the algorithm with standard PCG
      ! -------------------------------------------------
      !zgcr = 0._wp
      gcr  = 0._wp

      CALL mpp_lnk_2d( gcx, c_solver_pt, 1. )   ! lateral boundary condition

      ! p1 = r_1 = b - A*x_1
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
               &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
               &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
               &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
               &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
            gcr  (ji,jj) = zgcad
            gcdes(ji,jj) = zgcad
         END DO
      END DO

      DO k = 1, nn_nmax                               ! Iterative loop

         gammak = glob_sum_2d( gcr(:,:) * gcr(:,:) * bmask(:,:) ) 

         IF ( mpprank == 0 ) PRINT *, 'pcg iter',k,'gammak',gammak

         IF ( k == 1) THEN
            beta = 0.
            gcdes(:,:) = gcr(:,:)
         ELSE
            beta = gammak / gammakm1 
            ! p_k = r_k + beta_k * p_k-1 
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
               gcdes(ji,jj) = gcr (ji,jj) + beta * gcdes(ji,jj) 
               END DO
            END DO
         END IF    

         CALL mpp_lnk_2d( gcdes, c_solver_pt, 1. )   ! lateral boundary condition

         ! d_k = A*p_k
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
                  &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
                  &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
            END DO
         END DO 

         ! alph = <r_k-1|r_k-1>/<p_k-1|d_k-1>
         sigmak = glob_sum_2d(  gcdes(:,:) * gccd(:,:) * bmask (:,:) )
         alph = gammak / sigmak

         ! x_k+1 = x_k + alph_k * p_k
         ! r_k+1 = r_k - alph_k * d_k
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
               gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
            END DO
         END DO

         gammakm1 = gammak

         ! test of convergence
         IF( gammak < epsr .OR. k == nn_nmax ) THEN
            IF ( mpprank == 0 ) PRINT *, k,'pcg iterations' 
            res = SQRT( gammak )
            niter = k
            ncut = 999
         ENDIF

         ! indicator of non-convergence or explosion
         IF( k == nn_nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999

      END DO

999   CONTINUE
          
      CALL mpp_lnk_2d( gcx, c_solver_pt, 1. )      ! Output in gcx with lateral b.c. applied
      ! 
      !CALL wrk_dealloc( jpi, jpj, zgcr )

   END SUBROUTINE standard_cg2

   SUBROUTINE sol_pcg_sstep( kindic )
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solpcg
      !!
      INTEGER  ::   ji, jj, jn, jn2 , inum  ! dummy loop indices
      INTEGER  ::   itercount, yy
      INTEGER  ::   jimin, jimax, jjmin, jjmax
      REAL(wp) ::   zgcad, norme_util,local_rnorme        ! temporary scalars
      REAL(wp) ::   local_gcr, local_gcdmat               ! temporary scalars
      REAL(wp) ::   local_gcdes, local_gccd, local_zgcad  ! temporary scalars
      REAL(wp) ::   local_zgcr                            ! temporary scalars
      REAL(wp), DIMENSION(2) ::   zsum
      !REAL(wp), POINTER, DIMENSION(:,:) ::   zgcr
      !!----------------------------------------------------------------------
      !
      !CALL wrk_alloc( jpi, jpj, zgcr )
      !
      ! Initialization of the algorithm with standard PCG
      ! -------------------------------------------------
      !zgcr = 0._wp
      gcr  = 0._wp
      yy   = 250

      CALL mpp_lnk_2d_e( gcx, c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary condition

      ! gcdes = gcr = gcb - A*gcx
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
               &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
               &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
               &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
               &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
            gcr  (ji,jj) = zgcad
            gcdes(ji,jj) = zgcad
         END DO
      END DO

      ! rnorme = <gcr|gcr>
      !rnorme = glob_sum_2d(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
      !rnorme = glob_sum_2d( gcr(1:jpi,1:jpj) * gcdmat(1:jpi,1:jpj) * gcr(1:jpi,1:jpj) )
      rnorme = glob_sum_2d( gcr(2:jpim1,2:jpjm1) * gcdmat(2:jpim1,2:jpjm1) * gcr(2:jpim1,2:jpjm1) )
      
      local_gcr = SUM(gcr(2:jpim1,2:jpjm1))  
      IF (mpprank == rank_print) PRINT '(a,f23.16)','local_sum gcr   ',local_gcr                      
      IF (mpprank == rank_print) PRINT '(a,f23.16)','1st rnorme      ',rnorme                      

      CALL mpp_lnk_2d_e( gcdes, c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary condition

      ! gccd = A*gcdes
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
               &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
               &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
         END DO
      END DO 

      CALL mpp_lnk_2d_e( gccd , c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary condition

      !IF ( mpprank == rank_print ) THEN
      !   CALL iom_open('gccd_sstep',inum,.TRUE.)
      !   CALL iom_rp2d(1,1,inum,'gccd',gccd(1:jpi,1:jpj))
      !   CALL iom_rp2d(1,1,inum,'gcdes',gcdes(1:jpi,1:jpj))
      !   CALL iom_rp2d(1,1,inum,'bmask',bmask(1:jpi,1:jpj))
      !   CALL iom_close(inum)
      !END IF

      local_gcdes = SUM(gcdes(2:jpim1,2:jpjm1))
      local_gccd  = SUM(gccd (2:jpim1,2:jpjm1))
      IF (mpprank == rank_print) PRINT '(a,f23.16)','local_sum gcdes ',local_gcdes                      
      IF (mpprank == rank_print) PRINT '(a,f23.16)','local_sum gccd  ',local_gccd

      ! alph = <gcr|gcr>/<gcdes|gccd>
      radd = glob_sum_2d(  gcdes(2:jpim1,2:jpjm1) * gcdmat(2:jpim1,2:jpjm1) * gccd(2:jpim1,2:jpjm1)  )
      alph = rnorme /radd

      IF (mpprank == rank_print) PRINT '(a,f23.15)','1st radd',radd
      IF (mpprank == rank_print) PRINT '(a,f23.15)','1st alph',alph

      ! gcx = gcx + alph * gcdes
      ! gcr = gcr - alph * gccd
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
            gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
         END DO
      END DO

      ! Algorithm with Eijkhout rearrangement
      ! and s-step approach
      ! -------------------------------------
      itercount = 0        
      !                                                !================
      DO jn = 1, nn_nmax / sstep                       ! Iterative loop
         !                                             !================
         CALL mpp_lnk_2d_e( gcr, c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary condition  with sstep-1 extra-halo
         !CALL mpp_lnk_2d_e( gccd, c_solver_pt, 1., jpr2di, jpr2dj )   ! lateral boundary condition  with sstep-1 extra-halo

         DO jn2 = 1, sstep                             ! Inner Loop 
            
            IF ( jn2 == 1) THEN
            CALL mpp_lnk_2d( gccd, c_solver_pt, 1. )   ! lateral boundary condition  with sstep-1 extra-halo
            END IF 
            itercount = itercount + 1 
            !IF ( mpprank == 0 ) PRINT *, 'pcg iter',itercount,'rnorme',rnorme
            !IF ( mpprank == 0 ) PRINT '(a,i3,a,e17.11,a,e14.8)', 'pcg iter',itercount,' rnorme ',rnorme,' beta ',beta
            ! zgcr = matrix . gcr
            jjmin = 2-sstep+jn2         
            jjmax = jpjm1+sstep-jn2  
            jimin = 2-sstep+jn2      
            jimax = jpim1+sstep-jn2  
            !IF (mpprank == rank_print) PRINT'(a,i4,a,i4,a,i4,a,i4)','jj1:',jjmin,' jj2:',jjmax,' ji1:',jimin,' ji2:',jimax 
            !IF (jn == 1 .and. mpprank == 0) PRINT *,jjmin,jjmax,jimin,jimax
            DO jj = jjmin, jjmax
               DO ji = jimin,jimax   ! vector opt.
                  !zzgcr(ji,jj) = bmask(ji,jj)*( gcr(ji,jj)   &
                  zzgcr(ji,jj) = ( gcr(ji,jj)   &
                     &        +gcp(ji,jj,1)*gcr(ji,jj-1)+gcp(ji,jj,2)*gcr(ji-1,jj)   &
                     &        +gcp(ji,jj,4)*gcr(ji,jj+1)+gcp(ji,jj,3)*gcr(ji+1,jj)   )
               END DO
            END DO
 
            ! rr = rnorme = <gcr|gcr>
            rr = rnorme

            ! MPI_ALLREDUCE CALLS
            ! rnorme = <gcr|gcr> 
            ! zgcad  = <zgcr|gcr> 
            !rnorme = glob_sum_2d(gcr(1:jpi,1:jpj)*gcdmat(1:jpi,1:jpj)*   gcr(1:jpi,1:jpj))
            !zgcad  = glob_sum_2d(gcr(1:jpi,1:jpj)*gcdmat(1:jpi,1:jpj)* zzgcr(1:jpi,1:jpj) * bmask(1:jpi,1:jpj))
            !rnorme = glob_sum_2d(gcr(:,:)*gcdmat(:,:)*   gcr(:,:))
            !zgcad  = glob_sum_2d(gcr(:,:)*gcdmat(:,:)* zzgcr(:,:) * bmask(:,:))
            !rnorme = glob_sum_2d(gcr(jimin-1:jimax+1,jjmin-1:jjmax+1)*gcdmat(jimin-1:jimax+1,jjmin-1:jjmax+1)* & 
            !                        & gcr(jimin-1:jimax+1,jjmin-1:jjmax+1))
            !zgcad  = glob_sum_2d(gcr(jimin-1:jimax+1,jjmin-1:jjmax+1)*gcdmat(jimin-1:jimax+1,jjmin-1:jjmax+1)* & 
            !                        & zzgcr(jimin-1:jimax+1,jjmin-1:jjmax+1) * bmask(jimin-1:jimax+1,jjmin-1:jjmax+1))
            rnorme = glob_sum_2d(gcr(2:jpim1,2:jpjm1)*gcdmat(2:jpim1,2:jpjm1)*   gcr(2:jpim1,2:jpjm1))
            zgcad  = glob_sum_2d(gcr(2:jpim1,2:jpjm1)*gcdmat(2:jpim1,2:jpjm1)* zzgcr(2:jpim1,2:jpjm1) & 
                                                                           & * bmask(2:jpim1,2:jpjm1))

            local_zgcr   = SUM(zzgcr (2:jpim1,2:jpjm1))
            local_gcr    = SUM(gcr   (2:jpim1,2:jpjm1))
            local_gcdmat = SUM(gcdmat(2:jpim1,2:jpjm1))  
            local_rnorme = SUM(gcr(2:jpim1,2:jpjm1)*gcdmat(2:jpim1,2:jpjm1) * gcr(2:jpim1,2:jpjm1))
            local_zgcad  = SUM(gcr(2:jpim1,2:jpjm1)*gcdmat(2:jpim1,2:jpjm1) * zzgcr(2:jpim1,2:jpjm1) * bmask(2:jpim1,2:jpjm1))
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11)','pcg iter',itercount,' local_sum_zgcr   ',local_zgcr
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11)','pcg iter',itercount,' local_sum_gcr    ',local_gcr
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11,a,e14.8)','pcg iter',itercount,' local_gcdmat ',local_gcdmat
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11,a,e14.8)','pcg iter',itercount,' local_rnorme   ',local_rnorme
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11,a,e14.8)','pcg iter',itercount,' local_zgcad    ',local_zgcad
            223 FORMAT(a,2(e13.7,1x),a,4(e13.7,1x),a,2(e13.7,1x))
            IF (.TRUE.) THEN 
            IF ( mpprank == rank_print ) THEN
            PRINT*,'_____________________________________________________' 
            PRINT '(4(a,i3,15x))','pcg iter ',itercount,' pcg iter ',itercount,' pcg iter ',itercount,' pcg iter', itercount
            PRINT '(11a)'              ,'         ','       0      ','      1       ','|','       2      ','       3      ', &
                                   &                '    jpi-2     ','     jpi-1    ','|','     jpi      ','    jpi+1     '   
            PRINT 223                  ,' gcr     ',   gcr(0,yy)    ,   gcr(1,yy)    ,'|',   gcr(2,yy)    ,      gcr(3,yy) , &
                                                   &   gcr(jpi-2,yy),  gcr(jpi-1,yy) ,'|',  gcr(jpi,yy)   ,  gcr(jpi+1,yy)   
            PRINT 223                  ,'zgcr     ', zzgcr(0,yy)    , zzgcr(1,yy)    ,'|', zzgcr(2,yy)    ,    zzgcr(3,yy) , &
                                                   & zzgcr(jpi-2,yy),zzgcr(jpi-1,yy) ,'|',zzgcr(jpi,yy)   ,zzgcr(jpi+1,yy) 
            PRINT 223                  ,'gcdes    ', gcdes(0,yy)    , gcdes(1,yy)    ,'|', gcdes(2,yy)    ,    gcdes(3,yy) , &
                                                   & gcdes(jpi-2,yy),gcdes(jpi-1,yy) ,'|',gcdes(jpi,yy)   ,gcdes(jpi+1,yy) 
            PRINT 223                  ,'gccd     ',  gccd(0,yy)    ,  gccd(1,yy)    ,'|',  gccd(2,yy)    ,     gccd(3,yy) , &
                                                   &  gccd(jpi-2,yy), gccd(jpi-1,yy) ,'|', gccd(jpi,yy)   , gccd(jpi+1,yy) 
            END IF 
            END IF
            ! test of convergence
            ! test of convergence
            IF( rnorme < epsr .OR. itercount == nn_nmax ) THEN
               !IF ( mpprank == 0 ) PRINT *, jn,'pcg iterations' 
               res = SQRT( rnorme )
               niter = jn
               ncut = 999
            ENDIF

            ! beta = <rk+1|rk+1>/<rk|rk>
            beta = rnorme / rr
            radd = zgcad - beta*beta*radd
            alph = rnorme / radd

            IF (.TRUE. ) THEN
            IF ( mpprank == rank_print) THEN
               PRINT '(a,e25.15)', '  - rr - ', rr
               PRINT '(a,e25.15)', '-rnorme- ', rnorme
               PRINT '(a,e25.15)', ' -zgcad- ', zgcad
               PRINT '(a,e25.15)', ' -radd - ', radd
               PRINT '(a,e25.15)', '  -beta- ', beta
               PRINT '(a,e25.15)', '  -alph- ', alph
            END IF          
            END IF          

            ! gcx = gcx + alph * gcdes
            ! gcr = gcr - alph * gccd
            DO jj = jjmin, jjmax
               DO ji = jimin,jimax   ! vector opt.
                  gcdes(ji,jj) = gcr  (ji,jj)  + beta * gcdes(ji,jj) 
                  gccd (ji,jj) = zzgcr(ji,jj)  + beta * gccd (ji,jj) 
                  gcx  (ji,jj) = gcx  (ji,jj)  + alph * gcdes(ji,jj) 
                  gcr  (ji,jj) = gcr  (ji,jj)  - alph * gccd (ji,jj) 
               END DO
            END DO
           
            IF (.TRUE.) THEN
            IF ( mpprank == rank_print ) THEN
            PRINT 223                  ,'gcdes aft', gcdes(0,yy)    , gcdes(1,yy)    ,'|', gcdes(2,yy)    ,    gcdes(3,yy) , &
                                                   & gcdes(jpi-2,yy),gcdes(jpi-1,yy) ,'|',gcdes(jpi,yy)   ,gcdes(jpi+1,yy) 
            PRINT 223                  ,'gccd  aft',  gccd(0,yy)    ,  gccd(1,yy)    ,'|',  gccd(2,yy)    ,     gccd(3,yy) , &
                                                   &  gccd(jpi-2,yy), gccd(jpi-1,yy) ,'|', gccd(jpi,yy)   , gccd(jpi+1,yy) 
            PRINT 223                  ,' gcr  aft',   gcr(0,yy)    ,   gcr(1,yy)    ,'|',   gcr(2,yy)    ,      gcr(3,yy) , &
                                                   &   gcr(jpi-2,yy),  gcr(jpi-1,yy) ,'|',  gcr(jpi,yy)   ,  gcr(jpi+1,yy)   
            END IF 
            END IF                            
 
            local_gcr    = SUM(gcr   (2:jpim1,2:jpjm1))
            !IF ( mpprank == rank_print ) PRINT '(a,i3,a,e17.11,a,e14.8)','pcg iter',itercount,' local_gcr_after ',local_gcr

                                                       !===============
         END DO                                        ! End SSTEP    Loop  
                                                       !================

         ! indicator of non-convergence or explosion
         IF( jn == nn_nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999

                                                       !================
      END DO                                           ! End External Loop
      !                                                !================
999   CONTINUE
          
      !CALL lbc_lnk( gcx, c_solver_pt, 1. )      ! Output in gcx with lateral b.c. applied
      CALL mpp_lnk_2d_e( gcx, c_solver_pt, 1., jpr2di, jpr2dj )      ! Output in gcx with lateral b.c. applied
      CALL mpp_lnk_2d_e( gcr, c_solver_pt, 1., jpr2di, jpr2dj )      ! Output in gcr with boundary condition
      ! 
      !CALL wrk_dealloc( jpi, jpj, zgcr )
      !
      !IF( nn_timing == 1 )  CALL timing_stop('sol_pcg')
      !
   END SUBROUTINE sol_pcg_sstep
