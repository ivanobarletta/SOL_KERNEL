
   SUBROUTINE sstep_pcg( kindic )
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solpcg
      !!
      INTEGER  ::   ji, jj, jn, jn2 , inum  ! dummy loop indices
      INTEGER  ::   itercount, yy
      INTEGER  ::   jimin, jimax, jjmin, jjmax
      REAL(wp) ::   zgcscal,W,B,C
      !
      ! gcr   = gcb-a.gcx
      DO jj =    2,jpjm1
         DO ji = 2,jpim1
            zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
               &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
               &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
               &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
               &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
            gcr  (ji,jj) = gcdprc(ji,jj) * zgcad
         END DO
      END DO
   
      itercount = 0
      ! OUTER LOOP 
      DO jn = 1, nn_nmax / sstep

         CALL mpp_lnk_2d_e( gcr, c_solver_pt, 1. , jpr2di, jpr2dj )  

         zaus(:,:,1) = gcr(:,:)

         ! INNER LOOP
         DO jn2 = 1, sstep

            itercount = itercount + 1

            ! MV-product 
            DO jj =    2-sstep+jn2, jpjm1+sstep-jn2
               DO ji = 2-sstep+jn2, jpim1+sstep-jn2
                   zgcscal    =           ( gcr(ji,jj)   &
                     &        +gcp(ji,jj,1)*gcr(ji,jj-1)+gcp(ji,jj,2)*gcr(ji-1,jj)   &
                     &        +gcp(ji,jj,4)*gcr(ji,jj+1)+gcp(ji,jj,3)*gcr(ji+1,jj)   )
                   IF ( jn2 .LT. sstep ) THEN
                      zaus(ji,jj,jn2+1) = zgcscal
                   END IF
                   qaus(ji,jj,jn2  ) = zgcscal
               END DO
            END DO 

         ! INNER LOOP
         END DO
               
         IF ( jn == 1 ) THEN
            paus (:,:,:) = zaus (:,:,:)
         ELSE
            ! calculate W,C and update P   
               
         END IF            
          
         ! Calculate W
         W = 0._wp 
         DO jn2 = 1,sstep
            zgcscal  = SUM( qaus(:,:,jn2) * paus(:,:,jn2) )
            W = W + zgcscal  
         END DO                  
         ! call mpi_allreduce for W  

         ! Calculate ggg         
         DO jn2 = 1,sstep
            ggg(jn2) = SUM( paus(:,:,jn2) * gcr(:,:)      ) 
         END DO                  


      ! OUTER LOOP
      END DO

   END SUBROUTINE sstep_pcg
