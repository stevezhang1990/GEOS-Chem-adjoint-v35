      MODULE HIPPO_MOD
      
      IMPLICIT NONE

      INTEGER :: N_HIPPO
      REAL*8, ALLOCATABLE :: YEAR_HIPPO(:), DAY_HIPPO(:), TIME_HIPPO(:)
      REAL*8, ALLOCATABLE :: LON_HIPPO(:), LAT_HIPPO(:), ALT_HIPPO(:)
      REAL*8, ALLOCATABLE :: O3_HIPPO(:), O3_HIPPO_GC(:), MONTH_HIPPO(:)

      

      CONTAINS
      
      SUBROUTINE CLEANUP_HIPPO

      !! Clean up HIPPO arrays
      
      IF (ALLOCATED(YEAR_HIPPO)) DEALLOCATE(YEAR_HIPPO)
      IF (ALLOCATED(MONTH_HIPPO)) DEALLOCATE(MONTH_HIPPO)
      IF (ALLOCATED(DAY_HIPPO)) DEALLOCATE(DAY_HIPPO)
      IF (ALLOCATED(TIME_HIPPO)) DEALLOCATE(TIME_HIPPO)
      IF (ALLOCATED(LON_HIPPO)) DEALLOCATE(LON_HIPPO)
      IF (ALLOCATED(LAT_HIPPO)) DEALLOCATE(LAT_HIPPO)
      IF (ALLOCATED(ALT_HIPPO)) DEALLOCATE(ALT_HIPPO)
      IF (ALLOCATED(O3_HIPPO)) DEALLOCATE(O3_HIPPO)
      IF (ALLOCATED(O3_HIPPO_GC)) DEALLOCATE(O3_HIPPO_GC)

      END SUBROUTINE CLEANUP_HIPPO

      SUBROUTINE INIT_HIPPO
      
      !! Initialize HIPPO arrays

      !N_HIPPO = 18951
      !N_HIPPO = 18640
      !N_HIPPO = 21812
      N_HIPPO = 136
      ALLOCATE(YEAR_HIPPO(N_HIPPO))
      ALLOCATE(MONTH_HIPPO(N_HIPPO))
      ALLOCATE(DAY_HIPPO(N_HIPPO))
      ALLOCATE(TIME_HIPPO(N_HIPPO))
      ALLOCATE(LON_HIPPO(N_HIPPO))
      ALLOCATE(LAT_HIPPO(N_HIPPO))
      ALLOCATE(ALT_HIPPO(N_HIPPO))
      ALLOCATE(O3_HIPPO(N_HIPPO))
      ALLOCATE(O3_HIPPO_GC(N_HIPPO))

!     Initialize GC-O3 values to missing value
      O3_HIPPO_GC = -999.99

      END SUBROUTINE INIT_HIPPO

      SUBROUTINE READ_HIPPO_DATA

      !! Read HIPPO data from disk

      INTEGER :: I

      OPEN(UNIT=24,
     & FILE="/users/jk/07/xzhang/ACE_FTS_NO2/ace_fts_no2.csv")

      DO I=1,N_HIPPO
     
         READ(24, *) YEAR_HIPPO(I), MONTH_HIPPO(I), DAY_HIPPO(I), 
     &               TIME_HIPPO(I), LON_HIPPO(I), LAT_HIPPO(I), 
     &               ALT_HIPPO(I), O3_HIPPO(I)
     
      END DO

      CLOSE(UNIT=24)

      END SUBROUTINE READ_HIPPO_DATA

      SUBROUTINE WRITE_HIPPO_DATA

      !! Write HIPPO data to disk

      INTEGER :: I

      OPEN(UNIT=25,FILE="no2_gc.csv")

      DO I=1,N_HIPPO
     
         WRITE(25, *) YEAR_HIPPO(I), MONTH_HIPPO(I), DAY_HIPPO(I), 
     &                TIME_HIPPO(I), LON_HIPPO(I), LAT_HIPPO(I), 
     &                ALT_HIPPO(I), O3_HIPPO_GC(I)
     
      END DO

      CLOSE(UNIT=25)

      END SUBROUTINE WRITE_HIPPO_DATA

      SUBROUTINE COMPARE_HIPPO_DATA
      
      USE GRID_MOD,           ONLY : GET_IJ
      USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
      USE TIME_MOD,           ONLY : GET_MINUTE, GET_HOUR, GET_DAY
      USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      USE COMODE_MOD,         ONLY : JLOP
      USE COMODE_MOD,         ONLY : CSPEC
      USE ADJ_ARRAYS_MOD,     ONLY : ID2C
      USE TRACERID_MOD,       ONLY : IDNO2
      USE TRACERID_MOD,       ONLY : IDTOX
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      !USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AIRDEN
      USE DAO_MOD,            ONLY : AD
      USE TRACER_MOD,         ONLY : TCVV
      USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
      USE ADJ_ARRAYS_MOD,     ONLY : ID2C
      USE DIAG_MOD,           ONLY : AD43, LTNO2, CTNO2

#include "CMN_SIZE"             ! size parameters
#     include "CMN"             ! IFLX, LPAUSE
#     include "CMN_DIAG"        ! Diagnostic switches & arrays
!#     include "CMN_O3"          ! FMOL, XNUMOL
!#     include "comode.h"        ! IDEMS

      INTEGER                   :: I,J,L,K,M,N,NN
      INTEGER                   :: IIJJ(2), L_HIPPO, JLOOP
      REAL*8                    :: YEAR, MONTH, DAY, HOUR, MINUTE, TIME 
      REAL*8                    :: HEIGHT
      REAL*8, PARAMETER         :: ADJ_TCVVOX = 28.97d0/48.d0
      
      MINUTE = REAL(GET_MINUTE(),8)
      HOUR = REAL(GET_HOUR(),8)
      DAY = REAL(GET_DAY(),8)
      YEAR = REAL(GET_YEAR(),8)
      MONTH = REAL(GET_MONTH(),8)
      
      TIME = MINUTE + 100*HOUR
      DO K=1,N_HIPPO
         !PRINT *, "MONTH_HIPPO(K)", MONTH_HIPPO(K)
         
         IF( (DAY == DAY_HIPPO(K)) .AND.
     &       (YEAR == YEAR_HIPPO(K)) .AND.
     &       (MONTH == MONTH_HIPPO(K)) .AND. 
     &       ((TIME_HIPPO(K)-TIME) >= 0 ) .AND. 
     &       ((TIME_HIPPO(K)-TIME) < 30 ) ) THEN

            IIJJ = GET_IJ(REAL(LON_HIPPO(K),4),REAL(LAT_HIPPO(K),4))
            
            !PRINT *,"K_HIPPO",K
            !PRINT *,"TIME_HIPPO(K) - TIME", TIME_HIPPO(K) - TIME
            !PRINT *,"ALT_HIPPO(K)", ALT_HIPPO(K)
            

            I = IIJJ(1)
            J = IIJJ(2)
            
            HEIGHT = 0d0
            
            DO L=1,LLPAR
            
               !PRINT *,"ITS_IN_THE_TROP",ITS_IN_THE_TROP(I,J,L)

               HEIGHT = HEIGHT + BXHEIGHT(I,J,L)

               !PRINT *,"HEIGHT-ALT_HIPPO(K)",  HEIGHT - ALT_HIPPO(K)

               IF ( (HEIGHT .GE. ALT_HIPPO(K)) ) THEN
                 
               IF ( ITS_IN_THE_TROP(I,J,L) ) THEN
                  JLOOP = JLOP(I,J,L)     
                  O3_HIPPO_GC(K) = CSPEC(JLOOP,IDNO2)  
                  O3_HIPPO_GC(K) =  O3_HIPPO_GC(K) * 1d6 /
     &                      ( AIRDEN(L,I,J)   * XNUMOLAIR )
                     !O3_HIPPO_GC(K) = CHK_STT(I,J,L,IDTOX)* TCVV(IDTOX)/
     &                                 !AD(I,J,L) * 1d9
                  !PRINT *,"O3_HIPPO_GC(K):",O3_HIPPO_GC(K)
               ENDIF
                
               EXIT
               
               ENDIF

            ENDDO
            
         ENDIF

      ENDDO

      END SUBROUTINE COMPARE_HIPPO_DATA

      END MODULE HIPPO_MOD
