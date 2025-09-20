C $Header: /u/gcmpack/MITgcm/pkg/iceberg/iceberg_init_fixed.F,v 1.5 2011/12/28 20:45:05 jmc Exp $
C $Name: checkpoint65m $

#include "ICEBERG_OPTIONS.h"

      SUBROUTINE ICEBERG_INIT_FIXED( myThid )
C     *============================================================*
C     | SUBROUTINE ICEBERG_INIT_FIXED
C     | o Routine to initialize ICEBERG parameters and variables.
C     *============================================================*
C     | Initialize ICEBERG parameters and variables.
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#endif /* ALLOW_PTRACERS */
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"
#ifdef NONLIN_FRSURF
#include "SURFACE.h"
#endif /* NONLIN_FRSURF */

C     === Routine arguments ===
C     myThid -  Number of this instance of ICEBERG_INIT_FIXED
      INTEGER myThid

#ifdef ALLOW_ICEBERG
C     === Local variables ===
C     I,J,K,bi,bj - Loop counters, helpers for hFactors
      INTEGER I, J, K, bi, bj, Im1, Jm1, iTracer

C      Read iceberg mask
       IF (ICEBERGmaskFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGmaskFile,icebergMask,1,0,myThid)
       ENDIF

C      Read iceberg numbers file
       IF (ICEBERGmaskNumsFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGmaskNumsFile,icebergMaskNums,1,0,myThid)
       ENDIF
   
C      Read iceberg num per cell file
       IF (ICEBERGnumPerCellFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGnumPerCellFile,icebergNumBergs,1,0,myThid)
       ENDIF

C      Read melt mask
       IF (ICEBERGmeltFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGmeltFile,meltMask,1,0,myThid)
       ENDIF

C      Read drift mask
       IF (ICEBERGdriftFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGdriftFile,driftMask,1,0,myThid)
       ENDIF
   
C      Read barrier mask
       IF (ICEBERGbarrierFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGbarrierFile,barrierMask,1,0,myThid)
       ENDIF
   
C      Read iceberg area
       IF (ICEBERGareaFile .NE. ' ' ) THEN
        CALL READ_REC_XYZ_RL
     &     (ICEBERGareaFile,icebergArea3D,1,0,myThid)
       ENDIF

C      Read iceberg lengths
       IF (ICEBERGlengthFile .NE. ' ' ) THEN
        CALL READ_REC_3D_RL
     &     (ICEBERGlengthFile,precFloat64,
     &         maxBergCt,icebergLength,1,0,myThid) 
       ENDIF

C      Read iceberg widths
       IF (ICEBERGwidthsFile .NE. ' ' ) THEN
        CALL READ_REC_3D_RL
     &     (ICEBERGwidthsFile,precFloat64,
     &         maxBergCt,icebergWidths,1,0,myThid)
       ENDIF

C      Read iceberg depths
       IF (ICEBERGdepthsFile .NE. ' ' ) THEN
        CALL READ_REC_3D_RL
     &     (ICEBERGdepthsFile,precFloat64,
     &         maxBergCt,icebergDepths,1,0,myThid)
       ENDIF

C      Read open fraction file
      IF ( ICEBERGopenFracFile .NE. ' ' ) THEN
        CALL READ_REC_XYZ_RL
     &     (ICEBERGopenFracFile,openFraction,1,0,myThid)  
 
C          Partially block cells filled by icebergs, hFacC cell center
           DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
             DO J = 1, sNy
              DO I = 1, sNx
  
C             If we want to make icebergs physical barriers,
C             set both hFacC and recip_hFacC.
C               REQUIRES CPP_OPTIONS.h: #undef NONLIN_FRSURF
C             otherwise this is overwritten everytimestep.
               IF ( barrierMask(I,J,bi,bj) .GT. 0. _d 0 ) THEN
C              Loop through water column
                DO K = 1, Nr !1 is the surface
                  hFacC(I,J,K,bi,bj) = openFraction(I,J,K,bi,bj)
                 IF(hFacC(I,J,K,bi,bj).LT.hFacMin) THEN
                    hFacC(I,J,K,bi,bj) = 0. _d 0
                    recip_hFacC(I,J,K,bi,bj) = 0. _d 0
                 ELSE
                    recip_hFacC(I,J,K,bi,bj) = 
     &                             1. _d 0 / openFraction(I,J,K,bi,bj)
                 ENDIF !hFacC.LT.hFacMin
                 IF(hFacC(I,J,K,bi,bj).LT. 1. _d 0 .AND. 
     &              kBergBotC(I,J,bi,bj) .LT. K) THEN
                    kBergBotC(I,J,bi,bj) = K !finding deepest berg
                 ENDIF ! hFacC.LT.1
                ENDDO ! end k
               ELSE !BarrierMask.GT.0
                kBergBotC(i,j,bi,bj) = 0  
               ENDIF !BarrierMask.GT.0
              ENDDO ! end I
             ENDDO ! end J
            ENDDO ! end bi
           ENDDO !end bj

C     Partially block cells filled by icebergs, hFacS,hFacW, cell edges
C     largely copied from UPDATE_MASKS_ETC.F in CTRL pkg.
C     Smoothing might be desired, but not implemented here, see CTRL pkg 
C     for a way to do that.

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO K=1, Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           Im1=MAX(I-1,1-OLx)
           Jm1=MAX(J-1,1-OLy)
C          only update if you are in barrierMask or a neighbor,
C            and not boundary tile or face
           IF ((barrierMask(I,J,bi,bj) 
     &      + barrierMask(Im1,J,bi,bj) + barrierMask(I,Jm1,bi,bj)) 
     &      .GT. 0. _d 0 ) THEN
            IF(maskInW(I,J,bi,bj).NE.0) THEN !only if not boundary
             IF (DYG(I,J,bi,bj).EQ.0.) THEN
C              Dealing with thin walls, cribbed from CTRL pkg
               hFacW(I,J,K,bi,bj)=0.
               recip_hFacW(I,J,K,bi,bj)=0.
             ELSE
               hFacW(I,J,K,bi,bj)=
     &                   MIN(hFacC(I,J,K,bi,bj),hFacC(Im1,J,K,bi,bj))
C      &               0.5 * (hFacC(I,J,K,bi,bj) + hFacC(Im1,J,K,bi,bj))
               IF(hFacW(I,J,K,bi,bj).EQ.0.) THEN
                recip_hFacW(I,J,K,bi,bj) = 0. _d 0
              ELSE
                recip_hFacW(I,J,K,bi,bj) = 1. _d 0 / hFacW(I,J,K,bi,bj)
              ENDIF !hFacW.EQ.0
             ENDIF !GYF.EQ.0
            ENDIF !MaskInW.EQ.0
            IF(maskInS(I,J,bi,bj).NE.0) THEN 
             IF (DXG(I,J,bi,bj).EQ.0.) THEN
              hFacS(I,J,K,bi,bj)=0.
              recip_hFacS(I,J,K,bi,bj)=0.
             ELSE
              hFacS(I,J,K,bi,bj)= 
     &                 MIN(hFacC(I,J,K,bi,bj),hFacC(I,Jm1,K,bi,bj))
C      &                0.5 * (hFacC(I,J,K,bi,bj) + hFacC(I,Jm1,K,bi,bj))
              IF(hFacS(I,J,K,bi,bj).EQ.0.) THEN
                recip_hFacS(I,J,K,bi,bj) = 0. _d 0
              ELSE
                recip_hFacS(I,J,K,bi,bj) = 1. _d 0 / hFacS(I,J,K,bi,bj)
              ENDIF !hFacS.EQ.0
             ENDIF !GXG.EQ.0
            ENDIF !maskInW.NE.0
           ENDIF !barrierMask.GT.0
          ENDDO !I 
         ENDDO !J
        ENDDO !K
       ENDDO !bi
      ENDDO !bj        

#ifdef NONLIN_FRSURF
C--   Save initial geometrical hFac factor into h0Fac (fixed in time):
C     cribbed from ini_masks_etc.F
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO K=1, Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
            h0FacC(I,J,K,bi,bj) = _hFacC(I,J,K,bi,bj)
            h0FacW(I,J,K,bi,bj) = _hFacW(I,J,K,bi,bj)
            h0FacS(I,J,K,bi,bj) = _hFacS(I,J,K,bi,bj)
           ENDDO !I 
         ENDDO !J
        ENDDO !K
       ENDDO !bi
      ENDDO !bj        
#endif /* NONLIN_FRSURF */


#ifdef ALLOW_PTRACERS
C     Initialize tracer mask
      DO bi = 1,nSx
       DO bj = 1,nSy
        DO I = 1-Olx,sNx+Olx
         DO J = 1-Oly,sNy+Oly
           DO iTracer = 1,PTRACERS_num
            brg_ptracerMask (I,J,iTracer,bi,bj) = 0. _d 0
           ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C      Read ptracer mask data
       IF ( brg_ptracerMaskFile .NE. ' ' ) THEN
        CALL READ_REC_3D_RL
     &     (brg_ptracerMaskFile,precFloat64,PTRACERS_num,
     &      brg_ptracerMask,1,0,myThid)
       ENDIF
#endif /* ALLOW_PTRACERS */

      ENDIF !endif Berg OpenFrac file exists
   
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       CALL ICEBERG_DIAGNOSTICS_INIT(myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_ICEBERG */

      RETURN
      END
