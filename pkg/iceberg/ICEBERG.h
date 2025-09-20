C $Header: /u/gcmpack/MITgcm/pkg/iceberg/ICEBERG.h,v 1.20 2015/04/22 21:33:58 bdavison Exp $
C $Name: checkpoint65m $

#ifdef ALLOW_ICEBERG

CBOP
C !ROUTINE: ICEBERG.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | ICEBERG.h
C     | o Basic header thermodnynamic iceberg ice package.
C     |   Contains all ICEBERG field declarations.
C     *==========================================================*

C==============================================================================
C     FILES AND FILE VARIABLES
C     ICEBERGmaskFile                 :: File containing iceberg mask and flag for orientation
C     ICEBERGmaskNumsFile             :: File containing arbitrary numbers for each cell containing icebergs. Links to text files containing iceberg dimensions.
C     ICEBERGnumPerCellFile           :: File containing number of icebergs per cell.
C     ICEBERGdriftFile                :: File containing mask of where effect of iceberg drift on melting will be calculated (logical)
C     ICEBERGbarrierFile              :: File containing mask for where icebergs act as physical barrier to water flow
C     ICEBERGopenFracFile             :: File containing proportion of cell volume that is open (i.e. not taken up by icebergs)
C     ICEBERGareaFile                 :: File containing total submerged iceberg surface area in each cell.
C     ICEBERGlengthFile               :: File containing lengths of all icebergs.
C     ICEBERGwidthsFile               :: File containing widths of all icebergs.
C     ICEBERGdepthsFile               :: File containing depths of all icebergs.
C     icebergMask                     :: XY Mask for iceberg cells and iceberg orientation (1 = long axis oriented east-west)
C     icebergMaskNums                 :: XY field containing numbers corresponding to each column with icebergs
C     icebergNumBergs                 :: XY field containing number of icebergs per cell
C     driftMask                       :: XY mask of where to calculate iceberg drift velocity
C     meltMask                        :: XY mask of where to calculate iceberg melt
C     barrierMask                     :: XY mask of where to make icebergs physical barrier to water flow
C     openFraction                    :: XYZ field specifying proportion of cell that is open
C     icebergArea3D                   :: XYZ field of iceberg submerged surface area
C     icebergLength                   :: XY[maxBerg] containing lengths of all icebergs.
C     icebergWidths                   :: XY[maxBerg] containing widths of all icebergs.
C     icebergDepths                   :: XY[maxBerg] containing depths of all icebergs.
C
C===============================================================================
C-   CONSTANTS SET IN data.iceberg
C     icebergRho                      :: Iceberg density (def: 917 kg m^-3 s^-1)
C     brg_iceTemp                     :: Surface temperature on the top of icefront (def: 0 degC). Interior temperature of the changes linearly from ICEFRONTthetaSurface at surface to 0 oC at the bottom
C     brg_L                           :: Latent heat of fusion (def: 334.0*10^3 J kg^-1)
C     brg_c_i                         :: heat capacity of icebergs (def: 2000 J K^-1 kg^-1)
C     brg_Cd                          :: quadratic drag coefficient (def: 0.0025)
C     icebergBGvel                    :: Constant minimum background velocity applied to iceberg faces (m s^-1)
C     lambda1                         :: Freezing point slope (def: -0.0573 degC psu^-1)
C     lambda2                         :: Freezing point offset (def: 0.0832 degC)
C     lambda3                         :: Freezing point depth (def: -7.61*10^-4)
C     brg_GamT                        :: Thermal turbulent transfer coeffcient (def: 0.022)
C     brg_GamS                        :: Salt turbulent transfer coefficient (def: 0.00062)
C     brg_c_w                         :: Heat capacity of water (def: 3974 J kg^-1 degC^-1)
C     brg_SelectDrag                  :: select how drag is computed from velocity (def: 3)
C                                        (1:n = 2, 2:n = 2, 3: n = 1 + .75*hFacC, 4: 1 + .75*(hFacC)**3
C     brg_SelectFill                  :: select how frontal area scales with hFacC (def: 3)
C                                        (1:linear, 2:quad, 3:quartic)
C     brg_DragForm                    :: form drag across iceberg (default = 0.0025)
C     brg_useInputPtracers            :: select if using ptracers input

C=============================================================================
C     FIELDS
C     kBergBotC               :: index of the bottom "berg cell" (2D), if 0 no berg drag
C     icebergHeatFlux3D       :: upward heat flux (W/m^2)
C     icebergFWFlux3D         :: upward fresh water flux (virt. salt flux) (kg/m^2/s)
C     icebergMeltRate3D       :: Melt rate (m/d)
C     icebergTendT3D          :: Temperature tendency (Kelvin/s)
C     icebergTendS3D          :: Salinity tendency (psu/s)
#ifdef ALLOW_DIAGNOSTICS
C     icebergDragU          :: iceberg stress (for diagnostics), Zonal comp. 
C                               Units are N/m^2 ;   > 0 increase top uVel
C     icebergDragV          :: iceberg stress (for diagnostics), Merid. comp. 
C                               Units are N/m^2 ;   > 0 increase top vVel
#endif /* ALLOW_DIAGNOSTICS */
C==============================================================================
C \ev
CEOP

      COMMON /ICEBERG_PARMS_I/
     &     brg_SelectDrag,
     &     brg_SelectFill,
     &     brg_useInputPtracers
      INTEGER brg_SelectDrag
      INTEGER brg_SelectFill
      INTEGER brg_useInputPtracers

      COMMON /ICEBERG_PARMS_R/
     &     icebergRho,
     &     brg_iceTemp,
     &     icebergBGvel,
     &     brg_lambda1,
     &     brg_lambda2,
     &     brg_lambda3,
     &     brg_GamT,
     &     brg_GamS,
     &     brg_c_w,
     &     brg_c_i,
     &     brg_L,
     &     brg_Cd,
     &     brg_DragForm
      _RL icebergRho
      _RL brg_iceTemp
      _RL icebergBGvel
      _RL brg_lambda1
      _RL brg_lambda2
      _RL brg_lambda3
      _RL brg_GamT
      _RL brg_GamS
      _RL brg_c_w
      _RL brg_c_i
      _RL brg_L
      _RL brg_Cd
      _RL brg_DragForm

      COMMON /ICEBERG_FIELDS_I/ kBergBotC
      INTEGER kBergBotC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      COMMON /ICEBERG_FIELDS_RL/
     &     icebergHeatFlux3D,
     &     icebergFWFlux3D,
     &     icebergMeltRate3D,
     &     icebergTendT3D,
     &     icebergTendS3D,
     &     openFraction, 
     &     driftMask, 
     &     meltMask,
     &     icebergMask,
     &     barrierMask, 
     &     icebergMaskNums, 
     &     icebergNumBergs,
     &     icebergArea3D,
     &     icebergLength,
     &     icebergWidths,
     &     icebergDepths
      _RL icebergHeatFlux3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  
      _RL icebergFWFlux3D (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  
      _RL icebergMeltRate3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)  
      _RL icebergTendT3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL icebergTendS3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL openFraction(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL driftMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL meltMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL icebergMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL barrierMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL icebergMaskNums(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL icebergNumBergs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL icebergArea3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL icebergLength(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       maxBergCt,nSx,nSy) 
      _RL icebergWidths(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       maxBergCt,nSx,nSy)
      _RL icebergDepths(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       maxBergCt,nSx,nSy)

      COMMON /ICEBERG_PARM_C/
     &     ICEBERGmaskFile,
     &     ICEBERGmaskNumsFile,
     &     ICEBERGnumPerCellFile,
     &     ICEBERGdriftFile,
     &     ICEBERGmeltFile,
     &     ICEBERGopenFracFile,
     &     ICEBERGbarrierFile,
     &     ICEBERGareaFile,
     &     ICEBERGlengthFile,
     &     ICEBERGwidthsFile,
     &     ICEBERGdepthsFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGmaskFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGmaskNumsFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGnumPerCellFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGmeltFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGdriftFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGopenFracFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGbarrierFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGareaFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGlengthFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGwidthsFile
      CHARACTER*(MAX_LEN_FNAM) ICEBERGdepthsFile

#ifdef ALLOW_DIAGNOSTICS
      COMMON /ICEBERG_DIAG_DRAG/ icebergDragU, icebergDragV
      _RS icebergDragU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS icebergDragV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_PTRACERS
C Parameters relating to PTRACERS

      COMMON /ICEBERG_PTRACERS_RL/
     &     brg_ptr_addMass3D,
     &     brg_ptracerMask
      _RL brg_ptr_addMass3d
     &     (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy,PTRACERS_num)
      _RL brg_ptracerMask
     &     (1-Olx:sNx+Olx,1-Oly:sNy+Oly,PTRACERS_num,nSx,nSy)
   
      COMMON /ICEBERG_PTRACERS_FILES/
     &     brg_ptracerMaskFile
      CHARACTER*(512) 
     &     brg_ptracerMaskFile

#endif /* ALLOW_PTRACERS */

#endif /* ALLOW_ICEBERG */
