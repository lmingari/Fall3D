!***************************************************************
!*
!*    Module for CAL operations
!*
!***************************************************************
MODULE CAL_nc
  use KindType, ONLY : ip,rp
  use InpOut,   ONLY : ludat
  implicit none
  save
  !
  !*** CALMET Headers (NOT KIND-TYPE DEPENDENT FOR CALMET)
  !
  logical :: lcalgrd
  !
  integer(4), parameter :: mxncom = 1000
  !
  integer(4) :: nssta,nusta,npsta,irtype
  integer(4) :: icbyr,icbmo,icbdy,icbhr,icbsec
  integer(4) :: iceyr,icemo,icedy,icehr,icesec
  integer(4) :: nzc,irlg,iwfcod,nowsta,nlu,ibtz,iwat1,iwat2,iutmzn,ncom,icom
  integer(4) :: nxx,nyy,nzl,datec
  integer(4) :: idatc0,idatc1,isec0,isec1
  integer(4) :: idum = 0
  !
  real(4) :: dgrid4,xorigr4,yorigr4,feast,fnorth,rnlat0,relon0,xlat1,xlat2
  !
  character(len=132) :: comment(mxncom)
  character(len=64)  :: datamod
  character(len=16)  :: dataset,dataver
  character(len=12)  :: daten
  character(len=8 )  :: pmap,datum,axtz
  character(len=4 )  :: utmhem
  character(len=14)  :: ndatcb,ndatce,ndatc0,ndatc1
  character(len=8)   :: clabdu,clab1,clab10,clabtk
  character(len=8)   :: clabu,clabv,clabw,clabt,clabsc,clabus,clabzi,clabrl
  !
  !*** Rest of (normal) variables
  !
  integer(ip) :: ibyr_CAL,ibmo_CAL,ibdy_CAL,ibhr_CAL,ibsec_CAL
  integer(ip) :: nx_CAL,ny_CAL,nz_CAL,nt_CAL,it_CAL
  real(rp)    :: dt_CAL
  real(rp)    :: xo_CAL,yo_CAL,dx_CAL,dy_CAL
  !
  !*** Interpolation stuff
  !
  integer(ip) :: nelem_CAL
  integer(ip) :: npoin_DAT
  !
  integer(ip), allocatable :: lnods_CAL(:,:)
  integer(ip), allocatable :: lelpo_CAL(:)
  !
  real   (rp), allocatable :: s_po_CAL(:)
  real   (rp), allocatable :: t_po_CAL(:)
  !
  !*** Allocatable variables
  !
  real(rp), allocatable :: time_CAL(:)
  real(rp), allocatable :: timesec_CAL(:)
  real(rp), allocatable :: xutm_CAL(:,:)
  real(rp), allocatable :: yutm_CAL(:,:)
  !
  integer(4),allocatable :: LDU_CAL(:,:)
  real(4), allocatable   :: work4(:)       ! real4 work arrary. Allocated to (nx*ny*nz)
  !
  real(rp), allocatable :: work8(:)
  real(rp), allocatable :: Z0_CAL (:,:)
  real(rp), allocatable :: topg_CAL (:,:)
  !
  real(rp), allocatable :: z_CAL   (:,:,:)
  real(rp), allocatable :: UST_CAL (:,:,:)
  real(rp), allocatable :: PBL_CAL (:,:,:)
  real(rp), allocatable :: L_CAL   (:,:,:)
  real(rp), allocatable :: TSFC_CAL(:,:,:)
  !
  real(rp), allocatable :: u_CAL   (:,:,:,:)
  real(rp), allocatable :: v_CAL   (:,:,:,:)
  real(rp), allocatable :: w_CAL   (:,:,:,:)
  real(rp), allocatable :: T_CAL   (:,:,:,:)
  !
CONTAINS
  !
  !
  !
  subroutine CAL_getmem
    !***************************************************************
    !*
    !*   Gets memory for CALMET variables
    !*
    !***************************************************************
    implicit none
    !
    allocate(time_CAL(nt_CAL))
    allocate(timesec_CAL(nt_CAL))
    !
    !*** Variables
    !
    allocate(work4(nx_CAL*ny_CAL*nz_CAL))
    allocate(work8(nx_CAL*ny_CAL*nz_CAL))
    allocate(xutm_CAL (nx_CAL,ny_CAL))
    allocate(yutm_CAL (nx_CAL,ny_CAL))
    !
    allocate(Z0_CAL   (nx_CAL,ny_CAL       ))
    allocate(LDU_CAL  (nx_CAL,ny_CAL       ))
    allocate(topg_CAL (nx_CAL,ny_CAL       ))
    !
    allocate(z_CAL    (nx_CAL,ny_CAL,nz_CAL))
    allocate(UST_CAL  (nx_CAL,ny_CAL,nt_CAL))
    allocate(PBL_CAL  (nx_CAL,ny_CAL,nt_CAL))
    allocate(L_CAL    (nx_CAL,ny_CAL,nt_CAL))
    allocate(TSFC_CAL (nx_CAL,ny_CAL,nt_CAL))
    !
    allocate(u_CAL    (nx_CAL,ny_CAL,nz_CAL,nt_CAL))
    allocate(v_CAL    (nx_CAL,ny_CAL,nz_CAL,nt_CAL))
    allocate(w_CAL    (nx_CAL,ny_CAL,nz_CAL,nt_CAL))
    allocate(T_CAL    (nx_CAL,ny_CAL,nz_CAL,nt_CAL))
    !
    !*** Interpolation
    !
    allocate(lnods_CAL(4,nelem_CAL))
    allocate(lelpo_CAL(npoin_DAT))
    allocate(s_po_CAL(npoin_DAT))
    allocate(t_po_CAL(npoin_DAT))
    !
  end subroutine CAL_getmem
  !
END MODULE CAL_nc
