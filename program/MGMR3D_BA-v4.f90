! V17: Include radial dependence of pancake thickness, radial loop at the deepest level
! V18: Include radial dependence of pancake thickness, radial loop at a rather high level
!
  module BigArrays
    use constants, only : dp
  implicit none
  integer, save :: ObsDist_dim
  integer, save :: tTrace_dim_b, tTrace_dim_o, tTrace_Offset
  real(dp), save :: tTrace_Step, nuTrace_step
  integer, save :: nuTrace_dim
!  integer, save :: i_nu_ini,i_nu_max,Padding
!  integer, parameter :: AtmHei_dim=2000
  integer, save :: AtmHei_dim
  integer, parameter :: LamInt_dim=7000
  integer, parameter :: CoreDist_Dim=20  ! number of line-showers per shower
  real(dp), save ::  ObsDist_Step,  ObsDist,  AtmHei_step, IndRefOrho, TopAtmExpon
    real*8, allocatable, save :: Z_ch(:),T_ch(:)
    !real*8, allocatable, save :: Ex(:),AxD(:),Ey(:),AyD(:),Ar(:),Erh(:)
    complex*16, allocatable, save :: Ex_nu(:),Ey_nu(:),Er_nu(:)
    complex*16, allocatable, save :: CEx(:),CEy(:),CEr(:)
    !complex*16, allocatable, save :: Ex_nub(:,:,:),AxD_nub(:,:,:),Ey_nub(:,:,:),AyD_nub(:,:,:),Ar_nub(:,:,:),Erh_nub(:,:,:)
    complex*16, allocatable, save :: Ex_nu_dwn(:,:),Ey_nu_dwn(:,:),Er_nu_dwn(:,:),filt(:)
    real(dp), allocatable, save :: Ex_tb(:,:,:) ,AxD_tb(:,:,:),Ey_tb(:,:,:),AyD_tb(:,:,:),Ar_tb(:,:,:),Erh_tb(:,:,:)
    real(dp), allocatable, save :: Ex_spld(:,:,:) ,AxD_spld(:,:,:),Ey_spld(:,:,:),AyD_spld(:,:,:),Ar_spld(:,:,:),Erh_spld(:,:,:)
    real(dp), allocatable, save :: t_tb(:), t_to(:), Line2Obsrv(:)
    real(dp), allocatable, save :: Ex_to(:,:) ,Ey_to(:,:),Er_to(:,:),ObsPlsTime(:)
    real*8, allocatable, save :: L_I(:),JL_I(:),dL_I(:),zeta_I(:,:)
    real*8, allocatable, save :: T_obs(:),lambda(:),PenDepth(:),xi(:),dxi(:),ddxi(:)
    real*8, allocatable, save :: Ix(:),Iy(:),IQ(:),alpha_tr(:),Ix_int(:),Iy_Int(:)
    real(dp), save :: PancakeThi(0:CoreDist_Dim), CoreDist_A(0:CoreDist_Dim), Line2Core(0:CoreDist_Dim*2+3)
    integer, save :: LamGrid_Nrt, LamGrid_max
  contains
  subroutine AssignDim()
    allocate(Z_ch(ObsDist_dim))
    allocate(T_ch(ObsDist_dim))
!    allocate(Ex(1:tTrace_dim)) ;  allocate(AxD(1:tTrace_dim))
!    allocate(Ey(1:tTrace_dim))
!    allocate(AyD(1:tTrace_dim_o)) ;  allocate(Ar(1:tTrace_dim_o)) ;  allocate(Erh(1:tTrace_dim_o))
    allocate(t_tb(1:tTrace_dim_b)) ;     allocate(t_to(1:tTrace_dim_o))
    allocate(L_I(0:LamInt_dim)) ; allocate(JL_I(0:LamInt_dim)) ; allocate(dL_I(0:LamInt_dim))
    allocate(zeta_I(2,0:LamInt_dim))
    allocate(T_obs(0:AtmHei_dim)) ; allocate(lambda(0:AtmHei_dim))
    allocate(PenDepth(0:AtmHei_dim) , Line2Obsrv(1:ObsDist_dim))
    allocate(xi(0:AtmHei_dim)) ; allocate(dxi(0:AtmHei_dim)) ; allocate(ddxi(0:AtmHei_dim))
    allocate(Ix(0:AtmHei_dim)) ; allocate(Iy(0:AtmHei_dim)) ; allocate(IQ(0:AtmHei_dim))
    allocate(alpha_tr(0:AtmHei_dim)) ; allocate(Ix_int(0:AtmHei_dim)) ; allocate(Iy_Int(0:AtmHei_dim))
  end subroutine AssignDim
  subroutine DAssignArr()
    deallocate(Ex_nu, Ey_nu, Er_nu)
    deallocate(CEx, CEy, CEr) ; deallocate(ObsPlsTime)
    deallocate(Ex_nu_dwn, Ey_nu_dwn, Er_nu_dwn, filt)
! DwnSample_E:      deallocate(t_to, Ex_to, Ey_to, Er_to)
! REPE, pacakeint
!    deallocate(L_I, JL_I) ; deallocate(dL_I) !; deallocate(End_I) ;
!    deallocate(zeta_I) ;  deallocate(PenDepth , Line2Obsrv)
!    deallocate(T_obs) ; deallocate(lambda)
!    deallocate(xi) ; deallocate(dxi) ; deallocate(ddxi)
!    deallocate(Ix) ; deallocate(Iy) ; deallocate(IQ)
!    deallocate(alpha_tr) ; deallocate(Ix_int) ; deallocate(Iy_Int)
! REPE, radialint
!    deallocate(Z_ch, T_ch) !;
!    deallocate(t_tb, Ex_tb, AxD_tb, Ey_tb, AyD_tb, Ar_tb, Erh_tb)
!    deallocate(Ex_spld, AxD_spld, Ey_spld, AyD_spld, Ar_spld, Erh_spld)
!
  end subroutine DAssignArr
  end module BigArrays
! ==========================================================
