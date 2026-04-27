!===========================================================================
! Ice.f90
! Module and subroutines for extending pyMGMR3D to include in-ice showers.
!
! Physics summary:
!   When a cosmic-ray air shower enters ice it creates a secondary cascade.
!   This module provides:
!     (1) Tabulated ice density/refractivity profile from current_ice.dat
!     (2) In-ice shower current arrays (Ix_ice, Iy_ice, IQ_ice)
!     (3) xi_ice = n(z)-1 for retarded-time calculations in CurrField
!     (4) Cumulative slant-depth IcePenDepth, continuing from X_surface
!
! Coordinate convention (same sign convention as the rest of MGMR3D):
!   zeta > 0  : height above ground along shower axis  (atmosphere)
!   zeta < 0  : depth below ground along shower axis   (ice)
!   Ice index i=1 -> zeta = -1 * IceHei_step  (just below surface)
!
! Call order (mirrors the atmosphere setup):
!   AssignDim()         -> call AssignIceDim(...)
!   Initialize_shower() -> call Initialize_Ice_Layer(PenDepth(0), Cos_Zenith)
!                       -> call Compute_Ice_Shower_Currents(...)
!   CurrField()         -> second loop using Ix_ice, IceRefrac, etc.
!   DAssignArr()        -> call DAssignIceArr()
!===========================================================================

module IceArrays
    use constants, only : dp
    implicit none
    save

    ! --- Dimensions (set by AssignIceDim) ---
    integer  :: Ice_dim      = 0        ! number of ice depth steps; 0 = disabled
    real(dp) :: IceHei_step  = 5.0_dp   ! step size along shower axis in ice [m]
    integer  :: ice_model_id = 0        ! 0=none 1=constant 2=S.Pole 3=Greenland

    ! --- Profile arrays: index 0 = ice surface, i = i*IceHei_step below ---
    real(dp), allocatable :: IceDepth    (:)   ! depth below surface [m]
    real(dp), allocatable :: IceDensity  (:)   ! density [g/cm3]
    real(dp), allocatable :: IceRefrac   (:)   ! xi_ice = n(z)-1  (0.31 at surface, 0.78 deep)
    real(dp), allocatable :: dIceRefrac  (:)   ! d(xi_ice)/dz [1/m], for retarded distance
    real(dp), allocatable :: IcePenDepth (:)   ! unified slant-depth from top-of-atm [g/cm2]

    ! --- In-ice shower current arrays (same role as Ix,Iy,IQ in BigArrays) ---
    real(dp), allocatable :: Ix_ice      (:)   ! transverse current x
    real(dp), allocatable :: Iy_ice      (:)   ! transverse current y
    real(dp), allocatable :: IQ_ice      (:)   ! charge excess current (Askaryan)
    real(dp), allocatable :: Ix_int_ice  (:)   ! integrated Ix (moving-dipole term)
    real(dp), allocatable :: Iy_int_ice  (:)   ! integrated Iy (moving-dipole term)
    real(dp), allocatable :: alpha_tr_ice(:)   ! pancake thickness index

    ! --- Observer position ---
    real(dp) :: z_observer = 0.0_dp    ! antenna depth [m]; 0=surface, <0=in ice

    ! --- Derived scalars ---
    real(dp) :: X_surface     = 0.0_dp    ! total atmospheric slant-depth at surface [g/cm2]
    real(dp) :: ice_depth_max = 500.0_dp  ! max ice depth to simulate [m]; set by SetParams

end module IceArrays


!===========================================================================
subroutine AssignIceDim(Ice_dim_in, IceHei_step_in)
!   Allocate all ice arrays.
!   Called from AssignDim() in MGMR3D_BA-v4.f90.
!
!   Ice_dim_in    = NINT(ice_depth_max / IceHei_step_in)
!   IceHei_step_in = step size in metres (recommend = AtmHei_step)
!
    use IceArrays
    use constants, only : dp
    implicit none
    integer,  intent(in) :: Ice_dim_in
    real(dp), intent(in) :: IceHei_step_in

    Ice_dim     = Ice_dim_in
    IceHei_step = IceHei_step_in

    if (Ice_dim <= 0) return

    allocate( IceDepth    (0:Ice_dim) )
    allocate( IceDensity  (0:Ice_dim) )
    allocate( IceRefrac   (0:Ice_dim) )
    allocate( dIceRefrac  (0:Ice_dim) )
    allocate( IcePenDepth (0:Ice_dim) )
    allocate( Ix_ice      (0:Ice_dim) )
    allocate( Iy_ice      (0:Ice_dim) )
    allocate( IQ_ice      (0:Ice_dim) )
    allocate( Ix_int_ice  (0:Ice_dim) )
    allocate( Iy_int_ice  (0:Ice_dim) )
    allocate( alpha_tr_ice(0:Ice_dim) )

    ! Sensible defaults (overwritten by Initialize_Ice_Layer)
    IceDepth     = 0.0_dp
    IceDensity   = 0.917_dp   ! bulk ice
    IceRefrac    = 0.788_dp   ! n=1.788 for rho=0.917
    dIceRefrac   = 0.0_dp
    IcePenDepth  = 0.0_dp
    Ix_ice       = 0.0_dp
    Iy_ice       = 0.0_dp
    IQ_ice       = 0.0_dp
    Ix_int_ice   = 0.0_dp
    Iy_int_ice   = 0.0_dp
    alpha_tr_ice = 1.0_dp

end subroutine AssignIceDim


!===========================================================================
subroutine DAssignIceArr()
!   Deallocate ice arrays.
!   Called from DAssignArr() in MGMR3D_BA-v4.f90.
!
    use IceArrays
    implicit none
    if (Ice_dim <= 0) return
    deallocate( IceDepth, IceDensity, IceRefrac, dIceRefrac, IcePenDepth )
    deallocate( Ix_ice, Iy_ice, IQ_ice, Ix_int_ice, Iy_int_ice, alpha_tr_ice )
end subroutine DAssignIceArr


!===========================================================================
subroutine Initialize_Ice_Layer(X_surface_in, Cos_Zenith)
!   Read current_ice.dat (written by ice_models.py) and populate IceArrays.
!   MUST be called AFTER Initialize_shower so that PenDepth(0) is available.
!
!   current_ice.dat column format:
!     Lines starting with '!' are comments and are skipped.
!     depth[m]  density[g/cm3]  xi=n-1  X_cumul_vertical[g/cm2]
!
    use IceArrays
    use constants, only : dp
    implicit none
    real(dp), intent(in) :: X_surface_in   ! = PenDepth(0) from Initialize_shower
    real(dp), intent(in) :: Cos_Zenith     ! cos(zenith angle of shower)

    integer  :: i, ios, n_read
    real(dp) :: depth_in, rho_in, xi_in, Xcumul_in
    character(len=200) :: line_buf

    X_surface = X_surface_in

    if (Ice_dim <= 0) then
        write(2,*) 'Initialize_Ice_Layer: ice disabled (Ice_dim=0), skipping.'
        return
    endif

    open(unit=8, file='current_ice.dat', status='old', iostat=ios)
    if (ios /= 0) then
        write(2,*) 'FATAL: cannot open current_ice.dat'
        stop 'current_ice.dat not found - run ice_models.py first'
    endif

    n_read = 0
    do
        read(8, '(A)', iostat=ios) line_buf
        if (ios /= 0) exit
        if (len_trim(line_buf) == 0) cycle
        if (line_buf(1:1) == '!') cycle
        read(line_buf, *, iostat=ios) depth_in, rho_in, xi_in, Xcumul_in
        if (ios /= 0) exit
        i = n_read
        if (i > Ice_dim) then
            write(2,'(A,I6)') 'WARNING: current_ice.dat longer than Ice_dim=', Ice_dim
            exit
        endif
        IceDepth(i)    = depth_in
        IceDensity(i)  = rho_in
        IceRefrac(i)   = xi_in
        ! Unified slant-depth = atmospheric slant-depth at ice surface
        !                      + vertical ice grammage / cos(zenith)
        IcePenDepth(i) = X_surface + Xcumul_in / Cos_Zenith
        n_read = n_read + 1
    enddo
    close(unit=8)

    if (n_read < 2) then
        write(2,*) 'FATAL: current_ice.dat has fewer than 2 valid data rows.'
        stop 'current_ice.dat too short'
    endif

    if (n_read - 1 < Ice_dim) then
        write(2,'(A,I5,A,I5)') 'WARNING: current_ice.dat has', n_read, &
            ' rows, trimming Ice_dim from', Ice_dim
        Ice_dim = n_read - 1
    endif

    ! Central-difference derivative of refractivity (used in calD in CurrField)
    dIceRefrac(0) = (IceRefrac(1) - IceRefrac(0)) / IceHei_step
    do i = 1, Ice_dim-1
        dIceRefrac(i) = (IceRefrac(i+1) - IceRefrac(i-1)) / (2.0_dp * IceHei_step)
    enddo
    dIceRefrac(Ice_dim) = dIceRefrac(Ice_dim-1)

    write(2,'(/,A)') '=== Ice Layer Initialized ==='
    write(2,'(A,I5,A,F7.1,A)') &
        '  Steps:', Ice_dim+1, ',  step size:', IceHei_step, ' [m]'
    write(2,'(A,F8.4,A,F8.6)') &
        '  Surface: rho=', IceDensity(0), ' g/cm3,  xi=', IceRefrac(0)
    write(2,'(A,F8.4,A,F8.6)') &
        '  Deep:    rho=', IceDensity(Ice_dim), ' g/cm3,  xi=', IceRefrac(Ice_dim)
    write(2,'(A,F10.2,A,F10.2,A)') &
        '  X_surface=', X_surface, &
        ' g/cm2,  X_ice_max=', IcePenDepth(Ice_dim), ' g/cm2'
    write(2,'(A,F6.2,A)') &
        '  Cherenkov angle in deep ice: ', &
        acos(1.0_dp/(1.0_dp+IceRefrac(Ice_dim)))*180.0_dp/acos(-1.0_dp), ' deg'
    write(2,'(A,/)') '============================='

end subroutine Initialize_Ice_Layer


!===========================================================================
subroutine Compute_Ice_Shower_Currents( &
        Energy_sh2, X_max2, X_02, lamx2, R_02, L_02, RL_param, &
        J0Q, u0, F_over_beta, &
        Xb_0, Xc_0, a_ChX, &
        Force0_ice, alpha_frc0_ice)
!   Compute Ix_ice, Iy_ice, IQ_ice for the in-ice secondary shower.
!   Called from Initialize_shower in MGMR3D_shower-v5.f90 immediately
!   after the atmospheric currents are computed.
!
!   Parameters mirror those used in Initialize_shower for the air shower:
!     Energy_sh2      : energy normalization of the in-ice shower
!     X_max2          : Xmax in UNIFIED slant-depth [g/cm2]
!                       (= X_surface + X_max_ice_vertical/cos(zenith))
!                       computed by ice_models.py ice_depth_to_unified_slantdepth()
!     X_02            : shower start in unified slant-depth [g/cm2]
!                       = X_surface (shower starts at ice entry)
!     lamx2/R_02/L_02 : longitudinal profile shape parameters (same as air)
!     RL_param        : .true. = R-L profile, .false. = G-H profile
!     J0Q             : charge excess normalization constant
!     u0, F_over_beta : drift velocity parameters
!     Xb_0, Xc_0, a_ChX : charge excess shape parameters
!     Force0_ice      : geomagnetic force magnitude in ice [keV/m]
!     alpha_frc0_ice  : direction of geomagnetic force in ice [rad]
!
!   Key physics in ice vs. air:
!     - BoF uses sqrt(rho_Xmx/rho_i): in deep ice this ratio ~ 1 (correct:
!       suppresses geomagnetic emission relative to charge excess)
!     - IQ scales as sqrt(rho_i/0.06) * rho_Xmx/0.06
!
    use IceArrays
    use constants, only : dp, pi
    implicit none

    real(dp), intent(in) :: Energy_sh2, X_max2, X_02, lamx2, R_02, L_02
    logical,  intent(in) :: RL_param
    real(dp), intent(in) :: J0Q, u0, F_over_beta
    real(dp), intent(in) :: Xb_0, Xc_0, a_ChX
    real(dp), intent(in) :: Force0_ice, alpha_frc0_ice

    integer  :: i, iXmx_ice
    real(dp) :: X_rh, NPart2, BoF, ux, uy, s_vel
    real(dp) :: Force_x_ice, Force_y_ice
    real(dp) :: rho_Xmx, rho_i, b_ce
    real(dp) :: dX_step, decay_fac
    real(dp), parameter :: X_EM = 36.7_dp   ! radiation length in ice [g/cm2]

    if (Ice_dim <= 0)        return
    if (Energy_sh2 <= 0.0_dp) return

    Force_x_ice = Force0_ice * cos(alpha_frc0_ice)
    Force_y_ice = Force0_ice * sin(alpha_frc0_ice)

    ! Find ice index closest to Xmax for the density ratio in BoF
    iXmx_ice = 0
    do i = 0, Ice_dim
        if (IcePenDepth(i) <= X_max2) iXmx_ice = i
    enddo
    rho_Xmx = IceDensity(iXmx_ice)

    Ix_ice(:) = 0.0_dp
    Iy_ice(:) = 0.0_dp
    IQ_ice(:) = 0.0_dp

    do i = 0, Ice_dim
        X_rh = IcePenDepth(i)
        if (X_rh <= X_02) cycle      ! shower not started at this depth
        if (X_rh > X_02 + 6.0_dp * lamx2 .and. X_rh > X_max2) cycle   ! shower over

        ! ---- Longitudinal profile (identical formula to air shower) ----
        NPart2 = 0.0_dp
        if (RL_param) then
            if ((1.0_dp - R_02*(X_max2-X_rh)/L_02) > 0.0_dp) then
                NPart2 = Energy_sh2 &
                       * (1.0_dp - R_02*(X_max2-X_rh)/L_02)**(1.0_dp/(R_02*R_02)) &
                       * exp((X_max2-X_rh)/(L_02*R_02))
            endif
        else
            ! Gaisser-Hillas
            if ((X_rh-X_02) > 0.0_dp .and. (X_max2-X_02) > 0.0_dp) then
                NPart2 = Energy_sh2 &
                       * ((X_rh-X_02)/(X_max2-X_02))**((X_max2-X_02)/lamx2) &
                       * exp((X_max2-X_rh)/lamx2)
                if (NPart2 < 0.0_dp) NPart2 = 0.0_dp
            endif
        endif
        if (NPart2 <= 0.0_dp) cycle

        ! ---- Drift velocity in ice ----
        rho_i = IceDensity(i)
        BoF   = 4.0_dp / ((X_max2 - Xb_0)/(X_rh - Xb_0) + 3.0_dp) &
              * sqrt(rho_Xmx / rho_i) / F_over_beta

        ux    = Force_x_ice * BoF
        uy    = Force_y_ice * BoF
        s_vel = u0 / sqrt(u0*u0 + ux*ux + uy*uy)   ! saturation factor

        Ix_ice(i) = NPart2 * ux * s_vel
        Iy_ice(i) = NPart2 * uy * s_vel

        ! ---- Charge excess (Askaryan term) ----
        b_ce      = (X_rh - X_max2) / (500.0_dp + X_rh - Xc_0)
        IQ_ice(i) = J0Q * NPart2 &
                  * ((b_ce/a_ChX) + 1.0_dp) &
                  * (1.0_dp - exp(-5.0_dp*((X_rh-Xc_0)/(X_max2-Xc_0)))) &
                  * sqrt(rho_i/0.06_dp) * rho_Xmx/0.06_dp

        alpha_tr_ice(i) = 1.0_dp   ! pancake width parameter; extendable later

    enddo

    ! ---- Moving-dipole integrated currents (same decay logic as air shower) ----
    Ix_int_ice(Ice_dim) = 0.0_dp
    Iy_int_ice(Ice_dim) = 0.0_dp
    do i = Ice_dim-1, 0, -1
        ! grammage of one step: rho[g/cm3] * IceHei_step[m] * 100[cm/m]
        dX_step   = IceDensity(i) * IceHei_step * 100.0_dp
        decay_fac = exp(dX_step / X_EM)
        Ix_int_ice(i) = Ix_int_ice(i+1)*decay_fac &
                      - Ix_ice(i)*IceHei_step*(decay_fac - 1.0_dp)*X_EM/dX_step
        Iy_int_ice(i) = Iy_int_ice(i+1)*decay_fac &
                      - Iy_ice(i)*IceHei_step*(decay_fac - 1.0_dp)*X_EM/dX_step
    enddo

    write(2,'(A)') 'Compute_Ice_Shower_Currents:'
    write(2,'(A,E12.4,A,F8.2,A,F8.2,A)') &
        '  E_sh2=', Energy_sh2, &
        ',  X_max2=', X_max2, ' g/cm2,  X_02=', X_02, ' g/cm2'
    write(2,'(A,F7.4,A,E12.4,A,E12.4)') &
        '  rho@Xmax_ice=', rho_Xmx, &
        ',  max|IQ_ice|=', maxval(abs(IQ_ice)), &
        ',  max|Ix_ice|=', maxval(abs(Ix_ice))

end subroutine Compute_Ice_Shower_Currents
