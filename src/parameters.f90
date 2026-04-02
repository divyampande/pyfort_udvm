! =============================================================
! MODULE: parameters_mod
! Purpose: All physical and numerical parameters for the simulation.
!          Holds module-level variables that are initialised by
!          reading the namelist (via io_mod) and then calling
!          init_parameters() to derive secondary quantities.
!
! SIGN CONVENTIONS (read this before touching anything):
!   - Plunge h: positive DOWNWARD  (Theodorsen / Katz&Plotkin convention)
!   - Pitch α:  positive NOSE-UP
!   - Vortex Γ: positive COUNTER-CLOCKWISE
!   - Normal n̂: (0,1) pointing UP from the flat plate
!   - Kinematic wash w_j = hdot + U∞α + αdot*(xc_j - xpivot)
!       (positive = flow entering from below, requiring upwash from bound Γ)
!   - Boundary condition: ΣA(j,i)Γ_i + wake = -w_j
!   - Biot-Savart: u = +Γ(yp-yv)/(2π r²),  v = -Γ(xp-xv)/(2π r²)
!       A CCW vortex at origin induces DOWNWARD flow at (1,0). ✓
! =============================================================
module parameters_mod
  use iso_fortran_env, only: real64
  implicit none

  integer,  parameter :: wp     = real64
  real(wp), parameter :: PI     = 4.0_wp * atan(1.0_wp)
  real(wp), parameter :: TWO_PI = 2.0_wp * PI

  ! ---- Flow parameters ----------------------------------------
  real(wp) :: U_inf    ! Freestream velocity
  real(wp) :: rho      ! Fluid density
  real(wp) :: chord    ! Chord length
  real(wp) :: b_semi   ! Semi-chord  b = c/2  (Theodorsen convention)

  ! ---- Motion parameters --------------------------------------
  character(len=16) :: motion_type  ! 'PLUNGE', 'PITCH', 'COMBINED'
  real(wp) :: h0          ! Plunge amplitude [same units as chord]
  real(wp) :: alpha0      ! Pitch amplitude [radians]
  real(wp) :: phi_pitch   ! Phase of pitch relative to plunge [radians]
  real(wp) :: x_pivot     ! Pivot point as fraction of chord from LE
  real(wp) :: reduced_freq ! k = omega*c/(2*U_inf)   (input from namelist)
  real(wp) :: k_red        ! Same as reduced_freq (internal alias for clarity)
  real(wp) :: omega        ! Angular frequency [rad/s]
  real(wp) :: kh           ! Non-dim plunge velocity = k * (h0/chord)
                           !   = omega*h0 / (2*U_inf)
                           !   kh characterises the wake pattern and forces

  ! ---- Numerical parameters -----------------------------------
  integer            :: N_panels        ! Number of bound vortex panels
  integer            :: N_cycles        ! Number of oscillation cycles
  integer            :: steps_per_cycle ! Steps per cycle
  integer            :: N_steps         ! Total steps (derived)
  real(wp)           :: dt              ! Time step (derived)
  character(len=16)  :: wake_model      ! 'FROZEN' or 'ROLLUP'
  real(wp)           :: eps_core        ! Vortex core regularisation radius

  ! ---- Output parameters --------------------------------------
  character(len=128) :: output_dir
  integer            :: snap_interval
  character(len=64)  :: force_file
  character(len=64)  :: wake_snap_base

contains

  ! ------------------------------------------------------------
  ! Derive secondary parameters after the namelist has been read.
  ! Call this AFTER read_config() in io_mod.
  ! ------------------------------------------------------------
  subroutine init_parameters()
    b_semi  = chord / 2.0_wp

    ! Reduced frequency: k = omega*c/(2*U_inf)  =>  omega = 2*k*U_inf/c
    k_red   = reduced_freq
    omega   = 2.0_wp * k_red * U_inf / chord

    ! Non-dimensional plunge velocity amplitude
    kh      = k_red * (h0 / chord)

    ! Time discretisation: one oscillation period T = 2π/omega
    N_steps = N_cycles * steps_per_cycle
    dt      = TWO_PI / omega / real(steps_per_cycle, wp)
  end subroutine init_parameters

  ! ------------------------------------------------------------
  ! Pretty-print all relevant parameters to stdout.
  ! ------------------------------------------------------------
  subroutine print_parameters()
    write(*,'(/,A)') "  =============================================="
    write(*,'(A)')   "    Unsteady DVM  –  Simulation Parameters"
    write(*,'(A)')   "  =============================================="
    write(*,'(A,A)')    "    Motion type      : ", trim(motion_type)
    write(*,'(A,A)')    "    Wake model       : ", trim(wake_model)
    write(*,'(A,F10.5)') "    Reduced freq  k  : ", k_red
    write(*,'(A,F10.5)') "    Plunge amp h0/c  : ", h0 / chord
    write(*,'(A,F10.5)') "    kh = k*(h0/c)    : ", kh
    write(*,'(A,F10.5)') "    omega            : ", omega
    write(*,'(A,F10.5)') "    Period T         : ", TWO_PI / omega
    write(*,'(A,F10.5)') "    dt               : ", dt
    write(*,'(A,I8)')    "    N_panels         : ", N_panels
    write(*,'(A,I8)')    "    N_steps          : ", N_steps
    write(*,'(A,F10.5)') "    eps_core         : ", eps_core
    write(*,'(A,/)')   "  =============================================="
  end subroutine print_parameters

end module parameters_mod