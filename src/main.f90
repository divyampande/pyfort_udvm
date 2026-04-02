! =============================================================
! PROGRAM: unsteady_dvm
! Purpose: Main driver for the unsteady Discrete Vortex Method.
!
! Usage:
!   ./udvm [config_file]      (default: config.nml)
!
! Execution flow:
!   1. Read namelist config → init parameters
!   2. Initialise geometry (flat plate panels) and wake storage
!   3. Build constant AIC matrix A (only once for flat plate)
!   4. Time integration loop:
!        a. Compute kinematic normal wash w(j) at collocation pts
!        b. DVM step: solve for Γ_bound, shed new wake vortex, convect
!        c. Compute forces (Cl, Ct, Cm)
!        d. Write force output; write wake snapshots at intervals
!   5. Cleanup and exit
! =============================================================
program unsteady_dvm
  use parameters_mod
  use geometry_mod
  use kinematics_mod
  use influence_mod
  use wake_mod
  use solver_mod
  use forces_mod
  use io_mod
  implicit none

  ! ---- Working arrays ----------------------------------------
  real(wp), allocatable :: A_aic(:,:)           ! (N×N) constant AIC matrix
  real(wp), allocatable :: Gamma_bound(:)        ! Bound Γ at current  step
  real(wp), allocatable :: Gamma_bound_prev(:)   ! Bound Γ at previous step
  real(wp), allocatable :: Gamma_bound_new(:)    ! Bound Γ after solve
  real(wp), allocatable :: w_arr(:)              ! Normal wash vector

  ! ---- Scalars -----------------------------------------------
  real(wp) :: Cl, Ct, Cm_qc
  real(wp) :: t
  integer  :: step

  ! ---- Config file -------------------------------------------
  character(len=256) :: config_file
  integer :: nargs

  ! ============================================================
  ! 0.  Get config file name from command line (or use default)
  ! ============================================================
  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, config_file)
  else
    config_file = 'config.nml'
  end if

  ! ============================================================
  ! 1.  Read configuration and initialise parameters
  ! ============================================================
  call read_config(trim(config_file))
  call init_parameters()
  call print_parameters()

  ! Create output directory (POSIX shell; works on Linux/macOS)
  call execute_command_line('mkdir -p ' // trim(output_dir), wait=.true.)

  ! ============================================================
  ! 2.  Initialise geometry and wake storage
  ! ============================================================
  call init_geometry()
  call init_wake()

  ! ============================================================
  ! 3.  Allocate working arrays
  ! ============================================================
  allocate(A_aic(N_panels, N_panels))
  allocate(Gamma_bound(N_panels))
  allocate(Gamma_bound_prev(N_panels))
  allocate(Gamma_bound_new(N_panels))
  allocate(w_arr(N_panels))

  ! ============================================================
  ! 4.  Build constant AIC matrix
  !     A(j,i) = v at collocation j from unit Γ at bound vortex i
  !     Flat plate geometry is fixed → compute once.
  ! ============================================================
  call build_aic_matrix(xb, yb, xc, yc, N_panels, A_aic)

  ! ============================================================
  ! 5.  Initialise circulations to zero (impulsive start)
  ! ============================================================
  Gamma_bound      = 0.0_wp
  Gamma_bound_prev = 0.0_wp

  ! ============================================================
  ! 6.  Open output files
  ! ============================================================
  call open_force_file()

  ! ============================================================
  ! 7.  Time integration
  ! ============================================================
  write(*, '(A)') "  Starting time integration ..."

  do step = 1, N_steps

    ! Current time at the START of this step
    t = real(step - 1, wp) * dt

    ! (a) Kinematic normal wash
    call compute_normalwash(t, xc, N_panels, w_arr)

    ! (b) DVM step: modifies Gamma_bound_new and wake arrays
    Gamma_bound_prev = Gamma_bound
    call solver_step(A_aic, xb, yb, xc, yc, N_panels, w_arr, &
                     x_TE_body, y_TE_body, &
                     Gamma_bound, Gamma_bound_new)
    Gamma_bound = Gamma_bound_new

    ! (c) Forces
    call compute_forces(Gamma_bound, Gamma_bound_prev, xb, N_panels, &
                        Cl, Cm_qc, Ct)

    ! (d) Write force line
    call write_forces(t, Cl, Ct, Cm_qc, sum(Gamma_bound))

    ! (e) Wake snapshots
    if (mod(step, snap_interval) == 0) then
      call write_wake_snapshot(step, t, &
             xb, yb, Gamma_bound,   N_panels, &
             x_wake, y_wake, Gamma_wake, n_wake)
    end if

    ! (f) Progress report: one line per cycle
    if (mod(step, steps_per_cycle) == 0) then
      write(*, '(A,I4,A,I4,A,F7.3,A,F8.4,A,I6,A)') &
        "    Cycle ", step / steps_per_cycle, "/", N_cycles, &
        "  t=", t, "  Cl=", Cl, &
        "  n_wake=", n_wake, " vortices"
    end if

  end do

  write(*, '(A)') "  Time integration complete."

  ! ============================================================
  ! 8.  Cleanup
  ! ============================================================
  call close_force_file()
  call cleanup_geometry()
  call cleanup_wake()

  deallocate(A_aic, Gamma_bound, Gamma_bound_prev, Gamma_bound_new, w_arr)

  write(*, '(/,A,A)') "  Results written to: ", trim(output_dir)
  write(*, '(A,/)')   "  Done."

end program unsteady_dvm