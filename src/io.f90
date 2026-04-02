! =============================================================
! MODULE: io_mod
! Purpose: Namelist reading and file-based output.
!
! File formats:
!   forces.dat   –  ASCII columns: t  Cl  Ct  Cm_qc  Gamma_total
!   wake_NNNNNN.dat – ASCII columns: type  x  y  Gamma
!                     type=0 → bound vortex, type=1 → wake vortex
! =============================================================
module io_mod
  use parameters_mod
  implicit none

  integer, parameter :: FORCE_UNIT = 20
  integer, parameter :: SNAP_UNIT  = 21

contains

  ! ============================================================
  ! Read all namelist groups from the config file.
  ! Sets all variables in parameters_mod.
  ! Call init_parameters() AFTER this to derive secondary params.
  ! ============================================================
  subroutine read_config(config_file)
    character(len=*), intent(in) :: config_file

    ! Namelist groups (variable names match parameters_mod exactly)
    namelist /motion/    motion_type, h0, alpha0, phi_pitch, x_pivot, reduced_freq
    namelist /flow/      U_inf, rho, chord
    namelist /numerical/ N_panels, N_cycles, steps_per_cycle, wake_model, eps_core
    namelist /output/    output_dir, snap_interval, force_file, wake_snap_base

    integer :: unit_nml, ios

    ! ---- Defaults (in case namelist groups are missing) ------
    motion_type     = 'PLUNGE'
    h0              = 0.1_wp
    alpha0          = 0.0_wp
    phi_pitch       = 0.0_wp
    x_pivot         = 0.25_wp
    reduced_freq    = 0.5_wp
    U_inf           = 1.0_wp
    rho             = 1.0_wp
    chord           = 1.0_wp
    N_panels        = 20
    N_cycles        = 6
    steps_per_cycle = 100
    wake_model      = 'FROZEN'
    eps_core        = 1.0e-3_wp
    output_dir      = 'output'
    snap_interval   = 10
    force_file      = 'forces.dat'
    wake_snap_base  = 'wake'

    ! ---- Open namelist file -----------------------------------
    open(newunit=unit_nml, file=trim(config_file), status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') "  WARNING: Could not open config file: ", trim(config_file)
      write(*,'(A)')   "  Using built-in defaults."
      return
    end if

    ! Read each group; rewind between groups (Fortran namelist
    ! search scans forward from current file position).
    read(unit_nml, nml=motion,    iostat=ios); if (ios /= 0) write(*,*) "  INFO: /motion/ group not found"
    rewind(unit_nml)
    read(unit_nml, nml=flow,      iostat=ios); if (ios /= 0) write(*,*) "  INFO: /flow/ group not found"
    rewind(unit_nml)
    read(unit_nml, nml=numerical, iostat=ios); if (ios /= 0) write(*,*) "  INFO: /numerical/ group not found"
    rewind(unit_nml)
    read(unit_nml, nml=output,    iostat=ios); if (ios /= 0) write(*,*) "  INFO: /output/ group not found"
    close(unit_nml)

  end subroutine read_config

  ! ============================================================
  ! Open the force output file and write a header.
  ! ============================================================
  subroutine open_force_file()
    character(len=300) :: fpath
    integer :: ios

    fpath = trim(output_dir) // '/' // trim(force_file)
    open(unit=FORCE_UNIT, file=trim(fpath), status='replace', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') "  ERROR: Cannot open force file: ", trim(fpath)
      stop
    end if

    ! Header
    write(FORCE_UNIT,'(A)') '# Unsteady DVM – Force History'
    write(FORCE_UNIT,'(A,F8.4,A,F8.4,A,F8.4)') &
      '# k=', k_red, '  h0/c=', h0/chord, '  kh=', kh
    write(FORCE_UNIT,'(A)') '# Columns: t   Cl   Ct   Cm_qc   Gamma_total'
    write(FORCE_UNIT,'(A)') '# All quantities non-dimensional.'
  end subroutine open_force_file

  ! ============================================================
  subroutine write_forces(t, Cl, Ct, Cm_qc, Gamma_total)
    real(wp), intent(in) :: t, Cl, Ct, Cm_qc, Gamma_total
    write(FORCE_UNIT, '(5ES18.8)') t, Cl, Ct, Cm_qc, Gamma_total
  end subroutine write_forces

  ! ============================================================
  subroutine close_force_file()
    close(FORCE_UNIT)
  end subroutine close_force_file

  ! ============================================================
  ! Write a wake snapshot.
  ! Columns: type(0=bound,1=wake)   x   y   Gamma
  ! ============================================================
  subroutine write_wake_snapshot(step, t, &
                                  xb_arr, yb_arr, Gam_b, N_b, &
                                  xw_arr, yw_arr, Gam_w, N_w)
    integer,  intent(in) :: step, N_b, N_w
    real(wp), intent(in) :: t
    real(wp), intent(in) :: xb_arr(N_b), yb_arr(N_b), Gam_b(N_b)
    real(wp), intent(in) :: xw_arr(N_w), yw_arr(N_w), Gam_w(N_w)

    character(len=300) :: fpath
    character(len=8)   :: snum
    integer :: ios, i

    write(snum, '(I6.6)') step
    fpath = trim(output_dir) // '/' // trim(wake_snap_base) // '_' // &
            trim(snum) // '.dat'

    open(unit=SNAP_UNIT, file=trim(fpath), status='replace', iostat=ios)
    if (ios /= 0) return

    write(SNAP_UNIT, '(A,I6,A,ES12.5)') '# step=', step, '  t=', t
    write(SNAP_UNIT, '(A)') '# type  x  y  Gamma   (0=bound, 1=wake)'

    do i = 1, N_b
      write(SNAP_UNIT, '(I2,3ES18.8)') 0, xb_arr(i), yb_arr(i), Gam_b(i)
    end do
    do i = 1, N_w
      write(SNAP_UNIT, '(I2,3ES18.8)') 1, xw_arr(i), yw_arr(i), Gam_w(i)
    end do

    close(SNAP_UNIT)
  end subroutine write_wake_snapshot

end module io_mod