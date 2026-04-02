! =============================================================
! MODULE: kinematics_mod
! Purpose: Defines airfoil motion and computes the kinematic
!          normal wash (boundary condition RHS) at each panel's
!          collocation point.
!
! IMPLEMENTED:
!   - Pure sinusoidal plunge  (h positive downward)
!
! TODO:
!   - Sinusoidal pitch
!   - Combined plunge + pitch (with arbitrary phase)
!   - Non-sinusoidal motions (ramp, arbitrary time series)
! =============================================================
module kinematics_mod
  use parameters_mod
  implicit none

contains

  ! ------------------------------------------------------------
  ! PLUNGE STATE  (fully implemented)
  !
  ! h(t)    = h0 * sin(omega*t)       [positive = downward]
  ! hdot(t) = h0 * omega * cos(omega*t)
  ! ------------------------------------------------------------
  subroutine plunge_state(t, h, hdot)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: h, hdot
    h    = h0    * sin(omega * t)
    hdot = h0    * omega * cos(omega * t)
  end subroutine plunge_state

  ! ------------------------------------------------------------
  ! PITCH STATE  (TODO – currently returns zero)
  !
  ! alpha(t)    = alpha0 * sin(omega*t + phi_pitch)
  ! alphadot(t) = alpha0 * omega * cos(omega*t + phi_pitch)
  ! ------------------------------------------------------------
  subroutine pitch_state(t, alpha, alphadot)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: alpha, alphadot

    ! TODO: Implement pitch kinematics.
    ! Uncomment and modify the lines below:
    !   alpha    = alpha0 * sin(omega * t + phi_pitch)
    !   alphadot = alpha0 * omega * cos(omega * t + phi_pitch)
    !
    ! For now, pitch is disabled:
    alpha    = 0.0_wp
    alphadot = 0.0_wp

    ! Suppress unused-variable warning
    if (.false.) then
      alpha = real(t, wp)
    end if
  end subroutine pitch_state

  ! ------------------------------------------------------------
  ! COMBINED STATE  (TODO)
  !
  ! When motion_type == 'COMBINED', call both plunge_state and
  ! pitch_state and let compute_normalwash combine them.
  ! No separate subroutine is strictly needed, but you could add
  ! one here for clarity or to handle special coupling laws.
  !
  ! TODO: Add any non-trivial coupling between pitch and plunge.
  ! ------------------------------------------------------------

  ! ------------------------------------------------------------
  ! COMPUTE KINEMATIC NORMAL WASH  (fully implemented)
  !
  ! The kinematic normal wash at collocation point j is the
  ! velocity of the plate surface (in the plate-normal direction)
  ! that the bound + wake vortices must cancel to satisfy the
  ! no-penetration condition.
  !
  !   w_j = hdot + U_inf*alpha + alphadot*(xc_j - x_pivot*chord)
  !
  ! Sign convention:
  !   Positive w_j means the plate pushes INTO the fluid from below
  !   (requires upwash from the vortices, i.e. positive Gamma).
  !   Boundary condition: ΣA(j,i)Γ_i + wake_terms = -w_j
  !
  ! Arguments:
  !   t      – current time
  !   xc_arr – collocation x-positions (size N)
  !   N      – number of panels
  !   w_arr  – output: normal wash at each collocation point
  ! ------------------------------------------------------------
  subroutine compute_normalwash(t, xc_arr, N, w_arr)
    real(wp), intent(in)  :: t
    integer,  intent(in)  :: N
    real(wp), intent(in)  :: xc_arr(N)
    real(wp), intent(out) :: w_arr(N)

    real(wp) :: h, hdot
    real(wp) :: alpha, alphadot
    integer  :: j

    call plunge_state(t, h, hdot)
    call pitch_state(t, alpha, alphadot)

    do j = 1, N
      w_arr(j) = hdot &
               + U_inf * alpha &
               + alphadot * (xc_arr(j) - x_pivot * chord)
    end do
  end subroutine compute_normalwash

end module kinematics_mod