! =============================================================
! MODULE: wake_mod
! Purpose: Storage and management of shed wake vortices.
!
! IMPLEMENTED:
!   - Frozen wake: vortices convect at U_inf in x only; y = 0.
!     This is the simplest physically correct model.  It ignores
!     the self-induced velocity of the wake and the actual vertical
!     position of the trailing edge during plunge.  Adequate for
!     small-amplitude motions and validation against Theodorsen.
!
! TODO (in order of increasing complexity):
!   1. Place new vortex at the ACTUAL trailing-edge y-position
!      (y_TE_inertial = -h(t) for frozen-y but correct x/y shed point)
!   2. ROLLUP wake: each vortex convects with freestream PLUS the
!      velocity induced on it by ALL other vortices (and bound).
!      Requires:
!        - A Biot-Savart sum over all other vortices per step
!        - Sub-cycling or adaptive dt for numerical stability
!        - Potentially vortex amalgamation for long-wake cases
!   3. Free-wake with vortex stretching (3-D extension)
! =============================================================
module wake_mod
  use parameters_mod
  implicit none

  integer  :: n_wake = 0           ! Current number of shed wake vortices

  ! Dynamic arrays – pre-allocated to maximum possible size
  real(wp), allocatable :: x_wake(:)      ! x-position of each wake vortex
  real(wp), allocatable :: y_wake(:)      ! y-position of each wake vortex
  real(wp), allocatable :: Gamma_wake(:)  ! Circulation of each wake vortex

contains

  ! ------------------------------------------------------------
  subroutine init_wake()
    integer :: max_wake
    max_wake = N_steps + 10   ! One vortex shed per time step
    allocate(x_wake(max_wake), y_wake(max_wake), Gamma_wake(max_wake))
    x_wake    = 0.0_wp
    y_wake    = 0.0_wp
    Gamma_wake = 0.0_wp
    n_wake    = 0
  end subroutine init_wake

  ! ------------------------------------------------------------
  ! Append one new wake vortex.
  ! ------------------------------------------------------------
  subroutine add_wake_vortex(x, y, Gamma)
    real(wp), intent(in) :: x, y, Gamma
    n_wake = n_wake + 1
    x_wake(n_wake)     = x
    y_wake(n_wake)     = y
    Gamma_wake(n_wake) = Gamma
  end subroutine add_wake_vortex

  ! ------------------------------------------------------------
  ! Convect all wake vortices by one time step.
  !
  ! FROZEN wake: x_w += U_inf * dt,  y_w unchanged (stays at 0).
  ! The "frozen" assumption means the wake lies along the x-axis
  ! regardless of the plate's vertical motion.  This overestimates
  ! the distance between successive wake vortices when the TE moves
  ! vertically but is accurate for small kh.
  !
  ! ROLLUP wake (TODO):
  !   For each vortex i, compute velocity induced by ALL other
  !   wake vortices AND all bound vortices, then advect:
  !     x_w(i) += (U_inf + u_induced(i)) * dt
  !     y_w(i) +=          v_induced(i)  * dt
  !   Use the regularised Biot-Savart from influence_mod.
  !   Call biot_savart() for every pair → O(n²) per step.
  !   Use sub-cycling (e.g. 4 sub-steps) for stability.
  ! ------------------------------------------------------------
  subroutine convect_wake()
    integer :: i

    select case (trim(wake_model))
    case ('FROZEN')
      do i = 1, n_wake
        x_wake(i) = x_wake(i) + U_inf * dt
        ! y_wake(i) is intentionally kept at 0
      end do

    case ('ROLLUP')
      ! TODO: Implement rollup wake convection.
      ! Steps:
      !   1. Allocate temporary u_ind(n_wake), v_ind(n_wake).
      !   2. For each vortex i, loop over all j ≠ i:
      !        call biot_savart(x_w(i),y_w(i), x_w(j),y_w(j), Gam_w(j), u,v)
      !        u_ind(i) += u;  v_ind(i) += v
      !      Also include influence of bound vortices (pass as argument or use
      !      module variables once geometry_mod is accessible here).
      !   3. x_w += (U_inf + u_ind)*dt
      !      y_w +=             v_ind *dt
      stop "ERROR: ROLLUP wake not yet implemented"

    case default
      write(*,*) "WARNING: Unknown wake_model '", trim(wake_model), "'. Using FROZEN."
      do i = 1, n_wake
        x_wake(i) = x_wake(i) + U_inf * dt
      end do
    end select
  end subroutine convect_wake

  ! ------------------------------------------------------------
  ! Return the sum of all wake vortex circulations.
  ! Used in the Kelvin constraint: Γ_w_new = -ΣΓ_bound - Γ_wake_sum
  ! ------------------------------------------------------------
  pure function total_wake_circulation() result(Gamma_sum)
    real(wp) :: Gamma_sum
    if (n_wake == 0) then
      Gamma_sum = 0.0_wp
    else
      Gamma_sum = sum(Gamma_wake(1:n_wake))
    end if
  end function total_wake_circulation

  ! ------------------------------------------------------------
  subroutine cleanup_wake()
    if (allocated(x_wake))     deallocate(x_wake)
    if (allocated(y_wake))     deallocate(y_wake)
    if (allocated(Gamma_wake)) deallocate(Gamma_wake)
    n_wake = 0
  end subroutine cleanup_wake

end module wake_mod