! =============================================================
! MODULE: influence_mod
! Purpose: Biot-Savart velocity calculations and aerodynamic
!          influence coefficient (AIC) matrix assembly.
!
! KEY FORMULA (derived from complex potential W = iΓ/(2π)*log(z-z₀)):
!
!   u_ind = +Γ * (yp - yv) / (2π r²_reg)
!   v_ind = -Γ * (xp - xv) / (2π r²_reg)
!
!   where r²_reg = (xp-xv)² + (yp-yv)² + ε²   (Lamb-Oseen regularisation)
!
! SIGN CHECK (critical to get right):
!   A CCW vortex (Γ>0) at the ORIGIN induces at point (1, 0):
!     v = -Γ*(1-0)/(2π*1) < 0   →  DOWNWARD flow to the right. ✓
!   A CCW vortex induces downward flow to its right and upward to its left.
!   This is consistent with the angular velocity vector pointing OUT OF the page.
!
! AIC COEFFICIENT A(j,i):
!   = v-component at collocation j from a UNIT vortex at bound position i
!   = -(xc_j - xb_i) / (2π r²_reg)
!   For flat plate (yc = yb = 0):  r² = (xc_j - xb_i)² + ε²
!   Since xc_j > xb_i always (3/4 > 1/4 within same panel, and LE panels
!   affect TE collocation points downstream):
!   → A(j,i) < 0 when j >= i  (downward v to the right of the vortex)
!   → A(j,i) > 0 when j < i   (upward v to the left of the vortex)
!
!   The geometry is fully described in geometry_mod.
! =============================================================
module influence_mod
  use parameters_mod
  implicit none

contains

  ! ------------------------------------------------------------
  ! Full 2D Biot-Savart: both velocity components induced by a
  ! single point vortex of strength Gamma at (xv, yv),
  ! evaluated at (xp, yp).
  ! ------------------------------------------------------------
  pure subroutine biot_savart(xp, yp, xv, yv, Gamma, u_ind, v_ind)
    real(wp), intent(in)  :: xp, yp, xv, yv, Gamma
    real(wp), intent(out) :: u_ind, v_ind

    real(wp) :: dx, dy, r2_reg

    dx     = xp - xv
    dy     = yp - yv
    r2_reg = dx*dx + dy*dy + eps_core*eps_core

    u_ind =  Gamma * dy / (TWO_PI * r2_reg)
    v_ind = -Gamma * dx / (TWO_PI * r2_reg)
  end subroutine biot_savart

  ! ------------------------------------------------------------
  ! Scalar AIC coefficient:
  !   A_coeff = v-component at (xp,yp) from unit vortex at (xv,yv)
  ! Equivalent to biot_savart with Gamma=1, returning v_ind only.
  ! Keeping this as a separate function avoids computing u when
  ! we only need v (the normal velocity for the flat plate BC).
  ! ------------------------------------------------------------
  pure function aic_coeff(xp, yp, xv, yv) result(A_coeff)
    real(wp), intent(in) :: xp, yp, xv, yv
    real(wp)             :: A_coeff

    real(wp) :: dx, dy, r2_reg

    dx     = xp - xv
    dy     = yp - yv
    r2_reg = dx*dx + dy*dy + eps_core*eps_core

    A_coeff = -dx / (TWO_PI * r2_reg)
  end function aic_coeff

  ! ------------------------------------------------------------
  ! Build the (N × N) aerodynamic influence coefficient matrix.
  !
  !   A_mat(j, i) = v at collocation j from unit Γ at bound vortex i
  !
  ! This matrix is CONSTANT for a flat plate (geometry does not
  ! change with airfoil motion) and is therefore computed ONCE
  ! before the time loop.
  !
  ! For cambered / NACA geometries (TODO):
  !   Replace v_ind with the normal-direction component:
  !     A(j,i) = u_ind * (-sin β_j) + v_ind * cos β_j
  !   where β_j is the local panel inclination angle.
  ! ------------------------------------------------------------
  subroutine build_aic_matrix(xb_arr, yb_arr, xc_arr, yc_arr, N, A_mat)
    integer,  intent(in)  :: N
    real(wp), intent(in)  :: xb_arr(N), yb_arr(N)
    real(wp), intent(in)  :: xc_arr(N), yc_arr(N)
    real(wp), intent(out) :: A_mat(N, N)

    integer :: i, j

    do i = 1, N          ! Column: source panel i (bound vortex)
      do j = 1, N        ! Row:    evaluation at collocation j
        A_mat(j, i) = aic_coeff(xc_arr(j), yc_arr(j), xb_arr(i), yb_arr(i))
      end do
    end do
  end subroutine build_aic_matrix

  ! ------------------------------------------------------------
  ! Compute v-component at point (xp,yp) induced by the ENTIRE
  ! wake (all n_w wake vortices).  Returns a scalar sum.
  ! Used in the RHS assembly of the linear system each step.
  ! ------------------------------------------------------------
  pure function wake_v_influence(xp, yp, n_w, x_w, y_w, Gam_w) result(v_sum)
    real(wp), intent(in) :: xp, yp
    integer,  intent(in) :: n_w
    real(wp), intent(in) :: x_w(n_w), y_w(n_w), Gam_w(n_w)

    real(wp) :: v_sum
    integer  :: k
    real(wp) :: u_tmp, v_tmp

    v_sum = 0.0_wp
    do k = 1, n_w
      call biot_savart(xp, yp, x_w(k), y_w(k), Gam_w(k), u_tmp, v_tmp)
      v_sum = v_sum + v_tmp
    end do
  end function wake_v_influence

end module influence_mod