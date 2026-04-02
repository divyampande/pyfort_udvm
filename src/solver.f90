! =============================================================
! MODULE: solver_mod
! Purpose: Linear system solver and the core DVM time-step routine.
!
! THEORY: At each time step m the unknowns are the N bound vortex
! strengths Γ_i (i=1…N).  The boundary condition at collocation
! point j is:
!
!   ΣᵢA(j,i)·Γᵢ  +  B_new(j)·Γ_w_new  +  Σ_k B_old(j,k)·Γ_wk  =  -w_j
!
! where A is the pre-built AIC matrix, B_new / B_old are the
! influences of the NEW and OLD wake vortices, and w_j is the
! kinematic wash from kinematics_mod.
!
! KELVIN SUBSTITUTION:
! Kelvin's circulation theorem (total circulation = 0 always):
!   Σᵢ Γ_i  +  Γ_w_new  +  Σ_k Γ_wk  =  0
!   → Γ_w_new = -Σᵢ Γ_i - Σ_k Γ_wk   ← substitute into BC
!
! After substitution the system becomes:
!   Σᵢ [A(j,i) - B_new(j)] · Γ_i  =  -w_j + B_new(j)·Γ_wake_sum
!                                           - Σ_k B_old(j,k)·Γ_wk
!
! Effective matrix:  M(j,i) = A(j,i) - B_new(j)
!   (B_new(j) is CONSTANT across all columns i — rank-1 modification)
! RHS vector:        rhs(j) = -w_j + B_new(j)·Γ_wake_sum - Σ_k B_old(j,k)·Γ_wk
!
! After solving for Γ_i, recover:
!   Γ_w_new = -Σᵢ Γ_i - Γ_wake_sum
!
! NEW WAKE VORTEX PLACEMENT:
!   x_w_new = x_TE + 0.3·U_inf·dt   (standard in literature)
!   y_w_new = 0                       (frozen wake)
!   After the solve, convect ALL wake vortices (including new one)
!   by U_inf·dt before the next step.
! =============================================================
module solver_mod
  use parameters_mod
  use influence_mod
  use wake_mod
  implicit none

contains

  ! ============================================================
  ! Gaussian elimination with partial pivoting.
  ! Solves A_mat · x = b_vec.
  ! Both A_mat and b_vec are modified in place; on exit b_vec = x.
  ! ============================================================
  subroutine gauss_solve(A_mat, b_vec, N)
    integer,  intent(in)    :: N
    real(wp), intent(inout) :: A_mat(N, N)
    real(wp), intent(inout) :: b_vec(N)

    integer  :: i, j, k, pivot_row
    real(wp) :: pivot_val, factor, tmp
    real(wp) :: tmp_row(N)

    ! --- Forward elimination ---
    do k = 1, N

      ! Find partial pivot
      pivot_row = k
      pivot_val = abs(A_mat(k, k))
      do i = k + 1, N
        if (abs(A_mat(i, k)) > pivot_val) then
          pivot_val = abs(A_mat(i, k))
          pivot_row = i
        end if
      end do

      ! Swap rows k and pivot_row
      if (pivot_row /= k) then
        tmp_row          = A_mat(k, :)
        A_mat(k, :)      = A_mat(pivot_row, :)
        A_mat(pivot_row, :) = tmp_row
        tmp              = b_vec(k)
        b_vec(k)         = b_vec(pivot_row)
        b_vec(pivot_row) = tmp
      end if

      if (abs(A_mat(k, k)) < 1.0e-14_wp) then
        write(*, '(A,I4,A)') "  WARNING: Near-singular matrix at pivot row ", k, &
                              ". Check N_panels or eps_core."
      end if

      ! Eliminate rows below k
      do i = k + 1, N
        factor = A_mat(i, k) / A_mat(k, k)
        do j = k, N
          A_mat(i, j) = A_mat(i, j) - factor * A_mat(k, j)
        end do
        b_vec(i) = b_vec(i) - factor * b_vec(k)
      end do

    end do

    ! --- Back substitution ---
    do i = N, 1, -1
      do j = i + 1, N
        b_vec(i) = b_vec(i) - A_mat(i, j) * b_vec(j)
      end do
      b_vec(i) = b_vec(i) / A_mat(i, i)
    end do

  end subroutine gauss_solve

  ! ============================================================
  ! ONE DVM TIME STEP
  !
  ! Input:
  !   A_aic        – (N×N) constant AIC matrix (from influence_mod)
  !   xb,yb,xc,yc  – panel geometry (from geometry_mod)
  !   N            – number of panels
  !   w_arr        – kinematic normal wash at each collocation pt
  !   x_TE, y_TE   – trailing edge position (body frame; y_TE=0 here)
  !   Gamma_bound  – bound circulations from PREVIOUS step
  !
  ! Output:
  !   Gamma_bound_new  – updated bound circulations
  !
  ! Side effects (via wake_mod):
  !   - One new wake vortex is appended to wake arrays
  !   - All wake vortices (including new) are convected by U_inf*dt
  ! ============================================================
  subroutine solver_step(A_aic, xb, yb, xc, yc, N, w_arr, &
                         x_TE, y_TE, Gamma_bound, Gamma_bound_new)
    integer,  intent(in)  :: N
    real(wp), intent(in)  :: A_aic(N, N)
    real(wp), intent(in)  :: xb(N), yb(N)
    real(wp), intent(in)  :: xc(N), yc(N)
    real(wp), intent(in)  :: w_arr(N)
    real(wp), intent(in)  :: x_TE, y_TE
    real(wp), intent(in)  :: Gamma_bound(N)   ! Previous step
    real(wp), intent(out) :: Gamma_bound_new(N)

    ! Local workspace
    real(wp) :: M_work(N, N)   ! Effective system matrix (will be modified by gauss_solve)
    real(wp) :: rhs(N)         ! Right-hand side

    real(wp) :: B_new(N)       ! Influence of new wake vortex at collocation pts
    real(wp) :: x_w_new, y_w_new
    real(wp) :: Gamma_wake_sum
    real(wp) :: Gamma_w_new
    integer  :: i, j, k

    ! ---- 1. Position of the new wake vortex ----------------
    !   Placed 0.3·U_inf·dt downstream of the TE.
    !   This fractional offset (0 < f < 1) is a model parameter;
    !   0.3 is standard (Katz & Plotkin §13).
    !   For frozen wake, y is set to 0 (ignoring TE vertical motion).
    !   For rollup wake (TODO), y_w_new = y_TE_inertial = -h(t).
    x_w_new = x_TE + 0.3_wp * U_inf * dt
    y_w_new = 0.0_wp    ! Frozen wake: y = 0

    ! ---- 2. Influence of new wake vortex at each collocation pt ----
    do j = 1, N
      B_new(j) = aic_coeff(xc(j), yc(j), x_w_new, y_w_new)
    end do

    ! ---- 3. Sum of OLD wake circulations -------------------
    Gamma_wake_sum = total_wake_circulation()

    ! ---- 4. Effective matrix after Kelvin substitution -----
    !   M(j,i) = A(j,i) - B_new(j)
    !   B_new(j) is independent of i → same value subtracted from
    !   every column in row j.
    do i = 1, N
      do j = 1, N
        M_work(j, i) = A_aic(j, i) - B_new(j)
      end do
    end do

    ! ---- 5. Right-hand side --------------------------------
    !   rhs(j) = -w(j)
    !           + B_new(j) * Gamma_wake_sum          [Kelvin term]
    !           - Σ_k B_old(j,k) * Gamma_wake(k)    [old wake influence]
    do j = 1, N
      rhs(j) = -w_arr(j) + B_new(j) * Gamma_wake_sum

      ! Subtract the influence of each already-shed wake vortex
      do k = 1, n_wake
        rhs(j) = rhs(j) &
               - aic_coeff(xc(j), yc(j), x_wake(k), y_wake(k)) &
                 * Gamma_wake(k)
      end do
    end do

    ! ---- 6. Solve M_work · Gamma_bound_new = rhs ----------
    Gamma_bound_new = rhs
    call gauss_solve(M_work, Gamma_bound_new, N)

    ! ---- 7. New wake vortex strength (Kelvin) --------------
    !   Total circulation = ΣΓ_bound + Γ_w_new + Γ_wake_sum = 0
    Gamma_w_new = -sum(Gamma_bound_new) - Gamma_wake_sum

    ! ---- 8. Append new vortex to wake arrays ---------------
    call add_wake_vortex(x_w_new, y_w_new, Gamma_w_new)

    ! ---- 9. Convect ALL wake vortices (frozen or rollup) ---
    call convect_wake()

    ! Suppress unused-variable warning for yb (not needed for flat plate,
    ! but present for future cambered geometry)
    if (.false.) then
      i = int(yb(1))
    end if

  end subroutine solver_step

end module solver_mod