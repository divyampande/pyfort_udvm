! =============================================================
! MODULE: forces_mod
! Purpose: Compute aerodynamic force coefficients from the
!          bound vortex distribution.
!
! LIFT (Cl):
!   Uses the unsteady Bernoulli formulation, which combines:
!
!   (a) Quasi-steady (Kutta-Joukowski) term:
!         L_qs = ρ · U_inf · Σᵢ Γᵢ
!         Physical meaning: lift from steady circulation.
!
!   (b) Unsteady (added-mass) term from ∂Φ/∂t:
!         The perturbation potential at x is:
!           Φ(x,t) = Σ_{xb_i < x} Γᵢ(t)  [sum of bound Γ upstream of x]
!         Pressure difference across panel j:
!           ΔP_j = ρ [U_inf · γ_j + ∂Φ_j/∂t]
!         where Φ_j = Σ_{i=1}^{j} Γᵢ  (cumulative from LE)
!         Integrating over the chord:
!           L_uns = ρ · Σ_j (∂Φ_j/∂t) · Δx
!                 = ρ · (Δx/Δt) · Σ_j [Σᵢ₌₁ʲ ΔΓᵢ]
!         Rearranging the double sum (each ΔΓᵢ appears in panels j≥i):
!           L_uns = ρ · (Δx/Δt) · Σᵢ (N−i+1) · ΔΓᵢ
!         where ΔΓᵢ = Γᵢ_new − Γᵢ_old.
!
!   Total:  Cl = (L_qs + L_uns) / (½ ρ U_inf² c)
!
! THRUST (Ct) – TODO:
!   For a flat plate, thrust arises from LEADING EDGE SUCTION.
!   The discrete suction force is:
!     T = π ρ C_s²
!   where C_s is the LE singularity strength, approximated as:
!     C_s ≈ Γ₁ / sqrt(2 · xb₁)   [from thin airfoil theory]
!   Ct = T / (½ ρ U_inf² c)
!   NOTE: Leading edge suction is a delicate quantity; validate
!   the lift implementation first before adding thrust.
!
! MOMENT (Cm_qc) about quarter-chord:
!   Quasi-steady only for now.
!   TODO: add unsteady ∂Φ/∂t contribution to moment.
! =============================================================
module forces_mod
  use parameters_mod
  implicit none

contains

  ! ============================================================
  ! Compute Cl, Cm (about quarter-chord), and Ct.
  !
  ! Arguments:
  !   Gamma_new – bound circulations at current step (size N)
  !   Gamma_old – bound circulations at previous step (size N)
  !   xb_arr    – bound vortex x-positions (size N)
  !   N         – number of panels
  ! ============================================================
  subroutine compute_forces(Gamma_new, Gamma_old, xb_arr, N, Cl, Cm_qc, Ct)
    integer,  intent(in)  :: N
    real(wp), intent(in)  :: Gamma_new(N), Gamma_old(N)
    real(wp), intent(in)  :: xb_arr(N)
    real(wp), intent(out) :: Cl, Cm_qc, Ct

    real(wp) :: dx
    real(wp) :: Gamma_total
    real(wp) :: L_qs, L_uns, L_total
    real(wp) :: dGamma(N), unsteady_sum
    real(wp) :: Cm_sum
    real(wp) :: dyn_press_c    ! 0.5 * rho * U_inf^2 * chord
    integer  :: i

    dx         = chord / real(N, wp)
    dyn_press_c = 0.5_wp * rho * U_inf**2 * chord

    ! ---- (a) Quasi-steady lift --------------------------------
    Gamma_total = sum(Gamma_new)
    L_qs        = rho * U_inf * Gamma_total

    ! ---- (b) Unsteady (added-mass) lift ----------------------
    ! ΔΓᵢ for each panel
    dGamma = Gamma_new - Gamma_old

    ! Weighted sum:  Σᵢ (N−i+1) · ΔΓᵢ
    unsteady_sum = 0.0_wp
    do i = 1, N
      unsteady_sum = unsteady_sum + real(N - i + 1, wp) * dGamma(i)
    end do

    ! L_uns = ρ · (Δx/Δt) · Σᵢ (N−i+1)·ΔΓᵢ
    L_uns = rho * (dx / dt) * unsteady_sum

    ! ---- Total lift coefficient ------------------------------
    L_total = L_qs + L_uns
    Cl      = L_total / dyn_press_c

    ! ---- Pitching moment about quarter-chord (quasi-steady) --
    ! M_qc = ρ U_inf Σᵢ Γᵢ · (xb_i − c/4)   [nose-up positive]
    ! Cm_qc = M_qc / (½ ρ U²  c²)
    Cm_sum = 0.0_wp
    do i = 1, N
      Cm_sum = Cm_sum + Gamma_new(i) * (xb_arr(i) - 0.25_wp * chord)
    end do
    Cm_qc = rho * U_inf * Cm_sum / (dyn_press_c * chord)

    ! TODO: Add unsteady moment term (∂Φ/∂t weighted by moment arm).
    ! The structure is the same as L_uns but with (xb_i − c/4) as weight.

    ! ---- Thrust coefficient (Leading Edge Suction) -----------
    ! TODO: Implement Ct from leading edge suction.
    ! Temporary quasi-steady approximation:
    !   For a plunging flat plate with effective AoA ≈ ḣ/U_inf,
    !   the thrust in the quasi-steady limit is:
    !     Ct ≈ Cl · (ḣ_now / U_inf)
    !   but this requires passing ḣ(t) here.  For now return 0.
    !
    ! Full formula (once you are ready):
    !   C_s  = Gamma_new(1) / sqrt(2.0_wp * xb_arr(1))
    !   T    = PI * rho * C_s**2
    !   Ct   = T / dyn_press_c
    Ct = 0.0_wp

  end subroutine compute_forces

end module forces_mod