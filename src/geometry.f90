! =============================================================
! MODULE: geometry_mod
! Purpose: Defines the airfoil discretisation.
!
! IMPLEMENTED:
!   - Flat plate using the classical 1/4–3/4 panel rule
!       * Bound vortex at 1/4-chord of each panel (xb)
!       * Collocation point at 3/4-chord of each panel (xc)
!       * Both on y = 0  (flat plate)
!     This is the standard discrete thin-airfoil stagger that
!     recovers Cl = 2π α in the single-panel limit and satisfies
!     the Kutta condition automatically through the panel equations.
!
! TODO:
!   - Cambered plate (give xb, xc the correct y-offsets)
!   - NACA 4-digit profile  (use init_naca4)
!   - NACA 5-digit profile  (use init_naca5)
!   - Source + vortex panel method for thick airfoils
!     (requires solving a larger system; vortex-only is fine for
!      thin-airfoil / Cl / Ct purposes)
! =============================================================
module geometry_mod
  use parameters_mod
  implicit none

  ! Panel geometry arrays (allocated by init_geometry)
  real(wp), allocatable :: xb(:)          ! Bound vortex x-positions
  real(wp), allocatable :: yb(:)          ! Bound vortex y-positions
  real(wp), allocatable :: xc(:)          ! Collocation x-positions
  real(wp), allocatable :: yc(:)          ! Collocation y-positions
  real(wp), allocatable :: panel_dx(:)    ! Panel chord-wise length Δx

  ! Trailing edge position in the body frame
  ! (In the inertial frame, y_TE changes during plunge; the wake
  !  module uses this for rollup wake placement.)
  real(wp) :: x_TE_body   ! Always = chord (TE at the end of the plate)
  real(wp) :: y_TE_body   ! Always = 0     (flat plate)

contains

  ! ------------------------------------------------------------
  ! Flat-plate discretisation using the 1/4-3/4 rule.
  !
  ! Panels are numbered i = 1 … N_panels, ordered LE→TE.
  ! Panel i spans from (i-1)*Δx to i*Δx.
  !   xb(i) = (i - 3/4) * Δx   ← quarter-chord of panel i
  !   xc(i) = (i - 1/4) * Δx   ← three-quarter-chord of panel i
  !
  ! WHY THIS WORKS:
  !   The influence coefficient A(j,i) = v at xc(j) from unit Γ at xb(i)
  !   = -1/(2π*(xc(j)-xb(i))) for the flat plate.
  !   For j=i: xc(i)-xb(i) = Δx/2 > 0  →  A(i,i) < 0.
  !   All xc(j) > xb(i) when j >= i  →  diagonal and lower-left block < 0.
  !   Upper-right block (j < i): xc(j) < xb(i)  →  coefficient > 0.
  !   After Kelvin substitution the effective matrix is well-conditioned.
  ! ------------------------------------------------------------
  subroutine init_geometry()
    integer  :: i
    real(wp) :: dx

    allocate(xb(N_panels), yb(N_panels))
    allocate(xc(N_panels), yc(N_panels))
    allocate(panel_dx(N_panels))

    dx = chord / real(N_panels, wp)

    do i = 1, N_panels
      xb(i)       = (real(i, wp) - 0.75_wp) * dx
      yb(i)       = 0.0_wp
      xc(i)       = (real(i, wp) - 0.25_wp) * dx
      yc(i)       = 0.0_wp
      panel_dx(i) = dx
    end do

    x_TE_body = chord
    y_TE_body = 0.0_wp
  end subroutine init_geometry

  ! ------------------------------------------------------------
  subroutine cleanup_geometry()
    if (allocated(xb))       deallocate(xb)
    if (allocated(yb))       deallocate(yb)
    if (allocated(xc))       deallocate(xc)
    if (allocated(yc))       deallocate(yc)
    if (allocated(panel_dx)) deallocate(panel_dx)
  end subroutine cleanup_geometry

  ! ------------------------------------------------------------
  ! TODO: NACA 4-digit geometry
  ! Description:
  !   Build xb, xc, yb, yc for a cambered/thick NACA airfoil.
  !   Use the standard NACA 4-digit formulae to get the camber
  !   line and thickness distribution.  Place bound vortices at
  !   the quarter-chord of each panel along the camber line and
  !   collocation points at 3/4-chord.  Normal vectors must point
  !   perpendicular to the local panel tangent.
  !   The boundary condition becomes:
  !     w_j = (U_inf * tan(β_j) + hdot)*cos(β_j) + ...
  !   where β_j is the local camber angle.
  !
  ! subroutine init_naca4(digits)
  !   character(len=4), intent(in) :: digits  ! e.g. '2412'
  !   ...
  ! end subroutine
  !
  ! TODO: NACA 5-digit geometry (similar, different camber formula)
  ! ------------------------------------------------------------

end module geometry_mod