module module_neigh_local
  USE module_kind_variables, ONLY: kind_double

  ! ja atom
  integer  :: ja_atom, iw2
  ! the number of neighbours of ja_number
  integer  :: max_neigh_local

  integer, dimension(:), allocatable     :: i_type, i_central
  ! r_central   - all the distance with respect central ja atom
  ! ur_central  - the same thing as previuos but with transformated distances
  ! r_fcut      - the fcut function

  real(kind_double), dimension(:), allocatable :: r_central, ur_central, r_fcut
  ! tmp_xp(3,:) = xp(3,ia_n)
  ! tmp_dxp(3,:) = xp(3,ia_n) - xp(3,ja)
  ! d_r_central(3,:)  - the derivatives
  ! d_ur_central(3,:) - the derivatives
  ! d_r_fcut(3,:) - the derivatives
  real(kind_double), dimension(:, :), allocatable    :: &
    d_r_central, d_ur_central, d_r_fcut, tmpcos_dxp, tmp_dxp, tmp_xp
end module module_neigh_local



subroutine local_neighbours(ja, nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
  use T_kind_param_m, only: kind_double, double
  use module_neigh_local, only: r_central, i_central, i_type, tmp_dxp, tmp_xp, &
                                max_neigh_local, iw2
  use ml_in_lammps_module, only: imm_neigh
  use milady_data, only: r_cut_desc, k_desc_type, n_desc_type

  implicit none

  integer, intent(in)  :: ja
  integer, intent (in) ::  nmax, ntype, numneigh_full
  integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
  integer, dimension(nmax), intent(in)  :: a_type
  real(double), dimension(ntype), intent(in) :: map_mass
  real(double), dimension(3,nmax), intent(in) :: xp

  integer  :: iw, iw1, ia_n, ia, max_neigh, ia_type
  real(kind_double)    :: r_ji, r2_ji
  real(kind_double), dimension(3) :: dxp_ji


  max_neigh = imm_neigh
  ! begin small box or not 1/
  if (allocated(r_central)) deallocate (r_central); allocate (r_central(max_neigh))
  if (allocated(i_central)) deallocate (i_central); allocate (i_central(0:max_neigh))
  if (allocated(i_type)) deallocate (i_type); allocate (i_type(0:max_neigh))
  if (allocated(tmp_dxp)) deallocate (tmp_dxp); allocate (tmp_dxp(3, 0:max_neigh))
  if (allocated(tmp_xp)) deallocate (tmp_xp); allocate (tmp_xp(3, 0:max_neigh))

  ia_n = 0
  i_central(0) = ja
  i_type(0) = a_type(ja)
  do iw = 1, numneigh_full
    ia = firstneigh_full(iw)
    ia_type = a_type(ia)

    dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
    r2_ji = Sum( dxp_ji(1:3)**2 )
    r_ji = dsqrt(r2_ji)

    if (r_ji >= r_cut_desc) cycle
    ia_n = ia_n + 1
    r_central(ia_n) = r_ji
    i_central(ia_n) = ia
    i_type(ia_n) = ia_type
    tmp_dxp(:, ia_n) = dxp_ji(:)
    tmp_xp(:, ia_n) = xp(:, ia)
    k_desc_type(ia_n) = ia
  end do
  max_neigh_local = ia_n
  n_desc_type = ia_n

end subroutine local_neighbours
