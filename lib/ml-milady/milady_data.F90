module milady_data
      use T_kind_param_m, only: double
      implicit none
      real(double), dimension(:,:), allocatable :: weights
      real(double), dimension(:), allocatable :: details, rcutvec
      !real(double)   :: weight_factor
      real(double), dimension(:), allocatable :: weight_factor
      real(double), dimension(:), allocatable :: weight_per_type, weight_per_type_3ch, ref_energy_per_atom

      integer :: nelements, nweights, ninteractions, ndetails, ext_descriptor_type

      real(double), dimension(:),  allocatable   :: local_desc
      real(double), dimension(:, :, :), allocatable  :: local_desc_deriv
      real(double), dimension(:),  allocatable   :: q_local_desc
      real(double), dimension(:,:),  allocatable   :: q_local_desc_deriv

      real(double), dimension(:),  allocatable   :: local_desc_kernel
      real(double), dimension(:, :, :), allocatable  :: local_desc_kernel_deriv


      real(double), dimension(:,:), allocatable :: st_temp

      real(double) :: r_cut_desc
      integer, dimension(:), allocatable :: k_desc_type
      integer :: n_desc_type

    contains

   subroutine allocate_local_desc
      use temporary_data_cov, only : dim_xdesc, dim_xdesc_quadratic, dim_xdesc_kernel, &
                                     dim_xdesc_quadratic_only, dim_xdesc_kernel_only
      use ml_in_lammps_module, only: imm_neigh, snap_order, snap_quadratic, snap_linear, snap_kernel
      implicit none

      if (allocated(local_desc)) deallocate(local_desc) ; allocate(local_desc(dim_xdesc))
      if (allocated(local_desc_deriv)) deallocate(local_desc_deriv)
                                       allocate(local_desc_deriv(dim_xdesc, 0:imm_neigh, 3))
      if (snap_order == snap_quadratic) then
            if (allocated(q_local_desc)) deallocate(q_local_desc) ; allocate(q_local_desc(dim_xdesc_quadratic_only))
            if(allocated(q_local_desc_deriv)) deallocate(q_local_desc_deriv)
                                              allocate(q_local_desc_deriv(dim_xdesc_quadratic_only, 3))
      end if

      if (snap_order == snap_kernel) then
            if (allocated(local_desc_kernel)) deallocate(local_desc_kernel)
                                              allocate(local_desc_kernel(dim_xdesc_kernel_only))
            if (allocated(local_desc_kernel_deriv)) deallocate(local_desc_kernel_deriv)
                                              allocate(local_desc_kernel_deriv(dim_xdesc_kernel_only, 0:imm_neigh, 3))
      end if


      if (snap_order == snap_linear) then
            if (allocated(st_temp)) deallocate(st_temp) ; allocate(st_temp(dim_xdesc, 6))
      end if

      if (snap_order == snap_quadratic) then
            if (allocated(st_temp)) deallocate(st_temp) ; allocate(st_temp(dim_xdesc_quadratic-1, 6))
      end if


      if (snap_order == snap_kernel) then
            if (allocated(st_temp)) deallocate(st_temp) ; allocate(st_temp(dim_xdesc_kernel-1, 6))
      end if


  end subroutine allocate_local_desc

end module milady_data


module  milady_interfaces
!
  interface


  subroutine  compute_g2(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_g2 => local_desc , local_g2_deriv => local_desc_deriv, &
                            r_cut_desc, n_desc_type, k_desc_type
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine compute_g2


  subroutine  compute_afs(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_afs => local_desc , local_afs_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine compute_afs

  subroutine  compute_pow_so4(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_pso4 => local_desc , local_pso4_deriv => local_desc_deriv, r_cut_desc, &
                            n_desc_type, k_desc_type, weight_per_type, weight_per_type_3ch
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine  compute_pow_so4


  subroutine  compute_bispectrum_so4(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_bso4 => local_desc , local_bso4_deriv => local_desc_deriv, r_cut_desc, &
                            n_desc_type, k_desc_type, weight_per_type, weight_per_type_3ch
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine  compute_bispectrum_so4


  subroutine  compute_pow_so3(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_pso3 => local_desc , local_pso3_deriv => local_desc_deriv, r_cut_desc, &
                            n_desc_type, k_desc_type, weight_per_type, weight_per_type_3ch
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine  compute_pow_so3

  subroutine  compute_pow_so3_3body(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass,  numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use milady_data, only:  local_psb3 => local_desc , local_psb3_deriv => local_desc_deriv, r_cut_desc, &
                            n_desc_type, k_desc_type, weight_per_type, weight_per_type_3ch
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine  compute_pow_so3_3body

  subroutine  compute_soap(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, numneigh_full, &
                        firstneigh_full)
    USE T_kind_param_m, ONLY:  double
    use angular_functions
    use ml_in_lammps_module, ONLY: rangml,l_max, j_max, imm_neigh,    &
                            weighted, soap_dim, desc_forces, factor_weight_mass, &
                            r_cut_width_soap, alpha_soap, one_pi, sqrt_two, rb_soap, W_soap, &
                            S_factor_matrix, ns_soap_index, n_soap,  &
                            lsoap_diag, lsoap_norm, lsoap_lnorm, ns_soap_index, nspecies_soap, c_zero, vec_3d_zero
    use milady_data, only: local_soap => local_desc, local_soap_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type
    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp
  end subroutine compute_soap

end interface
!
end module  milady_interfaces
