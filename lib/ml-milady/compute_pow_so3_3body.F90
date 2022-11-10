
  subroutine compute_pow_so3_3body(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, &
                             numneigh_full, firstneigh_full)

    use T_kind_param_m, ONLY: double, kind_double, kind_double_complex
    use angular_functions, only: spherical_harm, grad_spherical_harm
    use ml_in_lammps_module, ONLY: rangml, debug, imm_neigh, one_pi,  n_rbf_so3, l_max,  &
                                pow_so3_dim, desc_forces, weighted, weighted_3ch, desc_forces

    !use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
    use module_so3, only: clmn, dclmn, W_pow_so3, coeff_rbf_so3, ini_rbf_so3, &
                          clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch
    use milady_data, only: local_psb3 => local_desc, local_psb3_deriv => local_desc_deriv, r_cut_desc, &
                           n_desc_type, k_desc_type, weight_per_type, weight_per_type_3ch

    implicit none
    integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
    integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
    integer, dimension(nmax), intent(in)  :: a_type
    real(double), dimension(ntype), intent(in) :: map_mass
    real(double), dimension(3,nmax), intent(in) :: xp

    integer   :: d_n_neigh
    !integer, dimension(imm, imm_neigh)    :: d_kind_neigh

    real(double), dimension(3) :: dxp_ji
    integer  :: ia, ja, iw, iw1, iw2, ia_n, i_count, ia_type, ja_type

    integer  :: p, l, m, n1, n2, i1, i2
    real(double)   :: r2_ji, r_ji, tmp_fact

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w_3ch
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w_3ch

    real(double), dimension(n_rbf_so3)     :: phi_ji, dphi_ji, rbf_ji, drbf_ji
    real(double)   :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch



    if (allocated(k_desc_type)) deallocate(k_desc_type) ; allocate(k_desc_type(numneigh_full))
    k_desc_type(:)=0
    local_psb3(:)=0.d0
    local_psb3_deriv(:,:,:)=0.d0

    do ja = i_start_at, i_final_at

      ja_type = a_type(ja)
      factor_ja=1.d0
      factor_ja_3ch=1.d0
      if (weighted) then
        factor_ja= weight_per_type(ja_type)
        if (weighted_3ch) factor_ja_3ch=  weight_per_type_3ch(ja_type)
      endif


      tmp_pow_so3_out(:) = 0.d0
      tmp_pow_so3_deriv_out(:, :, :) = 0.d0
      clmn(:, :) = (0.d0, 0.d0)
      if (weighted) then
        clmn_w(:, :) = (0.d0, 0.d0)
        tmp_pow_so3_out_w(:) = 0.d0
        tmp_pow_so3_deriv_out_w(:, :, :) = 0.d0
        if (weighted_3ch) then
          clmn_w_3ch(:, :) = (0.d0, 0.d0)
          tmp_pow_so3_out_w_3ch(:) = 0.d0
          tmp_pow_so3_deriv_out_w_3ch(:, :, :) = 0.d0
        end if
      end if

      ia_n = 0
      do iw = 1, numneigh_full
        ia = firstneigh_full(iw)
        ia_type = a_type(ia)

        dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
        r2_ji = Sum( dxp_ji(1:3)**2 )
        r_ji = dsqrt(r2_ji)

        if (r_ji >= r_cut_desc) cycle
        ia_n = ia_n + 1

        phi_ji(:) = 0.d0

        factor_ia = 1.d0
        factor_ia_3ch = 1.d0
        if (weighted) then
          factor_ia =  weight_per_type(ia_type)
          if (weighted_3ch) factor_ia_3ch = weight_per_type_3ch(ia_type)
        endif

        ! this gives clmn and d_clmn ...
        call spherical_3d(ia_n, factor_ia, factor_ia_3ch, desc_forces, dxp_ji, r_ji)
        k_desc_type(ia_n)=ia

      end do  ! ia end of neighbours iterations ...

      n_desc_type = ia_n
      ! update the componenet from the central atom ....
      i_count = 0
      do p = ini_rbf_so3, n_rbf_so3
        do l = 0, l_max
          i_count = i_count + 1
          tmp_fact =  dsqrt( (2.d0*dble(l)+1) / (4.d0 * one_pi) )
          clmn(0, i_count) = clmn(0, i_count) + tmp_fact
          if (weighted) then
            clmn_w(0, i_count) = clmn_w(0, i_count) + tmp_fact*factor_ja
            if (weighted_3ch) clmn_w_3ch(0, i_count) = clmn_w_3ch(0, i_count) + tmp_fact*factor_ja_3ch
          end if

        end do
      end do

      d_n_neigh = ia_n
      i_count = 0
      do n1 = ini_rbf_so3, n_rbf_so3
        do n2 = ini_rbf_so3, n_rbf_so3
          do l = 0, l_max
            i_count = i_count + 1
            i1 = (n1 - ini_rbf_so3)*(l_max + 1) + l + 1
            i2 = (n2 - ini_rbf_so3)*(l_max + 1) + l + 1
            do m = -l, l
              tmp_pow_so3_out(i_count) = tmp_pow_so3_out(i_count) + real(dconjg(clmn(m, i1))*clmn(m, i2), kind(0.d0))
              do ia = 1, d_n_neigh
                if (desc_forces) tmp_pow_so3_deriv_out(i_count, ia, 1:3) = &
                                   tmp_pow_so3_deriv_out(i_count, ia, 1:3) + &
                                   real( dconjg(clmn(m, i1))*dclmn(m, i2, ia, 1:3) + &
                                         clmn(m, i2)*dconjg( dclmn(m, i1, ia, 1:3)), kind(0.d0) )
              end do
              if (weighted) then
                tmp_pow_so3_out_w(i_count) = tmp_pow_so3_out_w(i_count) + &
                                             real(dconjg(clmn_w(m, i1))*clmn_w(m, i2), kind(0.d0))
                do ia = 1, d_n_neigh
                  if (desc_forces) tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) = &
                                      tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) + &
                                      real( dconjg(clmn_w(m, i1))*dclmn_w(m, i2, ia, 1:3) + &
                                            clmn_w(m, i2)*dconjg(dclmn_w(m, i1, ia, 1:3)), kind(0.d0) )
                end do
                if (weighted_3ch) then
                  tmp_pow_so3_out_w_3ch(i_count) = tmp_pow_so3_out_w_3ch(i_count) + &
                                                   real(dconjg(clmn_w_3ch(m, i1))*clmn_w_3ch(m, i2), kind(0.d0))
                  do ia = 1, d_n_neigh
                    if (desc_forces) tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) = &
                                        tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) + &
                                        real( dconjg(clmn_w_3ch(m, i1))*dclmn_w_3ch(m, i2, ia, 1:3) + &
                                              clmn_w_3ch(m, i2)*dconjg(dclmn_w_3ch(m, i1, ia, 1:3)), kind(0.d0) )
                  end do
                end if
              end if

            end do
          end do                  ! m
        end do                  ! l
      end do                  ! p


      ! forces .......
      if (desc_forces) then
        do ia = 1, d_n_neigh
          local_psb3_deriv(1:pow_so3_dim, ia, 1:3) = tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
          local_psb3_deriv(1:pow_so3_dim, 0, 1:3) = local_psb3_deriv(1:pow_so3_dim, 0, 1:3) - &
                                                    tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
        end do
        if (weighted) then
          do ia = 1, d_n_neigh
            local_psb3_deriv(pow_so3_dim+1:2*pow_so3_dim,  ia, 1:3) = &
                                  tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
            local_psb3_deriv(pow_so3_dim+1:2*pow_so3_dim,  0, 1:3) = &
                                  local_psb3_deriv(pow_so3_dim+1:2*pow_so3_dim, 0, 1:3) - &
                                  tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
          end do
          if (weighted_3ch) then
            do ia = 1, d_n_neigh
              local_psb3_deriv(2*pow_so3_dim+1:3*pow_so3_dim,  ia, 1:3) = &
                                tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
              local_psb3_deriv(2*pow_so3_dim+1:3*pow_so3_dim,  0, 1:3) = &
                                  local_psb3_deriv(2*pow_so3_dim+1:3*pow_so3_dim, 0, 1:3) - &
                                  tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
            end do
          end if
        end if
      end if

      ! energy .......
      local_psb3(1:pow_so3_dim) = tmp_pow_so3_out(1:pow_so3_dim)!*factor_ja
      if (weighted) then
        local_psb3(pow_so3_dim+1:2*pow_so3_dim) = tmp_pow_so3_out_w(1:pow_so3_dim)!*factor_ja
        if (weighted_3ch) then
          local_psb3(2*pow_so3_dim+1:3*pow_so3_dim) = tmp_pow_so3_out_w_3ch(1:pow_so3_dim)!*factor_ja_3ch
        end if
      end if

    end do   !ja main loop

  end subroutine compute_pow_so3_3body

  subroutine init_pow_so3_3body
  use ml_in_lammps_module, only: l_max, n_rbf_so3, pow_so3_dim, imm_neigh, &
                              weighted, weighted_3ch
  use module_so3, only: radial_pow_so3, radial_bartok, radial_sgg, &
                        clmn, dclmn, ini_rbf_so3, clmn_w, dclmn_w, &
                        clmn_w_3ch, dclmn_w_3ch

  implicit none

  if (radial_pow_so3 == radial_bartok) then
    ini_rbf_so3 = 1
    pow_so3_dim = int((1 + l_max))*n_rbf_so3**2
  end if


  if (radial_pow_so3 == radial_sgg) then
    ini_rbf_so3 = 0
    pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)**2
  end if

  if (allocated(clmn)) deallocate (clmn);    allocate (clmn(-l_max:l_max, pow_so3_dim))
  if (allocated(dclmn)) deallocate (dclmn); allocate (dclmn(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
  if (weighted) then
    if (allocated(clmn_w)) deallocate (clmn_w);    allocate (clmn_w(-l_max:l_max, pow_so3_dim))
    if (allocated(dclmn_w)) deallocate (dclmn_w); allocate (dclmn_w(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
    if (weighted_3ch) then
     if (allocated(clmn_w_3ch)) deallocate (clmn_w_3ch);    allocate (clmn_w_3ch(-l_max:l_max, pow_so3_dim))
     if (allocated(dclmn_w_3ch)) deallocate (dclmn_w_3ch); allocate (dclmn_w_3ch(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
    end if
  end if

end subroutine init_pow_so3_3body
