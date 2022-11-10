
  subroutine compute_pow_so3(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, &
                             numneigh_full, firstneigh_full)

    use T_kind_param_m, ONLY: double, kind_double, kind_double_complex
    use angular_functions, only: spherical_harm, grad_spherical_harm
    use ml_in_lammps_module, ONLY: rangml, debug, imm_neigh, one_pi,  n_rbf_so3, l_max,  &
                                pow_so3_dim, desc_forces, weighted, weighted_3ch, desc_forces

    !use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
    use module_so3, only: clmn, dclmn, W_pow_so3, coeff_rbf_so3, ini_rbf_so3, &
                          clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch
    use milady_data, only: local_pso3 => local_desc, local_pso3_deriv => local_desc_deriv, r_cut_desc, &
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

    integer  :: p, l, m
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
    local_pso3(:)=0.d0
    local_pso3_deriv(:,:,:)=0.d0

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
      do p = ini_rbf_so3, n_rbf_so3
        do l = 0, l_max
          i_count = i_count + 1
          do m = -l, l
            tmp_pow_so3_out(i_count) = tmp_pow_so3_out(i_count) + real(dconjg(clmn(m, i_count))*clmn(m, i_count), kind(0.d0))
            if (weighted) then
              tmp_pow_so3_out_w(i_count) = tmp_pow_so3_out_w(i_count) + &
                                           real(dconjg(clmn_w(m, i_count))*clmn_w(m, i_count), kind(0.d0))
              if (weighted_3ch) then
                tmp_pow_so3_out_w_3ch(i_count) = &
                                        tmp_pow_so3_out_w_3ch(i_count) + &
                                        real(dconjg(clmn_w_3ch(m, i_count))*clmn_w_3ch(m, i_count), kind(0.d0))
              end if
            end if
            do ia = 1, d_n_neigh
              if (desc_forces) tmp_pow_so3_deriv_out(i_count, ia, 1:3) = &
                                tmp_pow_so3_deriv_out(i_count, ia, 1:3) + &
                                2.d0*real(dconjg(clmn(m, i_count))*dclmn(m, i_count, ia, 1:3), kind(0.d0))
            end do
            if (weighted) then
              do ia = 1, d_n_neigh
                if (desc_forces) tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) = &
                                 tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) + &
                                 2.d0*real(dconjg(clmn_w(m, i_count))*dclmn_w(m, i_count, ia, 1:3), kind(0.d0))
              end do
              if (weighted_3ch) then
                do ia = 1, d_n_neigh
                  if (desc_forces) tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) = &
                                  tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) + &
                                  2.d0*real(dconjg(clmn_w_3ch(m, i_count))*dclmn_w_3ch(m, i_count, ia, 1:3), kind(0.d0))
                end do
              end if
            end if
          end do                  ! m
        end do                  ! p
      end do                  ! l
      ! forces .......
      if (desc_forces) then
        do ia = 1, d_n_neigh
          local_pso3_deriv(1:pow_so3_dim, ia, 1:3) = tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
          local_pso3_deriv(1:pow_so3_dim, 0, 1:3) = local_pso3_deriv(1:pow_so3_dim, 0, 1:3) - &
                                                    tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
        end do
        if (weighted) then
          do ia = 1, d_n_neigh
            local_pso3_deriv(pow_so3_dim+1:2*pow_so3_dim,  ia, 1:3) = &
                  tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
            local_pso3_deriv(pow_so3_dim+1:2*pow_so3_dim,  0, 1:3) = &
                   local_pso3_deriv(pow_so3_dim+1:2*pow_so3_dim, 0, 1:3) - &
                   tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
          end do
          if (weighted_3ch) then
            do ia = 1, d_n_neigh
              local_pso3_deriv(2*pow_so3_dim+1:3*pow_so3_dim,  ia, 1:3) = &
                    tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
              local_pso3_deriv(2*pow_so3_dim+1:3*pow_so3_dim,  0, 1:3) = &
                    local_pso3_deriv(2*pow_so3_dim+1:3*pow_so3_dim, 0, 1:3) - &
                    tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
            end do
          end if
        end if
      end if

      ! energy .......
      local_pso3(1:pow_so3_dim) = tmp_pow_so3_out(1:pow_so3_dim)!*factor_ja
      if (weighted) then
        local_pso3(pow_so3_dim+1:2*pow_so3_dim) = tmp_pow_so3_out_w(1:pow_so3_dim)!*factor_ja
        if (weighted_3ch) then
          local_pso3(2*pow_so3_dim+1:3*pow_so3_dim) = tmp_pow_so3_out_w_3ch(1:pow_so3_dim)!*factor_ja_3ch
        end if
      end if

    end do   !ja main loop

  end subroutine compute_pow_so3


  subroutine init_pow_so3()

    use ml_in_lammps_module, only: l_max, n_rbf_so3, pow_so3_dim, imm_neigh, &
                                weighted, weighted_3ch
    use module_so3, only: radial_pow_so3, radial_bartok, radial_sgg, clmn, dclmn, &
                          clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch,  ini_rbf_so3

    implicit none

    if (radial_pow_so3 == radial_bartok) then
      ini_rbf_so3 = 1
      pow_so3_dim = int((1 + l_max))*n_rbf_so3
    end if
    if (radial_pow_so3 == radial_sgg) then
      ini_rbf_so3 = 0
      pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)
    end if
    if (allocated(clmn)) deallocate (clmn); allocate (clmn(-l_max:l_max, pow_so3_dim))
    if (allocated(dclmn)) deallocate (dclmn); allocate (dclmn(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
    if (weighted) then
      if (allocated(clmn_w)) deallocate (clmn_w); allocate (clmn_w(-l_max:l_max, pow_so3_dim))
      if (allocated(dclmn_w)) deallocate (dclmn_w); allocate (dclmn_w(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
      if (weighted_3ch) then
        if (allocated(clmn_w_3ch)) deallocate (clmn_w_3ch); allocate (clmn_w_3ch(-l_max:l_max, pow_so3_dim))
        if (allocated(dclmn_w_3ch)) deallocate (dclmn_w_3ch); allocate (dclmn_w_3ch(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
      end if
    end if

  end subroutine init_pow_so3


  subroutine init_pow_so3_rbf()
    use ml_in_lammps_module, only: n_rbf_so3
    use module_so3, only: coeff_rbf_so3, W_pow_so3, &
                          chebT, chebU, d_chebT, &
                          scgg, d_scgg, radial_pow_so3, radial_bartok, &
                          radial_pow_so3, radial_sgg
    use milady_data, only: r_cut_desc
    !use compute_afs_mod
    implicit none


    real(kind(0.d0)), dimension(n_rbf_so3, n_rbf_so3)  :: S, V
    real(kind(0.d0)), dimension(n_rbf_so3) :: L
    !local
    integer  :: p, q, nb, ilaenv, lwork


    if (radial_pow_so3 == radial_bartok) then
      if (allocated(W_pow_so3)) deallocate (W_pow_so3); allocate (W_pow_so3(n_rbf_so3, n_rbf_so3))
      if (allocated(coeff_rbf_so3)) deallocate (coeff_rbf_so3); allocate (coeff_rbf_so3(n_rbf_so3))

      do q = 1, n_rbf_so3
        do p = 1, n_rbf_so3
          !Bartok paper
          S(p, q) = dsqrt((2.d0*dble(p) + 5.d0)*(2.d0*dble(q) + 5.d0))/(dble(p + q) + 5.d0)
          !trail version:
          ! S(p,q)=dsqrt((2.d0*dble(p)+5.d0)*(2.d0*dble(p)+6.d0)*(2.d0*dble(p)+7.d0)*  &
          !              (2.d0*dble(q)+5.d0)*(2.d0*dble(q)+6.d0)*(2.d0*dble(q)+7.d0) ) &
          !              /                                                             &
          !              ( (5.d0+dble(p)+dble(q)) * (6.d0+dble(p)+dble(q)) * (7.d0+dble(p)+dble(q)) )
        end do
      end do

      nb = ilaenv(1, 'DSYTRD', 'L', n_rbf_so3, -1, -1, -1)
      lwork = (nb + 2)*n_rbf_so3
      call diagsym(S, n_rbf_so3, lwork, L)
      V(:, :) = 0.d0
      do p = 1, n_rbf_so3
        if (L(p) == 0.d0) then
          V(p, p) = 0.d0
        else
          !V(p, p) = L(p)/dabs(L(p))*dabs(L(p))**(-0.5)
          V(p, p) = sign(1.d0, L(p)) * dabs(L(p))**(-0.5)
        end if
        ! Bartok paper
        coeff_rbf_so3(p) = dsqrt((2.d0*p + 5.d0)*r_cut_desc**(-2.d0*p - 5.d0))
        !trial version:
        !coeff_rbf_so3(p) =  dsqrt( 0.5d0 * (2.d0*p+5.d0) * (2.d0*p+6.d0) * (2.d0*p+7.d0) * r_cut**(-2*p-7))

      end do

      W_pow_so3(:, :) = matmul(S(:, :), matmul(V(:, :), transpose(S(:, :))))
    end if

    if (radial_pow_so3 == radial_sgg) then
      if (allocated(chebT)) deallocate (chebT); allocate (chebT(0:n_rbf_so3))
      if (allocated(d_chebT)) deallocate (d_chebT); allocate (d_chebT(0:n_rbf_so3))
      if (allocated(chebU)) deallocate (chebU); allocate (chebU(0:n_rbf_so3))

      if (allocated(scgg)) deallocate (scgg); allocate (scgg(0:n_rbf_so3))
      if (allocated(d_scgg)) deallocate (d_scgg); allocate (d_scgg(0:n_rbf_so3))
    end if

  end subroutine init_pow_so3_rbf



  subroutine spherical_3d(ia_n, factor_ia, factor_ia_3ch, desc_forces_local, dxp, rr)

    use T_kind_param_m, ONLY: double, kind_double, kind_double_complex
    use ml_in_lammps_module, only: l_max, pow_so3_dim, n_rbf_so3, imm_neigh, &
                                weighted, weighted_3ch
    use module_so3, only: dclmn, clmn, coeff_rbf_so3, W_pow_so3, radial_pow_so3, &
                          radial_bartok, radial_sgg, scgg, d_scgg, ini_rbf_so3, &
                          dclmn_w, clmn_w, dclmn_w_3ch, clmn_w_3ch
    use milady_data, only: r_cut_desc

    use angular_functions, only: spherical_harm, grad_spherical_harm
    !use time_check_general, only: debug_time, time, tot_time, debug_time, MY_MPI_WTIME

    implicit none
    real(double), intent(in)    :: rr, factor_ia, factor_ia_3ch
    real(double), dimension(3), intent(in)  :: dxp
    integer, intent(in)  :: ia_n
    logical, intent(in)  :: desc_forces_local
    !complex(kind_double_complex), dimension(-l_max:l_max,pow_so3_dim) :: cmm
    real(double), dimension(:), allocatable :: phi_ji, dphi_ji, rbf_ji, drbf_ji
    real(double)    :: tmp_rr, tmp_rrp, tmp_rbf, tmp_rbf_w, tmp_rbf_w_3ch, &
                        dtmp_rbf, dtmp_rbf_w, dtmp_rbf_w_3ch, tmp_cos(3)
    !------------
    real(double), parameter     :: alpha_gemv = 1.d0, beta_gemv = 0.d0
    integer, parameter   :: incx_gemv = 1, incy_gemv = 1
    !-------------
    integer  :: p, l, i_count, m
    double complex, dimension(0:l_max, -l_max:l_max)   :: sharm
    double complex, dimension(0:l_max, -l_max:l_max, 3)      :: grad_sharm

    if (allocated(phi_ji)) deallocate (phi_ji)
    allocate (phi_ji(n_rbf_so3))
    if (allocated(dphi_ji)) deallocate (dphi_ji)
    allocate (dphi_ji(n_rbf_so3))
    if (allocated(rbf_ji)) deallocate (rbf_ji)
    allocate (rbf_ji(n_rbf_so3))
    if (allocated(drbf_ji)) deallocate (drbf_ji)
    allocate (drbf_ji(n_rbf_so3))

    tmp_cos(1:3) = dxp(1:3)/rr

    if (radial_pow_so3 == radial_bartok) then
      rbf_ji(:) = 0.d0
      drbf_ji(:) = 0.d0
      tmp_rr = r_cut_desc - rr
      do p = 1, n_rbf_so3
        tmp_rrp = tmp_rr**(p + 1)
        phi_ji(p) = coeff_rbf_so3(p)*tmp_rr*tmp_rrp
        if (desc_forces_local) dphi_ji(p) = -coeff_rbf_so3(p)*(p + 2)*tmp_rrp
      end do
      ! y = alpha*A*x + beta * y
      !call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      ! y = \alpha *W phi_ji + \beta y
      !rbf_ji(:) =matmul(W_pow_so3(:,:) , phi_ji(:))
      call dgemv('N', n_rbf_so3, n_rbf_so3, alpha_gemv, W_pow_so3, n_rbf_so3, phi_ji, incx_gemv, beta_gemv, rbf_ji, incy_gemv)
      if (desc_forces_local) then
        !drbf_ji(:)=matmul(W_pow_so3(:,:) ,dphi_ji(:))
        call dgemv('N', n_rbf_so3, n_rbf_so3, alpha_gemv, W_pow_so3, n_rbf_so3, dphi_ji, incx_gemv, beta_gemv, drbf_ji, incy_gemv)
      end if
    end if

    if (radial_pow_so3 == radial_sgg) then
      tmp_rr = 2.d0*rr/r_cut_desc - 1.d0
      call paftouny_so3(desc_forces_local, n_rbf_so3, tmp_rr)
      call sc_stefanodg(desc_forces_local, n_rbf_so3, rr)
      if (allocated(rbf_ji)) deallocate (rbf_ji); allocate (rbf_ji(0:n_rbf_so3))
      if (allocated(drbf_ji)) deallocate (drbf_ji); allocate (drbf_ji(0:n_rbf_so3))
    end if

    ! compute spherical functions
    do l = 0, l_max
      do m = -l, l
        sharm(l, m) = spherical_harm(l, m, dxp(:))
        if (desc_forces_local) grad_sharm(l, m, :) = grad_spherical_harm(l, m, dxp(:))
      end do
    end do


    i_count = 0
    if (radial_pow_so3 == radial_bartok) then
    do p = ini_rbf_so3, n_rbf_so3
      tmp_rbf = rbf_ji(p)
      tmp_rbf_w = factor_ia * tmp_rbf
      tmp_rbf_w_3ch = factor_ia_3ch * tmp_rbf
      dtmp_rbf = drbf_ji(p)
      dtmp_rbf_w = factor_ia * dtmp_rbf
      dtmp_rbf_w_3ch = factor_ia_3ch * dtmp_rbf
      do l = 0, l_max
        i_count = i_count + 1
        do m = -l, l
          clmn(m, i_count) = clmn(m, i_count) + tmp_rbf*sharm(l, m)
          if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = dtmp_rbf*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf*grad_sharm(l, m, 1:3)
          if (weighted) then
            clmn_w(m, i_count) = clmn_w(m, i_count) + tmp_rbf_w*sharm(l, m)
            if (desc_forces_local) dclmn_w(m, i_count, ia_n, 1:3) = dtmp_rbf_w*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w*grad_sharm(l, m, 1:3)
            if (weighted_3ch) then
              clmn_w_3ch(m, i_count) = clmn_w_3ch(m, i_count) + tmp_rbf_w_3ch*sharm(l, m)
              if (desc_forces_local) dclmn_w_3ch(m, i_count, ia_n, 1:3) = dtmp_rbf_w_3ch*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w_3ch*grad_sharm(l, m, 1:3)

            end if
          end if
        end do                  ! m
      end do                  ! l
    end do                  ! p
    end if


    if (radial_pow_so3 == radial_sgg) then
    do p = ini_rbf_so3, n_rbf_so3
      tmp_rbf = scgg(p)
      tmp_rbf_w = tmp_rbf * factor_ia
      tmp_rbf_w_3ch =  tmp_rbf * factor_ia_3ch
      dtmp_rbf = d_scgg(p)
      dtmp_rbf_w = dtmp_rbf * factor_ia
      dtmp_rbf_w_3ch = dtmp_rbf * factor_ia_3ch
      do l = 0, l_max
        i_count = i_count + 1
        do m = -l, l
          clmn(m, i_count) = clmn(m, i_count) + tmp_rbf*sharm(l, m)
          if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = dtmp_rbf*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf*grad_sharm(l, m, 1:3)
          if (weighted) then
            clmn_w(m, i_count) = clmn_w(m, i_count) + tmp_rbf_w*sharm(l, m)
            if (desc_forces_local) dclmn_w(m, i_count, ia_n, 1:3) = dtmp_rbf_w*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w*grad_sharm(l, m, 1:3)
            if (weighted_3ch) then
              clmn_w_3ch(m, i_count) = clmn_w_3ch(m, i_count) + tmp_rbf_w_3ch*sharm(l, m)
              if (desc_forces_local) dclmn_w_3ch(m, i_count, ia_n, 1:3) = dtmp_rbf_w_3ch*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w_3ch*grad_sharm(l, m, 1:3)

            end if
          end if
        end do                  ! m
      end do                  ! l
    end do                  ! p
    end if


  end subroutine spherical_3d



  subroutine paftouny_so3(desc_forces_local, n_cheb, xx)

    use T_kind_param_m, ONLY: kind_double
    use module_so3, only: chebT, chebU, d_chebT
    implicit none

    logical, intent(in)  :: desc_forces_local
    real(kind_double), intent(in)    :: xx
    integer, intent(in)  :: n_cheb
    integer  :: nn

    chebT(0) = 1.d0
    chebT(1) = xx
    chebU(0) = 1.d0
    chebU(1) = 2.d0*xx

    d_chebT(0) = 0.d0
    d_chebT(1) = 1.d0
    if (n_cheb .ge. 2) then
      do nn = 2, n_cheb
        chebT(nn) = xx*chebT(nn - 1) - (1.d0 - xx**2)*chebU(nn - 2)
        chebU(nn) = xx*chebU(nn - 1) + chebT(nn)
      end do

      if (desc_forces_local) then
        do nn = 2, n_cheb
          d_chebT(nn) = dble(nn)*chebU(nn - 1)
        end do
      end if
    end if
  end subroutine paftouny_so3



  subroutine sc_stefanodg(desc_forces_local, n_cheb, rr)

    use T_kind_param_m, ONLY: kind_double
    use ml_in_lammps_module, only: pi
    use module_so3, only: chebT, d_chebT, scgg, d_scgg
    use milady_data, only: r_cut_desc
    implicit none

    integer, intent(in)  :: n_cheb
    logical, intent(in)  :: desc_forces_local
    real(kind_double), intent(in)    :: rr
    real(kind_double)    :: xx
    integer  :: nn

    xx = pi*(rr/r_cut_desc)
    scgg(0) = 0.d0
    scgg(1) = 0.5d0*(1.d0 + dcos(xx))

    d_scgg(0) = 0.d0
    d_scgg(1) = -0.5d0*dsin(xx)*pi/r_cut_desc

    if (n_cheb > 2) then
      do nn = 2, n_cheb
        scgg(nn) = 0.5d0*(1.d0 - chebT(nn - 1))*scgg(1)
        if (desc_forces_local) then
          ! THIS SHOULD BE x 2.D0/rcut
          d_scgg(nn) = -1.d0*scgg(1)*d_chebT(nn - 1)/r_cut_desc + 0.5d0*(1.d0 - chebT(nn - 1))*d_scgg(1)
        end if
      end do
    end if
  end subroutine sc_stefanodg
