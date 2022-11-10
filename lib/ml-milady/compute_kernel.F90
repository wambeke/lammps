
module module_compute_kernel
  !use mpi
  use T_kind_param_m, ONLY:  double
  use ml_in_lammps_module, ONLY: rangml,timetot, timetot00, timetot01, timetot02, timetot03, one_pi,  imm_neigh,    &
                               weighted, weighted_3ch,  desc_forces

  use milady_data, only: local_desc, local_desc_deriv, local_desc_kernel, local_desc_kernel_deriv, &
                         r_cut_desc, n_desc_type, k_desc_type
  use module_kernel, only: dim_kernel, global_kernel, draft_kernel, &
                           length_kernel, sigma_kernel, &
                           kernel_type, kernel_se, kernel_po, kernel_random, kernel_random_po, kernel_mcd, kernel_maha, &
                           min_ker, max_ker, var_ker, mean_ker, kernel_power, &
                           kernel_phase_random, &
                           maha_kernel, maha_norm_kernel, basis_random_po
  use temporary_data_cov, only: dim_xdesc, dim_xdesc_kernel
  implicit none

  contains


  subroutine  compute_kernel(i_start_at, i_final_at)
    implicit none
    integer, intent (in) :: i_start_at,i_final_at
    real(double), dimension(:,:), allocatable :: vfk_trans1, vfk_trans2, vfk_trans3
    real(double), dimension(:,:), allocatable :: force_tmp1, force_tmp2, force_tmp3
    real(double), dimension(:), allocatable :: vtemp1, vtemp2, vtemp3, vtemp4, &
                                           vtemp_ja1, vtemp_ja2, vtemp_ja3, &
                                           energy_ja, kernel_ik
    real(double), dimension(3) ::  dktemp, dktemp_ja, dktemp_norm, &
                                   vtemp, vtemp_ja, vtemp_ik, vtemp0
    real(double) :: dtmp_jaik, ktemp_norm
    real(double) :: l2kse, sigma_kse2, soverl
    integer :: j, iw, ja, inn, ja_neigh, ik, ii, iid, dim_omega
    logical :: debug_time
    real(double) :: dtmp, dtmp_ja, ktemp, k_ja, k_ja2, k_ja32, dk_ja, k_ik, k_ik2, tmp_f1
    real(double) :: time01, time01i, time02, time03, time04,  time03i
    real(double), dimension(:), allocatable :: ene_tmp

    if (allocated(energy_ja)) deallocate(energy_ja) ; allocate(energy_ja(dim_xdesc))
    if (allocated(kernel_ik)) deallocate(kernel_ik) ; allocate(kernel_ik(dim_xdesc))

    local_desc_kernel=0.d0
    if (desc_forces) local_desc_kernel_deriv(:,:,:)=0.d0
    if ((i_start_at==0) .and. (i_final_at==0)) return

    l2kse = 2.d0*length_kernel**2
    sigma_kse2 = sigma_kernel**2
    soverl=sigma_kernel**2/length_kernel**2
    do ja=i_start_at,i_final_at
      energy_ja(:) =  local_desc(:)

      if(kernel_type == kernel_po) then
        dtmp_ja=dot_product(energy_ja(:), energy_ja(:))
        k_ja = (sigma_kse2 + dtmp_ja/l2kse)**kernel_power
        k_ja2 = dsqrt(k_ja)
        k_ja32 = k_ja2**3
        dk_ja = 2.d0*dble(kernel_power)*(sigma_kse2 + dtmp_ja/l2kse)**(kernel_power-1)/l2kse
      end if

      if (desc_forces) then
        ja_neigh = n_desc_type
        if (allocated(vtemp1)) deallocate(vtemp1) ; allocate(vtemp1(ja_neigh))
        if (allocated(vtemp2)) deallocate(vtemp2) ; allocate(vtemp2(ja_neigh))
        if (allocated(vtemp3)) deallocate(vtemp3) ; allocate(vtemp3(ja_neigh))
        if (allocated(vtemp4)) deallocate(vtemp4) ; allocate(vtemp4(ja_neigh))
        !
        if (allocated(vfk_trans1)) deallocate(vfk_trans1) ; allocate(vfk_trans1(ja_neigh, dim_xdesc))
        if (allocated(vfk_trans2)) deallocate(vfk_trans2) ; allocate(vfk_trans2(ja_neigh, dim_xdesc))
        if (allocated(vfk_trans3)) deallocate(vfk_trans3) ; allocate(vfk_trans3(ja_neigh, dim_xdesc))
        do inn =1, ja_neigh
          !TODO not forget There are those alternatives ...
          ! fk_trans(:) =(config_desc(iconf)%force(:,ja,inn,ix) - min_ker(:))/(max_ker(:)-min_ker(:)) + min_ker(:)
          !  fk_trans(:) =(config_desc(iconf)%force(:,ja,inn,ix) - mean_ker(:))/var_ker(:)
          vfk_trans1(inn,:) = local_desc_deriv(:,inn,1)
          vfk_trans2(inn,:) = local_desc_deriv(:,inn,2)
          vfk_trans3(inn,:) = local_desc_deriv(:,inn,3)
        end do


        if (kernel_type == kernel_po) then
          if (allocated(vtemp_ja1)) deallocate(vtemp_ja1) ; allocate(vtemp_ja1(ja_neigh))
          if (allocated(vtemp_ja2)) deallocate(vtemp_ja2) ; allocate(vtemp_ja2(ja_neigh))
          if (allocated(vtemp_ja3)) deallocate(vtemp_ja3) ; allocate(vtemp_ja3(ja_neigh))
          call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, energy_ja,1, 0.d0, vtemp_ja1,1)
          call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, energy_ja,1, 0.d0, vtemp_ja2,1)
          call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, energy_ja,1, 0.d0, vtemp_ja3,1)
        end if

      end if


      do ik=1, dim_kernel
        kernel_ik(:) = global_kernel(:,ik)

        select case (kernel_type)

          !TODO - speed of matrix multiplication ...
          case (kernel_random)

            dtmp_jaik = dot_product(energy_ja(:),kernel_ik(:)) + kernel_phase_random(ik)
            ktemp_norm = cos (dtmp_jaik)
            if (desc_forces) then
              vtemp0(:)=0.d0
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik,1, 0.d0, vtemp1,1)
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik,1, 0.d0, vtemp2,1)
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik,1, 0.d0, vtemp3,1)
              do inn=1, ja_neigh
                dktemp_norm(1) = - vtemp1(inn)* sin (dtmp_jaik)
                dktemp_norm(2) = - vtemp2(inn)* sin (dtmp_jaik)
                dktemp_norm(3) = - vtemp3(inn)* sin (dtmp_jaik)
                vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
                !config_desc(iconf)%force_kernel(ik,ja,inn,1:3) local_desc_kernel_deriv (ik,inn,1:3)= dktemp_norm(1:3)
                local_desc_kernel_deriv (ik,inn,1:3) =    dktemp_norm(1:3)
              end do
              !config_desc(iconf)%force_kernel(ik,ja,0,1:3) = vtemp0(1:3)
              local_desc_kernel_deriv(ik,0,1:3) =    vtemp0(1:3)
            end if

          case(kernel_po)

            dtmp_jaik = dot_product(energy_ja(:), kernel_ik(:))
            dtmp=sigma_kse2+dtmp_jaik/l2kse
            ktemp=dtmp**kernel_power
            k_ik=(sigma_kse2 + dot_product(kernel_ik(:), kernel_ik(:))/l2kse)**kernel_power
            k_ik2=dsqrt(k_ik)
            ktemp_norm = ktemp / (k_ik2*k_ja2)
            if (desc_forces) then
              vtemp0(:)=0.d0
              !
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik,1, 0.d0, vtemp1,1)
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik,1, 0.d0, vtemp2,1)
              call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik,1, 0.d0, vtemp3,1)
              !
              tmp_f1= dble(kernel_power)* dtmp**(kernel_power-1)/l2kse
              do inn=1, ja_neigh
                dktemp(1)= tmp_f1 * vtemp1(inn)
                dktemp(2)= tmp_f1 * vtemp2(inn)
                dktemp(3)= tmp_f1 * vtemp3(inn)
                dktemp_ja(1)= vtemp_ja1(inn) * dk_ja
                dktemp_ja(2)= vtemp_ja2(inn) * dk_ja
                dktemp_ja(3)= vtemp_ja3(inn) * dk_ja
                dktemp_norm(1:3)= (dktemp(1:3) - ktemp * dktemp_ja(1:3) / (2.d0*k_ja))/(k_ja2*k_ik2)
                local_desc_kernel_deriv(ik,inn,1:3)  = dktemp_norm(1:3)
                vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
              end do
              local_desc_kernel_deriv(ik,0,1:3) = vtemp0(1:3)
            end if ! desc_forces_local

          case (kernel_random_po)

            dim_omega =  basis_random_po(ik)%dim_omega
            if (allocated(ene_tmp)) deallocate(ene_tmp) ; allocate(ene_tmp(dim_omega))
            dtmp = 1.d0
            do ii = 1, dim_omega
              ene_tmp(ii) =  dot_product(energy_ja(:),basis_random_po(ik)%omega(:,ii))
              dtmp = dtmp * ene_tmp(ii)
            end do
            ktemp_norm = dtmp / sqrt(dble(dim_kernel))
            if (desc_forces) then
              vtemp0(1:3) = 0.d0
              vtemp1(1:ja_neigh) = 0.d0
              vtemp2(1:ja_neigh) = 0.d0
              vtemp3(1:ja_neigh) = 0.d0
              if (allocated(force_tmp1)) deallocate(force_tmp1) ; allocate(force_tmp1(ja_neigh, dim_omega))
              if (allocated(force_tmp2)) deallocate(force_tmp2) ; allocate(force_tmp2(ja_neigh, dim_omega))
              if (allocated(force_tmp3)) deallocate(force_tmp3) ; allocate(force_tmp3(ja_neigh, dim_omega))
              do ii = 1, dim_omega
                kernel_ik(:) =  basis_random_po(ik)%omega(:,ii)
                call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans1, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
                force_tmp1(1:ja_neigh, ii) = vtemp4(1:ja_neigh)
                call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans2, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
                force_tmp2(1:ja_neigh,ii) = vtemp4(1:ja_neigh)
                call dgemv('N', ja_neigh, dim_xdesc, 1.d0, vfk_trans3, ja_neigh, kernel_ik, 1, 0.d0, vtemp4, 1)
                force_tmp3(1:ja_neigh,ii) = vtemp4(1:ja_neigh)
              end do
              do iid = 1, dim_omega
                dtmp = 1.d0
                do ii = 1, dim_omega
                  if (ii == iid ) cycle
                  dtmp = dtmp * ene_tmp(ii)
                end do
                vtemp1(1:ja_neigh) = vtemp1(1:ja_neigh) + dtmp * force_tmp1(1:ja_neigh, iid)
                vtemp2(1:ja_neigh) = vtemp2(1:ja_neigh) + dtmp * force_tmp2(1:ja_neigh, iid)
                vtemp3(1:ja_neigh) = vtemp3(1:ja_neigh) + dtmp * force_tmp3(1:ja_neigh, iid)
              end do

              do inn = 1, ja_neigh
                dktemp_norm(1) = vtemp1(inn)
                dktemp_norm(2) = vtemp2(inn)
                dktemp_norm(3) = vtemp3(inn)
                vtemp0(:) = vtemp0(:) - dktemp_norm(1:3)
                local_desc_kernel_deriv (ik,inn,1:3) = dktemp_norm(1:3)/ sqrt(dble(dim_kernel))
              end do

              local_desc_kernel_deriv(ik,0,1:3) = vtemp0(1:3)/ sqrt(dble(dim_kernel))
            end if


          end select
          local_desc_kernel(ik) = ktemp_norm
        end do

    end do ! ja

  return
  end subroutine compute_kernel
end module module_compute_kernel
