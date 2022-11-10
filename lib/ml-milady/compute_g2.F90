!subroutine compute_g2(i_start_at,i_final_at, local_d_n_neigh, local_d_kind_neigh,  local_g2_out,local_g2_deriv, iconf)
subroutine compute_g2(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
USE T_kind_param_m, ONLY:  double
!use gen_com_m, ONLY: imm,A2cm,lperiod,bg,at,indi2
!use tab_imm_m, ONLY : iwmax2,xp
use ml_in_lammps_module, ONLY: one_pi, weighted, &
                            g2_eta, g2_rs, &
                            g2_dim, factor_weight_mass, desc_forces
use milady_data, only: local_g2 => local_desc, local_g2_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type
!use derived_types, only: config_real
implicit none
integer, intent (in) :: i_start_at,i_final_at,nmax, ntype, numneigh_full
integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
integer, dimension(nmax), intent(in)  :: a_type
real(double), dimension(ntype), intent(in) :: map_mass
real(double), dimension(3,nmax), intent(in) :: xp

real(double), dimension(3) :: dxp_ji
integer :: iw

integer :: p
real(double) :: fcut,r_ji,r2_ji, dfcut
real(double) :: factor_ia, factor_ja, g2_func

integer :: ia,ja,ia_n


if (allocated(k_desc_type)) deallocate(k_desc_type) ; allocate(k_desc_type(numneigh_full))
local_g2(:)=0.d0
local_g2_deriv(:,:,:)=0.d0
k_desc_type(:)=0

if ( (i_start_at==0) .and. (i_final_at==0) ) return


  do ja=i_start_at,i_final_at

    if (weighted) then
         factor_ja = map_mass(a_type(ja)) / factor_weight_mass
    else
        factor_ja=1.d0
    endif
    ia_n = 0

    do iw=1, numneigh_full
        ia = firstneigh_full(iw)
        if (weighted) then
           factor_ia = map_mass(a_type(ia)) / factor_weight_mass
        else
           factor_ia=1.d0
        endif

        dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
        r2_ji = Sum( dxp_ji(1:3)**2 )
        r_ji = dsqrt(r2_ji)

        if (r_ji >= r_cut_desc) cycle
        ia_n = ia_n + 1
        k_desc_type(ia_n)=ia

        fcut=0.5d0*(cos(one_pi*r_ji/r_cut_desc)+1.d0)*factor_ia
        dfcut=-0.5d0*one_pi*sin(one_pi*r_ji/r_cut_desc)/r_cut_desc*factor_ia
        !debug write(*,*) fcut, dfcut, factor_ia
        do p=1,g2_dim
              g2_func = dexp(-g2_eta(p)*(r_ji-g2_rs(p))**2)
              local_g2(p) = local_g2(p) + fcut * g2_func
              if (desc_forces) local_g2_deriv(p, ia_n, 1:3) = &
                local_g2_deriv(p, ia_n, 1:3) + (-2d0*g2_eta(p)*(r_ji-g2_rs(p))*fcut + dfcut) * g2_func*dxp_ji(1:3)/r_ji
              if (desc_forces) local_g2_deriv(p, 0, 1:3) = local_g2_deriv(p, 0, 1:3) - local_g2_deriv(p, ia_n, 1:3)
        enddo
     end do  !ia
  local_g2(:) = local_g2(:) * factor_ja!/dble(ia_n)
  if (desc_forces) local_g2_deriv(:,:,:) = local_g2_deriv(:,:,:) * factor_ja!/dble(ia_n)
  n_desc_type = ia_n
  end do  !ja

return
end subroutine compute_g2
