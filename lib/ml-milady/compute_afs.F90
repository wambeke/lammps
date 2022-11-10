!     AFS computing the trinagles centered on the atom j and the angle ijk.
!
!                 j--------i
!                  \
!                   \
!                    \
!                     k
!used quantities ji, kj and jik



subroutine  compute_afs(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
USE T_kind_param_m, ONLY:  double
use ml_in_lammps_module, ONLY: n_rbf_afs,n_cheb,  weighted,  &
                               W_afs, desc_forces, factor_weight_mass, coeff_rbf_afs, afs_type, afs_type_bartok,  &
                               afs_type_homemade, weighted, weighted_3ch, afs_dim
!mi use derived_types, only : config_real
use milady_data, only: local_afs => local_desc, local_afs_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type, &
                                    weight_per_type, weight_per_type_3ch
 implicit none
integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
integer, dimension(numneigh_full), intent(in) :: firstneigh_full(numneigh_full)
integer, dimension(nmax), intent(in)  :: a_type
real(double), dimension(ntype), intent(in) :: map_mass
real(double), dimension(3,nmax), intent(in) :: xp



real(double), dimension(3) :: dxp_ji,dxp_jk, dxp_ik
real(double), dimension(3) :: cdxp_ji,cdxp_jk, cdxp_ik
real(double), dimension(n_rbf_afs) :: phi_ji, dphi_ji, phi_jk, dphi_jk, rbf_ji, rbf_jk, drbf_ji, drbf_jk
integer :: ia,ja,ka, ia_n, ka_n, iw, i_desc, i_desc2, i_desc3

integer :: iz,p1,p2, p0, pmin, pmax
integer :: ia_type, ja_type, ka_type
double precision :: r2_ji,r2_jk,r_ji,r_jk, r2_ik, r_ik, cos_kji, tmpr, term_local
double precision,dimension(0:n_cheb) :: chebT_kji,chebU_kji
double precision,dimension(3) :: dcos_kji_ji,dcos_kji_jk,dchebT_kji_ji,dchebT_kji_jk, term_local_3d
double precision :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch, factor_ka, factor_ka_3ch



if (allocated(k_desc_type)) deallocate(k_desc_type) ; allocate(k_desc_type(numneigh_full))
local_afs(:)=0.d0
local_afs_deriv(:,:,:)=0.d0
k_desc_type(:)=0

if ( (i_start_at==0) .and. (i_final_at==0) ) return


do ja=i_start_at,i_final_at

   ja_type = a_type(ja)
   if (weighted) then
      factor_ja= weight_per_type(ja_type)          !CHANGE config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
      if (weighted_3ch) factor_ja_3ch =  weight_per_type_3ch(ja_type)    !CHANGE config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
   else
      factor_ja=1.d0
   endif


   ia_n = 0

   do iw=1, numneigh_full
      ia = firstneigh_full(iw)
      ia_type = a_type(ia)
      if (weighted) then
         factor_ia =  weight_per_type(ia_type) !map_mass(a_type(ia)) / factor_weight_mass
         if (weighted_3ch) factor_ia_3ch = weight_per_type_3ch(ia_type)
      else
         factor_ia=1.d0
      endif

      dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
      r2_ji = Sum( dxp_ji(1:3)**2 )
      r_ji = dsqrt(r2_ji)

      if (r_ji >= r_cut_desc) cycle
      phi_ji(:)=0.d0
      ia_n = ia_n + 1

      do p1=1,n_rbf_afs
         !phi_ji(p1) = dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(r_cut-r_ji)**(p1+2.d0)
         tmpr= (r_cut_desc-r_ji)**(p1+1)
         phi_ji(p1) = coeff_rbf_afs(p1)*tmpr*(r_cut_desc-r_ji)
         !phi_ji(p1) = coeff_rbf_afs(p1)*(r_cut-r_ji)**(p1+2.d0)
         !if (desc_forces_local) dphi_ji(p1)=-coeff_rbf_afs(p1)*(p1+2.d0)*(r_cut-r_ji)**(p1+1.d0)
         dphi_ji(p1)=-coeff_rbf_afs(p1)*(p1+2)*tmpr
      enddo
      rbf_ji(:) =matmul(W_afs(:,:) , phi_ji(:))
      drbf_ji(:)=matmul(W_afs(:,:) ,dphi_ji(:))
      cdxp_ji(:)=dxp_ji(:)/r_ji

      k_desc_type(ia_n)=ia
      ka_n=0
      do iz=1, numneigh_full
         ka = firstneigh_full(iz)
         ka_type = a_type(ka)
         if (weighted) then
            factor_ka = weight_per_type(ka_type)
            if (weighted_3ch) factor_ka_3ch = weight_per_type_3ch(ka_type)
         else
            factor_ka=1.d0
         endif

         dxp_jk(1:3) = xp(1:3,ka) - xp(1:3,ja)
         r2_jk = Sum( dxp_jk(1:3)**2 )
         r_jk = dsqrt(r2_jk)

         dxp_ik(1:3) = xp(1:3,ka) - xp(1:3,ia)
         r2_ik = Sum( dxp_ik(1:3)**2 )
         r_ik = dsqrt(r2_ik)

         if (r_jk >= r_cut_desc) cycle
         ka_n = ka_n + 1
         if (r_ik <= 1d-20) cycle

         do p1=1,n_rbf_afs
            tmpr= (r_cut_desc-r_jk)**(p1+1)
            !phi_ji(p1) = coeff_rbf_afs(p1)*(r_cut_desc-r_ji)**(p1+2.d0)
            phi_jk(p1) = coeff_rbf_afs(p1)*(r_cut_desc-r_jk)*tmpr
            !if (desc_forces) dphi_ji(p1)=-coeff_rbf_afs(p1)*(p1+2.d0)*(r_cut_desc-r_ji)**(p1+1.d0)
            dphi_jk(p1)=-coeff_rbf_afs(p1)*dble(p1+2)*tmpr
         enddo

         rbf_jk(:) =matmul(W_afs(1:n_rbf_afs,1:n_rbf_afs) , phi_jk(1:n_rbf_afs))
         drbf_jk(:)=matmul(W_afs(1:n_rbf_afs,1:n_rbf_afs) ,dphi_jk(1:n_rbf_afs))

         cdxp_jk(:)=dxp_jk(:)/r_jk
         cdxp_ik(:)=dxp_ik(:)/r_ik
         cos_kji=dot_product(cdxp_ji(1:3),cdxp_jk(1:3))

         if (desc_forces) dcos_kji_ji(1:3)=(cdxp_jk(1:3) - cos_kji*cdxp_ji(1:3))/r_ji
         if (desc_forces) dcos_kji_jk(1:3)=(cdxp_ji(1:3) - cos_kji*cdxp_jk(1:3))/r_jk

         chebT_kji(0)=1.d0
         chebT_kji(1)=cos_kji
         chebU_kji(0)=1.d0
         chebU_kji(1)=2.d0*cos_kji
         if (n_cheb.ge.2) then
              do p2=2,n_cheb
                 chebT_kji(p2)=2d0*cos_kji*chebT_kji(p2-1)-chebT_kji(p2-2)
                 chebU_kji(p2)=2d0*cos_kji*chebU_kji(p2-1)-chebU_kji(p2-2)
              enddo
         endif

         i_desc=0
         do p0=1,n_rbf_afs
            if (afs_type == afs_type_bartok) then
               pmin=p0
               pmax=p0
            else
               pmin=1
               pmax=n_rbf_afs
            end if

            do p1=pmin, pmax
               do p2=0,n_cheb

                  if (p2==0) then
                    dchebT_kji_ji(1:3)=0d0
                    dchebT_kji_jk(1:3)=0d0
                  else
                    dchebT_kji_ji(1:3)=p2*chebU_kji(p2-1)*dcos_kji_ji(1:3)
                    dchebT_kji_jk(1:3)=p2*chebU_kji(p2-1)*dcos_kji_jk(1:3)
                  endif

                  i_desc = i_desc + 1
                  tmpr = rbf_ji(p0) * rbf_jk(p1)
                  !local_afs(i_desc) = local_afs(i_desc) + factor_ia*factor_ka* rbf_ji(p0) * rbf_jk(p1) * chebT_kji(p2)/2.d0
                  term_local= tmpr * chebT_kji(p2)/2.d0
                  local_afs(i_desc) = local_afs(i_desc) + term_local
                  if (weighted) then
                     i_desc2 = i_desc + afs_dim
                     i_desc3 = i_desc2 + afs_dim
                     local_afs(i_desc2) = local_afs(i_desc2) +  factor_ia*factor_ka*term_local
                     if (weighted_3ch) local_afs(i_desc3) = local_afs(i_desc3) +  factor_ia_3ch*factor_ka_3ch*term_local
                  end if

                  !if (desc_forces) local_afs_deriv(i_desc,ia_n,1:3) = local_afs_deriv(i_desc, ia_n, 1:3) + factor_ia* factor_ka* ( rbf_ji(p0)*rbf_jk(p1)*dchebT_kji_ji(1:3) + &
                  !                                        drbf_ji(p0)*cdxp_ji(1:3)*rbf_jk(p1)*chebT_kji(p2) )
                  if (desc_forces) then
                        term_local_3d(1:3)= ( rbf_ji(p0)*rbf_jk(p1)*dchebT_kji_ji(1:3) + &
                                              drbf_ji(p0)*cdxp_ji(1:3)*rbf_jk(p1)*chebT_kji(p2) ) !/2.d0
                        local_afs_deriv(i_desc,ia_n,1:3) = local_afs_deriv(i_desc, ia_n, 1:3) + term_local_3d(1:3)
                        if (weighted) then
                           local_afs_deriv(i_desc2,ia_n,1:3) = local_afs_deriv(i_desc2,  ia_n, 1:3) + &
                                                               factor_ia*factor_ka*term_local_3d(1:3)
                           if (weighted_3ch) local_afs_deriv(i_desc3,ia_n,1:3) = &
                                                  local_afs_deriv(i_desc3, ia_n, 1:3) + &
                                                  factor_ia_3ch*factor_ka_3ch*term_local_3d(1:3)
                        end if
                  end if
                  !
               enddo   !p2 cos
            enddo      !p1
         end do        !p0
      end do ! ka of the  neigh of ja

      if (desc_forces) then
         local_afs_deriv(:, 0,:) =    local_afs_deriv(:,0,:) -  local_afs_deriv(:,ia_n,:)
      end if
   end do  ! ia_n loop over neigh of ja

   if (weighted) then
      local_afs(1+afs_dim:2*afs_dim) = local_afs(1+afs_dim:2*afs_dim)*factor_ja
      if (weighted_3ch) local_afs(2*afs_dim+1:3*afs_dim) = local_afs(2*afs_dim+1:3*afs_dim)*factor_ja_3ch
   end if
   if (desc_forces)  then
      local_afs_deriv(1:afs_dim,0:ia_n,:) =  local_afs_deriv(1:afs_dim,0:ia_n,:)
      if (weighted) then
         local_afs_deriv(1+afs_dim:2*afs_dim,0:ia_n,:) = local_afs_deriv(1+afs_dim:2*afs_dim,0:ia_n,:)*factor_ja
         if (weighted_3ch) local_afs_deriv(2*afs_dim+1:3*afs_dim,0:ia_n,:) = &
                              local_afs_deriv(2*afs_dim+1:3*afs_dim,0:ia_n,:)*factor_ja_3ch
      end if
   end if

   n_desc_type = ia_n
 end do   !ja main loop


 !mi deallocate(xpnp)

return
end subroutine compute_afs



subroutine init_afs_rbf()
use ml_in_lammps_module, only: afs_dim, n_rbf_afs, n_cheb, W_afs, coeff_rbf_afs, &
                               afs_type, afs_type_bartok, afs_type_homemade
use milady_data, only : r_cut_desc
implicit none
real(kind(0.d0)),dimension(n_rbf_afs,n_rbf_afs) :: S,V
real(kind(0.d0)),dimension(n_rbf_afs) :: L
!local
integer :: p, q, nb, ilaenv,lwork

if (allocated(W_afs)) deallocate(W_afs) ; allocate(W_afs(n_rbf_afs,n_rbf_afs))
if (allocated(coeff_rbf_afs)) deallocate(coeff_rbf_afs) ; allocate(coeff_rbf_afs(n_rbf_afs))

!write(*,*) 'in inti afs ...', afs_type, n_rbf, n_cheb
select case(afs_type)
    case (afs_type_bartok)
        afs_dim = n_rbf_afs*(n_cheb+1)
    case (afs_type_homemade)
        afs_dim = n_rbf_afs**2*(n_cheb+1)
    case default
        write(*,*) 'Milady: Fatal no implementation for this afs_type. ONLY 1 and 2'
        stop
end select


do q=1,n_rbf_afs
   do p=1,n_rbf_afs
      S(p,q)=dsqrt((2.d0*dble(p)+5.d0)*(2.d0*dble(q)+5.d0))/(dble(p+q)+5.d0)
   enddo
enddo

 nb = ilaenv( 1, 'DSYTRD', 'L', n_rbf_afs, -1, -1, -1 )
 lwork=(nb+2)*n_rbf_afs
 call diagsym(S,n_rbf_afs,lwork,L)
 V(:,:)=0d0
 do p=1,n_rbf_afs
    !if (abs(L(p) - 0.d0).lt.1.d-15)  then
    if (L(p)==0.d0) then
       V(p,p) = 0.d0
    else
       V(p,p) = L(p)/dabs(L(p)) * dabs(L(p))**(-0.5)
    end if

    coeff_rbf_afs(p) = dsqrt((2.d0*p+5.d0)*r_cut_desc**(-2.d0*p-5.d0))
 enddo

W_afs(:,:)=matmul(S(:,:),matmul(V(:,:),transpose(S(:,:))))
!W_pow_so3(:,:)=matmul( matmul(S(:,:), V(:,:)),TRANSPOSE(S(:,:)))

return
end subroutine init_afs_rbf


subroutine diagsym(A,n_rbf,lwork,L)
! input - the symmetric matrix A
! output - the ortho-normalized vectors that diaonalize A and the eigenvalues L
use ml_in_lammps_module, only: rangml
implicit none
integer,intent(in) :: n_rbf,lwork
double precision,dimension(n_rbf,n_rbf),intent(inout) :: A
double precision,dimension(n_rbf),intent(out) :: L
double precision,dimension(lwork) :: work
integer :: info

call dsyev('V','L',n_rbf,A,n_rbf,L,work,lwork,info)

if (rangml==0) then
  if (.not.(info == 0)) then
   write(6,*) 'ML: WARNING the nomalization of the overlap matrix is WRONG in diagsym compute_afs.F90'
  end if
end if

return
end subroutine diagsym
