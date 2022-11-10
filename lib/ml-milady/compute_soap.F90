subroutine  compute_soap(i_start_at, i_final_at,  nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
use T_kind_param_m, ONLY:  double
!mi use gen_com_m, ONLY: imm
use angular_functions, only : spherical_harm, grad_spherical_harm
use ml_in_lammps_module, ONLY: rangml,l_max,  imm_neigh,    &
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


complex(kind(0.d0)), dimension(-l_max:l_max,0:l_max,n_soap) :: c
complex(kind(0.d0)), dimension(-l_max:l_max,0:l_max, n_soap, imm_neigh,3) :: dc
real(double), dimension(soap_dim) :: tmp_soap
real(double), dimension(soap_dim, 0:imm_neigh,3):: tmp_soap_deriv
real(double), dimension(soap_dim, 3):: d_tmp_soap
real(double), dimension (0:l_max, n_soap) :: phi_ji, dphi_ji, rbf_ji, drbf_ji
integer :: p, l, m
complex(kind(0.d0)) :: czero
real(double), dimension(3) :: dxp_ji
integer :: iw

integer ::  p1, p2, i_s1, i_s2, i_tmp, j_tmp
integer :: ia, ja, ix, ia_n, i_desc
double precision :: r2_ji,r_ji, norm_descriptor
double precision :: factor_ia, factor_ja, fcut, dfcut, x_temp
real(double) :: x_b !mi , ms_bessel(0:l_max), d_ms_bessel(0:l_max)
double precision :: bessel_l, bessel_lm, bessel_lp, bessel_lmm, exp_m, exp_p

if (rangml ==0) then
   write(*,*) '...in soap'
end if


if (allocated(k_desc_type)) deallocate(k_desc_type) ; allocate(k_desc_type(numneigh_full))
local_soap(:)=0.d0
local_soap_deriv(:,:,:)=0.d0
k_desc_type(:)=0

if ((i_start_at==0) .and. (i_final_at==0)) return

czero=cmplx(0.d0,0.d0, kind=kind(1.d0))


do ja=i_start_at,i_final_at
        !debug write(*,*) 'in the core 001 ', ja

    if (weighted) then
         factor_ja = map_mass(a_type(ja)) / factor_weight_mass
    else
        factor_ja=1.d0
    endif
    ia_n = 0



    tmp_soap(:)=0.d0
    tmp_soap_deriv(:,:,:)=0.d0
    c(:,:,:)=(0.d0,0.d0)
    dc(:,:,:,:,:)=(0.d0,0.d0)

    phi_ji(0,:) = 0.0d0
    phi_ji(0,1) = 1.0d0
    rbf_ji(0,:) = matmul( phi_ji(0,:), S_factor_matrix(:,:))

    do p = 1, n_soap
      c(0,0,p) = rbf_ji(0,p)*spherical_harm(0,0,vec_3d_zero(:))
      do l = 1, l_max
        c(:,l,p)=c_zero
      enddo
    enddo


   do iw=1, numneigh_full
        ia = firstneigh_full(iw)
        !debug write(*,*) 'in the core xyz ', ja, ia
        if (weighted) then
           factor_ia = map_mass(a_type(ia)) / factor_weight_mass
        else
           factor_ia=1.d0
        endif

        dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
        r2_ji = Sum( dxp_ji(1:3)**2 )
        r_ji = dsqrt(r2_ji)
!unitl here ... next for other soap ...

        if (r_ji .gt. r_cut_desc) cycle
        ia_n = ia_n + 1
        k_desc_type(ia_n)=ia
!myBesselNoErase         do p=1,n_soap
!myBesselNoErase            x_b = 2.d0*alpha_soap*r_ji*rb_soap(p)
!myBesselNoErase            do l=0,l_max
!myBesselNoErase               if (x_b==0) then
!myBesselNoErase                  ms_bessel(0)=1d0
!myBesselNoErase                  ms_bessel(1:l_max)=0d0
!myBesselNoErase               else
!myBesselNoErase                  if (l==0) then
!myBesselNoErase                     ms_bessel(0)=sin(x_b)/x_b
!myBesselNoErase                     d_ms_bessel(0)=cos(x_b)/x_b - sin(x_b)/(x_b)**2
!myBesselNoErase                  elseif (l==1) then
!myBesselNoErase                     ms_bessel(1)=-cos(x_b)/x_b + sin(x_b)/(x_b)**2
!myBesselNoErase                     d_ms_bessel(1)= sin(x_b)/x_b + cos(x_b)/x_b**2  + cos(x_b)/(x_b)**2 - 0.5d0*sin(x_b)/(x_b)**3
!myBesselNoErase                  else
!myBesselNoErase                     ms_bessel(l)=-ms_bessel(l-2)+(2*l-1)*ms_bessel(l-1)/x_b
!myBesselNoErase                     d_ms_bessel(l)=ms_bessel(l-1)-(l+1)*ms_bessel(l)/x_b
!myBesselNoErase                  endif
!myBesselNoErase                  phi_ji(l,p) = ms_bessel(l)
!myBesselNoErase                  if (desc_forces) dphi_ji(l,p)=-2.d0*alpha_soap*r_ji*ms_bessel(l) + l*ms_bessel(l)/r_ji + 2.d0*d_ms_bessel(l)*rb_soap(p)
!myBesselNoErase               endif
!myBesselNoErase            enddo ! l = 1, l_max
!myBesselNoErase         end do ! n_soap

!Quip Bessel implementation .........................................
          do p = 1, n_soap
            x_b = 2.d0*alpha_soap * r_ji * rb_soap(p)
            exp_p = exp( -alpha_soap*( r_ji + rb_soap(p) )**2 )
            exp_m = exp( -alpha_soap*( r_ji - rb_soap(p) )**2 )
            do l = 0, l_max
              if( l == 0 ) then
                if( abs(x_b  - 0.d0).lt.1.d-14) then
                  !bessel_l = 1.0_dp
                  bessel_l = exp( -alpha_soap * (rb_soap(p)**2 + r_ji**2) )
                  if(desc_forces) bessel_lp = 0.0d0
                else
                  !bessel_lm = cosh(x_b)/x_b
                  !bessel_l = sinh(x_b)/x_b
                  bessel_lm = 0.5d0 * (exp_m + exp_p) / x_b
                  bessel_l  = 0.5d0 * (exp_m - exp_p) / x_b
                  if(desc_forces) bessel_lp = bessel_lm - (2*l+1)*bessel_l / x_b
                endif
              else
                if ( abs(x_b  - 0.d0).lt.1.d-14) then
                  bessel_l = 0.0d0
                  if(desc_forces) bessel_lp = 0.0d0
                else
                  bessel_lmm = bessel_lm
                  bessel_lm = bessel_l
                  if(desc_forces) then
                    bessel_l = bessel_lp
                    bessel_lp = bessel_lm - (2*l+1)*bessel_l / x_b
                  else
                    bessel_l = bessel_lmm - (2*l-1)*bessel_lm / x_b
                  endif
                endif
              endif
!Quip Bessel implementation .............................................
              !radial_fun(l,a) = exp( -alpha_soap * (rb_soap(p)**2 + r_ji**2) ) * bessel_l !* rb_soap(p)
              phi_ji(l,p) = bessel_l !* rb_soap(p)
              if(desc_forces) dphi_ji(l,p) = -2.0d0 * alpha_soap * r_ji * bessel_l + &
                     l*bessel_l / r_ji + bessel_lp * 2.0d0 * alpha_soap * rb_soap(p)

            enddo !l
          enddo   !p


         !debug write(6,*) 'r_ji=', r_ji
         !call r8mat_print(l_max+1,n_soap,phi_ji(0:l_max,1:n_soap), 'Bessel')

         !fcut and derivatives ... the cos_directions are not included ...
         if (r_ji .gt. (r_cut_desc-r_cut_width_soap)) then
           x_temp = one_pi*(r_ji - r_cut_desc + r_cut_width_soap)/r_cut_width_soap
           fcut=0.5d0*(cos(x_temp)+1.d0)*factor_ia
           if (desc_forces) dfcut= -0.5d0*one_pi/r_cut_width_soap*sin(x_temp)*factor_ia
         else
           fcut=1.d0
           dfcut=0.d0
         end if
         !radial functions multiplied with the good matrix ....
         rbf_ji(:,:) =matmul(phi_ji(:,:), W_soap(:,:))*fcut
         if (desc_forces) drbf_ji(:,:) =  matmul(phi_ji(:,:), W_soap(:,:))*dfcut + matmul(dphi_ji(:,:), W_soap(:,:))*fcut
         ! call r8mat_print(l_max+1,n_soap,rbf_ji(0:l_max,1:n_soap), 'Bessel+transform')


         do p=1,n_soap
            do l=0,l_max
               do m=-l,l
                  c(m,l,p) = c(m,l,p) + rbf_ji(l,p) * spherical_harm(l,m,dxp_ji(:))
                  if (desc_forces) dc(m,l,p,ia_n,1:3) =  dc(m,l,p,ia_n,1:3) +  &
                                                         drbf_ji(l,p)*(dxp_ji(1:3)/r_ji) * spherical_harm(l,m,dxp_ji(:)) + &
                                                         rbf_ji(l,p) * grad_spherical_harm(l,m,dxp_ji(:))
                  !debug write(*,*) spherical_harm(l,m,dxp_ji(:))
                  !debug stop
               enddo !m
            enddo    !l
         enddo       !p
   end do  !iw
!  just the descriptor without forces ....
   i_desc = 0
   do i_tmp = 1, nspecies_soap*n_soap
      p1 = ns_soap_index(1,i_tmp)
      i_s1 = ns_soap_index(2,i_tmp)
      do j_tmp = 1, i_tmp
         p2 =  ns_soap_index(1,j_tmp)
         i_s2 = ns_soap_index(2,j_tmp)

         if(lsoap_diag .and. p1 /= p2) cycle

         do l = 0, l_max
            i_desc = i_desc + 1
            !do m=-l,l
                 tmp_soap(i_desc)=real(dot_product(dconjg(c(:,l,p1)),c(:,l,p2)), kind(0.d0))
                 !tmp_soap(i_desc)=real(dot_product(c(:,l,p1),c(:,l,p2)))
                 !tmp_soap(i_desc)=tmp_soap(i_desc) + real(c(m,l,p1))*real(c(m,l,p2)) + conjg(c(m,l,p1))*conjg(c(m,l,p2))
            !end do
            !write(73,'(i5,f30.15)') i_desc, tmp_soap(i_desc)
            if(lsoap_lnorm) tmp_soap(i_desc) = tmp_soap(i_desc) / sqrt(2.d0 * dble(l) + 1.d0)
            if( i_tmp /= j_tmp )    tmp_soap(i_desc) = tmp_soap(i_desc) * sqrt_two
         enddo !l
      enddo !j_tmp
   enddo !i_tmp
   !version+1 quip
   tmp_soap(soap_dim) = 0.0d0
   norm_descriptor = sqrt(dot_product(tmp_soap(:),tmp_soap(:)))

   if(lsoap_norm) then
      local_soap(1:soap_dim)= tmp_soap(1:soap_dim) / norm_descriptor
   else
      local_soap(1:soap_dim)= tmp_soap(1:soap_dim)
   endif
   !version+1 quip
   local_soap(soap_dim)= 0.5d0

! the derivatives of the descrptors ...
   if (desc_forces) then
      do ia=1,ia_n
      i_desc = 0
         do i_tmp = 1, nspecies_soap*n_soap
            p1 = ns_soap_index(1,i_tmp)
            i_s1 = ns_soap_index(2,i_tmp)
            do j_tmp = 1, i_tmp
               p2 =  ns_soap_index(1,j_tmp)
               i_s2 = ns_soap_index(2,j_tmp)

               if(lsoap_diag .and. p1 /= p2) cycle

               do l = 0, l_max
                  i_desc = i_desc + 1
                  !do m=-l,l
                  do ix=1,3
                     !d_tmp_soap(i_desc,ix)=real( matmul(conjg(dc(:,l,p1,ia,ix)),c(:,l,p2)) + matmul( dc(:,l,p2,ia, ix), conjg(c(:,l,p1))) )
                     d_tmp_soap(i_desc,ix)=real( dot_product( dconjg(dc(:,l,p1,ia,ix)) , c(:,l,p2)  ) + &
                                                 dot_product( dc(:,l,p2,ia, ix) , dconjg(c(:,l,p1)) ), kind(0.d0) )
                  end do
                  !end do
                  if(lsoap_lnorm) d_tmp_soap(i_desc,1:3) = d_tmp_soap(i_desc,1:3) / sqrt(2.d0 * dble(l) + 1.d0)
                  if( i_tmp /= j_tmp )    d_tmp_soap(i_desc,1:3) = d_tmp_soap(i_desc,1:3) * sqrt_two
               enddo !l
            enddo !j_tmp
         enddo !i_tmp


         !version+1 quip
         d_tmp_soap(soap_dim,1:3) = 0.0d0
         if (lsoap_norm) then
           d_tmp_soap(:,:) = d_tmp_soap(:,:)/norm_descriptor
           do ix=1,3
               local_soap_deriv(:,ia,ix)= d_tmp_soap(:,ix) - tmp_soap(:) * &
                                          DOT_PRODUCT( tmp_soap(:) , d_tmp_soap(:,ix) ) / norm_descriptor**2
           end do
         else
            local_soap_deriv(:,ia,:)=d_tmp_soap(:,:)
         end if
         local_soap_deriv(:, 0,:) =  local_soap_deriv(:,0,:) - local_soap_deriv(:, ia,:)

      end do !ia from ia_n
   end if  !from desc ...
    !NON NORMALIZED VERSION
!mi    local_soap_out(:,ja) = local_soap_out(:,ja)*factor_ja
!mi    if (desc_forces) local_soap_deriv_out(:,ja, :,:) =  local_soap_deriv_out(:,ja, :,:)*factor_ja
   local_soap(:) = local_soap(:)*factor_ja
   if (desc_forces) local_soap_deriv(:, :,:) =  local_soap_deriv(:,:,:)*factor_ja

   n_desc_type = ia_n
enddo !ja
!debug write(6,*) 'Bye soap'
return
end subroutine compute_soap


subroutine  gen_dimension_for_soap
use ml_in_lammps_module, ONLY: l_max, n_soap, W_soap, lsoap_diag, soap_dim, nspecies_soap, &
                               rb_soap, ns_soap_index, one_pi, alpha_soap, &
                               rangml,  debug, atom_sigma_soap, S_factor_matrix
use milady_data, only : r_cut_desc
implicit none
integer :: i1, i2, i, is
integer :: i_bi, info
double precision :: r_cut_basis
double precision :: S_soap(n_soap, n_soap), C_soap(n_soap,n_soap)
double precision :: S_scale_factor(n_soap)
double precision :: C_scale_factor(n_soap), C_factor_matrix(n_soap,n_soap)
logical :: lfactorise
if (allocated(rb_soap)) deallocate(rb_soap) ; allocate(rb_soap(n_soap))
if (allocated(W_soap))  deallocate(W_soap)  ; allocate(W_soap(n_soap,n_soap))
if (allocated(S_factor_matrix))  deallocate(S_factor_matrix)  ; allocate(S_factor_matrix(n_soap,n_soap))

atom_sigma_soap = sqrt(0.5d0/alpha_soap)
!if (debug)
if (debug) then
  if (rangml==0) write(6, '("ML: init in gen_dimension_for_soap :",  2i5, 2f20.10)') &
                          n_soap, l_max, alpha_soap, atom_sigma_soap
end if
!debug if (rangml==0) write(6,'("ML: init in gen_dimension_for_soap :" 2i5, 2f20.10)') n_soap, l_max, alpha_soap, atom_sigma_soap


if (lsoap_diag) then
    ! +1 is the QUIP version ...
    soap_dim=(l_max+1)*n_soap*nspecies_soap*(nspecies_soap+1)/2+1
else
    ! +1 is the QUIP version ...
    soap_dim=(l_max+1)*n_soap*nspecies_soap*(n_soap*nspecies_soap+1)/2+1
end if

if (allocated(ns_soap_index)) deallocate(ns_soap_index); allocate(ns_soap_index(2,n_soap*nspecies_soap))
i_bi=0
do is=1,nspecies_soap
   do i=1,n_soap
      i_bi = i_bi +1
      ns_soap_index(1,i_bi)=i
      ns_soap_index(2,i_bi)=is
   end do
end do

r_cut_basis=r_cut_desc + atom_sigma_soap * sqrt(2.d0 * 10.d0 * log(10.d0))

do i=1,n_soap
  rb_soap(i)=dble(i-1)*r_cut_basis/dble(n_soap)
end do
!debug write (6,'("debug rb", 4f20.10)') (rb_soap(i),i=1,n_soap)


do i1 = 1, n_soap
   do i2 = 1, n_soap
      C_soap(i2,i1) = exp(-alpha_soap * (rb_soap(i1) - rb_soap(i2))**2)

      S_soap(i2,i1) = ( exp( -alpha_soap*(rb_soap(i1)**2+rb_soap(i2)**2) ) * &
                        sqrt(2.0d0) * alpha_soap**1.5d0 * (rb_soap(i1) + rb_soap(i2)) + &
                       alpha_soap*exp(-0.5d0 * alpha_soap * (rb_soap(i1) - &
                       rb_soap(i2))**2)*sqrt(one_pi)*(1.0d0 + alpha_soap*(rb_soap(i1) + rb_soap(i2))**2 ) * &
                       ( 1.0d0 + erf( sqrt(alpha_soap/2.0d0) * (rb_soap(i1) + rb_soap(i2)) ) ) )
   enddo
enddo
S_soap(:,:) = S_soap(:,:) / sqrt(128.d0*alpha_soap**5)
!debug write(6,*) 'C_soap'
!debug call r8mat_print (n_soap, n_soap, C_soap)


!debug write(6,*) 'S_soap'
!debug call r8mat_print (n_soap, n_soap, S_soap)

!factorize S_soap -> Cholesky_S_soap...

!first step of factorization of S_soap ....
S_scale_factor = 1.0d0
do i1 = 1, n_soap
   S_scale_factor(i1) = 1.0d0 / sqrt(S_soap(i1,i1))
enddo

lfactorise = ( maxval(S_scale_factor(:)) / minval(S_scale_factor(:)) < 0.1d0 )

if( lfactorise ) then
   do i1 = 1, n_soap
      S_factor_matrix(:,i1) = S_soap(:,i1)*S_scale_factor(:)*S_scale_factor(i1)
   enddo
else
   S_factor_matrix(:,:) = S_soap(:,:)
endif

!second step of factorization ... the factorization itself
!fill the loawer part ....
call dpotrf('L', n_soap, S_factor_matrix, n_soap, info)

if( info /= 0 ) then
   if (rangml==0) write(6, &
     '("ML: ERROR the matrix S_soap cannot be factorized in gen_dimension_for_soap, info for dptrof is ", i6)') info
   stop '<stop> in gen_dimension_for_soap for S_soap factorization'
endif

!put the matrix in order ....also fill the upper band.
do i1 = 2, n_soap
   do i2 = 1, i1
      S_factor_matrix(i2,i1) = S_factor_matrix(i1,i2)
   enddo
enddo


if( lfactorise ) then
   S_factor_matrix = 0.0d0
   do i1 = 1, n_soap
      do i2 = 1, n_soap
         S_factor_matrix(i2,i1) = S_factor_matrix(i2,i1) / S_scale_factor(i1)
      enddo
   enddo
end if
!end of the factorization ....

 !debug call r8mat_print (n_soap, n_soap, S_factor_matrix, 'Cholesky S out')


!factorize C_soap -> Cholesky_X_soap...

!first step of factorization of C_soap ....
C_scale_factor = 1.0d0
do i1 = 1, n_soap
   C_scale_factor(i1) = 1.0d0 / sqrt(C_soap(i1,i1))
enddo

lfactorise = ( maxval(C_scale_factor(:)) / minval(C_scale_factor(:)) < 0.1d0 )

if( lfactorise ) then
   do i1 = 1, n_soap
      C_factor_matrix(:,i1) = C_soap(:,i1)*C_scale_factor(:)*C_scale_factor(i1)
   enddo
else
   C_factor_matrix(:,:) = C_soap(:,:)
endif

!second step of factorization ... the factorization itself
!fill the loawer part ....
call dpotrf('L', n_soap, C_factor_matrix, n_soap, info)

if( info /= 0 ) then
   if (rangml==0) write(6, &
     '("ML: ERROR the matrix C_soap cannot be factorized in gen_dimension_for_soap, info for dptrof is ", i6)') info
   stop '<stop> in gen_dimension_for_soap C_soap for factorization'
endif

!put the matrix in order ....also fill the upper band.
do i1 = 2, n_soap
   do i2 = 1, i1
      C_factor_matrix(i2,i1) = C_factor_matrix(i1,i2)
   enddo
enddo


if( lfactorise ) then
   C_factor_matrix = 0.0d0
   do i1 = 1, n_soap
      do i2 = 1, n_soap
         C_factor_matrix(i2,i1) = C_factor_matrix(i2,i1) / C_scale_factor(i1)
      enddo
   enddo
end if
!end of the factorization ....




!solving lieanr system of equations .....
do i1 = 1, n_soap
   do i2 = 1, i1-1 !i + 1, this%n_max
      S_factor_matrix(i2,i1) = 0.0d0
   enddo
enddo



if( lfactorise ) then
   do i1 = 1, n_soap
      W_soap(:,i1) = S_factor_matrix(:,i1)*S_scale_factor(i1)
   enddo
else
   W_soap(:,:)=S_factor_matrix
endif

call dpotrs( 'U', n_soap, n_soap, C_factor_matrix, n_soap, W_soap, n_soap, info )
!this covariance
!b -> cholesky overlap, factor_matrix
!x -> trnasform

!debug write(6,*) 'gen_dimension_for_soap after dpotrs'
if( lfactorise ) then
   do i = 1, n_soap
      W_soap(:,i) = W_soap(:,i)*S_scale_factor(i)
   enddo
endif

!debug call r8mat_print (n_soap, n_soap, W_soap)
!debug stop

return
end subroutine gen_dimension_for_soap
