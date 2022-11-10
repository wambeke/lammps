

subroutine  compute_bispectrum_so4(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
!use mpi
use T_kind_param_m, ONLY:  double
use angular_functions, only : spherical_4d
use ml_in_lammps_module, ONLY: rangml,timetot, timetot00, timetot01, timetot02, timetot03, one_pi, jj_max, imm_neigh,    &
                               weighted, weighted_3ch,  bisso4_dim, lbso4_diag, desc_forces

use module_bispectrum_so4, only: czero, ZAcmm, ZBcmm, ZCcmm, ZAcmm2, ZBcmm2, ZCcmm2, ZAcmm3, ZBcmm3, ZCcmm3, &
                                 cmm, cmm2, cmm3, &
                                 class_mml, &
                                 fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch
use milady_data, only: local_bso4 => local_desc, local_bso4_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type , &
                       weight_per_type, weight_per_type_3ch
implicit none
integer, intent (in) :: i_start_at,i_final_at, nmax, ntype, numneigh_full
integer, dimension(numneigh_full), intent(in) :: firstneigh_full
!integer, dimension(numneigh_full), intent(in) :: fiirstneigh_full(numneigh_full)
integer, dimension(nmax), intent(in)  :: a_type
real(double), dimension(ntype), intent(in) :: map_mass
real(double), dimension(3,nmax), intent(in) :: xp

!----------------------------------------------------!
integer :: n_axpy
integer, parameter :: incx=1, incy=1
!----------------------------------------------------!

double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max) :: Umm
double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,1:3):: dUmm

real(double), dimension(3) :: dxp_ji
integer :: j,iw

integer :: l,l1,l2,m1,m2
integer :: ia, ia_n, ja,  i_bi, i_bi2, i_bi3
double precision :: r2_ji,r_ji, etmp, etmp2, etmp3
double complex :: sum_bi, sum_bi_w, sum_bi_w_3ch

double complex, allocatable, dimension(:) :: dsum_bi_ia, dsum_bi_ia_w, dsum_bi_ia_w_3ch
real(double)  :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch, &
                    tmpll

double complex, dimension(:), pointer ::  read_vecall
real(double), dimension(:,:), allocatable :: tmp_sum
integer :: l1_min, l1_max,  ia_type, ja_type
logical :: debug_time
real(double) :: time01, time01i, time02, time03, time04,  time03i

if (allocated(dsum_bi_ia)) deallocate(dsum_bi_ia) ; allocate(dsum_bi_ia(3*imm_neigh))
if (weighted) then
    if (allocated(dsum_bi_ia_w)) deallocate(dsum_bi_ia_w) ; allocate(dsum_bi_ia_w(3*imm_neigh))
  if (weighted_3ch) then
    if (allocated(dsum_bi_ia_w_3ch)) deallocate(dsum_bi_ia_w_3ch) ; allocate(dsum_bi_ia_w_3ch(3*imm_neigh))
  end if
end if
if (allocated(k_desc_type)) deallocate(k_desc_type) ; allocate(k_desc_type(numneigh_full))
!...i_bi...
local_bso4(:)=0.d0
local_bso4_deriv(:,:,:)=0.d0
k_desc_type(:)=0
if ( (i_start_at==0) .and. (i_final_at==0) ) return

debug_time=.false.
if (debug_time) call cpu_time(time01) !=MPI_WTIME()

do ja=i_start_at,i_final_at

    !if (weighted) then
    !     factor_ja = map_mass(a_type(ja)) / factor_weight_mass
    !else
    !    factor_ja=1.d0
    !endif
   ja_type = a_type(ja)
   if (weighted) then
      factor_ja= weight_per_type(ja_type)          !CHANGE config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
      if (weighted_3ch) factor_ja_3ch=  weight_per_type_3ch(ja_type)    !CHANGE config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
   else
     factor_ja=1.d0
   endif

    if (debug_time) call cpu_time(time01i) !=MPI_WTIME()
!...INIT cmm, cmm2, cmm3 ...............!
    !cmm(:,:,:) = 0.d0
    do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
      cmm(m2,m1,j)= 0.d0
    end do ; end do ; end do
    if (weighted)  then
      !cmm2(:,:,:)=0.d0
      do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
        cmm2(m2,m1,j)= 0.d0
      end do ; end do ; end do

      if (weighted_3ch) then
        !cmm3(:,:,:)=0.d0
        do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
          cmm3(m2,m1,j)= 0.d0
        end do ; end do ; end do
      end if
    end if

!...END INIT cmm, cmm2, cmm3 ...............!



   ia_n=0

    if (debug_time) then
         call cpu_time(time03i) !=MPI_WTIME()
      timetot00 = timetot00 + time03i - time01i
    end if
    !call neighbours(numneigh_full, firstneigh_full)
    !write(*,*) 'inn ', ja, ja_type, numneigh_full, r_cut_desc
    do iw=1, numneigh_full
        ia = firstneigh_full(iw)

        dxp_ji(1:3) = xp(1:3,ia) - xp(1:3,ja)
        r2_ji = dxp_ji(1)**2 + dxp_ji(2)**2 + dxp_ji(3)**2
        r_ji = dsqrt(r2_ji)

        if (r_ji .gt. r_cut_desc) cycle
        ia_n = ia_n + 1
        ia_type = a_type(ia)
        if (weighted) then
          factor_ia=   weight_per_type(ia_type) !CHANGE config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
          if (weighted_3ch) factor_ia_3ch=   weight_per_type_3ch(ia_type)   !CHANGE config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
        else
          factor_ia=1.d0
        endif

        !write(*,*) 'cc ', ia_type, weight_per_type(ia_type), weight_per_type_3ch(ia_type),  weighted, weighted_3ch
        fcut=0.5d0*(cos(one_pi*r_ji/r_cut_desc)+1.d0)
        dfcut= -0.5d0*one_pi/r_cut_desc*sin(one_pi*r_ji/r_cut_desc)
        if (weighted) then
          fcut_w=fcut*factor_ia
          dfcut_w=dfcut*factor_ia
          if (weighted_3ch) then
            fcut_w_3ch=fcut*factor_ia_3ch
            dfcut_w_3ch=dfcut*factor_ia_3ch
          end if
        end if

        call spherical_4d(ia_n, dxp_ji, r_ji, r_cut_desc, Umm, dUmm)

        do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
          cmm(m2,m1,j)= cmm(m2,m1,j)+Umm(m2,m1,j)*fcut
        end do ; end do ; end do
        if (weighted) then
          do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
            cmm2(m2,m1,j)= cmm2(m2,m1,j)+Umm(m2,m1,j)*fcut_w
          end do ; end do ; end do
          if (weighted_3ch) then
            do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
              cmm3(m2,m1,j)= cmm3(m2,m1,j)+Umm(m2,m1,j)*fcut_w_3ch
            end do ; end do ; end do
          end if
        end if

        k_desc_type(ia_n)=ia
    end do  !iw


    do j = 0,jj_max
      do m1=-j,j,2
        cmm(m1,m1,j) = cmm(m1,m1,j) +  1.d0
        if (weighted)  then
          cmm2(m1,m1,j) = cmm2(m1,m1,j) +  1.d0 ! factor_ja
          if (weighted_3ch) cmm3(m1,m1,j) =cmm3(m1,m1,j) + 1.d0 ! factor_ja_3ch
        end if
      end do
    end do


    if (debug_time) then
     call cpu_time(time03) !=MPI_WTIME()
     timetot01 = timetot01 + time03 - time03i
    end if

    n_desc_type = ia_n
    call pack_cmm_into_ZcmmABC(jj_max)

    i_bi=0
    n_axpy = 3*ia_n

    do l1=0,jj_max

      if (lbso4_diag) then
        !GaborLike l2=l1
        l1_min=l1
        l1_max=l1
      else
        !Tlike do l2=0,l1
        l1_min=0
        l1_max=l1
      end if

      do l2=l1_min, l1_max
        do l=abs(l1-l2),min(jj_max,l1+l2)

            if (mod(l1+l2+l,2)==1) cycle
            if (.not.(lbso4_diag)) then
              if (l < l1) cycle  ! this comes from Aidan Thompson and SNAP
            end if
            i_bi=i_bi+1
            etmp = 0.d0
            dsum_bi_ia(:)=czero
            if (weighted) then
              etmp2 = 0.d0
              dsum_bi_ia_w(:)=czero
              if (weighted_3ch) then
                etmp3 = 0.d0
                dsum_bi_ia_w_3ch(:)=czero
              end if

            end if

            do m1=-l,l,2
              do m2=-l,l,2
                  !energy
                  sum_bi = ZAcmm(m1,m2,i_bi)
                  etmp  = etmp + real ( conjg(cmm(m1,m2,l)) * sum_bi, kind=kind(0.d0))
                  if (weighted) then
                    sum_bi_w = ZAcmm2(m1,m2,i_bi)
                    etmp2  = etmp2 + real ( conjg(cmm2(m1,m2,l)) * sum_bi_w, kind=kind(0.d0))
                    if (weighted_3ch) then
                      sum_bi_w_3ch = ZAcmm3(m1,m2,i_bi)
                      etmp3  = etmp3 + real ( conjg(cmm3(m1,m2,l)) * sum_bi_w_3ch, kind=kind(0.d0))
                    end if
                  end if
                  if (desc_forces) then
                    read_vecall => class_mml(m1,m2,l)%vecall
                    call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                    nullify(read_vecall)
                    !do ia=1,ia_n
                    !do ix=1,3
                    !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                    !end do
                    !end do
                    !forces
                    if (weighted) then
                      read_vecall => class_mml(m1,m2,l)%vecall_w
                      call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                      nullify(read_vecall)
                      !do ix=1,3
                      !do ia=1,ia_n
                      !   dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l,ia,ix))*sum_bi_w
                      !end do
                      !end do
                      if (weighted_3ch) then
                        read_vecall => class_mml(m1,m2,l)%vecall_w_3ch
                        call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                        nullify(read_vecall)
                        !do ix=1,3
                        !do ia=1,ia_n
                        !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l,ia,ix))*sum_bi_w_3ch
                        !end do
                        !end do
                      end if
                    end if

                  end if
              end do  !m2
            end do    !m1

            tmpll = dble(1+l)/dble(1+l1)
            do m1=-l1,l1,2
              do m2=-l1,l1,2
                  sum_bi = ZBcmm(m1,m2,i_bi)*tmpll
                  if (weighted) then
                    sum_bi_w = ZBcmm2(m1,m2,i_bi)*tmpll
                    if (weighted_3ch) sum_bi_w_3ch = ZBcmm3(m1,m2,i_bi)*tmpll
                  end if

                  if (desc_forces) then
                    read_vecall => class_mml(m1,m2,l1)%vecall
                    call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                    nullify(read_vecall)
                    !do ia=1,ia_n
                    !do ix=1,3
                        !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                    !end do
                    !end do
                    if (weighted) then
                      read_vecall => class_mml(m1,m2,l1)%vecall_w
                      call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                      nullify(read_vecall)
                      !do ix=1,3
                      !do ia=1,ia_n
                      !  !TOC dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l1,ia,ix))*sum_bi_w
                      !end do
                      !end do

                      if (weighted_3ch) then
                        read_vecall => class_mml(m1,m2,l1)%vecall_w_3ch
                        call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                        nullify(read_vecall)
                        !do ix=1,3
                        !do ia=1,ia_n
                        !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l1,ia,ix))*sum_bi_w_3ch
                        !end do
                        !end do
                      end if
                    end if
                  end if

              end do  !m2
            end do    !m1

            tmpll = dble(1+l)/dble(1+l2)
            do m1=-l2,l2,2
              do m2=-l2,l2,2
                  !sum_bi = ZCcmm(m1,m2,i_bi)*dble(1+l)/dble(1+l2)
                  sum_bi = ZCcmm(m1,m2,i_bi)*tmpll
                  if (weighted) then
                    sum_bi_w = ZCcmm2(m1,m2,i_bi)*tmpll
                    if (weighted_3ch) sum_bi_w_3ch = ZCcmm3(m1,m2,i_bi)*tmpll
                  end if

                  if (desc_forces) then
                    read_vecall => class_mml(m1,m2,l2)%vecall
                    call zaxpy(n_axpy,-sum_bi, read_vecall, incx, dsum_bi_ia,incy)
                    nullify(read_vecall)
                    !do ia=1,ia_n
                    !do ix=1,3
                        !dsum_bi_ia_all(1:3*ia_n) = dsum_bi_ia_all(1:3*ia_n)   - read_vecall(1:3*ia_n)*sum_bi
                    !end do
                    !end do
                    if (weighted) then
                      read_vecall => class_mml(m1,m2,l2)%vecall_w
                      call zaxpy(n_axpy,-sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                      nullify(read_vecall)
                      !do ix=1,3
                      !do ia=1,ia_n
                      !  !TOC dsum_bi_ia_w(ia,ix) = dsum_bi_ia_w(ia,ix)   - conjg(dcmm2(m1,m2,l2,ia,ix))*sum_bi_w
                      !end do
                      !end do
                      if (weighted_3ch) then
                        read_vecall => class_mml(m1,m2,l2)%vecall_w_3ch
                        call zaxpy(n_axpy,-sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                        nullify(read_vecall)
                        !do ix=1,3
                        !do ia=1,ia_n
                        !  !TOC dsum_bi_ia_w_3ch(ia,ix) = dsum_bi_ia_w_3ch(ia,ix)   - conjg(dcmm3(m1,m2,l2,ia,ix))*sum_bi_w_3ch
                        !end do
                        !end do
                      end if
                    end if
                  end if

              end do  !m2
            end do    !m1


            local_bso4(i_bi) = etmp
            if (weighted) then
              i_bi2 = i_bi+bisso4_dim
              i_bi3 = i_bi2+bisso4_dim
              local_bso4(i_bi2) = etmp2
              if (weighted_3ch) local_bso4(i_bi3) = etmp3
            end if

            if (desc_forces) then
              do ia=1,ia_n
                local_bso4_deriv(i_bi, ia,1) = real ( dsum_bi_ia(3*ia-2) , kind=kind(0.d0))
                local_bso4_deriv(i_bi, ia,2) = real ( dsum_bi_ia(3*ia-1) , kind=kind(0.d0))
                local_bso4_deriv(i_bi, ia,3) = real ( dsum_bi_ia(3*ia)   , kind=kind(0.d0))
              end do
              if (weighted) then
                do ia=1,ia_n
                  local_bso4_deriv(i_bi2, ia,1) = real ( dsum_bi_ia_w(3*ia-2) , kind=kind(0.d0))
                  local_bso4_deriv(i_bi2, ia,2) = real ( dsum_bi_ia_w(3*ia-1) , kind=kind(0.d0))
                  local_bso4_deriv(i_bi2, ia,3) = real ( dsum_bi_ia_w(3*ia)   , kind=kind(0.d0))
                end do
                if (weighted_3ch) then
                  do ia=1,ia_n
                    local_bso4_deriv(i_bi3, ia,1) = real ( dsum_bi_ia_w_3ch(3*ia-2) , kind=kind(0.d0))
                    local_bso4_deriv(i_bi3, ia,2) = real ( dsum_bi_ia_w_3ch(3*ia-1) , kind=kind(0.d0))
                    local_bso4_deriv(i_bi3, ia,3) = real ( dsum_bi_ia_w_3ch(3*ia)   , kind=kind(0.d0))
                  end do
                end if
              end if
            end if

        end do !l
      end do !l2
    end do !l1
    if (allocated(tmp_sum)) deallocate(tmp_sum) ; allocate(tmp_sum(size(local_bso4,1),3))
    !do ia=1,ia_n
       !if (desc_forces) local_bso4_deriv(:, 0,:) =  local_bso4_deriv(:, 0,:) -  local_bso4_deriv(:, ia,:)
       !if (desc_forces) local_bso4_deriv(:, 0,:) =   -  sum(local_bso4_deriv(:, :,:), dim=2)
       tmp_sum(:,:) = - sum(local_bso4_deriv(:, 1:ia_n,:), dim=2)
    !enddo  !ia from ia_m
    local_bso4_deriv(:, 0,:) =  tmp_sum(:,:)
    if (debug_time) then
     call cpu_time(time04) !=MPI_WTIME()
     timetot02 = timetot02 + time04 - time03
    end if


    if (desc_forces) local_bso4_deriv(:, 0:ia_n,1:3) =   - local_bso4_deriv(:,0:ia_n,1:3)
    if (weighted) then
      local_bso4(bisso4_dim+1:2*bisso4_dim) = local_bso4(bisso4_dim+1:2*bisso4_dim) !*factor_ja
      if (weighted_3ch) local_bso4(2*bisso4_dim+1:3*bisso4_dim) =    local_bso4(2*bisso4_dim+1:3*bisso4_dim) !*factor_ja_3ch
      !if (desc_forces) then
        !local_bso4_deriv(bisso4_dim+1:2*bisso4_dim, :,:) =  - local_bso4_deriv(bisso4_dim+1:2*bisso4_dim,:,:)!*factor_ja
        !if (weighted_3ch) local_bso4_deriv(2*bisso4_dim+1:3*bisso4_dim, :,:) =  - local_bso4_deriv(2*bisso4_dim+1:3*bisso4_dim, :,:)!*factor_ja_3ch
      !end if
    end if
  enddo !ja

if (debug_time) then
   call cpu_time(time02) !=MPI_WTIME()
   timetot = timetot + time02 - time01
   timetot03 = timetot03 + time02 - time04
  write(*,*) 'timetot  ', timetot
  write(*,*) 'timetot00', timetot00
  write(*,*) 'timetot01', timetot01
  write(*,*) 'timetot02', timetot02
  write(*,*) 'timetot03', timetot03
end if


return
end subroutine compute_bispectrum_so4


subroutine gen_dimension_for_bispectrum_so4
use module_bispectrum_so4, only: bisso4_cg_dim, bisso4_cg_full_dim, &
                                 bisso4_cg_A_dim,bisso4_cg_B_dim, bisso4_cg_C_dim, &
                                 cg_A, cg_B, cg_C, &
                                 cmm, cmm2, cmm3, &
                                 ZAcmm, ZAcmm2, ZAcmm3, &
                                 ZBcmm, ZBcmm2, ZBcmm3, &
                                 ZCcmm, ZCcmm2, ZCcmm3, &
                                 class_mml

use ml_in_lammps_module, ONLY:  imm_neigh, weighted, weighted_3ch, &
                                lbso4_diag, jj_max, bisso4_l, bisso4_l1, bisso4_l2, bisso4_dim, cg_vector
implicit none
integer :: l1, l2, l,j, l1_min, l1_max, m1,m2,ma,mb,m1ma, m2mb
integer :: i_bi, it_cg, it_cgA, it_cgB, it_cgC, it_umm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            pre-compute dimension of BSO4                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


i_bi=0
it_cg=0
do l1=0,jj_max
  !l2=l1 ! following Gabor  only the diagonal elements are important
    if (lbso4_diag) then
      !GaborLike l2=l1
      l1_min=l1
      l1_max=l1
    else
      !Tlike do l2=0,l1
      l1_min=0
      l1_max=l1
    end if
    do l2=l1_min, l1_max
      !do l2=0, l1  ! in the end we want only the componenets with l1 <= l2 <= l
      do l=abs(l1-l2),min(jj_max,l1+l2)
      if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
      if (.not.(lbso4_diag)) then
        if (l < l1) cycle  ! this comes from Thompson
      end if
      i_bi = i_bi+1
      do m1=-l,l,2
        do m2=-l,l,2
          do ma=max(-l1,m1-l2),min(l1,m1+l2),2
            do mb=max(-l1,m2-l2),min(l1,m2+l2),2
              it_cg = it_cg+1
            end do
          end do
        end do
      end do
     enddo
  end do
enddo

bisso4_cg_dim=it_cg
bisso4_dim=i_bi
if (allocated(bisso4_l  )) deallocate(bisso4_l  ) ; allocate (bisso4_l(bisso4_dim))
if (allocated(bisso4_l1 )) deallocate(bisso4_l1 ) ; allocate (bisso4_l1(bisso4_dim))
if (allocated(bisso4_l2 )) deallocate(bisso4_l2 ) ; allocate (bisso4_l2(bisso4_dim))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       pre-compute cgA,cgB,cgC in order to compute           !
!                      i) the dimensions 2) computation                       !
!                      used for packing ZAcmm ZBcmm and ZCcmm                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



i_bi=0
it_cgA=0
it_cgB=0
it_cgC=0
it_umm=0
do l1=0,jj_max
  !l2=l1 ! following Gabor  only the diagonal elements are important
    if (lbso4_diag) then
      !GaborLike l2=l1
      l1_min=l1
      l1_max=l1
    else
      !Tlike do l2=0,l1
      l1_min=0
      l1_max=l1
    end if
    do l2=l1_min, l1_max
     !do l2=0, l1
        do l=abs(l1-l2),min(jj_max,l1+l2)
          if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
          if (.not.(lbso4_diag)) then
          if (l < l1) cycle           ! in the end we want only the componenets with l1 <= l2 <= l
          end if
          i_bi = i_bi+1
          bisso4_l (i_bi) = l
          bisso4_l1(i_bi) = l1
          bisso4_l2(i_bi) = l2

          do m1=-l,l,2
            do m2=-l,l,2
              it_umm = it_umm+1
              do ma=max(-l1,m1-l2),min(l1,m1+l2),2
                do mb=max(-l1,m2-l2),min(l1,m2+l2),2
                  it_cgA = it_cgA+1
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


          do m1=-l1,l1,2
            do m2=-l1,l1,2
              it_umm = it_umm+1
              do ma=max(-l,m1-l2),min(l ,m1+l2),2
                do mb=max(-l,m2-l2),min(l ,m2+l2),2
                  it_cgB = it_cgB+1
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


          do m1=-l2,l2,2
            do m2=-l2,l2,2
              it_umm = it_umm+1
              do ma=max(-l1,m1-l),min(l1 ,m1+l),2
                do mb=max(-l1,m2-l),min(l1 ,m2+l),2
                  it_cgC = it_cgC+1
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


        end do
     end do
enddo


bisso4_cg_A_dim = it_cgA
bisso4_cg_B_dim = it_cgB
bisso4_cg_C_dim = it_cgC

if (allocated(cg_A)) deallocate(cg_A) ; allocate(cg_A(bisso4_cg_A_dim))
if (allocated(cg_B)) deallocate(cg_B) ; allocate(cg_B(bisso4_cg_B_dim))
if (allocated(cg_C)) deallocate(cg_C) ; allocate(cg_C(bisso4_cg_C_dim))
cg_A(:)= 0.d0
cg_B(:)= 0.d0
cg_C(:)= 0.d0

i_bi=0
it_cgA=0
it_cgB=0
it_cgC=0
do l1=0,jj_max
  !l2=l1 ! following Gabor  only the diagonal elements are important
    if (lbso4_diag) then
      !GaborLike l2=l1
      l1_min=l1
      l1_max=l1
    else
      !Tlike do l2=0,l1
      l1_min=0
      l1_max=l1
    end if
    do l2=l1_min, l1_max
        do l=abs(l1-l2),min(jj_max,l1+l2)
          if (mod(l1+l2+l,2)==1) cycle ! assure invariance par réflexion
          if (.not.(lbso4_diag)) then
          if (l < l1) cycle           ! in the end we want only the componenets with l1 <= l2 <= l
          end if
          i_bi = i_bi+1

          do m1=-l,l,2
            do m2=-l,l,2
              do ma=max(-l1,m1-l2),min(l1,m1+l2),2
                do mb=max(-l1,m2-l2),min(l1,m2+l2),2
                  m1ma= m1-ma
                  m2mb =m2-mb
                  it_cgA = it_cgA+1
                  cg_A(it_cgA) =cg_vector(l1,ma,l2,m1ma,l,m1)*cg_vector(l1,mb,l2,m2mb,l,m2)
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


          do m1=-l1,l1,2
            do m2=-l1,l1,2
              do ma=max(-l,m1-l2),min(l,m1+l2),2
                do mb=max(-l,m2-l2),min(l,m2+l2),2
                  m1ma= m1-ma
                  m2mb =m2-mb
                  it_cgB = it_cgB+1
                  cg_B(it_cgB) =cg_vector(l,ma,l2,m1ma,l1,m1)*cg_vector(l,mb,l2,m2mb,l1,m2)
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


          do m1=-l2,l2,2
            do m2=-l2,l2,2
              do ma=max(-l1,m1-l),min(l1,m1+l),2
                do mb=max(-l1,m2-l),min(l1,m2+l),2
                  m1ma= m1-ma
                  m2mb =m2-mb
                  it_cgC = it_cgC+1
                  cg_C(it_cgC) =cg_vector(l1,ma,l,m1ma,l2,m1)*cg_vector(l1,mb,l,m2mb,l2,m2)
                end do !ma
              end do   !mb
            end do     !m1
          end do       !m2


        end do
    end do
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (allocated(class_mml)) deallocate(class_mml) ; allocate(class_mml(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max))
do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
  if (allocated(class_mml(m2,m1,j)%vecall)) deallocate(class_mml(m2,m1,j)%vecall)
                                            allocate(class_mml(m2,m1,j)%vecall(1:3*imm_neigh))
end do ; end do ; end do

if (allocated(cmm))   deallocate(cmm);      allocate(  cmm(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )
if (allocated(ZAcmm)) deallocate(ZAcmm);    allocate(ZAcmm(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
if (allocated(ZBcmm)) deallocate(ZBcmm);    allocate(ZBcmm(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
if (allocated(ZCcmm)) deallocate(ZCcmm);    allocate(ZCcmm(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )

if (weighted) then

  if (allocated(cmm2))   deallocate(cmm2);      allocate(  cmm2(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )
  do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
    if (allocated(class_mml(m2,m1,j)%vecall_w)) deallocate(class_mml(m2,m1,j)%vecall_w)
                                                allocate(class_mml(m2,m1,j)%vecall_w(1:3*imm_neigh))
  end do ; end do ; end do
  if (allocated(ZAcmm2)) deallocate(ZAcmm2);    allocate(ZAcmm2(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
  if (allocated(ZBcmm2)) deallocate(ZBcmm2);    allocate(ZBcmm2(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
  if (allocated(ZCcmm2)) deallocate(ZCcmm2);    allocate(ZCcmm2(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )

  if (weighted_3ch) then

    if (allocated(cmm3))   deallocate(cmm3);      allocate(  cmm3(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )
    do j= 0, jj_max ; do m1 = -j,j,2; do m2 = -j,j,2
      if (allocated(class_mml(m2,m1,j)%vecall_w_3ch)) deallocate(class_mml(m2,m1,j)%vecall_w_3ch)
                                                      allocate(class_mml(m2,m1,j)%vecall_w_3ch(1:3*imm_neigh))
    end do ; end do ; end do
    if (allocated(ZAcmm3)) deallocate(ZAcmm3);    allocate(ZAcmm3(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
    if (allocated(ZBcmm3)) deallocate(ZBcmm3);    allocate(ZBcmm3(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )
    if (allocated(ZCcmm3)) deallocate(ZCcmm3);    allocate(ZCcmm3(-jj_max:jj_max,-jj_max:jj_max, bisso4_dim) )

  end if

end if

return
end subroutine gen_dimension_for_bispectrum_so4


subroutine test_parameters_bso4
use ml_in_lammps_module, only: rangml, j_max, inv_r0_input, inv_r0
implicit none
    if (rangml==0) then
         if (j_max<0) then
             write(6,'("ML: bi_so4 parametrization j_max should be larger or at least equals 0")')
             stop 'read_ml_file bi_so4 j_max negative'
         end if
         if (.not.((inv_r0_input > 0.d0).and.(inv_r0_input < 1.d0)) ) then
            write(6,*) 'ML:', inv_r0, inv_r0_input
            write(6,'("ML: bi_so4 parametrization, inv_r0_input should be larger than 0 and lower than 1")')
             stop 'read_ml_file bi_so4 inv_r0_input beyound the limits'
         endif
    endif


end subroutine test_parameters_bso4


subroutine pack_cmm_into_ZcmmABC(jj_max)
  use module_kind_variables, ONLY:  double
  use ml_in_lammps_module, only: lbso4_diag, weighted, weighted_3ch
  use module_bispectrum_so4, only: czero, cmm, cmm2, cmm3,&
                                   cg_a, cg_B, cg_C, &
                                   ZAcmm, ZBcmm, ZCcmm, &
                                   ZAcmm2, ZBcmm2, ZCcmm2, &
                                   ZAcmm3, ZBcmm3, ZCcmm3


  implicit none
  integer, intent(in) :: jj_max

  integer :: i_bi,  l1,l2,l, l1_min, l1_max, m1ma, m2mb
  integer :: m1,m2,ma,mb
  integer :: itfA, itfB, itfC
  double complex :: sum_bi, sum_bi_w, sum_bi_w_3ch
  real(double) :: cg_local

    !ZAcmm(:,:,:) = czero
    !ZBcmm(:,:,:) = czero
    !ZCcmm(:,:,:) = czero

    !ZAcmm2(:,:,:) = czero
    !ZBcmm2(:,:,:) = czero
    !ZCcmm2(:,:,:) = czero

    !ZAcmm3(:,:,:) = czero
    !ZBcmm3(:,:,:) = czero
    !ZCcmm3(:,:,:) = czero

    i_bi=0
    itfA =0
    itfB =0
    itfC =0
    do l1=0,jj_max

      if (lbso4_diag) then
        !GaborLike l2=l1
        l1_min=l1
        l1_max=l1
      else
        !Tlike do l2=0,l1
        l1_min=0
        l1_max=l1
      end if

      do l2=l1_min, l1_max
        do l=abs(l1-l2),min(jj_max,l1+l2)

            if (mod(l1+l2+l,2)==1) cycle
            if (.not.(lbso4_diag)) then
              if (l < l1) cycle  ! this comes from Aidan Thompson and SNAP
            end if
            i_bi=i_bi+1

            do m1=-l,l,2
              do m2=-l,l,2
                sum_bi=czero
                sum_bi_w=czero
                sum_bi_w_3ch=czero
                  do ma=max(-l1,m1-l2),min(l1,m1+l2),2
                    do mb=max(-l1,m2-l2),min(l1,m2+l2),2
                      itfA = itfA + 1
                      m1ma= m1-ma
                      m2mb =m2-mb
                      cg_local = cg_A(itfA)
                      sum_bi = sum_bi +  cg_local*cmm(ma,mb,l1)*cmm(m1ma, m2mb, l2)
                      if (weighted) then
                        sum_bi_w = sum_bi_w + cg_local*cmm2(ma,mb,l1)*cmm2(m1ma, m2mb, l2)
                        if (weighted_3ch)   sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma,mb,l1)*cmm3(m1ma, m2mb, l2)
                      end if
                    end do !mb
                  end do ! ma
                  ZAcmm(m1,m2, i_bi) = sum_bi
                  if (weighted) then
                    ZAcmm2(m1,m2,i_bi) = sum_bi_w
                    if (weighted_3ch) ZAcmm3(m1,m2,i_bi) = sum_bi_w_3ch
                  end if
              enddo  !m2
            enddo   !m1

            do m1=-l1,l1,2
              do m2=-l1,l1,2
                  sum_bi=czero
                  sum_bi_w=czero
                  sum_bi_w_3ch=czero
                  do ma=max(-l,m1-l2),min(l,m1+l2),2
                    do mb=max(-l,m2-l2),min(l,m2+l2),2
                      itfB = itfB + 1
                      m1ma= m1-ma
                      m2mb =m2-mb
                      cg_local = cg_B(itfB)
                      sum_bi = sum_bi +  cg_local*cmm(ma,mb,l)*cmm(m1ma, m2mb, l2)
                      if (weighted) then
                        sum_bi_w = sum_bi_w + cg_local*cmm2(ma,mb,l)*cmm2(m1ma, m2mb, l2)
                       if(weighted_3ch)   sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma,mb,l)*cmm3(m1ma, m2mb, l2)
                      end if
                    end do !mb
                  end do ! ma
                  ZBcmm(m1,m2, i_bi) = sum_bi
                  if (weighted) then
                    ZBcmm2(m1,m2,i_bi) = sum_bi_w
                    if (weighted_3ch) ZBcmm3(m1,m2,i_bi) = sum_bi_w_3ch
                  end if
              enddo  !m2
            enddo   !m1

            do m1=-l2,l2,2
              do m2=-l2,l2,2
                  sum_bi=czero
                  sum_bi_w=czero
                  sum_bi_w_3ch=czero
                  do ma=max(-l1,m1-l),min(l1,m1+l),2
                    do mb=max(-l1,m2-l),min(l1,m2+l),2
                      itfC = itfC + 1
                      m1ma= m1-ma
                      m2mb =m2-mb
                      cg_local = cg_C(itfC)
                      sum_bi = sum_bi +  cg_local*cmm(ma,mb,l1)*cmm(m1ma, m2mb, l)
                      if (weighted) then
                        sum_bi_w = sum_bi_w + cg_local*cmm2(ma,mb,l1)*cmm2(m1ma, m2mb, l)
                        if(weighted_3ch) sum_bi_w_3ch = sum_bi_w_3ch + cg_local*cmm3(ma,mb,l1)*cmm3(m1ma, m2mb, l)
                      end if
                    end do !mb
                  end do ! ma
                  ZCcmm(m1,m2, i_bi) = sum_bi
                  if (weighted) then
                       ZCcmm2(m1,m2,i_bi) = sum_bi_w
                       if (weighted_3ch) ZCcmm3(m1,m2,i_bi) = sum_bi_w_3ch
                  end if
              enddo  !m2
            enddo    !m1


          end do  !l
        end do  !l2  Tlike enddo
      enddo !l1

return
end subroutine pack_cmm_into_ZcmmABC
