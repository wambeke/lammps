subroutine compute_pow_so4(i_start_at, i_final_at, nmax, xp, ntype, a_type, map_mass, numneigh_full, firstneigh_full)
!use mpi
use T_kind_param_m, ONLY:  double
use angular_functions, only : spherical_4d
use ml_in_lammps_module, ONLY: rangml,timetot, timetot00, timetot01, timetot02, timetot03, one_pi, jj_max, imm_neigh,    &
                               weighted, weighted_3ch,  pow_so4_dim, desc_forces

use module_bispectrum_so4, only: cmm, cmm2, cmm3, class_mml, &
                                 czero, fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch

use milady_data, only: local_pso4 => local_desc, local_pso4_deriv => local_desc_deriv, r_cut_desc, n_desc_type, k_desc_type , &
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
integer :: ia, ia_n, ja,  icnt, icnt2, icnt3
double precision :: r2_ji,r_ji, etmp, etmp2, etmp3
double complex :: sum_bi, sum_bi_w, sum_bi_w_3ch, ctmp, ctmp2, ctmp3

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
local_pso4(:)=0.d0
local_pso4_deriv(:,:,:)=0.d0
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
      factor_ja = weight_per_type(ja_type)          !CHANGE config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
      if (weighted_3ch) factor_ja_3ch = weight_per_type_3ch(ja_type)    !CHANGE config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
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

    do j = 0,jj_max
      do m1=-j,j,2
        cmm(m1,m1,j) = 1.d0
        if (weighted)  then
          cmm2(m1,m1,j) = factor_ja
          if (weighted_3ch) cmm3(m1,m1,j) = factor_ja_3ch
        end if
      end do
    end do
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
          factor_ia = weight_per_type(ia_type) !CHANGE config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
          if (weighted_3ch) factor_ia_3ch = weight_per_type_3ch(ia_type)   !CHANGE config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
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

    if (debug_time) then
     call cpu_time(time03) !=MPI_WTIME()
     timetot01 = timetot01 + time03 - time03i
    end if

    n_desc_type = ia_n

    icnt=0
    n_axpy = 3*ia_n

    do l=0,jj_max
        icnt = icnt + 1
        etmp=0.d0
        etmp2=0.d0
        etmp3=0.d0
        if (desc_forces) then
          dsum_bi_ia(1:3*ia_n)=czero
          if (weighted) then
            dsum_bi_ia_w(1:3*ia_n)=czero
            if (weighted_3ch) then
              dsum_bi_ia_w_3ch(1:3*ia_n)=czero
            end if
          end if
        end if

        do m1=-l,l,2
          do m2=-l,l,2
            ctmp = cmm(m2,m1,l)
            !config_desc(iconf)%energy(icnt,ja) = config_desc(iconf)%energy(icnt,ja) + real(conjg(ctmp) * ctmp)
            etmp = etmp + real(conjg(ctmp) * ctmp)
            if (weighted) then
              ctmp2 = cmm2(m2,m1,l)
              etmp2 = etmp2 + real(conjg(ctmp2) * ctmp2)
              if (weighted_3ch) then
                ctmp3 = cmm3(m2,m1,l)
                etmp3 = etmp3 + real(conjg(ctmp3) * ctmp3)
              end if
            end if
            if (desc_forces) then
              sum_bi = 2.d0*ctmp
              read_vecall => class_mml(m2,m1,l)%vecall
              !do ia=1,ia_n
              !  config_desc(iconf)%force(icnt,ja,ia,1:3) = config_desc(iconf)%force(icnt,ja,ia,1:3) + conjg(dcmm(m2,m1,l,ia,1:3))*sum_bi
              !end do !ia
              call zaxpy(n_axpy, sum_bi, read_vecall, incx, dsum_bi_ia,incy)
              nullify(read_vecall)
              if (weighted) then
                sum_bi_w = 2.d0*ctmp2
                read_vecall => class_mml(m2,m1,l)%vecall_w
                call zaxpy(n_axpy, sum_bi_w, read_vecall, incx, dsum_bi_ia_w,incy)
                nullify(read_vecall)
              end if
              if (weighted_3ch) then
                sum_bi_w_3ch = 2.d0*ctmp3
                read_vecall => class_mml(m2,m1,l)%vecall_w_3ch
                call zaxpy(n_axpy, sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch,incy)
                nullify(read_vecall)
              end if
            end if ! desc_forces
          enddo  !m2
        enddo    !m1


        !fix energy descriptor  ...
        local_pso4(icnt) = etmp
        if (weighted) then
          icnt2 = pow_so4_dim + icnt
          icnt3 = pow_so4_dim + icnt2
          local_pso4(icnt2) = etmp2
          if (weighted_3ch) then
            local_pso4(icnt3) = etmp3
          end if
        end if

        ! fix forces descriptor ...
        if (desc_forces) then
          do ia = 1, ia_n
            local_pso4_deriv(icnt,ia,1) = real(dsum_bi_ia(3*ia-2), kind=kind(1.d0))
            local_pso4_deriv(icnt,ia,2) = real(dsum_bi_ia(3*ia-1), kind=kind(1.d0))
            local_pso4_deriv(icnt,ia,3) = real(dsum_bi_ia(3*ia  ), kind=kind(1.d0))
          end do
          if (weighted) then
            do ia = 1, ia_n
              local_pso4_deriv(icnt2,ia,1) = real(dsum_bi_ia_w(3*ia-2), kind=kind(1.d0))
              local_pso4_deriv(icnt2,ia,2) = real(dsum_bi_ia_w(3*ia-1), kind=kind(1.d0))
              local_pso4_deriv(icnt2,ia,3) = real(dsum_bi_ia_w(3*ia  ), kind=kind(1.d0))
            end do
            if (weighted_3ch) then
              do ia = 1, ia_n
                local_pso4_deriv(icnt3,ia,1) = real(dsum_bi_ia_w_3ch(3*ia-2), kind=kind(1.d0))
                local_pso4_deriv(icnt3,ia,2) = real(dsum_bi_ia_w_3ch(3*ia-1), kind=kind(1.d0))
                local_pso4_deriv(icnt3,ia,3) = real(dsum_bi_ia_w_3ch(3*ia  ), kind=kind(1.d0))
              end do
            end if
          end if
        end if

    enddo ! l


    if (allocated(tmp_sum)) deallocate(tmp_sum) ; allocate(tmp_sum(size(local_pso4,1),3))
    !do ia=1,ia_n
       !if (desc_forces) local_pso4_deriv(:, 0,:) =  local_pso4_deriv(:, 0,:) -  local_pso4_deriv(:, ia,:)
       tmp_sum(:,:) = - sum(local_pso4_deriv(:, 1:ia_n,:), dim=2)
    !enddo  !ia from ia_m
    local_pso4_deriv(:, 0,:) =  tmp_sum(:,:)


    !config_desc(iconf)%energy(1:pow_so4_dim,ja) = config_desc(iconf)%energy(1:pow_so4_dim,ja)
    !if (weighted) then
    !  local_pso4(pow_so4_dim+1:2*pow_so4_dim) = local_pso4(pow_so4_dim+1:2*pow_so4_dim) !*factor_ja
    !  if (weighted_3ch) local_pso4(2*pow_so4_dim+1:3*pow_so4_dim) = local_pso4(2*pow_so4_dim+1:3*pow_so4_dim) !*factor_ja_3ch
    !  if (desc_forces_local) then
    !    local_pso4_deriv(pow_so4_dim+1:2*pow_so4_dim, :,:) =  local_pso4_deriv(pow_so4_dim+1:2*pow_so4_dim, :,:) !*factor_ja
    !    if (weighted_3ch) local_pso4_deriv(2*pow_so4_dim+1:3*pow_so4_dim, :,:) =  local_pso4_deriv(2*pow_so4_dim+1:3*pow_so4_dim, :,:) !*factor_ja_3ch
    !  end if
    !end if
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
end subroutine compute_pow_so4

subroutine init_pow_so4
use ml_in_lammps_module, only: jj_max, imm_neigh, weighted, weighted_3ch
use module_bispectrum_so4, only: class_mml, cmm, cmm2, cmm3
implicit none
integer :: j,m1,m2

if (allocated(class_mml)) deallocate(class_mml) ; allocate(class_mml(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max))
do j= 0, jj_max
  do m1 = -j,j,2
    do m2 = -j,j,2
      if (allocated(class_mml(m2,m1,j)%vecall)) deallocate(class_mml(m2,m1,j)%vecall)
                                                allocate(class_mml(m2,m1,j)%vecall(1:3*imm_neigh))
    end do
  end do
end do
if (allocated(cmm)) deallocate(cmm)
                    allocate(  cmm(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )

if (weighted) then
  if (allocated(cmm2)) deallocate(cmm2)
                       allocate(  cmm2(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )
  do j= 0, jj_max
    do m1 = -j,j,2
      do m2 = -j,j,2
        if (allocated(class_mml(m2,m1,j)%vecall_w)) deallocate(class_mml(m2,m1,j)%vecall_w)
                                                    allocate(class_mml(m2,m1,j)%vecall_w(1:3*imm_neigh))
      end do
    end do
  end do

  if (weighted_3ch) then
    if (allocated(cmm3)) deallocate(cmm3)
                         allocate( cmm3(-jj_max:jj_max,-jj_max:jj_max, 0:jj_max) )
    do j= 0, jj_max
      do m1 = -j,j,2
        do m2 = -j,j,2
          if (allocated(class_mml(m2,m1,j)%vecall_w_3ch)) deallocate(class_mml(m2,m1,j)%vecall_w_3ch)
                                                          allocate(class_mml(m2,m1,j)%vecall_w_3ch(1:3*imm_neigh))
        end do
      end do
    end do
  end if

end if


end subroutine init_pow_so4
