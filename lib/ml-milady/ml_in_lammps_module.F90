
module module_so3
  use T_kind_param_m, only: kind_double, kind_double_complex

  integer, parameter   :: radial_bartok = 1, &
                          radial_sgg = 2
  integer  :: radial_pow_so3, ini_rbf_so3

  complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn
  complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn
  complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn_w
  complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn_w
  complex(kind_double_complex), allocatable, dimension(:, :, :, :)     :: dclmn_w_3ch
  complex(kind_double_complex), allocatable, dimension(:, :)     :: clmn_w_3ch
  real(kind_double), dimension(:), allocatable :: coeff_rbf_so3
  real(kind_double), dimension(:, :), allocatable    :: W_pow_so3

  real(8), dimension(:), allocatable     :: chebT, d_chebT                  ! T Pftouny first kind
  real(8), dimension(:), allocatable     :: chebU  ! T Pftouny second kind
  real(8), dimension(:), allocatable     :: scgg, d_scgg
end module module_so3



module module_bispectrum_so4
  use T_kind_param_m, ONLY:  double
  implicit none
  double complex, parameter :: czero=cmplx(0.d0,0.d0, kind=kind(1.d0))
  integer :: bisso4_cg_dim, bisso4_cg_full_dim, bisso4_cg_A_dim,bisso4_cg_B_dim, bisso4_cg_C_dim
  double complex, allocatable, dimension(:,:,:) :: cmm, cmm2, cmm3, &
                                                   ZAcmm, ZAcmm2, ZAcmm3, &
                                                   ZBcmm, ZBcmm2, ZBcmm3, &
                                                   ZCcmm, ZCcmm2, ZCcmm3


  double precision, allocatable, dimension(:) :: cg_A, cg_B, cg_C
  real(double) :: fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch

  type mml4D
    double complex, dimension(:), allocatable :: vecall, vecall_w, vecall_w_3ch
  end type

  type(mml4D), dimension(:,:,:), allocatable, target  :: class_mml

end module module_bispectrum_so4

 module module_kernel
      USE T_kind_param_m, ONLY:  double
      integer :: dim_kernel

      !building kernel du pauvre. Get evrything fron one class.
      logical :: write_kernel_matrix
      character(len=2), dimension(:), allocatable::  db_classes_for_kernel
      character(len=80) :: classes_for_kernel
      character(len=2), dimension(:), allocatable::  classes_train_for_kernel, classes_full_for_kernel

      real(double) :: kernel_power

      integer :: kernel_type ! kernel type, one of those below ...
      integer, parameter ::  kernel_se = 1,  &   ! square-exp - Gaussian
                             kernel_ou = 2,  &   !  ornstein_uhlenbeck
                             kernel_mc = 3,  &   ! matern_class
                             kernel_po = 4,  &   ! polynomial
                             kernel_mcd = 5, &   ! by MCD
                             kernel_maha = 6,&   ! by MAHA
                             kernel_random = 7, & ! random Gaussian
                             kernel_random_po = 44 ! random Poly

      integer :: kernel_dump, np_kernel_ref, np_kernel_full, np_omega
      integer, parameter :: kernel_dump_by_class=1, &
                            kernel_dump_by_mcd=2, &
                            kernel_dump_by_mahalanobis=3

      ! parameters of the kernel ...
      real(double)            :: time_for_desc_kernel
      real(double)            :: train_time_for_desc_kernel
      real(double)            :: test_time_for_desc_kernel
      real(double)            :: length_kernel
      real(double)            :: sigma_kernel
      real(double), parameter :: l_kou=1.d0
      !where the selected vectors are stores. To be ScaLpack in near future ...
      !if of dimension dim_descriptor x dim_kernel = (dim_desc, dim_kernel)

      real(double) :: power_mcd
      real(double), dimension(:,:), allocatable :: global_kernel, draft_kernel
      real(double), dimension(:,:), allocatable :: maha_kernel
      real(double), dimension(:), allocatable :: maha_norm_kernel

      !tenporary, vectors ...
      real(double), dimension(:), allocatable :: tmp_kernel, min_ker, max_ker, mean_ker, var_ker

      real(double), dimension(:), allocatable :: kernel_phase_random

      type basis_random
        ! desc dimension
        integer :: dim_xdesc
        ! dim of omega sampling
        integer :: dim_omega
        ! omega's matrix dim_xdesc x dim_omega
        real(double), dimension(:,:), allocatable :: omega
      end type basis_random
      type (basis_random), dimension(:), allocatable :: basis_random_po


 end module module_kernel


module ml_in_lammps_module
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
!      use tab_imm_m
!      use var_pot
!      use jqmod
!#if(PARAML)
!      use mpi
!      use mod_mpi_ml
!#endif
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none

      integer :: snap_order
      integer, parameter :: snap_linear = 1 , snap_quadratic = 2, snap_kernel = 7
      integer :: dim_xdesc

      real(double), parameter :: two_pi=8.d0*atan(1.d0)
      real(double), parameter :: one_pi=4.d0*atan(1.d0)
      real(double), dimension(3), parameter :: vec_3d_zero=(/0.d0, 0.d0, 0.d0/)
      real(double), parameter :: sqrt_two=sqrt(2.d0)
      double complex, parameter :: c_zero=(0.d0,0.d0)

      integer :: i_start_at,i_final_at
      real(double), dimension(:,:), allocatable :: design_mat  ! (D,M) D size of the database D values of x,
                                                               !  M the dimension of fingerprint ...
      real(double) :: eta_max_g2
      real(double), allocatable, dimension(:,:) :: g2, g3
      real(double), allocatable, dimension(:,:,:,:) :: g2_deriv, g3_deriv

      integer :: pow_so4_dim
      real(double), allocatable, dimension(:,:) :: pow_so4
      real(double), allocatable, dimension(:,:,:,:) :: pow_so4_deriv


      integer :: pow_so3_dim, l_max, n_rbf_so3,  n_rbf
      real(double), allocatable, dimension(:,:) :: pow_so3
      real(double), allocatable, dimension(:,:,:,:) :: pow_so3_deriv
      real(double), allocatable, dimension(:,:) :: W_pow_so3
      real(double), allocatable, dimension(:) :: coeff_rbf


      integer :: mtp_dim
      real(double), allocatable, dimension(:,:) ::  mtp
      real(double), allocatable, dimension(:,:,:,:) :: mtp_deriv


      integer :: soap_dim, nspecies_soap
      real(double), allocatable, dimension(:) ::  rb_soap
      real(double), allocatable, dimension(:,:) ::  soap, W_soap, S_factor_matrix
      real(double), allocatable, dimension(:,:,:,:) :: soap_deriv
      integer, allocatable, dimension(:,:) :: ns_soap_index
      integer :: n_soap
      logical :: lsoap,lsoap_fcut
      logical :: lsoap_diag, lsoap_norm, lsoap_lnorm
      double precision :: alpha_soap, atom_sigma_soap, r_cut_width_soap
      real(double), allocatable, dimension(:,:) :: k_soap
      real(double), allocatable, dimension(:,:,:) :: r_soap,at_soap
      integer, allocatable, dimension(:,:) :: iwmax2_soap,indi2_soap

      integer :: bisso4_dim
      integer, allocatable, dimension(:)  :: bisso4_l, bisso4_l1, bisso4_l2
      real(double), allocatable, dimension(:,:) :: bispectrum_so4
      real(double), allocatable, dimension(:,:,:,:) :: bispectrum_so4_deriv
      logical :: lbso4_diag
      real(double) :: j_max
      real(double) :: inv_r0, inv_r0_input
      integer :: jj_max

      integer:: afs_dim, n_cheb, n_rbf_afs
      integer :: afs_type
      integer, parameter :: afs_type_bartok=1, afs_type_homemade=2
      real(double), allocatable, dimension(:,:) :: afs
      real(double), allocatable, dimension(:,:,:,:) :: afs_deriv
      real(double), allocatable, dimension(:,:) :: W_afs
      real(double), allocatable, dimension(:) :: coeff_rbf_afs


      integer :: bisso3_dim
      logical :: lbso3_diag
      integer, allocatable, dimension(:)  :: bisso3_l, bisso3_l1, bisso3_l2
      real(double), allocatable, dimension(:,:) :: bispectrum_so3
      real(double), allocatable, dimension(:,:,:,:) :: bispectrum_so3_deriv


      double precision, allocatable, dimension(:,:,:,:,:,:) :: cg_vector


      real(double), allocatable, dimension(:,:) :: k_acd,distance_acd
      integer, allocatable, dimension(:) :: data_im
      integer, allocatable, dimension(:,:) :: data_natm
      real(double), allocatable, dimension(:,:,:) :: r_acd,at_acd
      integer, allocatable, dimension(:,:) :: iwmax2_acd,indi2_acd
      integer,dimension(:),allocatable :: natm

      real(double)  :: r_factorial(0:167)
      ! new version for g2 and g3
      integer :: n_g2_eta, n_g2_rs, g2_dim
      real(double) ,dimension(:),allocatable :: g2_eta, g2_rs

      integer :: n_g3_eta,  n_g3_zeta, n_g3_lambda, g3_dim
      real(double) ,dimension(:),allocatable :: g3_eta,  g3_zeta, g3_lambda
      integer :: behler_dim

      integer, parameter :: imm_neigh=100
      integer :: nd_data         ! dimension of the database M, internal
      integer :: nd_fingerprint  ! dimension of the feature space D, internal
      integer :: ml_type         ! input
      integer :: iread_ml        !input
      integer :: isave_ml        !input
      integer :: dim_kxx
      integer :: descriptor_type ! input
      integer, parameter ::  descriptor_g2=1, &
                             descriptor_g3=2, &
                             descriptor_behler=3, &
                             descriptor_afs=4, &
                             descriptor_soap=5, &
                             descriptor_pow_so3=6, &
                             descriptor_pow_so3_3body=603, &
                             descriptor_bispectrum_so3=7, &
                             descriptor_pow_so4=8, &
                             descriptor_g2_pow_so4=18, &
                             descriptor_bispectrum_so4=9, &
                             descriptor_g2_afs=14, &
                             descriptor_g2_bispectrum_so4=19, &
                             descriptor_mtp=100

      integer, parameter :: ml_type_krr=1, &
                            ml_type_basis=0, &
                            ml_type_gap=2,  &
                            ml_type_descriptors=-1





      integer :: target_type,force_comp
      integer, parameter ::  target_energy=1, &
                             target_force=2

      double precision,dimension(2) :: massat
      real(double) :: factor_weight_mass

      integer :: kernel_type, selection_type, ns_data,i_begin,kelem,seed,max_data
      real(double) :: lambda_krr, min_lambda_krr, max_lambda_krr
      integer :: n_values_lambda_krr
      character(len=4) :: regularization_name
      real(double), dimension(:), allocatable :: vector_lambda_krr
      logical :: toy_model, krr_error, debug, build_subdata,write_desc,weighted, weighted_3ch, rescale,search_hyp
      logical :: marginal_likelihood, desc_forces
      double precision :: n_frac,r0,s_max_r,s_min_r,s_max_i,s_min_i
      integer :: lpath
      character(len=60) :: path
      character(len=2) :: pref
      integer :: i_poscar
      integer :: kappa_acd,mc_step
      double precision :: alpha_acd,acd_threshold,ksi_ini,temp_ini,tau
      logical ::acd,acd_fcut,acd_weighted,sparsification,sparsification_by_acd,sparsification_by_entropy
      logical,dimension(:),allocatable :: reject
      integer :: size_indi2
      !database related
      ! iconf_data = iconf_data_train + iconf_data_test
      integer :: iconf_data, iconf_data_train, iconf_data_test
      integer, dimension(:), allocatable :: db_train, db_test

      character(len=80) :: db_file, db_path
      character(len=4)  :: char_desc
      integer, parameter ::  selection_type_first=1, &
                             selection_type_last=2, &
                             selection_type_random=3
      logical :: strict_behler
      integer :: dim_energy, dim_force
      !md_iconf is the corresponding iconf for molecular dynamics and is 1
      integer :: md_iconf
      real(kind(0.d0)) :: sign_stress, sign_stress_big_box
      integer :: snap_fit_type
      character(len=2) :: snap_class_constraints

      integer :: mtp_poly_min, mtp_poly_max, mtp_rad_order
      real(double) :: tmp_val_desc_max, val_desc_max

      ! weights optimization
      character(len=2), dimension(:), allocatable :: no_class_weights
      real(double) :: factor_force_error, factor_energy_error, factor_stress_error, timetot, timetot00, &
                      timetot01, timetot02, timetot03


contains


subroutine gen_param_behler()

   implicit none
  integer :: p_e,p_r, p_l, p_z
  integer :: icount
 ! G2 related ........
  if ((descriptor_type==descriptor_g2) .or. (descriptor_type==descriptor_g2_bispectrum_so4) &
                                       .or.(descriptor_type==descriptor_g2_afs) &
                                       .or.(descriptor_type==descriptor_g2_pow_so4)) then
      g2_dim=n_g2_eta*n_g2_rs
      if (allocated(g2_eta)) deallocate(g2_eta) ; allocate(g2_eta(g2_dim))
      if (allocated(g2_rs )) deallocate(g2_rs ) ; allocate(g2_rs (g2_dim))

       icount=0
       do p_e=1,n_g2_eta
         do p_r=1,n_g2_rs
            icount = icount + 1
            g2_eta(icount)=1.d-2+dble(p_e-1)*(eta_max_g2-1.d-2)/dble(n_g2_eta-1)
            g2_rs(icount)=0.d0+dble(p_r-1)
         enddo
       enddo
       if (icount /= g2_dim) then
         if (rangml==0) write(6,*) 'Error in ml_in_ndm_module.F90, subroutine gen_param_behler'
         stop 'dimension of g2 is incorrect'
       end if

      end if


  ! G3 related ........
  if (descriptor_type==descriptor_g3) then
      g3_dim=n_g3_eta*n_g3_zeta*n_g3_lambda
      if (strict_behler) then
         g3_dim=43
      end if
      if (allocated(g3_eta)   ) deallocate(g3_eta   ) ; allocate(g3_eta(g3_dim))
      if (allocated(g3_zeta)  ) deallocate(g3_zeta  ) ; allocate(g3_zeta(g3_dim))
      if (allocated(g3_lambda)) deallocate(g3_lambda) ; allocate(g3_lambda(g3_dim))


       icount=0
       do p_e=1,n_g3_eta
         do p_l=1,n_g3_lambda
            do p_z=1,n_g3_zeta
                icount=icount+1
                g3_eta(icount)=1.d-2+dble(p_e-1)*(0.80d0-1.d-2)/dble(n_g3_eta-1)
                g3_lambda(icount)=-1.d0+dble(p_l-1)*2.d0
                g3_zeta(icount)=2.d0**(p_z-1)
            enddo
         enddo
       enddo


       if (icount /= g3_dim) then
         if (rangml==0) write(6,*) 'Error in ml_in_ndm_module.F90, subroutine gen_param_behler'
         stop 'dimension of g2 is incorrect'
       end if

  end if
  !Behler related ....

return
end subroutine gen_param_behler







subroutine prepare_factorial()
implicit none

  r_factorial(0:167)=(/          &
  1.d0,                          &
  1.d0,                          &
  2.d0,                          &
  6.d0,                          &
  24.d0,                         &
  120.d0,                        &
  720.d0,                        &
  5040.d0,                       &
  40320.d0,                      &
  362880.d0,                     &
  3628800.d0,                    &
  39916800.d0,                   &
  479001600.d0,                  &
  6227020800.d0,                 &
  87178291200.d0,                &
  1307674368000.d0,              &
  20922789888000.d0,             &
  355687428096000.d0,            &
  6.402373705728d+15,            &
  1.21645100408832d+17,          &
  2.43290200817664d+18,          &
  5.10909421717094d+19,          &
  1.12400072777761d+21,          &
  2.5852016738885d+22,           &
  6.20448401733239d+23,          &
  1.5511210043331d+25,           &
  4.03291461126606d+26,          &
  1.08888694504184d+28,          &
  3.04888344611714d+29,          &
  8.8417619937397d+30,           &
  2.65252859812191d+32,          &
  8.22283865417792d+33,          &
  2.63130836933694d+35,          &
  8.68331761881189d+36,          &
  2.95232799039604d+38,          &
  1.03331479663861d+40,          &
  3.71993326789901d+41,          &
  1.37637530912263d+43,          &
  5.23022617466601d+44,          &
  2.03978820811974d+46,          &
  8.15915283247898d+47,          &
  3.34525266131638d+49,          &
  1.40500611775288d+51,          &
  6.04152630633738d+52,          &
  2.65827157478845d+54,          &
  1.1962222086548d+56,           &
  5.50262215981209d+57,          &
  2.58623241511168d+59,          &
  1.24139155925361d+61,          &
  6.08281864034268d+62,          &
  3.04140932017134d+64,          &
  1.55111875328738d+66,          &
  8.06581751709439d+67,          &
  4.27488328406003d+69,          &
  2.30843697339241d+71,          &
  1.26964033536583d+73,          &
  7.10998587804863d+74,          &
  4.05269195048772d+76,          &
  2.35056133128288d+78,          &
  1.3868311854569d+80,           &
  8.32098711274139d+81,          &
  5.07580213877225d+83,          &
  3.14699732603879d+85,          &
  1.98260831540444d+87,          &
  1.26886932185884d+89,          &
  8.24765059208247d+90,          &
  5.44344939077443d+92,          &
  3.64711109181887d+94,          &
  2.48003554243683d+96,          &
  1.71122452428141d+98,          &
  1.19785716699699d+100,         &
  8.50478588567862d+101,         &
  6.12344583768861d+103,         &
  4.47011546151268d+105,         &
  3.30788544151939d+107,         &
  2.48091408113954d+109,         &
  1.88549470166605d+111,         &
  1.45183092028286d+113,         &
  1.13242811782063d+115,         &
  8.94618213078297d+116,         &
  7.15694570462638d+118,         &
  5.79712602074737d+120,         &
  4.75364333701284d+122,         &
  3.94552396972066d+124,         &
  3.31424013456535d+126,         &
  2.81710411438055d+128,         &
  2.42270953836727d+130,         &
  2.10775729837953d+132,         &
  1.85482642257398d+134,         &
  1.65079551609085d+136,         &
  1.48571596448176d+138,         &
  1.3520015276784d+140,          &
  1.24384140546413d+142,         &
  1.15677250708164d+144,         &
  1.08736615665674d+146,         &
  1.03299784882391d+148,         &
  9.91677934870949d+149,         &
  9.61927596824821d+151,         &
  9.42689044888324d+153,         &
  9.33262154439441d+155,         &
  9.33262154439441d+157,         &
  9.42594775983835d+159,         &
  9.61446671503512d+161,         &
  9.90290071648618d+163,         &
  1.02990167451456d+166,         &
  1.08139675824029d+168,         &
  1.14628056373471d+170,         &
  1.22652020319614d+172,         &
  1.32464181945183d+174,         &
  1.44385958320249d+176,         &
  1.58824554152274d+178,         &
  1.76295255109024d+180,         &
  1.97450685722107d+182,         &
  2.23119274865981d+184,         &
  2.54355973347219d+186,         &
  2.92509369349301d+188,         &
  3.3931086844519d+190,          &
  3.96993716080872d+192,         &
  4.68452584975429d+194,         &
  5.5745857612076d+196,          &
  6.68950291344912d+198,         &
  8.09429852527344d+200,         &
  9.8750442008336d+202,          &
  1.21463043670253d+205,         &
  1.50614174151114d+207,         &
  1.88267717688893d+209,         &
  2.37217324288005d+211,         &
  3.01266001845766d+213,         &
  3.8562048236258d+215,          &
  4.97450422247729d+217,         &
  6.46685548922047d+219,         &
  8.47158069087882d+221,         &
  1.118248651196d+224,           &
  1.48727070609069d+226,         &
  1.99294274616152d+228,         &
  2.69047270731805d+230,         &
  3.65904288195255d+232,         &
  5.01288874827499d+234,         &
  6.91778647261949d+236,         &
  9.61572319694109d+238,         &
  1.34620124757175d+241,         &
  1.89814375907617d+243,         &
  2.69536413788816d+245,         &
  3.85437071718007d+247,         &
  5.5502938327393d+249,          &
  8.04792605747199d+251,         &
  1.17499720439091d+254,         &
  1.72724589045464d+256,         &
  2.55632391787286d+258,         &
  3.80892263763057d+260,         &
  5.71338395644585d+262,         &
  8.62720977423323d+264,         &
  1.31133588568345d+267,         &
  2.00634390509568d+269,         &
  3.08976961384735d+271,         &
  4.78914290146339d+273,         &
  7.47106292628289d+275,         &
  1.17295687942641d+278,         &
  1.85327186949373d+280,         &
  2.94670227249504d+282,         &
  4.71472363599206d+284,         &
  7.59070505394721d+286,         &
  1.22969421873945d+289,         &
  2.0044015765453d+291,          &
  3.28721858553429d+293,         &
  5.42391066613159d+295,         &
  9.00369170577843d+297,         &
  1.503616514865d+300 /)

return
end subroutine prepare_factorial

end module ml_in_lammps_module




module   temporary_data_cov
implicit none

integer :: dim_xdesc, dim_xdesc_linear, dim_xdesc_quadratic, dim_xdesc_quadratic_only
integer :: dim_xdesc_kernel, dim_xdesc_kernel_only
integer :: dim_extra,  dim_train, dim_valid

end module temporary_data_cov
