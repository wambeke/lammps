
! check my fortran please!

subroutine milady_setup(nchannels,    nel,    n_nint, nwei,  &
                        extweights, extrcutvec, ndet, &
                        extdetails, extweight_factor, extweight_factor_3ch, ext_ref_energy_per_atom,  desc_type)
      use T_kind_param_m, ONLY:  double
      use milady_data, only: nelements, nweights, ninteractions, ndetails, ext_descriptor_type, weight_factor, &
                             weights,details, rcutvec, r_cut_desc,  weight_per_type, weight_per_type_3ch, ref_energy_per_atom
      use ml_in_lammps_module, only: descriptor_type, descriptor_afs, descriptor_soap, descriptor_g2, &
                                     descriptor_pow_so4, descriptor_bispectrum_so4, descriptor_pow_so3, &
                                     descriptor_pow_so3_3body, &
                                     n_rbf_afs, n_cheb, &
                                     n_soap, l_max, alpha_soap, r_cut_width_soap, lsoap_diag, &
                                     lsoap_norm, lsoap_lnorm, nspecies_soap, &
                                     r0, inv_r0, inv_r0_input, one_pi, lbso4_diag, j_max, jj_max, &
                                     eta_max_g2, n_g2_eta, n_g2_rs, &
                                     rangml, &
                                     desc_forces, weighted, weighted_3ch, factor_weight_mass, &
                                     prepare_factorial, snap_order, snap_kernel, afs_type, &
                                     l_max, n_rbf_so3, &
                                     timetot, timetot00, timetot01, timetot02, timetot03
      use temporary_data_cov, only: dim_xdesc, dim_xdesc_kernel
      use module_kernel, only: dim_kernel, kernel_type, global_kernel, kernel_phase_random, kernel_random, &
                               kernel_po, sigma_kernel, length_kernel, kernel_power, &
                               kernel_random_po, basis_random_po
      use module_so3, only: radial_pow_so3
      implicit none
      integer, intent(in) :: nchannels, nel, n_nint, nwei, ndet, desc_type
      !real(double),dimension(1),  intent(in) ::extweight_factor
      real(double),  dimension(nel), intent(in) ::extweight_factor
      real(double),  dimension(nel), intent(in) ::extweight_factor_3ch
      real(double),  dimension(nel), intent(in) ::ext_ref_energy_per_atom
      real(double), dimension(nwei*n_nint), intent(in)::extweights
      real(double), dimension(n_nint), intent(in) ::extrcutvec
      real(double), dimension(ndet), intent(in) ::extdetails
      real(double), dimension(:), allocatable :: tmp_weights
      ! local variables ...
      integer :: i, i_inter, i_desc, read_dim_xdesc, ik, ii, id, icount, itmp
      rangml=1
      if (rangml==0) write(6,*) 'enter in setup ...'
      timetot = 0.d0
      timetot00 =0.d0
      timetot01 =0.d0
      timetot02 =0.d0
      timetot03 = 0.d0
! set global values
      !number of elements and weight_per_type ...
      if (allocated(weight_per_type))     deallocate(weight_per_type)     ; allocate(weight_per_type(nel))
      if (allocated(weight_per_type_3ch)) deallocate(weight_per_type_3ch) ; allocate(weight_per_type_3ch(nel))
      if (allocated(ref_energy_per_atom)) deallocate(ref_energy_per_atom) ; allocate(ref_energy_per_atom(nel))
      weight_per_type(1:nel) = extweight_factor(1:nel)
      weight_per_type_3ch(1:nel) = extweight_factor_3ch(1:nel)
      ref_energy_per_atom(1:nel) =  ext_ref_energy_per_atom(1:nel)
      if (nchannels == 1) then
            weighted=.false.
            weighted_3ch=.false.
      end if
      if (nchannels == 2) then
            weighted=.true.
            weighted_3ch=.false.
      end if
      if (nchannels == 3) then
            weighted=.true.
            weighted_3ch=.true.
      end if
      !number of elements
      nelements = nel
      !dimension of descriptor for energy
      nweights = nwei
      !number of interactions
      ninteractions = n_nint
      !number of details for that descriptor
      ndetails = ndet
      ext_descriptor_type = desc_type
      !descriptor details n_radial, n_angular ....

      ! the weigth_factor on each element.
      if (allocated(weight_factor)) deallocate(weight_factor) ; allocate(weight_factor(nel))
      weight_factor = extweight_factor
      factor_weight_mass = 1.0

      !weighted =.false.
      !do i = 1, size(weight_factor(:))
      !    if (dabs(weight_factor(i) - 1.d0) > 1.d-3) weighted=.true.
      !end do


! safe filling of weight array from C. This will be changed if kernel is read.
      if(allocated(weights)) deallocate(weights) ; allocate(weights(nweights,ninteractions))
      do i = 1,(ninteractions*nweights)
            i_inter = (i-1)/nweights
            i_desc = modulo(i-1,nweights)
            weights(1+i_desc,1+i_inter) = extweights(i)
      end do

      if(allocated(details)) deallocate(details) ; allocate(details(ndetails))
      do i = 1,ndetails
            details(i) = extdetails(i)
      end do

      if(allocated(rcutvec)) deallocate(rcutvec) ; allocate(rcutvec(ninteractions))
      do i =1, ninteractions
            rcutvec(i) = extrcutvec(i)
      end do

!cutoff of the descriptor ...
      !old_version  r_cut_desc=maxval(rcutvec(:))/2.d0
      r_cut_desc=maxval(rcutvec(:))

! fake rang of the processor ... waiting the good ones.


!fill the details of the descriptor and init them ...

      descriptor_type = ext_descriptor_type
      select case(descriptor_type)

      case(descriptor_afs)

        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel ) then
          if (size(extdetails) /= 7) then
            if (rangml==0) write(6,*) 'Milady: fatal error for AFS descriptor, the number of the details (9th line) is 4'
            if (rangml==0) write(6,*) 'Milady: snap order  afs_type n_rbf_afs n_cheb'
            stop 'milady_setup afs kernel'
          end if
        else
          if (size(extdetails) /= 4) then
            if (rangml==0) write(6,*) 'Milady: fatal error for AFS descriptor, the number of the details (9th line) is 4'
            if (rangml==0) write(6,*) 'Milady: snap order  afs_type n_rbf_afs n_cheb'
            stop 'milady_setup afs'
          end if
        end if

        afs_type = nint(extdetails(2))
        n_rbf_afs = nint (extdetails(3))
        n_cheb = nint (extdetails(4))

        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(5))
              read_dim_xdesc = nint (extdetails(6))
              dim_kernel = nint (extdetails(7))
        end if
        if (rangml==0) write(*,*) 'afs_type', snap_order, afs_type, n_rbf_afs , n_cheb

      case (descriptor_bispectrum_so4)

        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel ) then
          if (size(extdetails) /= 3+3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for BSO4 descriptor, the number of the details (9th line) is 4'
            stop 'milady_setup bso4 kernel'
          end if
        else
          if (size(extdetails) /= 3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for BSO4 descriptor, the number of the details (9th line) is 4'
            stop 'milady_setup bso4'
          end if
        end if

        inv_r0_input=(1.d0 - 0.02d0/one_pi)
        inv_r0=one_pi*inv_r0_input
        r0 = r_cut_desc / inv_r0
        j_max = extdetails(2)
        jj_max=nint(2*j_max)
        !edbug write(*,*) 'rrrrrrrrrr', j_max, jj_max
        lbso4_diag = .false. ; if (nint(extdetails(3))==1) lbso4_diag=.true.
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(4))
              read_dim_xdesc = nint (extdetails(5))
              dim_kernel = nint (extdetails(6))

        end if
        write(*,*) kernel_type, read_dim_xdesc, dim_kernel

      case (descriptor_pow_so4)
        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel) then
          if (size(extdetails) /= 2+3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO4 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso4 kernel'
          end if
        else
          if (size(extdetails) /= 2) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO4 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso4'
          end if
        end if
        inv_r0_input=(1.d0 - 0.02d0/one_pi)
        inv_r0=one_pi*inv_r0_input
        r0 = r_cut_desc / inv_r0
        j_max = extdetails(2)
        jj_max=nint(2*j_max)
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(3))
              read_dim_xdesc = nint (extdetails(4))
              dim_kernel = nint (extdetails(5))
        end if


      case (descriptor_pow_so3)
        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel) then
          if (size(extdetails) /= 4+3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO3 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso3 kernel'
          end if
        else
          if (size(extdetails) /= 4) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO3 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso3'
          end if
        end if
        l_max = nint(extdetails(2))
        n_rbf_so3 = nint(extdetails(3))
        radial_pow_so3 = nint(extdetails(4))
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(5))
              read_dim_xdesc = nint (extdetails(6))
              dim_kernel = nint (extdetails(7))
        end if



      case (descriptor_pow_so3_3body)
        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel) then
          if (size(extdetails) /= 4+3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO3 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso3 3body kernel'
          end if
        else
          if (size(extdetails) /= 4) then
            if (rangml==0) write(6,*) 'Milady: fatal error for PSO3 descriptor, the number of the details (9th line) is 2'
            stop 'milady_setup pso3 3body'
          end if
        end if
        l_max = nint(extdetails(2))
        n_rbf_so3 = nint(extdetails(3))
        radial_pow_so3 = nint(extdetails(4))
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(5))
              read_dim_xdesc = nint (extdetails(6))
              dim_kernel = nint (extdetails(7))
        end if

      case (descriptor_soap)
        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel) then
          if (size(extdetails) /= 8) then
            if (rangml==0) write(6,*) 'Milady: fatal error for SOAP descriptor, the number of the details (9th line) is 8'
            stop 'milady_setup soap kernel'
          end if
        else
          if (size(extdetails) /= 8) then
            if (rangml==0) write(6,*) 'Milady: fatal error for SOAP descriptor, the number of the details (9th line) is 8'
            stop 'milady_setup soap'
          end if
        end if
        n_soap = nint(extdetails(2))
        l_max = nint(extdetails(3))
        alpha_soap = extdetails(4)
        r_cut_width_soap = extdetails(5)
        lsoap_diag = .false. ; if (nint(extdetails(6))==1) lsoap_diag=.true.
        lsoap_norm = .false. ; if (nint(extdetails(7))==1) lsoap_norm=.true.
        lsoap_lnorm = .false. ; if (nint(extdetails(8))==1) lsoap_lnorm=.true.
        !debug write(*,*) n_soap, l_max, alpha_soap, r_cut_width_soap, lsoap_diag, lsoap_norm, lsoap_lnorm
        !extra-stuff not yet included in params.pot
        nspecies_soap=1
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(9))
              read_dim_xdesc = nint (extdetails(10))
              dim_kernel = nint (extdetails(11))
        end if


      case(descriptor_g2)
        snap_order = nint(extdetails(1))
        if (snap_order == snap_kernel) then
          if (size(extdetails) /= 4+3) then
            if (rangml==0) write(6,*) 'Milady: fatal error for G2 descriptor, the number of the details (9th line) is 4'
            stop 'milady_setup g2'
          end if
        else
          if (size(extdetails) /= 4) then
            if (rangml==0) write(6,*) 'Milady: fatal error for G2 descriptor, the number of the details (9th line) is 4'
            stop 'milady_setup g2'
          end if
        end if
        eta_max_g2 = extdetails(2)
        n_g2_eta = nint (extdetails(3))
        n_g2_rs = nint (extdetails(4))
        if (snap_order == snap_kernel) then
              kernel_type = nint (extdetails(5))
              read_dim_xdesc = nint (extdetails(6))
              dim_kernel = nint (extdetails(7))
        end if


      case default
        if (rangml==0) write(6,*) 'No implementation for descriptor_type...', ext_descriptor_type
        stop "fatal in lib/milady/milady_setup.F90 descriptor"
      end select


      if (rangml==0) write(6,*) 'desc read  ...'
      desc_forces=.true.
      call prepare_factorial
      call init_descriptors
      if (rangml==0) write(6,*) 'desc init  ...'

      if (snap_order == snap_kernel) then
            if (dim_xdesc /= read_dim_xdesc) then
               write(5,*) 'Milady: fatal error dimension  descriptor is not the same as the from potential', &
                          dim_xdesc, read_dim_xdesc
               stop "fatal read_dim_xdesc"
            end if
            dim_xdesc_kernel = 1 + dim_xdesc + dim_kernel

            if (allocated(tmp_weights)) deallocate(tmp_weights) ; allocate (tmp_weights(dim_xdesc_kernel))
            if (ninteractions > 1 ) stop 'ninteractions should be 1 ... and only 1 in this version'
            if (allocated(global_kernel)) deallocate(global_kernel) ; allocate(global_kernel(dim_xdesc, dim_kernel))
            ninteractions = 1
            !write(*,*) 'HHHEERREEE', dim_xdesc, dim_kernel, dim_xdesc_kernel, size(weights,1), size(weights,2)
            tmp_weights(1:dim_xdesc_kernel) = weights(1:dim_xdesc_kernel,ninteractions)
            icount = dim_xdesc_kernel

            if (kernel_type /= kernel_random_po) then

              do ik = 1, dim_kernel
                do id = 1, dim_xdesc
                  icount = icount + 1
                  global_kernel(id,ik) = weights(icount,ninteractions)
                end do
              end do
            else
               global_kernel(:,:) = 0.d0
               if (allocated(basis_random_po)) deallocate(basis_random_po) ; allocate(basis_random_po(dim_kernel))
               do ik = 1, dim_kernel
                 icount = icount +1
                 itmp  = nint(weights(icount,ninteractions))
                 basis_random_po(ik)%dim_omega = itmp
                 if (allocated(basis_random_po(ik)%omega)) deallocate(basis_random_po(ik)%omega)
                 allocate(basis_random_po(ik)%omega(dim_xdesc, itmp ))
                 do ii = 1, itmp
                    do id = 1, dim_xdesc
                      icount = icount + 1
                      basis_random_po(ik)%omega(id,ii) = weights(icount,ninteractions)
                    end do
                 end do

               end do
            end if


            if (kernel_type == kernel_random) then
              if ( icount /=  (dim_xdesc_kernel + dim_xdesc*dim_kernel)) then
                write(6,*) 'Problems in reading kernel file ', icount, dim_xdesc_kernel + dim_xdesc*dim_kernel
                stop 'kernel kernel_random'
              end if
              if (allocated(kernel_phase_random)) deallocate(kernel_phase_random) ; allocate(kernel_phase_random(dim_kernel))
                do ik = 1, dim_kernel
                  kernel_phase_random(ik) = weights(icount + ik, ninteractions)
                end do
            end if

           if (kernel_type == kernel_po) then
              if ( icount /=  (dim_xdesc_kernel + dim_xdesc*dim_kernel)) then
                write(6,*) 'Problems in reading kernel file ', icount, dim_xdesc_kernel + dim_xdesc*dim_kernel
                stop 'kernel kernel_po'
              end if
              sigma_kernel = weights(icount + 1, ninteractions)
              length_kernel = weights(icount + 2, ninteractions)
              kernel_power = weights(icount + 3, ninteractions)
           end if

           if (kernel_type == kernel_random_po ) then
              sigma_kernel = weights(icount + 1, ninteractions)
              length_kernel = weights(icount + 2, ninteractions)
              kernel_power = weights(icount + 3, ninteractions)
           end if

            if(allocated(weights)) deallocate(weights) ; allocate(weights(1:dim_xdesc_kernel,ninteractions))
            weights(1:dim_xdesc_kernel,ninteractions) = tmp_weights(1:dim_xdesc_kernel)
            deallocate(tmp_weights)
      end if


      if (rangml==0) write(6,*) 'out of  setup ...'
!      stop 'setup stop'
      !rangml=1
!reallocate the parameters for kernel case as well as some tests:
! safe filling of weight array from C. This will be changed if kernel is read.
      return
end subroutine milady_setup


subroutine deallocate_milady
use milady_data
implicit none
!fill the details of the descriptor

if (allocated(k_desc_type)) deallocate(k_desc_type)
if (allocated(local_desc)) deallocate(local_desc)
if (allocated(local_desc_deriv)) deallocate(local_desc_deriv)

end subroutine deallocate_milady



subroutine  init_descriptors

use   ml_in_lammps_module, ONLY : rangml,  &
                               descriptor_type, descriptor_g2, descriptor_g3, descriptor_behler, &
                               descriptor_afs, descriptor_g2_afs, descriptor_soap, descriptor_pow_so3, descriptor_bispectrum_so3, &
                               descriptor_pow_so4, &
                               descriptor_bispectrum_so4, jj_max, one_pi, inv_r0_input, inv_r0, &
                               descriptor_g2_bispectrum_so4, descriptor_pow_so3_3body,  &
                               descriptor_mtp, &
                               gen_param_behler, j_max, &
                               g2_dim, g3_dim, bisso4_dim, char_desc, descriptor_g2_pow_so4,  &
                               pow_so4_dim, mtp_dim, afs_dim, lbso4_diag, n_rbf, &
                               pow_so3_dim, l_max, lbso3_diag, bisso3_dim, &
                               soap_dim, snap_order, snap_kernel, snap_linear, snap_quadratic, &
                               weighted, weighted_3ch, n_rbf_so3
use temporary_data_cov, only : dim_xdesc,  dim_xdesc_linear, dim_xdesc_quadratic, dim_xdesc_quadratic_only, &
                               dim_xdesc_kernel, dim_xdesc_kernel_only
use milady_data, only: allocate_local_desc
use module_kernel, only: dim_kernel
use module_so3, only: radial_pow_so3, radial_sgg, radial_bartok
implicit none

  if (rangml==0) write(6,'("ML: descriptor was initialized in init_descriptors...")')


select case(descriptor_type)

  case (descriptor_g2)
     call gen_param_behler
     dim_xdesc = g2_dim
     char_desc = 'bhg2'

  case (descriptor_g3)
     call gen_param_behler
     dim_xdesc = g3_dim
     char_desc = 'bhg3'

  case (descriptor_behler)
     call gen_param_behler
     dim_xdesc = g2_dim+g3_dim
     char_desc = 'bhlr'

  case (descriptor_soap)
     call compute_cg_vector(dble(l_max),1)
     call gen_dimension_for_soap()
     dim_xdesc=soap_dim
     char_desc='soap'

  case (descriptor_pow_so3)
    pow_so3_dim = int((1 + l_max))*n_rbf_so3
    if (radial_pow_so3 == radial_bartok) then
      pow_so3_dim = int((1 + l_max))*n_rbf_so3
    end if
    if (radial_pow_so3 == radial_sgg) then
      pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)
    end if
    call init_pow_so3
    call init_pow_so3_rbf()
    call compute_cg_vector(dble(l_max), 1)
    if (weighted) then
      dim_xdesc = 2*pow_so3_dim
      if (weighted_3ch) dim_xdesc = 3*pow_so3_dim
    else
      dim_xdesc = pow_so3_dim
    end if
    char_desc = 'pso3'

  case (descriptor_pow_so3_3body)
    pow_so3_dim = int((1 + l_max))*n_rbf_so3**2
    if (radial_pow_so3 == radial_bartok) then
      pow_so3_dim = int((1 + l_max))*n_rbf_so3**2
    end if
    if (radial_pow_so3 == radial_sgg) then
      pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)**2
    end if
    call init_pow_so3_3body
    call init_pow_so3_rbf()
    call compute_cg_vector(dble(l_max), 1)
    if (weighted) then
      dim_xdesc = 2*pow_so3_dim
      if (weighted_3ch) dim_xdesc = 3*pow_so3_dim
    else
      dim_xdesc = pow_so3_dim
    end if
    char_desc = 'psb3'

  case (descriptor_bispectrum_so3)
 !future    call init_pow_so3_rbf()
     call compute_cg_vector(dble(l_max),1)
     if (lbso3_diag) then
!future       call gen_dimension_for_bispectrum_so3_diagonal()
     else
!future       call gen_dimension_for_bispectrum_so3_all()
     end if
     dim_xdesc=bisso3_dim
     char_desc = 'bso3'


  case (descriptor_afs)
     call init_afs_rbf()
     if (weighted) then
       dim_xdesc=2*afs_dim
       if (weighted_3ch) dim_xdesc=3*afs_dim
     else
       dim_xdesc=afs_dim
     end if
     char_desc = 'afsr'

  case (descriptor_g2_afs)
     call gen_param_behler
     call init_afs_rbf()
     dim_xdesc = g2_dim+afs_dim
     char_desc = 'g2af'

  case (descriptor_g2_pow_so4)
     call gen_param_behler
     call compute_cg_vector(dble(j_max),2)
     dim_xdesc = g2_dim+int(2*j_max)+1
     pow_so4_dim = int(2*j_max)+1
     char_desc = 'g2p4'

  case (descriptor_pow_so4)
     jj_max=int(2*j_max)
     call compute_cg_vector(dble(j_max),2)
     pow_so4_dim = int(2*j_max)+1
     if (weighted) then
       dim_xdesc=2*pow_so4_dim
       if (weighted_3ch) dim_xdesc=3*pow_so4_dim
     else
       dim_xdesc=pow_so4_dim
     end if
     call init_pow_so4
     char_desc = 'pso4'


  case (descriptor_bispectrum_so4)

     jj_max = int(2*j_max)
     inv_r0=one_pi*inv_r0_input
     call compute_cg_vector(dble(j_max),2)
     call gen_dimension_for_bispectrum_so4()
     if (weighted) then
       dim_xdesc=2*bisso4_dim
       if (weighted_3ch) dim_xdesc=3*bisso4_dim
     else
       dim_xdesc=bisso4_dim
     end if
     char_desc = 'bso4'
     call test_parameters_bso4

  case (descriptor_mtp)
 !future    call gen_dimension_for_mtp()
     dim_xdesc=mtp_dim
     char_desc = 'mtp3'

  case (descriptor_g2_bispectrum_so4)
     call gen_param_behler
     call compute_cg_vector(dble(j_max),2)
     !snap call gen_dimension_for_bispectrum_so4_all()
     ! Gabor version for which are taken only (J J_1 J_1) componenets
     if (lbso4_diag) then
!future       call gen_dimension_for_bispectrum_so4_diagonal()
     else
!future       call gen_dimension_for_bispectrum_so4_all()
     end if
     dim_xdesc=bisso4_dim+g2_dim
     if (rangml==0)  write(6,'("ML: hybrid descritor in init desc  bso4 + g2:  ", 2i5)') bisso4_dim, g2_dim
     char_desc = 'g2b4'
  case default
       if (rangml==0) write(6,*) 'No implementation for descriptor_type...', descriptor_type
       stop "fatal in File: descriptor.f90, Subroutine: init_descriptor"
end select

  if (rangml==0) write(6,'("ML: descriptor ",a," has the dimension ",i6)') char_desc, dim_xdesc
  if (snap_order == snap_linear) then
      dim_xdesc_linear = 1 + dim_xdesc
      if (rangml==0) write(6,'("ML: The LML based on the descriptor ",a," has the dimension ",i6)') char_desc, dim_xdesc_linear
  end if

  if (snap_order == snap_quadratic) then
      dim_xdesc_linear = 1 + dim_xdesc
      dim_xdesc_quadratic = 1 + dim_xdesc + dim_xdesc**2
      dim_xdesc_quadratic_only = dim_xdesc**2
      if (rangml==0) write(6,'("ML: The quadratic ML based on ",a," has the dimension ",i6)') char_desc, dim_xdesc_quadratic
  end if

  if (snap_order == snap_kernel) then
      dim_xdesc_linear = 1 + dim_xdesc
      dim_xdesc_kernel = 1 + dim_xdesc + dim_kernel
      dim_xdesc_kernel_only = dim_kernel
      if (rangml==0) write(6,'("ML: The kernel ML based on ",a," has the dimension ",i6, i8)') char_desc, dim_xdesc, dim_kernel

  end if

  call allocate_local_desc
return
end subroutine  init_descriptors
