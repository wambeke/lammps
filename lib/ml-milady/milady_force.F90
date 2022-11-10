!     Extern "C" declaration has the form:
!
!  void milady_force_(int *, int *,
!    int *, double *, int *, int *,
!    int *, double *, int *, int *, int *, int *,
!    double *, double *, double *, double *,
!		  double *, double *, int *);
!
!     Call from pair_milady.cpp has the form:
!
!    milady_force_(&i,&nmax,&eflag_either,&eflag_global,&eflag_atom,&vflag_atom,
!        &eng_vdwl,eatom,&ntype,type,fmap,&x[0][0],
!	       &numneigh[i],firstneigh[i],&numneigh_full[i],firstneigh_full[i],
!        &f[0][0],&vatom[0][0],&errorflag);
!
   ! a_type(1:ntype) = the array of species
   ! ntype = the number of species
   ! i = the atom
   !  eflag_either, eflag_global,
   ! eflag_atom   -> if the energy is computed
   ! vflag_atom   -> if the virial is computed
   ! eflag_global ->

subroutine milady_force(i, nmax,   &
          eflag_either, eflag_global, eflag_atom, vflag_atom, &
          eng_vdwl, eatom, ntype, a_type, fmap, map_mass, x, &
          numneigh, firstneigh, numneigh_full, firstneigh_full, &
          f, vatom, errorflag, escale, fscale)
      use T_kind_param_m, only: double
      use milady_data, only:  n_desc_type, k_desc_type, ninteractions, &
                              weights, local_desc, q_local_desc, local_desc_kernel, &
                              local_desc_deriv, q_local_desc_deriv, local_desc_kernel_deriv, &
                              ref_energy_per_atom
      use ml_in_lammps_module, only : i_start_at, i_final_at,  &
                                      descriptor_type, descriptor_afs, descriptor_soap, descriptor_g2, &
                                      descriptor_pow_so4, descriptor_bispectrum_so4, descriptor_pow_so3, &
                                      descriptor_pow_so3_3body, &
                                      afs_dim, &
                                      soap_dim, &
                                      g2_dim, &
                                      bisso4_dim, &
                                      rangml, imm, snap_order, snap_linear, snap_quadratic, snap_kernel, &
                                      weighted, weighted_3ch
      use temporary_data_cov, only : dim_xdesc, dim_xdesc_linear, dim_xdesc_quadratic, &
                                     dim_xdesc_quadratic_only, dim_xdesc_kernel, dim_xdesc_kernel_only
      use module_compute_kernel, only: compute_kernel, dim_kernel
      use milady_interfaces, only : compute_afs, compute_soap, compute_bispectrum_so4, compute_g2, &
                                    compute_pow_so3, compute_pow_so3_3body
      !use compute_pow_so3_mod, only: compute_pow_so3
      implicit none
! ----------in out variables------------------------------------!
      integer, intent(in) :: i, nmax, ntype, eflag_either, eflag_global, eflag_atom, vflag_atom
      real(kind=kind(1.d0)), dimension(3,nmax), intent(inout) :: x, f
      real(kind=kind(1.d0)), dimension(6,nmax), intent(inout) :: vatom
      real(kind=kind(1.d0)), dimension(6,nmax) :: st_temp
      real(kind=kind(1.d0)), intent(inout) ::   eng_vdwl
      real(kind=kind(1.d0)), dimension(nmax) ,  intent(inout) :: eatom
      integer, dimension(nmax), intent(in)  ::  a_type
      integer, dimension(ntype), intent(in) ::  fmap
      real(double), dimension(ntype), intent(in)  :: map_mass
      integer, intent(in)  :: numneigh, numneigh_full
      integer, dimension(numneigh), intent(in)  :: firstneigh
      integer, dimension(numneigh_full), intent(in)  ::  firstneigh_full
      real(double), intent(in) ::  fscale, escale
      integer, intent(inout) ::  errorflag
      real(double), dimension(:), allocatable :: tmp_ld, tmp_wd
! ----------------------local variables-----------------------!
      real(kind=kind(1.d0)) :: uu(3), rr, tmptmp
      integer j,elti,jn,cn, ci ! , eltj
      real(kind=kind(1.d0)) :: rtmp,etmp, q_etmp, k_etmp, tmp_ref_energy

      integer :: id1, id2, indx, indx_last, iitest
!--------------------------------------------------------------!
      if (rangml==0) write(6,*) 'enter in force ...', i
      errorflag = 0

      i_start_at=i
      i_final_at=i
      imm = nmax

      !call allocate_local_desc
      select case (descriptor_type)
            case (descriptor_afs)
                  if (weighted.and.(.not.(weighted_3ch))) then
                     if (2*afs_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors bisso4_dim vs. dim_xdesc",2i5)') bisso4_dim, dim_xdesc
                       stop 'milady_force bisso4_dim'
                     end if
                 end if
                  if ((weighted_3ch)) then
                    if (3*afs_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors bisso4_dim vs. dim_xdesc",2i5)') bisso4_dim, dim_xdesc
                       stop 'milady_force bisso4_dim'
                    end if
                  end if
                  call compute_afs(i_start_at, i_final_at,  nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)
            case (descriptor_soap)
                  if (soap_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors afs_dim vs. dim_xdesc",2i5)') afs_dim, dim_xdesc
                       stop 'milady_force soap_dim'
                  end if


                  call compute_soap(i_start_at, i_final_at,  nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)

            case (descriptor_g2)
                  if (g2_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors g2_dim vs. dim_xdesc",2i5)') afs_dim, dim_xdesc
                       stop 'milady_force g2_dim'
                  end if


                  call compute_g2(i_start_at, i_final_at,  nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)


            case (descriptor_bispectrum_so4)
                  if (weighted.and.(.not.(weighted_3ch))) then
                     if (2*bisso4_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors bisso4_dim vs. dim_xdesc",2i5)') bisso4_dim, dim_xdesc
                       stop 'milady_force bisso4_dim'
                     end if
                 end if
                  if ((weighted_3ch)) then
                    if (3*bisso4_dim /= dim_xdesc) then
                       write(*,'("ML: serious problems in milady descriptors bisso4_dim vs. dim_xdesc",2i5)') bisso4_dim, dim_xdesc
                       stop 'milady_force bisso4_dim'
                    end if
                  end if
                  call compute_bispectrum_so4(i_start_at, i_final_at, nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)
            case (descriptor_pow_so4)
                  call compute_pow_so4(i_start_at, i_final_at, nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)
            case (descriptor_pow_so3)
                  call compute_pow_so3(i_start_at, i_final_at, nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)

            case (descriptor_pow_so3_3body)
                  call compute_pow_so3_3body(i_start_at, i_final_at, nmax, x, ntype, a_type, map_mass, &
                                     numneigh_full, firstneigh_full)



      end select
      tmp_ref_energy = ref_energy_per_atom(a_type(i_start_at))
!     Compute forces atom i
      elti = fmap(a_type(i))
      if (elti.gt.0) then
        rtmp = 0.0
        cn = 0
        !energy flag  per atom ... total energy lower
        !etmp = weights(1,ninteractions) * float(nmax) + dot_product(local_desc(1:dim_xdesc),weights(2:1+dim_xdesc, ninteractions))
        if (snap_order >= snap_linear) then
            etmp = weights(1,ninteractions) + dot_product(local_desc(1:dim_xdesc),weights(2:1+dim_xdesc, ninteractions))
            tmptmp = 0.d0
            do iitest=1, dim_xdesc
              tmptmp=tmptmp + local_desc(iitest)*weights(1+iitest, ninteractions)
            !write(*,'(a2, 2i7,4f20.10)') 'f ', i, iitest, local_desc(iitest), weights(1+iitest, ninteractions),  &
            !            local_desc(iitest)*weights(1+iitest, ninteractions), tmptmp
            end do
        end if
        !write(*,'(a5,i8,2f20.10)') ' ftot', i, etmp, weights(1, ninteractions)

        q_etmp=0.d0
        if (snap_order == snap_quadratic) then
            indx_last = 1 + dim_xdesc
            do id1 = 1, dim_xdesc
                  do id2 = 1, dim_xdesc
                        indx =dim_xdesc*(id1-1) + id2
                        q_local_desc(indx) = local_desc(id1)*local_desc(id2)
                  end do
            end do
            q_etmp = dot_product( q_local_desc(1:dim_xdesc_quadratic_only), &
                                  weights(dim_xdesc_linear+1:dim_xdesc_quadratic , ninteractions) )
        end if

        k_etmp=0.d0
        if (snap_order == snap_kernel) then
            call compute_kernel(i_start_at, i_final_at)
            k_etmp = dot_product(local_desc_kernel(1:dim_xdesc_kernel_only), &
                        weights(dim_xdesc_linear+1:dim_xdesc_kernel, ninteractions))
        end if


        etmp = (etmp + q_etmp + k_etmp + tmp_ref_energy) * escale
        !write(*,*) i, etmp,  weights(1,ninteractions),  nmax , dot_product(local_desc(1:dim_xdesc),weights(2:1+dim_xdesc, ninteractions))
        if (eflag_atom.ne.0) then
          eatom(i) =  etmp
        endif

        if (eflag_either.ne.0) then
          if (eflag_global.ne.0) eng_vdwl = eng_vdwl + etmp
          if (eflag_atom.ne.0) then
            eatom(i) = etmp
          endif
        endif

      if (snap_order >= snap_linear) then
        ! those are the forces on the atom i ... the first part ... those in first  r_cut shell
        do ci=1,3
          f(ci,i) =  f(ci,i) - dot_product(local_desc_deriv(1:dim_xdesc,0,ci), weights(2:1+dim_xdesc, ninteractions)) * fscale
        end do

        !the other part of forces .... on atom j that was in the list of atom i ... not complete but lammps will do the job to complete tgrough mpi.
        do jn = 1, n_desc_type
          j = k_desc_type(jn)
          do ci=1,3
            f(ci, j)=f(ci,j) -  dot_product(local_desc_deriv(1:dim_xdesc,jn,ci), weights(2:1+dim_xdesc, ninteractions)) * fscale
          end do
        end do
      end if

      if (snap_order == snap_kernel) then
        ! those are the forces on the atom i ... the first part ... those in first  r_cut shell
        if (allocated(tmp_ld)) deallocate(tmp_ld) ; allocate(tmp_ld(dim_kernel))
        if (allocated(tmp_wd)) deallocate(tmp_wd) ; allocate(tmp_wd(dim_kernel))
        tmp_wd(1:dim_kernel) = weights(dim_xdesc_linear+1:dim_xdesc_kernel, ninteractions)
        do ci=1,3
          tmp_ld(1:dim_kernel) = local_desc_kernel_deriv(1:dim_xdesc_kernel_only,0,ci)
          !f(ci,i) =  f(ci,i) - dot_product(local_desc_kernel_deriv(:,0,ci), &
          !                    weights(dim_xdesc+1:dim_xdesc_kernel, ninteractions)) * fscale
          f(ci,i) = f(ci,i) - dot_product(tmp_ld(:),tmp_wd(:)) * fscale
        end do

        !the other part of forces .... on atom j that was in the list of atom i ... not complete but lammps will do the job to complete tgrough mpi.

        do jn = 1, n_desc_type
          j = k_desc_type(jn)
          do ci=1,3
            f(ci, j)=f(ci,j) -  dot_product(local_desc_kernel_deriv(1:dim_xdesc_kernel_only,jn,ci), &
                              tmp_wd(1:dim_kernel)) * fscale
          end do
        end do
      end if

      if (snap_order == snap_quadratic) then
            indx_last = 1 + dim_xdesc
            indx = 0
            do id1=1, dim_xdesc
              do id2 = 1, dim_xdesc
                indx =dim_xdesc*(id1-1) + id2
                q_local_desc_deriv(indx,:) = local_desc_deriv(id1,0,:)*local_desc(id2) + &
                                             local_desc_deriv(id2,0,:)*local_desc(id1)
              end do
            end do
! those are the forces on the atom i ... the first part ... those in first  r_cut shell
            do ci=1,3
              f(ci,i) = f(ci,i) - dot_product( q_local_desc_deriv(1:dim_xdesc_quadratic_only,ci), &
                                               weights(dim_xdesc_linear+1:dim_xdesc_quadratic, ninteractions) ) * fscale
            end do

!the other part of forces .... on atom j that was in the list of atom i ... not complete but lammps will do the job to complete tgrough mpi.
            do jn = 1, n_desc_type
              j = k_desc_type(jn)
              indx_last = 1 + dim_xdesc
              indx = 0
              do id1=1, dim_xdesc
                do id2 = 1, dim_xdesc
                  indx =dim_xdesc*(id1-1) + id2
                  q_local_desc_deriv(indx,:) = local_desc_deriv(id1,jn,:)*local_desc(id2) + &
                                               local_desc_deriv(id2,jn,:)*local_desc(id1)
                end do
              end do

              do ci=1,3
                f(ci, j) = f(ci,j) - dot_product( q_local_desc_deriv(1:dim_xdesc_quadratic_only,ci) , &
                                                 weights(dim_xdesc_linear+1:dim_xdesc_quadratic, ninteractions) ) * fscale
              end do
            end do
      end if

!put the stress ... one day
        if (vflag_atom.ne.0) then
            st_temp(:,:) = 0.d0
            do jn = 1, n_desc_type
                  j = k_desc_type(jn)

                  uu(1:3) = x(1:3,i) - x(1:3,j)
                  rr = dsqrt( Sum( uu(1:3)**2 ) )


                  if (snap_order == snap_quadratic) then
                    indx_last = 0
                    do id1 = 1 , dim_xdesc
                      do id2 = 1 , dim_xdesc
                        indx =  dim_xdesc*(id1-1) + id2
                        q_local_desc_deriv(indx,1:3) = local_desc_deriv(id1,jn,1:3)*local_desc(id2) + &
                                                             local_desc(id1)*local_desc_deriv(id2,jn,1:3)
                      end do
                    end do
                  end if


      if (snap_order >= snap_linear) then
        st_temp(1:dim_xdesc, 1) = st_temp(1:dim_xdesc, 1) + local_desc_deriv(1:dim_xdesc,jn,1)*uu(1)/rr
        st_temp(1:dim_xdesc, 2) = st_temp(1:dim_xdesc, 2) + local_desc_deriv(1:dim_xdesc,jn,2)*uu(2)/rr
        st_temp(1:dim_xdesc, 3) = st_temp(1:dim_xdesc, 3) + local_desc_deriv(1:dim_xdesc,jn,3)*uu(3)/rr
        st_temp(1:dim_xdesc, 4) = st_temp(1:dim_xdesc, 4) + local_desc_deriv(1:dim_xdesc,jn,2)*uu(3)/rr
        st_temp(1:dim_xdesc, 5) = st_temp(1:dim_xdesc, 5) + local_desc_deriv(1:dim_xdesc,jn,1)*uu(3)/rr
        st_temp(1:dim_xdesc, 6) = st_temp(1:dim_xdesc, 6) + local_desc_deriv(1:dim_xdesc,jn,1)*uu(2)/rr
      end if


      if (snap_order == snap_kernel) then
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 1) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 1) + local_desc_kernel_deriv(1:dim_kernel,jn,1)*uu(1)/rr
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 2) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 2) + local_desc_kernel_deriv(1:dim_kernel,jn,2)*uu(2)/rr
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 3) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 3) + local_desc_kernel_deriv(1:dim_kernel,jn,3)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 4) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 4) + local_desc_kernel_deriv(1:dim_kernel,jn,2)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 5) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 5) + local_desc_kernel_deriv(1:dim_kernel,jn,1)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 6) = &
          st_temp(dim_xdesc+1:dim_xdesc_kernel-1 , 6) + local_desc_kernel_deriv(1:dim_kernel,jn,1)*uu(2)/rr
      end if


      if (snap_order == snap_quadratic) then
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 1) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,1) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,1)*uu(1)/rr
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 2) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,2) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,2)*uu(2)/rr
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 3) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,3) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,3)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 4) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,4) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,2)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 5) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,5) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,1)*uu(3)/rr
        st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, 6) = &
          st_temp(dim_xdesc+1:dim_xdesc_quadratic-1,6) + q_local_desc_deriv(dim_xdesc+1:dim_xdesc_quadratic-1,1)*uu(2)/rr
      end if


            end do   ! jn over the neighbours ...

            if (snap_order >= snap_linear) then
              do ci = 1,6
                vatom(ci,i) = vatom(ci,i) + &
                    dot_product( st_temp(1:dim_xdesc, ci), weights(2:dim_xdesc+1, ninteractions) )/3.d0 * fscale
              enddo
            end if


            if (snap_order == snap_kernel) then
              do ci = 1,6
                vatom(ci,i) = vatom(ci,i) + &
                    dot_product( st_temp(dim_xdesc+1:dim_xdesc_kernel-1, ci), &
                                 weights(dim_xdesc+2:dim_xdesc_kernel, ninteractions) )/3.d0 * fscale
              enddo
            end if

            if (snap_order == snap_quadratic) then
              do ci = 1,6
                vatom(ci,i) = vatom(ci,i) + &
                    dot_product( st_temp(dim_xdesc+1:dim_xdesc_quadratic-1, ci), &
                                 weights(dim_xdesc+2:dim_xdesc_quadratic, ninteractions) )/3.d0 * fscale
              enddo
            end if
        endif ! vflas for stress

      end if  !elti ....


      if (rangml==0) write(6,*)'ML: after pair'
      return
end subroutine milady_force
