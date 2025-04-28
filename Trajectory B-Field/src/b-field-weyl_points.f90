!=============================================================
! This program perturbs the energy with a magnetic field
! and calculates the fermi surface.
!===============================================================
program generate_fermi
    use mpi
    use reading_module
    use fourier_module
    use perturbation
    use create_mesh
   ! #include <mpif.h>
    Implicit None
!--------Presets to be changed by User
    character(len=80):: prefix="BiTeI"
    !Adjust these parameters to obtain better resolution around alphacrit and see points closer to an effectively closed gap
    integer,parameter::np=5,npartitions=1,dim=3
         ! Flags
    logical :: useOptimized = .true.  ! Set to false for 4x4 case
    
    real*8,parameter::B_x = 0d0, B_y = 0.01d0, B_z = 0d0, dk = 0.000001d0!,& !Run again at B_y=0.05-0.06 to see the gap close
                    !   alpha_min = 0.50d0, alpha_max = 0.6d0 !4x4 range perturbed
                    !   alpha_min = 0d0, alpha_max = 1d0

                    !   alpha_min = 0.805d0, alpha_max = 0.81d0 !18x18 range  perturbed range

!---------MPI variables
    integer :: ierr, nprocs, rank, local_start, local_end, local_count
    real*8 :: mpi_start_time, mpi_end_time

!---------Variable allocation
    character(len=80) :: hamil_file_trivial, hamil_file_topological, nnkp, line, partnumber
    character(len=200) :: hamil_dir = '../Hamiltonians 18x18/'

    integer ik, ipart, ib, is,&
            i,j,k,kx,ky,kz,&
            n,nr_trivial, nr_topological,nb,&
            i1,i2,&
            lwork,info,&
            o,p,j1,j2,&
            total_pairs,&
            temp_index
    
    real*8 phase, dx, dy, dz,&
           bandgap,bandgapp,&
           twopi,jk,a,b,a1,b1,&
           spin_x(1,1),spin_y(1,1),spin_z(1,1),&
           spin_xp(1,1),spin_yp(1,1),spin_zp(1,1),&
           alpha,ef(npartitions),efp(npartitions),&
           gap(npartitions),gap_perturbed(npartitions),gap_unperturbed(npartitions),gap_min(npartitions),gapp_min(npartitions),&
           write_values(11:14),&
           bvec(3,3),avec(3,3),rvecs(3),&
           ktemp1(3),ktemp2(3),&
           kxtemp1(3),kxtemp2(3),&
           kytemp1(3),kytemp2(3),&
           kmesh(np,np),&
           mesh_kx(3,np, np), mesh_ky(3,np, np),&
           mesh_gap(3, np**2),&
           delkx, delky, delkz,&
           part_time,part_time2


    complex*16 chi(2,1),chip(2,1),&
               phi(3),phase_factor,overlap,&
							 overlap_x, overlap_y, overlap_z,&
							 phase_x, phase_y, phase_z

    real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
    
    integer,allocatable:: ndeg_trivial(:),ndeg_topological(:)

    real*8,allocatable:: phases(:,:),rvec_trivial(:,:),rvec_topological(:,:),&
                         ene(:,:),enep(:,:),eneBerry(:,:,:,:),&
												 ene_step(:,:), meshBerry(:,:,:,:),&
                         rwork(:),rvec(:,:),&
                         spin(:,:,:),spinp(:,:,:), magnitude_field(:,:,:)
                         
    complex*16,allocatable:: H(:,:), Hk(:,:), Hm(:,:), Hmag(:,:),&
														 H_step(:,:),Hk_step(:,:),&
														 u1(:), u0(:,:,:,:),&
														 connection(:,:,:,:),curvature(:,:,:,:),&
                             HK_trivial(:,:), HK_topological(:,:),&
                             Hamr_trivial(:,:,:), Hamr_topological(:,:,:),&
                             work(:), &
                             gauge_fixed_U_x(:), gauge_fixed_U_y(:), gauge_fixed_U_z(:)

    ! real*8, parameter :: dk=0.000001d0
    !----Box for weyl point----!
    real*8, parameter :: kbox_x=0.12d0,&!x_min = -0.06d0, x_max = 0.06d0,&
                         kbox_y=0.12d0,&!y_min = -0.06d0, y_max = 0.06d0,&
                         kbox_z=0.055d0

    ! real*8, parameter :: kbox_x=0.15d0,&!x_min = -0.06d0, x_max = 0.06d0,&
    !                      kbox_y=0.15d0,&!y_min = -0.06d0, y_max = 0.06d0,&
    !                      kbox_z=0.05d0!z_min = -0.03d0, z_max = 0.03d0

    integer, dimension(:), allocatable:: indices

    real*8, dimension(:,:), allocatable :: mesh, stepmesh
    !complex*16,dimension(2,2) :: sigx(2,2),sigy(2,2),sigz(2,2)


!--------- Initialize MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
    ! Start timing
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    mpi_start_time = MPI_WTIME()
    
    ! Only rank 0 prints informational messages
    if (rank == 0) then
        write(*,*) "Starting calculations with ", nprocs, " MPI processes"
    endif


    ! Modify the file path construction
    write(hamil_file_trivial,'(2a,a,a)') trim(adjustl(hamil_dir)), trim(adjustl(prefix)), "_hr_trivial.dat"
    write(hamil_file_topological,'(2a,a,a)') trim(adjustl(hamil_dir)), trim(adjustl(prefix)), "_hr_topological.dat"
    write(nnkp,'(2a,a,a)') trim(adjustl(hamil_dir)), trim(adjustl(prefix)), ".nnkp"
!-----Define pi
    twopi=4.0d0*atan(1.0d0)*2.0d0
! !--------- Define Pauli Matrices sigma x,y,z
!   data sigx /(0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
!   data sigy /(0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
!   data sigz /(1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/

! Only process 0 reads the input files
  if (rank == 0) then
!---------------  Read the vectors
		call read_lattice(nnkp, rank, avec, bvec, ierr)
!------Broadcast necessary values to all processes
		call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
		call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call read_header(hamil_file_trivial, nb, nr_trivial)
      write(*,*) nb,nr_trivial
    call read_header(hamil_file_topological, nb, nr_topological)
      write(*,*) nb,nr_topological
  endif

	call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(nr_trivial, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(nr_topological, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

	allocate(rvec_trivial(3, nr_trivial), rvec_topological(3, nr_topological))
	allocate(Hamr_trivial(nb, nb, nr_trivial), Hamr_topological(nb, nb, nr_topological))
	allocate(Hk_topological(nb,nb),Hk_trivial(nb,nb), u1(nb))
	allocate(H(nb,nb),Hk(nb,nb),H_step(nb,nb),Hk_step(nb,nb))
	allocate(ndeg_trivial(nr_trivial),ndeg_topological(nr_topological))

! Only rank 0 reads the Hamiltonian data
  if (rank == 0) then
    if (useOptimized) then
      if (rank == 0) then
        call read_optimized_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr_trivial, &
                                       	 Hamr_trivial, Hamr_topological, rvec_trivial, ndeg_trivial, avec)
        rvec_topological = rvec_trivial
      end if
    else
      call read_general_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr_trivial, nr_topological, &
                                    Hamr_trivial, Hamr_topological, rvec_trivial, rvec_topological, &
                                    ndeg_trivial, ndeg_topological, avec)
    endif
  endif

	allocate(enep(nb, (2*np+1)**dim), ene(nb, (2*np+1)**dim), ene_step(nb, (2*np+1)**dim), eneBerry(nb, 2*np+1, 2*np+1, 2*np+1))
! Broadcast the read data to all processes
    call MPI_BCAST(ndeg_trivial, nr_trivial, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ndeg_topological, nr_topological, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rvec_trivial, 3*nr_trivial, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rvec_topological, 3*nr_topological, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(hamr_trivial, nb*nb*nr_trivial, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(hamr_topological, nb*nb*nr_topological, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
!------ LAPACK-related array allocations
	lwork=max(1,2*nb-1)
	allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
! allocate(spin(3,nb,(np+1)**dim),spinp(3,nb,(np+1)**dim))

!------ open gap file
    ! open(777,file='gap.dat',status='replace', position='append', action='write')
    ! open(100,file='trajectory_4_0.001.dat',status='replace', position='append', action='write', iostat=ierr)
    ! if (ierr /= 0) then
    !     open(100, file='trajectory_4_0.001.dat', status='new', action='write')
    ! end if
    ! open(110,file='perturbed.dat',status='old', position='append', action='write')
    ! open(120,file='FERMISURFACE.dat', status='new', position='append',action='write')

!-------Generate Mesh
  ! allocate(mesh(3, (2*np+1)**dim))
	allocate(meshBerry(3, 2*np+1, 2*np+1, 2*np+1))
	! allocate(stepmesh(3, (2*np+1)**dim))
	allocate(U0(nb, 2*np+1, 2*np+1, 2*np+1))
	allocate(connection(3, 2*np+1, 2*np+1, 2*np+1))
	allocate(curvature(3, 2*np+1, 2*np+1, 2*np+1))
	allocate(magnitude_field(2*np+1, 2*np+1, 2*np+1))
  
	call berryLattice3D(np, dim, kbox_x, kbox_y, kbox_z, delkx, delky, delkz, meshBerry, bvec, j)
	! call Lattice3D(np, dim, kbox_x, kbox_y, kbox_z, delkx, delky, delkz, mesh, bvec, j)
  ! call SpecificLattice(mesh, np, -0.05d0, 0.05d0, 0.05d0, 0.05d0, 0.048d0, 0.050d0, bvec, j)


  if (rank == 0) then
  	write(*,*) "Total mesh points generated:", j
    write(*,*) "Expected mesh points:", (2*np+1)**dim
  endif
!------ Magnetic Field
  allocate(Hm(2,2), Hmag(nb,nb))
    ! call magnetic_field(nb, B_x, B_y, B_z, Hmag)

!----------- BEGIN INTERPOLATION -----------! 
!------ Fourrier transform H(R) to H(k)
    ! allocate(phases(nr, (np+1)**dim))

    local_count = npartitions / nprocs
    local_start = rank * local_count + 1
    local_end = (rank + 1) * local_count
    
    if (rank == nprocs - 1) then
        local_end = npartitions  ! Last process takes any remainder
    end if
    
    ! Each process handles its own partitions
    
    do ipart = local_start, local_end!1, npartitions
        if (rank == 0) then
            write(*,'(a,i5,a,i5)') 'Partition=', ipart, ' of ', npartitions
        endif
       !alpha=float(ipart-1)/float(npartitions-1)

        alpha=0.77966
        ! alpha = alpha_min + float(ipart - 1)*(alpha_max - alpha_min)/float(npartitions - 1)
        if (rank == 0) then
            write(*,'(A,I5,A,F12.6)') 'Partition ', ipart, ' alpha = ', alpha
        endif

        ! Initialize Hamiltonians for the current partition
    !    ene=0d0

! !----- FOURIER TRANSFORMS

!---4x4---!   
    !    call fourier_transform_general(2*np, dim, nr_trivial, nr_topological, nb, ndeg_trivial, ndeg_topological, mesh,&
    !                                   rvec_trivial, rvec_topological, Hamr_trivial, Hamr_topological,&
    !                                   Hmag, alpha, enep, ene, work, lwork, rwork, rank, ierr)
    !    gap_min(ipart)=abs(minval(ene(3,:))-maxval(ene(2,:))) !Unperturbed Case
    !    gapp_min(ipart)=abs(minval(enep(3,:))-maxval(enep(2,:))) !Perturbed Case
    !    ef(ipart)=(minval(ene(3,:))+maxval(ene(2,:)))/2d0 

!---18x18---!
    !    call fourier_transform_optimized(2*np, dim, nr_trivial, nb, ndeg_trivial, mesh, rvec_trivial, &
    !                                     Hamr_trivial, Hamr_topological, Hmag, alpha, &
    !                                     enep, ene, work, lwork, rwork, rank, ierr)

				!---Extra Fourier transform method for BC---!
							! do i=1,3
							! 	stepmesh = mesh
							! 	! Make sure mesh step sizes are correct

							! 	stepmesh(i, :) = stepmesh(i, :) + dk
							! 	call inner_ft_optimized(k, nr_trivial, nb, ndeg_trivial, stepmesh, rvec_trivial, &
							! 													Hamr_trivial, Hamr_topological, Hmag, alpha, &
							! 													Hk_trivial, Hk_topological, H_step, Hk_step, rank, ierr)
							! 	call zheev('V', 'U', nb, Hk_step, nb, ene_step(:, k), work, lwork, rwork, info)
							! 	u1=Hk_step(:,13)
							! 	overlap = dot_product(conjg(u0), u1)
							! 	u1=u1*(conjg(overlap)/abs(overlap))

							! 	connection(i, kx, ky, kz) = -aimag(dot_product(conjg(u0), (u1)/dk))
							! enddo
!--- Berry Transforms ---!
        do kx=1, 2*np+1
					do ky=1, 2*np+1
						do kz=1, 2*np+1

							! k = (kx-1) + (ky-1)*(2*np+1) + (kz-1)*(2*np+1)**2 + 1

              call inner_ft_berry(kx, ky, kz, nr_trivial, nb, ndeg_trivial, meshBerry, rvec_trivial, &
                                  Hamr_trivial, Hamr_topological, Hmag, alpha, &
                                  Hk_trivial, Hk_topological, H, Hk, rank, ierr)
							call zheev('V', 'U', nb, Hk, nb, eneBerry(:, kx, ky, kz), work, lwork, rwork, info)

							U0(:, kx, ky, kz)=Hk(:,13)
						enddo ! Exit k-loop
					enddo
				enddo 

				
        !---Finite Difference method
				do kx=1, 2*np
					do ky=1, 2*np
						do kz=1, 2*np
							! X direction
							overlap_x = dot_product((U0(:, kx, ky, kz)), U0(:, kx+1, ky, kz))
							phase_x = (overlap_x)/abs(overlap_x)
							gauge_fixed_U_x = U0(:, kx+1, ky, kz) * conjg(phase_x)
							connection(1, kx, ky, kz) = -aimag((dot_product(((U0(:, kx, ky, kz))),((gauge_fixed_U_x-U0(:, kx, ky, kz))/dk))))!delkx)))
              
							! Y direction
							overlap_y = dot_product((U0(:, kx, ky, kz)), U0(:, kx, ky+1, kz))
							phase_y = (overlap_y)/abs(overlap_y)
							gauge_fixed_U_y = U0(:, kx, ky+1, kz) * conjg(phase_y)
							connection(2, kx, ky, kz) = -aimag((dot_product(((U0(:, kx, ky, kz))),((gauge_fixed_U_y-U0(:, kx, ky, kz))/dk))))!delky)))
				
							! Z direction
							overlap_z = dot_product((U0(:, kx, ky, kz)), U0(:, kx, ky, kz+1))
							phase_z = (overlap_z)/abs(overlap_z)
							gauge_fixed_U_z = U0(:, kx, ky, kz+1) * conjg(phase_z)
							connection(3, kx, ky, kz) = -aimag((dot_product(((U0(:, kx, ky, kz))),((gauge_fixed_U_z-U0(:, kx, ky, kz))/dk))))!delkz)))
						enddo
					enddo
				enddo        
        do kx=1, 2*np-1
					do ky=1, 2*np-1
						do kz=1, 2*np-1
							curvature(1, kx, ky, kz) = (connection(3, kx, ky+1, kz)-connection(3, kx, ky, kz))/dk - &!delky
													   						 (connection(2, kx, ky, kz+1)-connection(2, kx, ky, kz))/dk    !delkz
							curvature(2, kx, ky, kz) = (connection(1, kx, ky, kz+1)-connection(1, kx, ky, kz))/dk - & !delkz
													   						 (connection(3, kx+1, ky, kz)-connection(3, kx, ky, kz))/dk   !delkx
							curvature(3, kx, ky, kz) = (connection(2, kx+1, ky, kz)-connection(2, kx, ky, kz))/dk - & !delkx
													   						 (connection(1, kx, ky+1, kz)-connection(1, kx, ky, kz))/dk   !delky
							magnitude_field(kx, ky, kz) = sqrt((curvature(1,kx,ky,kz))**2 + &
																					   		 (curvature(2,kx,ky,kz))**2 + &
																								 (curvature(3,kx,ky,kz))**2)
						enddo
					enddo
				enddo

        ! do kx=1, 2*np
				! 	do ky=1, 2*np
				! 		do kz=1, 2*np
				! 			connection(1, kx, ky, kz)=-aimag(dot_product(U0(:, kx, ky, kz),((U0(:, kx+1, ky, kz)-U0(:, kx, ky, kz))/dk)))
				! 			connection(2, kx, ky, kz)=-aimag(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky+1, kz)-U0(:, kx, ky, kz))/dk)))
				! 			connection(3, kx, ky, kz)=-aimag(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky, kz+1)-U0(:, kx, ky, kz))/dk)))
				! 		enddo ! Exit k-loop
				! 	enddo
				! enddo 
        !----Plaquette method
        ! do kx=1, 2*np
				! 	do ky=1, 2*np
				! 		do kz=1, 2*np
				! 			connection(1, kx, ky, kz)=(dot_product(U0(:, kx, ky, kz),((U0(:, kx+1, ky, kz)))/(abs(dot_product(U0(:, kx, ky, kz),((U0(:, kx+1, ky, kz))))))))
				! 			connection(2, kx, ky, kz)=(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky+1, kz)))/(abs(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky+1, kz))))))))
				! 			connection(3, kx, ky, kz)=(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky, kz+1)))/(abs(dot_product(U0(:, kx, ky, kz),((U0(:, kx, ky, kz+1))))))))
				! 		enddo ! Exit k-loop
				! 	enddo
				! enddo 
        do kx=1, 2*np
          do ky=1, 2*np
            do kz=1, 2*np
              ! X direction
              connection(1, kx, ky, kz) = dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx+1, ky, kz)) / &
                                          abs(dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx+1, ky, kz)))
              ! Y direction
              connection(2, kx, ky, kz) = dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx, ky+1, kz)) / &
                                          abs(dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx, ky+1, kz)))
              ! Z direction
              connection(3, kx, ky, kz) = dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx, ky, kz+1)) / &
                                          abs(dot_product(conjg(U0(:, kx, ky, kz)), U0(:, kx, ky, kz+1)))
            enddo
          enddo
        enddo
        !---

        ! do kx = 1, 2*np
        !   do ky = 1, 2*np
        !     do kz = 1, 2*np
        !       ! x-y plane
        !       curvature(1, kx, ky, kz) = aimag(log( connection(1, kx,   ky,   kz) * &
        !                                             connection(2, kx+1, ky,   kz) / &
        !                                             (connection(1, kx,   ky+1, kz) * &
        !                                              connection(2, kx,   ky,   kz)) ))
        
        !       ! y-z plane
        !       curvature(2, kx, ky, kz) = aimag(log( connection(2, kx,   ky,   kz) * &
        !                                             connection(3, kx,   ky+1, kz) / &
        !                                             (connection(2, kx,   ky,   kz+1) * &
        !                                              connection(3, kx,   ky,   kz)) ))
        
        !       ! z-x plane
        !       curvature(3, kx, ky, kz) = aimag(log( connection(3, kx,   ky,   kz) * &
        !                                             connection(1, kx,   ky,   kz+1) / &
        !                                             (connection(3, kx+1, ky,   kz) * &
        !                                              connection(1, kx,   ky,   kz)) ))
        
        !       magnitude_field(kx, ky, kz) = sqrt( real(curvature(1, kx, ky, kz))**2 + &
        !                                           real(curvature(2, kx, ky, kz))**2 + &
        !                                           real(curvature(3, kx, ky, kz))**2 )
        !     end do
        !   end do
        ! end do
        do kx=1, 2*np
					do ky=1, 2*np
						do kz=1, 2*np
							curvature(1, kx, ky, kz) = log(connection(2, kx, ky, kz)*connection(3, kx, ky+1, kz)/ &
													   						    (connection(2, kx, ky, kz+1)*connection(3, kx, ky, kz)))
							curvature(2, kx, ky, kz) = log(connection(3, kx, ky, kz)*connection(1, kx, ky, kz+1)/ &
                                            (connection(3, kx+1, ky, kz)*connection(1, kx, ky, kz)))
							curvature(3, kx, ky, kz) = log(connection(1, kx, ky, kz)*connection(2, kx+1, ky, kz)/ &
                                            (connection(1, kx, ky+1, kz)*connection(2, kx, ky, kz)))
							magnitude_field(kx, ky, kz) = sqrt((curvature(1,kx,ky,kz))**2 + &
																					   		 (curvature(2,kx,ky,kz))**2 + &
																								 (curvature(3,kx,ky,kz))**2)
						enddo
					enddo
				enddo
       gap_min(ipart)=abs(minval(ene(13,:))-maxval(ene(12,:))) !Unperturbed Case
    !    gapp_min(ipart)=abs(minval(enep(13,:))-maxval(enep(12,:))) !Perturbed Case
      !  ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0
       ef(ipart)=(minval(eneBerry(13,:,:,:))+maxval(eneBerry(12,:,:,:)))/2d0
!----- END FOURIER TRANSFORM

!------Band gap and Fermi energy calculations
    
   

    



!------Export data-------!

!------Only rank 0 writes output files
    if (rank == 0 .and. ipart >= local_start .and. ipart <= local_end) then

        write(partnumber,'(i5)') ipart
!-------Uncomment to write Energies to file
!-------These energies are needed for the fermi surface script
        !   do k=1,(np+1)**dim
        !         write(120, '(6(x,f12.6))') mesh(1:3,k), enep(13,k), ene(13,k), ef(ipart)
        !   end do


!-------Write Weyl points to file
!-------Uncoment depending on the Hamiltonian used, or if perturbation is added
            ! do k=1,(2*np+1)**dim
                !----4x4----!
                ! bandgap=ene(3,k)-ene(2,k) ! Unperturbed
                ! bandgapp=enep(3,k)-enep(2,k) ! Perturbed

                !-------18x18-------!
                ! bandgap=ene(13,k)-ene(12,k) ! Unperturbed
                ! bandgapp=enep(13,k)-enep(12,k) ! Perturbed

                !---------UNPERTURBED--------!
                ! if (abs(bandgap-(gap_min(ipart))) .lt. 0.001d0) then
                !     print*, abs(bandgap-(gap_min(ipart)))
                !     write(110, '(5(x,f12.6))') mesh(1:3,k), alpha, abs(bandgap-(gap_min(ipart)))!, gap_min(ipart)
                ! endif
                !---------PERTURBED--------!
                ! if (abs(bandgapp-(gapp_min(ipart))) .lt. 0.001d0) then
                !     print*, abs(bandgapp-(gapp_min(ipart)))
                !     write(110, '(5(x,f12.6))') mesh(1:3,k), alpha, abs(bandgapp-(gapp_min(ipart)))!, gap_min(ipart)
                ! endif
            ! end do

            ! write(100,*)
            ! write(100,*)
            ! write(110,*)
            ! write(110,*)
            
    ! endif

    ! if (rank == 0 .and. ipart >= local_start .and. ipart <= local_end .and. gap_perturbed(ipart) < 0.08) then
        write(partnumber,'(i5)') ipart
        write(line,'(3a)') 'curvatureweyl.dat'!'fermi_surface_energies_By_WSM_unfiltered.dat' 
        open(200,file=trim(line))
         
          ! do k=1,(2*np+1)**dim
          !       !  write(200, '(7(x,f12.6))') mesh(1:3,k), enep(12:13,k)!, ef(ipart)
          !        !Need the fermi energy of BCB again at the QCP, so need to plot only
          !        !a small selection of points at less than 0.05 eV around the node
          !       ! if (ene(13,k)-minval(ene(13,:)) .lt. 0.5d0) then
          !   ! if (abs(gap_min(ipart)) .lt. 0.05d0) then

          !       !  write(200, '(5(x,f12.6))') mesh(1:3,k), ene(12,k)-ef(ipart),ene(13,k)-ef(ipart)
					! 	write(200, '(5(x,f12.6))') mesh(1:3,k), ene(13,k), ef(ipart), curvature(1:3,)
          !       ! endif
          ! end do
          ! write(200,*)
          ! write(200,*)
          ! close(200)
    !  endif
            ! kz = np
						do kx=2, 2*np-1
							do ky=2, 2*np-1
								do kz=2, 2*np-1
						 			! k = (kx-1) + (ky-1)*(2*np+1) + (kz-1)*(2*np+1)**2 + 1
									write(200, '(13(x,f12.6))') meshBerry(1:3, kx, ky, kz), &
                                              real(curvature(1:3,kx,ky,kz)), &
                                              magnitude_field(kx, ky, kz)
								enddo
							enddo
						enddo
						write(200,*)
						write(200,*)
						close(200)
		endif
!------- Check time taken to calculate

       call cpu_time(part_time)
       part_time2 = part_time/60
       if (rank == 0) then
        print '(A, I3, A, F6.2)', "Partition ", ipart, " runtime (minutes): ", part_time2
       endif
    !    write(777, '(2(x,f12.6))') alpha, gap_min(ipart)
    end do ! Interpolation loop
    close(100)
    close(110)

    close(777)
    ! Gather results from all processes to rank 0
if (rank == 0) then
    ! We already have results for our local partitions
    ! Receive results from other processes
    do i = 1, nprocs-1
        local_start = i * local_count + 1
        local_end = (i + 1) * local_count
        if (i == nprocs - 1) then
            local_end = npartitions
        end if
        local_count = local_end - local_start + 1
        
        if (local_count > 0) then
            call MPI_RECV(gap_unperturbed(local_start), local_count, MPI_DOUBLE_PRECISION, i, 0, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            call MPI_RECV(ef(local_start), local_count, MPI_DOUBLE_PRECISION, i, 1, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        endif
    end do
    
    ! Write out the gap data for all partitions
    open(777, file='gap.dat')
    do ipart = 1, npartitions
        write(777, '(2(x,f12.6))') float(ipart-1)/float(npartitions-1), gap_unperturbed(ipart)
    end do
    close(777)
    
    ! Find the critical alpha with minimum gap
    temp_index = minloc(gap_unperturbed, dim=1)
    write(*,*) "Critical alpha with minimum gap: ", float(temp_index-1)/float(npartitions-1)
    write(*,*) "Minimum gap value: ", gap_unperturbed(temp_index)
    
else
    ! Send results to master process
    call MPI_SEND(gap(local_start), local_end-local_start+1, MPI_DOUBLE_PRECISION, &
                 0, 0, MPI_COMM_WORLD, ierr)
    call MPI_SEND(ef(local_start), local_end-local_start+1, MPI_DOUBLE_PRECISION, &
                 0, 1, MPI_COMM_WORLD, ierr)
endif

! Measure total runtime
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
mpi_end_time = MPI_WTIME()

if (rank == 0) then
    write(*,'(A, F10.4, A)') "Total MPI runtime: ", (mpi_end_time - mpi_start_time) / 60.0, " minutes"
endif

! Finalize MPI
call MPI_FINALIZE(ierr)

! Skip error handling in non-root processes
! if (rank /= 0) then
!     goto 999
! endif
! !-------Errors
! 333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),'" not found'
!     stop
! 444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),'" not found'
!     stop
! 445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),'" not found'
!     stop      

! 999 continue
!---------END
end program generate_fermi
