!==============================================================
! This program adiabatically interpolates the energy between two
! topological phases of BiTeI.
!==============================================================
program interpolate_topology
    use mpi
    use reading_module
    use fourier_module
    implicit none
    
    !--------Presets
    character(len=80):: prefix="BiTeI"
    integer, parameter :: np=30, npartitions=10, &  ! Adjust these parameters for resolution
                         dim=2                      ! dimensions in k-space
    real*8, parameter :: B_x = 0d0, B_y = 0.01d0, B_z = 0d0 ! Magnetic field components
    
    !---------MPI variables
    integer :: ierr, nprocs, rank, local_start, local_end, local_count
    real*8 :: mpi_start_time, mpi_end_time
    
    !---------Variable declarations
    character(len=80) :: hamil_file_trivial, hamil_file_topological, nnkp, line, partnumber
    character(len=200) :: hamil_dir = '/home/aleks/MPYS-PROJECT/Hamiltonians 18x18/'  ! 4x4 Hamiltonian directory
    
    integer :: ik, ipart, ib, is, i, j, k, n, nr, nb, i1, i2, lwork, info, &
               o, p, j1, j2, total_pairs, temp_index
    
    real*8 :: phase, dx, dy, dz, twopi, jk, a, b, a1, b1, &
              spin_x(1,1), spin_y(1,1), spin_z(1,1), &
              spin_xp(1,1), spin_yp(1,1), spin_zp(1,1), &
              alpha, ef(npartitions), gap(npartitions),efp(npartitions), &
              gapp(npartitions), &
              write_values(11:14), &
              bvec(3,3), avec(3,3), rvecs(3), &
              ktemp1(3), ktemp2(3), &
              kxtemp1(3), kxtemp2(3), &
              kytemp1(3), kytemp2(3), &
              kmesh(np,np), &
              mesh_kx(3,np,np), mesh_ky(3,np,np), &
              mesh_gap(3,np**2), &
              part_time, part_time2
    
    complex*16 :: sigx(2,2), sigy(2,2), sigz(2,2), &
                 chi(2,1), chip(2,1), &
                 phi(3), phase_factor
    
    real*8, dimension(npartitions) :: min_eigenvalue, alpha_values
    
    ! Allocatable arrays
    integer, allocatable :: ndeg(:),ndeg_trivial(:),ndeg_topological(:)
    real*8, allocatable :: mesh(:,:), phases(:,:), rvec_trivial(:,:), rvec_topological(:,:), &
                           ene(:,:), enep(:,:), rwork(:), rvec(:,:), &
                           spin(:,:,:), spinp(:,:,:)
    complex*16, allocatable :: H(:,:), Hk(:,:), Hm(:,:), Hmag(:,:), &
                              HK_trivial(:,:), HK_topological(:,:), &
                              Hamr_trivial(:,:,:), Hamr_topological(:,:,:), &
                              work(:)
    
    real*8, parameter :: x_min = -0.1d0, x_max = 0.1d0, &
                         y_min = -0.1d0, y_max = 0.1d0
    
    integer, allocatable :: indices(:)
    
    ! Flags
    logical :: useOptimized = .true.  ! Set to false for 4x4 case

!--------- Initialize MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
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
!--------- Define Pauli Matrices sigma x,y,z
    data sigx /(0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
    data sigy /(0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
    data sigz /(1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/

! Only process 0 reads the input files
    if (rank == 0) then
!---------------  Read the vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)', iostat=info)line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 111
      read(98,*)avec
      read(98,'(a)')line
      read(98,'(a)')line
      read(98,'(a)')line
      read(98,*)bvec

      call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!------read trivial and topological H(R)
    ! open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
    ! open(97,file=trim(adjustl(hamil_file_topological)),err=445)
    ! read(99,*)
    ! read(97,*)
    ! read(97,*)!Skip nb,nr same as other file
    ! read(99,*)nb,nr
      call read_header(hamil_file_trivial, nb, nr)
      call read_header(hamil_file_topological, nb, nr)
    endif
    ! Broadcast necessary values to all processes

  call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        allocate(rvec(3, nr),rvec_trivial(3, nr), rvec_topological(3, nr))
        allocate(Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr))
        allocate(Hk_topological(nb,nb),Hk_trivial(nb,nb))
        allocate(H(nb,nb),Hk(nb,nb))
        allocate(ndeg(nr),ndeg_trivial(nr),ndeg_topological(nr))
        allocate(enep(nb, (np+1)**dim), ene(nb, (np+1)**dim))
  ! Only rank 0 reads the Hamiltonian data
    if (rank == 0) then

    if (useOptimized) then
        if (rank == 0) then
            call read_optimized_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr, &
            Hamr_trivial, Hamr_topological, rvec, ndeg, avec)
            rvec_topological = rvec_trivial
        end if
        else
            call read_general_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr, &
            Hamr_trivial, Hamr_topological, rvec_trivial, rvec_topological, ndeg_trivial, ndeg_topological, avec)
        endif
    end if


! Broadcast the read data to all processes
    call MPI_BCAST(ndeg_trivial, nr, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ndeg_topological, nr, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rvec_trivial, 3*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rvec_topological, 3*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(hamr_trivial, nb*nb*nr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(hamr_topological, nb*nb*nr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
!------ LAPACK-related array allocations
      lwork=max(1,2*nb-1)
      allocate(work(max(1,lwork)),&
               rwork(max(1,3*nb-2)))
      !allocate(weightx(nb,(np+1)**2),weighty(nb,(np+1)**2),weightz(nb,(np+1)**2))
      allocate(spin(3,nb,(np+1)**dim),spinp(3,nb,(np+1)**dim))
!------ open gap file
      open(777,file='gap.dat')

!-------- Calculate the total number of pairs
      allocate(mesh(3, (np+1)**dim))
      !----- Loop through all x values and all y values
      !----- Store the pairs in the result array
      !-------Generate Mesh
      
          dx = (x_max - x_min)/np
          dy = (y_max - y_min)/np
         ! dz = (z_max - z_min)/np
          
          j=0
          do i = 0, np
              do n = 0, np
                      ! Calculate coordinates directly
                      j = j + 1
                      mesh(1, j) = (x_min + i * dx )* bvec(1,1)
                      mesh(2, j) = (y_min + n * dy )* (bvec(1,2)+bvec(2,2))
                      mesh(3, j) = 0.5d0 *bvec(3,3)! Exactly on the BZ boundary
                  
              end do
          end do

!------ Magnetic Field
      allocate(Hm(2,2),Hmag(nb,nb))
      Hm = B_x*sigx + B_y*sigy + B_z*sigz
!------ Turn Hm into an 18x18 to match Hk      
      Hmag = (0d0, 0d0)
      do i=1, nb/2
            Hmag(i,i)=Hm(1,1)
            Hmag(i,i+nb/2)=Hm(1,2)
            Hmag(i+nb/2,i)=Hm(2,1)
            Hmag(i+nb/2,i+nb/2)=Hm(2,2)
      end do
!------ Fourrier transform H(R) to H(k)
    !   allocate(phases(nr, (np+1)**dim))

      local_count = npartitions / nprocs
      local_start = rank * local_count + 1
      local_end = (rank + 1) * local_count
      
      if (rank == nprocs - 1) then
          local_end = npartitions  ! Last process takes any remainder
      end if
      
      ! Each process handles its own partitions
      do ipart = local_start, local_end
          if (rank == 0) then
              write(*,'(a,i5,a,i5)') 'Partition=', ipart, ' of ', npartitions
          endif
         alpha=float(ipart-1)/float(npartitions-1)
         !alpha=0d0
         ! Initialize Hamiltonians for the current partition
         ene=0d0
! !----- FOURIER TRANSFORM 
!          call fourier_transform_general(np, dim, nr, nb, ndeg_trivial, ndeg_topological, mesh, bvec, avec, &
!           rvec_trivial, rvec_topological, Hamr_trivial, Hamr_topological, &
!           Hmag, alpha, enep, ene, work, lwork, rwork, rank, ierr)

         call fourier_transform_optimized(np, dim, nr, nb, ndeg, mesh, rvec, &
         Hamr_trivial, Hamr_topological, Hmag, alpha, &
         enep, ene, work, lwork, rwork, rank, ierr)

!----- END FOURIER TRANSFORM
            
!----------Find spins by applying pauli             
           ! do ib=1,nb
           !       do is=1,nb/2
           !             chi(1,1) = Hk(is     ,ib) 
           !             chi(2,1) = Hk(is+nb/2,ib)
           !             chip(1,1) = H(is     ,ib) 
           !             chip(2,1) = H(is+nb/2,ib)
                        
                        !spin_x = matmul(conjg(transpose(chi)),matmul(sigx, chi))
                        !spin_y = matmul(conjg(transpose(chi)),matmul(sigy, chi))
                        !spin_z = matmul(conjg(transpose(chi)),matmul(sigz, chi))
                        !spin(1,ib,k)=spin(1,ib,k)+spin_x(1,1)
                        !spin(2,ib,k)=spin(2,ib,k)+spin_y(1,1)
                        !spin(3,ib,k)=spin(3,ib,k)+spin_z(1,1)
                        !Calculate spins for Perturbed Hamiltonian
                        !spin_xp = matmul(conjg(transpose(chip)),matmul(sigx, chip))
                        !spin_yp = matmul(conjg(transpose(chip)),matmul(sigy, chip))
                        !spin_zp = matmul(conjg(transpose(chip)),matmul(sigz, chip))
      
                        !spinp(1,ib,k)=spinp(1,ib,k)+spin_xp(1,1)
                        !spinp(2,ib,k)=spinp(2,ib,k)+spin_yp(1,1)
                        !spinp(3,ib,k)=spinp(3,ib,k)+spin_zp(1,1)
            !      enddo
            !enddo

         
!------calcualte gap and Fermi level
         gapp(ipart)= minval(enep(3,:))-maxval(enep(2,:))
         gap(ipart)= minval(ene(3,:))-maxval(ene(2,:))
         ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0
        !  efp(ipart)=(minval(enep(13,:))+maxval(enep(12,:)))/2d0

!------Export data
           ! Only rank 0 writes output files
      if (rank == 0 .and. ipart >= local_start .and. ipart <= local_end) then
          write(partnumber,'(i5)') ipart
          write(line,'(3a)') 'k_surface_fermi_energies_By_0.05_part_',trim(adjustl(partnumber)),'.dat' 
          open(100,file=trim(line))
             
            do k=1,(np+1)**dim
                  write(100,'(6(x,f12.6))') mesh(1:2,k), (ene(i,k), i=11,14)!,&
                                            ! spinp(1:3,i,k),&!need to minimize the energy wrt fermi energy
                                            ! sqrt(spinp(1,i,k)**2 +spinp(2,i,k)**2 +spinp(3,i,k)**2)) !This now writes into the files the coordinates as a function of the TCB and BCB energy difference
                  !write(200,'(3(x,f12.6))') mesh(1:2,k),ene(i,k)
            end do
            write(100,*)
            write(100,*)
            close(100)
      endif
       
!------- Check time taken to calculate

      call cpu_time(part_time)
      part_time2 = part_time/60
      if (rank == 0) then 
       print '(A, I3, A, F6.2)', "Partition ", ipart, " runtime (minutes): ", part_time2
      endif
   end do

   ! Gather results from all processes to rank 0
if (rank == 0) then
   ! We already have results for our local partitions
   ! Receive results from other processes
   do i = 1, nprocs-1
       local_start = i * local_count + 1
       local_end = (i + 1) * local_count
       if (i == nprocs - 1) then
           local_end = npartitions
       endif
       local_count = local_end - local_start + 1
       
       if (local_count > 0) then
           call MPI_RECV(gap(local_start), local_count, MPI_DOUBLE_PRECISION, i, 0, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
           call MPI_RECV(ef(local_start), local_count, MPI_DOUBLE_PRECISION, i, 1, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       endif
   end do
   
   ! Write out the gap data for all partitions
   open(777, file='gap.dat')
   do ipart = 1, npartitions
       write(777, '(3(x,f12.6))') float(ipart-1)/float(npartitions-1), gap(ipart), gapp(ipart)
   end do
   close(777)
   
   ! Find the critical alpha with minimum gap
   temp_index = minloc(gap, dim=1)
   write(*,*) "Critical alpha with minimum gap: ", float(temp_index-1)/float(npartitions-1)
   write(*,*) "Minimum gap value: ", gap(temp_index), gapp(temp_index)
   
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
if (rank /= 0) then
   goto 999
endif
!-------Errors
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),'" not found'
   stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),'" not found'
   stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),'" not found'
   stop      

999 continue 
!---------END
end program interpolate_topology
