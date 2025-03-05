!=============================================================
! This program perturbs the energy with a magnetic field
! and calculates the fermi surface.
!===============================================================
Program generate_fermi
    use mpi
   ! #include <mpif.h>
    Implicit None
!--------Presets
    character(len=80):: prefix="BiTeI"
    integer,parameter::np=35,npartitions=40!Adjust these parameters to obtain better resolution around alphacrit and see points closer to an effectively closed gap
    
    real*8,parameter::B_x = 0d0, B_y = 0.01d0, B_z = 0d0,& !Run again at B_y=0.05-0.06 to see the gap close
                      alpha_min = 0.76d0, alpha_max = 0.79d0

!---------MPI variables
    integer :: ierr, nprocs, rank, local_start, local_end, local_count
    real*8 :: mpi_start_time, mpi_end_time

!---------Variable allocation

    character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber

    integer ik, ipart, ib, is,&
            i,j,k,&
            n,nr,nb,&
            i1,i2,&
            lwork,info&
            ,o,p,j1,j2,&
            total_pairs,&
            temp_index
    
    real*8 phase, dx, dy, dz, bandgap,&
           twopi,jk,a,b,a1,b1,&
           spin_x(1,1),spin_y(1,1),spin_z(1,1),&
           spin_xp(1,1),spin_yp(1,1),spin_zp(1,1),&
           alpha,ef(npartitions),efp(npartitions),&
           gap(npartitions),gap_perturbed(npartitions),gap_unperturbed(npartitions),gap_min(npartitions),&
           write_values(11:14),&
           bvec(3,3),avec(3,3),rvecs(3),&
           ktemp1(3),ktemp2(3),&
           kxtemp1(3),kxtemp2(3),&
           kytemp1(3),kytemp2(3),&
           kmesh(np,np),&
           mesh_kx(3,np, np), mesh_ky(3,np, np),&
           mesh_gap(3, np**2),&
           part_time,part_time2


    complex*16 sigx(2, 2), sigy(2, 2), sigz(2, 2),&
               chi(2,1),chip(2,1),&
               phi(3),phase_factor

    real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
    
    integer,allocatable:: ndeg(:),ndeg_topological(:)

    real*8,allocatable:: phases(:,:),rvec_trivial(:,:),rvec_topological(:,:),&
                         ene(:,:),enep(:,:),&
                         rwork(:),rvec(:,:),&
                         spin(:,:,:),spinp(:,:,:)
                         
    complex*16,allocatable:: H(:,:), Hk(:,:), Hm(:,:), Hmag(:,:),&
                             HK_trivial(:,:), HK_topological(:,:),&
                             Hamr_trivial(:,:,:), Hamr_topological(:,:,:),&
                             work(:)

    real*8, parameter :: x_min = -0.07d0, x_max = 0.07d0,&
                         y_min = -0.07d0, y_max = 0.07d0,&
                         z_min = -0.03d0, z_max = 0.03d0

    integer, dimension(:), allocatable:: indices

    real*8, dimension(:,:), allocatable :: mesh
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


!------------------------------------------------------
  write(hamil_file_trivial,'(a,a)')trim(adjustl(prefix)), "_hr_trivial.dat"!Why were these (a,a,a)?
  write(hamil_file_topological,'(a,a)')trim(adjustl(prefix)), "_hr_topological.dat"
  write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
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
111   read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,*)bvec

!------read trivial and topological H(R)
  open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
  open(97,file=trim(adjustl(hamil_file_topological)),err=445)
  read(99,*)
  read(99,*)nb,nr
  endif
!------Broadcast necessary values to all processes
call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  allocate(rvec(3,nr))
  allocate(Hk(nb,nb),Hamr_trivial(nb,nb,nr),&
           Hamr_topological(nb,nb,nr),Hk_topological(nb,nb),&
           H(nb,nb), Hk_trivial(nb,nb),&
           ndeg(nr),ene(nb,(np+1)**3),enep(nb,(np+1)**3))

! Only rank 0 reads the Hamiltonian data
           if (rank == 0) then
  read(99,*)ndeg
  do i = 1, 80
    read(97,*)
  end do
  do k=1,nr
     do i=1,nb
        do j=1,nb
           read(99,*)rvecs(1),rvecs(2),rvecs(3),i1,i2,a,b
           hamr_trivial(i1,i2,k)=dcmplx(a,b)
           read(97,*)rvecs(1),rvecs(2),rvecs(3),i1,i2,a1,b1
           hamr_topological(i1,i2,k)=dcmplx(a1,b1)
        end do
     end do
     rvec(:,k) = rvecs(1)*avec(:,1) + rvecs(2)*avec(:,2) + rvecs(3)*avec(:,3)
  end do
endif

! Broadcast the read data to all processes
call MPI_BCAST(ndeg, nr, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(rvec, 3*nr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(hamr_trivial, nb*nb*nr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(hamr_topological, nb*nb*nr, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
!------ LAPACK-related array allocations
lwork=max(1,2*nb-1)
allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
allocate(spin(3,nb,(np+1)**2),spinp(3,nb,(np+1)**2))

    allocate(mesh(3, (np+1)**3))
!------ open gap file
    open(777,file='gap.dat',status='replace', position='append', action='write')
    open(100,file='k_points_for_alpha_0.01.dat',status='replace', position='append', action='write')

!-------Generate Mesh

    dx = (x_max - x_min)/np
    dy = (y_max - y_min)/np
    dz = (z_max - z_min)/np
    
    j=0
    do i = 0, np
        do n = 0, np
            do k = 0, np
                j = j+1  ! Direct index calculation
                ! Calculate coordinates directly
                mesh(1, j) = (x_min + i * dx )!* bvec(1,1)
                mesh(2, j) = (y_min + n * dy )!* (bvec(1,2)+bvec(2,2))
                mesh(3, j) = (z_min + k * dz ) + 0.5*bvec(3,3)
            end do
        end do
    end do
!------ Magnetic Field
!     allocate(Hm(2,2),Hmag(nb,nb))
!     Hm = B_x*sigx + B_y*sigy + B_z*sigz
! !------ Turn Hm into an 18x18 to match Hk      
!     Hmag = (0d0, 0d0)
!     do i=1, nb/2
!           Hmag(i,i)=Hm(1,1)
!           Hmag(i,i+nb/2)=Hm(1,2)
!           Hmag(i+nb/2,i)=Hm(2,1)
!           Hmag(i+nb/2,i+nb/2)=Hm(2,2)
!     end do
!------ Fourrier transform H(R) to H(k)
    allocate(phases(nr, (np+1)**3))

    ! local_count = npartitions / nprocs
    ! local_start = rank * local_count + 1
    ! local_end = (rank + 1) * local_count
    
    if (rank == nprocs - 1) then
        local_end = npartitions  ! Last process takes any remainder
    end if
    
    ! Each process handles its own partitions
    ! critical_alpha=0.789473712
    
    do ipart = 1, npartitions
        if (rank == 0) then
            write(*,'(a,i5,a,i5)') 'Partition=', ipart, ' of ', npartitions
        endif
       !alpha=float(ipart-1)/float(npartitions-1)
        alpha = alpha_min + float(ipart - 1)*(alpha_max - alpha_min)/float(npartitions - 1)
        if (rank == 0) then
            write(*,'(A,I5,A,F12.6)') 'Partition ', ipart, ' alpha = ', alpha
        endif

        ! Initialize Hamiltonians for the current partition
       ene=0d0

       do k=1,(np+1)**3
        do j=1,nr
          phases(j,k)=dot_product(mesh(:,k),rvec(:,j))
        end do
       end do
          ! Fourier transform terms
           !I AM CURRENTLY OPTIMIZING STUFF
             ! Compute phase factors
       do k=1,(np+1)**3
          HK_trivial=(0d0,0d0)
          HK_topological=(0d0,0d0)

          do j=1,nr
            phase = phases(j,k)
            phase_factor = dcmplx(cos(phase), -sin(phase))/float(ndeg(j))
             ! Sum over H(R) contributions for each Hamiltonian
            Hk_trivial = Hk_trivial + Hamr_trivial(:,:,j) * phase_factor
            Hk_topological = Hk_topological + Hamr_topological(:,:,j) * phase_factor
                
          end do
!----------Interpolate between the trivial and topological states 
          Hk=Hk_trivial*(1-alpha)+Hk_topological*alpha
!----------Perturb Hamiltonian
        !   H = Hk+Hmag
       !enddo
!----------Compute eigenvalues and eigenvectors
         ! call zheev('V','U',nb,H,nb,enep(:,k),work,lwork,rwork,info)
          call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
          if (info /= 0) then
            if (rank == 0) write(*,*) "ZHEEV failed with info =", info
            call MPI_ABORT(MPI_COMM_WORLD, info, ierr)
        endif
    end do      !Close k-loop 
!------calcualte gap and Fermi level
       gap_min(ipart)=minval(ene(13,:))-maxval(ene(12,:))
       print *, "MINGAP: ", gap_min(ipart)
       ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0
    !    efp(ipart)=(minval(enep(13,:))+maxval(enep(12,:)))/2d0
    !    gap_perturbed(ipart)=minval(enep(13,:))-maxval(enep(12,:))
    !    gap()=abs(ene(13,i)-ene(12,i))
!------Export data
!------Only rank 0 writes output files
    if (rank == 0 .and. ipart >= local_start .and. ipart <= local_end) then
!------Export alphas with minimum gap value
        write(partnumber,'(i5)') ipart
        ! open(101, file="gap_vs_alpha.dat")
        
        ! close(101)
!------Write k-points for small gap
    !    write(line,'(3a)') 'k_surface_fermi_energies_By_WSM_trajectory.dat' 
    !    write(line,'(3a)') 'k_points_for_alpha.dat'
        !   do k=1,(np+1)**3
        !         write(100, '(5(x,f12.6))') mesh(1:3,k), ene(12:13,k)- ef(ipart)
        !   end do
            do k=1,(np+1)**3
                bandgap=ene(13,k)-ene(12,k)
                ! print*, bandgap
                if (abs(bandgap-gap_min(ipart)) .lt. 0.01d0) then
                    ! print*, bandgap
                    write(100, '(5(x,f12.6))') mesh(1:3,k), alpha, abs(bandgap-gap_min(ipart))!, gap_min(ipart)
                endif
            end do

            write(100,*)
            write(100,*)
            
    endif
    ! if (rank == 0 .and. ipart >= local_start .and. ipart <= local_end .and. gap_perturbed(ipart) < 0.08) then
    !     write(partnumber,'(i5)') ipart
    !     write(line,'(3a)') 'k_surface_fermi_energies_By_WSM_perturbed.dat' 
    !     open(200,file=trim(line))
         
    !        do k=1,(np+1)**3
    !              write(200, '(7(x,f12.6))') mesh(1:3,k), enep(12:13,k)!, ef(ipart)
    !        end do
    !          write(200,*)
    !          write(200,*)
    !          close(200)
    !  endif
     
!------- Check time taken to calculate

       call cpu_time(part_time)
       part_time2 = part_time/60
       if (rank == 0) then
        print '(A, I3, A, F6.2)', "Partition ", ipart, " runtime (minutes): ", part_time2
       endif
       write(777, '(2(x,f12.6))') alpha, gap_min(ipart)
    end do
    close(100)
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
end program generate_fermi
