!==============================================================
! This program adibatically interpolates the energy between two
! topological phases of BiTeI.
!==============================================================
Program interpolate_topology
      Implicit None
!--------Presets
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=200,npartitions=1!Adjust these parameters to obtain better resolution around alphacrit and see points closer to an effectively closed gap
      
      real*8,parameter::B_x = 0d0, B_y = 0.05d0, B_z = 0d0
!---------Variable allocation
      character(len=30) :: klabel(nkpath),kxlabel(nkpath),kylabel(nkpath)
      character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber

      integer ik, ipart, ib, is,&
              i,j,k,&
              n,nr,nb,&
              i1,i2,&
              lwork,info&
              ,o,p,j1,j2,&
              total_pairs

      real*8 phase_trivial,phase_topological,&
             twopi,jk,a,b,a1,b1,&
             spin_x(1,1),spin_y(1,1),spin_z(1,1),&
             spin_xp(1,1),spin_yp(1,1),spin_zp(1,1),&
             alpha,ef(npartitions),gap(npartitions),&
             write_values(13:14),&
             bvec(3,3),avec(3,3),&
             ktemp1(3),ktemp2(3),&
             xkl(nkpath),ykl(nkpath),&
             kxtemp1(3),kxtemp2(3),&
             kytemp1(3),kytemp2(3),&
             kmesh(np,np),&
             mesh_kx(3,np, np), mesh_ky(3,np, np),&
             mesh_gap(3, np**2),&
             part_time,part_time2

      complex*16 sigx(2, 2), sigy(2, 2), sigz(2, 2),&
                 chi(2,1),chip(2,1),&
                 phi(3)

      real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
      !real*8,parameter::third=1d0/3d0
      integer,allocatable:: ndeg_trivial(:),ndeg_topological(:)
      real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),&
                           ene(:,:),enep(:,:),&
                           rwork(:),rvec(:,:),&
                           spin(:,:,:),spinp(:,:,:)
                           
      complex*16,allocatable:: H(:,:), Hk(:,:), Hm(:,:), Hmag(:,:),&
                               HK_trivial(:,:), HK_topological(:,:),&
                               Hamr_trivial(:,:,:), Hamr_topological(:,:,:),&
                               work(:)
      real*8, parameter :: x_min = -0.2, x_max = 0.2, y_min = -0.2, y_max = 0.2
      real*8, dimension(:), allocatable :: x_vals, y_vals
      real*8, dimension(:,:), allocatable :: mesh
      !complex*16,dimension(2,2) :: sigx(2,2),sigy(2,2),sigz(2,2)
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

!---------------  Read the vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 111
      read(98,*)avec
      read(98,'(a)')line
      read(98,'(a)')line
      read(98,'(a)')line
      read(98,*)bvec
!------- Generate x and y values
    allocate(x_vals(np+1), y_vals(np+1))
    do i=0,np
    	x_vals(i+1) = x_min + i * (x_max - x_min)/np
    	y_vals(i+1) = y_min + i * (y_max - y_min)/np
    enddo
!-------- Calculate the total number of pairs
    allocate(mesh(3, (np+1)**2))
!----- Loop through all x values and all y values
!----- Store the pairs in the result array
    j = 1
    do i = 1, np+1
        do n = 1, np+1
            mesh(1, j) = x_vals(i)*bvec(1,1)
            mesh(2, j) = y_vals(n)*(bvec(1,2)+bvec(2,2))
            mesh(3, j) = 0.5d0 *bvec(3,3) ! Exactly on the BZ boundary
            j = j + 1
        end do
    end do
!------read trivial H(R)
    open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec(3,nr))
    allocate(rvec_trivial(3,nr),Hk_trivial(nb,nb),Hamr_trivial(nb,nb,nr),&
             ndeg_trivial(nr),ene(nb,(np+1)**2),enep(nb,(np+1)**2))
    read(99,*)ndeg_trivial
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_trivial(1,k),rvec_trivial(2,k),rvec_trivial(3,k),i1,i2,a,b
             hamr_trivial(i1,i2,k)=dcmplx(a,b)
          enddo
       enddo
    enddo
    close(99)

!------read topolgical H(R)
    open(99,file=trim(adjustl(hamil_file_topological)),err=445)
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_topological(3,nr),Hk_topological(nb,nb),Hamr_topological(nb,nb,nr), &
           ndeg_topological(nr),HK(nb,nb),H(nb,nb))
    read(99,*)ndeg_topological
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_topological(1,k),rvec_topological(2,k),rvec_topological(3,k),i1,i2,a,b
             hamr_topological(i1,i2,k)=dcmplx(a,b)
          enddo
       enddo
       rvec(:,k) = rvec_topological(1,k)*avec(:,1) + rvec_topological(2,k)*avec(:,2) + rvec_topological(3,k)*avec(:,3)
    enddo
    close(99)
!------ LAPACK-related array allocations
      lwork=max(1,2*nb-1)
      allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
      !allocate(weightx(nb,(np+1)**2),weighty(nb,(np+1)**2),weightz(nb,(np+1)**2))
      allocate(spin(3,nb,(np+1)**2),spinp(3,nb,(np+1)**2))
!------ open gap file
      open(777,file='gap.dat')
!------ Magnetic Field
      allocate(Hm(2,2),Hmag(18,18))
      Hm = B_x*sigx + B_y*sigy + B_z*sigz
      
      do i=1, nb/2
            Hmag(i,i)=Hm(1,1)
            Hmag(i,i+nb/2)=Hm(1,2)
            Hmag(i+nb/2,i)=Hm(2,1)
            Hmag(i+nb/2,i+nb/2)=Hm(2,2)
      enddo
!------ Fourrier transform H(R) to H(k)
      ene=0d0
      do ipart=1,npartitions
         write(*,'(a,i5)') 'Partition=',ipart
         !alpha=float(ipart-1)/float(npartitions-1)
         alpha=0d0
         ! Initialize Hamiltonians for the current partition
         do k=1,(np+1)**2
                HK_trivial=(0d0,0d0)
            HK_topological=(0d0,0d0)

            ! Fourier transform terms
            do j=1,nr
   
                   phase_trivial=0.0d0
               phase_topological=0.0d0
   
               ! Compute phase factors
               do i=1,3
                      phase_trivial=phase_trivial    +mesh(i,k)*rvec(i,j)
                  phase_topological=phase_topological+mesh(i,k)*rvec(i,j)
                  
               enddo

               ! Sum over H(R) contributions for each Hamiltonian
               do i1=1,nb
                  do i2=1,nb
                      Hk_trivial(i1,i2)=Hk_trivial(i1,i2)+Hamr_trivial(i1,i2,j)* &
                                        dcmplx(cos(phase_trivial), &
                                              -sin(phase_trivial))/float(ndeg_trivial(j))

                     Hk_topological(i1,i2)=Hk_topological(i1,i2)+Hamr_topological(i1,i2,j)* &
                                          dcmplx(cos(phase_topological), &
                                                -sin(phase_topological))/float(ndeg_topological(j))
                  enddo
               enddo
            enddo
!----------Interpolate between the trivial and topological states 
            Hk=alpha*HK_trivial+(1d0-alpha)*HK_topological
!----------Perturb Hamiltonian
            H = HK+Hmag
         !enddo
!----------Compute eigenvalues and eigenvectors
            call zheev('V','U',nb,H,nb,enep(:,k),work,lwork,rwork,info)
            call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
      enddo      !Close k-loop
            
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
         gap(ipart)= minval(ene(13,:))-maxval(ene(12,:))
         ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0
!------Export data
         write(partnumber,'(i5)') ipart
         write(line,'(3a)') 'Energy_part_0.2B',trim(adjustl(partnumber)),'.dat' 
         open(100,file=trim(line))
         open(200,file='energy.dat')

            do k=1,(np+1)**2
                 ! write_values(13:14) = 0.0  ! Assuming i ranges from 11 to 14
                 ! do i = 13,14
                      ! Check if ene(i,k) - ef(ipart) is less than 0.01 and set it to 0 if true
                     ! if (abs(ene(i,k)-minval(ene(i,:))) .lt. 0.001) then
                     !     write_values(i) = 0.0
                     ! else
                     !     write_values(i) = ene(i,k)-minval(ene(i,:))
                     ! end if
                 ! enddo
                 ! write(*,*) 'Debug: ene(13,k) - ef(ipart) for k = ', k, ' = ', ene(13,k) - ef(ipart)
                  write(100,'(6(x,f12.6))') mesh(1:2,k), (enep(i,k)-ef(ipart), i=11,14)!,&
                                            ! spinp(1:3,i,k),&!need to minimize the energy wrt fermi energy
                                            ! sqrt(spinp(1,i,k)**2 +spinp(2,i,k)**2 +spinp(3,i,k)**2)) !This now writes into the files the coordinates as a function of the TCB and BCB energy difference
                  !write(200,'(3(x,f12.6))') mesh(1:2,k),ene(i,k)
            enddo
              write(100,*)
              write(100,*)
             ! write(200,*)
             ! write(200,*)
      
        
         close(100)
       
!------- Export Gap and Fermi energy
        ! write(777,'(3(x,f12.6))') alpha,gap(ipart),ef(ipart) !this had ipart in front of alpha
         call cpu_time(part_time)
         part_time2 = part_time/60
         print'(A, F6.2)', "Total runtime (minutes): ", part_time2
      enddo
!-------Errors
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),'" not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),'" not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),'" not found'
      stop      
!---------END
end program interpolate_topology
