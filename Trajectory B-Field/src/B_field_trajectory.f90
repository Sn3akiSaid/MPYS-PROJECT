!==============================================================
! This program adibatically interpolates the energy between two
! topological phases of BiTeI.
!==============================================================
Program trajectory_xyz
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=100,npartitions=20,nx=50,ny=50,nz=40
      real,parameter::B_x=0, B_y=0.1, B_z=0
!-----------------------------------------------------
      integer*4 ik,ipart
      real*8 alpha,ef(npartitions),gap(npartitions)
      character(len=30)::klabel(nkpath),kxlabel(nkpath),kylabel(nkpath)
      character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber
      integer*4 i,j,k,n,nr,i1,i2,nb,lwork,info,temp_index
      real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
      real*8,parameter::third=1d0/3d0
      real*8 phase_trivial,phase_topological,twopi,jk,a,b,thirdpi
      real*8 avec(3,3),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),ykl(nkpath)
      real*8 kxtemp1(3),kxtemp2(3),kytemp1(3),kytemp2(3),kmesh(np,np),min_k(3,npartitions)
      real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),ene(:,:),rwork(:),enep(:,:)
      real(8) :: mesh_kx(3,np, np), mesh_ky(3,np, np), mesh_gap(3, np**2),x2,y2,xtemp,ytemp
      integer*4,allocatable:: ndeg_trivial(:),ndeg_topological(:)
      complex*16,allocatable:: Hk(:,:),HK_trivial(:,:),HK_topological(:,:),Hamr_trivial(:,:,:)
      complex*16,allocatable:: Hamr_topological(:,:,:),work(:),Hm(:,:),H(:,:),Hmag(:,:)
      real(8), parameter :: x_min =0.016,x_max =0.051,y_min =-0.019,y_max = -0.053,z_min=0.498,z_max=0.502
      complex*8 sigx(2, 2), sigy(2, 2), sigz(2, 2)
    real(8), dimension(:), allocatable :: x_vals, y_vals, z_vals
    real(8), dimension(:,:), allocatable :: mesh,mesh2
    real(8), dimension(:,:,:), allocatable :: dpoint_traj
    real :: start_time, part_time, start_time2, part_time2
    integer, dimension(:), allocatable :: indices
    integer :: total_pairs
!------------------------------------------------------
      call cpu_time(start_time)
      write(hamil_file_trivial,'(a,a,a)') "hamiltonians/", trim(adjustl(prefix)), "_hr_trivial_4.dat"
      write(hamil_file_topological,'(a,a,a)') "hamiltonians/", trim(adjustl(prefix)), "_hr_topological_4.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
	
      thirdpi=4.0d0*atan(1.0d0)*third
      twopi=4.0d0*atan(1.0d0)*2.0d0
      data sigx /(0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
      data sigy /(0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
      data sigz /(1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/
      allocate(Hm(2,2))
      Hm = B_x*sigx + B_y*sigy + B_z*sigz
!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
222   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 222
      read(98,*)avec
      
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      read(98,*)bvec
	
	
	
	open(2000,file='trajectories.dat', status = 'replace', action = 'write')
	close(2000)
!---------------kpath
    ! Generate x and y values
    allocate(x_vals(nx+1), y_vals(ny+1), z_vals(nz+1))
    do i =0,nx
    	x_vals(i+1) = x_min + i * (x_max - x_min)/nx
    	y_vals(i+1) = y_min + i * (y_max - y_min)/ny
    	
enddo
do i = 0,nz
z_vals(i+1) = z_min + i * (z_max - z_min)/nz
enddo

    ! Calculate the total number of pairs
    allocate(mesh2(3, (nx+1)*(ny+1)*(nz+1)))

    ! Loop through all x values and all y values
    ! Store the pairs in the result array
    j = 1
    do i = 1, nx+1
        do n = 1,ny+1
        do k = 1,nz+1
            mesh2(1, j) = x_vals(i)
            mesh2(2, j) = y_vals(n)
            mesh2(3, j) = z_vals(k)
            j = j + 1
        end do
    end do
    enddo

      

!------read trivial H(R)
      open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec_trivial(3,nr),Hk_trivial(nb,nb),Hamr_trivial(nb,nb,nr),&
               ndeg_trivial(nr),ene(nb,(nx+1)*(ny+1)*(nz+1)))
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
	       ndeg_topological(nr),HK(nb,nb))
      read(99,*)ndeg_topological
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec_topological(1,k),rvec_topological(2,k),rvec_topological(3,k),i1,i2,a,b
               hamr_topological(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(99)
!------ LAPACK-related array allocations
      lwork=max(1,2*nb-1)
      allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
!------ open gap file
      open(777,file='gap.dat')


! Initialize the array to store minimum eigenvalue per alpha partition
      min_eigenvalue = 1d10  ! Start with a high value to find min
!---- Fourrier transform H(R) to H(k)
      ene=0d0
      
      do ipart=1,npartitions
      
         write(*,'(a,i5)') 'Partition=',ipart
         alpha= 0.39 + 0.11*float(ipart-1)/float(npartitions-1)
         ! Initialize Hamiltonians for the current partition
         do k=1,(nx+1)*(ny+1)*(nz+1)
                HK_trivial=(0d0,0d0)
            HK_topological=(0d0,0d0)

            ! Fourier transform terms
            do j=1,nr
   
                   phase_trivial=0.0d0
               phase_topological=0.0d0
   
               ! Compute phase factors
               do i=1,3
                 phase_trivial=phase_trivial    +dot_product(mesh2(i,k)*bvec(:,i),rvec_trivial(i,j)*avec(:,i))
                 phase_topological=phase_topological+dot_product(mesh2(i,k)*bvec(:,i),rvec_topological(i,j)*avec(:,i))
                  

               enddo
               ! Sum over H(R) contributions for each Hamiltonian
               do i1=1,nb
                  do i2=1,nb
                     Hk_trivial(i1,i2)=Hk_trivial(i1,i2)+Hamr_trivial(i1,i2,j)* &
                                       dcmplx(cos(phase_trivial),-sin(phase_trivial))/float(ndeg_trivial(j))
                     Hk_topological(i1,i2)=Hk_topological(i1,i2)+Hamr_topological(i1,i2,j)* &
                                       dcmplx(cos(phase_topological),-sin(phase_topological))/float(ndeg_topological(j))					

                  enddo
               enddo
            enddo
            
            HK=alpha*HK_trivial+(1d0-alpha)*HK_topological
            H = HK
            do i=1, nb/2
                  H(i,i)=H(i,i)+Hm(1,1)
                  H(i+nb/2,i)=H(i+nb/2,i)+Hm(2,1)
                  H(i,i+nb/2)=H(i,i+nb/2)+Hm(1,2)
                  H(i+nb/2,i+nb/2)=H(i+nb/2,i+nb/2)+Hm(2,2)
            enddo
            !call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
            call zheev('V','U',nb,H,nb,ene(:,k),work,lwork,rwork,info)
         enddo
         
!------calcualte gap and Fermi level
          temp_index=0
          do i = 1,(nx+1)*(ny+1)*(nz+1)
          if((ene(3,i)-ene(2,i)).lt.0.01) temp_index = temp_index+1
          enddo
	  if(allocated(indices)) deallocate(indices)
          allocate(indices(temp_index))
          temp_index = 0
          do i = 1,(nx+1)*(ny+1)*(nz+1)
          if((ene(3,i)-ene(2,i)).lt.0.01) temp_index = temp_index + 1 
          if((ene(3,i)-ene(2,i)).lt.0.01) indices(temp_index) = i
          enddo

!------export eigenvalus 
         open(2000,file='trajectories.dat', status = 'unknown', action = 'write', position = 'append')
         do i=1,size(indices)
         write(2000, '(5(x,f12.6))') alpha, mesh2(1,indices(i)),mesh2(2,indices(i)) & 
         ,mesh2(3,indices(i)),ene(3,indices(i))-ene(2,indices(i))
         enddo
         
         close(2000)
      call cpu_time(part_time)
      part_time2 = part_time/60
      print'(A, F6.2)', "Total runtime (minutes): ", part_time2
      enddo

333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),' not found'
      stop
end program trajectory_xyz


