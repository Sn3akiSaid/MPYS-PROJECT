!==============================================================
! This program adibatically interpolates the energy between two
! topological phases of BiTeI.
!==============================================================
Program interpolate_topology
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=20,npartitions=10 !Adjust these parameters to obtain better resolution around alphacrit and see points closer to an effectively closed gap
!------------------------------------------------------
      integer*4 ik,ipart
      real*8 alpha,ef(npartitions),gap(npartitions)
      character(len=30)::klabel(nkpath),kxlabel(nkpath),kylabel(nkpath)
      character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber
      integer*4,parameter::nk_LAL=(nkpath-1)*np+1, nk_HAH=(nkpath-1)*np+1,nk=(nkpath-1)*np+1
      integer*4 i,j,k,n,nr,i1,i2,nb,lwork,info
      real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
      real*8,parameter::third=1d0/3d0
      real*8 phase_trivial,phase_topological,twopi,jk,a,b
      real*8 klist(3,1:nk),kxlist(3,1:nk_LAL),kylist(3,1:nk_HAH)
      real*8 xk(np**2),yk(nk_HAH),kpath(3,np),kxpath(3,np),kypath(3,np)
      real*8 bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),ykl(nkpath)
      real*8 kxtemp1(3),kxtemp2(3),kytemp1(3),kytemp2(3),kmesh(np,np)
      real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),ene(:,:),rwork(:)
      real(8) :: mesh_kx(3,np, np), mesh_ky(3,np, np), mesh_gap(3, np**2)
      integer*4,allocatable:: ndeg_trivial(:),ndeg_topological(:)
      complex*16,allocatable:: Hk(:,:),HK_trivial(:,:),HK_topological(:,:),Hamr_trivial(:,:,:)
      complex*16,allocatable:: Hamr_topological(:,:,:),work(:)
      	real(8), parameter :: x_min = -0.1, x_max = 0.1, y_min = -0.1, y_max = 0.1
    real(8), dimension(:), allocatable :: x_vals, y_vals
    real(8), dimension(:,:), allocatable :: mesh
    integer :: total_pairs
!------------------------------------------------------
      write(hamil_file_trivial,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
      write(hamil_file_topological,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

      twopi=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(98,*)bvec
!---------------kpath
      data kpath(:,1) /    -0.1d0,      0.0d0,    0.5d0/  !L
      data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
      data kpath(:,3) /     0.1d0,      0.0d0,    0.5d0/  !L
      data kpath(:,4) /     0.0d0,      0.1d0,    0.5d0/  !H
      data kpath(:,5) /     0.0d0,     -0.1d0,    0.5d0/  !H

      data kxlabel     /'L','A','L'/
      data kylabel     /'H','A','H'/

      kxtemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+&
                 (kpath(2,1)-kpath(2,2))*bvec(:,2)+&
                 (kpath(3,1)-kpath(3,2))*bvec(:,3)

      

      kytemp1(:)= (kpath(1,5)-kpath(1,2))*bvec(:,1)+&
                  (kpath(2,5)-kpath(2,2))*bvec(:,2)+&
                  (kpath(3,5)-kpath(3,2))*bvec(:,3)

      yk(1)= -sqrt(dot_product(kytemp1,kytemp1))
      ykl(1)=yk(1)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Generate x and y values
    allocate(x_vals(np), y_vals(np))
    do i =1,np
    	x_vals(i) = x_min + i * (x_max - x_min)/np
    	y_vals(i) = y_min + i * (y_max - y_min)/np
end do

    ! Calculate the total number of pairs
    allocate(mesh(3, np**2))

    ! Loop through all x values and all y values
    ! Store the pairs in the result array
    j = 1
    do i = 1, np
        do n = 1, np
            mesh(1, j) = x_vals(i)
            mesh(2, j) = y_vals(n)
            mesh(3, j) = -0.5d0
            j = j + 1
        end do
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Creates a list of 100 points inbetween BZ points LAL
      k=0
      kxtemp1=0d0
       do j=1,np**2
        k=k+1
        kxtemp2=mesh(1,k)*bvec(:,1)+mesh(2,k)*bvec(:,2)+mesh(3,k)*bvec(:,3)
        kxtemp1=kxtemp2
       enddo



      kxtemp2=mesh(1,nk)*bvec(:,1)+kxlist(2,nk)*bvec(:,2)+kxlist(3,nk)*bvec(:,3)

      mesh=mesh*twopi
      

!------read trivial H(R)
      open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec_trivial(3,nr),Hk_trivial(nb,nb),Hamr_trivial(nb,nb,nr),&
               ndeg_trivial(nr),ene(nb,np**2))
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
         alpha=float(ipart-1)/float(npartitions-1)
         ! Initialize Hamiltonians for the current partition
         do k=1,np**2
                HK_trivial=(0d0,0d0)
            HK_topological=(0d0,0d0)

            ! Fourier transform terms
            do j=1,nr
   
                   phase_trivial=0.0d0
               phase_topological=0.0d0
   
               ! Compute phase factors
               do i=1,3
                      phase_trivial=phase_trivial    +mesh(i,k)*rvec_trivial(i,j)
                  phase_topological=phase_topological+mesh(i,k)*rvec_topological(i,j)
                  
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
            
            call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)

            ! Find minimum eigenvalue for the current k-point and update min_eigenvalue
            min_eigenvalue(ipart) = min(min_eigenvalue(ipart), minval(ene(:, k)))
            
         enddo
         
!------calcualte gap and Fermi level
         gap(ipart)= minval(ene(13,:))-maxval(ene(12,:))
          ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0
        

!------export eigenvalus 
         write(partnumber,'(i5)') ipart
         write(line,'(3a)') 'mesh_partition_',trim(adjustl(partnumber)),'.dat' 
         open(100,file=trim(line))
         open(1000, file='mesh_ene.dat', action = 'write')
            do k=1,np**2
              write(100,'(3(x,f12.6))') mesh(1,k)/twopi,mesh(2,k)/twopi,ene(13,k)-ene(12,k) !This now writes into the files the coordinates as a function of the TCB and BCB energy difference
            enddo
              write(100,*)
              write(100,*)
         do k=1, np**2
            if(ipart.eq.10) write(1000,'(3(x,f12.6))') mesh(1,k)/twopi,mesh(2,k)/twopi,ene(13,k)-ene(12,k)
         enddo
         close(100)
        close(1000)
!------- export gap
         write(777,'(f10.8,f8.5)') alpha,gap(ipart) !this had ipart in front of alpha
      enddo
    
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),' not found'
      stop
end program interpolate_topology
