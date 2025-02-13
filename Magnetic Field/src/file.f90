!==============================================================
! This program adibatically interpolates the energy between two
! topological phases of BiTeI.
! Calculates a magnetic perturbation
! Calculates spins for a given Hamiltonian
!==============================================================
Program interpolate_topology
    Implicit None
!--------to be midified by the user
    character(len=80):: prefix="BiTeI"
    integer,parameter::nkpath=3,np=100,npartitions=36
!---Magnetic Field to be modified by User
    real*8,parameter::B_x=0d0, B_y=0d0, B_z=0.1d0
!------------------------------------------------------
    integer ik, ipart, is, ib
    real*8 alpha,ef(npartitions),gap(npartitions)!,write_values(12:13)!,&
           !start_alpha,end_alpha,step_size!, B_field(3)
    character(len=30)::klabel(nkpath)
    character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber
    integer,parameter::nk=(nkpath-1)*np+1
    integer i,j,k,&
            nr,nb,&
            i1,i2,&
            lwork,info!,steps
    real*8,parameter::third=1d0/3d0,mthird=-1d0/3d0,&
                      twothird=2d0/3d0,mtwothird=-2d0/3d0

    real*8 phase_trivial,phase_topological,&
           twopi,jk,a,b,&
           spin_x(1,1),spin_y(1,1),spin_z(1,1),&
           spin_xp(1,1),spin_yp(1,1),spin_zp(1,1),&
           klist(3,1:nk),kpath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),&
           xk(nk),xkl(nkpath)

    real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),&
                         ene(:,:),enep(:,:),&
                         spin(:,:,:),spinp(:,:,:),&
                         rwork(:)

    integer,allocatable:: ndeg_trivial(:),ndeg_topological(:)

    complex*16 sigx(2, 2), sigy(2, 2), sigz(2, 2),&
               chi(2,1), chip(2,1)

    complex*16,allocatable::H(:,:), Hk(:,:), Hm(:,:), Hmag(:,:),&
                            HK_trivial(:,:), HK_topological(:,:),&
                            Hamr_trivial(:,:,:), Hamr_topological(:,:,:),&
                            work(:)
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
!ky
    !data kpath(:,1) /     0.1d0,  -0.2d0,        0.5d0/  !L
    !data kpath(:,2) /     0.0d0,  0.0d0,    0.5d0/  !A
    !data kpath(:,3) /     -0.1d0,  0.2d0,           0.5d0/  !H

!kx
   ! data kpath(:,1) /     -0.1d0,  0.0d0,        0.5d0/  !L
   ! data kpath(:,2) /     0.0d0,  0.0d0,    0.5d0/  !A
   ! data kpath(:,3) /     0.1d0,  0.0d0,           0.5d0/  !H

!LAH
    data kpath(:,1) /     -0.1d0,  0.1d0,        0.5d0/  !L
    data kpath(:,2) /     0.0d0,  0.0d0,    0.5d0/  !A
    data kpath(:,3) /     0.1d0,  0.1d0,           0.5d0/  !H

    data sigx /(0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
    data sigy /(0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
    data sigz /(1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/

    !data B_field / B_x, B_y, B_z/ !Magnetic field vector
    allocate(Hm(2,2),Hmag(18,18))
    Hm = B_x*sigx + B_y*sigy + B_z*sigz

    data klabel     /'L','A','H'/

    ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

    xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
    xkl(1)=xk(1)
    
!Creates a list of 100 points inbetween BZ points LAL
    k=0
    ktemp1=0d0
    do i=1,nkpath-1
     do j=1,np
      k=k+1
      jk=dfloat(j-1)/dfloat(np)
      klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
      ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
      if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      if(j.eq.1) xkl(i)=xk(k)
      ktemp1=ktemp2
     enddo
    enddo
    klist(:,nk)=kpath(:,nkpath)
    ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
    xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
    xkl(nkpath)=xk(nk)
    klist=klist*twopi

!------read trivial H(R)
    open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_trivial(3,nr),Hk_trivial(nb,nb),Hamr_trivial(nb,nb,nr),&
             ndeg_trivial(nr),ene(nb,nk),enep(nb,nk))
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
    enddo
    close(99)

!------ LAPACK-related array allocations
    lwork=max(1,2*nb-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
    allocate(spin(3,nb,(np+1)**2),spinp(3,nb,(np+1)**2))
!------ open gap file
    open(777,file='gap.dat')

    !start_alpha=0d0
   ! end_alpha=1d0
   ! step_size=0.1d0
    
   ! steps=int((end_alpha-start_alpha)/step_size)
!---- Fourrier transform H(R) to H(k)
   ! ene=0d0
   ! do ipart=1,steps+1
    !   write(*,'(a,i5)') 'Partition=',ipart
       !alpha=float(ipart-1)/float(npartitions-1)
     !  alpha=start_alpha+(ipart-1)*step_size
      ! print *, alpha
    !Hmag=(0d0,0d0)
    do i=1, nb/2
        !do j=1,nb/2
            Hmag(i,i)=Hm(1,1)
            Hmag(i,i+nb/2)=Hm(1,2)
            Hmag(i+nb/2,i)=Hm(2,1)
            Hmag(i+nb/2,i+nb/2)=Hm(2,2)
       ! enddo
    enddo
    !---- Fourrier transform H(R) to H(k)
    ene=0d0
    do ipart=30,npartitions
       write(*,'(a,i5)') 'Partition=',ipart
       !alpha=0.88888
       alpha=float(ipart-1)/float(npartitions-1)
       do k=1,nk
              HK_trivial=(0d0,0d0)
          HK_topological=(0d0,0d0)
          do j=1,nr
 
                 phase_trivial=0.0d0
             phase_topological=0.0d0
 
             do i=1,3
                    phase_trivial=phase_trivial    +klist(i,k)*rvec_trivial(i,j)
                phase_topological=phase_topological+klist(i,k)*rvec_topological(i,j)
             enddo
 
             do i1=1,nb
                do i2=1,nb
                   Hk_trivial(i1,i2)=Hk_trivial(i1,i2)+Hamr_trivial(i1,i2,j)* &
                                     dcmplx(cos(phase_trivial),-sin(phase_trivial))/float(ndeg_trivial(j))
                   Hk_topological(i1,i2)=Hk_topological(i1,i2)+Hamr_topological(i1,i2,j)* &
                                     dcmplx(cos(phase_topological),-sin(phase_topological))/float(ndeg_topological(j))
                enddo
             enddo
          enddo

          HK=alpha*HK_topological+(1d0-alpha)*HK_trivial

          H=Hk+Hmag
          

          call zheev('V','U',nb,H,nb,enep(:,k),work,lwork,rwork,info)
          call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)

          do ib=1,nb
            do is=1,nb/2
                  chi(1,1) = Hk(is     ,ib) 
                  chi(2,1) = Hk(is+nb/2,ib)
                  chip(1,1) = H(is     ,ib) 
                  chip(2,1) = H(is+nb/2,ib)
                  
                  spin_x = matmul(conjg(transpose(chi)),matmul(sigx, chi))
                  spin_y = matmul(conjg(transpose(chi)),matmul(sigy, chi))
                  spin_z = matmul(conjg(transpose(chi)),matmul(sigz, chi))
                  spin(1,ib,k)=spin(1,ib,k)+spin_x(1,1)
                  spin(2,ib,k)=spin(2,ib,k)+spin_y(1,1)
                  spin(3,ib,k)=spin(3,ib,k)+spin_z(1,1)
                  !Calculate spins for Perturbed Hamiltonian
                  spin_xp = matmul(conjg(transpose(chip)),matmul(sigx, chip))
                  spin_yp = matmul(conjg(transpose(chip)),matmul(sigy, chip))
                  spin_zp = matmul(conjg(transpose(chip)),matmul(sigz, chip))

                  spinp(1,ib,k)=spinp(1,ib,k)+spin_xp(1,1)
                  spinp(2,ib,k)=spinp(2,ib,k)+spin_yp(1,1)
                  spinp(3,ib,k)=spinp(3,ib,k)+spin_zp(1,1)
                  
            enddo
      enddo
       enddo
       
!------calcualte gap and Fermi level
       gap(ipart)= minval(enep(13,:))-maxval(enep(12,:))
        ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0

!------export eigenvalus 
       write(partnumber,'(i5)') ipart
       write(line,'(3a)') 'band_partition_',trim(adjustl(partnumber)),'.dat' 
       open(100,file=trim(line))
       !open(200,file="unperturbed_spins.dat")
      ! open(300,file="perturbed_spins.dat")
       do i = 11, 14
          do k=1,nk
         !       ! Check if ene(i,k) - ef(ipart) is less than 0.01 and set it to 0 if true
          !      if (abs(ene(i,k) - ef(ipart)) .lt. 0.01) then
           !     write_values = 0.0
            !    else
             !   write_values = ene(i,k) - ef(ipart)
              !  end if
                ! Write the values to the file
                write(100, '(5(x,f12.6))') xk(k), ene(i,k)-ef(ipart),&
                                           spin(2,i,k)/(sqrt(spin(1,i,k)**2 +spin(2,i,k)**2 +spin(3,i,k)**2)),&
                                           enep(i,k)-ef(ipart),&
                                           spinp(2,i,k)/(sqrt(spinp(1,i,k)**2 +spinp(2,i,k)**2 +spinp(3,i,k)**2))
                !write(200,'(4(x,f12.6))') spin(3,i,k),sqrt(spin(1,i,k)**2 +spin(2,i,k)**2 +spin(3,i,k)**2)
                !write(300,'(4(x,f12.6))') spinp(3,i,k),sqrt(spinp(1,i,k)**2 +spinp(2,i,k)**2 +spinp(3,i,k)**2)

             enddo
        write(100,*)
        write(100,*)
!        write(200,*)
 !       write(200,*)
  !      write(300,*)
   !     write(300,*)

        enddo
       close(100)
    !   close(200)
     !  close(300)
!------- export gap
       write(777,'(f10.8,f8.5)') alpha,gap(ipart) !this had ipart in front of alpha
    enddo
!------- call plt script of gnuploting
    !call write_plt(nkpath,xkl,klabel)
!------- stop
    stop
!------- Alarms go off!
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
    stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),' not found'
    stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),' not found'
    stop

    end program interpolate_topology