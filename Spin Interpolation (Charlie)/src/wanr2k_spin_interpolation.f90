program wannier_band_structure
    implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    integer,parameter::nkpath=8,np=100, npartitions = 10
!------------------------------------------------------ Assign values
    integer*4 ik,ikmax
    real*8 kz,alpha,ef(npartitions),gap(npartitions)
    character(len=30)::klabel(nkpath)
    character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber
    integer*4,parameter::nk=(nkpath-1)*np+1
    integer*4 i,j,k,nr,i1,i2,nb,lwork,info,ipart,o,p
    real*8,parameter::third=1d0/3d0!,kz=0d0
    real*8 phase_trivial,phase_topological,twopi,jk,a,b
    real*8 klist(3,1:nk),xk(nk),kpath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath)
    real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),ene(:,:),rwork(:)
    integer*4,allocatable:: ndeg_trivial(:),ndeg_topological(:)
    complex*16,allocatable:: Hk_trivial(:,:),Hk_topological(:,:),Hamr(:,:,:),work(:),Hk(:,:)
    complex*16,allocatable:: Hamr_trivial(:,:,:), Hamr_topological(:,:,:)
    real*16,dimension(2,2):: sigx(1:2,1:2),sigz(1:2,1:2)
    complex*16,dimension(2,2):: sigy(1:2,1:2)
    complex*16 temp1,temp2
    real*8, allocatable:: weightx(:,:),weighty(:,:),weightz(:,:)
!------------------------------------------------------
      write(hamil_file_trivial,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
      write(hamil_file_topological,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
    twopi=4.0d0*atan(1.0d0)*2.0d0
    sigx(1,1) = (0.0d0, 0.0d0)
    sigx(2,1) = (1.0d0, 0.0d0)
    sigx(1,2) = (1.0d0, 0.0d0)
    sigx(2,2) = (0.0d0, 0.0d0) 
    
    sigy(1,1) = (0.0d0, 0.0d0)
    sigy(2,1) = dcmplx(0.0d0, -1.0d0)
    sigy(1,2) = dcmplx(0.0d0, 1.0d0)
    sigy(2,2) = (0.0d0, 0.0d0)
    
    sigz(1,1) = (1.0d0, 0.0d0)
    sigz(2,1) = (0.0d0, 0.0d0)
    sigz(1,2) = (0.0d0, 0.0d0)
    sigz(2,2) = (-1.0d0, 0.0d0)


!---------------  reciprocal vectors
    open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
    
    read(98,*)bvec
!---------------kpath
    data kpath(:,1) /    -0.1d0,      0.1d0,    0.5d0/  !L
    data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
    data kpath(:,3) /     0.1d0,      0.1d0,    0.5d0/  !H

!      open(77,file='tmp')
!     read(77,*)ik,ikmax
!     kz=float(ik)*0.5d0/float(ikmax)
!     kpath(3,:)=kz

    data klabel     /'G','M','K','G','A','L','H','A'/

    ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

!      xk(1)= 0d0 !-sqrt(dot_product(ktemp1,ktemp1))
    xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
    xkl(1)=xk(1)
    data kpath(:,1) /    -0.1d0,      0.1d0,    0.5d0/  !L
    data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
    data kpath(:,3) /     0.1d0,      0.1d0,    0.5d0/  !H

    

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
!      write(*,*)klist
    klist=klist*twopi

!------read trivial H(R)
      open(99,file=trim(adjustl(hamil_file_trivial)),err=444)
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec_trivial(3,nr),Hk_trivial(nb,nb),Hamr_trivial(nb,nb,nr),&
               ndeg_trivial(nr),ene(nb,nk))
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

   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
   allocate(weightx(nb,nk),weighty(nb,nk),weightz(nb,nk))
   
   !---- Fourrier transform H(R) to H(k)
      ene=0d0
      do ipart=1,npartitions
      	weightx = (0d0,0d0)
        weighty = (0d0,0d0)
        weightz = (0d0,0d0)
         write(*,'(a,i5)') 'Partition=',ipart
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

            HK=alpha*HK_trivial+(1d0-alpha)*HK_topological
   
            call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
	
        do i = 1,7,3
        do j = 1,nb
        
        o = i+1
        p = i+2
	  weightx(j,k) = weightx(j,k) + conjg(HK(j,i))*(sigx(1,1)*Hk(j,i) + sigx(2,1)*Hk(j,i+9)) +&
	  	conjg(Hk(j,i+9))*(sigx(1,2)*Hk(j,i) + sigx(2,2)*Hk(j,i+9))
	  	
	  weighty(j,k) = weighty(j,k) + conjg(HK(j,o))*(sigy(1,1)*Hk(j,o) + sigy(2,1)*Hk(j,o+9)) &
	  	+ conjg(Hk(j,o+9))*(sigy(1,2)*Hk(j,o) + sigy(2,2)*Hk(j,o+9))
	  	
	  weightz(j,k) = weightz(j,k) + conjg(HK(j,p))*(sigz(1,1)*Hk(j,p) + sigz(2,1)*Hk(j,p+9)) &
	  	+ conjg(Hk(j,p+9))*(sigz(1,2)*Hk(j,p) + sigz(2,2)*Hk(j,p+9))
	  	
       enddo
       enddo

    enddo
            

     
!------calcualte gap and Fermi level
         gap(ipart)= minval(ene(13,:))-maxval(ene(12,:))
          ef(ipart)=(minval(ene(13,:))+maxval(ene(12,:)))/2d0

!------export eigenvalus 
         write(partnumber,'(i5)') ipart
         write(line,'(3a)') 'spin_partition_',trim(adjustl(partnumber)),'.dat' 
         open(100,file=trim(line))
         do i=1,nb
            do k=1,nk
              write(100,'(5(x,f12.6))') xk(k),ene(i,k)-ef(ipart), weightx(i,k), weighty(i,k), weightz(i,k)
            enddo
              write(100,*)
              write(100,*)
         enddo
         close(100)
!------- export gap
         write(777,'(f10.8,f8.5)')alpha,gap(ipart)
         enddo

333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_trivial)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_topological)),' not found'
      stop
      end
