!==============================================================
! This program adibatically interpolates the energy between two
! topological phases of BiTeI.
!==============================================================
Program interpolate_topology
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=100,npartitions=30
!------------------------------------------------------
      integer*4 ik,ipart
      real*8 alpha,ef(npartitions),gap(npartitions)
      character(len=30)::klabel(nkpath),kxlabel(nkpath),kylabel(nkpath)
      character(len=80) hamil_file_trivial,hamil_file_topological,nnkp,line,partnumber
      integer*4,parameter::nk_LAL=(nkpath-1)*np+1, nk_HAH=(nkpath-1)*np+1,nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,i1,i2,nb,lwork,info
      real*8,dimension(npartitions) :: min_eigenvalue, alpha_values
      real*8,parameter::third=1d0/3d0
      real*8 phase_trivial,phase_topological,twopi,jk,a,b
      real*8 klist(3,1:nk),kxlist(3,1:nk_LAL),kylist(3,1:nk_HAH),xk(nk_LAL),yk(nk_HAH),kpath(3,np),kxpath(3,np),kypath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),ykl(nkpath),kxtemp1(3),kxtemp2(3),kytemp1(3),kytemp2(3)
      real*8,allocatable:: rvec_trivial(:,:),rvec_topological(:,:),ene(:,:),rwork(:)
      integer*4,allocatable:: ndeg_trivial(:),ndeg_topological(:)
      complex*16,allocatable:: Hk(:,:),HK_trivial(:,:),HK_topological(:,:),Hamr_trivial(:,:,:)
      complex*16,allocatable:: Hamr_topological(:,:,:),work(:)
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
      data kpath(:,4) /     0.1d0,      0.1d0,    0.5d0/  !H
      data kpath(:,5) /     0.1d0,     -0.1d0,    0.5d0/  !H

      data kxlabel     /'L','A','L'/
      data kylabel     /'H','A','H'/

      kxtemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+&
                 (kpath(2,1)-kpath(2,2))*bvec(:,2)+&
                 (kpath(3,1)-kpath(3,2))*bvec(:,3)

      xk(1)= -sqrt(dot_product(kxtemp1,kxtemp1))
      xkl(1)=xk(1)

      kytemp1(:)= (kpath(1,5)-kpath(1,2))*bvec(:,1)+&
                  (kpath(2,5)-kpath(2,2))*bvec(:,2)+&
                  (kpath(3,5)-kpath(3,2))*bvec(:,3)

      yk(1)= -sqrt(dot_product(kytemp1,kytemp1))
      ykl(1)=yk(1)
      
!Creates a list of 100 points inbetween BZ points LAL
      k=0
      kxtemp1=0d0
      do i=1,nkpath-1
       do j=1,np
        k=k+1
        jk=dfloat(j-1)/dfloat(np)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        kxtemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(kxtemp2-kxtemp1,kxtemp2-kxtemp1))
        if(j.eq.1) xkl(i)=xk(k)
        kxtemp1=kxtemp2
       enddo
      enddo

!Creates a list of 100 points inbetween BZ points HAH
      k=0
      kytemp1=0d0
      do i=1,nkpath-1
       do j=1,np
        k=k+1
        jk=dfloat(j-1)/dfloat(np)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        kytemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) yk(k)=yk(k-1)+sqrt(dot_product(kytemp2-kytemp1,kytemp2-kytemp1))
        if(j.eq.1) ykl(i)=yk(k)
        kytemp1=kytemp2
       enddo
      enddo


      klist(:,nk)=kxpath(:,nkpath)+kypath(:,nkpath)
      kxtemp2=kxlist(1,nk)*bvec(:,1)+kxlist(2,nk)*bvec(:,2)+kxlist(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(kxtemp2-kxtemp1,kxtemp2-kxtemp1))
      xkl(nkpath)=xk(nk)
      kytemp2=kylist(1,nk)*bvec(:,1)+kylist(2,nk)*bvec(:,2)+kylist(3,nk)*bvec(:,3)
      yk(nk)=yk(nk-1)+sqrt(dot_product(kytemp2-kytemp1,kytemp2-kytemp1))
      ykl(nkpath)=yk(nk)
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
         do k=1,nk
                HK_trivial=(0d0,0d0)
            HK_topological=(0d0,0d0)

            ! Fourier transform terms
            do j=1,nr
   
                   phase_trivial=0.0d0
               phase_topological=0.0d0
   
               ! Compute phase factors
               do i=1,3
                      phase_trivial=phase_trivial    +klist(i,k)*rvec_trivial(i,j)
                  phase_topological=phase_topological+klist(i,k)*rvec_topological(i,j)
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
         write(line,'(3a)') 'band_partition_',trim(adjustl(partnumber)),'.dat' 
         open(100,file=trim(line))
         do i=1,nb
            do k=1,nk
              write(100,'(2(x,f12.6))') xk(k),ene(i,k)-ef(ipart)
            enddo
              write(100,*)
              write(100,*)
         enddo
         close(100)
!------- export gap
         write(777,'(f10.8,f8.5)') alpha,gap(ipart) !this had ipart in front of alpha
      enddo
    
     
!------- call plt script of gnuploting
      call write_plt(nkpath,xkl,klabel)
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

!=============================================
! This subroutine generates a gnuplot script
! to plot the band structures
!============================================!
      subroutine write_plt(nkp,xkl,kl)
      implicit none
      integer nkp,i
      real*8 xkl(nkp)
      character(len=30)kl(nkp)
      
      open(99,file='band.plt')
      write(99,'(a)') 'set xtics ( \'
      do i=1,nkp
         if(trim(adjustl(kl(i))).eq.'g'.or.trim(adjustl(kl(i))).eq.'G')kl(i)="{/Symbol \107}"
         if(i.ne.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i),", \"
         if(i.eq.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i)," )"
      enddo
      write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [',xkl(1),':',xkl(nkp),']'
      write(99,'(a)') &
           'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in'
      write(99,'(a,f4.2,a)')'set output "band.pdf"'
           write(99,'(9(a,/),a)') &
           'set encoding iso_8859_1',&
           'set size ratio 0 1.0,1.0',&
           'set ylabel "E-E_{VBM} (eV)"',&
           'set yrange [ -0.5 : 0.5 ]',&
           'unset key',&
           'set ytics 1.0 scale 1 nomirror out',&
           'set mytics 2',&
           'set terminal qt persist',&
           'set parametric',&
           'set trange [-10:10]',&
           'plot for [i=1:20]sprintf("band_partition_%d.dat", i) u 1:2 with l lt 1 lw 0.5,\'
     do i=2,nkp-1
       write(99,'(f12.6,a)') xkl(i),',t with l lt 2  lc -1,\'
     enddo
     write(99,'(a)') 't,0 with l lt 2  lc -1'
     end subroutine write_plt
