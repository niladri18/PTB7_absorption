program ppv

        use parameter
        use readfile

!open(unit=1,file="dipole_n.inp")
open(unit=10,file="./output/nonlinear.out")

print*,"Calculates the energy levels of PTB7 in the huckel limit"
print*,"then calculates the Hartree-Fock energy"
print*,"then calculates excited sates using SCI." 

print*, "HF calculations starts"

! ==== Reading Coordinates =====
   open(unit=90,file="./input/ptb7cord4.inp")
   open(unit=91,file="./input/unit4.inp")
   open(unit=25,file="./output/hf_ppt.out")
   read(90,*) n

   sulp=n/5
   !sulp=0
   nn=(n/2+sulp/2)*(n/2-sulp/2)
   lwork=3*n-1
   
allocate(arr(1:n,1:n))
allocate(f(1:n,1:n))
allocate(f0(1:n,1:n))
!allocate(v(n,n)) 
allocate(v_n(1:n,1:n))
allocate(charge(1:n))
allocate(p_on(1:n))
allocate(v_hf(1:n,1:n))
allocate(d_hf(1:n))
allocate(d(1:n))
allocate(d_n(1:n))
allocate(v_ij(1:n,1:n))
allocate(x_pos(1:n))
allocate(y_pos(1:n))
allocate(eps(1:n))
allocate(z(1:n))
allocate(z_pos(1:n))
allocate(U_c(1:n))
allocate(tmp(1:n))
allocate(vtmp(1:n,1:n))
allocate(h_s(1:nn,1:nn))
allocate(m_states(1:nn,1:nn))
allocate(E_tot(1:nn))
allocate(tsc(1:nn))
allocate(del(1:nn))
allocate(vsc(1:nn,1:nn))
allocate(dipole_x(1:nn))
allocate(dipole_y(1:nn))
allocate(dipole_z(1:nn))
allocate(dip_mag(1:nn))
allocate(s_index(1:n/2+sulp/2,n/2+sulp/2+1:n))
!allocate(s_index(n/2-sulp/2+1:n,n/2-sulp/2,1                                 ))
allocate(work(1:lwork))


 !do lv3=n/2-sulp/2+1,n
 !do lv4=n/2-sulp/2,1,-1
  dip_mag=0.0


   write(*,*) 'total number of atoms =',n,nn,sulp
   do i=1,n
   read(90,*) x_pos(i),y_pos(i),z_pos(i)
   end do
   print*,"Ended reading the coordinates"
! ==== End reading Coordinates ====

allocate(ilbond(2,200))

 do i=1,2
 do j=1,200
 ilbond(i,j)=0.0
 end do
 end do


! ====== Reading Main File ====

        call read_config
 
 !read(91,*) n
! read(91,*) nbonds1,nbonds2,nbonds3,nbonds4!,nbonds5 
! nbnd=nbonds1+nbonds2+nbonds3+nbonds4!+nbonds5
! print*,"total no of bonds",nbnd
 !nbonds1 = C-S; nbonds2 = C-C; nbonds3 = C=C; nbonds4 = benzene bonds; nbonds5 = interchain bonds
! read(91,*) t1,t2,t3,tb,kappa,eta
 !write(*,*) 't1,t2,t3',t1,t2,t3
 !Define this array variable: ilbnd(2,100)
! read(91,*) (ilbond(1,k),ilbond(2,k),k=1,nbnd)
 !write(*,*) 'links',(ilbond(1,k),ilbond(2,k),k=1,nbnd)
 !Define siten, Hub-U, and chem-pot arrays
! read(91,*) (eps(i),i=1,n)
! read(91,*) (U_c(i),i=1,n)
! read(91,*) (z(i),i=1,n)
 !write(*,*) 'z',(z(i),i=1,n)


! print*,"HOMO",n/2+sulp/2,"LUMO",n/2+sulp/2+1
! do i=1,n
! print*,x_pos(i),z(i),U(i)
! end do
 

! print*,"ended reading coordinates"
 print*, "Constructing Huckel matrix"
! ========= Constructing Huckel Matrix =======
 do i=1,n
 do j=1,n
 arr(i,j)=0.0
 enddo
 enddo
 do i=1,n
 arr(i,i)=eps(i)
 enddo
 do i=1,nbonds1
 i1=ilbond(1,i)
 i2=ilbond(2,i)
 arr(i1,i2)=t1
 arr(i2,i1)=t1
 !print*,i1,i2,arr(i1,i2),arr(i2,i1)
 enddo 
 do i=nbonds1+1,nbonds1+nbonds2
 i1=ilbond(1,i)
 i2=ilbond(2,i)
 arr(i1,i2)=t2
 arr(i2,i1)=t2
 !print*,i1,i2,arr(i1,i2),arr(i2,i1)
 enddo 
 do i=nbonds1+nbonds2+1,nbonds1+nbonds2+nbonds3
 i1=ilbond(1,i)
 i2=ilbond(2,i)
 arr(i1,i2)=t3
 arr(i2,i1)=t3
 !print*,i1,i2,arr(i1,i2),arr(i2,i1)
 enddo 
 
 do i=nbonds1+nbonds2+nbonds3+1,nbonds1+nbonds2+nbonds3+nbonds4
 i1=ilbond(1,i)
 i2=ilbond(2,i)
 arr(i1,i2)=tb
 arr(i2,i1)=tb
 !print*,i1,i2,arr(i1,i2),arr(i2,i1)
 enddo
 
 !do i=nbonds1+nbonds2+nbonds3+nbonds4+1,nbonds1+nbonds2+nbonds3+nbonds4+nbonds5
 !i1=ilbond(1,i)
 !i2=ilbond(2,i)
 !arr(i1,i2)=t3
 !arr(i2,i1)=t3
! print*,i1,i2,arr(i1,i2),arr(i2,i1)
 !enddo
! =========================================

 print*,"ended constructing huckel matrix"

! ========= Constructing Huckel Matrix =======
! do i=1,natoms
! do j=1,natoms
! hucmat(i,j)=0.0d0
! enddo
! enddo
! do i=1,natoms
! hucmat(i,i)=siten(i)
! enddo
! do i=1,nbonds1
! i1=ilbond(1,i)
! i2=ilbond(2,i)
! hucmat(i1,i2)=t1
! enddo 
! do i=1,nbonds2
! i1=ilbond(1,i)
! i2=ilbond(2,i)
! hucmat(i1,i2)=t2
! enddo 
! do i=1,nbonds3
! i1=ilbond(1,i)
! i2=ilbond(2,i)
! hucmat(i1,i2)=t3
! enddo 
! =========================================




	!do i=1,n
	!do j=1,n
		!if(arr(i,j).ne.0.0)then
		!print*,i,j,arr(i,j),arr(j,i)
		!end if
	!end do
	!end do



do i=1,n
do j=1,n
f(i,j)=0.0
f(j,i)=0.0
v_ij(i,j)=0.0
v_ij(j,i)=0.0
end do
end do



	f0=arr





 call DSYEV('V','U', N, arr, N, tmp, WORK, lwork, info )
 
 
 do i=1,n
 print*,i,tmp(i)
 end do

 print*,"Huckel calculation ends"
 
v_n=arr
d_n=tmp
tmp=0.0
vtmp=0.0


it=0


!  Calculating the coulomb interactions 
do i=1,n
 v_ij(i,i)=U_c(i)
 do j=1,n


r_ij=(x_pos(i)-x_pos(j))**2 + (y_pos(i)-y_pos(j))**2 + (z_pos(i)-z_pos(j))**2

if(i.ne.j)then
 v_ij(i,j)=14.397/(kappa*sqrt((28.794/(U_c(i)+U_c(j)))**2+r_ij))
! v_ij(i,j)=U_c(i)/(kappa*sqrt(1.0+0.611*r_ij))
end if
v_ij(j,i)=v_ij(i,j)
end do
end do


        ! end calculating the coulomb interactions
        
 print*,"Hartree-Fock calculation starts"   

do hf=1,400

!! self consistent calculation starts
        it=it+1


        ! calculating charge density at site i
 do i=1,n
        charge(i)=0.0
        !p_on(i)=0.0
        do k=1,n/2+sulp/2
                charge(i)=charge(i)+2.0*v_n(i,k)*v_n(i,k)
        end do
end do
f=0.0

! Assigning the matrix elements of the many body hamiltonian



 


!!  Assiging off-diagonal elements

	 do i=1,n

	 do j=1,n

	
	if(i.ne.j)then
	
		p_off=0.0
		do k=1,n/2+sulp/2
 			p_off=p_off+2.0*v_n(i,k)*v_n(j,k)
 			!print*,i,j,v_n(i,k)*v_n(j,k)
 		end do
		f(i,j)=f0(i,j)-0.5*p_off*v_ij(i,j)
		f(j,i)=f(i,j)
		
	end if
	
	end do
	end do
		
	


		
 	!! end of assigning the off diagonal elements
 	
 	!! Assigning the diagonal elements
 	
 	
 	
 	
 	
 	
 	do i=1,n
 	p_on(i)=0.0
 	do j=1,n
 		if(i.ne.j)then
 		p_on(i)=p_on(i) + (charge(j)-z(j))*v_ij(i,j)
 		!print*,"Chemical pot",i,j,z(j),charge(j),v_ij(i,j)
 		end if
 	end do
 	f(i,i)=p_on(i)+0.5*charge(i)*U_c(i)+eps(i)
 	end do
 	
 	!!!! end assigning the diagonal elements
 	
	


	
	call DSYEV('V','U', N, f , N, tmp, WORK, lwork, info )


 	
 	
 	err=0.0
 	do i=1,n
 	
 	err1=(tmp(i)-d_n(i))**2
 	err=err+err1
 	!print*,i,j,v_ij(i,j),v_ij(j,i)
 	end do
 	conv=sqrt(err)/n
 	

	v_n =f
 	d_n=tmp







 if(conv.le.0.0000001)exit
end do

 !call eigsrt(d_hf,v_hf)

print*,"converged in steps",it

do i=1,n
write(25,*)tmp(i),i,charge(i)
!print*,tmp(i),i,charge(i)
!print*,d_hf(i),i,charge(i)
end do



!!!! Starting SCI calculations!!!!
 !stop
 
print*, "HF calculations ends"
 
 h_s=0.0
 vtmp=f
 
 
 state1=1
 
 no_states=0
 
 do lv1=1,n/2+sulp/2
 do lv2=n/2+sulp/2+1,n
 
 state2=1
 
 
 s_index(lv1,lv2)=state1
 
 no_states=no_states+1
 
 do lv3=1,n/2+sulp/2
 do lv4=n/2+sulp/2+1,n
 	
 	
 	
 	
 	dir=0.0
 	xchange=0.0
 	
 	do mu=1,n
 	do nu=1,n
 		dir=dir + 2.0*vtmp(mu,lv1)*vtmp(mu,lv2)*vtmp(nu,lv3)*vtmp(nu,lv4)*v_ij(mu,nu)
 		xchange=xchange + 1.0*vtmp(mu,lv2)*vtmp(nu,lv1)*vtmp(nu,lv3)*vtmp(mu,lv4)*v_ij(mu,nu)
 	end do
 	end do
 	
	if(state1==state2)then
	h_s(state1,state2)=tmp(lv4)-tmp(lv3) + dir - xchange
	
	else 
	h_s(state1,state2) = dir - xchange
	
	
	end if

	state2=state2+1
	
 end do
 end do
 
 state2=1
 state1=state1+1
 end do
 end do
 print*,"No of states",no_states
 
 do i=1,nn
 do j=1,nn
 if(h_s(i,j).ne.0.0)then
 write(35,*)i,j,h_s(i,j),h_s(j,i)
 end if
 end do
 end do
 
 lwork=3*nn-1
 deallocate(work)
 allocate(work(lwork))

 	call DSYEV('V','L', nn, h_s , nn, tsc, WORK, lwork, info )
 
  vsc=h_s

 	
 	!do i=1,3
 	!if(i==2)cycle
 	!write(31,*)i,tsc(i)
 	do j=1,nn
 	write(31,*)h_s(j,2)
 	end do
 	
 	!end do
 
 
 !! Calculates dipole moment !!
 
 
 print*,"Calculating absorption from ground state"
 write(19,*)"1Bu state"
 write(19,*)"t_chain=",tb,"eps(F)=",eps(4),"eps(ROOC)=",eps(5),"eps(OR)=",eps(14)
 
 
 
 !state1=0
 tdip=0.0
 do ii=1,500
 !print*,"state",ii
 
 sumx=0.0
 sumy=0.0
 sumz=0.0
 
 do i=1,n
 
 rho=0.0
 
 do lv1=1,n/2+sulp/2
 do lv2=n/2+sulp/2+1,n
 
 
 
 
  
  
  !rho=rho+(vtmp(i,lv1)*vtmp(i,lv2))*vsc(s_index(lv1,lv2),ii)
  rho=rho+sqrt(2.0)*(vtmp(i,lv1)*vtmp(i,lv2))*vsc(s_index(lv1,lv2),ii)
  
 end do
 end do
  !print*,"rho=",rho
  sumx=sumx+rho*x_pos(i)
  sumy=sumy+rho*y_pos(i)
  sumz=sumz+rho*z_pos(i)
  
  !print*,"Dipole_x",dipole_x(state1)
  end do
  !write(15,*)"state",s_index(lv1,lv2),lv1,lv2
 
 dip=sqrt(sumx*sumx+sumy*sumy+sumz*sumz)
 if(dip>tdip)then
 tdip=dip
 tstate=ii
 end if
 write(19,*),ii,tsc(ii),dip!,"x-dipole",sumx,"y-dipole",sumy
 end do
 
      do i=1,nn
      del(i)=tsc(i)-tsc(1)
      end do
 
 
 
 print*,"1Bu=",tstate,tsc(tstate)
 
 !stop	
! Calculating excited states absorption
 print*,"Calculating excited states absorption"
!

 write(22,*)"mAg states"
 write(22,*)"t_chain=",tb,"eps(F)=",eps(4),"eps(ROOC)=",eps(5),"eps(OR)=",eps(14)
 
 do ii=1,500
!
 sum_x=0.d0
 sum_y=0.d0
 sum_z=0.d0
!
 do mu=1,n
!
  hsum = 0.d0
!
 do lv1=1,n/2+sulp/2 
 do lv2=n/2+sulp/2+1,n
!
 do lv3=1,n/2+sulp/2	 
!
   if(lv3/=lv1) &
 &  hsum = hsum - vsc(s_index(lv3,lv2),ii)*vsc(s_index(lv1,lv2),tstate)* &
 &  vtmp(mu,lv3)*vtmp(mu,lv1)
!
 end do
!
 do lv4=n/2+sulp/2+1,n

   if(lv4/=lv2) &
 &  hsum = hsum + vsc(s_index(lv1,lv4),ii)*vsc(s_index(lv1,lv2),tstate)* &
 &  vtmp(mu,lv2)*vtmp(mu,lv4)
!
 end do
 
 hsum1=0.0
 do i=1,n/2+sulp/2
 	if((i.ne.lv1).or.(i.ne.lv2))then
 	hsum=hsum + 2.0*vtmp(mu,i)*vtmp(mu,i)*vsc(s_index(lv1,lv2),ii)*vsc(s_index(lv1,lv2),tstate)
 	end if
 end do
 
 hsum=hsum + vsc(s_index(lv1,lv2),ii)*vsc(s_index(lv1,lv2),tstate)* &
 	(vtmp(mu,lv2)**2-vtmp(mu,lv1)**2) !- 2.0*vtmp(mu,lv1)**2- 2.0*vtmp(mu,lv2)**2
 
!
end do
end do
 	sum_x= sum_x + hsum*x_pos(mu)
 	sum_y= sum_y + hsum*y_pos(mu)
 	sum_z= sum_z + hsum*z_pos(mu)
! 
     end do
        dip_mag(ii)=sqrt(sum_x*sum_x+sum_y*sum_y+sum_z*sum_z)

	if(dip_mag(ii)>0.1)then
         write(21,*)ii,dip_mag(ii)
         write(22,*)tsc(ii),dip_mag(ii)
        end if
         write(10,*)ii,tsc(ii),"Dipole=",dip_mag(ii)

      end do
      
 
      
      print*,"Calculating absorption:"
      
      
      omega=0.0
      do i=1,800
      alpha=0.0
      do ii=2,300
      alpha=alpha+eta*omega*(dip_mag(ii))**2/((omega-del(ii))**2+eta**2)
      end do
      
      write(45,*)omega,alpha
      omega=omega+0.005
      end do
      
! 	
  print*,"1Bu=",tstate,"Energy=",tsc(tstate)

do i=1,nn
!write(28,*)i,E_tot(i)
end do
   
end program ppv
 