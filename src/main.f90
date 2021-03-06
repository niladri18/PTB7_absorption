program ppv

        use parameter
        use readfile
        use hartree
        use make_huckel
        use make_coulomb
        use allocation
        use sci

!open(unit=1,file="dipole_n.inp")
open(unit=10,file="./output/nonlinear.out")

print*,"Calculates the energy levels of PTB7 in the huckel limit"
print*,"then calculates the Hartree-Fock energy"
print*,"then calculates excited sates using SCI." 

print*, "HF calculations starts"

! ==== Open files =====
   open(unit=90,file="./input/ptb7cord4.inp")
   open(unit=91,file="./input/unit4.inp")
   open(unit=25,file="./output/hf_ppt.out")
   read(90,*) n

   sulp=n/5
   !sulp=0
   nn=(n/2+sulp/2)*(n/2-sulp/2)
   lwork=3*n-1
   

! ===== Allocate the arrays depending on the value of n ===== !

        call allocate_array

       !do lv3=n/2-sulp/2+1,n
       !do lv4=n/2-sulp/2,1,-1
       dip_mag=0.0
       ilbond = 0.0


! ====== Reading Main File ====== !

        call read_config
 
 

! ====== Constructing Huckel matrix ====== ! 

        call make_huckel_level


!  ===== Calculating the coulomb interactions ===== !

        call make_coulomb_interaction
        
       

! ====== Hartree-Fock calculation ====== !

         print*,"Hartree-Fock calculation starts"   

        v_n=arr
        d_n=tmp
        tmp=0.0
        vtmp=0.0

        it=0


        call hartree_fock

        !call eigsrt(d_hf,v_hf)

        print*,"converged in steps",it

        do i=1,n
                write(25,*)tmp(i),i,charge(i)
                !print*,tmp(i),i,charge(i)
                !print*,d_hf(i),i,charge(i)
        end do


print*, "HF calculations ends"



! ========== END HF CALCULATION  ==============
! ============================================= 
!====== Starting SCI calculations ======!!!!

        call sci_model
 
 
 
 
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
 
