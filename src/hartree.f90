module hartree

        use parameter
contains 

subroutine hartree_fock
do hf=1,400

        !! self consistent calculation starts
        it=it+1


         ! calculating charge density at site i
        do i=1,n
                charge(i)=0.0
                do k=1,n/2+sulp/2
                        charge(i)=charge(i)+2.0*v_n(i,k)*v_n(i,k)
                end do
        end do
        
        !!!  Assiging off-diagonal elements !!!

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
         
        ! end of assigning the off diagonal elements !




        !!! assigning the diagonal elements !!!

        do i=1,n
                p_on(i)=0.0
                do j=1,n
                        if(i.ne.j)then
                                p_on(i)=p_on(i) + (charge(j)-z(j))*v_ij(i,j)
                                !print*,"Chemical
                                !pot",i,j,z(j),charge(j),v_ij(i,j)
                        end if
                end do
                f(i,i)=p_on(i)+0.5*charge(i)*U_c(i)+eps(i)
        end do

        !!! end assigning the diagonal elements !!!





        call DSYEV('V','U', N, f , N, tmp, WORK, lwork, info )


        
        
        err=0.0
        do i=1,n
                err1=(tmp(i)-d_n(i))**2
                err=err+err1
        end do
        conv=sqrt(err)/n
        v_n =f
        d_n=tmp
        if(conv.le.0.0000001)exit

end do

end subroutine hartree_fock

end module hartree
