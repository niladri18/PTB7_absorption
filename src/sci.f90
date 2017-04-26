module sci

        use parameter


contains 

subroutine sci_model 
 
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

end subroutine sci_model


!! Next Calculates dipole moment !!

end module sci
 
