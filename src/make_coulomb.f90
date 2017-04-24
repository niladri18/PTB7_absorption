module make_coulomb
        use parameter

contains 

subroutine make_coulomb_interaction

!  ===== Calculating the coulomb interactions ===== !
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


! ===== end calculating the coulomb interactions ===== !



end subroutine make_coulomb_interaction

end module make_coulomb
