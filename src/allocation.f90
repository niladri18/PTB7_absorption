module allocation

        use parameter

contains

subroutine allocate_array

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
        allocate(ilbond(2,200))


end subroutine allocate_array

end module allocation

