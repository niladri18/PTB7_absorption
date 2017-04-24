module parameter

! define all the global variables

implicit none

        double precision,dimension(:,:),allocatable::arr,f,v_hf,v_n,f0,v_ij,h_s,m_states,vtmp,vsc
        integer,dimension(:,:),allocatable::s_index
        double precision,dimension(:),allocatable::d,charge,p_on,d_hf,x_pos,y_pos,z_pos,E_tot,dipole_x,dipole_y,eps,U_c,z,d_n,dipole_z,tmp,tsc,work,dip_mag,del
        double precision::t,en,smat,smatu,E_p,E_a,phi,theta,psi,r_s,r_b,r_d,rho,hsum,hsum1
        double precision::t2,t3,norm,t1,v_ii,kappa,r_ij,p_off,err,err1,conv,E_triplet,E_singlet,mid,dir,xchange,alpha,eta
        double precision::sumx,sumy,rnorm,sum_x,sum_y,tdip,dip,dip_x,dip_y,epsi,eps1,tcs,tb,t4,sumz,sum_z,dip_z,omega
        !real(SP)::r_s,r_b,r_d,pi
        integer::nrot,n,i,i1,i2,j,k,k1,l,x,xc,it,k2,i3,i4,nn,no_states,site,nbonds1,nbonds2,nbonds3,nbonds4,nbonds5,numst
        integer::j1,j2,j3,j4,nrot1,lv1,lv2,lv3,lv4,state1,state2,r_no,mu,nu,tstate,nbnd,hf,sulp,ii,lwork,info,lwork1
        integer,dimension(:,:),allocatable::ilbond

end module parameter

