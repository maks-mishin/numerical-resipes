! Author : Maks Mishin
! Date   : 1/27/2022

module ModTM
    use ModGaussComplex
    implicit none

    real(8),parameter::IE_pi=4.d0*atan(1.d0)
    real(8),parameter::IE_eps=1.d0/(36.d0*IE_pi)*1.d-9
    real(8),parameter::IE_mu=4.d0*IE_pi*1.d-7
    real(8),parameter::IE_k=2.d0*IE_pi
    real(8),parameter::IE_omega=IE_k/sqrt(IE_eps*IE_mu)
    real(8),parameter::IE_gamma=1.781d0
    contains

    complex function Hankel(x) result(res)
        real(8),intent(in)::x
        real(8)::bj,by,f0,t0

        if (x<=3.d0) then
            bj = (x/3.d0)**2
            bj = 1.d0+bj*(-2.2499997d0+bj*(1.2656208d0+bj*(-0.3163866d0+bj*(0.0444479d0+bj*(-0.0039444d0+bj*0.00021d0)))))
            by = (x/3.d0)**2
            by = 2.d0/3.1415926d0*log(x/2.d0)*bj+0.36746691d0+by*(0.60559366d0+by*(-0.74350384d0+by*(0.25300117d0+by*(-0.04261214d0+by*(0.00427916d0-by*0.00024846d0)))))
        else
            bj = 3.d0/x
            f0 = 0.79788456d0+bj*(-0.00000077d0+bj*(-0.00552740d0+bj*(-0.00009512d0+bj*(0.00137237d0+bj*(-0.00072805d0+bj*0.00014476d0)))))
            t0 = x-0.78539816d0+bj*(-0.04166397d0+bj*(-0.00003954d0+bj*(0.00262573d0+bj*(-0.00054125d0+bj*(-0.00029333d0+bj*0.00013558d0)))))
            by = SQRT(x)
            bj = f0*COS(t0)/by
            by = f0*SIN(t0)/by
        end if
        res=cmplx(bj,-by)
    end function

    subroutine GetMatrix_TM(nmax,length,z)
        integer,intent(in)::nmax
        real(8),intent(in)::length
        complex,intent(inout)::z(nmax,nmax)
        integer::n,m
        real(8)::xm,xn,re,im
        real(8)::delta,coeff
        complex::cdelta

        delta=length/nmax
        cdelta=cmplx(delta,0.d0)
        do m=1,nmax
            do n=1,nmax
                if (m==n) then
                    z(m,n)=cmplx(delta,-delta*2.d0/IE_pi*(log(IE_gamma*IE_k*delta/4.d0)-1.d0))
                elseif (m-n==1 .or. m-n==-1) then
                    z(m,n)=cmplx(delta,-delta/IE_pi*(3.d0*log(3.d0*IE_gamma*IE_k*delta/4.d0)-log(IE_gamma*IE_k*delta/4.d0)-2.d0))
                else
                    xm=(m-0.5d0)*delta
                    xn=(n-0.5d0)*delta
                    z(m,n)=Hankel(IE_k*abs(xm-xn))*cdelta
                end if
            end do
        end do
    end subroutine

    subroutine GetVector_TM(nmax,length,phi,b)
        integer,intent(in)::nmax
        real(8),intent(in)::length,phi
        complex,intent(inout)::b(nmax)
        integer::n
        real(8)::x,delta

        delta=length/nmax
        do n=1,nmax
            x=(n-0.5d0)*delta
            b(n)=cmplx(4.d0/(IE_omega*IE_mu)*cos(IE_k*x*cos(phi)),4.d0/(IE_omega*IE_mu)*sin(IE_k*x*cos(phi)))
        end do
    end subroutine

    real(8) function CrossSection_TM(nmax,length,phi,j) result(res)
        integer,intent(in)::nmax
        real(8),intent(in)::length,phi
        complex,intent(in)::j(:)
        integer::n
        real(8)::x,delta
        complex::cdelta,sum,e

        delta=length/nmax
        cdelta=cmplx(delta,0.d0)
        sum=(0.d0,0.d0)
        do n=1,nmax
            x=(n-0.5d0)*delta
            e=cmplx(cos(IE_k*x*cos(phi)),sin(IE_k*x*cos(phi)))
            sum=sum+j(n)*cdelta*e
        end do
        res=cabs(sum)**2*IE_k*IE_mu/IE_eps*0.25d0
    end function

    subroutine Solve_TM(nmax,length,phi)
        integer,intent(in)::nmax
        real(8),intent(in)::length,phi
        complex,allocatable::z(:,:),b(:),x(:)
        integer::ierr
        integer::i

        allocate(z(nmax,nmax),b(nmax),x(nmax))
        call GetMatrix_TM(nmax,length,z)
        call GetVector_TM(nmax,length,phi,b)
        call gauss_complex(nmax,z,b,x,ierr)
        deallocate(z,b,x)
    end subroutine

    subroutine SolveRCS_TM(nmax,length)
        integer,intent(in)::nmax
        real(8),intent(in)::length
        complex,allocatable::z(:,:),b(:),x(:)
        integer::ierr
        integer::i
        real(8)::phi

        allocate(z(nmax,nmax),b(nmax),x(nmax))
        do i=1,100
            phi=IE_pi/100.d0*(i-0.5d0)
            call GetMatrix_TM(nmax,length,z)
            call GetVector_TM(nmax,length,phi,b)
            call gauss_complex(nmax,z,b,x,ierr)
            write(200,*)phi,10.d0*log10(CrossSection_TM(nmax,length,phi,x))
        end do
        deallocate(z,b,x)
    end subroutine
end module ModTM