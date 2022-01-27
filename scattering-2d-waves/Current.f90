! Author : Maks Mishin
! Date   : 1/27/2022

module Current
    use ModGaussComplex
    use ModTM

    implicit none
    contains

    subroutine CalcCurrent()
        integer, parameter::nmax=100
        real(8),parameter::length=3
        complex,allocatable::z(:,:),b(:),x(:)
        integer::ierr
        integer::i
        real(8)::phi

        open(200, file = "TM_RCS.txt")

        allocate(z(nmax,nmax),b(nmax),x(nmax))
        do i=1,100
            phi=IE_pi/100.d0*(i-0.5d0)
            call GetMatrix_TM(nmax,length,z)
            call GetVector_TM(nmax,length,phi,b)
            call gauss_complex(nmax,z,b,x,ierr)
            write(200,*)phi,10.d0*log10(CrossSection_TM(nmax,length,phi,x))
        end do
        deallocate(z,b,x)
    end subroutine CalcCurrent
end module Current