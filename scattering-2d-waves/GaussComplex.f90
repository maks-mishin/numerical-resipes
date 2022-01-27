! Author : Maks Mishin
! Date   : 1/27/2022

module ModGaussComplex
    implicit none
    contains

    subroutine gauss_complex(n,a,b,x,ierr)
        integer,intent(in)::n
        complex,intent(inout)::a(n,n),b(n)
        complex,intent(out)::x(n)
        integer,intent(inout)::ierr
        integer::i,j,m
        complex::l,sum
        real(8)::eps

        eps=1.d-15
        ierr=0
        do i=1,n
            m=maxnum(a,n,i,i) !!!
            if (cabs(a(m,i))>eps) then
                call swap(a,b,n,i,m) !!!
            else
                ierr=1
                return
            end if
            do j=i+1,n
                l=-a(j,i)/a(i,i)
                call elem(a,b,n,i,j,l) !!!
            end do
        end do
        do i=n,1,-1
            sum=(0.d0,0.d0)
            do j=i+1,n
                sum=sum+a(i,j)*x(j)
            end do
            x(i)=(b(i)-sum)/a(i,i)
        end do
    end subroutine

    integer function maxnum(a,n,i,l) result(res)
        integer,intent(in)::n,i,l
        complex,intent(in)::a(n,n)
        integer::j,ind
        real(8)::max

        max=cabs(a(l,i))
        ind=l
        do j=l,n
            if (cabs(a(j,i))>max) then
                max=cabs(a(j,i))
                ind=j
            end if
        end do
        res=ind
    end function

    subroutine swap(a,b,n,i,m)
        integer,intent(in)::n,i,m
        complex,intent(inout)::a(n,n),b(n)
        integer::j
        complex::x

        do j=1,n
            x=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=x
        end do
        x=b(i)
        b(i)=b(m)
        b(m)=x
    end subroutine

    subroutine elem(a,b,n,i,j,l)
        integer,intent(in)::n,i,j
        complex,intent(inout)::a(n,n),b(n)
        complex,intent(in)::l
        integer::k

        do k=1,n
            a(j,k)=a(j,k)+l*a(i,k)
        end do
        b(j)=b(j)+l*b(i)
    end subroutine

    subroutine swaprow(a,n,i,m)
        integer,intent(in)::n,i,m
        real(8),intent(inout)::a(n,n)
        integer::j
        real(8)::x

        do j=1,n
            x=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=x
        end do
    end subroutine

    subroutine elemrow(a,n,i,j,l)
        integer,intent(in)::n,i,j
        real(8),intent(inout)::a(n,n)
        real(8),intent(in)::l
        integer::k

        do k=1,n
            a(j,k)=a(j,k)+l*a(i,k)
        end do
    end subroutine

end module ModGaussComplex