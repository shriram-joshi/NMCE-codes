Module library

    use kind
    use Data_
    use omp_lib
    implicit none

contains
!***********************************************************************************************************************
! Hermite Cubic functions
!***********************************************************************************************************************
    subroutine hermite_basis_phi1d()
        integer :: i, j

        ph=0.0_rk
        dph=0.0_rk
        do j=1,3 ! For the three gauss points
            ph(1,j) = 1.0_rk - 3.0_rk*xiarr(j)*xiarr(j) + 2.0_rk*xiarr(j)*xiarr(j)*xiarr(j)
            dph(1,j)= -6.0_rk*xiarr(j) + 6.0_rk*xiarr(j)*xiarr(j)

            ph(2,j) = xiarr(j) - 2*xiarr(j)*xiarr(j) + xiarr(j)*xiarr(j)*xiarr(j)
            dph(2,j)= 1.0_rk - 4.0_rk*xiarr(j) + 3*xiarr(j)*xiarr(j)

            ph(3,j) = 3*xiarr(j)*xiarr(j) - 2*xiarr(j)*xiarr(j)*xiarr(j)
            dph(3,j)= 6.0_rk*xiarr(j) - 6.0_rk*xiarr(j)*xiarr(j)

            ph(4,j) = xiarr(j)*xiarr(j)*xiarr(j) - xiarr(j)*xiarr(j)
            dph(4,j)= 3.0_rk*xiarr(j)*xiarr(j) - 2*xiarr(j)
        end do

        print*, "ph ="
        do i = 1, nl
            do j = 1, nl
                write(*, fmt='(F15.2)', advance='no') ph(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do

        print*, "dph ="
        do i = 1, nl
            do j = 1, nl
                write(*, fmt='(F15.2)', advance='no') dph(i, j)
            end do
            print *  ! Move to the next line after printing each row
        end do
    end subroutine hermite_basis_phi1d
!***********************************************************************************************************************
! Quadratic Basis functions End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Generate Variably Spaced Mesh
!***********************************************************************************************************************
    subroutine generate_varmesh()
            
        implicit none

        integer :: i
        real(kind=rk) :: A, B

        A = ((n+1)*xSpan(2)-xSpan(1))/n
        B = xSpan(1)-A
        ! Mesh generated according x = A+B/i
        do i = 1, n+1
            xMesh(i) = A + B/i
        end do

    end subroutine generate_varmesh
!***********************************************************************************************************************
! Generate Variably Spaced End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Generate Uniform Mesh
!***********************************************************************************************************************
    subroutine generate_unifmesh()
            
        implicit none

        integer :: i
        real(kind=rk) :: A, B

        A = ((n+1)*xSpan(1)-xSpan(2))/n
        B = xSpan(1)-A
        ! Mesh generated according x = A+B/i
        do i = 1, n+1
            xMesh(i) = A + B*i
        end do

    end subroutine generate_unifmesh
!***********************************************************************************************************************
! Generate Uniform Mesh End
!***********************************************************************************************************************

end module library
