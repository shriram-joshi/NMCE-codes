Module library

    use kind
    use Data_
    use omp_lib
    implicit none

contains
!***********************************************************************************************************************
! Linear Basis functions
!***********************************************************************************************************************
    subroutine lin_basis_phi1d()
        integer :: j

        ! if(3 /= size(ph1)/size(ph1(:,6)))then
        !     print*, (size(ph1)/size(ph1(:,6)))
        !     print*, "Cant use linear basis function subroutine. Change nl in data.f to match the required type of basis functions."
        !     stop
        ! end if

        ph1=0.0_rk
        dph1=0.0_rk
        do j=1,3
            ph1(1,j) = 1.0_rk - xiarr(j)
            dph1(1,j)= -1.0_rk
            ph1(2,j) = xiarr(j)
            dph1(2,j)= 1.0_rk
        end do
    end subroutine lin_basis_phi1d
!***********************************************************************************************************************
! Linear Basis functions End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Quadratic Basis functions
!***********************************************************************************************************************
    subroutine quad_basis_phi1d()
        integer :: j

        ! if(3 /= size(ph1)/size(ph1(:,6)))then
        !     print*, size(ph1)/size(ph1(:,6))
        !     print*, "Cant use quadratic basis function subroutine. Change nl in data.f to match the required type of basis functions."
        !     stop
        ! end if

        ph1=0.0_rk
        dph1=0.0_rk
        do j=1,3
            ph1(1,j) = 1.0_rk - 3.0_rk*xiarr(j) + 2.0_rk*xiarr(j)*xiarr(j)
            dph1(1,j)= -3.0_rk + 4.0_rk*xiarr(j)
            ph1(2,j) = 4.0_rk*(xiarr(j) - xiarr(j)*xiarr(j))
            dph1(2,j)= 4.0_rk - 8.0_rk*xiarr(j)
            ph1(3,j) = - xiarr(j) + 2.0_rk*xiarr(j)*xiarr(j)
            dph1(3,j)= -1.0_rk+4.0_rk*xiarr(j)
        end do
    end subroutine quad_basis_phi1d
!***********************************************************************************************************************
! Quadratic Basis functions End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Generate Variably Spaced Mesh
!***********************************************************************************************************************
    subroutine generate_mesh()
            
        implicit none

        integer :: i
        real(kind=rk) :: A, B

        A = ((n+1)*xSpan(2)-xSpan(1))/n
        B = xSpan(1)-A
        ! Mesh generated according x = A+B/i
        do i = 1, n+1
            xMesh(i) = A + B/i
        end do

        ! Uncomment to print out mesh points
        print*, "Variable mesh points "
        do i = 1, n+1
            write(*, fmt='(F15.4)', advance='no') xMesh(i)
        end do

    end subroutine generate_mesh
!***********************************************************************************************************************
! Generate Variably Spaced End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Generate Uniform Mesh
!***********************************************************************************************************************
    subroutine generate_uniform_mesh()
            
        implicit none

        integer :: i
        real(kind=rk) :: A, B

        A = ((n+1)*xSpan(1)-xSpan(2))/n
        B = xSpan(1)-A
        ! Mesh generated according x = A+B/i
        do i = 1, n+1
            xMesh(i) = A + B*i
        end do

        ! Uncomment to print out mesh points
        print*, "Uniform mesh points"
        do i = 1, n+1
            write(*, fmt='(F15.4)', advance='no') xMesh(i)
        end do

    end subroutine generate_uniform_mesh
!***********************************************************************************************************************
! Generate Uniform Mesh End
!***********************************************************************************************************************

end module library
