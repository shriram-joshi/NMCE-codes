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
        integer :: i

        ph1=0.0_rk
        dph1=0.0_rk
        do i=1,3
            ph1(1,i) = 1.0_rk - xiarr(i)
            dph1(1,i)= -1.0_rk
            ph1(2,i) = xiarr(i)
            dph1(2,i)= 1.0_rk
        end do
    end subroutine lin_basis_phi1d
!***********************************************************************************************************************
! Linear Basis functions End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Quadratic Basis functions
!***********************************************************************************************************************
    subroutine quad_basis_phi1d()
        implicit none
        integer :: i

        ph1=0.0_rk
        dph1=0.0_rk
        do i=1,3
            ph1(1,i) = 1.0_rk - 3.0_rk*xiarr(i) + 2.0_rk*xiarr(i)*xiarr(i)
            dph1(1,i)= -3.0_rk + 4.0_rk*xiarr(i)
            ph1(2,i) = 4.0_rk*(xiarr(i) - xiarr(i)*xiarr(i))
            dph1(2,i)= 4.0_rk - 8.0_rk*xiarr(i)
            ph1(3,i) = - xiarr(i) + 2.0_rk*xiarr(i)*xiarr(i)
            dph1(3,i)= -1.0_rk+4.0_rk*xiarr(i)
        end do
    end subroutine quad_basis_phi1d
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

        A = ((nMesh)*xSpan(2)-xSpan(1))/(nMesh-1)
        B = xSpan(1)-A

        ! Mesh generated according x = A+B/i
        do i = 1, nMesh
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

        A = ((nMesh)*xSpan(1)-xSpan(2))/(nMesh-1)
        B = xSpan(1)-A

        ! Mesh generated according x = A+B/i
        do i = 1, nMesh
            xMesh(i) = A + B*i
        end do

    end subroutine generate_unifmesh
!***********************************************************************************************************************
! Generate Uniform Mesh End
!***********************************************************************************************************************

!***********************************************************************************************************************
! Global Assembly 
!***********************************************************************************************************************
    subroutine global_assembly()
        implicit none

        integer:: i, j, k, l
        real(kind=rk):: dxdxi

        do k = 1,nE ! Change according to basis function used

            ! Calculate the local j matrix
            jL = 0
            rL = 0
            do l = 1,3 ! Sum over quadrature points       
                
                ! Calculate dx/dxi
                dxdxi = 0
                do i=1,nLP
                    dxdxi = dxdxi + xMesh(nop(i,k))*dph1(i,l) 
                end do

                do i = 1,nLVar
                    do j = 1,nLVar
                        jL(i,j) = jL(i,j) + wei(l)*(dph1(i,l)*dph1(j,l)/dxdxi + Pe*ph1(i,l)*dph1(j,l))
                    end do
                    rL(i) = 0
                end do
            end do
    
            ! Dump local j into the global J matrix
            do i = 1, nLP
                do j = 1, nLP
                    JG(nop(i,k),nop(j,k)) = JG(nop(i,k),nop(j,k)) + jL(i,j)
                end do
                RG(nop(i,k)) = RG(nop(i,k)) + rL(i)
            end do
        end do

        ! Changing J and R according to the BCs
        RG(1) = uBC(1)
        RG(nVar) = uBC(2)

        JG(1,1:nVar) = 0.0_rk
        JG(nVar,1:nVar) = 0.0_rk
        JG(1,1) = 1.0_rk
        JG(nVar,nVar) = 1.0_rk
    end subroutine global_assembly
!***********************************************************************************************************************
! Global Assembly Ends
!***********************************************************************************************************************

end module library
