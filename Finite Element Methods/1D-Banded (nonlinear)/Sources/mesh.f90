module meshing
    
    use kind
    use data
    use omp_lib
    implicit none

contains
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
end module meshing