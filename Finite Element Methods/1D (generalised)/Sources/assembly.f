module assembly

    use kind
    use data
    implicit none
    
contains
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
end module assembly