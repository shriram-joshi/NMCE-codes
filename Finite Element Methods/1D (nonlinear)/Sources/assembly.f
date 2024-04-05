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
        real(kind=rk):: dxdxi, uL, duLdxi

        do k = 1,nE ! Change according to basis function used

            ! Calculate the local j matrix
            jL = 0.0_rk
            rL = 0.0_rk
            do l = 1,3 ! Sum over quadrature points       
                
                ! Interpolate variables
                dxdxi = 0.0_rk
                uL = 0.0_rk
                duLdxi = 0.0_rk
                do i=1,nLP
                    dxdxi = dxdxi + xMesh(nop(i,k))*dph(i,l) 
                    uL = uL + UG(nop(i,k))*ph(i,l)
                    duLdxi = duLdxi + UG(nop(i,k))*dph(i,l) 
                end do

                do i = 1,nLVar
                    do j = 1,nLVar
                        jL(i,j) = jL(i,j) + wei(l)*(2.0_rk*dxdxi*KK*uL*ph(i,l)*ph(j,l) + dph(i,l)*dph(j,l)/dxdxi)
                    end do
                    rL(i) = rL(i) + wei(l)*(dxdxi*KK*uL*uL*ph(i,l) + dph(i,l)*duLdxi/dxdxi)
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
        RG(1) = 0.0_rk
        RG(nVar) = 0.0_rk

        ! For nonlinear code JJ*delU = -R
        RG = -RG

        JG(1,1:nVar) = 0.0_rk
        JG(nVar,1:nVar) = 0.0_rk
        JG(1,1) = 1.0_rk
        JG(nVar,nVar) = 1.0_rk
    end subroutine global_assembly
    !***********************************************************************************************************************
    ! Global Assembly Ends
    !***********************************************************************************************************************
end module assembly