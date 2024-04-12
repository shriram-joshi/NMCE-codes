module assembly
    use kind
    use data
    use library
    implicit none
    
contains
    
    subroutine global_assembly()
        implicit none

        integer:: k

        do k = 1,nE ! Change according to basis function used

            call local_assembly(k)

            ! Dump local j into the global J matrix
            call dump(k)

        end do

        call boundary_conditions()
    end subroutine global_assembly
    
    subroutine local_assembly(n)
        implicit none
        
        ! n = global element number 
        integer, intent(in) :: n

        integer:: l, i, j

        real(kind=rk):: dxdxi, uL, duLdxi
        
        ! Calculate the local j matrix
        jL = 0.0_rk
        rL = 0.0_rk
        do l = 1,3 ! Sum over quadrature points       
            
            ! Interpolate variables
            dxdxi = 0.0_rk
            uL = 0.0_rk
            duLdxi = 0.0_rk
            do i=1,nLP
                dxdxi = dxdxi + xMesh(nop(i,n))*dph(i,l) 
                uL = uL + UG(nop(i,n))*ph(i,l)
                duLdxi = duLdxi + UG(nop(i,n))*dph(i,l) 
            end do

            do i = 1,nLVar
                do j = 1,nLVar
                    jL(i,j) = jL(i,j) + wei(l)*(dph(i,l)*dph(j,l)/dxdxi + Pe*ph(i,l)*dph(j,l))
                end do
                rL(i) = rL(i) + wei(l)*(duLdxi*dph(i,l)/dxdxi + Pe*duLdxi*ph(i,l))
            end do
        end do
        
    end subroutine local_assembly
    
    subroutine dump(n)
        implicit none
        integer :: n, i, j, k, l, o, p
        
        ! Dump local j into the global J matrix
        do i = 1, nLVar
            ! k loop varies over eqautions for different variables at same node
            do k = 0, mdf(nop(i,n))-1
                o = nopp(nop(i,n)) + k
                do j = 1, nLVar
                    ! l loop varies over different variables at same node
                    do l = 0, mdf(nop(j,n))-1
                        p = nopp(nop(j,n)) + l
                        JG(bgbl(o,p,1),bgbl(o,p,2)) = JG(bgbl(o,p,1),bgbl(o,p,2)) + jL(i,j)
                    end do
                end do             
            end do
            RG(nop(i,n)) = RG(nop(i,n)) + rL(i)
        end do
    end subroutine dump

    subroutine boundary_conditions()
        implicit none
        
        ! Changing J and R according to the BCs
        RG(1) = 0.0_rk
        RG(nVar) = 0.0_rk

        ! For nonlinear code JJ*delU = -R
        RG = -RG

        JG(1,:) = 0.0_rk
        JG(bgbl(1,1,1),bgbl(1,1,2)) = 1.0_rk

        JG(nVar, :) = 0.0_rk
        JG(bgbl(nVar, nVar, 1), bgbl(nVar, nVar,2)) = 1.0_rk
        
    end subroutine boundary_conditions
    
end module assembly