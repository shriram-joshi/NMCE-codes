!-*- mode: f90;-*-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!FullGaussSolve is program used to solve full           !
!a matrix equation of the form:                         !
!                                                       !
!JJ*u=FF                                                !
!                                                       !
! where JJ and FF are the Jacobian and RHS matrix       !
! sent into the program with dimensions (n,n+1) and     !
! n, respectively.  The program solves the matrices     !
! by full gaussian elimination and the resulting        !
! solution is stored and passed out as SolutionMat(FF)  !
!                                                       ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Subroutine FullGaussSolverp(JJ,FF,Nodes)
        use omp_lib
        !Implicit Integer (i-n)
        !Implicit Real*8 (a-h,o-z)
        implicit none
        integer*4 :: i, j, m, imaxpos, imaxposl, Nodes
        Real*8 :: rowhold(1,Nodes), entrymax, c, b
        Real*8 :: JJ(Nodes,Nodes+1), FF(Nodes)



        !Append the response matrix to jacobian

        Do i = 1,Nodes
           JJ(i,Nodes+1) = FF(i)
        End Do

        !call progress_gauss(0.0d0) 
        Do j = 1,Nodes       !j = p Nodes = s

           !Find the largest entry in the row, the pivot
           !$omp parallel shared(JJ,j,Nodes,imaxpos,c) private(i,imaxposl,entrymax)
           entrymax = 0.0D0
           c = 0.0d0
           !$omp do schedule(dynamic)
           Do i = j,Nodes
              If (ABS(JJ(i,j)).gt.entrymax) then
                 entrymax = ABS(JJ(i,j))
                 imaxposl = i
              End If
           End Do
           !$omp end do
           !$omp critical
           if (ABS(entrymax).gt.c) then
              c = entrymax
              imaxpos = imaxposl

           end if
           !$omp end critical
           !$omp end parallel

           If (imaxpos.ne.j) then
              !$omp parallel shared(JJ,Nodes,imaxpos,j) private(c,m)
              !$omp do schedule(dynamic)
              Do m = 1,Nodes+1
                 c = JJ(j,m)

                 JJ(j,m) = JJ(imaxpos,m)
                 JJ(imaxpos,m) = c

              End Do
              !$omp end do
              !$omp end parallel
           End If



           !Divide the pivot row by the current entry value
           b = JJ(j,j) 

           !$omp parallel shared(JJ,Nodes,b,j) private(i,m,c)
           !$omp do schedule(dynamic)
           Do m = j,Nodes+1
              JJ(j,m) = JJ(j,m)/b
           End Do
           !$omp end do

           !Do Elimination
           !$omp do schedule(dynamic)
           Do i = j+1,Nodes
              c = JJ(i,j)
              Do m = j,Nodes+1
                 JJ(i,m) = JJ(i,m) - c*JJ(j,m)
              End Do
           End Do
           !$omp end do
           !$omp end parallel

           if (mod(j,10).eq.0) then
              !call progress_gauss(DBLE(j)/DBLE(Nodes))
           end if
        End Do
        !call progress_gauss(1.0d0)
        !call progress_back(0.0d0) 
        Do j = Nodes,1,-1

           !$omp parallel shared(j,JJ,Nodes) private(i,m,c)
           !$omp do schedule(dynamic)
           Do i = (j-1),1,-1
              c = JJ(i,j)

              Do m = j,Nodes+1
                 JJ(i,m) = JJ(i,m)-c*JJ(j,m)
              End Do

           End Do
           !$omp end do
           !$omp end parallel
           ! if (mod(j,10).eq.0) then
           !    call progress_back(DBLE(Nodes-j)/DBLE(Nodes))
           ! end if
        end do
        !call progress_back(1.0d0) 
        !$omp parallel shared(Nodes,JJ,FF) private(i)
        !$omp do schedule(dynamic)
        Do i = 1,Nodes
           FF(i) = JJ(i,Nodes+1)
        End Do
        !$omp end do
        !$omp end parallel
        Return
      End Subroutine FullGaussSolverp
