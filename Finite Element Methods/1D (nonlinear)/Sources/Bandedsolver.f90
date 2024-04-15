!***********************************************************************************************************************
! The Banded Solver
!***********************************************************************************************************************
!                                                                           !
! 	BANDEDSOLVER	Solves a system of equations that forms a                   !
! 			banded matrix; so that the diagonal lies in                     !
! 	the center of the banded matrix and the band width gives the            !
! 	sum of the maximum entries to the left and right of the                 !
! 	diagonal (not necessarly on the same row).                              !
!                                                                           !
!       Varaible meanings                                                   !
!    NUMVAR - # of equations                                                !
!	 Bandwidth - bandwidth                                                  !
!	 BMAT - global coefficient matrix                                       !
!	 sRHS - global residue                                                  !
!	 Take each band and straighten them and place into columns              !
!	 Add (3*bandwidth - 1)/2 columns of zeros to form sAMAT                 !
!                                                                           !
!       Original code                                                       !
! 	PROGRAMMER:    CHRIS POMMER                                             !
! 	DATE OF USE:   02/23/06                                                 !
!       LAST MODIFIED: 02/23/06                                             !
!                                                                           !
!       Parallelized code                                                   !
!       PROGRAMMER:    CHRIS ANTHONY                                        !
! 	DATE OF USE:   03/14/13                                                 !
!       LAST MODIFIED: 03/14/13                                             !
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BANDEDSOLVER(NUMVAR, Bandwidth, BMAT, sRHS)
  use OMP_LIB
  Implicit None

  Integer :: NUMVAR, Bandwidth, iDIAG
  Integer :: i, j, k, iROW, iROWPOS, iROWPOSL, iSTOP, iCOL
  Integer :: i_1, i_2, i_3, i_4, threads = 1

  Double Precision :: c_MAX, c_SAVE, c_MULT, c_MAXL
  Double Precision :: BMAT(NUMVAR,(3*Bandwidth-1)/2), sRHS(NUMVAR)


  call omp_set_num_threads(threads)

  !CONCATENATE LHS AND RHS MATRICES

  !$omp parallel shared(NUMVAR,BMAT,Bandwidth) private(i,j)
  !$omp do
  DO i = 1,NUMVAR
     BMAT(i,Bandwidth+1) = sRHS(i)
  END DO
  !$omp end do
  !$omp end parallel


  !SOLVE SYSTEM OF EQUATIONS USING GAUSSIAN-ELIMINATION
  iDIAG = (Bandwidth-1)/2

  !$omp parallel shared(BMAT,iDIAG) private(i,k,j)
  !$omp do
  DO i = 1,iDIAG

     DO j = iDIAG-i+2,Bandwidth
        k = j + (-1)*(iDIAG-i+1)
        BMAT(i,k) = BMAT(i,j)
        BMAT(i,j) = 0.0D0
     END DO
  END DO
  !$omp end do
  !$omp end parallel



     DO iROW = 1,NUMVAR

        IF ((iROW+iDIAG).GT.NUMVAR) THEN
           iSTOP = NUMVAR
        ELSE
           iSTOP = iROW+iDIAG
        END IF

        !FIND THE PIVOT

        !$omp parallel shared(c_MAX,BMAT,iSTOP) private(c_MAXL,i_1,iROWPOSL)
        c_MAX = 0.0D0
        c_MAXL = 0.0D0
        !$omp do
        DO  i_1 = iROW,iSTOP
           IF (DABS(BMAT(i_1,1)).GE.DABS(c_MAXL)) THEN
              c_MAXL = BMAT(i_1,1)
              iROWPOSL = i_1
           END IF
        END DO
        !$omp end do
        !$omp critical
        if (DABS(c_MAXL).gt.DABS(c_MAX)) then
           c_MAX = c_MAXL
           iROWPOS = iROWPOSL

        end if
        !$omp end critical
        !$omp end parallel


        !SWAP ROWS IF NECESSARY AND DIVIDE BY PIVOT
        IF (iROWPOS.NE.iROW) THEN
           !$omp parallel shared(BMAT,Bandwidth,iROWPOS,iROW,c_MAX) private(c_SAVE,i_2)
           !$omp do
           DO i_2 = 1,Bandwidth+1
              c_SAVE = BMAT(iROWPOS,i_2)
              BMAT(iROWPOS,i_2) = BMAT(iROW,i_2)
              BMAT(iROW,i_2) = c_SAVE
              BMAT(iROW,i_2) = BMAT(iROW,i_2)/c_MAX
           END DO
           !$omp end do
           !$omp end parallel
        ELSE
           !$omp parallel shared(c_MAX,BMAT,iROW,Bandwidth) private(i_2)
           !$omp do
           DO i_2 = 1,Bandwidth+1
              BMAT(iROW,i_2) = BMAT(iROW,i_2)/c_MAX
           END DO
           !$omp end do
           !$omp end parallel
        END IF




        !$omp parallel shared(BMAT,Bandwidth,iROW,c_MAX) private(i,j,i_3,i_4,c_MULT,k)
        !PERFORM FORWARD ELIMINATION
        !$omp do
        DO i_3 = iROW+1,iSTOP
           c_MULT = BMAT(i_3,1)
           DO i_4 = 1,Bandwidth+1
              BMAT(i_3,i_4) = BMAT(i_3,i_4) - c_MULT*BMAT(iROW,i_4)
           END DO
        END DO
        !$omp end do

        !SHIFT ELIMINATED ROWS
        !$omp do
        DO i = iROW+1,iSTOP
           k = 1
           DO j = 2,Bandwidth
              BMAT(i,k) = BMAT(i,j)
              BMAT(i,j) = 0.0D0
              k = k + 1
           END DO
        END DO
        !$omp end do
        !$omp end parallel



     END DO

     !PERFORM BACK SUBSTITUTION
     DO i_1 = (NUMVAR-1),1,-1
        IF (i_1.LT.Bandwidth) THEN
           iSTOP = i_1+1
        ELSE
           iSTOP = Bandwidth
        END IF
        !$omp parallel shared(BMAT,Bandwidth,iSTOP,i_1) private(iROW,iCOL,c_MULT)
        !$omp do
        DO iCOL = 2,iSTOP
           iROW = i_1-iCOL+2
           c_MULT = BMAT(iROW,iCOL)
           BMAT(iROW,Bandwidth+1) = BMAT(iROW,Bandwidth+1) - c_MULT*BMAT(i_1+1,Bandwidth+1)

        END DO
        !$omp end do
        !$omp end parallel

     END DO


     !UPDATE SOLUTION
     !$omp parallel shared(BMAT,Bandwidth,sRHS) private(i)
     !$omp do
     DO i = 1,NUMVAR
        sRHS(i)= BMAT(i,Bandwidth+1)
     END DO
     !$omp end do
     !$omp end parallel

END SUBROUTINE BANDEDSOLVER
!***********************************************************************************************************************
! The Banded Solver END
!***********************************************************************************************************************

!!      ! Pseudo code for global assembly: we solve Ax=b at each Newton iteration
!!      ! Zero out A and b
!!		vector_b=0.0_rk; Matrix_A=0.0_rk
!!
!!      ! Bandwidth local to global (Use this when you do global assembly)
!!      integer function bgbl(o,p,ij,iBw)
!!          implicit none
!!          integer, intent(in) :: o, ij
!!          integer, intent(in) :: p
!!          integer :: iBw
!!             if (ij==1) then
!!                  bgbl = o
!!             elseif (ij==2) then
!!                  bgbl = p - o + (iBW - 1)/2 + 1
!!             end if
!!          return
!!      end function bgbl
!!      iBw := bandwidth
!!      o := NOPP(NOP(n,i)) + k
!!      p := NOPP(NOP(n,j)) + z
!!      , where i,j,k,z are integers that ranges from
!!      i=[1,2,3], j=[1,2,3], k=[0,1,...,MDF(NOP(n,i))-1], z=[0,1,...,MDF(NOP(n,j))-1]
!!
!!		! Do Gaussian integration here
!!      ! (The subroutine below does both local and global assembly )
!!		do n=1,ne (n:just a local integer, ne:total number of elements)
!!			call gauss_int_Bridge(n)
!!		end do
!!
!!      ! Now that I have constructed A and b, I am now in a poisition to impose bc
!!
!!      ! (h,omega,v) <-- local ordering of variables
!!      ! In what follows, the last argument of bgbl is suppressed for notational simplicity
!!      ! nn := total number of nodes 3*nn stands for total number of unknowns
!!
!!		! At z=0: h=1 for Bridge dhdz=0 for Jet
!!		Matrix_A(1,:) = 0.0_rk; vector_b(1)=0.0_rk;
!!		Matrix_A(Bgbl(1,1,1),Bgbl(1,1,2)) = 1.0_rk;
!!		0vector_b(1) = solution(1) - 0.0_rk
!!
!!		! z=0: v=0
!!		Matrix_A(2,:) = 0.0_rk; vector_b(2)=0.0_rk;
!!		Matrix_A(Bgbl(2,2,1),Bgbl(2,2,2)) = 1.0_rk;
!!		vector_b(2) = solution(2) - 0.0_rk
!!
!!		! z=L: h=1 for Bridge dhdz=0 for Jet
!!		Matrix_A(3*nn-2,:) = 0.0_rk; vector_b(3*nn-2)=0.0_rk;
!!		Matrix_A(Bgbl(3*nn-2,3*nn-2,1),Bgbl(3*nn-2,3*nn-2,2)) = 1.0_rk;
!!		vector_b(3*nn-2) = solution(3*nn-2) - 0.0_rk
!!
!!		! z=L: v=uf (uf can be non zero for stretching bridge; uf=0 for jet)
!!		Matrix_A(3*nn-1,:) = 0.0_rk; vector_b(3*nn-1)=0.0_rk;
!!		Matrix_A(Bgbl(3*nn-1,3*nn-1,1),Bgbl(3*nn-1,3*nn-1,2)) = 1.0_rk;
!!		vector_b(3*nn-1) = solution(3*nn-1) - uf
!!
!!		! Call Bandedsolver (no constraints)
!!		! iBw: Band width
!!		! The solution will be returned to the vector, vector_b
!!		vector_b = -1.0_rk*vector_b
!!		call BANDEDSOLVER(nu, iBw, Matrix_A, vector_b)
!!		delu = vector_b
!!
!!      solution = solution + delu
