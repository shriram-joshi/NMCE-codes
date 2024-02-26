SUBROUTINE ArrowSolver(NOD, IBAND, IARROW, NRHS, TAU)
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !Use your_data_module_name, only: AMAT=>your_amat, ROWMAT=>your_rowmat, &
 !                                 COLMAT=>your_colmat, HEAD=>your_head, RHS=>your_rhs
 !you have to define your_amat, your_rowmat, your_colmat, your_head, your_rhs in your data module
 !see below for the information on the size and the rank for these arrays    
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 !************************************************************************************************
 
 use Data_
 Implicit None
 
 INTEGER, INTENT(IN) :: NOD, IBAND, IARROW, NRHS
 DOUBLE PRECISION, INTENT(IN) ::  TAU

 INTEGER :: ISIGN, NPIV, MID, MEX, NEXX
 INTEGER :: KKS, LLS, MS, IS1, JS1, KS1, JS2
 INTEGER :: IS,JS,N
 DOUBLE PRECISION :: DETLOG, BIG, TRANS, DIAG, PIV
 
 
 !-----------------------------------------------------------------------|
 !       NOD -           DIMENSION OF THE BANDED MATRIX A AND OF THE     |
 !                       COLUMNS AND ROWS OF THE ARROW                   |
 !       IBAND -         BAND WIDTH OF MATRIX A, MUST BE ODD             |
 !       IARROW -        WIDTH OF THE COLUMNS AND ROWS IN THE ARROW      |
 !       NRHS -          NUMBER OF B VECTORS                             |
 !                                                                       |
 !       TAU -           MINIMUM RATIO OF DIAGONAL TO MAXIMUM ELEMENT    |
 !                       FOR NO PIVOTING AND  L U, SET TAU=0.0           |
 !                       FOR FULL PIVOTING, SET TAU=1.0                  |
 !                       FOR SOLVING  L U X = B, SET TAU=-1.0            |
 !       DETLOG -        LOG 10 OF THE DETERMINANT                       |
 !       ISIGN -         SIGN OF THE DETERMINANT                         |
 !       NPIV -          NUMBER OF TIMES ROWS HAVE BEEN PIVOTED          |
 !                                                                       |
 !       RHS(NOD+IARROW,NRHS) -  HOLDS B VECTORS.  THE SOLUTION          |
 !                               VECTOR X IS RETURNED IN THIS VECTOR     |
 !       HEAD(IARROW,IARROW) -   HOLDS THE ELEMENTS IN THE SQUARE        |
 !                               SECTION OF THE ARROW MATRIX             |
 !       ROWMAT(IARROW,NOD) -       HOLDS THE ELEMENTS OF THE ROWS       |
 !                               OF THE ARROW MATRIX                     |
 !       COLMAT(NOD,IARROW) -       HOLDS THE ELEMENTS OF THE COLUMNS    |
 !                               OF THE ARROW MATRIX                     |
 !       AMAT(NOD, 3*IBAND-1 ) -    BANDED MATRIX A.  NOTE THAT THE      |
 !              ---------        DIMENSIONS ARE 50% GREATER THAN         |
 !                  2            THE BANDWIDTH.                          |
 !-----------------------------------------------------------------------|

 MID   = (IBAND+1)/2
 ISIGN = 1
 MEX   = 0
 NEXX  = 0
 NPIV  = 0
 DETLOG = 0.0D0

 ! NCODing !
 ! ------- !

!  $$$      Do N = 1, NOD
!  $$$
!  $$$         If (NCOD(N) .EQ. 1) Then
!  $$$
!  $$$            AMAT(N, :) = 0.0D0
!  $$$            COLMAT(N, :) = 0.0D0
!  $$$
!  $$$            AMAT(N, MID) = 1.0D0
!  $$$            RHS(N) = 0.0D0
!  $$$
!  $$$         End If
!  $$$
!  $$$      End Do
!  $$$
!  $$$      Do N = 1, IARROW
!  $$$
!  $$$         If (NCOD(NOD + N) .EQ. 1) Then
!  $$$
!  $$$            ROWMAT(N, :) = 0.0D0
!  $$$            HEAD(N, :) = 0.0D0
!  $$$
!  $$$            HEAD(N, N) = 1.0D0
!  $$$            RHS(NOD + N) = 0.0D0
!  $$$
!  $$$         End If
!  $$$
!  $$$      End Do

 !-----------------------------------------------------------------------|
 !       AMAT(I,MID)        CONTAINS THE DIAGONAL ELEMENTS OF MATRIX A   |
 !       MEX             KEEPS TRACK OF THE MAXIMUM ROW EXCHANGE         |
 !                       IN THE MATRIX A                                 |
 !       NEXX            KEEPS TRACK OF THE LAST NON-ZERO ELEMENT        |
 !                       IN THE ROW IN MATRIX A                          |
 !                                                                       |
 !       CLEAN OUT MATRIX A OUTSIDE THE BAND WIDTH                       |
 !-----------------------------------------------------------------------|

 DO IS = 1, NOD
    DO JS = IBAND + 1, IBAND + MID - 1
       AMAT(IS, JS) = 0.0D0
    END DO
 END DO

 !-----------------------------------------------------------------------|
 !       GAUSSIAN ELIMINATION OF MATRIX A.  CHECK THE COLUMN             |
 !       FOR THE MAXIMUM ELEMENT.  IF THE DIAGONAL IS LESS THAN A        |
 !       FRACTION, TAU, OF THE MAXIMUM ELEMENT THEN SWITCH ROWS.         |
 !       IF TAU = 0.0 THEN NO PIVOTING WILL TAKE PLACE AND THE  L U      |
 !       DECOMPOSITION WILL BE RETURNED.  IF TAU=1.0 THEN FULL ROW       |
 !       PIVOTING WILL BE USED.  IF TAU = -1.0 THEN THE  L U  MATRIX     |
 !       WILL BE FRONT-SUBSTITUTED, AND THEN BACK-SUBSTITUTED.           |
 !                                                                       |
 !-----------------------------------------------------------------------|

 IF (TAU.LT.0.0) THEN

    !-----------------------------------------------------------------------|
    !       IF TAU < 0 THE MATRIX IS NOT DECOMPOSED, AND THE PREVIOUS       |
    !       L U  MATRIX IS USED.  THE LOWER MATRIX  L  IS FORWARD           |
    !       SUBSTITUTED.                                                    |
    !-----------------------------------------------------------------------|
    DO IS = 1, NOD

       IS1 = MIN(IS + MID - 1, NOD)

       DO N = 1, NRHS

          PIV = RHS(IS)

          DO JS = IS + 1, IS1
             RHS(JS) = RHS(JS) - PIV*AMAT(JS, IS + MID - JS)
          END DO

          DO JS = 1, IARROW
             RHS(NOD + JS) = RHS(NOD + JS) - PIV*ROWMAT(JS, IS)
          END DO

       END DO

    END DO

    DO IS = 1, IARROW
       DO N = 1, NRHS

          PIV = RHS(NOD + IS)

          DO JS = IS + 1, IARROW
             RHS(NOD + JS) = RHS(NOD + JS) - PIV*HEAD(JS, IS)
          END DO

       END DO
    END DO

 ELSE

    !-----------------------------------------------------------------------|
    !       IF TAU > 0

    DO IS = 1, NOD - 1

       LLS  = IS
       BIG = DABS(AMAT(IS, MID))

       DO JS = IS + 1, MIN(LLS + MID - 1, NOD)

          IF(DABS(AMAT(JS, MID + IS - JS)).LT.BIG) CYCLE

          BIG = DABS(AMAT(JS, MID + IS - JS))
          LLS = JS

       END DO

       IF (DABS(AMAT(IS, MID)).LT.TAU*BIG) THEN
          !-----------------------------------------------------------------------|
          !       EXCHANGE ROW I WITH ROW LL.  THE VECTOR B AND ARRAY COL ARE     |
          !       EXCHANGED.  THE DETERMINANT OF THE MATRIX CHANGES SIGN.         |
          !-----------------------------------------------------------------------|
          MEX  = MAX(MEX,  LLS - IS)
          NEXX = MAX(NEXX, LLS - IS)
          KS1   = MIN(IBAND + NEXX, NOD - IS + MID)

          DO JS = MID, KS1
             TRANS   = AMAT(IS, JS)
             AMAT(IS, JS) = AMAT(LLS, JS + IS - LLS)
             AMAT(LLS, JS + IS - LLS) = TRANS
          END DO

          DO JS = 1, IARROW
             TRANS = COLMAT(IS, JS)
             COLMAT(IS, JS)  = COLMAT(LLS, JS)
             COLMAT(LLS, JS) = TRANS
          END DO

          DO N = 1, NRHS
             TRANS = RHS(IS)
             RHS(IS)  = RHS(LLS)
             RHS(LLS) = TRANS
          END DO

          ISIGN = -ISIGN
          NPIV  = NPIV + 1

       END IF
       !-----------------------------------------------------------------------|
       !       THE LOWER BAND MATRIX IS ELIMINATED AND REPLACED WITH THE       |
       !       PIVOTS FOR THE L U DECOMPOSITION.  THE B VECTOR AND THE         |
       !       ARRAY COL IS UPDATED.                                           |
       !-----------------------------------------------------------------------|

       DIAG = AMAT(IS, MID)
       NEXX = MAX(LLS - IS, NEXX - 1)
       JS1 = MAX(1, IS + MID - NOD)
       JS2 = MAX(1 - NEXX, IS + MID - NOD)

       DO JS = MID - 1, JS1, -1

          IS1  = MID + IS - JS
          PIV = AMAT(IS1, JS)/DIAG
          AMAT(IS1, JS) = PIV

          DO N = 1, NRHS
             RHS(IS1) = RHS(IS1) - PIV*RHS(IS)
          END DO

          DO MS = 1, IARROW
             COLMAT(IS1, MS) = COLMAT(IS1,MS) - PIV*COLMAT(IS, MS)
          END DO

          DO KKS = JS + 1, JS + MID - JS2
             AMAT(IS1, KKS) = AMAT(IS1, KKS) - PIV*AMAT(IS, MID + KKS - &
                  &           JS)
          END DO

       END DO

    END DO

    !-----------------------------------------------------------------------|
    !       THE ARRAY ROW IS ELIMINATED EXCEPT FOR THE LAST COLUMN AND      |
    !       REPLACED WITH THE PIVOT FOR THE L U DECOMPOSITION.  THE B       |
    !       VECTOR AND THE ARRAY HEAD ARE UPDATED.                          |
    !-----------------------------------------------------------------------|
    DO IS = 1, NOD - 1

       DIAG = AMAT(IS, MID)
       JS1 = MIN(IS + MID - 1 + MEX, NOD)

       DO JS = 1, IARROW

          PIV = ROWMAT(JS, IS)/DIAG
          ROWMAT(JS, IS)=PIV

          DO N = 1, NRHS
             RHS(NOD + JS) = RHS(NOD + JS) - PIV*RHS(IS)
          END DO

          DO KKS = IS + 1, JS1
             ROWMAT(JS, KKS) = ROWMAT(JS, KKS) - PIV*AMAT(IS, MID + KKS &
                  &           - IS)
          END DO

          DO LLS =1, IARROW
             HEAD(JS, LLS) = HEAD(JS, LLS) - PIV*COLMAT(IS, LLS)
          END DO

       END DO

    END DO

    !-----------------------------------------------------------------------|
    !       THE FINAL ELEMENT OF THE BANDED MATRIX A IS CHECKED TO SEE      |
    !       IF IT IS LARGE ENOUGH.  IF IT IS TOO SMALL, IT IS PIVOTED       |
    !       WITH AN ELEMENT FROM THE LAST COLUMN OF THE ARRAY ROW.          |
    !-----------------------------------------------------------------------|

    BIG = DABS(AMAT(NOD, MID))
    !        write(*,*) LLS
    DO JS = 1, IARROW

       IF (DABS(ROWMAT(JS, NOD)).LT.BIG) CYCLE

       BIG = DABS(ROWMAT(JS, NOD))
       LLS  = JS
       !          write(*,*) lls
    END DO

    IF (DABS(AMAT(NOD,MID)).LT.TAU*BIG.AND.LLS.LE.IARROW) THEN

       DO N = 1, NRHS
          TRANS = RHS(NOD)
          !            write(*,*) lls
          RHS(NOD) = RHS(NOD + LLS)
          RHS(NOD+LLS) = TRANS
       END DO

       DO JS = 1,IARROW
          TRANS = COLMAT(NOD, JS)
          COLMAT(NOD, JS) = HEAD(LLS, JS)
          HEAD(LLS, JS)    = TRANS
       END DO

       TRANS = AMAT(NOD, MID)
       AMAT(NOD, MID) = ROWMAT(LLS, NOD)
       ROWMAT(LLS, NOD) = TRANS
       ISIGN = -ISIGN
       NPIV  = NPIV + 1

    END IF

    !-----------------------------------------------------------------------|
    !       THE LAST COLUMN OF THE ARRAY ROW IS ELIMINATED.  THE B          |
    !       VECTOR AND THE ARRAY HEAD ARE UPDATED.                          |
    !-----------------------------------------------------------------------|

    DIAG = AMAT(NOD, MID)

    DO IS = 1, IARROW

       PIV = ROWMAT(IS, NOD)/DIAG
       ROWMAT(IS, NOD) = PIV

       DO N = 1, NRHS
          RHS(NOD + IS) = RHS(NOD + IS) - PIV*RHS(NOD)
       END DO

       DO JS = 1, IARROW
          HEAD(IS, JS) = HEAD(IS, JS) - PIV*COLMAT(NOD, JS)
       END DO

    END DO

    !-----------------------------------------------------------------------|
    !       THE ARRAY HEAD IS ELIMINATED USING THRESHOLD ROW PIVOTING.      |
    !       IF THE DIAGONAL ELEMENT IS LESS THAN A FRACTION, TAU, OF THE    |
    !       MAXIMUM COLUMN ELEMENT, THEN THE ROWS ARE EXCHANGED.  AS THE    |
    !       COLUMN IS ELIMINATED, THE B VECTOR AND THE ARRAY HEAD ARE       |
    !       UPDATED.                                                        |
    !-----------------------------------------------------------------------|

    DO IS = 1,IARROW - 1

       BIG = DABS(HEAD(IS, IS))

       DO JS = IS + 1, IARROW

          IF (DABS(HEAD(JS, IS)).LT.BIG) CYCLE

          BIG = DABS(HEAD(JS, IS))
          LLS  = JS

       END DO

       IF (DABS(HEAD(IS, IS)).LT.TAU*BIG) THEN

          DO N = 1, NRHS
             TRANS = RHS(NOD + IS)
             RHS(NOD + IS) = RHS(NOD + LLS)
             RHS(NOD + LLS) = TRANS
          END DO

          DO JS = 1, IARROW
             TRANS = HEAD(IS, JS)
             HEAD(IS,  JS) = HEAD(LLS, JS)
             HEAD(LLS, JS) = TRANS
          END DO

          ISIGN = -ISIGN
          NPIV  = NPIV + 1

       END IF

       DIAG = HEAD(IS, IS)

       DO JS = IS + 1, IARROW

          PIV = HEAD(JS, IS)/DIAG
          HEAD(JS, IS) = PIV

          DO N = 1, NRHS
             RHS(NOD + JS) = RHS(NOD + JS) - PIV*RHS(NOD + IS)
          END DO

          DO KKS = IS + 1,IARROW
             HEAD(JS, KKS) = HEAD(JS, KKS) - PIV*HEAD(IS, KKS)
          END DO

       END DO

    END DO

 END IF

 !-----------------------------------------------------------------------|
 !       THE UPPER MATRIX U CONTAINING THE ARROW IS BACK-SUBSTITUTED     |
 !       FOR THE SOLUTION VECTORS X.  THE MAGNITUDE AND SIGN OF THE      |
 !       MATRIX ARE CALCULATED.                                          |
 !-----------------------------------------------------------------------|
 DO IS = IARROW, 1, -1

    DIAG   = HEAD(IS, IS)
    DETLOG = DETLOG + DLOG10(DABS(DIAG))

    IF (DIAG.LT.0.0) ISIGN = -ISIGN

    DO N = 1, NRHS
       PIV = RHS(NOD + IS)/DIAG
       RHS(NOD + IS) = PIV

       DO JS = IS - 1, 1, -1
          RHS(NOD + JS) = RHS(NOD + JS) - PIV*HEAD(JS, IS)
       END DO

       DO KKS = NOD, 1, -1
          RHS(KKS) = RHS(KKS) - PIV*COLMAT(KKS, IS)
       END DO

    END DO

 END DO

 !-----------------------------------------------------------------------|
 !       THE UPPER MATRIX U CONTAINING THE BANDED MATRIX A IS BACK-      |
 !       SUBSTITUTED FOR THE SOLUTION VECTORS X.  THE MAGNITUDE AND      |
 !       SIGN OF THE DETERMINANT ARE CALCULATED.                         |
 !-----------------------------------------------------------------------|
 DO IS = NOD, 1, -1

    DIAG   = AMAT(IS, MID)
    DETLOG = DETLOG + DLOG10(DABS(DIAG))

    IF(DIAG.LT.0.0) ISIGN = -ISIGN

    JS1 = MAX(1, IS - MID - 1 - MEX)

    DO N = 1, NRHS

       PIV = RHS(IS)/DIAG
       RHS(IS) = PIV

       DO JS = IS - 1, JS1, -1
          RHS(JS) = RHS(JS) - AMAT(JS, MID + IS - JS)*PIV
       END DO

    END DO

 END DO


 RETURN

END SUBROUTINE ArrowSolver
!***********************************************************************************************************************
! The Arrow Solver END
!***********************************************************************************************************************
