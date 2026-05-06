  SUBROUTINE GAUSS_JORDAN(A,N,NP,B,call_writes)

    use definitions
    use healpix_types
    IMPLICIT NONE
    !      integer(kind=i4b), intent(in) :: coolant
    INTEGER(kind=i4b):: I,J,K,L,LL,IROW,ICOL
    INTEGER(kind=i4b), intent(in):: N,NP!,M,MPP
    integer(kind=i4b), PARAMETER :: NMAX=100
    INTEGER(kind=i4b):: IPIV(1:NMAX),INDXR(1:NMAX),INDXC(1:NMAX)
    real(kind=dp), intent(inout) :: A(1:NP,1:NP)
    real(kind=dp), intent(inout) :: B(1:NP)!,1:MPP)
    real(kind=dp) :: BIG,DUM,PIVINV
    logical,intent(out)::call_writes

    ICOL=0
    IROW=0
    IPIV=0
    DO I=1,N
      BIG=0.0D0
      DO J=1,N
        IF(IPIV(J).NE.1) THEN
          DO K=1,N
            IF(IPIV(K).EQ.0) THEN
              IF(ABS(A(J,K)).GE.BIG) THEN
                BIG=ABS(A(J,K))
                IROW=J
                ICOL=K
              ENDIF
            ELSE IF(IPIV(K).GT.1) THEN
              PRINT *,'ERROR! Singular matrix in GAUSS_JORDAN'
              call_writes=.true.
              return
              !                     write(6,*) 'Crashed in first loop'
              !                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
              !                     write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(p)
              !                     STOP
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IPIV(ICOL)=IPIV(ICOL)+1
      IF(IROW.NE.ICOL) THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        ENDDO
        !            DO L=1,M
        !               DUM=B(IROW,L)
        !               B(IROW,L)=B(ICOL,L)
        !               B(ICOL,L)=DUM
        !            ENDDO
        !================================
        DUM=B(IROW)
        B(IROW)=B(ICOL)
        B(ICOL)=DUM
        !================================
      ENDIF
      INDXR(I)=IROW
      INDXC(I)=ICOL
      IF(A(ICOL,ICOL).EQ.0.0D0) THEN
        PRINT *,'ERROR! Singular matrix found by GAUSS_JORDAN'
        call_writes=.true.
        return
        !            write(6,*) 'Crashed in second loop'
        !            write(6,*) 'grid point = ',p, ' coolant = ',coolant
        !            write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(p)
        !            STOP
      ENDIF
      PIVINV=1.0D0/A(ICOL,ICOL)
      A(ICOL,ICOL)=1.0D0
      DO L=1,N
        A(ICOL,L)=A(ICOL,L)*PIVINV
      ENDDO
      !         DO L=1,M
      !            B(ICOL,L)=B(ICOL,L)*PIVINV
      !         ENDDO
      !=======================================
      B(ICOL)=B(ICOL)*PIVINV
      !=======================================
      DO LL=1,N
        IF(LL.NE.ICOL) THEN
          DUM=A(LL,ICOL)
          A(LL,ICOL)=0.0D0
          DO L=1,N
            A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
          ENDDO
          !               DO L=1,M
          !                  B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
          !               ENDDO
          !=============================================
          B(LL)=B(LL)-B(ICOL)*DUM
          !=============================================
        ENDIF
      ENDDO
    ENDDO
    DO L=N,1,-1
      IF(INDXR(L).NE.INDXC(L)) THEN
        DO K=1,N
          DUM=A(K,INDXR(L))
          A(K,INDXR(L))=A(K,INDXC(L))
          A(K,INDXC(L))=DUM
        ENDDO
      ENDIF
    ENDDO
    RETURN
  END subroutine
