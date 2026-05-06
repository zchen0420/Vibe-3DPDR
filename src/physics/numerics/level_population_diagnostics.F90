  SUBROUTINE GAUSS_JORDAN_writes(A,N,NP,B,ill)

    use definitions
    use healpix_types
    use maincode_module, only : thermal

    IMPLICIT NONE
    integer(kind=i4b), intent(in) :: ill!,coolant
    INTEGER(kind=i4b):: I,J,K,L,LL,IROW,ICOL
    INTEGER(kind=i4b), intent(in):: N,NP!,M,MPP
    integer(kind=i4b), PARAMETER :: NMAX=100
    INTEGER(kind=i4b):: IPIV(1:NMAX),INDXR(1:NMAX),INDXC(1:NMAX)
    real(kind=dp), intent(inout) :: A(1:NP,1:NP)
    real(kind=dp), intent(inout) :: B(1:NP)!,1:MPP)
    real(kind=dp) :: BIG,DUM,PIVINV
    integer(kind=i4b) :: diagnostic_point

    diagnostic_point = 1

    write(6,*) 'b'
    do i=1,np
      write(6,*) b(i),i
    enddo

    write(6,*) 'a'
    do i=1,np
      do j=1,np
        write(6,*) a(i,j)
      enddo
    enddo


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
              ENDIF !ABS(A
            ELSE IF(IPIV(K).GT.1) THEN
              PRINT *,'ERROR! Singular matrix in GAUSS_JORDAN'
              write(6,*) 'Crashed in first loop'
              !                     write(6,*) 'grid point = ',p, ' coolant = ',coolant
              write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(diagnostic_point)
              STOP
            ENDIF !IPIV(K).EQ.0
          ENDDO !K=1,N
        ENDIF !IPIV(J).NE.1
      ENDDO !J=1,N
      IPIV(ICOL)=IPIV(ICOL)+1
      IF(IROW.NE.ICOL) THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        ENDDO !L=1,N
        !            DO L=1,M
        !               DUM=B(IROW,L)
        !               B(IROW,L)=B(ICOL,L)
        !               B(ICOL,L)=DUM
        !            ENDDO
        !================================
        DUM=B(IROW)
        if (i.eq.ill) write(6,*) 'dum=',dum,'A'
        !write(6,*) 'DUM=',DUM
        B(IROW)=B(ICOL)
        if (i.eq.ill) write(6,*) 'b(',irow,')=',b(irow),'B'
        !write(6,*) 'irow=',irow
        !write(6,*) 'B(irow)=',b(irow)
        B(ICOL)=DUM
        if (i.eq.ill) write(6,*) 'b(',icol,')=',b(icol),'C'
        !write(6,*) 'icol=',icol
        !write(6,*) 'b(icol)=',b(icol)
        !================================
      ENDIF !IROW.NE.ICOL
      INDXR(I)=IROW
      INDXC(I)=ICOL
      IF(A(ICOL,ICOL).EQ.0.0D0) THEN
        PRINT *,'ERROR! Singular matrix found by GAUSS_JORDAN'
        write(6,*) 'Crashed in second loop'
        !            write(6,*) 'grid point = ',p, ' coolant = ',coolant
        write(6,*) 'thermal%gas_temperature = ',thermal%gas_temperature(diagnostic_point)
        STOP
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
      if (i.eq.ill) write(6,*) b(icol),pivinv,'D'
      B(ICOL)=B(ICOL)*PIVINV
      if (i.eq.ill) write(6,*) b(icol),'E'
      !write(6,*) 'pivinv=',pivinv
      !write(6,*) 'b(icol)=',b(icol)
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
          if (i.eq.ill) then
            write(6,*) 'll=',ll,'F'
            write(6,*) 'b(ll)=',b(ll),'G'
            write(6,*) 'b(icol)=',b(icol),'H'
            write(6,*) 'dum=',dum,'I'
          endif
          B(LL)=B(LL)-B(ICOL)*DUM
          if (i.eq.ill) write(6,*) 'b(ll) after=',b(ll),'J'
          !=============================================
        ENDIF !LL.NE.ICOL
      ENDDO !LL=1,N
    ENDDO ! I=1,N
    DO L=N,1,-1
      IF(INDXR(L).NE.INDXC(L)) THEN
        DO K=1,N
          DUM=A(K,INDXR(L))
          A(K,INDXR(L))=A(K,INDXC(L))
          A(K,INDXC(L))=DUM
        ENDDO !K=1,N
      ENDIF !INDXR(L).NE.INDXC(L)
    ENDDO !L=N,1,-1
    do i=1,n
      write(6,*) b(i),i
    enddo
    RETURN
  END subroutine
