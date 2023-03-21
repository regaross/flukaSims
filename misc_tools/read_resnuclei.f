      PROGRAM READRN
      CHARACTER*125 LINE, FILINP, FILOUT
      PARAMETER (MAXZ = 72, MINNMZ = -4, MAXNMZ = 135, K = -5)
      DIMENSION RESULT(MAXZ, MINNMZ-K:MAXNMZ-K)

      WRITE(*,*) "Filename?"
      READ(*,'(A)') FILINP
      OPEN(UNIT=1, FILE=FILINP, STATUS='OLD')
      LQ = INDEX(FILINP,' ') - 1
      FILOUT = FILINP(1:LQ)//'.rn'
      OPEN(UNIT=2, FILE=FILOUT, STATUS='UNKNOWN')

      DO 1 I = 1, 14
         READ(1,'(A)') LINE    !  skip header lines
  1   CONTINUE

      READ(1,100,END=4) RESULT
  4   CONTINUE

      WRITE(2,'(A)')   '   Z   A    Residual nuclei'
      WRITE(2,'(A,/)') '         per cm**3 per primary'
      DO 2 I = 1, MAXZ
         DO 3 J = MINNMZ-K, MAXNMZ-K
            IF(RESULT(I,J) .GT. 0.D0)
     &       WRITE(2,'(2I4,1P, G15.6)') I, J+K+2*I, RESULT(I,J)
  3      CONTINUE
  2   CONTINUE
 100  FORMAT(1(5X,1P,10(1X,E11.4)))
      END
