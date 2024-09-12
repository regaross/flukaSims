*
*=== Usrrnc ===========================================================*
*
      SUBROUTINE USRRNC ( IZ, IA, IS, X, Y, Z, MREG, WEE, ICALL )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'resnuc.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*     USeR Residual NuClei:                                            *
*                                                                      *
*     Created on   06 April 2005   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*----------------------------------------------------------------------*
*
      LOGICAL LFIRST
      DATA LFIRST / .TRUE. /
      SAVE LFIRST

      IF ( MREG .EQ. 1 .OR. MREG .EQ. 2 ) THEN
        IF ( IA .GE. 5 .AND. IZ .GE. 2 ) THEN
          WRITE (96 ,'(A)') 'USRRNC'
          WRITE(96, '(4(a, I3))') 'Scored residual in region ', MREG, 
     & ' with A = ', IA, ' Z = ', IZ, ' IS = ', IS
        END IF
      END IF
      RETURN
*=== End of subroutine Usrrnc =========================================*
      END
