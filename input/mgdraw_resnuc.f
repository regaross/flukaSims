*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2003-2019:  CERN & INFN                            *
*     All Rights Reserved.                                             *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created on   01 March 1990   by        Alfredo Ferrari           *
*                                              INFN - Milan            *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE 'caslim.inc'
      INCLUDE 'comput.inc'
      INCLUDE 'sourcm.inc'
      INCLUDE 'fheavy.inc'
      INCLUDE 'flkstk.inc'
      INCLUDE 'genstk.inc'
      INCLUDE 'mgddcm.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'quemgd.inc'
      INCLUDE 'sumcou.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'resnuc.inc'
      INCLUDE 'rdpstk.inc'
      INCLUDE 'rdcycm.inc'
*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = 96, STATUS = 'UNKNOWN', FORM =
     &          'FORMATTED')
      END IF
      RETURN
*
*======================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
      RETURN
*
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*
      ENTRY EEDRAW ( ICODE )
      RETURN
*                                                                      *
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*             16: recoil from (heavy) bremsstrahlung                   *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )
      RETURN
*                                                                      *
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY SODRAW
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*             110: radioactive decay products                          *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*             237: mu pair production     secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*     Icode = 50x: call from synstp                                    *
*             500: synchrotron radiation photons from e-/e+            *
*             501: synchrotron radiation photons from other charged    *
*                  particles                                           *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=npflka          *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = 96, STATUS = 'UNKNOWN', FORM =
     &          'FORMATTED')
      END IF

      ! REDUCED OUTPUT
      IF ((ICODE .EQ. 100 .OR. ICODE .EQ. 101 .OR. ICODE .EQ. 102 .OR.
     &     ICODE .EQ. 300 ) .AND. 
     &    (MREG .EQ. 1 .OR. MREG .EQ. 2)) THEN
         IF ( NPHEAV .GT. 0 .OR. IBRES .GT. 1) THEN
             WRITE (96, *)
             WRITE (96, '(A)') 'MGDRAW'
             WRITE (96, '(A, I3, A, I3)') 'Reaction with ICODE ', ICODE, 
     &       ' in region ', MREG
             WRITE(96, '(A, I15)') 'Primary history: ', NCASE
             WRITE(96, '(A, 2I15)') 'JTRACK, LTRACK: ', JTRACK, LTRACK
             WRITE(96, '(A, I5, A, 100(I3))') 'Generated ', NP-NP0, 
     &      ' particles of ID: ', (KPART(I),I=NP0+1, NP)
             WRITE(96, '(A, I5, A)') 'Generated ', NPHEAV, 
     &       ' fragments with ' 
             WRITE(96, '(A, 10I5)') '     A = ', 
     &       (IBHEAV(KHEAVY(I)),I=1,NPHEAV)
             WRITE(96, '(A, 10I5)') '     Z = ', 
     &       (ICHEAV(KHEAVY(I)),I=1,NPHEAV)
             WRITE(96, '(A, 2(a, I3))') 'Generated a residual ',
     &       'with A = ', IBRES, ' , Z = ', ICRES
         END IF
      END IF

      IF (ICODE .EQ. 110 .AND. (MREG .EQ. 1 .OR. MREG .EQ. 2)) THEN
         WRITE (96, *)
         WRITE (96, '(A)') 'MGDRAW' 
         WRITE (96, '(A, I3, A, I3)') 'Reaction with ICODE ', ICODE, 
     &   ' (radioactive decay products) in region ', MREG
         WRITE(96, '(A, I15)') 'Primary history: ', NCASE
         WRITE(96, '(A, 2(a, I3))') 'Decay of nucleus ',
     &   'with A = ', IARDCY, ' , Z = ', IZRDCY

         WRITE(96, '(A)') 'Decay products: '
         DO I = 1, NP
            WRITE(96, '(I3)') KPART(I)
         END DO

      END IF
      RETURN
      END

