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
*
      DIMENSION DTQUEN ( MXTRCK, MAXQMG )
*
      CHARACTER*20 FILNAM
      LOGICAL NEUTRONS
      INTEGER ind
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /
*
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
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
      WRITE (IODRAW) NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
     &               SNGL (WTRACK)
      WRITE (IODRAW) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
     &                 SNGL (ZTRACK (I)), I = 0, NTRACK ),
     &               ( SNGL (DTRACK (I)), I = 1, MTRACK ),
     &                 SNGL (CTRACK)
*  +-------------------------------------------------------------------*
*  |  Quenching is activated
      IF ( LQEMGD ) THEN
         IF ( MTRACK .GT. 0 ) THEN
            RULLL  = ZERZER
            CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
            WRITE (IODRAW) ( ( SNGL (DTQUEN (I,JBK)), I = 1, MTRACK ),
     &                         JBK = 1, NQEMGD )
         END IF
      END IF
*  |  End of quenching
*  +-------------------------------------------------------------------*
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
!       IF ( .NOT. LFCOPE ) THEN
!             LFCOPE = .TRUE.
!             IF ( KOMPUT .EQ. 2 ) THEN
!                FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
!             ELSE
!                FILNAM = CFDRAW
!             END IF
!             OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
!      &          'FORMATTED' )
!          END IF

      ! IF (MREG .EQ.8 .AND. NEWREG.EQ.9 .AND. JTRACK .EQ. 8 ) THEN ! Neutron crossed into the TPC
      !       WRITE(99,*) 'Neutron Entered TPC!'
      !       WRITE(99,*) 'NCASE:', NCASE
      !       WRITE(99,*) 'ETRACK', ETRACK
      !       WRITE(99,*) 'LTRACK', LTRACK
      !       ! Tagging the particle for killing by setting weight = 0 in usrmed.f routine
      !       LLOUSE = 0
      ! END IF
      

      RETURN
*
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
*
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
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
      WRITE (IODRAW)  SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
*  +-------------------------------------------------------------------*
*  |  Quenching is activated : calculate quenching factor
*  |  and store quenched energy in DTQUEN(1, jbk)
      IF ( LQEMGD ) THEN
         RULLL = RULL
         CALL QUENMG ( ICODE, MREG, RULLL, DTQUEN )
         WRITE (IODRAW) ( SNGL (DTQUEN(1, JBK)), JBK = 1, NQEMGD )
      END IF
*  |  end quenching

      ! ! This is a neutron directly produced by a muon somewhere within the OD.
      IF (MREG .LE. 8 .AND. JTRACK .EQ. 8 .AND. LTRACK .EQ. 2 .AND. LLOUSE .EQ. 0) THEN
            WRITE(98, *) ICODE, JTRACK, MREG, LTRACK, ETRACK, 
     &       XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &      (ISPUSR(I),I=1,4), (SPAUSR(I),I=1,7) ! Parent data

            LLOUSE = 2
      END IF

      ! This is a neutron detected within the TPC Xenon
      IF (MREG .EQ. 1 .AND. JTRACK .EQ. 8 .AND. LLOUSE .NE. 1) THEN 
            WRITE(98, *) ICODE, JTRACK, MREG, LTRACK, ETRACK, 
     &       XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &      (ISPUSR(I),I=1,4), (SPAUSR(I),I=1,7) ! Parent data

            ! Tagging the particle for ignoring by setting LLOUSE = 1
            LLOUSE = 1
      END IF


      ! ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
      ! Setting parent values
      ! Integer values
      ISPUSR(1) = ICODE
      ISPUSR(2) = JTRACK
      ISPUSR(3) = MREG
      ISPUSR(4) = LTRACK
      ! Float values
      SPAUSR(1) = ETRACK
      SPAUSR(2) = XSCO
      SPAUSR(3) = YSCO
      SPAUSR(4) = ZSCO
      SPAUSR(5) = CXTRCK
      SPAUSR(6) = CYTRCK
      SPAUSR(7) = CZTRCK
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF
      WRITE (IODRAW) -NCASE, NPFLKA, NSTMAX, SNGL (TKESUM),
     &                SNGL (WEIPRI)
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope: it works only for 1 source particle on
*  |  the stack for the time being
      IF ( ILOFLK (NPFLKA) .GE. 100000 .AND. LRADDC (NPFLKA) ) THEN
         IARES  = MOD ( ILOFLK (NPFLKA), 100000  )  / 100
         IZRES  = MOD ( ILOFLK (NPFLKA), 10000000 ) / 100000
         IISRES = ILOFLK (NPFLKA) / 10000000
         IONID  = ILOFLK (NPFLKA)
         WRITE (IODRAW) ( IONID,SNGL(-TKEFLK(I)),
     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
      ELSE IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
         IONID = ILOFLK (NPFLKA)
         CALL DCDION ( IONID )
         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-IONID)),
     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: ???
      ELSE IF ( ILOFLK (NPFLKA) .LT. -6 ) THEN
         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-ILOFLK(NPFLKA))),
     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
         WRITE (IODRAW) ( ILOFLK(I), SNGL (TKEFLK(I)+AM(ILOFLK(I))),
     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
      END IF
*  |
*  +-------------------------------------------------------------------*
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
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
     &          'UNFORMATTED' )
      END IF

      ! ! This is a neutron produced by a muon somewhere within the OD.
      IF (MREG .LE. 8 .AND. JTRACK .EQ. 8 .AND. LLOUSE .EQ. 0) THEN
            WRITE(70, *) ICODE, NCASE, JTRACK, MREG, LTRACK, ETRACK, 
     &       XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &      (ISPUSR(I),I=1,4), (SPAUSR(I),I=1,7) ! Parent data

            LLOUSE = 2
      END IF

      ! This is a neutron detected within the TPC Xenon
      IF (MREG .EQ. 1 .AND. JTRACK .EQ. 8 .AND. LLOUSE .NE. 1) THEN 
            WRITE(72, *) ICODE, NCASE, JTRACK, MREG, LTRACK, ETRACK, 
     &       XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK,
     &      (ISPUSR(I),I=1,4), (SPAUSR(I),I=1,7) ! Parent data

            ! Tagging the particle for ignoring by setting LLOUSE = 1
            LLOUSE = 1
      END IF
      
      
      ! ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK
      ! Setting parent values
      ! Integer values
      ISPUSR(1) = ICODE
      ISPUSR(2) = JTRACK
      ISPUSR(3) = MREG
      ISPUSR(4) = LTRACK
      ! Float values
      SPAUSR(1) = ETRACK
      SPAUSR(2) = XSCO
      SPAUSR(3) = YSCO
      SPAUSR(4) = ZSCO
      SPAUSR(5) = CXTRCK
      SPAUSR(6) = CYTRCK
      SPAUSR(7) = CZTRCK
      

* No output by default:
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

