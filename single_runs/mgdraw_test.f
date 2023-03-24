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
*     Muon collider interaction region user defined scoring.           *
*     By D.Calzolari                                                   *
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
      INCLUDE 'flkmat.inc'
*
*
      PARAMETER (IPLANE1 = 81)
      LOGICAL LPRINT, LFCOPE, LPREVIOUS
      CHARACTER*8 NRGNAM, MRGNAM
      SAVE LFCOPE
      SAVE ZPREVIOUS
      DATA LFCOPE / .FALSE. /
      DATA LPREVIOUS / .TRUE. /
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

      
      ! Initialization
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         OPEN (UNIT = IPLANE1, FILE = "Part_list.dat")
         CALL WRITE_PREAMBLE(IPLANE1)
      END IF
*
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
      CALL GEOR2N ( MREG, MRGNAM, IERR1 ) 
      CALL GEOR2N ( NEWREG, NRGNAM, IERR1 )
      LPRINT = .FALSE.
      IOUTPUT = 0
      ! Hardcoded conditions for each planes
      ! Plane 1: particles leaving the slab and going in vacuum
      IDX_MATERIAL_NEW = MEDFLK ( NEWREG, 1 )
      IDX_MATERIAL_OLD = MEDFLK ( MREG, 1 )
      IF ( IDX_MATERIAL_NEW .EQ. 2 .AND. IDX_MATERIAL_OLD .NE. 2 ) THEN
         LPRINT = .TRUE.
         IOUTPUT = IPLANE1
      ENDIF
      ! Score the particle if needed
      IF ( LPRINT ) THEN
         WRITE( IOUTPUT,'(1(I4), 7(1PE16.8), 2(I16), 7(1PE16.8))' ) 
     &               JTRACK, ETRACK, 
     &               XTRACK(NTRACK), YTRACK(NTRACK), ZTRACK(NTRACK),
     &               CXTRCK, CYTRCK, CZTRCK,
     &  ISPUSR(1), ISPUSR(2), SPAUSR(1), SPAUSR(2), SPAUSR(3), SPAUSR(4),
     &  SPAUSR(5), SPAUSR(6), SPAUSR(7)
      ENDIF
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
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
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
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=npflka          *
*                                                                      *
*======================================================================*
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      ! Write out the event attributes, type, parent, position, region
      WRITE(72,*) ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO,
     & CXTRCK, CYTRCK, CZTRCK 
      ! Write out how many product particles were produced and their respective codes
      WRITE(72,*) NP, (KPART(I),I=1,NP)
      ! Write out fragments produce (pieces of nuclei) A followed by Z
      WRITE(72,*) NPHEAV, (IBHEAV(KHEAVY(I)), ICHEAV(KHEAVY(I)), I=1,NPHEAV)
      ! Write out the Residual A and Z (if there is one)
      WRITE(72,*) ICRES, IBRES
      ! Write out the parent attributes
      WRITE(72,*) (ISPUSR(I),I=1,4), (SPAUSR(I),I=1,7)


      ! ICODE, JTRACK, MREG, LTRACK, ETRACK, XSCO, YSCO, ZSCO, CXTRCK, CYTRCK, CZTRCK

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

      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END

      SUBROUTINE WRITE_PREAMBLE ( IUNIT )
         INTEGER IUNIT
         WRITE( IUNIT, '(1(A4), 17(A16))' ) ! Preamble
     &               "#ID ", "KIN",
     &               "X", "Y", "Z", "CX", "CY", "CZ",
     & "id_mo", "id_pr", "X_mo", "Y_mo", "Z_mo", "CX_mo", "CY_mo", 
     & "CZ_mo", "Kin_mo"
      END SUBROUTINE

