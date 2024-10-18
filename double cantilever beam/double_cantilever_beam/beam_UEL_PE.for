C ======================================================================
C User Subroutine UEL for Abaqus: two elements for the split scheme
C  operator damage gradient and displacement problem:
C ----- Two-dimensional solid element ------
C U1: C2DPH3 phase-field damage gradient 3-node linear
C U2: C2DPH4 phase-field damage gradient 4-node bilinear
C ----- Three-dimensional solid elements ----- 
C U3: C3DPH4 phase-field damage gradient 4-node linear tetrahedron
C U4: C3DPH8 phase-field damage gradient 8-node linear brick
C
C ======================================================================
C
C Copyright© Gergely Molnar, CNRS, INSA of Lyon, LaMCoS, UMR5259
C
C gergely.molnar@insa-lyon.fr
C http://www.molnar-research.com
C
C This software is a computer program whose purpose is to calculate
C the fracture resistance of elastic-brittle solids in Abaqus under
C static and dynamic conditions.
C
C This software is governed by the CeCILL-B license under French law and
C abiding by the rules of distribution of free software. You can  use, 
C modify and/or redistribute the software under the terms of the CeCILL-B
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info". 
C
C As a counterpart to the access to the source code and  rights to copy,
C modify and redistribute granted by the license, users are provided only
C with a limited warranty  and the software's author,  the holder of the
C economic rights,  and the successive licensors  have only  limited
C liability. 
C
C In this respect, the user's attention is drawn to the risks associated
C with loading,  using,  modifying and/or developing or reproducing the
C software by the user in light of its specific status of free software,
C that may mean  that it is complicated to manipulate,  and  that  also
C therefore means  that it is reserved for developers  and  experienced
C professionals having in-depth computer knowledge. Users are therefore
C encouraged to load and test the software's suitability as regards their
C requirements in conditions enabling the security of their systems and/or 
C data to be ensured and,  more generally, to use and operate it in the 
C same conditions as regards security. 
C
C The fact that you are presently reading this means that you have had
C knowledge of the CeCILL-B
C
C ======================================================================
C Material properties to be given through the input file (*.inp):
C
C For UEL (phase field):
C PROPS(1) = Length scale parameter (lc)
C PROPS(2) = Crack surface energy (gc)
C PROPS(3) = AT switch (1: AT1, 2: AT2)
C
C For UMAT (displacement):
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Length scale parameter (lc)
C PROPS(4) = Crack surface energy (gc)
C PROPS(5) = AT switch (1: AT1, 2: AT2)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C            by 2 - (N_phase+N_UMAT)/2 (to be changed for each model)
C
C NSTVTT - variables to transfer between elements
C NSTVTO - solution dependent variables for the phase-field element
C            (phase, energy history)
C NSTV - overall solution dependent variables (NSTVTO+5), where
C           the additional 4 variables are the: time and iteration numbers
C                          +1 connectivity table
C
C     ==================================================================
C     Comments on solution dependent variables
C     ==================================================================
C
C     Stress/strain element in UMAT form
C                               SDV(1) - tensile elastic energy
C                               SDV(2) - damage phase-field
C                               SDV(3-6) - El. strain at the beginning of the step
C
C     Phase-field element (2*INNODE)+NNODE
C     SDV(1-2)xINNODE:
C                               SDV(1) - damage phase-field
C                               SDV(2) - tensile elastic energy
C     SDV(...): RHS(t_{n-1}) components previous internal
C                                force vector for HHT (dynamic) for each DOF
C
C ======================================================================
C
C   Upload connectivity list for Lagrange multipliers
C
C ======================================================================
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER(N_ELEM=88872,NSTV=16,ZERO=0.D0,HALF=0.5D0)
      DIMENSION TIME(2)
      CHARACTER XOUTDIR*255, XFNAME*200, SKSTEP*10, EXT*255
      CHARACTER DMKNAME*255, FNAMEX*200, SKINC*10
C
      INTEGER CN(8)
      INTEGER I,J,IE,JE, IS,JS
C     
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
C      
      IF (LOP.EQ.0) THEN
C --------------------------------------------------------
C -------- Uploading connectivity table ------------------      
C --------------------------------------------------------
      CALL GETJOBNAME(XFNAME,LXFNAME)     !Input file name
      CALL GETOUTDIR(XOUTDIR,LXOUTDIR)    !output directory
C
      FNAMEX=DMKNAME(XFNAME(1:LXFNAME),XOUTDIR(1:LXOUTDIR),'_connec.dat')
      OPEN(UNIT=96,FILE=FNAMEX,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(7,*)'Opening external file:'
      WRITE(7,*) FNAMEX
      WRITE(7,*)'Reading connectivity structure...'
      WRITE(7,*)'...'
      DO I = 1,N_ELEM
       READ(96,*) CN
       DO J = 1,8
        USRVAR(I,NSTV,J)=CN(J)
       ENDDO
      END DO
      CLOSE(96)
C        
      WRITE(7,*)'Complete.'
      ELSEIF (LOP.EQ.2) THEN
C       
C ----------------------------------------------------------------
C -------- Writing COMMON block to restart file ------------------      
C ----------------------------------------------------------------
       IF ((LRESTART.EQ.1).OR.(LRESTART.EQ.2)) THEN
        CALL GETJOBNAME(XFNAME,LXFNAME)     !Input file name
        CALL GETOUTDIR(XOUTDIR,LXOUTDIR)    !output directory
        WRITE (SKSTEP,'(I10)') KSTEP
        WRITE (SKINC,'(I10)') KINC
        CALL TRIM(SKSTEP,IS) 
        CALL TRIM(SKINC,JS) 
        IE=LEN(SKSTEP)
        JE=LEN(SKINC)
        IF (LRESTART.EQ.1) THEN
         EXT='_'//SKSTEP(IS:IE)//'_'//SKINC(JS:JE)//'_common.res'
        ELSEIF (LRESTART.EQ.2) THEN
         EXT='_'//SKSTEP(IS:IE)//'_common.res'
        ENDIF
        FNAMEX=DMKNAME(XFNAME(1:LXFNAME),XOUTDIR(1:LXOUTDIR),EXT)
        WRITE(7,*)'Writing COMMON block to restart file...'
        OPEN(UNIT=98,FILE=FNAMEX,STATUS='UNKNOWN',FORM='FORMATTED')
        CLOSE (98, STATUS='DELETE')       
        OPEN(UNIT=98,FILE=FNAMEX,STATUS='UNKNOWN',FORM='FORMATTED')
        WRITE(98,*) USRVAR
        WRITE(7,*)'Complete.'
        CLOSE(98)
       ENDIF
      ELSEIF (LOP.EQ.4) THEN 
C ------------------------------------------------------------------
C -------- Reading COMMON block from restart file ------------------      
C ------------------------------------------------------------------
       CALL GETJOBNAME(XFNAME,LXFNAME)     !Input file name
       CALL GETOUTDIR(XOUTDIR,LXOUTDIR)    !output directory
      FNAMEX=DMKNAME(XFNAME(1:LXFNAME-4),XOUTDIR(1:LXOUTDIR),'_common.res')
       WRITE(7,*)'Reading COMMON block from restart file...'
       OPEN(UNIT=98,FILE=FNAMEX,STATUS='UNKNOWN',FORM='FORMATTED')
       READ(98,*) USRVAR
       WRITE(7,*)'Complete.'
       CLOSE(98)
      ENDIF
      RETURN
      END    
C
C ======================================================================
C
C   Phase-field UEL
C
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-12,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 EIGHT=8.D0,TEN=10.D0,
     3 N_ELEM=88872,NSTVTO=2,NSTVTT=9,NSTV=16)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
       INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,NALL,CNT,INNODE

       REAL*8 AINTW(8),XII(8,MCRD),XI(MCRD),dNdxi(NNODE,MCRD),
     1 VJACOB(MCRD,MCRD),dNdx(NNODE,MCRD),VJABOBINV(MCRD,MCRD),
     2 AN(NNODE),BP(MCRD,NNODE),
     2 DP(MCRD),SDV(NSTV),AMASS(NDOFEL,NDOFEL),RHSINI(NDOFEL),
     3 QSML(MCRD),CLOC(NNODE),AMATRXPH(NNODE,NNODE),RHSPH(NNODE),
     4 ULOC(NNODE),DULOC(NNODE),ULANG(NNODE),UOLD(NNODE)

       REAL*8 DTM,THCK,HIST,CLPAR,GCPAR,ENGN,
     1 ENGDG,ENGD,PLSWT,ANISOSWT,REITER,CA,GD,DGDD,DGDD2,
     2 AD,DADD,DADDT,SMUL,ATSWT

C
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
C
C     History variables
        ENGDG=ZERO
        NALL=NSTVTO+NSTVTT
C        
C     ==================================================================
C     Time an iteration variables
C     ==================================================================
       IF (TIME(2).EQ.ZERO) THEN
        TIMEZ=-999.D0
        TIMEZLOC=-99.D0
         DO K2=1,NSTV-1
          DO K3=1,8
           USRVAR(JELEM,K2,K3)=ZERO
          END DO
         END DO
       ELSE
        TIMEZ=USRVAR(JELEM,NALL+1,1)
        TIMEZLOC=TIME(2)-DTIME
       ENDIF
       DTZERO=USRVAR(JELEM,NALL+2,1)
       IF (TIMEZ.LT.TIMEZLOC) THEN
        USRVAR(JELEM,NALL+1,1)=TIMEZLOC
        USRVAR(JELEM,NALL+2,1)=DTIME
        USRVAR(JELEM,NALL+3,1)=ZERO
        USRVAR(JELEM,NALL+4,1)=ZERO
       ELSE
        IF (DTZERO.GT.DTIME*(ONE+TOLER)) THEN
C -----   New correcting iteration   -----
         USRVAR(JELEM,NALL+2,1)=DTIME
         USRVAR(JELEM,NALL+3,1)=USRVAR(JELEM,NALL+3,1)+ONE
         USRVAR(JELEM,NALL+4,1)=ZERO
        ELSE
C -----   New local step   -----
         USRVAR(JELEM,NALL+4,1)=USRVAR(JELEM,NALL+4,1)+ONE
        ENDIF
       ENDIF
C ---- Time integration control
       IF (JELEM.EQ.ONE) THEN
        IF ((LFLAGS(1).EQ.1).OR.(LFLAGS(1).EQ.11)) THEN
         USRVAR(1,NALL+1,2)=ONE
        ELSE
         USRVAR(1,NALL+1,2)=ZERO
        ENDIF
       ENDIF
       REITER=USRVAR(JELEM,NALL+3,1)
       STEPITER=USRVAR(JELEM,NALL+4,1)
C       
C     ==================================================================
C     Material parameters
C     ==================================================================
       CLPAR = PROPS(1)
       GCPAR = PROPS(2)
       ATSWT = PROPS(3)
       THCK = ONE
C     ==================================================================
C     Initial preparations
C     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
       DO K1 = 1, NNODE     
        RHSPH(K1) = ZERO
        DO K2 = 1, NNODE
         AMATRXPH(K2,K1) = ZERO
        END DO
       END DO
       DO K1 = 1, NNODE     
        ULOC(K1) = U((K1-1)*2+1)
        ULANG(K1) = U((K1-1)*2+2)
        DULOC(K1) = DU((K1-1)*2+1,1)
        UOLD(K1) = ULOC(K1) - DULOC(K1)
       END DO
C
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
       IF (JTYPE.EQ.1) THEN
C U1: C2DPH3 phase-field damage gradient 3-node linear
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = 1
        AINTW(1) = HALF
       ELSEIF (JTYPE.EQ.2) THEN
C U2: C2DPH4 phase-field damage gradient 4-node bilinear
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = -ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = 4
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ELSEIF (JTYPE.EQ.3) THEN
C U3: C3DPH4 phase-field damage gradient 4-node linear tetrahedron
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE= 1
        AINTW(1) = ONE/SIX
       ELSEIF (JTYPE.EQ.4) THEN
C U4: C3DPH8 phase-field damage gradient 8-node linear brick
        XII(1,1) = MONE/THREE**HALF
        XII(1,2) = MONE/THREE**HALF
        XII(1,3) = MONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = MONE/THREE**HALF
        XII(2,3) = MONE/THREE**HALF
        XII(3,1) = MONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = MONE/THREE**HALF
        XII(4,1) = ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = MONE/THREE**HALF
        XII(5,1) = MONE/THREE**HALF
        XII(5,2) = MONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = MONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = MONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = ONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
       DO INPT=1,INNODE
C     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTO
          SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
C
C     Local coordinates of the integration point
        DO I=1,MCRD
         XI(I) = XII(INPT,I)
        END DO
C        
C     Shape functions and local derivatives
        IF (JTYPE.EQ.1) THEN
         CALL SHAPEFUNTRI(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.2) THEN
         CALL SHAPEFUNQUAD(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.3) THEN
         CALL SHAPEFUNTET(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.4) THEN
         CALL SHAPEFUNBRICK(AN,dNdxi,XI)
        ENDIF
C        
C     Jacobian
        DO I = 1,MCRD
         DO J = 1,MCRD
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
C        
        DTM = ZERO
        IF (MCRD.EQ.2) THEN
C --------- 2D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        ELSE
C --------- 3D Elements ------------        
        DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)        
        ENDIF 
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        ENDIF
C
C     Inverse of Jacobian
        IF (MCRD.EQ.2) THEN
C --------- 2D Elements ------------        
         VJABOBINV(1,1)=VJACOB(2,2)/DTM
         VJABOBINV(1,2)=-VJACOB(1,2)/DTM
         VJABOBINV(2,1)=-VJACOB(2,1)/DTM
         VJABOBINV(2,2)=VJACOB(1,1)/DTM
        ELSE
C --------- 3D Elements ------------        
         VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,2))/DTM
         VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1    VJACOB(1,3))/DTM
         VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,2))/DTM
         VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,1))/DTM
         VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1    VJACOB(2,1))/DTM      
        ENDIF 
C        
C     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,MCRD
          dNdx(K,I) = ZERO
          DO J = 1,MCRD
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
C
C     Calculating B matrix (B=LN)
       DO INODE=1,NNODE
        DO K1=1,MCRD
         BP(K1,INODE)=dNdx(INODE,K1)
        END DO
       END DO
C
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        PHASE=ZERO
        DPHASE=ZERO
        DO I=1,NNODE
         PHASE=PHASE+AN(I)*ULOC(I)
        END DO
        DO I=1,NNODE
         DPHASE=DPHASE+AN(I)*DULOC(I)
        END DO
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         SDV(1)=PHASE-DPHASE
        ELSE
         SDV(1)=PHASE
        ENDIF
C
C     Gradient
        DO I=1,MCRD
         DP(I)=ZERO
        END DO
        DO I=1,MCRD
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*ULOC(J)
         END DO
        END DO
C 
C   Calculate geometric crack and energy degradation functions
        IF (ATSWT.EQ.ONE) THEN
         CALL AFUNAT1(PHASE,AD,DADD,DADDT)
         CA = EIGHT/THREE ! Normalisation constant 
        ELSEIF (ATSWT.EQ.TWO) THEN
         CALL AFUNAT2(PHASE,AD,DADD,DADDT)
         CA = TWO ! Normalisation constant
        ENDIF
        CALL GFUN(PHASE,GD,DGDD,DGDDT)
C
C     ==================================================================
C     Calculating elastic ENERGY history
C     ==================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         ENGN=USRVAR(JELEM,NSTVTO+1,INPT)
        ELSE
         ENGN=USRVAR(JELEM,2,INPT)
        ENDIF
        HIST=ENGN
        SDV(2)=HIST
C        
C     ==================================================================
C     Calculating fracture energy for history output
C     ==================================================================
C
        ENGD=ZERO
        ENGD=AD/CA/CLPAR*DTM*THCK*GCPAR
C
        DO J=1,MCRD
         ENGD=ENGD+DP(J)*DP(J)*CLPAR*DTM*THCK/CA*GCPAR
        END DO
        ENGDG=ENGDG+ENGD
C
C     ==================================================================
C     Constants for residue and tangent
C     ==================================================================
        QCAP = MONE*DGDD*HIST-GCPAR/CA/CLPAR*DADD
        DQDD = MONE*DGDDT*HIST-GCPAR/CA/CLPAR*DADDT
        DO I = 1,MCRD
         QSML(I) = TWO*DP(I)*GCPAR*CLPAR/CA
        END DO
C       
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
        DO I=1,NNODE
         DO K=1,NNODE
          DO J=1,MCRD
           AMATRXPH(I,K)=AMATRXPH(I,K)+BP(J,I)*BP(J,K)*DTM*
     1      THCK*GCPAR*CLPAR*TWO/CA*AINTW(INPT)
          END DO
          AMATRXPH(I,K)=AMATRXPH(I,K)-AN(I)*AN(K)*DTM*THCK*
     1     AINTW(INPT)*DQDD
         END DO
        END DO
C
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
        DO I=1,NNODE
         DO J=1,MCRD
           RHSPH(I)=RHSPH(I)-BP(J,I)*QSML(J)*AINTW(INPT)*DTM*THCK
         END DO
         RHSPH(I)=RHSPH(I)+AN(I)*AINTW(INPT)*DTM*THCK*QCAP
        END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         IF (LFLAGS(3).EQ.5) THEN
         ELSE
          USRVAR(JELEM,I,INPT)=SVARS(NSTVTO*(INPT-1)+I)
         ENDIF
        END DO
C        
       END DO
C
      ENERGY(7)=ENGDG
C     ==================================================================
C       ---------- END OF CALCULATION ON INTEGRATION POINTS ----------
C     ==================================================================
C       
C     ==================================================================
C     HHT
C     ==================================================================
C
       DO I=1,NNODE
        RHSINI(I)=SVARS(I+NSTVTO*INNODE)
        SVARS(I+NSTVTO*INNODE)=RHSPH(I)
       END DO
       IF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
        PARALPHA=PARAMS(1)
        DO I=1,NNODE
         RHSPH(I)=RHSPH(I)*(ONE+PARALPHA)+RHSINI(I)*PARALPHA
        END DO
        DO I=1,NNODE
         DO K=1,NNODE
           AMATRXPH(I,K)=AMATRXPH(I,K)*(ONE+PARALPHA)
         END DO
        END DO
       ENDIF
C
C     ==================================================================
C     Adding Lagrange multipliers to enforce damage irreversibility
C     ==================================================================
C
      DO I=1,NNODE
       RHS((I-1)*2+1,1)=RHSPH(I)
       DO J=1,NNODE
        AMATRX((I-1)*2+1,(J-1)*2+1)=AMATRXPH(I,J)
       ENDDO
      END DO      
C
      DO I=1,NNODE
       CLOC(I) = USRVAR(JELEM,NSTV,I)
      END DO
C      
      DO I=1,NNODE
       IF (CLOC(I).GT.TOLER) THEN
        IF ((ULANG(I).GT.-TOLER).OR.(DULOC(I).LT.-TOLER)) THEN
         RHS((I-1)*2+1,1)=RHS((I-1)*2+1,1)+ULANG(I)
         RHS((I-1)*2+2,1)=DULOC(I)
         AMATRX((I-1)*2+1,(I-1)*2+2)=MONE
         AMATRX((I-1)*2+2,(I-1)*2+1)=MONE
        ELSE
         AMATRX((I-1)*2+2,(I-1)*2+2)=ONE
        ENDIF
       ELSEIF (CLOC(I).LT.-TOLER) THEN
         AMATRX((I-1)*2+2,(I-1)*2+2)=ONE
       ENDIF
      END DO
C
      RETURN
      END
C
C --------------------------------------------------------------      
C Energy degradation function
C --------------------------------------------------------------      
C
      SUBROUTINE GFUN(PHASE,GD,DGDD,DGDDT)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0)
C      
       GD=(ONE-PHASE)**TWO
       DGDD=MONE*TWO*(ONE-PHASE)
       DGDDT=TWO
C
      RETURN
      END
C
C --------------------------------------------------------------      
C Geometric crack function
C --------------------------------------------------------------      
C
      SUBROUTINE AFUNAT1(PHASE,AD,DADD,DADDT)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,TWO=2.D0,MONE=-1.D0,ONE=1.D0)
C      
       AD=PHASE
       DADD=ONE
       DADDT=ZERO
C
      RETURN
      END
C      
      SUBROUTINE AFUNAT2(PHASE,AD,DADD,DADDT)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,TWO=2.D0,MONE=-1.D0,ONE=1.D0)
C      
       AD=PHASE*PHASE
       DADD=TWO*PHASE
       DADDT=TWO
C
      RETURN
      END
C
C --------------------------------------------------------------      
C Shape functions for U1: 3-node linear
C --------------------------------------------------------------      
      SUBROUTINE SHAPEFUNTRI(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)
C
C     Derivatives of shape functions respect to local ccordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
      RETURN
      END
C
C --------------------------------------------------------------      
C Shape functions U2: 4-node bilinear
C --------------------------------------------------------------      
      SUBROUTINE SHAPEFUNQUAD(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
C
C     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
C
C ---------
C U3: C2DPH6 phase-field damage gradient 6-node quadratic
C U4: C2DPH8 phase-field damage gradient 8-node biquadratic
C
C ---------
C 
C --------------------------------------------------------------      
C Shape functions U5: 4-node linear tetrahedron
C --------------------------------------------------------------   
      SUBROUTINE SHAPEFUNTET(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,3)
      Real*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = ONE-XI(1)-XI(2)-XI(3)
      AN(2) = XI(1)
      AN(3) = XI(2)
      AN(4) = XI(3)
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE
      dNdxi(1,2) =  MONE
      dNdxi(1,3) =  MONE
C
      dNdxi(2,1) =  ONE
      dNdxi(2,2) =  ZERO
      dNdxi(2,3) =  ZERO
C
      dNdxi(3,1) =  ZERO
      dNdxi(3,2) =  ONE
      dNdxi(3,3) =  ZERO
C
      dNdxi(4,1) =  ZERO
      dNdxi(4,2) =  ZERO
      dNdxi(4,3) =  ONE
C
      RETURN
      END
C
C --------------------------------------------------------------      
C Shape functions for U6: 8-node linear brick
C --------------------------------------------------------------      
      SUBROUTINE SHAPEFUNBRICK(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(8),dNdxi(8,3)
      REAL*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
C     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      
      RETURN
      END
C ==============================================================
C ==============================================================
C                       Subroutine UMAT 
C                      Elastic material
C ==============================================================
C ==============================================================
C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0,HALF=0.5,
     1 TEN=10.0,ZERO=0.0, 
     2 DENGMAX=0.5,DTMIN=1.0D-9,
     3 N_ELEM=88872,NSTVTO=2,NSTVTT=9,NSTV=16) 
       DATA NEWTON,TOLER/40,1.D-6/ 
C       
       REAL*8 EPSC(6),CMAT(6,6),EIGV(3),VECTI(3),ALPHAI(3)
C
       REAL*8 EMOD,ENU,CLPAR,GCPAR,ELAMEG,ELAMEL,PLSWT,
     1  ANISOSWT,REITER,STEPITER,ENGN,ENGP,PARK,PHASE,
     2  GFUN,DENG,DENGMAXV,DTIMENEXT,ATSWT
C
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
C 
C ----------------------------------------------------------- 
C          Material properties
C ----------------------------------------------------------- 
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Length scale parameter (lc)
C PROPS(4) = Crack surface energy (gc)
C PROPS(5) = AT switch (1: AT1, 2: AT2)
C
C ----------------------------------------------------------- 
C
       NALL=NSTVTO+NSTVTT
       REITER=USRVAR(NOEL-N_ELEM,NALL+3,1)
       STEPITER=USRVAR(NOEL-N_ELEM,NALL+4,1)
C
C	Elastic properties
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       CLPAR=PROPS(3)
       GCPAR=PROPS(4)
       ATSWT=PROPS(5)
       ELAMEG=EMOD/(TWO*(ONE+ENU))
       ELAMEL=ELAMEG*TWO*ENU/(ONE-TWO*ENU)
       PARK = TOLER
       ANISOSWT=ONE
C
C ======================================================================
C   Recovering damage phase-field from the phase-field element
C ======================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         PHASE=USRVAR(NOEL-N_ELEM,1,NPT)
        ELSE
         PHASE=USRVAR(NOEL-N_ELEM,NSTVTO+2,NPT)
        ENDIF
C
        STATEV(2)=PHASE
        GFUN=(ONE-PHASE)**TWO
C
C
C ======================================================================
C   Strain components
C ======================================================================
       DO K1=1,6
        EPSC(K1)=ZERO
       END DO
       DO K1=1,NTENS
        STRAN(K1)=STRAN(K1)+DSTRAN(K1)
       END DO
C
       IF ((STEPITER.LE.(3-REITER)).AND.(REITER.LT.4)) THEN
C ------- Normal case: the stiffness is refreshed in every iteration
C         based on the actual strain state
         DO K1=1,NDI
          EPSC(K1)=STRAN(K1)
          USRVAR(NOEL-N_ELEM,NSTVTO+2+K1,NPT)=EPSC(K1)
         END DO
         DO K1=1,NTENS-NDI
          EPSC(K1+3)=STRAN(NDI+K1)
          USRVAR(NOEL-N_ELEM,NSTVTO+2+NDI+K1,NPT)=EPSC(K1+3)
         END DO
C     
       ELSEIF ((REITER.GE.4).AND.(STEPITER.EQ.ZERO)) THEN
C ------- If maximum trials are exceeded, the strain state from
C         the original (converged) step is taken
         DO K1=1,NDI
          EPSC(K1)=STRAN(K1)-DSTRAN(K1)
          USRVAR(NOEL-N_ELEM,NSTVTO+2+K1,NPT)=EPSC(K1)
         END DO
         DO K1=1,NTENS-NDI
          EPSC(K1+3)=STRAN(NDI+K1)-DSTRAN(NDI+K1)
          USRVAR(NOEL-N_ELEM,NSTVTO+2+NDI+K1,NPT)=EPSC(K1+3)
         END DO
C
       ELSE
C ------- If the stiffness is not refreshed, the strain is recovered
C         from an older NR iteration
        DO K1=1,NTENS
         EPSC(K1)=USRVAR(NOEL-N_ELEM,NSTVTO+2+K1,NPT)
        END DO
       ENDIF
C
C ======================================================================
C	Stiffness tensor
C ======================================================================
C
C     ==================================================================
C     Calculating degrated elastic stiffness matrix
C     ==================================================================
       DO I=1,6
        DO J=1,6
        CMAT(I,J)=ZERO
        END DO
       END DO
       IF (PHASE.LT.TOLER) THEN
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=CMAT(1,1)
        CMAT(3,3)=CMAT(1,1)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=CMAT(1,2)
        CMAT(1,3)=CMAT(1,2)
        CMAT(3,1)=CMAT(1,3)
        CMAT(2,3)=CMAT(1,2)
        CMAT(3,2)=CMAT(2,3)
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=CMAT(4,4)
        CMAT(6,6)=CMAT(4,4)
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=CMAT(I,J)*GFUN
         END DO
        END DO        
       ELSE
        CALL HMAT(EPSC,CMAT,ENU,EMOD,ANISOSWT,PHASE,PARK)
       ENDIF
C       
        DO K1=1, NTENS
         DO K2=1, NTENS
          DDSDDE(K2, K1)=ZERO
         END DO
        END DO
        DO I=1,NDI
         DO J=1,NDI
          DDSDDE(I,J)=CMAT(I,J)
         END DO
        END DO
        DO I=1,NTENS-NDI
         DO J=1,NTENS
          DDSDDE(J,I+NDI)=CMAT(J,I+3)
          DDSDDE(I+NDI,J)=CMAT(I+3,J)
         END DO
        END DO
C
C ======================================================================
C	Calculate Stresses
C ======================================================================
C
       DO K1=1, NTENS
        STRESS(K1)=ZERO
        DO K2=1, NTENS
         STRESS(K1)=STRESS(K1)+DDSDDE(K1, K2)*STRAN(K2)
        END DO
       END DO
C
C ======================================================================
C	Energies
C ======================================================================
C
C ======================================================================
C	Energies
C ======================================================================
       DO K1=1,6
        EPSC(K1)=ZERO
       END DO
       DO K1=1,NDI
        EPSC(K1)=STRAN(K1)
       END DO
       DO K1=1,NTENS-NDI
        EPSC(K1+3)=STRAN(NDI+K1)
       END DO       
       CALL EIGOWN(EPSC,EIGV,ALPHA,ALPHAI,VECTI)
C       
       ENGP=(ELAMEL*(ALPHA*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1  TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2  ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
       ENGN=(ELAMEL*((ONE-ALPHA)*(EIGV(1)+EIGV(2)+EIGV(3)))**
     1  TWO)/TWO+ELAMEG*((EIGV(1)*(ONE-ALPHAI(1)))**TWO+(EIGV(2)*
     2  (ONE-ALPHAI(2)))**TWO+(EIGV(3)*(ONE-ALPHAI(3)))**TWO)
C
       SSE=ENGP*GFUN+ENGN
       DENG=ENGP-STATEV(1)
       STATEV(1)=ENGP
       USRVAR(NOEL-N_ELEM,NSTVTO+1,NPT)=ENGP
       USRVAR(NOEL-N_ELEM,NSTVTO+2,NPT)=PHASE
C       
C ======================================================================
C	Time step control
C ======================================================================
C 
        IF (USRVAR(1,NALL+1,2).GT.HALF) THEN
         IF (ATSWT.EQ.ONE) THEN
          DENGMAXV=3.0*GCPAR/8.0/CLPAR*DENGMAX
         ELSE
          DENGMAXV=GCPAR/2.0/CLPAR*DENGMAX
         ENDIF
C --------- If automatic time integration is used
         IF ((DENG.GT.DENGMAXV).AND.(REITER.LT.4)) THEN
C --------- If energy criterion is violated (dE>dE_max)
          IF (DTIME.GT.DTMIN*(ONE+TOLER)) THEN
C --------- If time-step still can be reduced
           PNEWDT=DENGMAXV/DENG/TEN
           DTIMENEXT=PNEWDT*DTIME
           IF (DTIMENEXT.LT.(DTMIN*(ONE+TOLER))) THEN
            PNEWDT=PNEWDT*DTMIN/DTIMENEXT
           ENDIF
          ELSE
C --------- If time-step is already too small
          PNEWDT=ONE
          ENDIF
         ENDIF
        ENDIF
C       
       RETURN
       END      
C ==============================================================
C ==============================================================
C
C       Miscellaneous code 
C           (eigenvalue calculator,
C              3rd order polynomial solver, etc.)
C
C ==============================================================
C ==============================================================
C       Stiffness matrix based on asymmetric energy degradation
C ==============================================================
C
       SUBROUTINE HMAT(EPSI,CMAT,ENU,EMOD,ANISOSWT,PHASE,PARK)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12,TOLD=1.0D-7)
       INTEGER I, J, K, NDIM
       REAL*8 EPS(6),VECTI(3),CLMAT(6,6),EIGV(6),ALPHAI(3),
     1  CMAT(6,6),EPSI(6),EIGVEC(3,3),EIGVAL(3,3),TEPS(6,6),
     2  TSIGINV(6,6)
C     
       REAL*8 VMAXE, ALPHA, ANISOSWT, ENU, EMOD,EMU
     1 PHASE, CNT
C
C Rescaling the strain vector
       VMAXE=MAXVAL(ABS(EPSI))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPSI(K1)/VMAXE
        END DO
       ELSE
        DO K1=1,6
         EPS(K1)=EPSI(K1)
        END DO
       ENDIF
C	Initialising CMAT
        DO I=1,6
         EIGV(I)=ZERO
         DO J=1,6
          CMAT(I,J)=ZERO
          CLMAT(I,J)=ZERO
         END DO
        END DO
C
C Obtaining eigenvalues and eigenvectors
        CALL JACOBIEIG(EPS,EIGVEC,EIGVAL)
        DO I=1,3
         EIGV(I)=EIGVAL(I,I)
        END DO
C        
C ************ Starting spectral decomposition ************************        
       IF ((MINVAL(EIGV).GE.-TOLER).OR.(ANISOSWT.LT.TOLER)) THEN
C All eigenvalues are in tension, therefore the stiffness matrix can be
C     degraded without the decomposition
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=CMAT(I,J)*((ONE-PHASE)**TWO+PARK)
         END DO
        END DO
C
       ELSEIF (MAXVAL(EIGV).LT.-TOLER) THEN
C All eigenvalues are in compression, therefore the stiffness
C     matrix does not need to be degraded 
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
C
       ELSE
C Degradation switches       
        ALPHA=ZERO
        IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
         ALPHA=ONE
        ENDIF
        DO K1=1,3
         ALPHAI(K1)=ZERO
         IF (EIGV(K1).GT.TOLER) THEN
          ALPHAI(K1)=ONE
         ENDIF
        END DO
C
C Elementary stiffness matrix in the principal directions
        GAMMA=ENU/(ONE-TWO*ENU)
        CLMAT(1,1)=((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,2)=((ONE-ALPHAI(2)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(3,3)=((ONE-ALPHAI(3)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,2)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,1)=CLMAT(1,2)
        CLMAT(3,1)=CLMAT(1,3)
        CLMAT(3,2)=CLMAT(2,3)
        DO I=1,3
         DO J=1,3
          CLMAT(I,J)=CLMAT(I,J)*EMOD/(ONE+ENU)
         END DO
        END DO
        EMU=EMOD/(TWO*(ONE+ENU))
        CLMAT(4,4)=(ABS(EIGV(1))*((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+
     1    ABS(EIGV(2))*((ONE-ALPHAI(2)*PHASE)**TWO+PARK))/
     2   (ABS(EIGV(1))+ABS(EIGV(2)))*EMU
        CLMAT(5,5)=(ABS(EIGV(1))*((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+
     1    ABS(EIGV(3))*((ONE-ALPHAI(3)*PHASE)**TWO+PARK))/
     2   (ABS(EIGV(1))+ABS(EIGV(3)))*EMU
        CLMAT(6,6)=(ABS(EIGV(3))*((ONE-ALPHAI(3)*PHASE)**TWO+PARK)+
     1    ABS(EIGV(2))*((ONE-ALPHAI(2)*PHASE)**TWO+PARK))/
     2   (ABS(EIGV(3))+ABS(EIGV(2)))*EMU
C
C   Rotating the stiffness matrix in the original system
C
        CALL ROTMATS(EIGVEC,TEPS,TSIGINV)
        DO I=1,6
            DO J=1,6
                DO K=1,6
                    DO L=1,6
                        CMAT(I,J)=CMAT(I,J)+TSIGINV(I,K)*CLMAT(K,L)*TEPS(L,J);
                    END DO
                END DO
            END DO
        END DO
       ENDIF
C       
       RETURN
       END
C
      SUBROUTINE ROTMATS(EIGVEC,TEPS,TSIGINV)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(TWO=2.D0)       
       REAL*8 EIGVEC(3,3), TEPS(6,6), TSIGINV(6,6)
C
C Rotation matrix for the strain tensor in Voigt notation
       TEPS(1,1)=EIGVEC(1,1)**TWO
       TEPS(1,2)=EIGVEC(2,1)**TWO
       TEPS(1,3)=EIGVEC(3,1)**TWO
       TEPS(1,4)=EIGVEC(2,1)*EIGVEC(1,1)
       TEPS(1,5)=EIGVEC(3,1)*EIGVEC(1,1)
       TEPS(1,6)=EIGVEC(2,1)*EIGVEC(3,1)
       TEPS(2,1)=EIGVEC(1,2)**TWO
       TEPS(2,2)=EIGVEC(2,2)**TWO
       TEPS(2,3)=EIGVEC(3,2)**TWO
       TEPS(2,4)=EIGVEC(2,2)*EIGVEC(1,2)
       TEPS(2,5)=EIGVEC(3,2)*EIGVEC(1,2)
       TEPS(2,6)=EIGVEC(2,2)*EIGVEC(3,2)
       TEPS(3,1)=EIGVEC(1,3)**TWO
       TEPS(3,2)=EIGVEC(2,3)**TWO
       TEPS(3,3)=EIGVEC(3,3)**TWO
       TEPS(3,4)=EIGVEC(2,3)*EIGVEC(1,3)
       TEPS(3,5)=EIGVEC(3,3)*EIGVEC(1,3)
       TEPS(3,6)=EIGVEC(2,3)*EIGVEC(3,3)
       TEPS(4,1)=TWO*EIGVEC(1,2)*EIGVEC(1,1)
       TEPS(4,2)=TWO*EIGVEC(2,2)*EIGVEC(2,1)
       TEPS(4,3)=TWO*EIGVEC(3,2)*EIGVEC(3,1)
       TEPS(4,4)=EIGVEC(2,2)*EIGVEC(1,1)+EIGVEC(1,2)*EIGVEC(2,1)
       TEPS(4,5)=EIGVEC(3,1)*EIGVEC(1,2)+EIGVEC(1,1)*EIGVEC(3,2)
       TEPS(4,6)=EIGVEC(2,2)*EIGVEC(3,1)+EIGVEC(3,2)*EIGVEC(2,1)
       TEPS(5,1)=TWO*EIGVEC(1,3)*EIGVEC(1,1)
       TEPS(5,2)=TWO*EIGVEC(2,3)*EIGVEC(2,1)
       TEPS(5,3)=TWO*EIGVEC(3,3)*EIGVEC(3,1)
       TEPS(5,4)=EIGVEC(1,1)*EIGVEC(2,3)+EIGVEC(2,1)*EIGVEC(1,3)
       TEPS(5,5)=EIGVEC(3,1)*EIGVEC(1,3)+EIGVEC(1,1)*EIGVEC(3,3)
       TEPS(5,6)=EIGVEC(2,3)*EIGVEC(3,1)+EIGVEC(3,3)*EIGVEC(2,1)
       TEPS(6,1)=TWO*EIGVEC(1,3)*EIGVEC(1,2)
       TEPS(6,2)=TWO*EIGVEC(2,3)*EIGVEC(2,2)
       TEPS(6,3)=TWO*EIGVEC(3,3)*EIGVEC(3,2)
       TEPS(6,4)=EIGVEC(2,2)*EIGVEC(1,3)+EIGVEC(2,3)*EIGVEC(1,2)
       TEPS(6,5)=EIGVEC(1,3)*EIGVEC(3,2)+EIGVEC(1,2)*EIGVEC(3,3)
       TEPS(6,6)=EIGVEC(2,3)*EIGVEC(3,2)+EIGVEC(3,3)*EIGVEC(2,2)
C
C Inverse of the rotation matrix for the stress tensor in Voigt notation
       TSIGINV(1,1)=EIGVEC(1,1)**TWO
       TSIGINV(1,2)=EIGVEC(1,2)**TWO
       TSIGINV(1,3)=EIGVEC(1,3)**TWO
       TSIGINV(1,4)=TWO*EIGVEC(1,2)*EIGVEC(1,1)
       TSIGINV(1,5)=TWO*EIGVEC(1,3)*EIGVEC(1,1)
       TSIGINV(1,6)=TWO*EIGVEC(1,2)*EIGVEC(1,3)
       TSIGINV(2,1)=EIGVEC(2,1)**TWO
       TSIGINV(2,2)=EIGVEC(2,2)**TWO
       TSIGINV(2,3)=EIGVEC(2,3)**TWO
       TSIGINV(2,4)=TWO*EIGVEC(2,2)*EIGVEC(2,1)
       TSIGINV(2,5)=TWO*EIGVEC(2,3)*EIGVEC(2,1)
       TSIGINV(2,6)=TWO*EIGVEC(2,2)*EIGVEC(2,3)
       TSIGINV(3,1)=EIGVEC(3,1)**TWO
       TSIGINV(3,2)=EIGVEC(3,2)**TWO
       TSIGINV(3,3)=EIGVEC(3,3)**TWO
       TSIGINV(3,4)=TWO*EIGVEC(3,2)*EIGVEC(3,1)
       TSIGINV(3,5)=TWO*EIGVEC(3,3)*EIGVEC(3,1)
       TSIGINV(3,6)=TWO*EIGVEC(3,2)*EIGVEC(3,3)
       TSIGINV(4,1)=EIGVEC(2,1)*EIGVEC(1,1)
       TSIGINV(4,2)=EIGVEC(2,2)*EIGVEC(1,2)
       TSIGINV(4,3)=EIGVEC(2,3)*EIGVEC(1,3)
       TSIGINV(4,4)=EIGVEC(2,2)*EIGVEC(1,1)+EIGVEC(2,1)*EIGVEC(1,2)
       TSIGINV(4,5)=EIGVEC(2,1)*EIGVEC(1,3)+EIGVEC(2,3)*EIGVEC(1,1)
       TSIGINV(4,6)=EIGVEC(2,3)*EIGVEC(1,2)+EIGVEC(2,2)*EIGVEC(1,3)
       TSIGINV(5,1)=EIGVEC(3,1)*EIGVEC(1,1)
       TSIGINV(5,2)=EIGVEC(3,2)*EIGVEC(1,2)
       TSIGINV(5,3)=EIGVEC(3,3)*EIGVEC(1,3)
       TSIGINV(5,4)=EIGVEC(3,2)*EIGVEC(1,1)+EIGVEC(3,1)*EIGVEC(1,2)
       TSIGINV(5,5)=EIGVEC(3,1)*EIGVEC(1,3)+EIGVEC(3,3)*EIGVEC(1,1)
       TSIGINV(5,6)=EIGVEC(3,3)*EIGVEC(1,2)+EIGVEC(3,2)*EIGVEC(1,3)
       TSIGINV(6,1)=EIGVEC(3,1)*EIGVEC(2,1)
       TSIGINV(6,2)=EIGVEC(3,2)*EIGVEC(2,2)
       TSIGINV(6,3)=EIGVEC(3,3)*EIGVEC(2,3)
       TSIGINV(6,4)=EIGVEC(2,2)*EIGVEC(3,1)+EIGVEC(2,1)*EIGVEC(3,2)
       TSIGINV(6,5)=EIGVEC(3,1)*EIGVEC(2,3)+EIGVEC(3,3)*EIGVEC(2,1)
       TSIGINV(6,6)=EIGVEC(2,3)*EIGVEC(3,2)+EIGVEC(2,2)*EIGVEC(3,3)       
       RETURN
       END
C
C ==============================================================
C Eigenstrains from Voigt notation (3rd order polynomial)
C ==============================================================
C
      SUBROUTINE EIGOWN(EPS,EIGV,ALPHA,ALPHAI,VECTI)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12)
       INTEGER I, J, K
       REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
       REAL*8 PC, QC, ALPHA, DISC, PI, CNT
C
       PI=FOUR*ATAN(ONE)
C Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
C    
C   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
C
C   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
C
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
C       
       ALPHA=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHA=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
C
C    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
C    
       RETURN
       END
C
C
C ==============================================================
C   Jacobi eigenvector and eigenvalue in case of instability
C ==============================================================
C
      SUBROUTINE JACOBIEIG(EPS,EIGVEC,EIGVAL)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12)
       INTEGER I, J, K, N, CNT
       REAL*8 EPS(6),EIGVEC(3,3),EIGVAL(3,3)
       REAL*8 B2, BAR, BETA, COEFF, S, C, CS, SC, VAG
C
       DO I=1,3
        DO J=1,3
         EIGVEC(I,J) = ZERO
         EIGVAL(I,J) = ZERO
        END DO
       END DO
C       
       DO I=1,3
         EIGVEC(I,I) = ONE
       END DO
C       
       EIGVAL(1,1) = EPS(1)
       EIGVAL(2,2) = EPS(2)
       EIGVAL(3,3) = EPS(3)
       EIGVAL(1,2) = EPS(4)/TWO
       EIGVAL(1,3) = EPS(5)/TWO
       EIGVAL(2,3) = EPS(6)/TWO
       EIGVAL(2,1) = EIGVAL(1,2)
       EIGVAL(3,1) = EIGVAL(1,3)
       EIGVAL(3,2) = EIGVAL(2,3)
C     
       B2 = ZERO
       DO I=1,3
         DO J=1,3
           IF (I.NE.J) THEN
               B2 = B2 + EIGVAL(I,J)**TWO
           ENDIF
         END DO
       END DO
       BAR = B2/(THREE*THREE)/TWO
C      
       CNT=ONE
       DO WHILE ((B2.GT.TOLER).AND.(CNT.LT.CNTM))
          DO I=1,2
            DO J=I+1,3
              IF ((EIGVAL(J,I)**TWO-BAR).LE.TOLER) THEN
              ELSE
               B2 = B2 - TWO*EIGVAL(J,I)**TWO
               BAR = HALF*B2/(THREE*THREE)
               BETA = (EIGVAL(J,J)-EIGVAL(I,I))/(TWO*EIGVAL(J,I))
               COEFF = HALF*BETA/SQRT(ONE+BETA**TWO)
               S = SQRT(MAX(HALF+COEFF,ZERO))
               C = SQRT(MAX(HALF-COEFF,ZERO))
               DO K=1,3
                CS =  C*EIGVAL(I,K)+S*EIGVAL(J,K)
                SC = -S*EIGVAL(I,K)+C*EIGVAL(J,K)
                EIGVAL(I,K) = CS
                EIGVAL(J,K) = SC
               END DO
               DO K=1,3
                CS =  C*EIGVAL(K,I)+S*EIGVAL(K,J)
                SC = -S*EIGVAL(K,I)+C*EIGVAL(K,J)
                EIGVAL(K,I) = CS
                EIGVAL(K,J) = SC
                CS =  C*EIGVEC(K,I)+S*EIGVEC(K,J)
                SC = -S*EIGVEC(K,I)+C*EIGVEC(K,J)
                EIGVEC(K,I) = CS
                EIGVEC(K,J) = SC
               END DO
              ENDIF
            END DO
          END DO
       CNT=CNT+ONE
       END DO
C
C ---------- Sorting eigenvalues and vectors ---------------
        VAG=ZERO
        IF (EIGVAL(1,1).GT.EIGVAL(2,2)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(2,2)
            EIGVAL(2,2)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,2)
                EIGVEC(I,2)=VAG
            END DO
        ENDIF
        IF (EIGVAL(1,1).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF 
        IF (EIGVAL(2,2).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(2,2)
            EIGVAL(2,2)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,2)
                EIGVEC(I,2)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF        
       RETURN
       END    
C
C  Compose a filename  directory/jobname.exten

      CHARACTER*(*) FUNCTION DMKNAME(FNAME,DNAME,EXTEN)
C
      CHARACTER*(*) FNAME,DNAME,EXTEN
C     FNAME  I   JOBNAME
C     DNAME  I   DIRECTORY
C     EXTEN  I   EXTENSION
C     DMKNAME O DIRECTORY/JOBNAME.EXTEN
C
      LTOT = LEN(FNAME)
      LF = 0
      DO K1 = LTOT,2,-1
        IF (LF.EQ.0.AND.FNAME(K1:K1).NE.' ')  LF = K1
      END DO
C 
      LTOT = LEN(DNAME)
      LD = 0
      DO K1 = LTOT,2,-1
        IF (LD.EQ.0.AND.DNAME(K1:K1).NE.' ')  LD = K1
      END DO
C
      LTOT = LEN(EXTEN)
      LE = 0
      DO K1 = LTOT,2,-1
        IF (LE.EQ.0.AND.EXTEN(K1:K1).NE.' ')  LE = K1
      END DO
C      
      IF ((LF + LD + LE) .LE. LEN(DMKNAME)) THEN
        DMKNAME = DNAME(1:LD)//'/'//FNAME(1:LF)
        LTOT = LD + LF + 1
        IF ( LE.GT.0) THEN
           DMKNAME = DMKNAME(1:LTOT)//EXTEN(1:LE)
        END IF
      END IF
C
      RETURN
      END

      SUBROUTINE TRIM(STR,LSTR) 
C
      CHARACTER*(*) STR
      INTEGER I, LSTR
C     STR string
C     TRIM: first character which is not empty
        LSTR=0
        DO I = 1,LEN(STR)
         IF (STR(I:I).NE.' ') THEN
          LSTR=I
          EXIT
         ENDIF 
        END DO
C
      RETURN
      END

      
