      PROGRAM calcNyeNorm

      IMPLICIT NONE

      INCLUDE 'PARDIS.H'
      

      REAL(8):: G0XYZ(MAXCRD)

      INTEGER :: IJK(MNELX)
      INTEGER :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      INTEGER :: NX,NELX,NEQ,NODE,MDIM,NGAUSS,NDF

      real(8), allocatable:: elemFp(:,:)
      real(8) :: Ep_Equiv, glStrain(6)

      real(8) :: Fp_t              ! time of the imported solution
      integer :: Fp_N              ! step number of the imported solution
      integer, parameter :: nDOF_Fp = 9
      
   ! some dummy variables
      integer :: iElem,nElem,iNode,iGloNode,iLocNode,nodeIdx
      integer :: istep,nstep,nstep_max,nstepInp,val_ind
      integer :: I,II,III,J,JJ,JJJ,K,iDim

      logical :: fileExists, successReadNodeConnect
      logical :: commentsInData

      integer :: error
      
      
      character(len=30) :: meshName, fileName, solnFileName
      
      !------------------------------------------------------------------------------      
      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      
      ! set element defaults
      MDIM=3
      NDF=3
      NGAUSS=1
      NODE=4
      NEQ=NX*NDF

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      allocate(elemFp(nDOF_Fp,NELX))
      
      write(*,*) 'file name? (e.g. fp)'
      read(*,*) fileName
      
      open(104,file=trim(fileName)//'.out')   ! import
      open(105,file='eqvPlStrain.out')   ! export
      open(106,file='creep_calc.out')   ! export
      
      write(*,*) 'comments in data file? (T/F)'
      read(*,*) commentsInData      
      
      error = 0
      nstep = 0
      do while (error==0)
 
         ! import Fp for this time step
         error = 0
         if (commentsInData) then
            read(104,*,iostat=error) Fp_N, Fp_t
            if(error.NE.0) then   !end of input file      
               write(*,*) 'End of input file.'
               exit 
            endif 
         endif

         do iElem = 1, NELX      
            read(104,*,iostat=error) elemFp(:,iElem)
            if(error.NE.0) exit
         enddo
         if(error.NE.0) then   
            write(*,*) 'Corrupt input file: ',trim(fileName)
            exit
         endif
         
         ! EXPORT equivalent plastic strain
         if (commentsInData) &
            write(105,*,iostat=error) Fp_N, Fp_t
         do iElem = 1, NELX      
            CALL calcE_Equiv(elemFp(:,iElem),Ep_Equiv)
            write(105,*,iostat=error) Ep_Equiv
         enddo
         
         ! EXPORT green lagrange strain
         if (commentsInData) &
            write(106,*,iostat=error) Fp_N, Fp_t
         do iElem = 1, NELX      
            CALL calcGLStrain(elemFp(:,iElem),glStrain(:))
            write(106,*,iostat=error) glStrain(1:6)
         enddo
         
         nstep = nstep + 1
         if (commentsInData) then
            write(*,*) 'time step processed:', Fp_N, Fp_t
         else
            write(*,*) 'STEP DONE:',nstep
         endif
      enddo
      
      deallocate(elemFp)
      
      write(*,*) 'files created: eqvPlStrain.out, creep_calc.out'
      
      close(104)      ! close output file for Fp
      close(105)      ! close output file for equivalent plastic strain
      close(106)      ! close output file for plastic GL strain

      
      END



      SUBROUTINE READSTR(LR)
      IMPLICIT REAL*8(A-H,O-Z)
      !old version:
      !DIMENSION FLAG(20)
      !READ(LR,505) (FLAG(I),I=1,20)
!505   FORMAT(20A4)

! **** modified by deniz
      ! this code skips 1 line of text in the input file.
      character(150)::dummyString
      read(LR,'(A)') dummyString

      RETURN
      END
!*************************************************

      SUBROUTINE SHAPE_TET4(XEL,IPT,SHP,DVOL,XC)
      IMPLICIT NONE

      REAL(8):: SHP1(3,4),SHP(3,4),XJAC(3,3),XJACINV(3,3)
      REAL(8):: XEL(3,4),XC(3)
      REAL(8):: X(4),Y(4),Z(4),XX(3),YY(3),ZZ(3)
      REAL(8):: DVOL
      INTEGER:: I,J, INODE, IPT

      DO I=1,3
      DO J=1,4
      SHP1(I,J)=0.D0
      SHP(I,J)=0.D0
      ENDDO
      ENDDO

      SHP1(1,1)=1.D0

      SHP1(2,2)=1.D0

      SHP1(3,3)=1.D0


      SHP1(1,4)=-1.D0

      SHP1(2,4)=-1.D0

      SHP1(3,4)=-1.D0


   !		WRITE(*,*) 'XEL   IS:  ', XEL

      DO I=1,4
         X(I)=XEL(1,I)
         Y(I)=XEL(2,I)
         Z(I)=XEL(3,I)
      ENDDO

      DO I=1,3
         XX(I)=XEL(1,I)-XEL(1,4)
         YY(I)=XEL(2,I)-XEL(2,4)
         ZZ(I)=XEL(3,I)-XEL(3,4)
      ENDDO

      XJAC(1:3,1)=XX(1:3)
      XJAC(1:3,2)=YY(1:3)
      XJAC(1:3,3)=ZZ(1:3)



      CALL MATINV3_UEL(XJAC,XJACINV,DVOL)
      IF(DVOL.LT.1.D-16)THEN
         WRITE(*,*)'NEGATIVE VOLUME'
         STOP
      ENDIF


      DVOL=1.D0/6.D0*DABS(DVOL)

      DO INODE=1,4
         SHP(:,INODE)=MATMUL(XJACINV,SHP1(:,INODE))
      ENDDO



      RETURN
      END


!------------------------------------------------
      SUBROUTINE MATINV3_UEL(A,AI,DET)
      IMPLICIT NONE
      REAL(8), intent(in) :: A(3,3)
      REAL(8), intent(out):: AI(3,3)
      REAL(8), intent(out):: DET

      DET=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)  &
            *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)  &
           *A(1,3)*A(2,2))

      AI(1,1) =  ( A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AI(1,2) = -( A(1,2)*A(3,3)-A(1,3)*A(3,2))/DET
      AI(1,3) = -(-A(1,2)*A(2,3)+A(1,3)*A(2,2))/DET
      AI(2,1) = -( A(2,1)*A(3,3)-A(2,3)*A(3,1))/DET
      AI(2,2) =  ( A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AI(2,3) = -( A(1,1)*A(2,3)-A(1,3)*A(2,1))/DET
      AI(3,1) =  ( A(2,1)*A(3,2)-A(2,2)*A(3,1))/DET
      AI(3,2) = -( A(1,1)*A(3,2)-A(1,2)*A(3,1))/DET
      AI(3,3) =  ( A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
      RETURN
      END

      subroutine READGM(G0XYZ,IJK,NX,NELX,meshName)
      implicit none
      include 'PARDIS.H'
      real(8), intent(out) :: G0XYZ(MAXCRD)
      integer, intent(out) :: IJK(MNELX)
      integer, intent(out) :: NX, NELX
      character*(*), intent(out) :: meshName
      !locals
      integer :: IJKNEL(MAXELN)
      real(8) :: xc(3)
      integer :: MDIM,NDFADD,NDF,NNODE,NGAUSS
      character(len=50) :: fileName
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ
      
      open(201,file='loadtime.inp')
      call READSTR(201)
      read(201,*) icyc
      call READSTR(201)
      read(201,*) meshName
      close(201)
      
      len=len_trim(meshName)
      fileName=meshName(1:len)//'.inp'

      open(unit=LR,file=fileName)

      call READSTR(LR)
      read(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD     
      
      call READSTR(LR)
      read(LR,*) NX
      NEQ=NDF*NX
      
      IF(NDF.GT.NDOFELX)THEN
         WRITE(*,*)'INSUFFICIENT MEM-NDOFELX'
         STOP
      ENDIF

      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         write(*,*) 'Increase MAXNODE to ', NX
         WRITE(*,*)' MAXNODE is',MAXNODE
         STOP
      ENDIF
      NEQ=NDF*NX
      IF(NDF.GT.MDOFX)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY-MDOFX'
         STOP
      ENDIF
      IF(NX*MDIM.GT.MAXCRD)THEN
         WRITE(*,*)"INSUFFICIENT MEM -MAXCRD"
         STOP
      ENDIF

      IF(MDIM.GT.MAXDIM)THEN
         WRITE(*,*)'INSUFFICIENT MEM-MAXDIM'
         STOP
      ENDIF
      
      ! READ INITIAL POSITIONS
      do iNode=1,NX
         read(LR,*) nnum,(xc(ii),ii=1,MDIM)
         do ii=1,MDIM
            G0XYZ(MDIM*(nnum-1)+ii)=xc(ii)
         enddo
      enddo
      
      call READSTR(LR)
      if(IBRICKORTET.EQ.1)then
         NNODE=8
         NGAUSS=8
      else
         NGAUSS=1
         NNODE=4
      endif

      read(LR,*) NELX

      ! READ IN THE IJK/CONNECTIVITY MATRIX
      jj=0
      do iElem=1,NELX
         read(LR,*) NEL,(IJKNEL(j),j=1,NNODE)
         do j=1,NNODE
            jj=jj+1
            IJK(jj)=IJKNEL(j)
         enddo
      enddo
      
      close(LR)

      return
      end
      
      function factorial(n) result(n_fact)
         integer, intent(in) :: n
         integer:: n_fact
         n_fact = 1
         do i=1,n
            n_fact=n_fact*i
         enddo
      end function factorial
     
      SUBROUTINE printIntMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      integer, intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN(TRIM(matName)).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE printVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      real(8), intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN(TRIM(vecName)).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE
      
      SUBROUTINE printIntVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      integer, intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN(TRIM(vecName)).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE

      SUBROUTINE printMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      real(8), intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN(TRIM(matName)).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
      
      
      SUBROUTINE readGrainTexture(grain_rot,grain_euler,elementGrains,nGrains,NELX)
      implicit none
      
      real(8), intent(out):: grain_rot(3,3,nGrains), grain_euler(3,nGrains)
      integer, intent(in) :: elementGrains(NELX)
      integer, intent(in) :: nGrains,NELX
      
      integer :: iElem,iGrain,nGrains_in, iGrain_in
      real(8) :: euler(3)

      OPEN(UNIT=201,FILE='grainTexture.inp')
      read(201,*) nGrains_in
      
      if(nGrains_in /= nGrains) then
         write(*,*) 'nGrains in grainTexture is inconsisten with the rest of the input'
         stop
      endif

      DO iGrain=1,nGrains
         read(201,*) iGrain_in,grain_euler(1:3,iGrain)
         call eulerRot(grain_euler(:,iGrain),grain_rot(:,:,iGrain))
      ENDDO
      
      CLOSE(201)
      END SUBROUTINE
      



! output, rot,:
!   an active transformation matrix corresponding to the actual rotation due to euler angles {v}_i = TLGT_ij {v0}_j
!   as a (passive) coordinate transformation matrix, it transforms coordinates from crystal to global spatial frame. {v}_i = TLGT_ij {v_cry}_j
! input - euler angles
!     phi   - euler(1)
!     theta - euler(2)
!     omega - euler(3)
!---------------------------------------------------
      subroutine eulerRot(euler,rot)   ! restored to original
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),rot(3,3)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)

  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct

      rot = tlg               ! Deka
      !rot = TRANSPOSE(tlg)   ! Kourosh

      return
      end   

      SUBROUTINE getNc_fromEuler(nx,ny,nz,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: nx,ny,nz
      real(8), intent(in) :: phi1,theta,phi2
      real(8), parameter :: PI=3.14159265359D0

      nx = sin(theta)*sin(phi1)
      ny = -sin(theta)*cos(phi1)
      nz = cos(theta)

      END SUBROUTINE
      
      SUBROUTINE findDispIndex(this_N,disp_ind,disp_Ind_last,nodalPositions_N,nStep,error)
      implicit none
      integer, intent(in) :: this_N
      integer, intent(out) :: disp_ind
      integer, intent(inout) :: disp_Ind_last
      integer, intent(in) :: nodalPositions_N(nStep),nStep
      logical, intent(out) :: error

      error = .FALSE.
      
      disp_ind = disp_Ind_last
      do while(nodalPositions_N(disp_ind).NE.this_N)
         disp_ind = disp_ind + 1
         if(disp_ind.GT.nStep) then
            error=.TRUE.
            return
         endif
      enddo
      
      disp_Ind_last = disp_ind
      
      END SUBROUTINE
      
!***********************************************************************
!C     READS LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE readNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX,NX,success)
      implicit none
      INCLUDE 'PARDIS.H'
      character(len=*), intent(in) :: strConnectFile      
      logical, intent(out) :: success
      integer, intent(out):: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      integer, intent(in) :: NX
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, nNodes, I, lastIndex

      success = .false.
      
      inquire(file=trim(strConnectFile),exist=fileExists)
      if(.NOT.fileExists) return
      
      open(716,file=trim(strConnectFile))
      read(716,*,iostat=error) nNodes
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      lastIndex = 1
      do iNode=1,nNodes
         read(716,*,iostat=error) nElem
         read(716,*,iostat=error) (IELEMNO(I),I=lastIndex,lastIndex+nElem-1)
         lastIndex = lastIndex + nElem
         NINDX(iNode) = lastIndex
      enddo
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      
      success = .true.
      close(716)
      
      END SUBROUTINE

!***********************************************************************
!C     CREATES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE
!        modified by Kourosh for faster calculation
!***********************************************************************
      SUBROUTINE NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
      implicit none
      integer, intent(in) :: IJK(MNELX),NX,NELX,MNELX,MAXNODELEM
      integer, intent(out) :: IELEMNO(MAXNODELEM*NX),NINDX(NX)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************
      integer, allocatable :: IELEMNO_TMP(:,:),NINDX_TMP(:)
      integer :: nodeID, iNode, elemID, iElem, KOUNT

      ALLOCATE(IELEMNO_TMP(NX,MAXNODELEM),NINDX_TMP(NX))
      NINDX_TMP=0
      DO elemID=1,NELX
         DO iNode=1,4
            nodeID = IJK((elemID-1)*4+iNode)
            NINDX_TMP(nodeID)=NINDX_TMP(nodeID)+1
            IELEMNO_TMP(nodeID,NINDX_TMP(nodeID))=elemID
         ENDDO
      ENDDO
      
      KOUNT=1
      DO nodeID=1,NX
         DO iElem=1,NINDX_TMP(nodeID)
            IELEMNO(KOUNT)=IELEMNO_TMP(nodeID,iElem)
            KOUNT=KOUNT+1
         ENDDO
         NINDX(nodeID)=KOUNT
      ENDDO
      
      DEALLOCATE(IELEMNO_TMP,NINDX_TMP)

      
      END SUBROUTINE
      
      ! Calculates the equivalent finite strain (scalar) from a deformation gradient tensor:
      ! E_GL = 1/2*[F^T*F - I]
      ! E_equiv = Sqrt(2/3*||E_GL||)
      SUBROUTINE calcE_Equiv(F,E_Equiv)
      implicit none
      real(8), intent(in) :: F(9)
      real(8), intent(out):: E_Equiv
      
      real(8) :: F_2D(3,3), E_GreenLagrange(3,3)
      
      F_2D(1,1:3) = F(1:3)
      F_2D(2,1:3) = F(4:6)
      F_2D(3,1:3) = F(7:9)
      
      E_GreenLagrange = MATMUL(TRANSPOSE(F_2D),F_2D)
      E_GreenLagrange(1,1) = E_GreenLagrange(1,1) - 1.D0
      E_GreenLagrange(2,2) = E_GreenLagrange(2,2) - 1.D0
      E_GreenLagrange(3,3) = E_GreenLagrange(3,3) - 1.D0
      E_GreenLagrange = 0.5D0 * E_GreenLagrange
      
      E_Equiv = DSqrt(2.0D0/3.0D0*(E_GreenLagrange(1,1)**2 + &
                                   E_GreenLagrange(1,2)**2 + &
                                   E_GreenLagrange(1,3)**2 + &
                                   E_GreenLagrange(2,1)**2 + &
                                   E_GreenLagrange(2,2)**2 + &
                                   E_GreenLagrange(2,3)**2 + &
                                   E_GreenLagrange(3,1)**2 + &
                                   E_GreenLagrange(3,2)**2 + &
                                   E_GreenLagrange(3,3)**2 ))
      END SUBROUTINE
      
      SUBROUTINE calcGLStrain(F,glStrain)
      implicit none
      real(8), intent(in) :: F(9)
      real(8), intent(out):: glStrain(6)
      ! locals
      real(8) :: F_2D(3,3), F_2D_TR(3,3), fp_t_fp(3,3), glStrain_2D(3,3)
      real(8) :: ckrone_2d(3,3)
      integer :: i,j
      
      ckrone_2d = 0.d0
      ckrone_2d(1,1) = 1.d0
      ckrone_2d(2,2) = 1.d0
      ckrone_2d(3,3) = 1.d0
      
      F_2D(1,1:3) = F(1:3)
      F_2D(2,1:3) = F(4:6)
      F_2D(3,1:3) = F(7:9)
!   Green-Lagrange strain for plastic deformation
      F_2D_TR = TRANSPOSE(F_2D)
      fp_t_fp = MATMUL(F_2D_TR,F_2D)
      glStrain_2D = 0.d0
      do i=1,3
         do j=1,3
            glStrain_2D(i,j)=0.5d0*(fp_t_fp(i,j)-ckrone_2d(i,j))
         enddo
      enddo
      glStrain(1)=glStrain_2D(1,1)
      glStrain(2)=glStrain_2D(2,2)
      glStrain(3)=glStrain_2D(3,3)
      glStrain(4)=glStrain_2D(1,2)
      glStrain(5)=glStrain_2D(1,3)
      glStrain(6)=glStrain_2D(2,3)
      END SUBROUTINE