      PROGRAM printGrainNeighbors
      
      use meshtools
      use options
      use grains
      
      implicit none
      
      include 'PARDIS.H'
      
      integer :: NX, NELX
      integer, parameter :: NODE = 4
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX)
      integer :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      character(len=30) :: meshName
      
      integer :: iElem,iGrain,jGrain
      real(8), allocatable :: grainRot(:,:,:)
      integer, allocatable :: elemFaceNeighbors(:,:)
      real(8), allocatable :: elemCenters(:,:)
      real(8), allocatable :: elemVol(:)
      real(8), allocatable :: f(:)
      real(8) :: Aij, deltaV, dGrain1, dGrain2
      real(8) :: volDomain, max_r, r, rGrains
      integer :: ri
      integer, parameter :: res_r = 30
      real(8) :: grainMisAngle, thisMisInd, thisMisAngle
      real(8) :: cAxis(3)
         
      character :: eulerDirection
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: distVec(3), dist, weight, normalizationConstant
      
      logical :: success, successReadNeighbors, successReadNodeConnect
      
      ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)


      integer :: staIndex, endIndex, elemIndex
      integer :: iElemID, jElemID, nElem
      integer :: grainIdx, iFace, iGrainID, jGrainID
      
      integer :: i,j
      
      ! read in the mesh and connectivity
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      
      write(*,'(A)',advance='no') ' identifying node-to-element connectivity... '
      CALL readNODE_ELEM_CONNEC('nodeElementConnect.dat',IELEMNO,NINDX,NX,successReadNodeConnect)
      if (.not.successReadNodeConnect) then
         write(*,'(A)',advance='no') 'calculating... '
         CALL NODE_ELEM_CONNEC(IJK,IELEMNO,NINDX,NX,NELX,NODE)
         write(*,*) 'done.'
         CALL writeNODE_ELEM_CONNEC('nodeElementConnect.dat', &
                                    IELEMNO,NINDX,NELX,NX,MAXNODELEM,MAXNODE)
         write(*,*) 'saved connectivity data file: nodeElementConnect.dat'
      else
         write(*,*) 'imported.'
      endif
      
      ! ---------------------- IMPORT GRAINS ------------------------!
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      CALL readGrains('grains.inp',NX,NELX)
      if(grainsImported) &
         write(*,'(A,I0,A)') ' ', nGrains,' grains imported from grains.inp'   

      if(.NOT.grainsImported) then
         CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
         if(grainsImported) &
            write(*,'(A,I0,A)') ' ', nGrains,' grains imported from elementGrains.inp'
      endif
      
      if(.NOT.grainsImported) then
         CALL createSingleGrain(NX,NELX)
         write(*,'(A,I0,A)') 'single grain assumed'
      endif

      ! CALCULATE ELEMENT-ELEMENT CONNECTIVITY
      allocate(elemFaceNeighbors(NELX,4))
      write(*,'(A)',advance='no') ' identifying element-element connectivity... '
      CALL readFaceNeighborsTET('elementNeighbors.dat',elemFaceNeighbors,NELX,successReadNeighbors)
      if (.not.successReadNeighbors) then
         CALL calcFaceNeighborsTET(elemFaceNeighbors(:,:),IJK,IELEMNO,NINDX, & 
                                   NELX,NODE,MAXNODELEM,MAXNODE)
         write(*,*) 'done'
      else
         write(*,*) 'imported'
      endif
      if(.not.successReadNeighbors) then
         CALL writeFaceNeighborsTET('elementNeighbors.dat',elemFaceNeighbors,NELX)
         write(*,*) 'saved connectivity data file: elementNeighbors.dat'
      endif

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      CALL identifyGrainNodes(NX,NELX,IJK)

      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)

      ! domainRange0: initial domain range
      ! domainRange: may evolve in time
      domainLen(:) = dabs(domainRange(:,2) - domainRange(:,1))
      
      !read grain sizes: grainSize(:)
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported) write(*,*) 'grain sizes imported'
            
      !read grain phases: grainPhase(:)
      CALL readGrainPhases('grainPhases.inp')
      if(grainPhasesImported) write(*,*) 'grain phases imported'
      
      !read texture: grainTexture(:)
      CALL readTexture('grainTexture.inp')
      if(textureImported) write(*,*) 'texture imported'
      

      eulerDirection=' '
      do while (eulerDirection.NE.'A'.AND.eulerDirection.NE.'P')
         write(*,*) 'Grain orientations are provided in form of Euler angles'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         write(*,*) '(A) Euler angles rotate the specimen coordinate axes onto the crystal coordinate axes '&
                   &'(OIM standard, Kourosh, Deniz)'
         write(*,*) '(P) Euler angles rotate the crystal coordinate axes onto the specimen coordinate axes '&
                   &'(was used by Deka)'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         read(*,'(A1)') eulerDirection
         if (eulerDirection.EQ.'a') eulerDirection='A'
         if (eulerDirection.EQ.'p') eulerDirection='P'
      enddo
      write(*,*)
      
      allocate(grainRot(3,3,nGrains))
      do iGrainID = 1, nGrains
         call euler_slip(grainTexture(1,iGrainID), &
                         grainTexture(2,iGrainID), &
                         grainTexture(3,iGrainID), &
                         grainRot(:,:,iGrainID))
         if (eulerDirection=='P') then
            grainRot(:,:,iGrainID) = TRANSPOSE(grainRot(:,:,iGrainID))
         endif
      enddo      
      
      write(*,*) 'calculating grain centroids...'
      allocate(elemVol(NELX))
      allocate(elemCenters(3,NELX))
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL getElemXCenter(iElem,IJK,G0XYZ,NODE,elemCenters(:,iElem))
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      CALL calcGrainCentriods(elemCenters,elemVol,NELX)
      
      volDomain = 0.d0
      do iGrain=1,nGrains
         volDomain = volDomain + grainVolumes(iGrain)
      enddo
      max_r = (domainLen(1)+domainLen(2)+domainLen(3)) / 3
      allocate(f(res_r))
      f = 0.d0
      do iGrain=1,nGrains
         do jGrain=1,nGrains
            dGrain1 = grainSize(iGrain)
            dGrain2 = grainSize(jGrain)
            Aij = dabs(DOT_PRODUCT(grainRot(:,3,iGrain),grainRot(:,3,jGrain))) / volDomain * (dGrain1**2+dGrain2**2)
            do ri = 1,res_r
               r = max_r / (res_r-1) * (ri-1)
               CALL calcIntersectingVolume(rGrains,dGrain1,dGrain2,r,deltaV)
               f(ri) = f(ri) + Aij * deltaV
            enddo
         enddo
      enddo
      
      
      deallocate(f)
      deallocate(elemFaceNeighbors)
      
      
      END PROGRAM

      

!*********************************************************
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
      integer :: error
      character(len=50) :: fileName
      character(len=150) :: lineStr
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ
      
      open(201,file='loadtime.inp')
      call READSTR(201,lineStr,error)
      read(201,*) icyc
      call READSTR(201,lineStr,error)
      read(201,*) meshName
      close(201)
      
      len=len_trim(meshName)
      fileName=meshName(1:len)//'.inp'

      open(unit=LR,file=fileName)

      call READSTR(LR,lineStr,error)
      read(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD     
      
      call READSTR(LR,lineStr,error)
      read(LR,*) NX
      NEQ=NDF*NX
      
      IF(NDF.GT.NDOFELX)THEN
         WRITE(*,*)'INSUFFICIENT MEM-NDOFELX'
         STOP
      ENDIF

      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         write(*,*) 'Increase MAXNODE to ', NX
         WRITE(*,*)'Increase MAXNODE to ', NX
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
      
      
      call READSTR(LR,lineStr,error)
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
      
      SUBROUTINE READSTR(fileNum,lineStr,error)
      IMPLICIT REAL*8(A-H,O-Z)
      ! this code skips 1 line of text in the input file.
      integer, intent(in) :: fileNum
      character(150), intent(out) :: lineStr
      integer, intent(out) :: error
      read(fileNum,'(A)',iostat=error) lineStr
      RETURN
      END
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
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
      SUBROUTINE NODE_ELEM_CONNEC(IJK,IELEMNO,NINDX,NX,NELX,NODE)
      implicit none
      INCLUDE 'PARDIS.H'
      integer, intent(in) :: NX,NELX,NODE
      integer, intent(in):: IJK(MNELX)
      integer, intent(out)::IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************
      integer, allocatable :: IELEMNO_TMP(:,:),NINDX_TMP(:)
      integer :: nodeID, I, J, KOUNT

! ----------- slow way of doing this
!      ITEMP=1
!      DO I=1,NX
!         DO J=1,NELX
!            DO K=1,NODE
!               IF(IJK((J-1)*NODE+K).EQ.I)THEN
!                  IELEMNO(ITEMP)=J
!                  ITEMP=ITEMP+1
!                  GOTO 20
!               ENDIF
!            ENDDO
!   20    ENDDO
!         NINDX(I)=ITEMP
!      ENDDO
! ----------- faster way of doing this:
      ALLOCATE(IELEMNO_TMP(NX,MAXNODELEM),NINDX_TMP(NX))
      NINDX_TMP=0
      DO I=1,NELX
         DO J=1,NODE
            nodeID = IJK((I-1)*NODE+J)
            NINDX_TMP(nodeID)=NINDX_TMP(nodeID)+1
            IELEMNO_TMP(nodeID,NINDX_TMP(nodeID))=I
         ENDDO
      ENDDO
      
      KOUNT=1
      DO I=1,NX
         DO J=1,NINDX_TMP(I)
            IELEMNO(KOUNT)=IELEMNO_TMP(I,J)
            KOUNT=KOUNT+1
         ENDDO
         NINDX(I)=KOUNT
      ENDDO
      
      DEALLOCATE(IELEMNO_TMP,NINDX_TMP)
!------------------------------------------
      
      
      END SUBROUTINE
      
      
!***********************************************************************
!C     WRITES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE writeNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX, &
                                       NELX,NX,MAXNODELEM,MAXNODE)
                                       
      implicit none
      character(len=*), intent(in) :: strConnectFile      

      integer, intent(in) :: NX,NELX,MAXNODELEM,MAXNODE
      integer, intent(in) :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, I, lastIndex
      integer :: elemStart, elemEnd

      
      open(716,file=trim(strConnectFile))
      write(716,*,iostat=error) NX
      lastIndex = 1
      do iNode=1,NX
         elemStart=1
         if(iNode.ne.1) elemStart=NINDX(iNode-1)
         elemEnd = NINDX(iNode)-1
         nElem = elemEnd - elemStart + 1
         
         write(716,*,iostat=error) nElem
         write(716,*,iostat=error) (IELEMNO(iElem),iElem=elemStart,elemEnd)
         lastIndex = lastIndex + nElem
      enddo
      
      close(716)
      
      END SUBROUTINE
      
      subroutine euler_slip(phi,theta,omega,tlgt)
      implicit double precision (a-h,o-z)  

      dimension tlg(3,3),tlgt(3,3)
      real(8) :: phi,theta,omega

      pi=4.d0*datan(1.d0)

  
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

      tlgt = TRANSPOSE(tlg)
      
      end subroutine      

      SUBROUTINE calcIntersectingVolume(rGrains,dGrain1,dGrain2,r,deltaV)
      implicit none
      real(8), intent(in) :: rGrains,dGrain1,dGrain2,r
      real(8), intent(out):: deltaV
      
      deltaV = r-rGrains+(dGrain1+dGrain2)/2.d0*(dGrain2**2)
      END SUBROUTINE