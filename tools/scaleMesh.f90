
      program scaleMesh
      
      implicit none
      
      integer, parameter :: MAXNODE=60000
      integer, parameter :: mesh_in=10, mesh_out=11
      integer, parameter :: grains_in=20, grains_out=21
      real(8), parameter :: PI=3.14159265359d0

      real(8) :: G0XYZ(3,MAXNODE)
      integer :: NX, NELX

      real(8) :: xc(3), domainLenZ, domainLenX,notchHeight, notchDepth, stretchX_new
      real(8) :: X,Y,Z,X_new,Y_new,Z_new,Z_notch
      integer :: MDIM,NDFADD,NDF,NNODE,NGAUSS
      integer :: error
      integer :: iPos
      integer :: iGrain,nGrains,grainID
      real(8) :: grainSize,meanGrainSize,targetMeanGrainSize
      logical :: scaleGrainSizesToo
      character(len=64) :: meshName
      character(len=50) :: fileName, fileName_out
      integer, parameter :: linesize = 16000
      character(len=linesize) :: lineStr
      real(8) :: factorScaling
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ

      open(201,file='loadtime.inp')
      call READSTR(201,lineStr,linesize,error)
      read(201,*) icyc
      call READSTR(201,lineStr,linesize,error)
      read(201,*) meshName
      close(201)

      len=len_trim(meshName)
      fileName=meshName(1:len)//'.inp'
      fileName_out=meshName(1:len)//'_scaled.inp'

      open(unit=mesh_in,file=fileName)
      open(unit=mesh_out,file=fileName_out)

      call READSTR(mesh_in,lineStr,linesize,error)
      write(mesh_out,'(A)') trim(lineStr)
      
      read(mesh_in,*) MDIM,NDFADD
      write(mesh_out,*) MDIM,NDFADD
      NDF=MDIM+NDFADD     

      call READSTR(mesh_in,lineStr,linesize,error)
      write(mesh_out,'(A)') trim(lineStr)
      read(mesh_in,*) NX
      NEQ=NDF*NX

      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         write(*,*) 'Increase MAXNODE to ', NX
         WRITE(*,*)'Increase MAXNODE to ', NX
         STOP
      ENDIF
      NEQ=NDF*NX

      ! READ INITIAL POSITIONS
      do iNode=1,NX
         read(mesh_in,*) nnum,(G0XYZ(ii,iNode),ii=1,MDIM)
      enddo
      
      write(*,*) 'enter the factor of scaling to be applied on the nodal coordinates,'
      write(*,*) 'or enter 0 to specify target mean grain size'
      read(*,*) factorScaling
      
      if (factorScaling == 0.d0) then  ! get scaling factor from target / current mean grain size
      
         write(*,*) 'enter target mean grain size'
         read(*,*) targetMeanGrainSize
      
         open(unit=grains_in,file='grainSizes.inp')
         read(grains_in,*) nGrains
         meanGrainSize = 0.d0
         do iGrain=1,nGrains
            read(grains_in,*) grainID,grainSize
            meanGrainSize = meanGrainSize + grainSize
         enddo
         meanGrainSize = meanGrainSize / nGrains
         close(grains_in)
         
         factorScaling = targetMeanGrainSize / meanGrainSize
         scaleGrainSizesToo = .true.
         
      else
         
         write(*,*) 'scale grainSizes.inp too? (T/F)'
         read(*,*) scaleGrainSizesToo
         
      endif

      ! scale nodes
      G0XYZ = factorScaling*G0XYZ
      
      ! write new mesh
      write(mesh_out,*) NX
      do iNode=1,NX
         write(mesh_out,*) iNode,G0XYZ(:,iNode)
      enddo
      
      error = 0
      do while(error == 0)
         call READSTR(mesh_in,lineStr,linesize,error)
         if (error /= 0 ) exit
         write(mesh_out,'(A)') trim(lineStr)
         do iPos=lineSize-20,linesize
            if (lineStr(iPos:iPos)/=' ') then
               write(*,*) 'line truncated! increase linesize beyond',linesize
               write(*,*) trim(lineStr(lineSize-20:lineSize))
               stop
            endif
         enddo
      enddo

      close(mesh_in)
      close(mesh_out)
      
      if (scaleGrainSizesToo) then
         open(unit=grains_in,file='grainSizes.inp')
         open(unit=grains_out,file='grainSizes_scaled.inp')
         read(grains_in,*) nGrains
         write(grains_out,*) nGrains
         do iGrain=1,nGrains
            read(grains_in,*) grainID,grainSize
            write(grains_out,*) grainID,grainSize*factorScaling
         enddo
      endif
      
      write(*,*) 'done'
      end program
      
      SUBROUTINE READSTR(fileNum,lineStr,linesize,error)
      IMPLICIT REAL*8(A-H,O-Z)
      ! this code skips 1 line of text in the input file.
      integer, intent(in) :: fileNum,linesize
      character(len=linesize), intent(out) :: lineStr
      integer, intent(out) :: error
      read(fileNum,'(A)',iostat=error) lineStr
      RETURN
      END
