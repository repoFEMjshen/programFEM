      PROGRAM shuffleCrystallographicAssignments
      
      use grains, only : nGrains, &
         readTexture_direct,grainTexture,textureImported
            
      implicit none

      real(8) :: rndNum
      integer :: nGrainsInFile
      real(8), allocatable :: grainTexture_Shuffled(:,:)
      real(8), allocatable :: listRandnGrains(:)
      integer, allocatable :: newOrderGrains(:)
      real(8), parameter :: PI = 3.14159265359d0
      
      logical :: success
            
      integer :: i,j
      
      ! initialize seed
      CALL init_random_seed()
      
      !read texture: grainTexture(:)
      CALL readTexture_direct('grainTexture.inp',nGrainsInFile)
      if(textureImported) write(*,*) 'texture imported for # grains:',nGrainsInFile
      if(.NOT.textureImported) then
         write(*,*) 'WARNING ! no texture information (grainTexture.inp) found. assuming single crystal, euler:',grainTexture(:,1)
      endif
      
      !allocate the orientation list
      allocate(grainTexture_Shuffled(3,nGrains))
      allocate(listRandnGrains(nGrains))
      allocate(newOrderGrains(nGrains))

      do i = 1,nGrains
            
         CALL RANDOM_NUMBER(rndNum)
         listRandnGrains(i) = rndNum
         
      enddo
      
      newOrderGrains = (/(i,i=1,nGrains)/)
      call SortBottom(listRandnGrains, newOrderGrains, nGrains, nGrains)
      
      grainTexture_Shuffled(:,:) = grainTexture(:,newOrderGrains(:))
      
      call saveGrainTexture('grainTexture_shuffled.inp',grainTexture_Shuffled,nGrains)
      
      END PROGRAM
