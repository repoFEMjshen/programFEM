      module MeshInfo
      
         implicit none
         integer:: nNodes,nElements,nNodesPerElement,nDimensions,nGrainsAttachedToNodesMax,&
                   nNodesVirtual,MAXADJGRAINELEM
         real(8):: elementSize
         integer, allocatable:: connectivity(:,:) , mappingElementToGrain(:), &
                   connectivityVirtual(:,:), INODGR(:,:), NNODGR(:), NATEL(:), &
                   IATEL(:,:),NODE_REPOSITORY(:,:)
         real(8), allocatable:: coordinateNodes(:,:) , coordinateNodesVirtual(:,:) , &
                   fieldValueAtGauss(:,:), WATEL(:,:)
         

      end module MeshInfo