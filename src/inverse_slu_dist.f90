subroutine inverse_slu_dist(n,iproc,nprocs,N_update,in_st,irowptr,ibindx,val,resid,isize)

! Seven basic steps are required:
!   1. Create C structures used in SuperLU_DIST
!   2. Initialize the MPI environment and the SuperLU process grid
!   3. Set up the input matrix and the right-hand side
!   4. Set the options argument
!   5. Call f_pdgssvx
!   6. Release the process grid and terminate the MPI environment
!   7. Release all structures

      use superlu_mod
      implicit real*8(a-h,o-z)

!#include "mafdecls.fh"
!#include "global.fh"
      include 'MPIF.H'
      include 'PARDIS.H'

      dimension ibindx(isize),val(isize),resid(N_update),irowptr(N_update+1)
      integer,dimension(:),allocatable::ibindx_loc,irowptr_loc

! From superlupara.f90 (module): integer, parameter :: superlu_ptr = 8
! (64-bit)
      integer(superlu_ptr) :: grid
      integer(superlu_ptr) :: options
      integer(superlu_ptr) :: ScalePermstruct
      integer(superlu_ptr) :: LUstruct
      integer(superlu_ptr) :: SOLVEstruct
      integer(superlu_ptr) :: A
      integer(superlu_ptr) :: stat
      


      call f_create_gridinfo_handle(grid)
      call f_create_options_handle(options)
      call f_create_ScalePerm_handle(ScalePermstruct)
      call f_create_LUstruct_handle(LUstruct)
      call f_create_SOLVEstruct_handle(SOLVEstruct)
      call f_create_SuperMatrix_handle(A)
      call f_create_SuperLUStat_handle(stat)


!!$      nprow = 1
!!$      npcol = nprocs

      open(102,file='slu_proc_distrib.inp')
      read(102,*)nprow,npcol
      close(102)
!      nprow = 8
!      npcol = 4

      if(nprow*npcol.ne.nprocs)then
         write(*,*)'nprow*npcol not equal to nprocs'
         stop
      endif
         
         

      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

      call f_superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid)
      call get_GridInfo(grid, iam=iam)

      if ( iam >= nprow * npcol ) then
         go to 100
      endif


!------Write---------------------------------
!!$      write(901+iproc,*)irowptr(N_update+1)
!!$      do i=1,N_update
!!$         is=irowptr(i)
!!$         ie=irowptr(i+1)-1
!!$         if(i.eq.N_update)ie=irowptr(i+1)
!!$         do j=is,ie
!!$            write(901+iproc,*)i+in_st,ibindx(j),val(j)
!!$         enddo
!!$         write(801+iproc,*)i+in_st,resid(i)
!!$      enddo
!---------------------------------------------

      allocate(ibindx_loc(irowptr(N_update+1)),irowptr_loc(N_update+1))
      
      m=n

      nnz_loc=0
      do i=1,N_update
         is=irowptr(i)
         ie=irowptr(i+1)-1
         if(i.eq.N_update)ie=irowptr(i+1)
         do j=is,ie
            ibindx_loc(j)=ibindx(j)-1
            nnz_loc=nnz_loc+1
         enddo
      enddo


      do i=1,N_update
         irowptr_loc(i)=irowptr(i)-1
      enddo
      irowptr_loc(N_update+1)=irowptr(N_update+1)

      call mpi_barrier(mpi_comm_world,ierror)

!      stop


      call f_dCreate_CompRowLoc_Mat_dist(A, m, n, nnz_loc, N_update, in_st, &
           val, ibindx_loc, irowptr_loc, SLU_NR_loc, SLU_D, SLU_GE)

      nrhs = 1


      call get_CompRowLoc_Matrix(A, nrow_loc=N_update)

      call f_set_default_options(options)
      call set_superlu_options(options,ColPerm=NATURAL)
      call set_superlu_options(options,RowPerm=NOROWPERM)
      call get_SuperMatrix(A,nrow=m,ncol=n)
      call f_ScalePermstructInit(m, n, ScalePermstruct)
      call f_LUstructInit(m, n, LUstruct)
      call f_PStatInit(stat)

      call f_pdgssvx(options, A, ScalePermstruct, resid, N_update, nrhs, &
                     grid, LUstruct, SOLVEstruct, berr, stat, info)


      deallocate(ibindx_loc,irowptr_loc)


! Deallocate the storage allocated by SuperLU_DIST
      call f_PStatFree(stat)
      call f_Destroy_SuperMat_Store_dist(A)
      call f_ScalePermstructFree(ScalePermstruct)
      call f_Destroy_LU(n, grid, LUstruct)
      call f_LUstructFree(LUstruct)
      call get_superlu_options(options, SolveInitialized=init)
      if (init == YES) then
         call f_dSolveFinalize(options, SOLVEstruct)
      endif

!      write(*,*)'dbg20',iam

100   call f_superlu_gridexit(grid)

! Deallocate the C structures pointed to by the Fortran handles
      call f_destroy_gridinfo_handle(grid)
      call f_destroy_options_handle(options)
      call f_destroy_ScalePerm_handle(ScalePermstruct)
      call f_destroy_LUstruct_handle(LUstruct)
      call f_destroy_SOLVEstruct_handle(SOLVEstruct)
      call f_destroy_SuperMatrix_handle(A)
      call f_destroy_SuperLUStat_handle(stat)

      write(*,*)'exit inv'
   
      return
      end

