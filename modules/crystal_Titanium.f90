
      module crystal_Titanium
      
      implicit none 
      
      integer, parameter :: nslipfam=8, nslip_max=100    
      integer, parameter :: HCP=1, BCC=2
      integer, parameter :: pri_A=1, trans_B=2
      
      real(8) :: cst_HCP(3,nslip_max)  ! slip directions for HCP
      real(8) :: cmt_HCP(3,nslip_max)  ! slip planes for HCP
      real(8) :: cst_BCC(3,nslip_max)  ! slip directions for HCP
      real(8) :: cmt_BCC(3,nslip_max)  ! slip planes for HCP
      
      real(8) :: sin_nm_HCP(nslip_max,nslip_max)
      real(8) :: sin_nt_HCP(nslip_max,nslip_max)
      real(8) :: sin_nn_HCP(nslip_max,nslip_max)
      real(8) :: cos_nm_HCP(nslip_max,nslip_max)
      real(8) :: cos_nt_HCP(nslip_max,nslip_max)
      real(8) :: cos_nn_HCP(nslip_max,nslip_max)


      integer :: nslips_famHCP(nslipfam)
      integer :: nslips_famBCC(nslipfam)
      integer :: nslip_HCP, nslip_BCC

      contains
      !-------------------------------------------------------------------------
      !    THIS SUBROUTINE INITIALISES THE CRYSTAL SLIP SYSTEMS, 
      !    CALCULATES THE TRANSFORMATION MATRIX AND CONTAINS 
      !    STATE VARIABLES
      !----------------------------------------------------------
      !     cmt - slip plane in iCrystal system
      !     cst - slip directions in iCrystal system
      !     tlg0 - transformation matrix (kalidindis paper)
      !     nslip - number of slip systems: HCP:30, bcc:48
      !     nslips_fam(ifam) - number of slip systems in each slip family. 
      !                       ifam: slip familty
      !------------------------------------------------------------------------
      !     Modified by Deniz: 
      !     (the signs of) the burgers vectors are chosen such that
      !     they are towards the soft directions 
      !     (see Roundy Krenn Cohen Morris 2001, for a discussion of 
      !     soft/hard directions along (221)<111> and (123)<111>)
      !-------------------------------------------------------------------
      SUBROUTINE init_crys_HCP()
      
      implicit none
      
      integer :: nslips_pra_HCP(nslipfam)
      
      real(8), parameter :: r1_2=1.d0/2.d0, &
                            rsq3=dsqrt(3.d0)
      integer :: ic,i
      real(8) :: bunb,bush,bunbo,bunshi,bunsh

      ! HCP -- 30 slip systems
      nslip_HCP = 30      

      !{0001}<11-20>(3) basal, a
      !{10-10}<11-20>(3) prismatic, a
      !{10-1-1}<11-20>(6)1st pyramidal, a
      !{10-11}<11-23>(12)1st pyramidal, c+a
      !{11-22}<11-23>(6) 2nd pyramidal c+a

      ! generate the slip system families of an HCP lattice within transformed-Beta.
      ! HCP symmetries are broken in trans-Beta grains, due to neighboring BCC laths:
      ! all <a>basals and <a>prismatics are now distinguishable
      nslips_famHCP(1)=1 ! 1. basal<a>
      nslips_famHCP(2)=1 ! 2. basal<a>
      nslips_famHCP(3)=1 ! 3. basal<a>

      nslips_famHCP(4)=1 ! prismatic<a>
      nslips_famHCP(5)=1 ! prismatic<a>
      nslips_famHCP(6)=1 ! prismatic<a>

      nslips_famHCP(7)=6 ! 6 1st order pyramidal<a>'s
      nslips_famHCP(8)=18! # of <c+a>'s         

      ! generate the slip system families for HCP
      nslips_pra_HCP(1)=3 ! basal<a>
      nslips_pra_HCP(2)=3 ! prismatic<a>
      nslips_pra_HCP(3)=6 ! 1st pyramidal<a>
      nslips_pra_HCP(4)=12! 1st pyramidal<c+a>
      nslips_pra_HCP(5)=6 ! 2nd pyramidal<c+a>
      
      !      e2|  a2
      !        | /
      !        |/
      ! a3-----*------e1
      !         \
      !          \
      !           a1
      ic=0
      !---basal slip system
      ! slip planes
      do i=1,3
      cmt_HCP(1,i)=  0.0
      cmt_HCP(2,i)=  0.0 
      cmt_HCP(3,i)=  1.0d0
      enddo

      ! slip directions
      cst_HCP(1,1)=  r1_2
      cst_HCP(2,1)= -r1_2*rsq3 
      cst_HCP(3,1)=  0.0

      cst_HCP(1,2)=  r1_2
      cst_HCP(2,2)=  rsq3*r1_2
      cst_HCP(3,2)=  0.0

      cst_HCP(1,3)= -1.0d0
      cst_HCP(2,3)=  0.0
      cst_HCP(3,3)=  0.0
      ic=ic+nslips_pra_HCP(1)

      !---    prismatic slip system
      cmt_HCP(1,1+ic)=  0.0
      cmt_HCP(2,1+ic)=  1.0d0 
      cmt_HCP(3,1+ic)=  0.0

      cmt_HCP(1,2+ic)= -rsq3*r1_2
      cmt_HCP(2,2+ic)=  r1_2
      cmt_HCP(3,2+ic)=  0.0
      
      cmt_HCP(1,3+ic)= -rsq3*r1_2
      cmt_HCP(2,3+ic)= -r1_2 
      cmt_HCP(3,3+ic)=  0.0

      cst_HCP(1,1+ic)=  -1.0d0
      cst_HCP(2,1+ic)=  0.0
      cst_HCP(3,1+ic)=  0.0

      cst_HCP(1,2+ic)=  r1_2
      cst_HCP(2,2+ic)=  rsq3*r1_2
      cst_HCP(3,2+ic)=  0.0

      cst_HCP(1,3+ic)=  r1_2
      cst_HCP(2,3+ic)= -rsq3*r1_2
      cst_HCP(3,3+ic)=  0.0
      
      ic=ic+nslips_pra_HCP(2)

      !       pyramidal a slip
      bunbo=sqrt(4*(1.587)**2+3)
      bunshi=1.587
      cmt_HCP(1,1+ic)=0.0
      cmt_HCP(2,1+ic)=-2*bunshi/bunbo
      cmt_HCP(3,1+ic)=rsq3/bunbo

      cmt_HCP(1,2+ic)=rsq3*bunshi/bunbo
      cmt_HCP(2,2+ic)=-bunshi/bunbo
      cmt_HCP(3,2+ic)=rsq3/bunbo

      cmt_HCP(1,3+ic)= rsq3*bunshi/bunbo
      cmt_HCP(2,3+ic)= bunshi/bunbo
      cmt_HCP(3,3+ic)= rsq3/bunbo 

      cmt_HCP(1,4+ic)= 0.0
      cmt_HCP(2,4+ic)= 2*bunshi/bunbo
      cmt_HCP(3,4+ic)= rsq3/bunbo

      cmt_HCP(1,5+ic)=-rsq3*bunshi/bunbo
      cmt_HCP(2,5+ic)= bunshi/bunbo
      cmt_HCP(3,5+ic)= rsq3/bunbo

      cmt_HCP(1,6+ic)=-rsq3*bunshi/bunbo
      cmt_HCP(2,6+ic)=-bunshi/bunbo
      cmt_HCP(3,6+ic)= rsq3/bunbo


      cst_HCP(1,1+ic)= 1.0d0
      cst_HCP(2,1+ic)= 0.0
      cst_HCP(3,1+ic)= 0.0 

      cst_HCP(1,2+ic)= r1_2
      cst_HCP(2,2+ic)= r1_2*rsq3
      cst_HCP(3,2+ic)= 0.0        

      cst_HCP(1,3+ic)=-r1_2
      cst_HCP(2,3+ic)= r1_2*rsq3
      cst_HCP(3,3+ic)= 0.0 

      cst_HCP(1,4+ic)=-1.0d0
      cst_HCP(2,4+ic)= 0.0
      cst_HCP(3,4+ic)= 0.0

      cst_HCP(1,5+ic)=-r1_2
      cst_HCP(2,5+ic)=-r1_2*rsq3
      cst_HCP(3,5+ic)= 0.0  

      cst_HCP(1,6+ic)= r1_2
      cst_HCP(2,6+ic)=-r1_2*rsq3
      cst_HCP(3,6+ic)= 0.0
      ic=ic+nslips_pra_HCP(3)

      !       1st order <c+a> slip
      bunbo=sqrt(4.0*(1.587)**2+3.0)
      bunshi=1.587

      cmt_HCP(1,1+ic)=  0.0
      cmt_HCP(2,1+ic)=  -2*bunshi/bunbo 
      cmt_HCP(3,1+ic)=  rsq3/bunbo

      cmt_HCP(1,2+ic)=  rsq3*bunshi/bunbo
      cmt_HCP(2,2+ic)=  -bunshi/bunbo 
      cmt_HCP(3,2+ic)=  rsq3/bunbo

      cmt_HCP(1,3+ic)=  rsq3*bunshi/bunbo
      cmt_HCP(2,3+ic)=  bunshi/bunbo 
      cmt_HCP(3,3+ic)=  rsq3/bunbo

      cmt_HCP(1,4+ic)=  0.0
      cmt_HCP(2,4+ic)=  2*bunshi/bunbo 
      cmt_HCP(3,4+ic)=  rsq3/bunbo

      cmt_HCP(1,5+ic)=  -rsq3*bunshi/bunbo
      cmt_HCP(2,5+ic)=  bunshi/bunbo 
      cmt_HCP(3,5+ic)=  rsq3/bunbo

      cmt_HCP(1,6+ic)=  -rsq3*bunshi/bunbo
      cmt_HCP(2,6+ic)=  -bunshi/bunbo 
      cmt_HCP(3,6+ic)=  rsq3/bunbo

      cmt_HCP(1,7+ic)=  0.0
      cmt_HCP(2,7+ic)=  -2*bunshi/bunbo 
      cmt_HCP(3,7+ic)=  rsq3/bunbo

      cmt_HCP(1,8+ic)=  rsq3*bunshi/bunbo
      cmt_HCP(2,8+ic)=  -bunshi/bunbo 
      cmt_HCP(3,8+ic)=  rsq3/bunbo

      cmt_HCP(1,9+ic)=  rsq3*bunshi/bunbo
      cmt_HCP(2,9+ic)=  bunshi/bunbo 
      cmt_HCP(3,9+ic)=  rsq3/bunbo

      cmt_HCP(1,10+ic)=  0.0
      cmt_HCP(2,10+ic)=  2*bunshi/bunbo 
      cmt_HCP(3,10+ic)=  rsq3/bunbo

      cmt_HCP(1,11+ic)= - rsq3*bunshi/bunbo
      cmt_HCP(2,11+ic)=  bunshi/bunbo 
      cmt_HCP(3,11+ic)=  rsq3/bunbo

      cmt_HCP(1,12+ic)=  -rsq3*bunshi/bunbo
      cmt_HCP(2,12+ic)=  -bunshi/bunbo 
      cmt_HCP(3,12+ic)=  rsq3/bunbo

      !-----------slip directions

      bunb=2.0*(sqrt((1.587)**2+1.0))
      bush=1.587
      cst_HCP(1,1+ic)=  1.0/bunb
      cst_HCP(2,1+ic)=  rsq3/bunb 
      cst_HCP(3,1+ic)=  2*bush/bunb

      cst_HCP(1,2+ic)=  -1.0/bunb
      cst_HCP(2,2+ic)=  rsq3/bunb 
      cst_HCP(3,2+ic)=  2*bush/bunb

      cst_HCP(1,3+ic)=  -2.0/bunb
      cst_HCP(2,3+ic)=  0.0 
      cst_HCP(3,3+ic)=  2*bush/bunb

      cst_HCP(1,4+ic)=  -1.0/bunb
      cst_HCP(2,4+ic)=  -rsq3/bunb 
      cst_HCP(3,4+ic)=  2*bush/bunb

      cst_HCP(1,5+ic)=  1.0/bunb
      cst_HCP(2,5+ic)=  -rsq3/bunb 
      cst_HCP(3,5+ic)=  2*bush/bunb

      cst_HCP(1,6+ic)=  2.0/bunb
      cst_HCP(2,6+ic)=  0.0 
      cst_HCP(3,6+ic)=  2*bush/bunb

      cst_HCP(1,7+ic)=  -1.0/bunb
      cst_HCP(2,7+ic)=  rsq3/bunb 
      cst_HCP(3,7+ic)=  2*bush/bunb

      cst_HCP(1,8+ic)=  -2.0/bunb
      cst_HCP(2,8+ic)=  0.0 
      cst_HCP(3,8+ic)=  2*bush/bunb

      cst_HCP(1,9+ic)=  -1.0/bunb
      cst_HCP(2,9+ic)=  -rsq3/bunb 
      cst_HCP(3,9+ic)=  2*bush/bunb

      cst_HCP(1,10+ic)=  1.0/bunb
      cst_HCP(2,10+ic)=  -rsq3/bunb 
      cst_HCP(3,10+ic)=  2*bush/bunb

      cst_HCP(1,11+ic)=  2.0/bunb
      cst_HCP(2,11+ic)=  0.0 
      cst_HCP(3,11+ic)=  2*bush/bunb

      cst_HCP(1,12+ic)=  1.0/bunb
      cst_HCP(2,12+ic)=  rsq3/bunb 
      cst_HCP(3,12+ic)=  2*bush/bunb

      ic=ic+nslips_pra_HCP(4)

      !       2nd order <c+a> slip
      bunb=2.0*(sqrt((1.587)**2+1.0))
      bush=1.587
      cmt_HCP(1,1+ic)=  bush/bunb
      cmt_HCP(2,1+ic)=  -rsq3*bush/bunb 
      cmt_HCP(3,1+ic)=  2.0/bunb

      cmt_HCP(1,2+ic)= 2*bush/bunb
      cmt_HCP(2,2+ic)=  0.0
      cmt_HCP(3,2+ic)=  2.0/bunb

      cmt_HCP(1,3+ic)=  bush/bunb
      cmt_HCP(2,3+ic)=  rsq3*bush/bunb 
      cmt_HCP(3,3+ic)=  2.0/bunb

      cmt_HCP(1,4+ic)=  -bush/bunb
      cmt_HCP(2,4+ic)=  rsq3*bush/bunb 
      cmt_HCP(3,4+ic)=  2.0/bunb

      cmt_HCP(1,5+ic)= -2.0*bush/bunb
      cmt_HCP(2,5+ic)= 0.0 
      cmt_HCP(3,5+ic)=  2.0/bunb

      cmt_HCP(1,6+ic)=  -bush/bunb
      cmt_HCP(2,6+ic)=  -rsq3*bush/bunb 
      cmt_HCP(3,6+ic)=  2.0/bunb

      !----------------slip directions------------

      cst_HCP(1,1+ic)=  -1.0/bunb
      cst_HCP(2,1+ic)=  rsq3/bunb 
      cst_HCP(3,1+ic)=  2.0*bush/bunb

      cst_HCP(1,2+ic)=  -2.0/bunb
      cst_HCP(2,2+ic)=  0.0 
      cst_HCP(3,2+ic)=  2.0*bush/bunb

      cst_HCP(1,3+ic)=  -1.0/bunb
      cst_HCP(2,3+ic)=  -rsq3/bunb 
      cst_HCP(3,3+ic)=  2.0*bush/bunb

      cst_HCP(1,4+ic)=  1.0/bunb
      cst_HCP(2,4+ic)=  -rsq3/bunb 
      cst_HCP(3,4+ic)=  2.0*bush/bunb

      cst_HCP(1,5+ic)=  2.0/bunb
      cst_HCP(2,5+ic)=  0.0 
      cst_HCP(3,5+ic)=  2.0*bush/bunb

      cst_HCP(1,6+ic)=  1.0/bunb
      cst_HCP(2,6+ic)=  rsq3/bunb 
      cst_HCP(3,6+ic)=  2.0*bush/bunb
      ic=ic+nslips_pra_HCP(5)
      
      ! calculate relative cosines between slip system normal and slip directions
      call init_cosines_HCP() ! assumes nslip_HCP is initialized

      END SUBROUTINE
      
      ! calculates the cosines between slip system normals and slip directions
      ! taken from Kourosh's code
      subroutine init_cosines_HCP()
      implicit none      
      
      real(8) :: CTT(3,nslip_max)
      integer :: I,J,ISYS
      ! CALCULATING THE COSINE OF ANGLES BETWEEN M,N AND T
      CTT=0.D0
      DO ISYS=1,nslip_HCP
            CTT(1,ISYS)=cmt_HCP(2,ISYS)*cst_HCP(3,ISYS)-cmt_HCP(3,ISYS)* &
                       	cst_HCP(2,ISYS)
            CTT(2,ISYS)=cmt_HCP(3,ISYS)*cst_HCP(1,ISYS)-cmt_HCP(1,ISYS)* &
                       	cst_HCP(3,ISYS)
            CTT(3,ISYS)=cmt_HCP(1,ISYS)*cst_HCP(2,ISYS)-cmt_HCP(2,ISYS)* &
                       	cst_HCP(1,ISYS)	
      ENDDO
      
      DO I=1,nslip_HCP
         DO J=1,nslip_HCP
            cos_nm_HCP(I,J)=cmt_HCP(1,I)*cst_HCP(1,J)+cmt_HCP(2,I)*cst_HCP(2,J)+ &
                                    cmt_HCP(3,I)*cst_HCP(3,J) 
            cos_nt_HCP(I,J)=cmt_HCP(1,I)*CTT(1,J)+cmt_HCP(2,I)*CTT(2,J)+ &
                                    cmt_HCP(3,I)*CTT(3,J)
            cos_nn_HCP(I,J)=cmt_HCP(1,I)*cmt_HCP(1,J)+cmt_HCP(2,I)*cmt_HCP(2,J)+ &
                                    cmt_HCP(3,I)*cmt_HCP(3,J)

            sin_nm_HCP(I,J)=DSQRT(DABS(1.-cos_nm_HCP(I,J)*cos_nm_HCP(I,J)))
            sin_nt_HCP(I,J)=DSQRT(DABS(1.-cos_nt_HCP(I,J)*cos_nt_HCP(I,J)))
            sin_nn_HCP(I,J)=DSQRT(DABS(1.-cos_nn_HCP(I,J)*cos_nn_HCP(I,J)))
            
         ENDDO
      ENDDO
      
      
      
      end subroutine
      

      SUBROUTINE init_crys_BCC_DEKA()
     
      implicit none
      
      real(8), parameter :: R2=1.0D0/DSQRT(2.0D0), &
                            R3=1.0D0/DSQRT(3.0D0), &
                            R6=1.0D0/DSQRT(6.0D0), &
                            R26=2.0D0/DSQRT(6.0D0), &
                            R14=1.0D0/DSQRT(14.0D0), &
                            R214=2.0D0/DSQRT(14.0D0), &
                            R314=3.0D0/DSQRT(14.0D0)
      integer :: ic,I
      !-------BCC 48 slip systems -------
      nslip_BCC = 48

!-- SLIP PLANE
	DO  I=1,2
	   cmt_BCC(1,I)= 0.0
	   cmt_BCC(2,I)= R2
	   cmt_BCC(3,I)= R2
	enddo
	DO  I=3,4
	   cmt_BCC(1,I)= R2
	   cmt_BCC(2,I)= R2
	   cmt_BCC(3,I)= 0.0
	enddo
	DO  I=5,6
	   cmt_BCC(1,I)= R2
	   cmt_BCC(2,I)= 0.0
	   cmt_BCC(3,I)= R2
	enddo
	DO  I=7,8
	   cmt_BCC(1,I)= 0.0
	   cmt_BCC(2,I)= R2
	   cmt_BCC(3,I)=-R2
	enddo
	DO  I=9,10
	   cmt_BCC(1,I)=-R2
	   cmt_BCC(2,I)= R2
	   cmt_BCC(3,I)= 0.0
	enddo
	DO  I=11,12
	   cmt_BCC(1,I)= R2
	   cmt_BCC(2,I)= 0.0
	   cmt_BCC(3,I)=-R2
	enddo
!--number:13
        cmt_BCC(1,13)= R26
        cmt_BCC(2,13)= R6
        cmt_BCC(3,13)= R6
!--number:14
        cmt_BCC(1,14)= R6
        cmt_BCC(2,14)= R6
        cmt_BCC(3,14)= R26
!--number:15
        cmt_BCC(1,15)= R6
        cmt_BCC(2,15)= R26
        cmt_BCC(3,15)= R6
!--number:16
        cmt_BCC(1,16)=-R26
        cmt_BCC(2,16)= R6
        cmt_BCC(3,16)= R6
!--number:17
        cmt_BCC(1,17)=-R6
        cmt_BCC(2,17)= R6
        cmt_BCC(3,17)= R26
!--number:18
        cmt_BCC(1,18)=-R6
        cmt_BCC(2,18)= R26
        cmt_BCC(3,18)= R6
!--number:19
        cmt_BCC(1,19)= R26
        cmt_BCC(2,19)= R6
        cmt_BCC(3,19)=-R6
!--number:20
        cmt_BCC(1,20)= R6
        cmt_BCC(2,20)= R6
        cmt_BCC(3,20)=-R26
!--number:21
        cmt_BCC(1,21)= R6
        cmt_BCC(2,21)= R26
        cmt_BCC(3,21)=-R6
!--number:22
        cmt_BCC(1,22)= R26
        cmt_BCC(2,22)=-R6
        cmt_BCC(3,22)= R6
!--number:23
        cmt_BCC(1,23)= R6
        cmt_BCC(2,23)=-R6
        cmt_BCC(3,23)= R26
!--number:24
        cmt_BCC(1,24)= R6
        cmt_BCC(2,24)=-R26
        cmt_BCC(3,24)= R6
!--number:25
        cmt_BCC(1,25)= R314
        cmt_BCC(2,25)= R14
        cmt_BCC(3,25)= R214
!--number:26
        cmt_BCC(1,26)= R14
        cmt_BCC(2,26)= R214
        cmt_BCC(3,26)= R314
!--number:27
        cmt_BCC(1,27)= R214
        cmt_BCC(2,27)= R314
        cmt_BCC(3,27)= R14
!--number:28
        cmt_BCC(1,28)= R214
        cmt_BCC(2,28)= R14
        cmt_BCC(3,28)= R314
!--number:29
        cmt_BCC(1,29)= R14
        cmt_BCC(2,29)= R314
        cmt_BCC(3,29)= R214
!--number:30
        cmt_BCC(1,30)= R314
        cmt_BCC(2,30)= R214
        cmt_BCC(3,30)= R14
!--number:31
        cmt_BCC(1,31)=-R314
        cmt_BCC(2,31)= R14
        cmt_BCC(3,31)= R214
!--number:32
        cmt_BCC(1,32)=-R14
        cmt_BCC(2,32)= R214
        cmt_BCC(3,32)= R314
!--number:33
        cmt_BCC(1,33)=-R214
        cmt_BCC(2,33)= R314
        cmt_BCC(3,33)= R14
!--number:34
        cmt_BCC(1,34)=-R214
        cmt_BCC(2,34)= R14
        cmt_BCC(3,34)= R314
!--number:35
        cmt_BCC(1,35)=-R14
        cmt_BCC(2,35)= R314
        cmt_BCC(3,35)= R214
!--number:36
        cmt_BCC(1,36)=-R314
        cmt_BCC(2,36)= R214
        cmt_BCC(3,36)= R14
!--number:37
        cmt_BCC(1,37)= R314
        cmt_BCC(2,37)= R14
        cmt_BCC(3,37)=-R214
!--number:38
        cmt_BCC(1,38)= R14
        cmt_BCC(2,38)= R214
        cmt_BCC(3,38)=-R314
!--number:39
        cmt_BCC(1,39)= R214
        cmt_BCC(2,39)= R314
        cmt_BCC(3,39)=-R14
!--number:40
        cmt_BCC(1,40)= R214
        cmt_BCC(2,40)= R14
        cmt_BCC(3,40)=-R314
!--number:41
        cmt_BCC(1,41)= R14
        cmt_BCC(2,41)= R314
        cmt_BCC(3,41)=-R214
!--number:42
        cmt_BCC(1,42)= R314
        cmt_BCC(2,42)= R214
        cmt_BCC(3,42)=-R14
!--number:43
        cmt_BCC(1,43)= R314
        cmt_BCC(2,43)=-R14
        cmt_BCC(3,43)= R214
!--number:44
        cmt_BCC(1,44)= R14
        cmt_BCC(2,44)=-R214
        cmt_BCC(3,44)= R314
!--number:45
        cmt_BCC(1,45)= R214
        cmt_BCC(2,45)=-R314
        cmt_BCC(3,45)= R14
!--number:46
        cmt_BCC(1,46)= R214
        cmt_BCC(2,46)=-R14
        cmt_BCC(3,46)= R314
!--number:47
        cmt_BCC(1,47)= R14
        cmt_BCC(2,47)=-R314
        cmt_BCC(3,47)= R214
!--number:48
        cmt_BCC(1,48)= R314
        cmt_BCC(2,48)=-R214
        cmt_BCC(3,48)= R14

!--SLIP DIRECTION
!--number:1
        cst_BCC(1,1)= R3
        cst_BCC(2,1)= R3
        cst_BCC(3,1)=-R3
!--number:2
        cst_BCC(1,2)= R3
        cst_BCC(2,2)=-R3
        cst_BCC(3,2)= R3
!--number:3
        cst_BCC(1,3)=-R3
        cst_BCC(2,3)= R3
        cst_BCC(3,3)= R3
!--number:4
        cst_BCC(1,4)= R3
        cst_BCC(2,4)=-R3
        cst_BCC(3,4)= R3
!--number:5
        cst_BCC(1,5)=-R3
        cst_BCC(2,5)= R3
        cst_BCC(3,5)= R3
!--number:6
        cst_BCC(1,6)= R3
        cst_BCC(2,6)= R3
        cst_BCC(3,6)=-R3
!--number:7
        cst_BCC(1,7)= R3
        cst_BCC(2,7)= R3
        cst_BCC(3,7)= R3
!--number:8
        cst_BCC(1,8)=-R3
        cst_BCC(2,8)= R3
        cst_BCC(3,8)= R3
!--number:9
        cst_BCC(1,9)= R3
        cst_BCC(2,9)= R3
        cst_BCC(3,9)= R3
!--number:10
        cst_BCC(1,10)= R3
        cst_BCC(2,10)= R3
        cst_BCC(3,10)=-R3
!--number:11
        cst_BCC(1,11)= R3
        cst_BCC(2,11)= R3
        cst_BCC(3,11)= R3
!--number:12
        cst_BCC(1,12)= R3
        cst_BCC(2,12)=-R3
        cst_BCC(3,12)= R3
!--number:13
        cst_BCC(1,13)=-R3
        cst_BCC(2,13)= R3
        cst_BCC(3,13)= R3
!--number:14
        cst_BCC(1,14)= R3
        cst_BCC(2,14)= R3
        cst_BCC(3,14)=-R3
!--number:15
        cst_BCC(1,15)= R3
        cst_BCC(2,15)=-R3
        cst_BCC(3,15)= R3
!--number:16
        cst_BCC(1,16)= R3
        cst_BCC(2,16)= R3
        cst_BCC(3,16)= R3
!--number:17
        cst_BCC(1,17)= R3
        cst_BCC(2,17)=-R3
        cst_BCC(3,17)= R3
!--number:18
        cst_BCC(1,18)= R3
        cst_BCC(2,18)= R3
        cst_BCC(3,18)=-R3
!--number:19
        cst_BCC(1,19)= R3
        cst_BCC(2,19)=-R3
        cst_BCC(3,19)= R3
!--number:20
        cst_BCC(1,20)= R3
        cst_BCC(2,20)= R3
        cst_BCC(3,20)= R3
!--number:21
        cst_BCC(1,21)=-R3
        cst_BCC(2,21)= R3
        cst_BCC(3,21)= R3
!--number:22
        cst_BCC(1,22)= R3
        cst_BCC(2,22)= R3
        cst_BCC(3,22)=-R3
!--number:23
        cst_BCC(1,23)=-R3
        cst_BCC(2,23)= R3
        cst_BCC(3,23)= R3
!--number:24
        cst_BCC(1,24)= R3
        cst_BCC(2,24)= R3
        cst_BCC(3,24)= R3
!--number:25
        cst_BCC(1,25)=-R3
        cst_BCC(2,25)= R3
        cst_BCC(3,25)= R3
!--number:26
        cst_BCC(1,26)= R3
        cst_BCC(2,26)= R3
        cst_BCC(3,26)=-R3
!--number:27
        cst_BCC(1,27)= R3
        cst_BCC(2,27)=-R3
        cst_BCC(3,27)= R3
!--number:28
        cst_BCC(1,28)= R3
        cst_BCC(2,28)= R3
        cst_BCC(3,28)=-R3
!--number:29
        cst_BCC(1,29)= R3
        cst_BCC(2,29)=-R3
        cst_BCC(3,29)= R3
!--number:30
        cst_BCC(1,30)=-R3
        cst_BCC(2,30)= R3
        cst_BCC(3,30)= R3
!--number:31
        cst_BCC(1,31)= R3
        cst_BCC(2,31)= R3
        cst_BCC(3,31)= R3
!--number:32
        cst_BCC(1,32)= R3
        cst_BCC(2,32)=-R3
        cst_BCC(3,32)= R3
!--number:33
        cst_BCC(1,33)= R3
        cst_BCC(2,33)= R3
        cst_BCC(3,33)=-R3
!--number:34
        cst_BCC(1,34)= R3
        cst_BCC(2,34)=-R3
        cst_BCC(3,34)= R3
!--number:35
        cst_BCC(1,35)= R3
        cst_BCC(2,35)= R3
        cst_BCC(3,35)=-R3
!--number:36
        cst_BCC(1,36)= R3
        cst_BCC(2,36)= R3
        cst_BCC(3,36)= R3
!--number:37
        cst_BCC(1,37)= R3
        cst_BCC(2,37)=-R3
        cst_BCC(3,37)= R3
!--number:38
        cst_BCC(1,38)= R3
        cst_BCC(2,38)= R3
        cst_BCC(3,38)= R3
!--number:39
        cst_BCC(1,39)=-R3
        cst_BCC(2,39)= R3
        cst_BCC(3,39)= R3
!--number:40
        cst_BCC(1,40)= R3
        cst_BCC(2,40)= R3
        cst_BCC(3,40)= R3
!--number:41
        cst_BCC(1,41)=-R3
        cst_BCC(2,41)= R3
        cst_BCC(3,41)= R3
!--number:42
        cst_BCC(1,42)= R3
        cst_BCC(2,42)=-R3
        cst_BCC(3,42)= R3
!--number:43
        cst_BCC(1,43)= R3
        cst_BCC(2,43)= R3
        cst_BCC(3,43)=-R3
!--number:44
        cst_BCC(1,44)=-R3
        cst_BCC(2,44)= R3
        cst_BCC(3,44)= R3
!--number:45
        cst_BCC(1,45)= R3
        cst_BCC(2,45)= R3
        cst_BCC(3,45)= R3
!--number:46
        cst_BCC(1,46)=-R3
        cst_BCC(2,46)= R3
        cst_BCC(3,46)= R3
!--number:47
        cst_BCC(1,47)= R3
        cst_BCC(2,47)= R3
        cst_BCC(3,47)= R3
!--number:48
        cst_BCC(1,48)= R3
        cst_BCC(2,48)= R3
        cst_BCC(3,48)=-R3
!-- 	


      end subroutine
      
      SUBROUTINE init_crys_BCC()
     
      implicit none
      
      real(8), parameter :: R2=1.0D0/DSQRT(2.0D0), &
                            R3=1.0D0/DSQRT(3.0D0), &
                            R6=1.0D0/DSQRT(6.0D0), &
                            R26=2.0D0/DSQRT(6.0D0), &
                            R14=1.0D0/DSQRT(14.0D0), &
                            R214=2.0D0/DSQRT(14.0D0), &
                            R314=3.0D0/DSQRT(14.0D0)
      integer :: ic,i
      !-------BCC 48 slip systems -------
      nslip_BCC = 48

      !----------------------------------
      ! BCC slip families (used in indexing of g0_fam,xh_0,xh_inf,xg_0,xg_inf)
      ! [1] = all n={101}, all b=<111> systems
      ! [2] = all n={112}, b=(soft) directions (-)
      ! [3] = all n={112}, b=(hard) directions (-)
      ! [4] = all n={123}, b=(soft) directions (+)      
      ! [5] = all n={123}, b=(hard) directions (-)
      nslips_famBCC = 0
      nslips_famBCC(1)=12  ! {110}<111> type slip systems
      nslips_famBCC(2)=12  ! {112}<111> type slip systems - soft
      nslips_famBCC(3)=12  ! {112}<111> type slip systems - hard
      nslips_famBCC(4)=24  ! {112}<111> type slip systems - soft
      nslips_famBCC(5)=24  ! {112}<111> type slip systems - hard
  

      ! BCC SLIP PLANES:
      ! BCC - {110} type planes
      DO  I=1,2
         cmt_BCC(1,I)= 0.0
         cmt_BCC(2,I)= R2
         cmt_BCC(3,I)= R2
      enddo
      DO  I=3,4
         cmt_BCC(1,I)= R2
         cmt_BCC(2,I)= R2
         cmt_BCC(3,I)= 0.0
      enddo
      DO  I=5,6
         cmt_BCC(1,I)= R2
         cmt_BCC(2,I)= 0.0
         cmt_BCC(3,I)= R2
      enddo
      DO  I=7,8
         cmt_BCC(1,I)= 0.0
         cmt_BCC(2,I)= R2
         cmt_BCC(3,I)=-R2
      enddo
      DO  I=9,10
         cmt_BCC(1,I)=-R2
         cmt_BCC(2,I)= R2
         cmt_BCC(3,I)= 0.0
      enddo
      DO  I=11,12
         cmt_BCC(1,I)= R2
         cmt_BCC(2,I)= 0.0
         cmt_BCC(3,I)=-R2
      enddo
      
      ! BCC - {112} type planes
      !  number:13
      cmt_BCC(1,13)= R26
      cmt_BCC(2,13)= R6
      cmt_BCC(3,13)= R6
      !  number:14
      cmt_BCC(1,14)= R6
      cmt_BCC(2,14)= R6
      cmt_BCC(3,14)= R26
      !  number:15
      cmt_BCC(1,15)= R6
      cmt_BCC(2,15)= R26
      cmt_BCC(3,15)= R6
      !  number:16
      cmt_BCC(1,16)=-R26
      cmt_BCC(2,16)= R6
      cmt_BCC(3,16)= R6
      !  number:17
      cmt_BCC(1,17)=-R6
      cmt_BCC(2,17)= R6
      cmt_BCC(3,17)= R26
      !  number:18
      cmt_BCC(1,18)=-R6
      cmt_BCC(2,18)= R26
      cmt_BCC(3,18)= R6
      !  number:19
      cmt_BCC(1,19)= R26
      cmt_BCC(2,19)= R6
      cmt_BCC(3,19)=-R6
      !  number:20
      cmt_BCC(1,20)= R6
      cmt_BCC(2,20)= R6
      cmt_BCC(3,20)=-R26
      !  number:21
      cmt_BCC(1,21)= R6
      cmt_BCC(2,21)= R26
      cmt_BCC(3,21)=-R6
      !  number:22
      cmt_BCC(1,22)= R26
      cmt_BCC(2,22)=-R6
      cmt_BCC(3,22)= R6
      !  number:23
      cmt_BCC(1,23)= R6
      cmt_BCC(2,23)=-R6
      cmt_BCC(3,23)= R26
      !  number:24
      cmt_BCC(1,24)= R6
      cmt_BCC(2,24)=-R26
      cmt_BCC(3,24)= R6
      
      ! BCC - {123} type of planes
      !  number:25
      cmt_BCC(1,25)= R314
      cmt_BCC(2,25)= R14
      cmt_BCC(3,25)= R214
      !  number:26
      cmt_BCC(1,26)= R14
      cmt_BCC(2,26)= R214
      cmt_BCC(3,26)= R314
      !  number:27
      cmt_BCC(1,27)= R214
      cmt_BCC(2,27)= R314
      cmt_BCC(3,27)= R14
      !  number:28
      cmt_BCC(1,28)= R214
      cmt_BCC(2,28)= R14
      cmt_BCC(3,28)= R314
      !  number:29
      cmt_BCC(1,29)= R14
      cmt_BCC(2,29)= R314
      cmt_BCC(3,29)= R214
      !  number:30
      cmt_BCC(1,30)= R314
      cmt_BCC(2,30)= R214
      cmt_BCC(3,30)= R14
      !  number:31
      cmt_BCC(1,31)=-R314
      cmt_BCC(2,31)= R14
      cmt_BCC(3,31)= R214
      !  number:32
      cmt_BCC(1,32)=-R14
      cmt_BCC(2,32)= R214
      cmt_BCC(3,32)= R314
      !  number:33
      cmt_BCC(1,33)=-R214
      cmt_BCC(2,33)= R314
      cmt_BCC(3,33)= R14
      !  number:34
      cmt_BCC(1,34)=-R214
      cmt_BCC(2,34)= R14
      cmt_BCC(3,34)= R314
      !  number:35
      cmt_BCC(1,35)=-R14
      cmt_BCC(2,35)= R314
      cmt_BCC(3,35)= R214
      !  number:36
      cmt_BCC(1,36)=-R314
      cmt_BCC(2,36)= R214
      cmt_BCC(3,36)= R14
      !  number:37
      cmt_BCC(1,37)= R314
      cmt_BCC(2,37)= R14
      cmt_BCC(3,37)=-R214
      !  number:38
      cmt_BCC(1,38)= R14
      cmt_BCC(2,38)= R214
      cmt_BCC(3,38)=-R314
      !  number:39
      cmt_BCC(1,39)= R214
      cmt_BCC(2,39)= R314
      cmt_BCC(3,39)=-R14
      !  number:40
      cmt_BCC(1,40)= R214
      cmt_BCC(2,40)= R14
      cmt_BCC(3,40)=-R314
      !  number:41
      cmt_BCC(1,41)= R14
      cmt_BCC(2,41)= R314
      cmt_BCC(3,41)=-R214
      !  number:42
      cmt_BCC(1,42)= R314
      cmt_BCC(2,42)= R214
      cmt_BCC(3,42)=-R14
      !  number:43
      cmt_BCC(1,43)= R314
      cmt_BCC(2,43)=-R14
      cmt_BCC(3,43)= R214
      !  number:44
      cmt_BCC(1,44)= R14
      cmt_BCC(2,44)=-R214
      cmt_BCC(3,44)= R314
      !  number:45
      cmt_BCC(1,45)= R214
      cmt_BCC(2,45)=-R314
      cmt_BCC(3,45)= R14
      !  number:46
      cmt_BCC(1,46)= R214
      cmt_BCC(2,46)=-R14
      cmt_BCC(3,46)= R314
      !  number:47
      cmt_BCC(1,47)= R14
      cmt_BCC(2,47)=-R314
      cmt_BCC(3,47)= R214
      !  number:48
      cmt_BCC(1,48)= R314
      cmt_BCC(2,48)=-R214
      cmt_BCC(3,48)= R14
      
      !  SLIP DIRECTIONS
      !  number:1
      cst_BCC(1,1)= R3
      cst_BCC(2,1)= R3
      cst_BCC(3,1)=-R3
      !  number:2
      cst_BCC(1,2)= R3
      cst_BCC(2,2)=-R3
      cst_BCC(3,2)= R3
      !  number:3
      cst_BCC(1,3)=-R3
      cst_BCC(2,3)= R3
      cst_BCC(3,3)= R3
      !  number:4
      cst_BCC(1,4)= R3
      cst_BCC(2,4)=-R3
      cst_BCC(3,4)= R3
      !  number:5
      cst_BCC(1,5)=-R3
      cst_BCC(2,5)= R3
      cst_BCC(3,5)= R3
      !  number:6
      cst_BCC(1,6)= R3
      cst_BCC(2,6)= R3
      cst_BCC(3,6)=-R3
      !  number:7
      cst_BCC(1,7)= R3
      cst_BCC(2,7)= R3
      cst_BCC(3,7)= R3
      !  number:8
      cst_BCC(1,8)=-R3
      cst_BCC(2,8)= R3
      cst_BCC(3,8)= R3
      !  number:9
      cst_BCC(1,9)= R3
      cst_BCC(2,9)= R3
      cst_BCC(3,9)= R3
      !  number:10
      cst_BCC(1,10)= R3
      cst_BCC(2,10)= R3
      cst_BCC(3,10)=-R3
      !  number:11
      cst_BCC(1,11)= R3
      cst_BCC(2,11)= R3
      cst_BCC(3,11)= R3
      !  number:12
      cst_BCC(1,12)= R3
      cst_BCC(2,12)=-R3
      cst_BCC(3,12)= R3
      !  number:13               *b direction modified to point towards soft dir
      cst_BCC(1,13)= R3
      cst_BCC(2,13)=-R3
      cst_BCC(3,13)=-R3
      !  number:14               *b direction modified to point towards soft dir
      cst_BCC(1,14)=-R3
      cst_BCC(2,14)=-R3
      cst_BCC(3,14)=+R3
      !  number:15               *b direction modified to point towards soft dir
      cst_BCC(1,15)=-R3
      cst_BCC(2,15)=+R3
      cst_BCC(3,15)=-R3
      !  number:16               *b direction modified to point towards soft dir
      cst_BCC(1,16)=-R3
      cst_BCC(2,16)=-R3
      cst_BCC(3,16)=-R3
      !  number:17
      cst_BCC(1,17)= R3
      cst_BCC(2,17)=-R3
      cst_BCC(3,17)= R3
      !  number:18
      cst_BCC(1,18)= R3
      cst_BCC(2,18)= R3
      cst_BCC(3,18)=-R3
      !  number:19
      cst_BCC(1,19)= R3
      cst_BCC(2,19)=-R3
      cst_BCC(3,19)= R3
      !  number:20               *b direction modified to point towards soft dir
      cst_BCC(1,20)=-R3
      cst_BCC(2,20)=-R3
      cst_BCC(3,20)=-R3
      !  number:21
      cst_BCC(1,21)=-R3
      cst_BCC(2,21)= R3
      cst_BCC(3,21)= R3
      !  number:22
      cst_BCC(1,22)= R3
      cst_BCC(2,22)= R3
      cst_BCC(3,22)=-R3
      !  number:23
      cst_BCC(1,23)=-R3
      cst_BCC(2,23)= R3
      cst_BCC(3,23)= R3
      !  number:24               *b direction modified to point towards soft dir
      cst_BCC(1,24)=-R3
      cst_BCC(2,24)=-R3
      cst_BCC(3,24)=-R3
      !  number:25               *b direction modified to point towards soft dir
      cst_BCC(1,25)=+R3
      cst_BCC(2,25)=-R3
      cst_BCC(3,25)=-R3
      !  number:26               *b direction modified to point towards soft dir
      cst_BCC(1,26)=-R3
      cst_BCC(2,26)=-R3
      cst_BCC(3,26)=+R3
      !  number:27               *b direction modified to point towards soft dir
      cst_BCC(1,27)=-R3
      cst_BCC(2,27)=+R3
      cst_BCC(3,27)=-R3
      !  number:28               *b direction modified to point towards soft dir
      cst_BCC(1,28)=-R3
      cst_BCC(2,28)=-R3
      cst_BCC(3,28)=+R3
      !  number:29               *b direction modified to point towards soft dir
      cst_BCC(1,29)=-R3
      cst_BCC(2,29)=+R3
      cst_BCC(3,29)=-R3
      !  number:30               *b direction modified to point towards soft dir
      cst_BCC(1,30)=+R3
      cst_BCC(2,30)=-R3
      cst_BCC(3,30)=-R3
      !  number:31               *b direction modified to point towards soft dir
      cst_BCC(1,31)=-R3
      cst_BCC(2,31)=-R3
      cst_BCC(3,31)=-R3
      !  number:32
      cst_BCC(1,32)= R3
      cst_BCC(2,32)=-R3
      cst_BCC(3,32)= R3
      !  number:33
      cst_BCC(1,33)= R3
      cst_BCC(2,33)= R3
      cst_BCC(3,33)=-R3
      !  number:34
      cst_BCC(1,34)= R3
      cst_BCC(2,34)=-R3
      cst_BCC(3,34)= R3
      !  number:35
      cst_BCC(1,35)= R3
      cst_BCC(2,35)= R3
      cst_BCC(3,35)=-R3
      !  number:36               *b direction modified to point towards soft dir
      cst_BCC(1,36)=-R3
      cst_BCC(2,36)=-R3
      cst_BCC(3,36)=-R3
      !  number:37
      cst_BCC(1,37)= R3
      cst_BCC(2,37)=-R3
      cst_BCC(3,37)= R3
      !  number:38               *b direction modified to point towards soft dir
      cst_BCC(1,38)=-R3
      cst_BCC(2,38)=-R3
      cst_BCC(3,38)=-R3
      !  number:39
      cst_BCC(1,39)=-R3
      cst_BCC(2,39)= R3
      cst_BCC(3,39)= R3
      !  number:40               *b direction modified to point towards soft dir
      cst_BCC(1,40)=-R3
      cst_BCC(2,40)=-R3
      cst_BCC(3,40)=-R3
      !  number:41
      cst_BCC(1,41)=-R3
      cst_BCC(2,41)= R3
      cst_BCC(3,41)= R3
      !  number:42
      cst_BCC(1,42)= R3
      cst_BCC(2,42)=-R3
      cst_BCC(3,42)= R3
      !  number:43
      cst_BCC(1,43)= R3
      cst_BCC(2,43)= R3
      cst_BCC(3,43)=-R3
      !  number:44
      cst_BCC(1,44)=-R3
      cst_BCC(2,44)= R3
      cst_BCC(3,44)= R3
      !  number:45               *b direction modified to point towards soft dir
      cst_BCC(1,45)=-R3
      cst_BCC(2,45)=-R3
      cst_BCC(3,45)=-R3
      !  number:46
      cst_BCC(1,46)=-R3
      cst_BCC(2,46)= R3
      cst_BCC(3,46)= R3
      !  number:47               *b direction modified to point towards soft dir
      cst_BCC(1,47)=-R3
      cst_BCC(2,47)=-R3
      cst_BCC(3,47)=-R3
      !  number:48
      cst_BCC(1,48)= R3
      cst_BCC(2,48)= R3
      cst_BCC(3,48)=-R3
      END SUBROUTINE
      
      SUBROUTINE crys_ElastMat_HCP(C_crys,c11,c12,c13,c33,c55)
      implicit none 
      real(8), intent(in) :: c11,c12,c13,c33,c55 ! HCP elasticity coeffs
      real(8), intent(out):: C_crys(6,6)         ! 2D elasticity matrix
      ! HCP
      C_crys = 0.D0
      C_crys(1,1)=c11   ! props(1)
      C_crys(2,2)=c11   ! props(1)
      C_crys(3,3)=c33   ! props(4)
      C_crys(4,4)=(c11-c12)      !props(1)-props(2)
      C_crys(5,5)=c55*2.d0       !props(5)*2
      C_crys(6,6)=c55*2.d0       !props(5)*2
      C_crys(1,2)=c12   ! props(2)
      C_crys(2,1)=c12   ! props(2)
      C_crys(1,3)=c13   ! props(3)
      C_crys(3,1)=c13   ! props(3)
      C_crys(2,3)=c13   ! props(3)
      C_crys(3,2)=c13   ! props(3)
      
      !c11_HCP=props(1)
      !c12_HCP=props(2)
      !c13_HCP=props(3)
      !c33_HCP=props(4)
      !c55_HCP=props(5)
      
      ! CMRL:
      ! C11 = CMRLprops(1)
      ! C22 = CMRLprops(2)
      ! C33 = CMRLprops(3)
      ! C44 = CMRLprops(4) x 2
      ! C55 = CMRLprops(5) x 2
      ! C66 = CMRLprops(6) x 2
      ! C12 = CMRLprops(7)
      ! C13 = CMRLprops(8)
      ! C23 = CMRLprops(9)
      
      END SUBROUTINE
      
      SUBROUTINE crys_ElastMat_BCC(C_crys,c1,c2,c3) 
      implicit none
      real(8), intent(in) :: c1,c2,c3 ! BCC elasticity coeffs
      real(8), intent(out):: C_crys(6,6)  ! 2D elasticity matrix
      ! BCC
      C_crys=0.D0
      C_crys(1,1)=c1
      C_crys(2,2)=c1
      C_crys(3,3)=c1
      C_crys(4,4)=c3
      C_crys(5,5)=c3
      C_crys(6,6)=c3
      C_crys(1,2)=c2
      C_crys(2,1)=c2
      C_crys(1,3)=c2
      C_crys(3,1)=c2
      C_crys(2,3)=c2
      C_crys(3,2)=c2
      END SUBROUTINE
      
      ! assigns a variable related to slip families to corresponding slip systems      
      SUBROUTINE slipFamily_to_SystemsHCP(var_fam,var_slipSys,nslips_fam)
      implicit none

      real(8), intent(out):: var_slipSys(nslip_HCP)
      real(8), intent(in) :: var_fam(nslipfam)
      integer, intent(in) :: nslips_fam(nslipfam)
      integer :: iFam,posFam,iSys
      posFam=0
      do iFam=1,nslipfam
         var_slipSys(posFam+1: &
                     posFam+nslips_fam(iFam)) = var_fam(iFam)   ! assign the slip-family variable to all corresponding slip systems
         posFam=posFam+nslips_fam(iFam)
      enddo
      END SUBROUTINE
         
      SUBROUTINE slipFamily_to_SystemsBCC(var_fam,var_slipSys,nslips_fam)
      implicit none

      real(8), intent(out):: var_slipSys(nslip_max,2)
      real(8), intent(in) :: var_fam(nslipfam)
      integer, intent(in) :: nslips_fam(nslipfam)
      integer :: iFam,iSys
      
      ! {110}<111> systems
      var_slipSys(1:12,1) = var_fam(1)
      var_slipSys(1:12,2) = var_fam(1)
      
      ! {112}<111> systems
      var_slipSys(13:24,1) = var_fam(2)   ! soft direction
      var_slipSys(13:24,2) = var_fam(3)   ! hard direction
      
      ! {123}<111> systems
      var_slipSys(25:48,1) = var_fam(4)   ! soft direction
      var_slipSys(25:48,2) = var_fam(5)   ! hard direction
      
      END SUBROUTINE
      
      SUBROUTINE slipFamily_to_SystemsBCC_b1(var_fam,var_slipSys,nslips_fam)
      implicit none

      real(8), intent(out):: var_slipSys(nslip_max,2)
      real(8), intent(in) :: var_fam(nslipfam)
      integer, intent(in) :: nslips_fam(nslipfam)
      integer :: iFam,iSys
      ! locals
      real(8), parameter :: R3=1.0D0/DSQRT(3.0D0)
      real(8), parameter :: b1(3) = (/-R3,R3,R3/) ! b1=[-1, 1, 1]
      real(8), parameter :: b2(3) = (/R3,R3,-R3/) ! b2=[ 1, 1,-1]
      real(8) :: dot_b1, dot_b2

      ! {110} systems
      do iSys=1,12
         dot_b1 = dabs(DOT_PRODUCT(cst_BCC(:,iSys),b1(:)))
         if(dot_b1 > 0.95 .and. dot_b1 < 1.05) then
            var_slipSys(iSys,1) = var_fam(1)
            var_slipSys(iSys,2) = var_fam(1)
         endif
      enddo
      ! {112} systems
      do iSys=13,24
         dot_b1 = dabs(DOT_PRODUCT(cst_BCC(:,iSys),b1(:)))
         if(dot_b1 > 0.95 .and. dot_b1 < 1.05) then
            var_slipSys(iSys,1) = var_fam(2)
            var_slipSys(iSys,2) = var_fam(3)
         endif
      enddo
      ! {123} systems
      do iSys=25,48
         dot_b1 = dabs(DOT_PRODUCT(cst_BCC(:,iSys),b1(:)))
         if(dot_b1 > 0.95 .and. dot_b1 < 1.05) then
            var_slipSys(iSys,1) = var_fam(4)
            var_slipSys(iSys,2) = var_fam(5)
         endif
      enddo
      
      END SUBROUTINE
      
      END MODULE
      