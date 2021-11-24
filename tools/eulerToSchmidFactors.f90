      program eul_rot
      use crystal_Titanium
      
      implicit none
      
      character :: eulerDirection
            
      real(8), parameter :: PI = 3.14159265359D0
      real(8) :: euler(3), rot_HCP(3,3), rot_BCCtoHCP(3,3), rot_BCC(3,3)
      integer :: i,j,dir,isys
      real(8) :: gst_HCP(3,nslip_max),gmt_HCP(3,nslip_max), sch(3,3,nslip_max)
      real(8) :: gst_BCC(3,nslip_max),gmt_BCC(3,nslip_max)
      real(8) :: detRot
            
      write(*,*) 'Enter euler angles [phi 1, theta, phi 2] (in radians):'
      read(*,*) euler(1:3)
!      euler(:) = euler(:)*2.D0*PI/360.D0
      
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

      call getRotationMatrix(euler,rot_HCP)
      if (eulerDirection=='P') then
         rot_HCP = TRANSPOSE(rot_HCP)
      endif
      CALL calcDeterminant(rot_HCP,detRot)
      write(*,*) 'determinant:',detRot
      
      
      write(*,*) 'Active rotation matrix for the euler angles:'
      do i=1,3
         write(*,'(3F10.3)') rot_HCP(i,1:3)
      enddo

      write(*,*) 'Passive rotation (coordinate transformation) matrix for the euler angles:'
      do i=1,3
         write(*,'(3F10.3)') rot_HCP(1:3,i)
      enddo
      
      ! calculate and print Schmidt factors      
      call init_crys_HCP()
      call init_crys_BCC()
      
      do isys=1,nslip_HCP
         gst_HCP(:,isys) = MATMUL(rot_HCP,cst_HCP(:,isys))  ! slip directions
         gmt_HCP(:,isys) = MATMUL(rot_HCP,cmt_HCP(:,isys))  ! slip plane normals
      enddo
      call create_schmidt(sch,gst_HCP,gmt_HCP,nslip_max)
      do dir=1,3
         call print_schmidt(sch,dir)
      enddo
      
      call create_BCCtoHCP(rot_BCCtoHCP)
      
      rot_BCC = MATMUL(rot_HCP,rot_BCCtoHCP)
      do isys=1,nslip_BCC
         gst_BCC(:,isys) = MATMUL(rot_BCC,cst_BCC(:,isys))  ! slip directions
         gmt_BCC(:,isys) = MATMUL(rot_BCC,cmt_BCC(:,isys))  ! slip plane normals
      enddo
      
      write(*,*) 'HCP basal plane','isys=1','basal'
      write(*,*) gmt_HCP(:,1)
      write(*,*)
      write(*,*) 'HCP slip direction','isys=1','a3-basal'
      write(*,*) gst_HCP(:,1)
      write(*,*) 'HCP slip direction','isys=2','a1-basal'
      write(*,*) gst_HCP(:,2)
      write(*,*) 'HCP slip direction','isys=3','a2-basal'
      write(*,*) gst_HCP(:,3)
      write(*,*) 'HCP slip direction','isys=4','a2-prism ??'
      write(*,*) gst_HCP(:,4)
      write(*,*) 'HCP slip direction','isys=5','a1-basal ??'
      write(*,*) gst_HCP(:,5)
      write(*,*) 'HCP slip direction','isys=6','a3-basal ??'
      write(*,*) gst_HCP(:,6)
      write(*,*)
      write(*,*) 'BCC (101) plane','isys=5','must be in line with HCP-basal'
      write(*,*) gmt_BCC(:,5)
      write(*,*) 'BCC slip direction','isys=5','b1 ?? '
      write(*,*) gst_BCC(:,5)
      write(*,*) 'BCC slip direction','isys=6','b2 ??'
      write(*,*) gst_BCC(:,6)
      

      
      
      
      end program

!    THIS SUBROUTINE CALCULATES THE TRANSFORMATION MATRIX
!----------------------------------------------------------
!     phi   - euler(1)
!     theta - euler(2)
!     omega - euler(3)
!
! deniz - comment
!
! output, tlgt,:
!   an active transformation matrix corresponding to the actual rotation due to euler angles {v}_i = TLGT_ij {v0}_j
!   as a (passive) coordinate transformation matrix, it transforms coordinates from iCrystal to global spatial frame. {v}_i = TLGT_ij {v_cry}_j
!---------------------------------------------------
      subroutine getRotationMatrix(euler,rot)   ! Active
      implicit double precision (a-h,o-z)  
      real(8), intent(out) :: rot(3,3)
      
      real(8) :: euler(3),tlg(3,3),tlgt(3,3)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)

  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      ! tlg is a passive rotation correspt to euler angles
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct
      
      ! tlgt is an active rotation correspt to euler angles
      tlgt = TRANSPOSE(tlg)
      
      rot = tlgt  ! Active
      
      end subroutine 

!---------------------------------------------------------
      
      subroutine create_schmidt(sch,cst,cmt,nslip)
      implicit none
      
      real(8),intent(out):: sch(3,3,nslip)
      real(8),intent(in) :: cmt(3,nslip),cst(3,nslip)
      integer,intent(in) :: nslip
      
      integer :: i,j,isys
      ! deniz - comment
      ! here the schmidt tensors are formed, using the slip direction and plane normal vectors
      sch=0.0d0      
      do isys=1,nslip
         do j=1,3
          do i=1,3
             sch(i,j,isys)=cst(i,isys)*cmt(j,isys)
          enddo
         enddo
      enddo
      end subroutine

      SUBROUTINE print_schmidt(sch,dir)
      implicit none

      integer, parameter :: nslip=30
      real(8), intent(in) :: sch(3,3,nslip)
      integer, intent(in) :: dir

      character(len=3) :: strDir
      strDir = 'XYZ'

      write(*,*) 'Schmidt Factors for loading along direction ',strDir(dir:dir)
      write(*,*) '      <a3>basal      <a1>basal      <a2>basal      <a2>prismatic  <a3>prismatic  <a1>prismatic'
      write(*,'(6F15.3)') abs(sch(dir,dir,1)),abs(sch(dir,dir,2)),abs(sch(dir,dir,3)),  &  ! 3 <a>basals
                                 abs(sch(dir,dir,4)),abs(sch(dir,dir,5)),abs(sch(dir,dir,6))      ! 3 <a>prismatics
      write(*,*) 'MAX among <a>basals:', MAXVAL(abs(sch(dir,dir,1:3)))
      write(*,*) 'MAX among <a>prismatics:', MAXVAL(abs(sch(dir,dir,4:6)))
      write(*,*) 'MAX SCHMIDT FACTOR among all:', MAXVAL(abs(sch(dir,dir,1:6)))
      write(*,*) ''

      END SUBROUTINE
      
      subroutine create_BCCtoHCP(rot_BCCtoHCP)
      implicit none
      
      real(8), intent(out) :: rot_BCCtoHCP(3,3)

      real(8) :: bc1,bc2,bc3,bc4,bc5
   
      bc1=dsqrt(2.d0)
      bc2=2.d0*dsqrt(3.d0)
      bc3=bc2**2.d0+6.d0
      bc4=dsqrt(bc3)
      bc5=(bc4-bc2)/12.d0             
      rot_BCCtoHCP(1,1)=bc5
      rot_BCCtoHCP(1,2)=2.d0*bc5+dsqrt(3.d0)/2.d0
      rot_BCCtoHCP(1,3)=-bc5
      rot_BCCtoHCP(3,1)=1.d0/bc1
      rot_BCCtoHCP(3,2)=0.0d0
      rot_BCCtoHCP(3,3)=1.d0/bc1
      rot_BCCtoHCP(2,1)=rot_BCCtoHCP(3,2)*rot_BCCtoHCP(1,3)-  &
        rot_BCCtoHCP(1,2)*rot_BCCtoHCP(3,3)
      rot_BCCtoHCP(2,2)=rot_BCCtoHCP(3,3)*rot_BCCtoHCP(1,1)-  &
        rot_BCCtoHCP(3,1)*rot_BCCtoHCP(1,3)
      rot_BCCtoHCP(2,3)=rot_BCCtoHCP(3,1)*rot_BCCtoHCP(1,2)-  &
        rot_BCCtoHCP(1,1)*rot_BCCtoHCP(3,2)       
                 
      end subroutine
      

      subroutine calcDeterminant(a,det)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      det=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(2,1)*a(1,2)  &
           *a(3,3)+a(2,1)*a(1,3)*a(3,2)+a(3,1)*a(1,2)*a(2,3)-a(3,1)  &
           *a(1,3)*a(2,2))
      end
