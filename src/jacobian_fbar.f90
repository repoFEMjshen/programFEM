      ! ************************************************************************************************ !
      !
      !  Description:   This subroutine computes the Material Jacobian defined as dP/dF where P and F
      !                 correspond to the 1st PK stress and total deformation gradient, respectively 
      !
      !       Inputs:       
      !
      !                    
      !      Outputs:        
      ! 
      !      Details:   
      !
      !   References:   
      !
      !     Creation:   Ahmad Shahba  | 2016/02/04 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      subroutine computeTensorElasticityFirstCP                   &  
      (                                                           &
         tensorElasticityFirst,                                  &
         stressSecondPkThermalVoigtNew,                           &
         deformationGradientPlasticNew,                           &
         deformationGradientNew,                                  &
         deformationGradientPlasticOld,                           &
         tensorElasticityGlobal,                                  &
         tensorSchmidGlobal,                                      &
         deformationGradientThermalNew,                           &
         deformationGradientElasticNew,                           &
         tensorResolverStress,                                    &
         dIncSlipPlastic_dStressResolved,                         &
         stressCorrectorSecondPk,                                 &
         incrementSlipPlastic,                                    &
         nSlipSystems                                             &
      )




      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


         implicit none

      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
         real(8), parameter :: &
            ZERO = 0.d0,                            &
            ONE = 1.d0,                             &
            HALF = 0.5d0
            
      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------


      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

            

      !  Inputs ------------------------------------------------------------------------------------------
         integer, intent(in) ::                    &
            nSlipSystems
            
         real(8), intent(in) ::                       &
            deformationGradientNew(3,3),                       &
            deformationGradientPlasticOld(3,3),                &
            tensorElasticityGlobal(3,3,3,3),                   &
            tensorSchmidGlobal(3,3,nSlipSystems),              &
            deformationGradientThermalNew(3,3),                &
            deformationGradientElasticNew(3,3),                &
            tensorResolverStress(3,3,nSlipSystems),            &
            dIncSlipPlastic_dStressResolved(nSlipSystems),                &
            stressCorrectorSecondPk(3,3,nSlipSystems),                    &
            incrementSlipPlastic(nSlipSystems),                           &
            stressSecondPkThermalVoigtNew(6),                  &
            deformationGradientPlasticNew(3,3)

      !  Outputs -----------------------------------------------------------------------------------------

         real(8), intent(out)  ::                     &
            tensorElasticityFirst(3,3,3,3)
       


      ! --------------------------------------------------------------------------------------------------
      !  Local Constants ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
      
         include 'jacobian_fbar.inc' ! function interfaces
         
         real(8) :: matrixIdentity3x3(3,3)            
            

         integer ::                                &
            ii,                                                &
            jj,                                                &
            kk,                                                &
            ll,                                                &
            mm,                                                &
            nn,                                                &
            iSlipSystem,                                       &
            indexRow,                                          &
            indexColumn
            
         integer :: iflaginv,ierr
            
         integer ::                                &
            map4dTo2dFormat(81, 4)
            
         real(8) ::                                   &
            tensorIndentityTranspose(3, 3, 3, 3),              &
            tensorFactorO(3, 3, 3, 3),                         &
            tensorFactorP(3, 3, 3, 3),                         &
            tensorFactorQ(3, 3),                               &
            tensorFactorR(3, 3),                               &
            dFtF_dF(3, 3, 3, 3),                               &
            dStress2PkThermal_dF(3, 3, 3, 3),                  &
            deformationGradPlasticOldInv(3, 3),                &
            dAtilda_dF(3, 3, 3, 3),                            &
            deformationGradThermalNewInv(3, 3),                &
            dStressCorrector_dF(3, 3, 3, 3, nSlipSystems),     &
            dStressTrial_dF(3, 3, 3, 3),                       &
            dIncrementSlip_dStress2PK(3, 3, nSlipSystems),     &
            tensorIndentity(3, 3, 3, 3),                       &
            vectorRHS(81),                                     &
            matrixCoefficients(81, 81),                        &
            matrixCoefficientsInv(81, 81),                     &  ! REMOVE THIS
            dStress2PkThermal_dFStack(81),                     &
            dIncrementSlip_dF(3, 3),                           &
            dDeformationGradPlasticInv_dF(3, 3, 3, 3),         &
            dStress2Pk_dF(3, 3, 3, 3),                         &
            stressSecondPkThermalNew(3, 3),                    &
            deformationGradPlasticNewInv(3, 3),                &
            Stress2Pk(3, 3)

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !  =================================== Initialization =============================================
         tensorElasticityFirst(:,:,:,:) = ZERO
         
         call cal_ckrone_2d(matrixIdentity3x3)  
         
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of A_tilda w.r.t. F  ------------------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         tensorIndentityTranspose(:,:,:,:) = ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     if ( ii == ll .and. jj == kk) tensorIndentityTranspose(ii, jj, kk, ll) = ONE
                  end do
               end do
            end do
         end do 
         
         
      !  Calculate derivative of (F^t*F) w.r.t. F
         tensorFactorP = computeProductTensorLower(matrixIdentity3x3,                                    &
                                                                     transpose (deformationGradientNew))+& 
                         computeProductTensorUpper(transpose(deformationGradientNew),                   &
                                                                                      matrixIdentity3x3)
         dFtF_dF = contractTensorsDouble4d4d(tensorFactorP, tensorIndentityTranspose)                                                                                
         
      !  Calculate Lower product of plastic deformation gradient at time t
         !deformationGradPlasticOldInv = computeMatrixInverse(deformationGradientPlasticOld)
         deformationGradPlasticOldInv = deformationGradientPlasticOld
         !call la_matinv(deformationGradPlasticOldInv,3,3,iflaginv)
         call lapack_invert(3,deformationGradPlasticOldInv,iflaginv)
         
         tensorFactorP = computeProductTensorLower(transpose(deformationGradPlasticOldInv),              &
                                                                  transpose(deformationGradPlasticOldInv))
                                                                  
      !  alculate derivative of A_tilda w.r.t. F
         dAtilda_dF = contractTensorsDouble4d4d( tensorFactorP, dFtF_dF)
       
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of stress correcter tensor w.r.t. F ---------------------------------------
      !  ------------------------------------------------------------------------------------------------
         !deformationGradThermalNewInv = computeMatrixInverse(deformationGradientThermalNew)
         deformationGradThermalNewInv = deformationGradientThermalNew
         !call la_matinv(deformationGradThermalNewInv,3,3,iflaginv)
         call lapack_invert(3,deformationGradThermalNewInv,iflaginv)
         
         do iSlipSystem = 1, nSlipSystems
         
            tensorFactorP = computeProductTensorLower(transpose(deformationGradThermalNewInv),           &
                                                   transpose(matmul(tensorSchmidGlobal(:,:, iSlipSystem),&   
                                                                    deformationGradThermalNewInv)))
            dStressCorrector_dF(:,:,:,:, iSlipSystem) =                                                  &
               contractTensorsDouble4d4d(tensorElasticityGlobal ,                                      &
                                         contractTensorsDouble4d4d(tensorFactorP, dAtilda_dF))                                                         
         
         end do

      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of trial stress w.r.t. F   ------------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         tensorFactorP = computeProductTensorLower(transpose(deformationGradThermalNewInv),              &
                                                                  transpose(deformationGradThermalNewInv))
         dStressTrial_dF = HALF *                                                                        &
            contractTensorsDouble4d4d(tensorElasticityGlobal,                                          &
                                      contractTensorsDouble4d4d(tensorFactorP, dAtilda_dF))
                                      
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of plastic shear strain increment w.r.t. 2ndPK Stress in Thermally-expanded
      !  configuration  ---------------------------------------------------------------------------------
      !  ------------------------------------------------------------------------------------------------                                
         do iSlipSystem = 1, nSlipSystems
            tensorFactorQ = matmul(matmul(transpose(deformationGradientElasticNew),                      &
                                                    deformationGradientElasticNew),                      &
                                   tensorResolverStress(:,:, iSlipSystem)) 
            dIncrementSlip_dStress2PK(:,:, iSlipSystem) = dIncSlipPlastic_dStressResolved(iSlipSystem) * &
                                                          HALF *(tensorFactorQ + transpose(tensorFactorQ))
         end do
         
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of 2nd PK stress in thermally-activated configuration w.r.t. F   ----------
      !  ------------------------------------------------------------------------------------------------
         tensorIndentity(:,:,:,:) = ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     if ( ii == kk .and. jj == ll) tensorIndentity(ii, jj, kk, ll) = ONE
                  end do
               end do
            end do
         end do

      !  Setting the 4th order coefficient and RHS tensors   
         tensorFactorO(:,:,:,:) = ZERO
         tensorFactorP(:,:,:,:) = ZERO
         do iSlipSystem = 1, nSlipSystems
            tensorFactorO(:,:,:,:) = tensorFactorO(:,:,:,:) + incrementSlipPlastic(iSlipSystem) *        &
                                                                 dStressCorrector_dF(:,:,:,:, iSlipSystem)
         
            tensorFactorP(:,:,:,:) = tensorFactorP(:,:,:,:) +                                            &
                              computeProductDyadic2d2d(stressCorrectorSecondPk(:,:, iSlipSystem),        &
                                                        dIncrementSlip_dStress2PK(:,:, iSlipSystem)) 
         end do

         tensorFactorP(:,:,:,:) = tensorFactorP(:,:,:,:) + tensorIndentity(:,:,:,:)
         tensorFactorO(:,:,:,:) = dStressTrial_dF(:,:,:,:) - tensorFactorO(:,:,:,:) 

      !  Converting the 4th order tensors into 2d format   
         map4dTo2dFormat(:,:) = 0
         vectorRHS(:) = ZERO
         matrixCoefficients(:,:)= ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     indexRow = ll + (kk - 1) * 3 + (jj - 1) * 9 + (ii - 1) * 27
                     map4dTo2dFormat(indexRow, 1) = ii
                     map4dTo2dFormat(indexRow, 2) = jj
                     map4dTo2dFormat(indexRow, 3) = kk
                     map4dTo2dFormat(indexRow, 4) = ll
                     
                     vectorRHS(indexRow) =  tensorFactorO(ii, jj, kk, ll)
                     
                     do mm = 1, 3
                        do nn = 1, 3
                           indexColumn = ll + (kk - 1) * 3 + (nn - 1) * 9 + (mm - 1) * 27
                           matrixCoefficients(indexRow, indexColumn) = tensorFactorP(ii, jj, mm, nn)
                        end do
                     end do
                     
                  end do
               end do
            end do
         end do

         !dStress2PkThermal_dFStack = solveMatrixEquation(matrixCoefficients, vectorRHS)
         dStress2PkThermal_dFStack = vectorRHS
         matrixCoefficientsInv = matrixCoefficients
         call lapack_invert(81,matrixCoefficientsInv,ierr)
         vectorRHS = matmul(matrixCoefficientsInv,vectorRHS)
         !CALL solveAXb(matrixCoefficients,81,vectorRHS,ierr)
         if (ierr /= 0) then
            OPEN(587,FILE='lapack_error.log')
            write(587,*) 'error code:',ierr
            !CALL printMatrix(matrixCoefficients,81,81, &
            !         'lapack solve error - matrix:',587)
            CALL printMatrix(matrixCoefficientsInv,81,81, &
                     'lapack solve error - matrix:',587)
            write(587,*) 'RHS after call:'
            write(587,*) vectorRHS
            close(587)
            write(*,*) 'lapack error -- see file: lapack_error.log -- error code: ierr'
            CALL MPI_ABORT(ierr)
         endif
         
         
      !  Converting Stack Version To 4th order format
         dStress2PkThermal_dF = ZERO
         do indexRow = 1, 81
            ii = map4dTo2dFormat(indexRow, 1)
            jj = map4dTo2dFormat(indexRow, 2)
            kk = map4dTo2dFormat(indexRow, 3)
            ll = map4dTo2dFormat(indexRow, 4)
            
            dStress2PkThermal_dF(ii, jj, kk, ll) = dStress2PkThermal_dFStack(indexRow)
         end do

      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of Inverse Plastic Deformation Griadent w.r.t. F   ------------------------
      !  ------------------------------------------------------------------------------------------------   
         dDeformationGradPlasticInv_dF(:,:,:,:) = ZERO
         do iSlipSystem = 1, nSlipSystems
            tensorFactorQ(:,:) = - matmul(deformationGradPlasticOldInv,                                  &
                                          tensorSchmidGlobal(:,:, iSlipSystem))
            
            dIncrementSlip_dF(:,:) = contractTensorsDouble2d4d(                                          &
                                                          dIncrementSlip_dStress2PK(:,:, iSlipSystem),   &
                                                          dStress2PkThermal_dF )
                                                          
            dDeformationGradPlasticInv_dF(:,:,:,:) = dDeformationGradPlasticInv_dF(:,:,:,:) +            &
                                                     computeProductDyadic2d2d(tensorFactorQ,             &
                                                                              dIncrementSlip_dF)
                                                      
         end do

      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of 2nd PK stress w.r.t. F   -----------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         call transform1to2(stressSecondPkThermalVoigtNew, stressSecondPkThermalNew)
         !deformationGradPlasticNewInv = computeMatrixInverse(deformationGradientPlasticNew)
         deformationGradPlasticNewInv = deformationGradientPlasticNew
         !call la_matinv(deformationGradPlasticNewInv,3,3,iflaginv)
         call lapack_invert(3,deformationGradPlasticNewInv,iflaginv)
         
         dStress2Pk_dF(:,:,:,:) = ZERO
         tensorFactorR = matmul(deformationGradThermalNewInv, matmul(stressSecondPkThermalNew,           &
                         matmul(transpose(deformationGradThermalNewInv),                                 &
                                transpose(deformationGradPlasticNewInv))))
         tensorFactorO = computeProductTensorLower(matrixIdentity3x3, transpose(tensorFactorR))
         dStress2Pk_dF = contractTensorsDouble4d4d(tensorFactorO, dDeformationGradPlasticInv_dF)
         
         tensorFactorQ = matmul(deformationGradPlasticNewInv, deformationGradThermalNewInv)
         tensorFactorO = computeProductTensorLower(tensorFactorQ, tensorFactorQ)
         dStress2Pk_dF = dStress2Pk_dF + contractTensorsDouble4d4d(tensorFactorO, dStress2PkThermal_dF)
         
         tensorFactorO = computeProductTensorUpper(transpose(tensorFactorR), matrixIdentity3x3)
         dStress2Pk_dF = dStress2Pk_dF + contractTensorsDouble4d4d(tensorFactorO,                        &
                                                                            dDeformationGradPlasticInv_dF)
                                                                           
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of 1st PK stress w.r.t. F   -----------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         Stress2Pk = computeDeterminant3x3(deformationGradientThermalNew) *                                 &
                     matmul(tensorFactorQ, matmul(stressSecondPkThermalNew, transpose(tensorFactorQ))) 
                     
         tensorFactorO = computeProductTensorLower(matrixIdentity3x3, Stress2Pk)
         tensorElasticityFirst = contractTensorsDouble4d4d(tensorFactorO, tensorIndentity)
         
         tensorFactorO = computeProductTensorUpper(deformationGradientNew, matrixIdentity3x3)
         tensorElasticityFirst(:,:,:,:) = tensorElasticityFirst(:,:,:,:) +                               &
                                          contractTensorsDouble4d4d(tensorFactorO, dStress2Pk_dF)
         
      end
      
      subroutine jacobian_new(dfgrd0,dfgrd1,fp_t,fp_tau,s_tan, & 
       c_alpha_1d,tau, & 
       delta_gamma,dGamma_dTau, & 
       dtime,de_g,tpk_1d_tau, & 
       ddsdde,ntens,islip) 
      
      implicit real*8(a-h,o-z)
      
      parameter (nslip=100)
      dimension de_g(3,3,3,3),dGamma_dTau(nslip),tpk_2d_g(3,3), & 
           ckrone(3,3,3,3),c_alpha(3,3,nslip),g_alpha_tau(nslip), & 
           dfgrd0(3,3),ddsdde(6,6),tau(nslip),dfgrd1(3,3)
     

      dimension c_alpha_1d(6,nslip),fe_t(3,3),fe_tau(3,3), & 
           fe_tau_t(3,3),tot_gamma(nslip), & 
           s_tan(3,3,nslip), & 
           fp_t(3,3),fp_t_inv(3,3), & 
           delta_gamma(nslip),fp_tau(3,3),fp_tau_inv(3,3), & 
           gamma_s_2d(3,3),ckrone_2d(3,3), & 
           tpk_1d_tau(6), & 
           dfgrdt(3,3),fp_tau_t(3,3)   


      
      real(8):: sqr_1(3,3,islip)
      real(8)::sqr_2(3,3,3,3),sqr_tmp(3,3,3,3)
      real(8)::sqr_3(3,3,3,3,islip),sqr_tmp1(3,3,3,3)
      real(8)::sqr_tmp2(3,3,3,3)
      real(8)::sqr_3_tmp(3,3),sqr_3_tmp_t(3,3)
      real(8):: fp_tau_inv_t(3,3)
      real(8)::cir_tmp1(3,3),cir_tmp2(3,3)
      real(8):: cir_1(3,3)
        real(8):: cir_2_tmp0(3,3,3,3),cir_2_tmp2_2d_inv(9,9)
        real(8):: cir_2_tmp1(3,3,3,3),cir_2_tmp2(3,3,3,3)
        real(8):: cir_2_tmp2_inv(3,3,3,3),cir_2_tmp3(3,3,3,3)
        real(8):: cir_2_tmp4(3,3,3,3),cir_2_tmp4_2d(9,9)
        real(8):: cir_2(3,3,3,3), cir_2_2d(9,9)
      real(8):: cir_3_tmp0(3,3),cir_3_tmp1(3,3)
      real(8):: cir_3(3,3,3,3)
      
      real(8):: final_0(3,3,3,3), final_0_2d_inv(9,9)
      real(8):: final_1(3,3,3,3),final_2(3,3,3,3)
      real(8):: final_3(3,3,3,3),final_4(3,3,3,3)
      real(8):: final_3_tmp(3,3),final_2_tmp(3,3)
      real(8):: final_5(3,3,3,3),final_6(3,3,3,3)
      real(8):: final_7(3,3,3,3),final_7_2d(9,9)
      real(8):: final_4d(3,3,3,3),final_2d_tmp(6,6)
      real(8):: final_2d(9,9),ddsdde_SE(6,6)
      real(8):: TS_ddsdde(3,3,3,3),dfgrd1_inv(3,3)
      real(8):: ddsdde_TK(6,6)
      
      integer:: a,b,m,n,i,j,k,l
      integer :: INFO
      


!__________________________________________________________
! initialization
!__________________________________________________________

      sqr_2=0.0
      sqr_3_tmp=0.0
      sqr_3=0.0
      cir_2_tmp1=0.0
      cir_2_tmp3=0.0
      cir_2_2d=0.0
      cir_3_tmp0=0.0
      cir_3_tmp1=0.0
      cir_3=0.0
      final_2_tmp=0.0
      final_3_tmp=0.0
      final_6=0.0
      final_2d=0.0



      do isys=1,islip

        do i=1,3
          do j=1,3
            sqr_1(i,j,isys)=0.5d0*dGamma_dTau(isys) & 
              *(s_tan(i,j,isys)+s_tan(j,i,isys))
          enddo
        enddo
      enddo
      
      
      
      
 !___________________________________________________________  
 !square_2
 !___________________________________________________________  

      
      
      call matinv3(fp_tau,fp_tau_inv,det_fp)
      do i=1,3
        do j=1,3
            do k=1,3
                do l=1,3
      sqr_tmp(i,j,k,l)=fp_tau_inv(i,k)*fp_tau_inv(j,l)
                enddo
            enddo
         enddo
      enddo
      
      
      do i=1,3
        do j=1,3
            do k=1,3
                do l=1,3  
                    do m=1,3
                        do n=1,3  
      sqr_2(i,j,k,l)=sqr_2(i,j,k,l)+de_g(i,j,m,n)*sqr_tmp(m,n,k,l)   
                        enddo
                    enddo
                enddo
            enddo
         enddo
      enddo
      
      
      
      
!___________________________________________________________       
!sqr_3
!___________________________________________________________ 

    
      
      do isys=1,islip   

      sqr_3_tmp=0.0
       do i=1,3
        do j=1,3
            do k=1,3
      
        sqr_3_tmp(i,j)=sqr_3_tmp(i,j)+fp_tau_inv(i,k)*s_tan(k,j,isys)
      
            enddo
        enddo
      enddo
      
      sqr_3_tmp_t = transpose(sqr_3_tmp)
      fp_tau_inv_t = transpose(fp_tau_inv)
      
       do i=1,3
        do j=1,3
            do k=1,3
                do l=1,3
      sqr_tmp1(i,j,k,l)=fp_tau_inv_t(i,k)*sqr_3_tmp_t(j,l)
      sqr_tmp2(i,j,k,l)=sqr_3_tmp_t(i,k)*fp_tau_inv_t(j,l)
      
                enddo
            enddo
         enddo
      enddo
      
      
      
       do i=1,3
        do j=1,3
            do k=1,3
                do l=1,3  
                    do m=1,3
                        do n=1,3 
      sqr_3(i,j,k,l,isys)=sqr_3(i,j,k,l,isys)+ & 
       de_g(i,j,m,n)*(sqr_tmp1(m,n,k,l)+sqr_tmp2(m,n,k,l)) 
                        enddo
                    enddo
                enddo
            enddo
         enddo
        enddo
      
      
      
      
      enddo
      
      
      
!___________________________________________________________        
! cir_1      
!___________________________________________________________       

      
      call sigmat(tpk_1d_tau,tpk_2d_g)
      fp_tau_t = transpose(fp_tau)
      
        call mat33(cir_tmp1,tpk_2d_g,fp_tau_t,3)
        call mat33(cir_tmp2,fp_tau,cir_tmp1,3)
      
        do i=1,3
        do j=1,3
        cir_1(i,j)=1.0/det_fp*cir_tmp2(i,j)
        enddo
        enddo
        
        
        
      
!___________________________________________________________        
! cir_2
!___________________________________________________________  

      
      
      call cal_ckrone_2d(ckrone_2d)  
      
      do i=1,3
         do j=1,3
            do k=1,3
                do l=1,3 
      cir_2_tmp0(i,j,k,l)=ckrone_2d(i,k)*ckrone_2d(j,l)
                enddo
            enddo
          enddo
        enddo

      call sigmat_crys(c_alpha_1d,c_alpha,islip)


      do isys=1,islip
        
       do i=1,3
         do j=1,3
            do k=1,3
                do l=1,3   
        cir_2_tmp1(i,j,k,l)=cir_2_tmp1(i,j,k,l)+ & 
       c_alpha(i,j,isys)*sqr_1(k,l,isys)
                enddo
             enddo
         enddo
       enddo
       
       enddo
       
       
     
       do i=1,3
         do j=1,3
            do k=1,3
                do l=1,3    
                
        cir_2_tmp2(i,j,k,l)= cir_2_tmp0(i,j,k,l)+cir_2_tmp1(i,j,k,l)   
        
                enddo
             enddo
         enddo
       enddo
      
       call tr4dto2d(cir_2_tmp2,cir_2_tmp2_2d_inv) 
     
       !call la_matinv(cir_2_tmp2_2d_inv,9,9,INFO)
       call lapack_invert(9,cir_2_tmp2_2d_inv,INFO)
       
      do isys=1,islip

       do i=1,3
         do j=1,3
            do k=1,3
                do l=1,3  
           cir_2_tmp3(i,j,k,l)=cir_2_tmp3(i,j,k,l)+ & 
       delta_gamma(isys)*sqr_3(i,j,k,l,isys)    
                enddo
             enddo
         enddo
       enddo
       
      enddo 
       
       do i=1,3
         do j=1,3
            do k=1,3
                do l=1,3
       cir_2_tmp4(i,j,k,l)=sqr_2(i,j,k,l)-cir_2_tmp3(i,j,k,l)
                enddo
             enddo
         enddo
       enddo


      call tr4dto2d(cir_2_tmp4,cir_2_tmp4_2d)
      
      do i=1,9
      do j=1,9
      do k=1,9
      
      cir_2_2d(i,j)=cir_2_2d(i,j)+ & 
        cir_2_tmp2_2d_inv(i,k)*cir_2_tmp4_2d(k,j)
     
      enddo
      enddo
      enddo

      call tr2dto4d(cir_2_2d, cir_2)
      
      
      

!___________________________________________________________       
!cir_3
!___________________________________________________________  

      
      
      do isys=1,islip

      cir_3_tmp0=0.0
      do i=1,3
      do j=1,3
      do k=1,3
      cir_3_tmp0(i,j)=cir_3_tmp0(i,j)+s_tan(i,k,isys)*fp_t(k,j)
      enddo
      enddo
      enddo
      
      cir_3_tmp1=0.0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      cir_3_tmp1(i,j)=cir_3_tmp1(i,j)+ &
      sqr_1(k,l,isys)*cir_2(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      cir_3(i,j,k,l)=cir_3(i,j,k,l)+ &
      cir_3_tmp0(i,j)*cir_3_tmp1(k,l)
      enddo  
      enddo
      enddo
      enddo
      
      enddo
      
      
!___________________________________________________________      
! final jacobian
!___________________________________________________________



      
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      
      final_0(i,j,k,l)=fp_tau(i,k)*fp_tau(j,l)
      
      enddo  
      enddo
      enddo
      enddo   
      
      
      call tr4dto2d(final_0, final_0_2d_inv)
      !call la_matinv(final_0_2d_inv,9,9,INFO)
      call lapack_invert(9,final_0_2d_inv,INFO)
      
      
      
      do i=1,9
      do j=1,9 
      
       final_0_2d_inv(i,j)=det_fp*final_0_2d_inv(i,j)     
      enddo  
      enddo
      
      
      
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
       final_1(i,j,k,l)=cir_1(i,j)*fp_tau_inv_t(k,l)       
      enddo  
      enddo
      enddo
      enddo 
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      final_2_tmp(j,i)=final_2_tmp(j,i)+ &
      tpk_2d_g(i,k)*fp_tau_t(k,j) 
      enddo  
      enddo
      enddo
      
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      final_2(i,j,k,l)=ckrone_2d(i,k)*final_2_tmp(j,l)
      enddo  
      enddo
      enddo
      enddo
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      final_3_tmp(i,j)=final_3_tmp(i,j)+ &
      fp_tau(i,k)*tpk_2d_g(k,j) 
      enddo  
      enddo
      enddo
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      final_3(i,j,k,l)=final_3_tmp(i,l)*ckrone_2d(j,k)
      enddo  
      enddo
      enddo
      enddo
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3

      final_4(i,j,k,l)=1.0/det_fp*(final_2(i,j,k,l)+final_3(i,j,k,l))
     
      enddo
      enddo
      enddo
      enddo
      
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      final_5(i,j,k,l)=final_1(i,j,k,l)-final_4(i,j,k,l)
      enddo  
      enddo
      enddo
      enddo
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
      final_6(i,j,k,l)=final_6(i,j,k,l)+final_5(i,j,m,n)*cir_3(m,n,k,l)
      enddo
      enddo
      enddo  
      enddo
      enddo
      enddo      
      
      
      
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      final_7(i,j,k,l)=cir_2(i,j,k,l)+final_6(i,j,k,l)
      enddo  
      enddo
      enddo
      enddo
      
      
      call tr4dto2d(final_7,final_7_2d)
      
      do i=1,9
      do j=1,9
      do k=1,9
      final_2d(i,j)=final_2d(i,j)+final_0_2d_inv(i,k)* & 
       final_7_2d(k,j) 
      enddo
      enddo
      enddo




      call tr2dto4d(final_2d,final_4d)
      ! final_4d is for:  T_ij=final_4d_ijkl*E_kl
      
      call tr4to2(final_4d,ddsdde)
       

    
      return
      end
      
      
      
      
      
      
      
      
      
      
      subroutine cvrt(a,b,c)
      
      implicit none
      integer:: a,b,c
      
      if (a==1) then
      b=1
      c=1
      elseif (a==2) then
      b=2
      c=2
      elseif (a==3) then
      b=3
      c=3
      elseif (a==4) then
      b=1
      c=2
      elseif (a==5) then
      b=1
      c=3
      elseif (a==6) then
      b=2
      c=3
      elseif (a==7) then
      b=2
      c=1
      elseif (a==8) then
      b=3
      c=1
      elseif (a==9) then
      b=3
      c=2
      endif
      
      return 
      end
      
      
      subroutine cvrt_inv(b,c,a)
      
      implicit none
      integer:: a,b,c
      
      if ((b==1).and.(c==1)) then
      a=1
      elseif ((b==2).and.(c==2)) then
      a=2
      elseif ((b==3).and.(c==3)) then
      a=3
      elseif ((b==1).and.(c==2)) then
      a=4
      elseif ((b==1).and.(c==3)) then
      a=5
      elseif ((b==2).and.(c==3))then
      a=6
      elseif ((b==2).and.(c==1)) then
      a=7
      elseif ((b==3).and.(c==1)) then
      a=8
      elseif ((b==3).and.(c==2)) then
      a=9
      endif
      
      return 
      end     
      
      
      
      
      subroutine tr4dto2d(A_4d, A_2d)
      
      implicit none
      real(8):: A_4d(3,3,3,3), A_2d(9,9)
      integer:: i,j,k,l,m,n
      
      
      do i=1,3
      do j=1,3
      call cvrt_inv(i,j,m)
      
      do k=1,3
      do l=1,3
      call cvrt_inv(k,l,n)
      
      A_2d(m,n)=A_4d(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
      
      end
      
      
      
      
      
      
      subroutine tr2dto4d(A_2d, A_4d)
      
      implicit none
      real(8):: A_4d(3,3,3,3), A_2d(9,9)
      integer:: i,j,k,l,m,n
      
      
      do m=1,9
      call cvrt(m,i,j)
      do n=1,9
      call cvrt(n,k,l)
      
      A_4d(i,j,k,l)=A_2d(m,n)
      enddo
      enddo
      
      end 
            
            
      ! ************************************************************************************************ !
      !
      !  Description:   Perform double contraction of inner indices on two 4d tensors
      !
      !       Inputs:   tensorA(:,:,:,:),          type = real(wp)
      !                 tensorB(:,:,:,:),          type = real(wp)
      !
      !       Result:   tensorContracted(:,:,:,:), type = real(wp)
      ! 
      !      Details:   C_ijkl = A_ijpq * B_pqkl
      !
      !   References:   
      !
      !     Creation:   Coleman Alleman | 2015/07/10 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function contractTensorsDouble4d4d     &
      (                                      &
         tensorA,                            &
         tensorB                             &
      )                                      &
         result(tensorContracted)


         implicit none


      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
            
         real(8), parameter :: &
            ZERO = 0.d0

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================

      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic :: size


      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         real(8), intent(in) ::     &
            tensorA(3,3,3,3),                &
            tensorB(3,3,3,3)

      !  Result ------------------------------------------------------------------------------------------

         real(8) ::                                                                          &
            tensorContracted(size(tensorA,1), size(tensorA,2), size(tensorB,3), size(tensorB,4))


      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Integers ----------------------------------------------------------------------------------------

         integer ::     &
            ii,                     &
            jj,                     &
            kk,                     &
            ll,                     &
            pp,                     &
            qq



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      ! ==================================================================================================
      !  Error Checks ------------------------------------------------------------------------------------
      ! ==================================================================================================

         if ((size(tensorA,3) /= size(tensorB,1)) .or. (size(tensorA,4) /= size(tensorB,2))) then
            stop "Error - array shapes do not match"
         end if


      ! ==================================================================================================
      !  Initializations ---------------------------------------------------------------------------------
      ! ==================================================================================================

         tensorContracted = ZERO


      ! ==================================================================================================
      !  Computation -------------------------------------------------------------------------------------
      ! ==================================================================================================

         do ii = 1, size(tensorA, 1)
            do jj = 1, size(tensorA, 2)
               do kk = 1, size(tensorB, 3)
                  do ll = 1, size(tensorB, 4)
                     do pp = 1, size(tensorA, 3)
                        do qq = 1, size(tensorA, 4)
                           tensorContracted(ii, jj, kk, ll) =     &
                                 tensorContracted(ii, jj, kk, ll)     &
                                 + tensorA(ii, jj, pp, qq) * tensorB(pp, qq, kk, ll)
                        end do
                     end do
                  end do
               end do
            end do
         end do

      end function contractTensorsDouble4d4d
            
      ! ************************************************************************************************ !
      !
      !  Description:   Perform double contraction of inner indices on a 2d tensor with a 2d tensor
      !
      !       Inputs:   tensorA(:,:,),          type = real(wp)
      !                 tensorB(:,:),           type = real(wp)
      !
      !       Result:   tensorContracted(:,:),  type = real(wp)
      ! 
      !      Details:   c = sum(A_ij * B_ij)
      !
      !   References:   
      !
      !     Creation:   Coleman Alleman | 2015/07/10 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function contractTensorsDouble2d2d     &
      (                                      &
         tensorA,                            &
         tensorB                             &
      )                                      &
         result(scalarContracted)



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================

         implicit none


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic :: size


      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         real(8), intent(in) ::     &
            tensorA(3,3),                    &
            tensorB(3,3)

      !  Result ------------------------------------------------------------------------------------------

         real(8) ::     &
            scalarContracted

         real(8) :: tensorMult(3,3)
         integer :: i

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! ==================================================================================================
      !  Error Checks ------------------------------------------------------------------------------------
      ! ==================================================================================================

         if ((size(tensorA,1) /= size(tensorB,1)) .or. (size(tensorA,2) /= size(tensorB,2))) then
            stop "Error - array shapes do not match"
         end if


      ! ==================================================================================================
      !  Computation -------------------------------------------------------------------------------------
      ! ==================================================================================================
         tensorMult = matmul(transpose(tensorA), tensorB)

         scalarContracted = 0.d0
         do i = 1,size(tensorMult,1)
            scalarContracted = scalarContracted + tensorMult(i,i)
         enddo


      end function contractTensorsDouble2d2d
                  
      ! ************************************************************************************************ !
      !
      !  Description:   Perform double contraction of inner indices on a 2d tensor with a 4d tensor
      !
      !       Inputs:   tensorA(:,:),       type = real(wp)
      !                 tensorB(:,:,:,:),   type = real(wp)
      !
      !       Result:   tensorContracted(:,:),  type = real(wp)
      ! 
      !      Details:   C_ij = A_kl * B_klij
      !
      !   References:   
      !
      !     Creation:   Ahmad Shahba | 2015/2/4 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function contractTensorsDouble2d4d     &
      (                                      &
         tensorA,                            &
         tensorB                             &
      )                                      &
         result(tensorContracted)



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         implicit none
         
      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================

         real(8), parameter :: &
            ZERO = 0.d0

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic :: size


      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         real(8), intent(in) ::     &
            tensorA(3,3),                    &
            tensorB(3,3,3,3)

      !  Result ------------------------------------------------------------------------------------------

         real(8) ::                                      &
            tensorContracted(size(tensorB,3), size(tensorB,4))


      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Integers ----------------------------------------------------------------------------------------

         integer ::     &
            ii,                     &
            jj,                     &
            kk,                     &
            ll



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      ! ==================================================================================================
      !  Error Checks ------------------------------------------------------------------------------------
      ! ==================================================================================================

         if ((size(tensorA,1) /= size(tensorB,1)) .or. (size(tensorA,2) /= size(tensorB,2))) then
            stop "Error - array shapes do not match"
         end if


      ! ==================================================================================================
      !  Initializations ---------------------------------------------------------------------------------
      ! ==================================================================================================

         tensorContracted = ZERO


      ! ==================================================================================================
      !  Computation -------------------------------------------------------------------------------------
      ! ==================================================================================================

         do ii = 1, size(tensorB, 3)
            do jj = 1, size(tensorB, 4)
               do kk = 1, size(tensorB, 1)
                  do ll = 1, size(tensorB, 2)
                     tensorContracted(ii, jj) = tensorContracted(ii, jj)      &
                                              + tensorA(kk, ll) * tensorB(kk, ll, ii, jj)
                  end do
               end do
            end do
         end do

      end function contractTensorsDouble2d4d
            
      ! ************************************************************************************************ !
      !
      !  Description:   Perform dyadic product of two 2d tensors
      !
      !       Inputs:   tensorA(:,:),           type = real(wp)
      !                 tensorB(:,:),           type = real(wp)
      !
      !       Result:   tensorContracted(:,:,:,:),  type = real(wp)
      ! 
      !      Details:   C_ijkl = A_ij * B_kl
      !
      !   References:   
      !
      !     Creation:   Ahmad Shahba | 2015/02/04 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function computeProductDyadic2d2d            &
      (                                            &
         tensorA,                                  &
         tensorB                                   &
      )                                            &
         result(tensorProductDyadic)



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         implicit none
      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================

         real(8), parameter :: &
            ZERO = 0.d0

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic :: size


      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         real(8), intent(in) ::     &
            tensorA(3,3),                    &
            tensorB(3,3)

      !  Result ------------------------------------------------------------------------------------------

         real(8) ::                                      &
            tensorProductDyadic(size(tensorA,1), size(tensorA,2),size(tensorB,1), size(tensorB,2))


      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Integers ----------------------------------------------------------------------------------------

         integer ::     &
            ii,                     &
            jj,                     &
            kk,                     &
            ll



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



      ! ==================================================================================================
      !  Initializations ---------------------------------------------------------------------------------
      ! ==================================================================================================

         tensorProductDyadic = ZERO


      ! ==================================================================================================
      !  Computation -------------------------------------------------------------------------------------
      ! ==================================================================================================

         do ii = 1, size(tensorA, 1)
            do jj = 1, size(tensorA, 2)
               do kk = 1, size(tensorB, 1)
                  do ll = 1, size(tensorB, 2)
                     tensorProductDyadic(ii, jj, kk, ll) = tensorA(ii, jj) * tensorB(kk, ll)
                  end do
               end do
            end do
         end do

      end function computeProductDyadic2d2d
      
      !
      ! ************************************************************************************************ !
      !
      !  Description:   This function computes the lower dyadic product of two 2nd order tensors 
      !
      !       Inputs:       
      !
      !                    
      !      Outputs:        
      ! 
      !      Details:   
      !
      !   References:   
      !
      !     Creation:   Ahmad Shahba  | 2016/2/1 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function computeProductTensorLower                    &
      (                                                     &
         A,                                                 &
         B                                                  &
      )                                                     &
      result(productTensorLower)




      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         implicit none
         
      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
             
      real(8), parameter :: &
            ZERO = 0.d0            

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------


      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
         
      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

            

      !  Inputs ------------------------------------------------------------------------------------------
         
         real(8), intent(in) ::                       &
            A(3, 3),                                           &
            B(3, 3)


      !  Outputs -----------------------------------------------------------------------------------------

         real(8)  ::                                  &
            productTensorLower(3, 3, 3, 3)


      ! --------------------------------------------------------------------------------------------------
      !  Local Constants ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         integer ::                                &
            ii,                                                &
            jj,                                                &
            kk,                                                &
            ll

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         productTensorLower = ZERO
         
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     productTensorLower(ii, jj, kk, ll) = A(ii, kk) * B(jj, ll)
                  end do
               end do
            end do
         end do

      end function computeProductTensorLower


      !
      ! ************************************************************************************************ !
      !
      !  Description:   This function computes the upper dyadic product of two 2nd order tensors 
      !
      !       Inputs:       
      !
      !                    
      !      Outputs:        
      ! 
      !      Details:   
      !
      !   References:   
      !
      !     Creation:   Ahmad Shahba  | 2016/2/1 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function computeProductTensorUpper                    &
      (                                                     &
         A,                                                 &
         B                                                  &
      )                                                     &
      result(productTensorUpper)




      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         implicit none

      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
             
      real(8), parameter :: &
            ZERO = 0.d0
      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------


      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
         
      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

            

      !  Inputs ------------------------------------------------------------------------------------------
         
         real(8), intent(in) ::                       &
            A(3, 3),                                           &
            B(3, 3)


      !  Outputs -----------------------------------------------------------------------------------------

         real(8)  ::                                  &
            productTensorUpper(3, 3, 3, 3)


      ! --------------------------------------------------------------------------------------------------
      !  Local Constants ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         integer ::                                &
            ii,                                                &
            jj,                                                &
            kk,                                                &
            ll

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         productTensorUpper = ZERO
         
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     productTensorUpper(ii, jj, kk, ll) = A(ii, ll) * B(jj, kk)
                  end do
               end do
            end do
         end do

      end function computeProductTensorUpper
            
      ! ************************************************************************************************ !
      !
      !  Description:   This subroutine computes the elemental tangent stiffness matrix
      !
      !       Inputs:
      !
      !      Outputs:
      !
      !      Details:
      !
      !   References:
      !
      !     Creation:   Ahmad Shahba | 2016/04/19 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      subroutine computeMatrixStiffnessTangent     &
      (                                            &
         stressCauchyNewVoigt,                     &
         tensorFirstElasticity,                    &
         deformationGradNew,                       &
         repositoryOperatorGradientNew,            &
         repositoryDetJacobianNew,                 &
         iElement,                                 &
         TotalElmsInPatch,                         &
         matrixStiffnessTangent,                   &
         matrixStiffnessTangentNonlocal            &
      )


      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         implicit none

      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================

         real(8), parameter :: &
            ZERO = 0.d0,                            &
            ONE = 1.d0,                             &
            TWO = 2.d0, &
            HALF = 0.5d0,                            &
            THIRD = 1.d0/3.d0
            
      real(8), parameter :: weightGaussPoint = 1.d0

      integer, parameter :: nNodesPerElement = 4, &
                            nDOFPerNode = 3, &
                            nDimensions = 3
      character(len=30), parameter :: typeElement = 'TET4F'

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================



      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic ::      &
            matmul,        &
            transpose,     &
            sum

      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
                  include 'jacobian_fbar.inc'

      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         integer, intent(in) ::           &
            iElement,                                 &
            TotalElmsInPatch

         real(8), intent(in) ::              &
            stressCauchyNewVoigt(6),                  &
            deformationGradNew(3,3),                  &
            tensorFirstElasticity(3,3,3,3),           &
            repositoryOperatorGradientNew(3,4,TotalElmsInPatch),   &
            repositoryDetJacobianNew(TotalElmsInPatch)

      !  Outputs -----------------------------------------------------------------------------------------

         real(8), intent(inout) ::        &
            matrixStiffnessTangent(12,12),           &
            matrixStiffnessTangentNonlocal(12,12,TotalElmsInPatch)


      ! --------------------------------------------------------------------------------------------------
      !  Local Constants ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
         real(8) :: matrixIdentity3x3(3,3)            
         
         integer ::                    &
            ii,                                    &
            jj,                                    &
            kk,                                    &
            mm,                                    &
            nn,                                    &
            pp,                                    &
            indexNode,                             &
            jndexNode,                             &
            iDof,                                  &
            jDof,                                  &
            numberRow,                             &
            numberColumn

         real(8) ::                       &
            detDeformationGradNew,                 &
            integrationConstant


         real(8) ::                             &
            stressCauchyNew(3, 3),                       &
            tensorFirstSpatialElasticity(3, 3, 3, 3),    &
            factorP(3, 3, 3, 3),                         &
            factorQ(3, 3, 3, 3),                         &
            factorR(3, 3, 3, 3),                         &
            factorS(3,1),                                &
            factorT(3,1),                                &
            matrixElasticity3x3(3,3),                    &
            factorU(3, 3),                               &
            factorV(1, 1)

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
         matrixStiffnessTangent = 0.d0
         matrixStiffnessTangentNonlocal = 0.d0
         
         call cal_ckrone_2d(matrixIdentity3x3)  
         
         call transform1To2(stressCauchyNewVoigt, stressCauchyNew)

         if ( typeElement == 'TET4F' ) then

            detDeformationGradNew = computeDeterminant3x3(deformationGradNew) ! this is the 'fbar-modified def grad'
            integrationConstant = HALF * THIRD

         else

               write(*, *) 'This type of element is not implemented', typeElement
               stop

         end if


         tensorFirstSpatialElasticity(:,:,:,:) = ZERO
         do ii = 1,3
            do jj = 1,3
               do kk = 1,3
                  do mm = 1,3
                     do nn = 1,3
                        do pp = 1,3
                           tensorFirstSpatialElasticity(ii, jj, kk, mm) =                             &
                              tensorFirstSpatialElasticity(ii, jj, kk, mm) + ONE /                    &    
                              detDeformationGradNew * deformationGradNew(jj, nn) *                    &
                              deformationGradNew(mm, pp) * tensorFirstElasticity(ii, nn, kk, pp)
                        end do
                     end do
                  end do
               end do
            end do
         end do

         factorP(:,:,:,:) = ZERO
         factorQ(:,:,:,:) = ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do mm = 1, 3

                     do nn = 1, 3
                        factorP(ii, jj, kk, mm) = factorP(ii, jj, kk, mm) +                           &
                              tensorFirstSpatialElasticity(ii, jj, nn, nn) * matrixIdentity3x3(kk, mm)
                     end do

                     factorQ(ii, jj, kk, mm) = stressCauchyNew(ii, jj) * matrixIdentity3x3(kk, mm)
                     
                  end do
               end do
            end do
         end do

         factorR = THIRD * (factorP - TWO * factorQ)      
         
         if ( typeElement == 'TET4F' ) then

            do indexNode = 1, nNodesPerElement
               factorS(1 : 3, 1) = repositoryOperatorGradientNew(1 : 3, indexNode, iElement)
               do iDof = 1, nDOFPerNode
               
                  do jndexNode = 1, nNodesPerElement
                     factorT(1 : 3, 1) = repositoryOperatorGradientNew(1:3, jndexNode, iElement)
                     do jDof = 1, nDOFPerNode
                                
                        do ii = 1, 3
                           do jj = 1, 3
                              matrixElasticity3x3(ii, jj) = tensorFirstSpatialElasticity(iDof, ii, jDof, jj)
                              factorU(ii,jj) = factorR(iDof, ii, jDof, jj)
                           end do
                        end do
                        
                        factorV =                                                                        &
                           matmul(matmul(transpose(factorS), matrixElasticity3x3), factorT) *            &
                           integrationConstant * repositoryDetJacobianNew(iElement) *        &
                           weightGaussPoint + matmul(matmul(transpose(factorS), factorU), factorT)*      &
                           integrationConstant * repositoryDetJacobianNew(iElement) *       &
                           weightGaussPoint * (repositoryDetJacobianNew(iElement) /         &
                           sum(repositoryDetJacobianNew(:)) - ONE)
              
                        matrixStiffnessTangent((indexNode -1) * nDOFPerNode + iDof,             &
                                               (jndexNode-1) * nDOFPerNode + jDof) =            &
                           matrixStiffnessTangent((indexNode -1) * nDOFPerNode + iDof,          &
                                                  (jndexNode-1) * nDOFPerNode + jDof) + factorV(1, 1)

                     end do
                  end do
                  
                  
                  do ii = 1, TotalElmsInPatch
                      
                     if (ii /= iElement) then

                        do jndexNode = 1, nNodesPerElement
                              
                           factorT(1:3, 1)= repositoryOperatorGradientNew(1:3, jndexNode, ii)
                          
                           do jDof = 1, nDOFPerNode

                              do jj = 1, 3
                                 do kk = 1, 3
                                    factorU(jj, kk) = factorR(iDof, jj, jDof, kk)
                                 end do
                              end do
                      
                              factorV =                                                                  &
                                 matmul(matmul(transpose(factorS), factorU), factorT) *                  &
                                 integrationConstant * repositoryDetJacobianNew(iElement) * &
                                 weightGaussPoint * (repositoryDetJacobianNew(ii) /         &
                                 sum(repositoryDetJacobianNew(:)))
        
                              matrixStiffnessTangentNonlocal((indexNode-1) * nDOFPerNode + iDof,         &
                                 (jndexNode-1) * nDOFPerNode + jDof, ii) =                               &
                                    matrixStiffnessTangentNonlocal((indexNode-1) * nDOFPerNode + iDof,   &
                                    (jndexNode - 1) * nDOFPerNode + jDof, ii) + factorV(1, 1)
                                  
                           end do
                        end do
                     end if
                  end do            
               end do
            end do

         else

            write(*, *) 'This type of element is not implemented', typeElement
            stop         
         end if 

      end subroutine computeMatrixStiffnessTangent
      
      ! ************************************************************************************************ !
      !
      !  Description:    This subroutine computes the Material Jacobian defined as dP/dF where P and F
      !                 correspond to the 1st PK stress and total deformation gradient, respectively.
      !                 This subroutine should only be used for the thermo-elastic constitutive model 
      !
      !
      !       Inputs:
      !
      !      Outputs:
      !
      !      Details:
      !
      !   References:
      !
      !     Creation:   Ahmad Shahba | 2016/02/25 | Baltimore, MD
      !
      ! ************************************************************************************************ !
      subroutine computeTensorElasticityFirstThermElastCm         &
      (                                                           &
         tensorElasticityFirst,                                   &
         stressSecondPkThermal,                                   &
         deformationGradientNew,                                  &
         tensorElasticity,                                        &
         deformationGradientThermalNew,                           &
         deformationGradientElasticNew                            &
      )

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         implicit none

      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
         real(8), parameter :: &
            ZERO = 0.d0,                            &
            ONE = 1.d0,                             &
            TWO = 2.d0, &
            HALF = 0.5d0

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================



      ! --------------------------------------------------------------------------------------------------
      !  Intrinsic Functions -----------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         intrinsic ::      &
            matmul,        &
            transpose

      ! --------------------------------------------------------------------------------------------------
      !  Include Files -----------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
                  include 'jacobian_fbar.inc'

      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------
         real(8), intent(in) ::        &
            stressSecondPkThermal(3,3),         &
            deformationGradientNew(3,3),        &
            tensorElasticity(3,3,3,3),          &
            deformationGradientThermalNew(3,3), &
            deformationGradientElasticNew(3,3)

      !  Outputs -----------------------------------------------------------------------------------------
         real(8), intent(out) ::       &
            tensorElasticityFirst(3,3,3,3)

      ! --------------------------------------------------------------------------------------------------
      !  Local Constants ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------
         real(8) :: matrixIdentity3x3(3,3)            
         
         integer :: iflaginv
         
         integer ::                             &
            ii,                                             &
            jj,                                             &
            kk,                                             &
            ll

         real(8) ::                                &
            tensorIndentityTranspose(3, 3, 3, 3),           &
            tensorFactorP(3, 3, 3, 3),                      &
            dEe_dFe(3, 3, 3, 3),                            &
            deformationGradThermalNewInv(3, 3),             &
            tensorIndentity(3, 3, 3, 3),                    &
            dFe_dF(3, 3, 3, 3),                             &
            dStress2PkThermal_dEe(3, 3, 3, 3),              &
            dStress2PkThermal_dF(3, 3, 3, 3),               &
            dStress2Pk_dF(3, 3, 3, 3),                      &
            Stress2Pk(3, 3)


      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      !  =================================== Initialization =============================================
         tensorElasticityFirst(:,:,:,:) = ZERO
         
         call cal_ckrone_2d(matrixIdentity3x3)  
         
         !  ------------------------------------------------------------------------------------------------
      !  Calculate identity tensor  ------------------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         tensorIndentityTranspose(:,:,:,:) = ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     if ( ii == ll .and. jj == kk) tensorIndentityTranspose(ii, jj, kk, ll) = ONE
                  end do
               end do
            end do
         end do 
         
         tensorIndentity(:,:,:,:) = ZERO
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     if ( ii == kk .and. jj == ll) tensorIndentity(ii, jj, kk, ll) = ONE
                  end do
               end do
            end do
         end do


         
      !  Calculate derivative of Ee w.r.t. Fe
         tensorFactorP = computeProductTensorLower(matrixIdentity3x3,                                    &
                                                                     transpose (deformationGradientNew))+& 
                         computeProductTensorUpper(transpose (deformationGradientNew),                   &
                                                                                      matrixIdentity3x3)
         dEe_dFe = HALF * contractTensorsDouble4d4d(tensorFactorP, tensorIndentityTranspose)                                                                                
         
         
      !  Calculate derivative of Fe w.r.t. F  
         !deformationGradThermalNewInv = computeMatrixInverse(deformationGradientThermalNew)
         deformationGradThermalNewInv = deformationGradientThermalNew
         !call la_matinv(deformationGradThermalNewInv,3,3,iflaginv)
         call lapack_invert(3,deformationGradThermalNewInv,iflaginv)
         
         tensorFactorP = computeProductTensorLower(matrixIdentity3x3,                                    &
                                                               transpose(deformationGradThermalNewInv)) 
         dFe_dF = contractTensorsDouble4d4d(tensorFactorP, tensorIndentity)   


      !  Calculate derivative of 2nd PK in theramlly-expanded confguration w.r.t. Ee
         do ii = 1, 3
            do jj = 1, 3
               do kk = 1, 3
                  do ll = 1, 3
                     dStress2PkThermal_dEe(ii, jj, kk, ll) = HALF * (tensorElasticity(ii, jj, kk, ll) +  &
                                                                     tensorElasticity(ii, jj, ll, kk) )
                  end do
               end do
            end do
         end do

      !  Calculate derivative of 2nd PK in theramlly-expanded confguration w.r.t. F
         dStress2PkThermal_dF = contractTensorsDouble4d4d(contractTensorsDouble4d4d(                     &
                                                                  dStress2PkThermal_dEe, dEe_dFe), dFe_dF)   

      !  Calculate derivative of 2nd PK w.r.t. F
         dStress2Pk_dF = contractTensorsDouble4d4d(                                                      &
                                 computeProductTensorLower(deformationGradThermalNewInv,                 &
                                                         deformationGradThermalNewInv),                  &
                                 dStress2PkThermal_dF)

                                                            
      !  ------------------------------------------------------------------------------------------------
      !  Calculate derivative of 1st PK stress w.r.t. F   -----------------------------------------------
      !  ------------------------------------------------------------------------------------------------
         Stress2Pk = computeDeterminant3x3(deformationGradientThermalNew) *                                 &
                     matmul(deformationGradThermalNewInv, matmul(stressSecondPkThermal,                  &
                                                               transpose(deformationGradThermalNewInv))) 

                   
         tensorFactorP = computeProductTensorLower(matrixIdentity3x3, Stress2Pk)
         tensorElasticityFirst = contractTensorsDouble4d4d(tensorFactorP, tensorIndentity)
         
         tensorFactorP = computeProductTensorUpper(deformationGradientNew, matrixIdentity3x3)
         tensorElasticityFirst(:,:,:,:) = tensorElasticityFirst(:,:,:,:) +                               &
                                          contractTensorsDouble4d4d(tensorFactorP, dStress2Pk_dF)

      end subroutine computeTensorElasticityFirstThermElastCm
      !                                                                     
      !  arranges the six components of a vector into a symmetric 3x3 array
      !          
                                                                 
      subroutine transform1To2(x,xx)    

      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================
             

      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================

         implicit none

         real(8) :: x(6), xx(3,3)
         
         xx(1,1)=x(1)                                                    
         xx(2,2)=x(2)                                                    
         xx(3,3)=x(3)                                                    
         xx(1,2)=x(4)                                                    
         xx(2,1)=x(4)
         xx(3,1)=x(5)                                                    
         xx(1,3)=x(5)                                                    
         xx(2,3)=x(6)                                                    
         xx(3,2)=x(6)                                                    
                                                          
      end
      

      ! ************************************************************************************************ !
      !
      !  Description:   Compute the determinant of a 3x3 matrix
      !
      !       Inputs:   matrix(3,3), type = real(wp)
      !
      !       Result:   determinant, type = real(wp)
      ! 
      !      Details: 
      !
      !   References:   
      !
      !     Creation:   Coleman Alleman | 2015/06/11 | Baltimore, MD
      !
      ! ************************************************************************************************ !

      function computeDeterminant3x3      &
      (                                   &
         matrix                           &
      )                                   &
         result(determinant)



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ----------------------------------- Specification Section ----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      ! ==================================================================================================
      !  Modules -----------------------------------------------------------------------------------------
      ! ==================================================================================================


      ! ==================================================================================================
      !  Declarations ------------------------------------------------------------------------------------
      ! ==================================================================================================

         implicit none


      ! --------------------------------------------------------------------------------------------------
      !  Argument Variables ------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

      !  Inputs ------------------------------------------------------------------------------------------

         real(8), intent(in) ::     &
            matrix(3,3)

      !  Result ------------------------------------------------------------------------------------------

         real(8) ::     &
            determinant



      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --------------------------------------------------------------------------------------------------
      ! ------------------------------------- Executable Section -----------------------------------------
      ! --------------------------------------------------------------------------------------------------
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         determinant = matrix(1,1) * matrix(2,2) * matrix(3,3)    &
                     + matrix(1,2) * matrix(2,3) * matrix(3,1)    &
                     + matrix(1,3) * matrix(2,1) * matrix(3,2)    &
                     - matrix(3,1) * matrix(2,2) * matrix(1,3)    &
                     - matrix(3,2) * matrix(2,3) * matrix(1,1)    &
                     - matrix(3,3) * matrix(2,1) * matrix(1,2)


      end function computeDeterminant3x3
            

! in FEM_UEL
!
!      SUBROUTINE solveAXb(A,N,b,err)
!      integer, intent(in):: N
!      real(8), intent(in):: A(N,N)
!      real(8), intent(inout):: b(N)
!      integer, intent(out):: err
!      integer*4 :: err_4
!         CALL DPOSV('U',N,1,A,N,b,N,err_4)
!         err = err_4
!      END SUBROUTINE