module FatesParameterDerivedMod

  ! -------------------------------------------------------------------------------------
  ! This module contains all procedures types and settings for any quantities that are
  ! statically derived from static model parameters.  These are unchanging quantities
  ! and are based off of simple relationships from parameters that the user can
  ! vary.  This should be called once, and early in the model initialization call
  ! sequence immediately after FATES parameters are read in.
  !
  ! -------------------------------------------------------------------------------------

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesConstantsMod,     only : umolC_to_kgC
  use FatesConstantsMod,     only : g_per_kg
  use FatesInterfaceTypesMod,     only : nleafage
  use FatesInterfaceTypesMod,     only : ncrowndamage
  use FatesGlobals     ,     only : fates_log
  
  implicit none
  private

  type, public :: param_derived_type

     real(r8), allocatable :: jmax25top(:,:)  ! canopy top: maximum electron transport 
                                              ! rate at 25C (umol electrons/m**2/s)
     real(r8), allocatable :: tpu25top(:,:)   ! canopy top: triose phosphate utilization
                                              ! rate at 25C (umol CO2/m**2/s)
     real(r8), allocatable :: kp25top(:,:)    ! canopy top: initial slope of CO2 response
                                              ! curve (C4 plants) at 25C

     real(r8), allocatable :: branch_frac(:)  ! fraction of aboveground biomass in branches (as
                                              ! oppose to stems) - for use in damage allometries

     real(r8), allocatable :: damage_transitions(:,:,:) ! matrix of transition probabilities between
                                                      ! damage classes - one per PFT
     
   contains
     
     procedure :: Init
     procedure :: InitDamageTransitions
     procedure :: InitAllocate
     procedure :: InitAllocateDamageTransitions
     
  end type param_derived_type
  
  type(param_derived_type), public :: param_derived
  
contains

  real(r8) function integral_exp(x, y)

    real(r8), intent(in) :: x
    real(r8), intent(in) :: y

    integral_exp = (1/y) * exp(y*x)
    return
  end function integral_exp

  real(r8) function damage_integral(l, u, y)

    real(r8), intent(in) :: l
    real(r8), intent(in) :: u
    real(r8), intent(in) :: y

    real(r8) :: lower
    real(r8) :: upper

    lower = integral_exp(l, y)
    upper = integral_exp(u, y)
    damage_integral = upper - lower
    return
  end function damage_integral
  
  ! ===================================================================================
  subroutine InitAllocate(this,numpft)
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    allocate(this%jmax25top(numpft,nleafage))
    allocate(this%tpu25top(numpft,nleafage))
    allocate(this%kp25top(numpft,nleafage))

    allocate(this%branch_frac(numpft))
    
    
    return
  end subroutine InitAllocate

  ! =====================================================================================

 ! ===================================================================================
  subroutine InitAllocateDamageTransitions(this,ncrowndamage, numpft)
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: ncrowndamage
    integer, intent(in)                      :: numpft

    allocate(this%damage_transitions(ncrowndamage,ncrowndamage, numpft))
    
    return
  end subroutine InitAllocateDamageTransitions

  ! =====================================================================================
 
  subroutine Init(this,numpft)

    use EDPftvarcon, only: EDPftvarcon_inst
    use SFParamsMod, only: SF_val_CWD_frac
    use FatesLitterMod, only : ncwd
    
    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: numpft
    
    ! local variables
    integer  :: ft                 ! pft index
    integer  :: iage               ! leaf age class index
    integer  :: c                  ! cwd index

    associate( vcmax25top => EDPftvarcon_inst%vcmax25top ) 
    
      call this%InitAllocate(numpft)
      
      do ft = 1,numpft
         
         do iage = 1, nleafage

            ! Parameters derived from vcmax25top. 
            ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
            ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of 
            ! Experimental Botany 44:907-920.  Here use a factor "1.67", from 
            ! Medlyn et al (2002) Plant, Cell and Environment 25:1167-1179
            
            ! RF - copied this from the CLM trunk code, but where did it come from, 
            ! and how can we make these consistant? 
            ! jmax25top(ft) =  &
            ! (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrzc),11._r8),35._r8)) * vcmax25top(ft)
            
            this%jmax25top(ft,iage) = 1.67_r8   * vcmax25top(ft,iage)
            this%tpu25top(ft,iage)  = 0.167_r8  * vcmax25top(ft,iage)
            this%kp25top(ft,iage)   = 20000._r8 * vcmax25top(ft,iage)
         
         end do

         ! Allocate fraction of biomass in branches
         this%branch_frac(ft) = sum(SF_val_CWD_frac(1:3))
         
      end do !ft

    end associate
    return
  end subroutine Init

!=========================================================================
  
  subroutine InitDamageTransitions(this, ncrowndamage, numpft)

    use FatesConstantsMod,      only : days_per_year
    use EDPftvarcon, only: EDPftvarcon_inst


    class(param_derived_type), intent(inout) :: this
    integer, intent(in)                      :: ncrowndamage
    integer, intent(in)                      :: numpft

    ! local variables
    integer  :: ft                ! pft index
    integer  :: i                 ! crowndamage index (rows)
    integer  :: j                 ! crowndamage index (columns)
    real(r8) :: exponent          ! in the function for damage transitions 
    real(r8) :: cd_real           ! for in functions
    real(r8), allocatable :: transition_vec(:)

    allocate(transition_vec(ncrowndamage))

    call this%InitAllocateDamageTransitions(ncrowndamage, numpft)

    do ft = 1, numpft

        exponent = EDPftvarcon_inst%damage_exponent(ft)

       ! make a look up matrix that holds transition rates between
       ! damage classes. 
       ! Transition rates of moving from one class to the next are integral of
       ! a negative exponential with exponent a parameter from param file
!       do i = 1, ncrowndamage
!          cd_real = real(i)
!          transition_vec(i) = damage_integral(cd_real-1.0_r8, cd_real, exponent)
!       end do
!        normalise so they sum to 1
!       transition_vec = transition_vec/sum(transition_vec)

       ! ! populate a matrix
       ! do i = 1, ncrowndamage     ! new
       !    do j = 1, ncrowndamage  ! current
       !       if(j > i) then       ! can't move to less damaged
       !          this%damage_transitions(j,i,ft) = 0.0_r8
       !       else
       !          this%damage_transitions(j,i,ft) = transition_vec(i - j + 1)
       !       end if
       !    end do
       ! end do
       
       ! ! The above is for background small damage - it results in gradual build up of damage
       ! ! Below we account for catastrophic damage - which shifts trees to high damage classes
       ! ! regardless of starting damage level
       ! this%damage_transitions(1:ncrowndamage-2,ncrowndamage-1,ft) = 0.05/days_per_year
       ! this%damage_transitions(1:ncrowndamage-1,ncrowndamage,ft)   = 0.025/days_per_year
       
       ! do i = 1, ncrowndamage
       !    this%damage_transitions(i,:,ft) = this%damage_transitions(i,:,ft)/ &
       !         sum(this%damage_transitions(i,:,ft))
       ! end do

       do i = 1, ncrowndamage
          do j = 1, ncrowndamage
             if (j .eq. i) then
                this%damage_transitions(j,i,ft) = 1._r8
             else
                this%damage_transitions(j,i,ft) = 0._r8
             end if
          end do
       end do

       ! This is now annual - so each year 10% of trees will get some amount of damage
       this%damage_transitions(1,:,ft) = 0.025_r8
       

       write(fates_log(),'(a/,5(F12.6,1x))') 'JN transition matrix : ', this%damage_transitions(:,:,ft)
    end do

    return
  end subroutine InitDamageTransitions

end module FatesParameterDerivedMod
