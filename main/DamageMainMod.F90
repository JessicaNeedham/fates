module DamageMainMod

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : i4 => fates_int
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : years_per_day
!  use FatesInterfaceMod     , only : hlm_freq_day
  use FatesGlobals          , only : fates_log

  use EDPftvarcon           , only : EDPftvarcon_inst

  use EDTypesMod            , only : element_pos
  use EDtypesMod            , only : ed_site_type
  use EDtypesMod            , only : ed_patch_type
  use EDtypesMod            , only : ed_cohort_type
  use EDtypesMod            , only : AREA

  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : SetState

  implicit none
  private

  !public :: get_disturbance_collateral_damage_frac
  !public :: get_disturbance_canopy_damage_frac
  public :: get_crown_reduction
  public :: get_crown_damage
  public :: adjust_bdead
  public :: damage_function
  public :: get_damage_frac
  
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains

  real(r8) function damage_function(crowndamage)

    real(r8), intent(in) :: crowndamage ! plant crowndamage level

    ! hard coded here - could have different cases and parameters?
    damage_function = exp(-0.2*crowndamage)
    return
  end function damage_function

  !---------------------------------------------------------------------------

  subroutine get_damage_frac(cc_cd, nc_cd, dist_frac)

    ! There is a transition matrix describing the probabilities of cohorts
    ! increasing in damage. Movement from one damage class to a higher
    ! damage class is described by the damage_function above - this form is yet
    ! to be decided on or properly parameterised.
    ! movement from i to i+1, i+2 etc is given by the vector transition_probs
    ! this would populate the lower diagonal of a  matrix of dim(ncrowndamage,ncrowndamage)
    ! the matrix columns need to be normalised so that transitions from any damage class
    ! (including staying in the same damage class) sum to 1.
    ! since this function takes current cohort damage class and retuns the probability of
    ! moving to a new damage class, we don't actually make the whole matrix here
    
    ! USES
    use FatesInterfaceTypesMod, only : ncrowndamage

    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    real(r8), intent(out) :: dist_frac             ! probability of current cohort moving to new damage level

    ! LOCAL VARIABLES
    real(r8), allocatable :: transition_probs(:)   ! vector of transition probabilities between damage classes
    integer                :: i                    ! counter
    real(r8)               :: i_real               ! real version of counter
    integer                :: n                    ! how long vector of transition probs needs to be

    n = ncrowndamage + 1 - cc_cd
    allocate(transition_probs(n))

    do i = 1,n
       i_real = real(i)
       transition_probs(i) = damage_function(i_real)
    end do

    ! normalise it
    transition_probs = transition_probs/sum(transition_probs)

    i = nc_cd - cc_cd + 1
    dist_frac = transition_probs(i)

    return
  end subroutine get_damage_frac
  
  !-------------------------------------------------------------------------------

    
  
  

! subroutine get_disturbance_collateral_damage_frac(crowndamage, dist_frac)

!     use FatesInterfaceTypesMod, only : ncrowndamage

!     integer, intent(in) :: crowndamage
!     real(r8), intent(out) :: dist_frac

!     real(r8) :: t1         ! lower bound for integration                                 
!     real(r8) :: t2         ! upper bound for integration
!     real(r8) :: rl_cd      ! crowndamage as a real
!     real(r8) :: max         ! arbitrary large number                                                                   

!     max = 1000000000.0_r8
!     rl_cd = real(crowndamage)

!     t1 = negative_exponential(rl_cd)
!     t2 = negative_exponential(rl_cd-1.0_r8)

!     if(crowndamage == ncrowndamage) then
!        t1 = negative_exponential(max)
!     end if

!     dist_frac = t2-t1

!     return
!   end subroutine get_disturbance_collateral_damage_frac


  
!   !------------------------------------------------------------------------------------
!   ! This subroutine calculates damage of canopy trees
  
!   subroutine get_disturbance_canopy_damage_frac(crowndamage, ft, dist_frac)

!     use FatesInterfaceTypesMod, only : ncrowndamage
!     use FatesConstantsMod, only : years_per_day
    
!     integer, intent(in) :: crowndamage
!     real(r8), intent(out) :: dist_frac
!     integer, intent(in) :: ft
    
!     ! local variables
!     real(r8) :: damage_fracs(ncrowndamage)

!     damage_fracs = EDPftvarcon_inst%crowndamage_fracs(ft,:)
!     damage_fracs = damage_fracs * years_per_day
!     damage_fracs(1) = 1.0_r8 - sum(damage_fracs(2:ncrowndamage))
    
!     dist_frac = damage_fracs(crowndamage)
    
    
!     return
!   end subroutine get_disturbance_canopy_damage_frac

!-----------------------------------------------------------------------------------
  
  subroutine get_crown_reduction(crowndamage, crown_reduction)

    !------------------------------------------------------------------                                                                     
    ! This function takes the crown damage class of a cohort (integer)                                                                      
    ! and returns the fraction of the crown that is lost                                                                                    
    ! Since crowndamage class = 1 means no damage, we subtract one                                                                          
    ! before multiplying by 0.2                                                                                                             
    ! Therefore, first damage class is 20% loss of crown, second 40% etc.                                                                   
    !-------------------------------------------------------------------                                                                    
    use FatesInterfaceTypesMod     , only : ncrowndamage

    integer(i4), intent(in)   :: crowndamage
    real(r8),    intent(out)  :: crown_reduction

    crown_reduction = min(1.0_r8, (real(crowndamage) - 1.0_r8) * 0.2_r8)

    return
  end subroutine get_crown_reduction

  !----------------------------------------------------------------------------------------                                                 
  subroutine get_crown_damage(leaf_c, target_leaf_c, crowndamage)

    ! This subroutine calculates which damage class a cohort should be in
    ! based on leaf biomass relative to target leaf biomass - i.e. it allows
    ! for cohorts to recover

    use FatesInterfaceTypesMod          , only : ncrowndamage

    real(r8),    intent(in)   :: leaf_c
    real(r8),    intent(in)   :: target_leaf_c
    integer(i4), intent(out)  :: crowndamage

    ! Local variables
    real(r8) :: frac

    frac = min(leaf_c/target_leaf_c, 1.0_r8)

    crowndamage = max(1.0_r8, real(ceiling((1.0_r8-frac)/0.2_r8)))

    return
  end subroutine get_crown_damage
  

  !----------------------------------------------------------------------------------------

  subroutine adjust_bdead(bt_sap, dbt_sapdd, bt_agb, dbt_agbdd, agb_frac, branch_frac, &
    crown_reduction)

    ! This subroutine takes structural biomass from pft%GetState
    ! and scales it up to give what it should be for an undamaged tree

    real(r8), intent(inout) :: bt_sap
    real(r8), intent(inout) :: dbt_sapdd
    real(r8), intent(inout) :: bt_agb
    real(r8), intent(inout) :: dbt_agbdd
    real(r8), intent(in) :: agb_frac
    real(r8), intent(in) :: branch_frac
    real(r8), intent(in) :: crown_reduction
   
    
    bt_sap = bt_sap * agb_frac * branch_frac * (1.0_r8 - crown_reduction)
    dbt_sapdd = dbt_sapdd * agb_frac * branch_frac * (1.0_r8 - crown_reduction)

    bt_agb = bt_agb * branch_frac * (1.0_r8 - crown_reduction)
    dbt_agbdd = dbt_agbdd * branch_frac * (1.0_r8 - crown_reduction)

    return
  end subroutine adjust_bdead
  

  !----------------------------------------------------------------------------------------

end module DamageMainMod

