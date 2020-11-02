module DamageMainMod

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : i4 => fates_int
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : years_per_day
  use FatesGlobals          , only : fates_log

  use EDPftvarcon           , only : EDPftvarcon_inst

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

  public :: get_crown_reduction
  public :: get_crown_damage
  public :: adjust_bdead
  public :: get_damage_frac
  
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains


  subroutine get_damage_frac(cc_cd, nc_cd, pft, dist_frac)


    ! There is a transition matrix describing the rates of cohorts
    ! increasing in damage. Movement from one damage class to a higher
    ! damage class is described by the negative exponential. Both functions are in
    ! FatesParameterDerivedMod.
    
    ! USES
    use FatesInterfaceTypesMod, only : ncrowndamage
    use FatesConstantsMod, only : years_per_day
    use FatesParameterDerivedMod, only : param_derived

       
    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    integer, intent(in) :: pft
    real(r8), intent(out) :: dist_frac             ! probability of current cohort moving to new damage level

    dist_frac = param_derived%damage_transitions(cc_cd, nc_cd, pft) * years_per_day
    
    
  end subroutine get_damage_frac

  !-------------------------------------------------------    
  
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

