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

  logical, protected :: damage_time  ! if true then damage occurs during current time step

  public :: get_crown_reduction
  public :: get_crown_damage
  public :: adjust_bdead
  public :: get_damage_frac
  public :: is_it_damage_time
  public :: damage_time
  public :: get_damage_mortality
  
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains


  subroutine is_it_damage_time(is_master, currentSite)

    !----------------------------------------------------------------------------
    ! This subroutine determines whether damage should occur (it is called daily)
    !-----------------------------------------------------------------------------

    use FatesInterfaceTypesMod , only : hlm_day_of_year

    integer, intent(in) :: is_master
    type(ed_site_type), intent(inout), target :: currentSite
    
    
    damage_time = .false.

    if (hlm_day_of_year .eq. 1) then
       damage_time = .true.
    end if
   
  end subroutine is_it_damage_time
  
  !----------------------------------------------------------------------------

  subroutine get_damage_frac(cc_cd, nc_cd, pft, dist_frac)


    ! given current cohort damage class find the fraction of individuals
    ! going to the new damage class.
    ! Consults a look up table of transitions from param derived. 
    
    ! USES
    use FatesInterfaceTypesMod, only : ncrowndamage
    use FatesConstantsMod, only : years_per_day
    use FatesParameterDerivedMod, only : param_derived

       
    ! ARGUMENTS
    integer, intent(in) :: cc_cd                   ! current cohort crown damage
    integer, intent(in) :: nc_cd                   ! new cohort crown damage
    integer, intent(in) :: pft
    real(r8), intent(out) :: dist_frac             ! probability of current cohort moving to new damage level

    dist_frac = param_derived%damage_transitions(cc_cd, nc_cd, pft) !* years_per_day (if damage is occuring annually don't do this)
    
    
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

    ! local variables
    real(r8) :: class_width

    class_width = 1.0_r8/ncrowndamage
    crown_reduction = min(1.0_r8, (real(crowndamage) - 1.0_r8) * class_width)

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
    real(r8) :: class_width

    frac = min(leaf_c/target_leaf_c, 1.0_r8)
    class_width = 1.0_r8/ncrowndamage

    crowndamage = max(1.0_r8, real(ceiling((1.0_r8-frac)/class_width)))

    return
  end subroutine get_crown_damage
  

  !----------------------------------------------------------------------------------------

  subroutine adjust_bdead(bt_sap, dbt_sapdd, bt_agb, dbt_agbdd, agb_frac, branch_frac, &
    crown_reduction)

    ! This subroutine scales target allometries to the damaged state -
    ! this is only for use in ForceDBH where we compare actual and target structural biomass
    ! to find a dbh. Since actual structural biomass has been reduced target biomass needs to also
    ! be reduced or the search algorimthm in ForceDBH will fail. This is only called when crown damage
    ! is greater than 1. 

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

  subroutine get_damage_mortality(crowndamage,pft, dgmort)

    use FatesInterfaceTypesMod     , only : ncrowndamage
    use EDPftvarcon                , only : EDPftvarcon_inst
    
    integer(i4), intent(in) :: crowndamage
    integer(i4), intent(in) :: pft
    real(r8),    intent(out) :: dgmort

    ! local variables
    integer(i4) :: i
    real(r8), allocatable :: dgmort_vec(:)
    real(r8) :: damage_mort_p1
    real(r8) :: damage_mort_p2

    ! parameter to determine slope of exponential
    damage_mort_p1 = EDPftvarcon_inst%damage_mort_p1(pft)
    damage_mort_p2 = EDPftvarcon_inst%damage_mort_p2(pft)
    
    ! JN - set up a vector of damage mortality values 
    allocate(dgmort_vec(1:ncrowndamage))
    do i = 1,ncrowndamage
       dgmort_vec(i) = real(i)
    end do

    ! JN - could make these proper cases?
    
    ! 1.  JN power function
   ! dgmort_vec = dgmort_vec**damage_mort_p1
    ! JN but must be bounded between 0 and 1
   ! dgmort_vec = dgmort_vec/sum(dgmort_vec)

    ! 2. JN logistic function
    dgmort_vec = 1.0_r8 / (1.0_r8 + exp(-1.0_r8 * damage_mort_p2 * &
         (dgmort_vec - damage_mort_p1) ) )

    
    if (crowndamage .eq. 1 ) then
       dgmort = 0.0_r8
    else
       dgmort = dgmort_vec(crowndamage)
    end if

    return
  end subroutine get_damage_mortality
  !----------------------------------------------------------------------------------------
  
end module DamageMainMod

