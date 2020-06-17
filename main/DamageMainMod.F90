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

  public :: get_disturbance_collateral_damage_frac
  public :: get_disturbance_canopy_damage_frac
  public :: get_crown_reduction
 
  logical :: debug = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains


    subroutine get_disturbance_collateral_damage_frac(crowndamage, dist_frac)

    use EDTypesMod, only : ncrowndamagemax

    integer, intent(in) :: crowndamage
    real(r8), intent(out) :: dist_frac

    real(r8) :: t1         ! lower bound for integration
    real(r8) :: t2         ! upper bound for integration
    real(r8) :: rl_cd      ! crowndamage as a real
    real(r8) :: max         ! arbitrary large number 

    max = 1000000000.0_r8

   
    rl_cd = real(crowndamage)

    t1 = -1.0_r8*exp(-1.0_r8*rl_cd)
    t2 = -1.0_r8*exp(-1.0_r8*(rl_cd-1.0_r8))

    if(crowndamage == ncrowndamagemax) then
       t1 = -1.0_r8*exp(-1.0*max)
    end if

    dist_frac = t1-t2

    return
  end subroutine get_disturbance_collateral_damage_frac

  
  !------------------------------------------------------------------------------------
  ! This subroutine calculates damage of canopy trees
  
  subroutine get_disturbance_canopy_damage_frac(crowndamage, dist_frac)

    use EDTypesMod, only : ncrowndamagemax
    use FatesConstantsMod, only : years_per_day
    
    integer, intent(in) :: crowndamage
    real(r8), intent(out) :: dist_frac

    ! local variables
    real(r8) :: damage_fracs(ncrowndamagemax)

    damage_fracs = (/0.0_r8, 0.18_r8, 0.07_r8, 0.04_r8, 0.03_r8/)
    damage_fracs = damage_fracs * 0.001_r8
    
    damage_fracs(1) = 1.0_r8 - sum(damage_fracs(2:ncrowndamagemax))

    dist_frac = damage_fracs(crowndamage)
    
    
    return
  end subroutine get_disturbance_canopy_damage_frac

!-----------------------------------------------------------------------------------
  
  subroutine get_crown_reduction(crowndamage, crown_reduction)

    !------------------------------------------------------------------                                                                     
    ! This function takes the crown damage class of a cohort (integer)                                                                      
    ! and returns the fraction of the crown that is lost                                                                                    
    ! Since crowndamage class = 1 means no damage, we subtract one                                                                          
    ! before multiplying by 0.2                                                                                                             
    ! Therefore, first damage class is 20% loss of crown, second 40% etc.                                                                   
    !-------------------------------------------------------------------                                                                    
    use EDTypesMod     , only : ncrowndamagemax

    integer(i4), intent(in)   :: crowndamage
    real(r8),    intent(out)  :: crown_reduction

    crown_reduction = min(1.0_r8, (real(crowndamage) - 1.0_r8) * 0.2_r8)

    return
  end subroutine get_crown_reduction

  !----------------------------------------------------------------------------------------                                                 

  ! subroutine get_crown_damage_bounds(crowndamage, lower, upper)

  !   use EDTypesMod       , only : ncrowndamagemax

  !   integer(i4), intent(inout)  :: crowndamage
  !   real(r8)   , intent(out) :: lower
  !   real(r8)   , intent(out) :: upper

  !   integer(i4)              :: cdplus

  !   cdplus = crowndamage+1

  !   call get_crown_reduction(crowndamage, lower)
  !   call get_crown_reduction(cdplus, upper)

  !   return

  ! end subroutine get_crown_damage_bounds


  !----------------------------------------------------------------------------------------


end module DamageMainMod

