SUBROUTINE VinducedFW3(n1, k1, n2, k2)

  USE NWTC_Num,        Only: TwoPi_D
  USE FVW_Parm,       Only: delta_psi_Est
  USE MultTurb_params, Only: FWake!GCoord

  IMPLICIT NONE

  INTEGER, INTENT( IN ) :: n1, n2, k1, k2
  FWake%VinducedFarWakeRj( :, k2, n2, k1, n1 ) = FWake%VinducedFarWakeRjm1( :, k2, n2, k1-1, n1 ) + &
     & FWake%VinducedFarWakeRj( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est ), n1 ) - &
     & FWake%VinducedFarWakeRjm1( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est )-1, n1 )

END SUBROUTINE VinducedFW3
