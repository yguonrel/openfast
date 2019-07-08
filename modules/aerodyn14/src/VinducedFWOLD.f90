SUBROUTINE VinducedFWOLD( n1, k1, n2, k2 )

  USE NWTC_Num, Only: TwoPi_D
  USE FVW_Parm, Only: delta_psi_Est
  USE MultTurb_params, Only: FWake!GCoord

  IMPLICIT NONE

  INTEGER n1, n2, k1, k2

  FWake%VinducedFarWakej( :, k2, n2, k1, n1 ) = FWake%VinducedFarWakejm1( :, k2, n2, k1-1, n1 ) + &
     & FWake%VinducedFarWakej( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est ), n1 )- &
     & FWake%VinducedFarWakejm1( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est )-1, n1 )


END SUBROUTINE VinducedFWOLD
