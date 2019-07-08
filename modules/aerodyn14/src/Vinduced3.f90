SUBROUTINE Vinduced3( rblade, Gamma, rp, rblade2, n, jold, kold )

 !******************************************************************
 !* Computes induced velocity on NEAR WAKE due to the Far Wake
 !*
 !*  
 !******************************************************************
		
  USE FVW_Parm
  USE Precision
  USE NWTC_Num,         Only: Pi_D, TwoPi_D, D2R_D
  USE MathOps,          Only: Norm, Dot, Cross
  USE MultTurb_Params,  Only: NumWakes, NTurb, FWake!GCoord

  IMPLICIT NONE

  INTEGER, INTENT( IN ) :: n, jold, kold	!=Blade #
  REAL( ReKi ), DIMENSION( 3                            ), INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: rblade, rblade2



  INTEGER i, kx, indx, limit, BC_loc, a
  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2, delta, strain
  REAL( ReKi ) :: len1, len2, zeta, rc0, denom, rc, INTEGRAL
  REAL( ReKi ), DIMENSION( 3 ) :: r1, r2, crossr1r2
  REAL( ReKi ), ALLOCATABLE, DIMENSION( : ) :: rbladetemp

  REAL :: WhichTurb

  ALLOCATE( rbladetemp( 3 ))
  INTEGRAL = 0.0_ReKi; rc=0.0_ReKi
  DO i = 1, NumWakes
     WhichTurb = REAL(i-0.01)/REAL(NumBl)
        limit = CUTOFF_up( CEILING( WhichTurb ))
        BC_loc = BC( CEILING( WhichTurb ) )

     DO kx = NnearMax, WakeAgeLimit-1
        IF ( kx .LT. limit ) THEN
           r2( : ) = rp( : ) - rblade( :, kx,   i )
           r1( : ) = rp( : ) - rblade( :, kx+1, i )
           
           CALL NORM( r1, mag_r1 )
           CALL NORM( r2, mag_r2 )
           CALL DOT(  r1, r2, dotr1r2 )

           delta = 1.0_ReKi + a1 * ( abs( Gamma( kx, i ))) / nu
           rc0 = sqrt( 4.0_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )
           rbladetemp = 0.0_ReKi
           rbladetemp = rblade(  :, kx+1, i ) - rblade(  :, kx, i )
           CALL NORM( rbladetemp, len2 )
           rbladetemp = 0.0_ReKi
           rbladetemp = rblade2( :, kx,   i ) - rblade2( :, kx-1, i )
           CALL NORM( rbladetemp, len1 )
              
           zeta = dble( kx ) * delta_psi(1)

           strain = ( len2 - len1 ) / len1
           INTEGRAL = delta_psi(1) / ( 1.0_ReKi + strain )
           rc = sqrt(( rc0 * rc0 + 4.0_ReKi * alpha_param * delta * nu * zeta / &
              & Omega ) * INTEGRAL )

           IF ( kx .GE. 2 ) THEN
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           CALL CROSS( r1, r2, crossr1r2 )
           
           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.0_ReKi + &
              & close_to_zero **4.0_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.0_ReKi * &
              & dotr1r2 ) **2.0_ReKi
           
           IF ( mag_r1 .GT. 0.0_ReKi .AND. mag_r2 .GT. 0.0_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.0_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.0_ReKi .OR. crossr1r2( 3 ) .NE. 0.0_ReKi ) .AND. &
              & ( denom .NE. 0.0_ReKi )) THEN
           		
              DO indx=1, 3
                 FWake%VinducedFarWakeRj( indx, kold, n, kx+1, i ) = Gamma( kx, i ) / ( 4.0_ReKi * Pi_D ) * &
                    & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.0_ReKi - dotr1r2 / &
                    & ( mag_r1 * mag_r2 )) / sqrt( denom )
              END DO
           ELSE
              FWake%VinducedFarWakeRj( :, kold, n, kx+1, i ) = 0.0_ReKi
           END IF         
        ELSE IF ( jold .GE. ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc ) .AND. &
           & kx .LE. ( jold - ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc - limit)) .AND. &
           & kx .LT. ( limit + NINT( TwoPi_D / delta_psi_Est ))) THEN
           CALL VinducedFW3( i, kx+1, n, kold )
        ELSE
           FWake%VinducedFarWakeRj( :, kold, n, kx+1, i ) = 0.0_ReKi
        END IF
     END DO
  END DO

  DEALLOCATE( rbladetemp)
END SUBROUTINE Vinduced3
