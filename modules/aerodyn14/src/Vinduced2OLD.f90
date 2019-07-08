SUBROUTINE Vinduced2OLD( rblade, Gamma, rp, rblade2, n, jold, kold )

  USE NWTC_Num,         Only: Pi_D, TwoPi_D, D2R_D
  USE FVW_Parm
  USE MathOps,          Only: Norm, Dot, Cross
  USE MultTurb_Params, 	Only: NumWakes, NTurb, FWake!GCoord

  IMPLICIT NONE

  INTEGER,                                           INTENT( IN    ) :: n, jold, kold
  REAL( ReKi ), DIMENSION( 3 ),                            INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: rblade, rblade2


  INTEGER :: i, kx, indx, q, limit, BC_loc

  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2
  REAL( ReKi ) :: delta, strain, len1, len2, zeta, rc0, Denom
  REAL( ReKi ), ALLOCATABLE, DIMENSION( : ) :: r1, r2, crossr1r2, rbladetemp
  REAL( ReKi ) :: INTEGRAL, rc

  REAL :: WhichTurb

  ALLOCATE( rbladetemp( 3 ), r1(3), r2(3), crossr1r2(3))

  INTEGRAL = 0.00_ReKi; rc=0.00_ReKi

  WhichTurb = REAL(n-0.01)/REAL(NumBl)
  q = n - FLOOR( WhichTurb )*NumBl

  DO i = 1, NumWakes
     WhichTurb = REAL(i-0.01)/REAL(NumBl)

     limit = CUTOFF_up( CEILING( WhichTurb ))
     BC_loc = BC( CEILING( WhichTurb ) )
     DO kx = NnearMax, WakeAgeLimit-1
        IF ( kx .LT. limit ) THEN
           zeta = dble( kx ) * delta_psi(1)

           DO indx = 1, 3
              r2( indx ) = rp( indx ) - rblade( indx, kx,   i )
              r1( indx ) = rp( indx ) - rblade( indx, kx+1, i )
           END DO

           CALL NORM( r1, mag_r1 )
           CALL NORM( r2, mag_r2 )
           CALL DOT( r1, r2, dotr1r2 )

           delta = 1.00_ReKi + a1 * ( abs( Gamma( kx, i ))) / nu
           rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade(  :, kx+1, i ) - rblade(  :, kx, i )
           CALL NORM( rbladetemp, len2 )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade2( :, kx,   i ) - rblade2( :, kx-1, i )
           CALL NORM( rbladetemp, len1 )

           strain = ( len2 - len1 ) / len1

           INTEGRAL = delta_psi(1) / ( 1.00_ReKi + strain )
           rc = sqrt(( rc0 * rc0 + 4.00_ReKi * alpha_param * delta * nu * zeta / &
              & Omega ) * INTEGRAL )

           IF ( kx .GE. 2 ) THEN
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           CALL CROSS( r1, r2, crossr1r2 )

           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
              & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.00_ReKi * &
              & dotr1r2 ) **2.00_ReKi

           IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
              & (denom .NE. 0.00_ReKi)) THEN
              DO indx = 1, 3
                 FWake%VinducedFarWakej( indx, kold, q, kx+1, i ) = Gamma( kx, i ) / &
                    & ( 4.00_ReKi * Pi_D ) * crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * &
                    & ( 1.00_ReKi - dotr1r2 / ( mag_r1 * mag_r2 )) / sqrt(denom)
              END DO
           ELSE
              FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
           END IF
        ELSE IF ( jold .GE. ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc ) .AND. &
               & kx .LE. ( jold - ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc - limit )) .AND. &
               & kx .LT. ( limit  + NINT( TwoPi_D / delta_psi_Est ))) THEN
           CALL VinducedFWOLD( i, kx+1, q, kold )
        ELSE IF (kx .LE. WakeAgeLimit ) THEN
           FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
        END IF
     END DO
  END DO

  DEALLOCATE( r1, r2, crossr1r2, rbladetemp)

END SUBROUTINE Vinduced2OLD
