SUBROUTINE VinducedBC( rblade, Gamma, rp, Vind )


  USE NWTC_Num,         Only: Pi_D, D2R_D
  USE FVW_Parm
  USE MathOps,          Only: Norm, Dot, Cross
  USE MultTurb_params,  Only: NumWakes

  IMPLICIT NONE

  REAL( ReKi ), DIMENSION( 3 ), INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    NumBS+1, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NumWakes ), INTENT( IN    ) :: rblade

  INTEGER  i, indx, nbs

  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2
  REAL( ReKi ) :: delta, rc0, denom
  REAL( ReKi ), DIMENSION( 3 ) :: Vind, r1, r2, crossr1r2


  Vind = 0.00_ReKi

  DO i = 1, NumWakes
     DO nbs = Num_start, NumBS
        DO indx = 1, 3
           r1( indx ) = rp( indx ) - rblade( indx, nbs,   i )
           r2( indx ) = rp( indx ) - rblade( indx, nbs+1, i )
        END DO

        CALL NORM( r1, mag_r1 )
        CALL NORM( r2, mag_r2 )
        CALL DOT( r1, r2, dotr1r2 )

        delta = 1.00_ReKi + a1 * ( abs( Gamma( nbs, i ))) / nu
        rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )

        close_to_zero = rc0

        CALL CROSS( r1, r2, crossr1r2 )
        
        denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
           & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - &
           & 2.00_ReKi * dotr1r2 ) **2.00_ReKi
        
        IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR.&
           & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
           & (denom .NE. 0.00_ReKi )) THEN

           DO indx = 1, 3
              Vind( indx ) = Vind( indx ) + Gamma( nbs, i ) / ( 4.00_ReKi * Pi_D ) * &
                 & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.00_ReKi - dotr1r2 / &
                 & ( mag_r1 * mag_r2 )) / sqrt( denom )
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE VinducedBC
