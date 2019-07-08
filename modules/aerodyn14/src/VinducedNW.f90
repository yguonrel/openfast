SUBROUTINE VinducedNW( rblade, Gamma, rp, Vind, rblade2, up )


 !******************************************************************
 !* Computes induced velocity on referance LAGRANGIAN MARKER # "m" due to all of the near wake points
 !*     --> include near wake away from blade up to NnearMax-1 and along blade from Num_Start to NumBS+1
 !*  *** When called in UpdateAeroVals, this instead computes the induced velocity ON THE BLADE
 !*                                     due to the near wake
 !*
 !* dV = F_c*(GAMMA_mu/4pi)*(r_1xr_2)*[1/mag(r_1)+1/mag(r_2)]*{1/[mag(r_1)*mag(r_2)+r_1DOTr_2]}
 !* F_c = mag(r)^2/sqrt[mag(r)^4+r_c^4]
 !* r_c( zeta, epsilon ) = sqrt{r_c0^2+(4*alpha*delta*mu*zeta/Omega)INTEGRAL[(1+epsilon)^-1*dzeta]}
 !*                        ^*                 diffusion          *^ ^*  effect of strain field  *^
 !* r_c0^2 = 4*alpha*delta*mu*zeta_0/Omega; zeta_0 = 30 degrees (location of near-wake roll-up)
 !*                                            ^^^^^hardcoded.
 !*
 !* where:
 !*      dV = Vind          <--- induced velocity in NEAR WAKE (Outputted to VinducedNWx)
 !*      r_blade1, r_blade2 <---distances between the point of interest and the end points of the vortex segment)
 !*                 and these points come from r_nearj, r_nearjm1, or r_nearjm2
 !*      GAMMA = Gamma; comes from Gamma_nearj or Gamma_nearjm1
 !*      INTEGRAL[1+epsilon)^-1*dzeta] = INTEGRAL
 !*            = rp <--- reference point
 !*                  from PREDICTOR: r_oldj(:,m,n), r_oldj(:,m-1,n), r_oldj(:,m,n), r_oldjm1(:,m-1,n)
 !*                       CORRECTOR: r_primej(:,m,n), r_primej(:,m-1,n), r_primejm1(:,m,n), r_primejm1(:,m-1,n)
 !*                       OTHER:     BladeThreeQuarterChord(:, nbs, n), 
 !*                             
 !*       
 !*       
 !* 
 !******************************************************************


  USE NWTC_Num,        Only: Pi_D, D2R_D
  USE MathOps,         Only: Norm, Dot, Cross
  USE MultTurb_Params, Only: NumWakes, NTurb
  USE FVW_Parm

  IMPLICIT NONE

  INTEGER,                                             INTENT( IN    ) :: up
  REAL( ReKi ), DIMENSION( 3 ),                              INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION( NnearMax, NumBS+1, NumWakes    ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NnearMax, NumWakes ), INTENT( IN    ) :: rblade, rblade2

  REAL( ReKi ), DIMENSION( 3 ),                              INTENT(   OUT ) :: Vind

  INTEGER :: i, kx, indx, nbs

  REAL( ReKi ) :: delta, strain, len1, len2, zeta, rc0, denom, close_to_zero, mag_r1, mag_r2, dotr1r2, rc, INTEGRAL

  REAL( ReKi ), DIMENSION( 3 ) :: r1, r2, crossr1r2
  REAL( ReKi ), ALLOCATABLE, DIMENSION( : ) :: rbladetemp

  ALLOCATE( rbladetemp( 3 ))
  INTEGRAL = 0.00_ReKi; rc = 0.00_ReKi; Vind = 0.00_ReKi

  DO i = 1, up      !   NumWakes for pred/corr and NumBl for UpdateAero & Circulation
     DO kx = 1, NnearMax - 1
        zeta = dble( kx ) * delta_psi(1)
        DO nbs = Num_start, NumBS+1
              r2( : ) = rp( : ) - rblade( :, nbs, kx,   i )
              r1( : ) = rp( : ) - rblade( :, nbs, kx+1, i )

           CALL NORM( r1, mag_r1 )
           CALL NORM( r2, mag_r2 )
           CALL DOT( r1, r2, dotr1r2 )
		   
           !Calculate the core radius for the Vind cut off distance
           delta = 1.00_ReKi + a1 * ( abs( Gamma( kx, nbs, i ))) / nu
           rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )

           IF ( kx .GE. 2 ) THEN
              rbladetemp = 0.00_ReKi
              rbladetemp = rblade(  :, nbs, kx+1, i ) - rblade(  :,  nbs,  kx, i )
              CALL NORM( rbladetemp, len2 )
              rbladetemp = 0.00_ReKi
              rbladetemp = rblade2( :, nbs, kx,   i ) - rblade2( :, nbs, kx-1, i )
              CALL NORM( rbladetemp, len1 )
              strain = ( len2 - len1 ) / len1
              INTEGRAL = delta_psi(1) / ( 1.0_ReKi + strain )
              rc = sqrt(( rc0 * rc0 + 4.00_ReKi * alpha_param * delta * nu * zeta / &
                 & Omega ) * INTEGRAL )
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF
		   
           CALL CROSS( r1, r2, crossr1r2 )
		   
           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
              & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.00_ReKi * dotr1r2 ) **2.00_ReKi

           IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
              & ( denom .NE. 0.00_ReKi )) THEN
           	
              DO indx = 1, 3
                 Vind( indx ) = Vind( indx ) + Gamma( kx, nbs, i ) / ( 4.00_ReKi * Pi_D ) * &
                    & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.00_ReKi - dotr1r2 / &
                    & (mag_r1 * mag_r2 )) / sqrt( denom )
              END DO
           END IF
        END DO
     END DO
  END DO
  DEALLOCATE( rbladetemp)
END SUBROUTINE VinducedNW
