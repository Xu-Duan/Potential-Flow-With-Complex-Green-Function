module MGreenFunction
    USE MMesh, ONLY:TMesh
    USE MConstants, ONLY:Omega2k
    implicit none
    private :: PL2
    public :: GF1, GF2
contains
    subroutine GF1(Mesh, G1, GV1)
        type(TMesh), intent(in) :: Mesh
        real, dimension(Mesh%nPanels, Mesh%nPanels), intent(out) :: G1
        real, dimension(Mesh%nPanels, Mesh%nPanels, 3), intent(out) :: GV1
        ! local variables
        integer :: i, j, k
        real  :: rho, rho1, pi = acos(-1.0)
        G1 = 0.
        GV1 = 0.
        do i = 1, Mesh%nPanels
            do j = 1, Mesh%nPanels
                do k = 1, 4
                    rho = norm2(Mesh%center(:, i) - Mesh%GQPoints(:, k, j))
                    rho1 = sqrt(norm2(Mesh%center(1:2, i) - Mesh%GQPoints(1:2, k, j))**2 + &
                        &(Mesh%center(3, i) + Mesh%GQPoints(3, k, j))**2)
                    G1(i, j) = G1(i, j) + (1/rho - 1/rho1)* Mesh%GQWeights(k, j)
                    GV1(i, j, :) = GV1(i, j, :) + (Mesh%GQPoints(:, k, j) -&
                        & Mesh%center(:, i))/rho**3 * Mesh%GQWeights(k, j)
                    GV1(i, j, 1:2) = GV1(i, j, 1:2) - (Mesh%GQPoints(1:2, k, j) -&
                        & Mesh%center(1:2, i))/rho1**3 * Mesh%GQWeights(k, j)
                    GV1(i, j, 3) = GV1(i, j, 3) + (Mesh%GQPoints(3, k, j) +&
                        & Mesh%center(3, i))/rho1**3 * Mesh%GQWeights(k, j)
                enddo
                if (i == j) then
                    GV1(i, j, :) = GV1(i, j, :) - 2 * pi * Mesh%normalVector(:, i)
                endif
            enddo
        enddo
        G1 = -G1 / (4 * pi)
        GV1 = -GV1 / (4 * pi)
    endsubroutine GF1

    subroutine GF2(Mesh, omega, D1, D2, Z1, Z2, G2, GV2)
        type(TMesh), intent(in) :: Mesh
        real, intent(in) :: omega
        real, dimension(328, 46) :: D1, D2, Z1, Z2
        complex, dimension(Mesh%nPanels, Mesh%nPanels), intent(out) :: G2
        complex, dimension(Mesh%nPanels, Mesh%nPanels, 3), intent(out) :: GV2
        
        ! local variables
        integer :: i, j, k
        real  :: originalR(328), originalZ(46)  ! the points table
        real :: R, Z ! R is the horizontal distance from source to field, while Z is vertical
        real :: KR, KZ ! non-dimensional variables of R and Z
        real :: wavenumber, RR, temp, xl(3), zl(3), pi = acos(-1.0)
        ! xl and zl are weights of lagrange interpolation
        real :: interZ1, interZ2, interD1, interD2 ! these are interpolated values
        integer :: RI, ZI ! the index of where R and Z lies
        real :: epz, SQ, CSK, SIK

        wavenumber = Omega2k(omega)
        originalR(1) = 0.0
        do i = 2, 31
            originalR(i) = 10**((i-1.) / 5. - 6.)
        enddo
        do i = 32, 328
            originalR(i) = 4./3. + (i-32.)/3.0
        enddo
        ! Initialize Z
        do i = 1, 20
            originalZ(i) = -10**(i/5.-6.)
        enddo
        do i = 21, 46
            originalZ(i) = -10**(i/8.-4.5)
        enddo
        G2 = complex(0., 0.)
        GV2 = complex(0., 0.)

        do i = 1, Mesh%nPanels
            do j = 1, Mesh%nPanels
                do k = 1, 4
                    R = norm2(Mesh%center(1:2, i) - Mesh%GQPoints(1:2, k, j))
                    Z = Mesh%center(3, i) + Mesh%GQPoints(3, k, j)
                    KR = wavenumber * R
                    KZ = wavenumber * Z
                    RR = sqrt(R**2 + Z**2)
                    temp = pi / (wavenumber * RR)**3
                    if (KZ > -16) then      ! kz is normal
                        if (KR < 99.7) then   ! kr is normal
                            ! following code is used to determine where R and Z lies
                            IF (KR < 1) THEN           !     0 < KR < 1
                                RI = INT(5*(log10(KR)+6)+1)
                            ELSE                       !    1 < KR < 99.7
                                RI = INT(3*KR+28)
                            ENDIF
                            RI = MAX(MIN(RI, 327), 2)   ! make sure RI between 2 and 327
                            IF (KZ < -1e-2) THEN       !   -16 < KZ < -1e-2
                                ZI = INT(8*(log10(-KZ)+4.5))
                            ELSE                        ! -1e-2 < KZ < -1.5e-6
                                ZI = INT(5*(log10(-KZ)+6))
                            ENDIF
                            ZI = MAX(MIN(ZI, 45), 2)   ! make sure ZI between 2 and 45
                            ! RI and ZI are determined
                            ! perform interpolation
                            XL(1) = PL2(originalR(RI-1),   originalR(RI), originalR(RI+1), KR)
                            XL(2) = PL2(originalR(RI), originalR(RI-1), originalR(RI+1),   KR)
                            XL(3) = PL2(originalR(RI+1), originalR(RI-1),   originalR(RI), KR)
                            ZL(1) = PL2(originalZ(ZI-1),   originalZ(ZI), originalZ(ZI+1), KZ)
                            ZL(2) = PL2(originalZ(ZI), originalZ(ZI-1), originalZ(ZI+1),   KZ)
                            ZL(3) = PL2(originalZ(ZI+1), originalZ(ZI-1),   originalZ(ZI), KZ)
                            interZ1 = DOT_PRODUCT(XL, MATMUL(Z1(RI-1:RI+1, ZI-1:ZI+1), ZL))
                            interZ2 = DOT_PRODUCT(XL, MATMUL(Z2(RI-1:RI+1, ZI-1:ZI+1), ZL))
                            interD1 = DOT_PRODUCT(XL, MATMUL(D1(RI-1:RI+1, ZI-1:ZI+1), ZL))
                            interD2 = DOT_PRODUCT(XL, MATMUL(D2(RI-1:RI+1, ZI-1:ZI+1), ZL))
                        else ! kr is large, use asympotic expression
                            EPZ  = EXP(KZ)
                            SQ   = SQRT(2*PI/KR)
                            CSK  = COS(KR-PI/4)
                            SIK  = SIN(KR-PI/4)
                            
                            interZ1 = temp*KZ - PI*EPZ*SQ*SIK
                            interZ2 =                EPZ*SQ*CSK
                            interD1 = PI*EPZ*SQ*(CSK - 0.5*SIK/KR) - temp*KR
                            interD2 =    EPZ*SQ*(SIK + 0.5*CSK/KR)
                        endif
                        G2(i, j) = G2(i, j) + (interZ1 * 2 * wavenumber / pi + &
                            & (0., 1.) * 2 * wavenumber * interZ2) * Mesh%GQWeights(k, j)
                        GV2(i, j, 1:2) = GV2(i, j, 1:2) + (Mesh%GQPoints(1:2, k, j) - Mesh%center(1:2, i)) / R * &
                            & ((interD1 + temp * KR) * 2 * wavenumber**2 / pi + &
                            & (0., 1.) * 2 * wavenumber**2 * interD2) * Mesh%GQWeights(k, j)
                        GV2(i, j, 3) = GV2(i, j, 3) + ((interZ1 - temp * KZ) * 2 * wavenumber**2 / pi + &
                            & (0., 1.) * 2 * wavenumber**2 * interZ2) * Mesh%GQWeights(k, j)
                    else ! kz is large(negative), can be calculated directly in this situation
                        G2(i, j) = G2(i, j) + complex(temp/pi*2*wavenumber*kz, 0.)
                        GV2(i, j, :) = GV2(i, j, :) + complex(0., 0.)
                    endif
                enddo
            enddo
        enddo
        G2 = G2 / (-4 * pi)
        GV2 = GV2 / (-4 * pi)
    endsubroutine GF2

    real function PL2(x1, x2, x3, x0)
        real, intent(in) :: x1, x2, x3, x0
        PL2 = (x0 - x2) * (x0 - x3) / (x1 - x2) / (x1 - x3)
    endfunction PL2
endmodule MGreenFunction