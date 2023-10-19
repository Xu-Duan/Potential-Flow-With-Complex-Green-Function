PROGRAM Main
    USE MMesh     ,     ONLY: TMesh, GenerateMesh
    USE MOutput   ,     ONLY: WriteFK, WriteSummaryMesh, WriteDiffractionVelocity, WriteRadiationVelocity
    USE MConstants,     ONLY: nw, wmax, wmin, nbeta, betamax, betamin, czero, rho
    USE MBodyCondition, ONLY: DiffractionCondition
    USE MGreenFunction, ONLY: GF1, GF2
    IMPLICIT NONE

    character(30) :: pointsFile = 'input_files/Points.txt'
    character(30) :: connectFile = 'input_files/Connect.txt'
    character(30) :: motionFile = 'input_files/Motion.txt'
    TYPE(TMesh)   :: Mesh
    Real, DIMENSION(:), allocatable          :: w    ! all omega
    REAL, DIMENSION(:),allocatable           :: beta ! all beta
    Complex, DIMENSION(:, :, :),allocatable  :: FK   ! FK Force
    Complex, DIMENSION(:), allocatable       :: pressure
    COMPLEX, DIMENSION(:), allocatable       :: velocity            ! velocity is temporary storation of normal velocity
    COMPLEX, DIMENSION(:, :, :), allocatable :: diffractionVelocity ! normal velocity in diffraction problem
    COMPLEX, DIMENSION(:, :, :), allocatable :: radiationVelocity   ! normal velocity in radiation problem
    real, DIMENSION(:, :),allocatable           :: G1
    complex, DIMENSION(:, :),allocatable        :: G2, G ! first and second component of Green function
    real, DIMENSION(:, :, :),allocatable        :: GV1
    complex, DIMENSION(:, :, :),allocatable     :: GV2, GV ! gradiant of first and second component of Green function
    complex, DIMENSION(:, :),allocatable        :: V, V2
    complex, DIMENSION(:),allocatable :: sigma, phi, solution  ! source strength and velocity potential
    real, DIMENSION(:),allocatable :: realsigma, imagsigma
    integer, DIMENSION(:),allocatable :: ipiv
    INTEGER       :: i, j, k, l, iflag
    real, dimension(328, 46) :: D1, D2, Z1, Z2
    real :: addms(6), dmp(6)
    complex :: msAndDmp(6), diffractionForce(6)
    allocate(w(nw), beta(nbeta), FK(nw, nbeta, 6))

    CALL GenerateMesh(pointsFile, connectFile, motionFile, Mesh)

    CALL WriteSummaryMesh('output_files/Mesh.txt', Mesh) ! output the summary of mesh

    allocate(pressure(Mesh%nPanels), velocity(Mesh%nPanels))
    allocate(diffractionVelocity(nw, nbeta, Mesh%nPanels), radiationVelocity(nw, 6, Mesh%nPanels))

    IF (nw .EQ. 1) THEN
        w(1) = wmin
    ELSE
        DO i = 1, nw
            w(i) = wmin + (wmax - wmin) / (nw - 1) * (i-1)
        END DO
    END IF
    IF (nbeta .EQ. 1) THEN
        beta(1) = betamin
    ELSE
        DO i = 1, nbeta
            beta(i) = betamin + (betamax - betamin) / (nbeta - 1) * (i-1)
        END DO
    END IF

    ! Calculate FK Force and Normal Velocity
    open(unit=1, file='output_files/Pressure.txt', action='write')
    DO i = 1, nw
        ! Diffraction Problem
        DO j = 1, nbeta
            CALL DiffractionCondition(w(i), beta(j), Mesh, pressure, velocity)
            DO l = 1, Mesh%nPanels
                write(1, '(2(E14.6), I6, 2(ES14.6))') w(i), beta(j), l, REAL(pressure(l)), IMAG(pressure(l))
            END DO
            DO k = 1, 6
                FK(i, j, k) = czero
                DO l = 1, Mesh%nPanels
                    FK(i, j, k) = FK(i, j, k) - pressure(l) * Mesh%UNVA(k, l)
                END DO
            END DO
            
            DO k = 1, Mesh%nPanels
                diffractionVelocity(i, j, k) = conjg(velocity(k))
            END DO
        END DO

        
        ! Radiation Problem
        DO j = 1, 6
            DO k = 1, Mesh%nPanels
                ! Note here that Mesh%UNV is real, while radiationVelocity is complex
                radiationVelocity(i, j, k) = Mesh%UNV(j, k)
            END DO
        END DO
    END DO
    close(1)
    CALL WriteFK('output_files/FKForce.txt', w, beta, FK) ! output the summary of FKForce
    CALL WriteDiffractionVelocity('output_files/DiffractionVelocity.txt', w, beta, Mesh%nPanels, diffractionVelocity)
    CALL WriteRadiationVelocity('output_files/RadiationVelocity.txt', w, Mesh%nPanels, radiationVelocity)
    deallocate(pressure, velocity)

    ! Solve the problem
    open(1, file = 'tableForDel.txt')
    do i = 1, 328
        read(1, '(46ES14.6)') (D1(i, j), j = 1,46)
    enddo
    read(1, *)
    do i = 1, 328
        read(1, '(46ES14.6)') (D2(i, j), j = 1,46)
    enddo
    read(1, *)
    do i = 1, 328
        read(1, '(46ES14.6)') (Z1(i, j), j = 1,46)
    enddo
    read(1, *)
    do i = 1, 328
        read(1, '(46ES14.6)') (Z2(i, j), j = 1,46)
    enddo
    close(1)
    allocate(G1(Mesh%nPanels, Mesh%nPanels), G2(Mesh%nPanels, Mesh%nPanels), G(Mesh%nPanels, Mesh%nPanels),&
            &GV1(Mesh%nPanels, Mesh%nPanels, 3), GV2(Mesh%nPanels, Mesh%nPanels, 3),&
            & GV(Mesh%nPanels, Mesh%nPanels, 3), V(Mesh%nPanels, Mesh%nPanels), V2(Mesh%nPanels, Mesh%nPanels)&
            ,sigma(Mesh%nPanels), solution(Mesh%nPanels), phi(Mesh%nPanels), ipiv(Mesh%nPanels),&
            & realsigma(Mesh%nPanels), imagsigma(Mesh%nPanels))
    CALL GF1(Mesh, G1, GV1)
    ! =============== validate ==========================
    open(1, file = 'Green1.txt')
    do i = 1, Mesh%nPanels
        do j = 1, Mesh%nPanels
            write(1, '(4ES14.6)') G1(i, j), GV1(i, j, 1), GV1(i, j, 2), GV1(i, j, 3)
        enddo
    enddo
    close(1)
    ! ===================================================
    DO i = 1, nw
        call GF2(Mesh, w(i), D1, D2, Z1, Z2, G2, GV2)
        G = G1 + G2
        GV = GV1 + GV2
        open(1, file = 'Green2.txt')
        do j = 1, Mesh%nPanels
            do k = 1, Mesh%nPanels
                write(1, '(8ES14.6)') real(G2(j, k)), imag(G2(j, k)), real(GV2(j, k, 1)),&
                & imag(GV2(j, k, 1)), real(GV2(j, k, 2)), imag(GV2(j, k, 2)), real(GV2(j, k, 3)), imag(GV2(j, k, 3))
            enddo
        enddo
        close(1)
        do j = 1, Mesh%nPanels
            do k = 1, Mesh%nPanels
                V(j, k) = dot_product(GV(j, k, :), Mesh%normalVector(:, j))
            enddo
        enddo
        open(1, file = 'Green.txt')
        do j = 1, Mesh%nPanels
            do k = 1, Mesh%nPanels
                write(1, '(4ES14.6)') real(G(j, k)), imag(G(j, k)), real(V(j, k)), imag(V(j, k))
            enddo
        enddo
        close(1)
        V2 = V
        solution = diffractionVelocity(i, 1, :)
        call cgesv(Mesh%nPanels, 1, V2, Mesh%nPanels, ipiv, solution, Mesh%nPanels, iflag)
        sigma = -solution
        !open(1, file='sourceStrength.txt')
        !do j = 1, Mesh%nPanels
        !    read(1, '(2(ES14.6))') realsigma(j), imagsigma(j)
        !enddo
        !close(1)
        !sigma = cmplx(realsigma, imagsigma)
        phi = matmul(G, sigma)
    enddo
    !msAndDmp = cmplx(0., 0.)
    !do j = 1, 6
    !    do i = 1, Mesh%nPanels
    !        msAndDmp(j) = rho * phi(i) * Mesh%UNV(j, i)
    !    enddo
    !    addms(j) = real(msAndDmp(j))
    !    dmp(j) = aimag(msAndDmp(j)) * w(1)
    !enddo
    !open(1, file = 'addedMassAndDamping.txt')
    !do j = 1, 6
    !    write(1, '(ES14.6, 2(I5), 2(ES14.6))') w(1), j, 1, addms(j), dmp(j)
    !enddo
    !close(1)
    open(1, file='sourceStrengthMy.txt')
        do j = 1, Mesh%nPanels
            write(1, '(4(ES14.6))') real(sigma(j)), imag(sigma(j)), real(phi(j)), imag(phi(j))
        enddo
        close(1)
    diffractionForce = cmplx(0., 0.)
    do j=1,6
        do i = 1, Mesh%nPanels
            diffractionForce(j) = diffractionForce(j) + rho * phi(i) * Mesh%UNVA(j, i)
        enddo
    enddo
    open(1, file = 'DiffractionForce.txt')
    do j = 1, 6
        write(1, '(ES14.6, I5, 2(E14.6))') w(1), j, abs(diffractionForce(j)), atan2(&
        &aimag(diffractionForce(j)), real(diffractionForce(j)))
    enddo
    close(1)
END PROGRAM Main