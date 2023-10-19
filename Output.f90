MODULE MOutput
    USE MMesh,      ONLY: TMesh
    use MConstants, ONLY: nw, nbeta
    IMPLICIT NONE
    PUBLIC :: WriteSummaryMesh, WriteFK, WriteDiffractionVelocity, WriteRadiationVelocity
CONTAINS
    SUBROUTINE WriteSummaryMesh(MeshFile, Mesh)
        ! summary of mesh, in order of Points connectivity, center, area, normalVector, UNV, UNVA
        character(*), intent(in) :: MeshFile
        Type(TMesh), intent(in)  :: Mesh

        !local variables
        INTEGER :: i, j

        open(unit=1, file=MeshFile, action='write')
        write(1, '(A, I5)') 'Number of Points', Mesh%nPoints
        write(1, *) 'Number     x      y      z' 
        DO i = 1, Mesh%nPoints
            write(1, *) i, (Mesh%x(j, i), j=1,3)
        END DO
        write(1, '(A, I5)') 'Number of Panels', Mesh%nPanels
        write(1, '(A)') 'Number   #PT1   #PT2   #PT3   #PT4'
        DO i = 1, Mesh%nPanels
            write(1, *) i, (Mesh%connect(j, i), j=1,4)
        END DO
        write(1, '(A)') 'Center of each panel'
        write(1, '(A)') 'Number   x   y   z'
        DO i = 1, Mesh%nPanels
            write(1, *) i, (Mesh%center(j, i), j=1,3)
        END DO
        write(1, '(A)') 'Area of each panel'
        write(1, '(A)') 'Number   area'
        DO i = 1, Mesh%nPanels
            write(1, *) i, Mesh%area(i)
        END DO
        write(1, '(A)') 'Normal vector of each panel'
        write(1, '(A)') 'Number   x   y   z'
        DO i = 1, Mesh%nPanels
            write(1, *) i, (Mesh%normalVector(j, i), j = 1, 3)
        END DO
        write(1, '(A)') 'UNV of each panel (That is unit velocity in each mode induced normal velocity to each panel)'
        write(1, '(A)') 'Number    mode1    mode2     mode3    mode4    mode5     mode6'
        DO i = 1, Mesh%nPanels
            write(1, *) i, (Mesh%UNV(j, i), j = 1, 6)
        END DO
        write(1, '(A)') 'UNVA of each panel (That is unit velocity in each mode induced normal velocity to&
        & each panel multiply area)'
        write(1, '(A)') 'Number    mode1    mode2     mode3    mode4    mode5     mode6'
        DO i = 1, Mesh%nPanels
            write(1, *) i, (Mesh%UNVA(j, i), j = 1, 6)
        END DO
        write(1, '(A)') 'Gauss quadrature points for Gauss quadrature'
        write(1, '(A)') '   Panel             x             y            z'
        DO i = 1, Mesh%nPanels
            do j=1,4
                write(1, '(I5, I2,3(ES14.6))') i, j, Mesh%GQPoints(:, j, i)
            end do
        END DO
        close(1)
    END SUBROUTINE WriteSummaryMesh

    SUBROUTINE WriteFK(FKFile, w, beta, FK)
        character(*),INTENT(IN) :: FKFile
        REAL   , INTENT(IN), DIMENSION(nw)           :: w
        REAL   , INTENT(IN), DIMENSION(nbeta)        :: beta
        COMPLEX, DIMENSION(nw, nbeta, 6), intent(in) :: FK

        ! local variables
        INTEGER :: i, j, l

        open(unit=1, file=FKFile, action='write')
        write(1, '(2(A14), A6, 4(A14))') 'Freq/(rad/s)', 'waveAngle/rad', &
                &'dof', 'abs(FKForce)', 'angle(FKForce)', 'real(FKForce)', 'imag(FKForce)'
        DO i = 1, nbeta
            DO j = 1, nw
                DO l = 1, 6
                    write(1, '(2(E14.6), I6, 4(ES14.6))') w(j), beta(j), l, abs(FK(j, i, l)),&
                        & atan2(imag(FK(j, i, l)), real(FK(j, i, l))), real(FK(j, i, l)), IMAG(FK(j, i, l))
                enddo
            END DO
        END DO
        close(1)
    END SUBROUTINE WriteFK

    SUBROUTINE WriteDiffractionVelocity(DVFile, w, beta, nPanels, diffractionVelocity)
        character(*), INTENT(IN)                           :: DVFile
        REAL   , INTENT(IN), DIMENSION(nw)                 :: w
        REAL   , INTENT(IN), DIMENSION(nbeta)              :: beta
        INTEGER , INTENT(IN)                               :: nPanels
        COMPLEX, DIMENSION(nw, nbeta, nPanels), INTENT(IN) :: diffractionVelocity

        ! local variables
        INTEGER :: i, j, k

        open(unit=1, file=DVFile, action='write')
        DO i = 1, nbeta
            DO j = 1, nw
                write(1, '(2(A, F5.3))') 'angle :', beta(i), '   Frequency (rad/s) :', w(j)
                write(1, '(A)') 'Number        Real()        IMAG()'
                DO k = 1, nPanels
                    write(1, *) k, REAL(diffractionVelocity(j, i, k)), IMAG(diffractionVelocity(j, i, k))
                END DO
            END DO
        END DO
        close(1)
    END SUBROUTINE WriteDiffractionVelocity

    SUBROUTINE WriteRadiationVelocity(RVFile, w, nPanels, radiationVelocity)
        character(*), INTENT(IN)                           :: RVFile
        REAL   , INTENT(IN), DIMENSION(nw)                 :: w
        INTEGER , INTENT(IN)                               :: nPanels
        COMPLEX, DIMENSION(nw, 6, nPanels), INTENT(IN)     :: radiationVelocity

        ! local variables
        INTEGER :: i, j, k

        open(unit=1, file=RVFile, action='write')
        DO i = 1, nw
            DO j = 1, 6
                write(1, '(A, F5.3, A, I5)') 'Frequency (rad/s) :', w(i), '   Mod:', j
                write(1, '(A)') 'Number        Real()        IMAG()'
                DO k = 1, nPanels
                    write(1, *) k, REAL(radiationVelocity(j, i, k)), IMAG(radiationVelocity(j, i, k))
                END DO
            END DO
        END DO
        close(1)
    END SUBROUTINE WriteRadiationVelocity
END MODULE MOutput