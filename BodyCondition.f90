MODULE MBodyCondition
    USE MMesh, ONLY: TMesh
    USE MConstants, ONLY: czero, ComputeWave
    IMPLICIT NONE
    public :: DiffractionCondition
CONTAINS
    SUBROUTINE DiffractionCondition(w, beta, Mesh, pressure, velocity)
        REAL, INTENT(IN)                                 :: w, beta
        TYPE(TMesh), INTENT(IN)                          :: Mesh
        COMPLEX, DIMENSION(Mesh%nPanels), INTENT(OUT)    :: pressure
        COMPLEX, DIMENSION(Mesh%nPanels), INTENT(oUT)    :: velocity

        ! local variables
        INTEGER :: i
        COMPLEX :: phi
        COMPLEX, DIMENSION(3) :: vector
        
        DO i = 1, Mesh%nPanels
            IF (Mesh%center(3, i) .lt. 0.) THEN 
                CALL ComputeWave(w, beta, Mesh%center(:, i), phi, pressure(i), vector)
            ELSE
                pressure(i) = czero
                vector = (/czero, czero, czero/)
            END IF
            velocity(i) = dot_product(vector, Mesh%normalVector(:, i))
        END DO
    END SUBROUTINE DiffractionCondition

END MODULE MBodyCondition