MODULE MMesh
    USE MConstants, ONLY:Norm, CrossProduct, Average, Norm
    IMPLICIT NONE
    Type TMesh
        INTEGER :: nPoints
        INTEGER :: nPanels
        REAL,   DIMENSION(:, :),allocatable :: x       ! Location of each point (3 x NPoints)
        INTEGER,DIMENSION(:, :),allocatable :: connect ! connectivity of each panel (4 x NPanel) quadrilateral
        REAL,   DIMENSION(:, :),allocatable :: center  ! location of each panel (3 x NPanels)
        REAL,   DIMENSION(:),   allocatable :: area    ! area of each panel (NPanels)
        REAL,   DIMENSION(:, :),allocatable :: normalVector ! normal vector of each panel (3 x NPanels)
        REAL,   DIMENSION(:, :),allocatable :: UNV     ! unit velocity in each mode induced normal velocity (6 x NPanels)
        REAL,   DIMENSION(:, :),allocatable :: UNVA    ! unit velocity in each mode induced normal velocity multiply area (6 x NPanels)
        REAL, dimension(:, :, :),allocatable:: GQPoints! points for Gauss quadrature, note that we choose 4 points in a plane
        REAL,   dimension(:, :),allocatable :: GQWeights! weights for Gauss quadrature, note that there are four weghts
    END TYPE TMesh
    TYPE TDATGQ
        REAL, DIMENSION(4):: XGQ
        REAL, DIMENSION(4):: YGQ
        REAL, DIMENSION(4):: WGQ
    END TYPE

    PUBLIC  :: GenerateMesh ! generate mesh for this specific program
CONTAINS
    SUBROUTINE GenerateMesh(pointsFile, connectFile, motionFile, Mesh)
        character(*),INTENT(IN)  :: pointsFile, connectFile, motionFile

        TYPE(TMesh), INTENT(OUT) :: Mesh
        
        ! local variables
        INTEGER :: i, id, j, l
        REAL,DIMENSION(3   ) :: pt1, pt2, pt3, pt4, vec1, vec2, w1, w2, Tangen1, Normal, Tangen2
        REAL,DIMENSION(4   ) :: xl, yl
        REAL,DIMENSION(3, 3) :: translation
        REAL,DIMENSION(3, 6) :: rotation
        REAL                 :: area1, area2, AA, BB, CC, DD, A, B, C, D, xg1, yg1
        type(TDATGQ)         :: DataGQ
        
        ! read points
        open(unit=1, file=pointsFile, action='read', status='old')
        read(1, *) Mesh%nPoints
        allocate(Mesh%x(3, Mesh%nPoints))
        DO i = 1, Mesh%nPoints
            read(1, *) id, (Mesh%x(j, i), j= 1, 3)
        END DO
        close(1)

        ! read connect
        open(unit=2, file=connectFile, action='read', status='old')
        read(2, *) Mesh%nPanels
        allocate(Mesh%connect(4, Mesh%nPanels), Mesh%center(3, Mesh%nPanels),&
                Mesh%area(Mesh%nPanels), Mesh%normalVector(3, Mesh%nPanels), &
                Mesh%UNV(6, Mesh%nPanels), Mesh%UNVA(6, Mesh%nPanels), &
                &Mesh%GQPoints(3, 4, Mesh%nPanels), Mesh%GQWeights(4, Mesh%nPanels))
        DO i = 1, Mesh%nPanels
            read(2, *) (Mesh%connect(j, i), j= 1, 4)
        END DO
        close(2)

        ! read motion
        open(unit=3, file=motionFile, action='read', status='old')
        DO i = 1, 3
            read(3, *) (translation(i, j), j= 1, 3)
        END DO
        DO i = 1, 3
            read(3, *) (rotation(i, j), j= 1, 6)
        END DO
        close(3)

        CALL GAUSS_QUADRATURE_DATA(DataGQ)
        ! perform extra calculation and prepare Gauss Quaduature
        DO i = 1, Mesh%nPanels
            pt1 = Mesh%x(:, Mesh%connect(1, i))
            pt2 = Mesh%x(:, Mesh%connect(2, i))
            pt3 = Mesh%x(:, Mesh%connect(3, i))
            pt4 = Mesh%x(:, Mesh%connect(4, i))

            ! Area of 1-2-4 triangle
            vec1 = pt2-pt1
            vec2 = pt4-pt2
            CALL CrossProduct(vec1, vec2, w1)
            area1 = 0.5 * Norm(w1)

            ! Area of 2-3-4 triangle
            vec1 = pt4-pt3
            vec2 = pt2-pt3
            CALL CrossProduct(vec1, vec2, w2)
            area2 = 0.5 * Norm(w2)

            Mesh%area(i) = area1 + area2

            Mesh%normalVector(:, i) = (w1 + w2) / Norm(w1 + w2)

            CALL Average(pt1, pt2, pt4, w1)
            CALL Average(pt2, pt3, pt4, w2)
            Mesh%center(:, i) = (area1 * w1 + area2 * w2) / (area1 + area2)

            DO j = 1, 3
                Mesh%UNV(j, i) = dot_product(Mesh%normalVector(:, i), translation(j, :))
                Mesh%UNVA(j, i) = dot_product(Mesh%normalVector(:, i), translation(j, :)) * Mesh%area(i)
            END DO
            DO j = 1, 3
                CALL CrossProduct(rotation(j, 1:3), Mesh%center(:, i) - rotation(j, 4:6), w1)
                Mesh%UNV(3+j, i) = dot_product(Mesh%normalVector(:, i), w1)
                Mesh%UNVA(3+j, i) = dot_product(Mesh%normalVector(:, i), w1) * Mesh%area(i)
            END DO
            
            ! prepare Gauss Quaduature
            Tangen1=pt3-pt1
            Normal=Mesh%normalVector(:, i) !Unit normal vector had been calculated before

            Tangen1=Tangen1(1:3)/SQRT(Tangen1(1)**2+Tangen1(2)**2+Tangen1(3)**2) !Unit Tangen
            call CrossProduct(Normal,Tangen1, Tangen2)

            XL(1)=dot_product(Tangen1, (pt1-Mesh%center(:, i)))
            XL(2)=dot_product(Tangen1, (pt2-Mesh%center(:, i)))
            XL(3)=dot_product(Tangen1, (pt3-Mesh%center(:, i)))
            XL(4)=dot_product(Tangen1, (pt4-Mesh%center(:, i)))
            YL(1)=-dot_product(Tangen2, (pt1-Mesh%center(:, i)))
            YL(2)=-dot_product(Tangen2, (pt2-Mesh%center(:, i)))
            YL(3)=-dot_product(Tangen2, (pt3-Mesh%center(:, i)))
            YL(4)=-dot_product(Tangen2, (pt4-Mesh%center(:, i)))
            DO L=1,4
                AA=.25*(1-DATAGQ%XGQ(L))*(1-DATAGQ%YGQ(L))
                BB=.25*(1-DATAGQ%XGQ(L))*(1+DATAGQ%YGQ(L))
                CC=.25*(1+DATAGQ%XGQ(L))*(1+DATAGQ%YGQ(L))
                DD=.25*(1+DATAGQ%XGQ(L))*(1-DATAGQ%YGQ(L))
        
                XG1=AA*XL(1)+BB*XL(2)+CC*XL(3)+DD*XL(4)
                YG1=AA*YL(1)+BB*YL(2)+CC*YL(3)+DD*YL(4)
        
                Mesh%GQPoints(1:3,L, i)=Mesh%center(1:3, i)+Tangen1(1:3)*XG1-Tangen2(1:3)*YG1
        
                A=(1.-DATAGQ%YGQ(L))*(XL(4)-XL(1))+(1+DATAGQ%YGQ(L))*(XL(3)-XL(2))
                B=(1.-DATAGQ%XGQ(L))*(YL(2)-YL(1))+(1+DATAGQ%XGQ(L))*(YL(3)-YL(4))
                C=(1.-DATAGQ%XGQ(L))*(XL(2)-XL(1))+(1+DATAGQ%XGQ(L))*(XL(3)-XL(4))
                D=(1.-DATAGQ%YGQ(L))*(YL(4)-YL(1))+(1+DATAGQ%YGQ(L))*(YL(3)-YL(2))
                Mesh%GQWeights(L, i)=(ABS(A*B-C*D)*DATAGQ%WGQ(L)*.25)
            ENDDO
        END DO
    END SUBROUTINE GenerateMesh


    SUBROUTINE GAUSS_QUADRATURE_DATA(DATAGQ)
        !This subroutine provide Gauss quadrature nodes and weigthing function data
      
        !INPUT/OUTPUT
        TYPE(TDATGQ),         INTENT(OUT):: DATAGQ
        !Local
        REAL,DIMENSION(4)  :: XGQ,YGQ,WGQ
      
              DATA XGQ /.57735027, .57735027,-.57735027,-.57735027/
              DATA YGQ /-.57735027,.57735027,-.57735027,.57735027/
              DATA WGQ /.25,.25,.25,.25/
      
        DATAGQ%XGQ=XGQ(:)
        DATAGQ%YGQ=YGQ(:)
        DATAGQ%WGQ=WGQ(:)
      
    END SUBROUTINE
END MODULE MMesh