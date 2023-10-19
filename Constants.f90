MODULE MConstants
    IMPLICIT NONE
    real    :: rho = 1025.0 ! density of fluid
    real    :: g = 9.81     ! gravity acceleration
    real    :: depth = -1.   ! depth of fluid
    complex :: ii = Cmplx(0., 1.)
    real    :: pi = acos(-1.0)
    INTEGER :: nw = 1     ! number of frequencies
    real    :: wmin = 1.0     ! minimum frequency (rad/s)
    real    :: wmax = 1.0    ! maximum frequency (rad/s)
    real    :: xref = 0     ! reference position of x
    real    :: yref = 0     ! reference position of y
    INTEGER :: nbeta = 1    ! number of beta, the angle of input wave
    real    :: betamin = 0  ! angle of direction of input wave
    real    :: betamax = 0  ! angle with positive of x axis
    COMPLEX :: czero = COMPLEX(0., 0.) !complex zero

    PUBLIC  :: ComputeWave, CrossProduct, Average, Norm, Omega2k
CONTAINS
    SUBROUTINE ComputeWave(w, beta, location, phi, pressure, velocity)
        ! Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)
        REAL,intent(in) :: w,beta,location(3)
        COMPLEX,intent(out) :: phi,pressure,velocity(3)
        
        ! local variable
        REAL :: k
        REAL :: wbar
        k = Omega2k(w)

        wbar=(location(1)-xref)*COS(Beta)+(location(2)-yref)*SIN(Beta)
        phi=-ii*g/w*CIH(k,location(3))*CEXP(II*k*wbar)
        pressure=rho*g*CIH(k,location(3))*CEXP(II*k*wbar)
        velocity(1)=g/w*k*COS(beta)*CIH(k,location(3))*CEXP(II*k*wbar)
        velocity(2)=g/w*k*SIN(beta)*CIH(k,location(3))*CEXP(II*k*wbar)
        velocity(3)=-ii*g/w*k*SIH(k,location(3))*CEXP(II*k*wbar)
      END SUBROUTINE ComputeWave
    
    REAL FUNCTION CIH(k, z)
        real, intent(in)  :: k, z
        if (depth < 0) Then 
            CIH = EXP(k * z)
        ELSE 
            CIH = COSH(k * (z + depth)) / COSH(k * depth)
        endif
    END FUNCTION CIH

    REAL FUNCTION SIH(k, z)
        real, intent(in)  :: k, z
        if (depth < 0) Then 
            SIH = exp(k * z)
        else 
            SIH = SINH(k * (z + depth)) / COSH(k * depth)
        endif
    END FUNCTION SIH

    REAL Function Omega2k(omega) result(k)
        real, intent(in) :: omega
        ! local variables
        REAL             :: k0, k1, diff, eps
        if (depth < 0.) Then 
            k = omega**2 / g
            return
        endif   
        k0 = 0;
        k1 = 0.1;
        DO while (omega**2 - g * k1 * tanh(k1 * depth) > 0)
            k1 = k1 + 0.1;
        END DO
        
        eps = 2e-5
        ! bisection to get root
        diff = 1.
        DO while (abs(diff) > eps)
            
            k = (k0 + k1) / 2
            diff = omega**2 - g * k * tanh(k * depth)
            IF (diff > 0) THEN
                k0 = k
            ELSE
                k1 = k
            END IF
        END DO
    END FUNCTION Omega2k

    SUBROUTINE CrossProduct(var1, var2, result)
        REAL, INTENT(IN), DIMENSION(3) :: var1, var2
        REAL, INTENT(OUT),DIMENSION(3) :: result

        result(1) = var1(2)*var2(3) - var1(3)*var2(2)
        result(2) = var1(3)*var2(1) - var1(1)*var2(3)
        result(3) = var1(1)*var2(2) - var1(2)*var2(1)
    END SUBROUTINE CrossProduct

    SUBROUTINE Average(var1, var2, var3, result)
        REAL, INTENT(IN), DIMENSION(3) :: var1, var2, var3
        REAL, INTENT(OUT),DIMENSION(3) :: result
        
        result = (var1 + var2 + var3)/3
    END SUBROUTINE Average

    REAL FUNCTION Norm(var)
        REAL, INTENT(IN), DIMENSION(3) :: var
        
        Norm = sqrt((var(1)**2 + var(2)**2 + var(3)**2))
    END FUNCTION Norm
END MODULE MConstants