module helper_m

    implicit none

    ! CONSTANTS
    real,parameter :: PI = 3.141592653589793
    real,parameter :: PI2 = PI*0.5
    real,parameter :: PI4 = PI*0.25
    real,parameter :: g_ssl = 9.80665 ! Standard Gravity at sea level 
    real,parameter :: r_ez = 6356766.0 ! US standard atmosphere earth radius in meters
    
contains
    
    !!! COORDINATE TRANSFORM !!!

    function quat_mult(q1, q2) result(q3)

        implicit none
        
        real, intent(in) :: q1(4), q2(4)
        real :: q3(4)

        q3(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
        q3(2) = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3)
        q3(3) = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2)
        q3(4) = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)

    end function quat_mult


    function quat_base_to_dependent(vec1, q1) result(vec2)

        implicit none

        real, intent(in) :: q1(4), vec1(3)
        real :: vec2(3)

        real :: t(4)

        t(1) = -vec1(1)*q1(2) - vec1(2)*q1(3) - vec1(3)*q1(4)
        t(2) = vec1(1)*q1(1) + vec1(2)*q1(4) - vec1(3)*q1(3)
        t(3) = vec1(1)*q1(4) + vec1(2)*q1(1) + vec1(3)*q1(2)
        t(4) = vec1(1)*q1(3) - vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) - q1(2)*t(1) - q1(3)*t(4) + q1(4)*t(3)
        vec2(2) = q1(1)*t(3) + q1(2)*t(4) - q1(3)*t(1) - q1(4)*t(2)
        vec2(3) = q1(1)*t(4) - q1(2)*t(3) + q1(3)*t(2) - q1(4)*t(1)

    end function quat_base_to_dependent


    function quat_dependent_to_base(vec1, q1) result(vec2)

        implicit none
        
        real, intent(in) :: vec1(3), q1(4)
        real :: vec2(3)

        real :: t(4)

        t(1) = vec1(1)*q1(2) + vec1(2)*q1(3) + vec1(3)*q1(4)
        t(2) = vec1(1)*q1(1) - vec1(2)*q1(4) + vec1(3)*q1(3)
        t(3) = vec1(1)*q1(4) + vec1(2)*q1(1) - vec1(3)*q1(2)
        t(4) = -vec1(1)*q1(3) + vec1(2)*q1(2) + vec1(3)*q1(1)

        vec2(1) = q1(1)*t(2) + q1(2)*t(1) + q1(3)*t(4) - q1(4)*t(3)
        vec2(2) = q1(1)*t(3) - q1(2)*t(4) + q1(3)*t(1) + q1(4)*t(2)
        vec2(3) = q1(1)*t(4) + q1(2)*t(3) - q1(3)*t(2) + q1(4)*t(1)
        
    end function quat_dependent_to_base


    subroutine quat_norm(q1)

        implicit none
        
        real, intent(inout) :: q1(4)

        q1 = q1/sqrt(q1(1)**2 + q1(2)**2 + q1(3)**2 + q1(4)**2)

    end subroutine quat_norm


    function euler_to_quat(e1) result(q1)
        implicit none
        
        real, intent(in) :: e1(3)
        real :: q1(4)

        !q1(1) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(2) = sin(e1(1)*0.5)*cos(e1(2)*0.5)*cos(e1(3)*0.5) - cos(e1(1)*0.5)*sin(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(3) = cos(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5) + sin(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5)
        !q1(4) = cos(e1(1)*0.5)*cos(e1(2)*0.5)*sin(e1(3)*0.5) - sin(e1(1)*0.5)*sin(e1(2)*0.5)*cos(e1(3)*0.5)

        real :: s(3), c(3)

        s(1) = sin(e1(1)*0.5)
        s(2) = sin(e1(2)*0.5)
        s(3) = sin(e1(3)*0.5)

        c(1) = cos(e1(1)*0.5)
        c(2) = cos(e1(2)*0.5)
        c(3) = cos(e1(3)*0.5)

        q1(1) = c(1)*c(2)*c(3) + s(1)*s(2)*s(3)    
        q1(2) = s(1)*c(2)*c(3) - c(1)*s(2)*s(3)    
        q1(3) = c(1)*s(2)*c(3) + s(1)*c(2)*s(3)    
        q1(4) = c(1)*c(2)*s(3) - s(1)*s(2)*c(3)    
        
    end function euler_to_quat


    function quat_to_euler(q1) result(e1)
        implicit none
        
        real, intent(in) :: q1(4)
        real :: e1(3)

        if (q1(1)*q1(3) - q1(2)*q1(4) == 0.5) then
            
            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = PI2
            e1(3) = 0.0

        else if (q1(1)*q1(3) - q1(2)*q1(4) == -0.5) then

            e1(1) = 2*asin(q1(2)/cos(PI4))
            e1(2) = -PI2
            e1(3) = 0.0

        else

            e1(1) = atan2(2*(q1(1)*q1(2) + q1(3)*q1(4)), (q1(1)**2 + q1(4)**2 - q1(2)**2 - q1(3)**2))
            e1(2) = asin(2*(q1(1)*q1(3) - q1(2)*q1(4)))
            e1(3) = atan2(2*(q1(1)*q1(4) + q1(2)*q1(3)), (q1(1)**2 + q1(2)**2 - q1(3)**2 - q1(4)**2))

        end if
    
        
    end function quat_to_euler

    !!! STD ATMOSPHERE !!!

    function gravity_SI(h) result(g)
        implicit none
        
        real, intent(in) :: h 
        real :: g

        g = g_ssl*(r_ez/(r_ez + h))**2
        
    end function gravity_SI

    function gravity_English(h_ft) result(g)
        implicit none
        
        real, intent(in) :: h_ft
        real :: g

        g = meters_to_feet(gravity_SI(feet_to_meters(h_ft)))
        
    end function gravity_English


    !!! UNIT CONVERSIONS !!!

    function feet_to_meters(l_ft) result(l_m)

        implicit none
        
        real, intent(in) :: l_ft
        real :: l_m

        l_m = 0.3048*l_ft

    end function feet_to_meters

    function meters_to_feet(l_m) result(l_ft)

        implicit none
        
        real, intent(in) :: l_m
        real :: l_ft

        l_ft = 3.280839895013123*l_m


    end function meters_to_feet



end module helper_m


!program test_gravity

    !use helper_m

    !implicit none

    !integer :: i, j
    !real :: g

    !write(*,*) "SI units"
    !do i = 0, 100000, 5000
        !g = gravity_SI(real(i))

        !write(*,*) "Altitude:", i, " meters. Gravity", g, "m/s2"

    !end do

    !write(*,*)
    !write(*,*) "ENGLISH UNITS"

    !do j = 0, 200000, 10000
        !g = gravity_English(real(j))

        !write(*,*) "Altitude:", j, " meters. Gravity", g, "ft/s2"

    !end do
    
!end program test_gravity