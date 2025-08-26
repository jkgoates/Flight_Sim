module helper_m
    implicit none
    real,parameter :: PI = 3.141592653589793
    real,parameter :: PI2 = PI*0.5
    real,parameter :: PI4 = PI*0.25
    
contains
    
    function quat_mult(q1, q2) result(q3)

        implicit none
        
        real, intent(in) :: q1(4), q2(4)
        real :: q3(4)

        ! 1.31250000 secs
        q3(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
        q3(2) = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3)
        q3(3) = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2)
        q3(4) = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)

        ! 1.65625000 secs
        !q3 = (/q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4),&
                !q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3),&
                !q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2), &
                !q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)/)

    end function quat_mult


    function quat_base_to_dependent(vec1, q1) result(vec2)

        implicit none

        real, intent(in) :: q1(4), vec1(3)
        real :: vec2(3)

        ! explicit cross mult approach FASTEST
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






end module helper_m