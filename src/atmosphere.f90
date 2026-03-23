module atmosphere_m
    use goates_m
    use jsonx_m
    implicit none

    type atmosphere_t
        real, allocatable :: wind(:)
        character(len=:), allocatable :: turb_model, turb_intensity
        real :: wingspan, hstab_dist, vstab_dist
        logical :: turb_repeatable, turb_on

        real :: light_hag(3) = (/2000., 8000., 17000./)
        real :: light_sig(3) = (/5., 5., 3./)
        real :: moderate_hag(3) = (/2000., 11000., 45000./)
        real :: moderate_sig(3) = (/10., 10., 3./)
        real :: severe_hag(4) = (/2000., 4000., 20000., 80000./)
        real :: severe_sig(4) = (/15., 21., 21., 3./)

        real, allocatable :: turb_hag(:), turb_sig(:)
        real :: prev_xyz(3), prev_f, prev_g
        real, allocatable :: prev_turb(:,:), prev_xff(:)
        real :: Lu, Lv, Lw, Lb
        real :: xff
    end type atmosphere_t
    
contains

    subroutine atmosphere_init(this, j_atmosphere)

        implicit none
        type(atmosphere_t) :: this
        type(json_value), pointer :: j_atmosphere, j_turb, j_sample
        logical :: found
        integer :: i, n, seed, length
        integer, allocatable :: seed_array(:)

        write(*,*) "Initializing Atmosphere Model..."
        call jsonx_get(j_atmosphere, 'constant_wind[ft/s]', this%wind, 0.0, 3)

        write(*,*) "          Constant Wind [ft/s] = ", this%wind
        call json_get(j_atmosphere, 'turbulence', j_turb, found)
        if (found) then
            call jsonx_get(j_turb, 'model', this%turb_model, "none")
            if (this%turb_model /= "none") then
                call jsonx_get(j_turb, 'wingspan[ft]', this%wingspan)
                call jsonx_get(j_turb, 'hstab_distance[ft]', this%hstab_dist)
                call jsonx_get(j_turb, 'vstab_distance[ft]', this%vstab_dist)
                call jsonx_get(j_turb, 'intensity', this%turb_intensity)
                call jsonx_get(j_turb, 'repeatable', this%turb_repeatable)
                !call jsonx_get(j_turb, 'history_length', length)
                !allocate(this%prev_turb(length,3))
                !allocate(this%prev_xff(length))
                !this%prev_turb(:,:) = 0.0
                !do i = 1, length
                    !this%prev_xff(i) = real(i-1)*this%vstab_dist
                !end do

                write(*,*) "          Turbulence Intensity = ", this%turb_intensity
                write(*,*) "          Turbulence Model = ", this%turb_model

                ! set up random number generator
                if (this%turb_repeatable) then
                    seed = 12345
                else
                    call system_clock(count=seed)
                end if
                call random_seed(size=n)
                allocate(seed_array(n))
                seed_array = seed + 7*[(i-1, i=1,n)]
                call random_seed(put=seed_array)
                deallocate(seed_array)

                select case (trim(this%turb_intensity))
                case ('light')
                    allocate(this%turb_hag(3))
                    allocate(this%turb_sig(3))
                    this%turb_hag = this%light_hag
                    this%turb_sig = this%light_sig
                case ('moderate')
                    allocate(this%turb_hag(3))
                    allocate(this%turb_sig(3))
                    this%turb_hag = this%moderate_hag
                    this%turb_sig = this%moderate_sig
                case ('severe')
                    allocate(this%turb_hag(4))
                    allocate(this%turb_sig(4))
                    this%turb_hag = this%severe_hag
                    this%turb_sig = this%severe_sig
                end select

                select case (trim(this%turb_model))
                case('dryden_beal')
                    this%Lu = 1750.0
                    this%Lv = 875.0
                    this%Lw = 875.0
                    this%Lb = 4*this%wingspan/pi
                case('dryden_8785')
                    this%Lu = 1750.0
                    this%Lv = 875.0
                    this%Lw = 875.0
                end select

                !n_hist = max(nint(this%vstab_dist/this%dx), nint(this%hstab_dist/this%dx))
                !if (n_hist < 1) n_hist = 1
                !write(*,*) "          Turbulence history length = ", n_hist
                !allocate(this%prev_turb(n_hist,6))

                !this%prev_turb(:) = 0.0
                this%prev_xyz(:) = 0.0
                this%prev_f = 0.0
                this%prev_g = 0.0
                this%turb_on = .true.

                call json_get(j_turb, 'sample', j_sample, found)
                if (found) call turbulence_sample(this, j_sample)
            end if
        end if

    end subroutine atmosphere_init


    subroutine turbulence_sample(this, j_sample)

        implicit none
        type(atmosphere_t) :: this
        type(json_value), pointer :: j_sample
        character(len=:), allocatable :: fn
        integer :: i, j, n, n_psd, iunit, psd_mean_unit
        real, allocatable :: vals(:,:), psd_mean(:,:), psd_temp(:,:)
        real :: hag, sigma, dx, turb(6)
        real :: mean, stdev, stdev_temp
        logical :: found


        call test_rand_normal()


        write(*,*) "   Sampling Atmospheric Turbulence..."
        call jsonx_get(j_sample, 'save_filename', fn)
        open(newunit=iunit, file=fn, status='replace')
        write(*,*) "     saving sample to ", fn
        call jsonx_get(j_sample, 'number_of_points', n)
        call jsonx_get(j_sample, 'dx[ft]', dx)
        call jsonx_get(j_sample, 'height_above_ground[ft]', hag)
        allocate(vals(n,6))

        sigma = interpolate_1D(this%turb_hag, this%turb_sig, hag)
        write(*,*) "     Altitude = ", hag
        write(*,*) "     Turbulence Standard Deviation = ", sigma

        write(iunit, *) 'distance[ft], uprime[ft/s], vprime[ft/s], wprime[ft/s], pprime[rad/s], qprime[rad/s], rprime[rad/s]'
        do i = 1, n
            turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
            write(iunit, *) dx*real(i-1), ',', turb(1), ',', turb(2), ',', turb(3), ',', turb(4), ',', turb(5), ',', turb(6)
            vals(i,:) = turb(:)
        end do
        close(iunit)


        call json_get(j_sample, 'psd_analyses', n_psd, found)
        if (found) then
            call jsonx_get(j_sample, 'psd_analyses', n_psd)
            allocate(psd_mean(n/2+1,2))
            allocate(psd_temp(n/2+1,2))
            psd_mean(:,:) = 0.0
            psd_temp(:,:) = 0.0
            stdev = 0.0
            stdev_temp = 0.0

            write(*,*) 'Turbulence PSD Analyses:'
            open(newunit=psd_mean_unit, file='psd_mean.csv', status='replace')
            write(*,*) '    saving mean PSD to psd_mean.csv'
            do j = 1, n_psd
                write(*,*) 'u PSD ', j, ' of ', n_psd
                do i = 1, n
                    turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
                    vals(i,:) = turb(:)
                end do
                call psd(vals(:,6), dx, psd_norm=psd_temp, stdev=stdev_temp) ! 1 for u, 2 for v, 3 for w
                stdev = stdev + stdev_temp/n_psd
                psd_mean(:,2) = psd_mean(:,2) + psd_temp(:,2)/n_psd
            end do
            
            do i = 1, n/2+1
                write(psd_mean_unit,*) psd_temp(i,1)*2*pi, ',', psd_mean(i,2)/2/pi
            end do
            write(*,*) "     Mean Standard Deviation = ", stdev
            close(psd_mean_unit)

        end if

        deallocate(vals)


    end subroutine turbulence_sample

    function atmosphere_get_turbulence(this, states) result(ans)
        implicit none
        
        type(atmosphere_t), intent(inout) :: this
        real, intent(in) :: states(21)
        real :: ans(6)
        real :: dx, sigma

        dx = sqrt((states(7)-this%prev_xyz(1))**2 + (states(8) - this%prev_xyz(2))**2 + (states(9) - this%prev_xyz(3))**2)
        sigma = interpolate_1D(this%turb_hag, this%turb_sig, -states(9))

        if (this%turb_on .and. dx > 1.e-12) then
            ans(:) = get_turbulence(this, dx, sigma, sigma, sigma)
        else
            ans(:) = 0.0
        end if
        
        this%prev_xyz(:) = states(7:9)

    end function atmosphere_get_turbulence

    function get_turbulence(this, dx, sigma_u, sigma_v, sigma_w) result(ans)

        implicit none
        type(atmosphere_t) :: this
        real :: dx, sigma_u, sigma_v, sigma_w
        real :: ans(6)

        select case (trim(this%turb_model))
        case('dryden_beal')
            ans(:) = dryden_beal(this, dx, sigma_u, sigma_v, sigma_w)
        case('dryden_8785')
            continue
            !ans(:) = dryden_8785(this, dx, sigma_u, sigma_v, sigma_w)
        end select

    end function get_turbulence

    function dryden_beal(this, dx, sigma_u, sigma_v, sigma_w) result(ans)

        implicit none
        type(atmosphere_t) :: this
        real, intent(in):: dx, sigma_u, sigma_v, sigma_w
        real :: sigma_np
        real :: ans(6)
        real :: xff
        real :: Au, Av, Aw, Ap, f, g, etau, etav, etaw
        integer :: n,i

        if (.not. allocated(this%prev_turb)) then
            n = max(nint(this%vstab_dist/dx), nint(this%hstab_dist/dx))
            if (n < 1) n = 1
            write(*,*) "          Turbulence history length = ", n*10
            write(*,*) "                                dx  = ", dx
            allocate(this%prev_turb(n*10,3))
            allocate(this%prev_xff(n*10))
            this%prev_turb(:,:) = 0.0
            do i = 1, n*10
                this%prev_xff(i) = real(i-1)*this%vstab_dist
            end do
        end if

        Au = 0.5*dx/this%Lu
        Av = 0.25*dx/this%Lv
        Aw = 0.25*dx/this%Lw

        etau = rand_normal()*sigma_u*sqrt(2.0*this%Lu/dx)
        etav = rand_normal()*sigma_v*sqrt(2.0*this%Lv/dx)
        etaw = rand_normal()*sigma_w*sqrt(2.0*this%Lw/dx)

        f = ((1.0-Av)*this%prev_f + 2.0*Av*etav)/(1.0+Av)
        g = ((1.0-Aw)*this%prev_g + 2.0*Aw*etaw)/(1.0+Aw)

        ans(1) = ((1.0-Au)*this%prev_turb(1,1) + 2.0*Au*etau)/(1.0+Au)
        ans(2) = ((1.0-Av)*this%prev_turb(1,2) + Av*(f + this%prev_f) + sqrt(3.0)*(f-this%prev_f))/(1.0+Av)
        ans(3) = ((1.0-Aw)*this%prev_turb(1,3) + Aw*(g + this%prev_g) + sqrt(3.0)*(g-this%prev_g))/(1.0+Aw)

        sigma_np = sigma_w*sqrt(0.8*pi*(this%Lw/this%Lb)**(1.0/3.0)/(this%Lw*dx))
        Ap = 0.5*dx/this%Lb
        ans(4) = ((1.0-Ap)*this%prev_turb(-1,4) + 2.0*Ap*rand_normal()*sigma_np)/(1.0+Ap)

        this%prev_turb(1,:) = ans(1:3)

        ans(5) =  (ans(3)-interpolate_1D( this%prev_xff, this%prev_turb(:,3),this%hstab_dist))/this%hstab_dist
        ans(6) = -(ans(2)-interpolate_1D( this%prev_xff, this%prev_turb(:,2),this%vstab_dist))/this%vstab_dist

        ! shift history and save new values
        this%prev_f = f
        this%prev_g = g
        do i = size(this%prev_turb,1), 2, -1
            this%prev_turb(i,:) = this%prev_turb(i-1,:)
            this%prev_xff(i) = this%prev_xff(i-1)
        end do
        !this%prev_turb(2:,:) = this%prev_turb(1:-1,:)
        !this%prev_xff(2:) = this%prev_xff(1:-1)
        this%prev_xff = this%prev_xff + dx
        this%prev_xff(1) = 0.0
        this%prev_turb(1,:) = ans(1:3)
    end function dryden_beal



    
end module atmosphere_m