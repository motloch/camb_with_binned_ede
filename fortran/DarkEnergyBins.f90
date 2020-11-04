    !<pavel>
    !DE model inspired by 1304.3724
    module DarkEnergyBins
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none


    type, extends(TDarkEnergyModel) :: TDarkEnergyBins
    contains
    procedure :: ReadParams =>  TDarkEnergyBins_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyBins_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyBins_SelfPointer
    procedure :: Init => TDarkEnergyBins_Init
    procedure :: w_de => TDarkEnergyBins_w_de
    procedure :: grho_de => TDarkEnergyBins_grho_de
    procedure :: PerturbedStressEnergy => TDarkEnergyBins_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyBins_PerturbationEvolve
    procedure :: BackgroundDensityAndPressure => TDarkEnergyBins_BackgroundDensityAndPressure
    end type TDarkEnergyBins
    !These are missing:
    !procedure :: diff_rhopi_Add_Term     -    this means there is no DE quadrupole
    !procedure :: Effective_w_wa !Used as approximate values for non-linear corrections
    !                   - we can keep this at the LCDM values

    contains

    subroutine TDarkEnergyBins_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyBins) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i

    call this%TDarkEnergyModel%ReadParams(Ini)
    !Read the number of DE bins and allocate the arrays
    this%de_n_bins  = Ini%Read_Double('DE_n_bins')
    this%de_step_type = Ini%Read_Int('DE_step_type', 1)
    this%de_use_perturbations = Ini%Read_Int('DE_use_perturbations', 1)
    this%de_overflow_cutoff = Ini%Read_Double('DE_OVERFLOW_CUTOFF', 7.d2)
    allocate(this%de_bin_ai(this%de_n_bins+1))
    allocate(this%de_bin_amplitudes(this%de_n_bins))
    !Read the DE bin boundaries and amplitudes in the bin
    do i = 1, this%de_n_bins+1
        this%de_bin_ai(i)  = Ini%Read_Double_Array('DE_bin_ai', i)
    enddo
    do i = 1, this%de_n_bins
        this%de_bin_amplitudes(i)  = Ini%Read_Double_Array('DE_bin_amplitudes', i)
    enddo
    !Read the smoothing scale and soundspeed
    this%de_tau  = Ini%Read_Double('DE_tau')
    this%de_cs2 = Ini%Read_Double('DE_cs2')

    end subroutine TDarkEnergyBins_ReadParams


    function TDarkEnergyBins_PythonClass()
    character(LEN=:), allocatable :: TDarkEnergyBins_PythonClass

    TDarkEnergyBins_PythonClass = 'DarkEnergyBins'
    end function TDarkEnergyBins_PythonClass

    subroutine TDarkEnergyBins_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyBins), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyBins_SelfPointer

    subroutine TDarkEnergyBins_Init(this, State)
    use classes
    class(TDarkEnergyBins), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State
    real(dl) :: grho_rad, F, p, mu, xc, n

    select type(State)
    class is (CAMBdata)
        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 2
    end select

    end subroutine TDarkEnergyBins_Init

    subroutine TDarkEnergyBins_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyBins), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TDarkEnergyBins_PerturbedStressEnergy

    !Top function (choice from two)
    function SmoothedTopFunction(a, ai, aip1, tau, step_type, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    integer, intent(in) :: step_type
    real(dl) :: SmoothedTopFunction
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedTopFunction = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedTopFunction = 0
        else
            !We are in the step region
            if(step_type == 1) then
                exparg = exp(arg)
                expargp1 = exp(argp1)
                SmoothedTopFunction = (exparg - expargp1)/(1. + expargp1)/(1. + exparg)
            else if(step_type == 2) then
                SmoothedTopFunction = (erf(arg/sqrt(2.)) - erf(argp1/sqrt(2.)))/2.
            else
               stop 'step'
            endif
        endif
    end function

    !Derivative of step function (choice from two) wrt ln a
    function SmoothedTopFunctionDer(a, ai, aip1, tau, step_type, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    integer, intent(in) :: step_type
    real(dl) :: SmoothedTopFunctionDer
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedTopFunctionDer = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedTopFunctionDer = 0
        else
            !We are in the step region
            if(step_type == 1) then
                exparg = exp(arg)
                expargp1 = exp(argp1)
                SmoothedTopFunctionDer = -expargp1/(1. + expargp1)**2/tau + &
                    exparg/(1. + exparg)**2/tau
            else if(step_type == 2) then
                SmoothedTopFunctionDer = (exp(-arg**2/2.) - exp(-argp1**2/2.))/sqrt(const_twopi)/tau
            else
                stop 'step'
            endif
        endif

    end function

    !Second derivative of step function (choice from two) wrt ln a
    function SmoothedTopFunctionDerDer(a, ai, aip1, tau, step_type, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    integer, intent(in) :: step_type
    real(dl) :: SmoothedTopFunctionDerDer
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedTopFunctionDerDer = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedTopFunctionDerDer = 0
        else
            !We are in the step region
            if(step_type == 1) then
                exparg = exp(arg)
                expargp1 = exp(argp1)
                SmoothedTopFunctionDerDer = 2*expargp1**2/(1. + expargp1)**3/tau**2&
                    - expargp1/(1. + expargp1)**2/tau**2 &
                    -2*exparg**2/(1. + exparg)**3/tau**2 &
                    + exparg/(1. + exparg)**2/tau**2
            else if(step_type == 2) then
                SmoothedTopFunctionDerDer = (-arg*exp(-arg**2/2.) + argp1*exp(-argp1**2/2.))/sqrt(const_twopi)/tau**2
            else
                stop 'step'
            endif
        endif
    end function

    !Step function
    function SmoothedStepFunction(a, ai, aip1, tau, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    real(dl) :: SmoothedStepFunction
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedStepFunction = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedStepFunction = 0
        else
            exparg = exp(arg)
            expargp1 = exp(argp1)
            SmoothedStepFunction = (exparg - expargp1)/(1. + expargp1)/(1. + exparg)
        endif
    end function

    !Derivative of step function wrt ln a
    function SmoothedStepFunctionDer(a, ai, aip1, tau, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    real(dl) :: SmoothedStepFunctionDer
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedStepFunctionDer = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedStepFunctionDer = 0
        else
            exparg = exp(arg)
            expargp1 = exp(argp1)
            SmoothedStepFunctionDer = -expargp1/(1. + expargp1)**2/tau + &
                exparg/(1. + exparg)**2/tau
        endif

    end function

    !Second derivative of step function wrt ln a
    function SmoothedStepFunctionDerDer(a, ai, aip1, tau, cutoff)
    real(dl), intent(in) :: a, ai, aip1, tau, cutoff
    real(dl) :: SmoothedStepFunctionDerDer
    real(dl) :: arg, argp1
    real(dl) :: exparg, expargp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -cutoff .and. argp1 < -cutoff) then
            SmoothedStepFunctionDerDer = 0
        else if(arg > cutoff .and. argp1 > cutoff) then
            SmoothedStepFunctionDerDer = 0
        else
            exparg = exp(arg)
            expargp1 = exp(argp1)
            SmoothedStepFunctionDerDer = 2*expargp1**2/(1. + expargp1)**3/tau**2&
                - expargp1/(1. + expargp1)**2/tau**2 &
                -2*exparg**2/(1. + exparg)**3/tau**2 &
                + exparg/(1. + exparg)**2/tau**2
        endif
    end function

    subroutine TDarkEnergyBins_BackgroundDensityAndPressure(this, grhov, a, grhov_t, grhoa2_noDE, w, w_bg)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyBins), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(in) :: grhoa2_noDE
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl), optional, intent(in) :: w_bg
    real(dl) :: grhov_t_lcdm, grho_t, grhov_t_beyond
    real(dl) :: arg1, arg2
    real(dl) :: delta, Q
    integer :: i

    grhov_t_lcdm = grhov * a * a
    delta = 0

    if(a > 0) then
        grho_t = grhoa2_noDE/a/a + grhov_t_lcdm

        !Get contributions from each of the dark energy bins
        do i = 1, this%de_n_bins

            delta = delta + this%de_bin_amplitudes(i)*SmoothedStepFunction(a, &
                                this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau, &
                                this%de_overflow_cutoff)
        enddo

        grhov_t_beyond = delta*grho_t
        grhov_t = grhov_t_lcdm + grhov_t_beyond

        !Factor that goes into the DE equation of state
        Q = grhoa2_noDE/a/a/grhov_t_lcdm
        if (present(w)) then
            if (.not. present(w_bg)) then
                stop 'w_bg'
            else
                w = this%w_de(a, delta, Q, w_bg)
            endif
        endif

    !At a = 0, just use LCDM
    else
        grho_t = grhov_t_lcdm
        grhov_t_beyond = 0
        grhov_t = grhov_t_lcdm
        if (present(w)) then
            w = -1
        endif
    endif

    end subroutine TDarkEnergyBins_BackgroundDensityAndPressure

    function TDarkEnergyBins_grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyBins) :: this
    real(dl) :: TDarkEnergyBins_grho_de, apow
    real(dl), intent(IN) :: a

    !If my reading is correct, by redefining BackgroundDensity we no longer need
    !this function
    stop 'grho'

    end function TDarkEnergyBins_grho_de

    function TDarkEnergyBins_w_de(this, a, delta, Q, w_bg)
    class(TDarkEnergyBins) :: this
    real(dl) :: TDarkEnergyBins_w_de
    real(dl), intent(IN) :: a
    real(dl), intent(in) :: delta, Q, w_bg
    real(dl) :: ddelta_dlna
    integer :: i

    ddelta_dlna = 0

    !Get contributions from each of the dark energy bins
    do i = 1, this%de_n_bins

        ddelta_dlna = ddelta_dlna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDer(a, &
                        this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau, &
                        this%de_overflow_cutoff)
    enddo

    !Eq 3 from 1304.3724
    TDarkEnergyBins_w_de = -1 + Q*delta/(1+delta+delta*Q)*(1+w_bg)&
        -(1+Q)/3./(1+delta+delta*Q)*ddelta_dlna

    end function TDarkEnergyBins_w_de

    subroutine TDarkEnergyBins_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y, Q, Q_dot, w_bg, w_bg_dot)
    class(TDarkEnergyBins), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    real(dl), intent(in) :: Q, Q_dot, w_bg, w_bg_dot
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv
    real(dl) delta, ddelta_dlna, d2delta_d2lna
    real(dl) :: dQ_dlna, dw_bg_dlna
    real(dl) :: denom, t1, t2, t3, t4
    real(dl) :: w_ratio
    integer :: i

    if(this%de_use_perturbations == 0) then
        ayprime(w_ix) = 0
        ayprime(w_ix + 1) = 0
        return
    endif

    dQ_dlna = Q_dot/adotoa
    dw_bg_dlna = w_bg_dot/adotoa

    delta = 0
    ddelta_dlna = 0
    d2delta_d2lna = 0

    !Get contributions from each of the dark energy bins
    do i = 1, this%de_n_bins

        delta = delta + this%de_bin_amplitudes(i)*SmoothedStepFunction(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau, &
            this%de_overflow_cutoff)

        ddelta_dlna = ddelta_dlna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDer(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau, &
            this%de_overflow_cutoff)

        d2delta_d2lna = d2delta_d2lna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDerDer(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau, &
            this%de_overflow_cutoff)

    enddo

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    ! dw/dlog a/(1+w) - derivative of Eq 3 from 1304.3724

    !We expect there is no phantom crossing
    if(1+w < 0) then
        stop 'w'
    !Allows us to run LCDM with the same code
    else if(1+w .eq. 0) then
        deriv = 0
    else
        denom = 1 + delta + delta*Q
        t1 = -(1+Q)/3.*d2delta_d2lna/denom
        t2 = delta*Q*dw_bg_dlna/denom
        t3 = ddelta_dlna/3./denom**2*(ddelta_dlna*(1+Q)**2 + 3*Q*(1+w_bg))
        t4 = -dQ_dlna/3./denom**2*(ddelta_dlna - 3*delta*(1+delta)*(1+w_bg))

        deriv  = (t1 + t2 + t3 + t4)/(1+w)
    endif

    !TEST 6
    !if(abs(k - 0.524) < 0.002) then
    !    write(*,'(48e19.6)') a, deriv
    !endif
    !</pavel>

    !density perturbation
    !Looks like there is a typo in 1806.10608: in eq 22 they have [delta] =
    ![theta/k^2], in eq 23 they have [delta] = [theta H / k^2]
    ayprime(w_ix) = -3 * adotoa * (this%de_cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ! The deriv term is from w' in [(1+w)v]' = w'v + (1+w)v'
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%de_cs2 - deriv) * y(w_ix + 1) + &
        k * this%de_cs2 * y(w_ix)

    end subroutine TDarkEnergyBins_PerturbationEvolve

    end module DarkEnergyBins
    !</pavel>
