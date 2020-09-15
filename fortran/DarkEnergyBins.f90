    !<pavel>
    module DarkEnergyBins
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none


    type, extends(TDarkEnergyModel) :: TDarkEnergyBins
        integer :: de_n_bins
        real(dl) :: de_cs2
        real(dl) :: de_tau
        real(dl), allocatable :: de_bin_ai(:)
        real(dl), allocatable :: de_bin_amplitudes(:)
    contains
    procedure :: ReadParams =>  TDarkEnergyBins_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyBins_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyBins_SelfPointer
    procedure :: Init => TDarkEnergyBins_Init
    procedure :: w_de => TDarkEnergyBins_w_de
    procedure :: grho_de => TDarkEnergyBins_grho_de
    procedure :: PerturbedStressEnergy => TDarkEnergyBins_PerturbedStressEnergy
    !procedure :: diff_rhopi_Add_Term => TDarkEnergyBins_diff_rhopi_Add_Term
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

    !Step function from 1/[1 + exp{ln(a/ai)/tau}]
    function SmoothedStepFunction(a, ai, aip1, tau)
    real(dl), intent(in) :: a, ai, aip1, tau
    real(dl) :: SmoothedStepFunction
    real(dl) :: arg, argp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -700 .and. argp1 < -700) then
            SmoothedStepFunction = 0
        else if(arg > 700 .and. argp1 > 700) then
            SmoothedStepFunction = 0
        else
            SmoothedStepFunction = (exp(arg) - exp(argp1))/(1. + exp(argp1))/(1. + exp(arg))
            SmoothedStepFunction = (erf(arg/sqrt(2.)) - erf(argp1/sqrt(2.)))/2.
        endif
    end function

    !Derivative of step function from 1/[1 + exp{ln(a/ai)/tau}] wrt ln a
    function SmoothedStepFunctionDer(a, ai, aip1, tau)
    real(dl), intent(in) :: a, ai, aip1, tau
    real(dl) :: SmoothedStepFunctionDer
    real(dl) :: arg, argp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -700 .and. argp1 < -700) then
            SmoothedStepFunctionDer = 0
        else if(arg > 700 .and. argp1 > 700) then
            SmoothedStepFunctionDer = 0
        else
            SmoothedStepFunctionDer = -exp(argp1)/(1. + exp(argp1))**2/tau + &
                exp(arg)/(1. + exp(arg))**2/tau
            SmoothedStepFunctionDer = (exp(-arg**2/2.) - exp(-argp1**2/2.))/sqrt(const_twopi)/tau
        endif

    end function

    !Second derivative of step function from 1/[1 + exp{ln(a/ai)/tau}] wrt ln a
    function SmoothedStepFunctionDerDer(a, ai, aip1, tau)
    real(dl), intent(in) :: a, ai, aip1, tau
    real(dl) :: SmoothedStepFunctionDerDer
    real(dl) :: arg, argp1

        arg = log(a/ai)/tau
        argp1 = log(a/aip1)/tau

        if(arg < -700 .and. argp1 < -700) then
            SmoothedStepFunctionDerDer = 0
        else if(arg > 700 .and. argp1 > 700) then
            SmoothedStepFunctionDerDer = 0
        else
            SmoothedStepFunctionDerDer = 2*exp(argp1)**2/(1. + exp(argp1))**3/tau**2&
                - exp(argp1)/(1. + exp(argp1))**2/tau**2 &
                                        -2*exp(arg)**2/(1. + exp(arg))**3/tau**2&
                + exp(arg)/(1. + exp(arg))**2/tau**2 
            SmoothedStepFunctionDerDer = (-arg*exp(-arg**2/2.) + argp1*exp(-argp1**2/2.))/sqrt(const_twopi)/tau**2
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

    !Checks
    !write(*,*) 'a = 1e-7'
    !write(*,*) SmoothedStepFunction(1d-7, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDer(1d-7, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDerDer(1d-7, 0.1d0, 0.08d0)
    !write(*,*) 'a = 0.05'
    !write(*,*) SmoothedStepFunction(0.05d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDer(0.05d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDerDer(0.05d0, 0.1d0, 0.08d0)
    !write(*,*) 'a = 0.1'
    !write(*,*) SmoothedStepFunction(0.1d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDer(0.1d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDerDer(0.1d0, 0.1d0, 0.08d0)
    !write(*,*) 'a = 0.15'
    !write(*,*) SmoothedStepFunction(0.15d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDer(0.15d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDerDer(0.15d0, 0.1d0, 0.08d0)
    !write(*,*) 'a = 15'
    !write(*,*) SmoothedStepFunction(15d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDer(15d0, 0.1d0, 0.08d0)
    !write(*,*) SmoothedStepFunctionDerDer(15d0, 0.1d0, 0.08d0)
    !stop

    !write(*,*) 'a = 1e-7'
    !write(*,*) this%w_de(1d-7, 0.1d0, 0.2d0, 0.5d0)
    !write(*,*) 'a = 0.05'
    !write(*,*) this%w_de(5d-2, 0.1d0, 0.2d0, 0.5d0)
    !write(*,*) 'a = 0.1'
    !write(*,*) this%w_de(1d-1, 0.1d0, 0.2d0, 0.5d0)
    !write(*,*) 'a = 0.15'
    !write(*,*) this%w_de(1.5d-1, 0.1d0, 0.2d0, 0.5d0)
    !write(*,*) 'a = 15'
    !write(*,*) this%w_de(15d0, 0.1d0, 0.2d0, 0.5d0)
    !stop

    grhov_t_lcdm = grhov * a * a
    delta = 0

    if(a > 0) then
        grho_t = grhoa2_noDE/a/a + grhov_t_lcdm

        !Get contributions from each of the dark energy bins
        do i = 1, this%de_n_bins

            delta = delta + this%de_bin_amplitudes(i)*SmoothedStepFunction(a, &
                                this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau)
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
                !write(*,'(5e20.8e3)') a, delta, Q, w_bg, w
                !if(abs(a-1.78e-4) < 3e-6) then
                !    write(*,'(17e14.5)') a, delta, Q, w_bg, w
                !endif
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
                        this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau)
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
    real(dl) :: w0
    real(dl) :: w_ratio
    !real(dl) :: da, dlna, am, ap, w0, wm, wp, dw_dlna
    integer :: i

    dQ_dlna = Q_dot/adotoa
    dw_bg_dlna = w_bg_dot/adotoa

    delta = 0
    ddelta_dlna = 0
    d2delta_d2lna = 0

    !Get contributions from each of the dark energy bins
    do i = 1, this%de_n_bins

        delta = delta + this%de_bin_amplitudes(i)*SmoothedStepFunction(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau)

        ddelta_dlna = ddelta_dlna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDer(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau)

        d2delta_d2lna = d2delta_d2lna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDerDer(a, &
            this%de_bin_ai(i), this%de_bin_ai(i+1), this%de_tau)

    enddo

    !write(*,'(7e20.8e3)') a, delta, ddelta_dlna, Q, dQ_dlna, w_bg, dw_bg_dlna

    !This was just a check
    !da = 0.002*a
    !dlna = da/a
    !wm = this%w_de(a - da, delta - ddelta_dlna*dlna, Q - dQ_dlna*dlna, w_bg-dw_bg_dlna*dlna)
    w0 = this%w_de(a, delta, Q, w_bg)
    !!write(*,'(5e20.8e3)') a, delta, Q, w_bg, w0
    !wp = this%w_de(a + da, delta + ddelta_dlna*dlna, Q + dQ_dlna*dlna, w_bg+dw_bg_dlna*dlna)
    !dw_dlna = (wp - wm)/2./dlna
    !if(abs(1+w) .ge. 1e-6) then
    !    dw_dlna = dw_dlna/(1+w)
    !else
    !    dw_dlna = 0
    !endif

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    ! dw/dlog a/(1+w) - derivative of Eq 3 from 1304.3724
    denom = 1 + delta + delta*Q
    t1 = -(1+Q)/3.*d2delta_d2lna/denom
    t2 = delta*Q*dw_bg_dlna/denom
    t3 = ddelta_dlna/3./denom**2*(ddelta_dlna*(1+Q)**2 + 3*Q*(1+w_bg))
    t4 = -dQ_dlna/3./denom**2*(ddelta_dlna - 3*delta*(1+delta)*(1+w_bg))
    !If we are far away from any DE bins, derivative is zero
    if(abs(1+w) .ge. W_CUTOFF) then
        deriv  = (t1 + t2 + t3 + t4)/(1+w)
    else if(delta .le. 1e-8) then
        deriv = 0
    else
        w_ratio = (1+w)/W_CUTOFF
        deriv = (t1 + t2 + t3 + t4)/W_CUTOFF*(2*w_ratio-w_ratio**3)
    endif
    !density perturbation
    !Looks like there is a typo in 1806.10608: in eq 22 they have [delta] =
    ![theta/k^2], in eq 23 they have [delta] = [theta H / k^2]
    ayprime(w_ix) = -3 * adotoa * (this%de_cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ! The deriv term is from w' in [(1+w)v]' = w'v + (1+w)v'
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%de_cs2 - deriv) * y(w_ix + 1) + &
        k * this%de_cs2 * y(w_ix)

    !write(*,'(6e20.8e3)') a, k, w, deriv, ayprime(w_ix), ayprime(w_ix+1)

    !if(abs(k-0.066886) < 0.002) then
    !    write(*,'(27e14.5)') &
    !    a, &
    !    k, &
    !    ayprime(w_ix), &
    !    ayprime(w_ix+1), &
    !    - 3 * adotoa * (this%de_cs2 - w) *  (y(w_ix) + Hv3_over_k), &
    !    -   k * y(w_ix + 1), &
    !    - (1 + w) * k * z, &
    !    - adotoa*deriv* Hv3_over_k,&
    !    - adotoa * (1 - 3 * this%de_cs2 - deriv) * y(w_ix + 1), &
    !      k * this%de_cs2 * y(w_ix),&
    !      adotoa, &
    !      Hv3_over_k, &
    !      w, &
    !      deriv, &
    !      w0
    !endif

    end subroutine TDarkEnergyBins_PerturbationEvolve

    !function TDarkEnergyBins_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
    !    Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    !!Get derivative of anisotropic stress
    !class(TDarkEnergyBins), intent(in) :: this
    !real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &
    !    k, grhov_t, z, k2, yprime(:), y(:), Kf1
    !integer, intent(in) :: w_ix
    !real(dl) :: ppiedot, hdotoh

    !if (this%is_cosmological_constant) then
    !    ppiedot = 0
    !else
    !    hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
    !    ppiedot = 3._dl * dgrhoe + dgqe * &
    !        (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
    !        grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
    !        (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
    !    ppiedot = ppiedot * adotoa / Kf1
    !end if

    !end function TDarkEnergyBins_diff_rhopi_Add_Term

    end module DarkEnergyBins
    !</pavel>
