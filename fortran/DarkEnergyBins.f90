    !<pavel>
    module DarkEnergyBinsModule
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    DE_CUTOFF = 20

    type, extends(TDarkEnergyModel) :: TDarkEnergyBins
        real(dl) :: de_n_bins
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
    procedure :: PerturbationEvolve => TDarkEnergyBins_PerturbationEvolve
    procedure :: BackgroundDensityAndPressure => TDarkEnergyBins_BackgroundDensityAndPressure
    procedure, private :: SmoothedStepFunction
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

    call this%TDarkEnergyModel%ReadParams(Ini)
    !Read the number of DE bins and allocate the arrays
    this%de_n_bins  = Ini%Read_Double('DE_n_bins')
    this%de_tau  = Ini%Read_Double('DE_tau')
    this%de_cs2 = Ini%Read_Double('DE_cs2')
    allocate(this%de_bin_ai(this%de_n_bins+1))
    allocate(this%de_bin_amplitudes(this%de_n_bins))
    !Read the DE bin boundaries and amplitudes in the bin
    this%de_bin_ai  = Ini%Read_Double_Array('DE_bin_ai')
    this%de_bin_amplitudes  = Ini%Read_Double('DE_bin_amplitudes')

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

    !1/[1 + exp{ln(a/ai)/tau}]
    function SmoothedStepFunction(this, a, ai)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: x
    real(dl) :: arg

        arg = log(a/ai)/this%de_tau

        !Make sure the argument of exponential sensible
        if(arg < -DE_CUTOFF) then
            return 1
        else if(arg < DE_CUTOFF) then
            return 1./(1. + exp(x))
        else
            return 0
        endif
    end function

    !Derivative of 1/[1 + exp{ln(a/ai)/tau}] wrt ln a
    function SmoothedStepFunctionDer(this, a, ai)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: x
    real(dl) :: arg

        arg = log(a/ai)/this%de_tau

        !Make sure the argument of exponential sensible
        if( (arg < -DE_CUTOFF) .or. (arg > DE_CUTOFF) ) then
            return 0
        else
            return -exp(x)/(1. + exp(x))**2/this%de_tau
        endif
    end function

    subroutine TDarkEnergyBins_BackgroundDensityAndPressure(this, grhov, a, grhov_t, grhoa2_noDE, w, w_bg)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(in) :: grhoa2_noDE
    real(dl), intent(in) :: w_bg
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl), optional, intent(in) :: w_bg
    real(dl) :: grhov_t_lcdm, grho_t, grhov_t_beyond
    real(dl) :: arg1, arg2
    real(dl) :: delta, Q
    integer :: i

    grhov_t_lcdm = grhov * a * a
    grho_t = grhoa2_noDE/a2 + grhov_t_lcdm
    delta = 0

    !Get contributions from each of the dark energy bins
    do i = 1, de_n_bins

        delta = delta + this%de_bin_amplitudes(i)*SmoothedStepFunction(a, this%de_bin_ai(i+1))
        delta = delta - this%de_bin_amplitudes(i)*SmoothedStepFunction(a, this%de_bin_ai(i))

    enddo

    grhov_t_beyond = delta*grho_t
    grhov_t = grhov_t_lcdm + grhov_t_beyond

    !Factor that goes into the DE equation of state
    Q = grhoa2_noDE/a2/grhov_t_lcdm
    if (present(w)) then
        if (not present(w_bg)) then
            stop 'w_bg'
        else
            w = this%w_de(a, delta, Q, w_bg)
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

    ddelta_dlna = 0

    !Get contributions from each of the dark energy bins
    do i = 1, de_n_bins

        ddelta_dlna = ddelta_dlna + this%de_bin_amplitudes(i)*SmoothedStepFunctionDer(a, this%de_bin_ai(i+1))
        ddelta_dlna = ddelta_dlna - this%de_bin_amplitudes(i)*SmoothedStepFunctionDer(a, this%de_bin_ai(i))

    enddo

    TDarkEnergyBins_w_de = -1 + Q*delta/(1+delta+delta*Q)*(1+w_bg)&
        -(1+Q)/3./(1+delta+delta*Q)*ddelta_dlna

    end function TDarkEnergyBins_w_de

    !I am here

    subroutine TDarkEnergyBins_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TDarkEnergyBins), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    ! dw/dlog a/(1+w)
    deriv  = 
    !density perturbation
    !Looks like there is a typo in 1806.10608: in eq 22 they have [delta] =
    ![theta/k^2], in eq 23 they have [delta] = [theta H / k^2]
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ! The deriv term is from w' in [(1+w)v]' = w'v + (1+w)v'
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + &
        k * cs2 * y(w_ix)

    end subroutine TDarkEnergyBins_PerturbationEvolve


    end module DarkEnergyBinsModule
    !</pavel>
