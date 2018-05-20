    ! covariance matrix from paper 1703.05297v1

    module priorwde
    use CosmologyTypes    !da cui pesco CAMBParams
    use Likelihood_Cosmology	!da cui pesco LogLike..
    use MatrixUtils		!per avere l'inversione della matrice di covarianza
    use precision
    implicit none
    private

    logical :: debugging=.true.

    !likelihood variables
    type, extends(TCosmoCalcLikelihood) :: wLikelihood
	real(dl), dimension(:),allocatable :: v_w, v_wf, v
	real(dl), dimension(:,:),allocatable :: cm, inv_cm
	integer :: n
    contains

    procedure :: LogLikeDataParams => w_LnLike
    end type wLikelihood

    public wLikelihood, wLikelihood_Add
    contains

!...................................................................

    subroutine wLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(wLikelihood), pointer :: this

    integer :: i, j

    character(LEN=:), allocatable :: data_file
    integer file_unit
	
	allocate(this)

	this%n=Ini%Read_Double('numbins')
!        this%n=nb
!        if (debugging) write(*,*) 'sono 1'

        if (allocated(this%v_wf) .eqv. .false.) allocate (this%v_w(this%n), this%v_wf(this%n), this%cm(this%n,this%n) , this%inv_cm(this%n,this%n))

!        if (debugging) write(*,*) 'sono 2'

        data_file=Ini%Read_String_Default('data_file',trim(DataDir)//'cov_mat_wde.txt')

        open(newunit=file_unit, file=trim(data_file), status='old')             
        read(file_unit,*) ((this%cm(i,j), j=1,this%n), i=1,this%n)
        close(file_unit)

	if (debugging) write(*,*) this%cm

        this%needs_background_functions = .true.
        call LikeList%Add(this)   !added to the list of likelihoods
    

    end subroutine wLikelihood_Add

!-------------------------------------------------------------------

    real(mcp) function w_LnLike(this, CMB, DataParams)
    Class(wLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) DataParams(:)                              !capire cos'è

    integer :: i,j
    real :: chi2

!	this%n=CMB%numbins


	do i=1,this%n
		this%v_w(i)=CMB%binw(i)
		this%v_wf(i)=-1._dl
	end do

	if (debugging) write(*,*) 'la dimensione dei vettori e', this%n
	if (debugging) write(*,*) 'il vettore v_w', this%v_w
        if (debugging) write(*,*) 'il vettore v_wfid', this%v_wf

!	if (debugging) write(*,*) 'la matrice e', this%cm

        if (allocated(this%v) .eqv. .false.) allocate (this%v(CMB%numbins))
	do i=1,CMB%numbins
		this%v(i) = this%v_w(i)-this%v_wf(i)
        end do

        chi2 = 0._dl

	!definisco la matrice inversa tramite la subroutine Matrix_Inverse(M) :: RICORDA utilizzabile se matrice simmetrica definita positiva
        this%inv_cm = this%cm
	call Matrix_Inverse(this%inv_cm)      

	do i=1,CMB%numbins
		do j=1,CMB%numbins
			chi2 = chi2 + this%v(i)*this%inv_cm(i,j)*this%v(j)
		end do
	end do

	w_Lnlike = chi2/2._dl 

    end function  w_LnLike


    end module priorwde
