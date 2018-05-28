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
	real(dl), dimension(:),allocatable :: v_w, v_wf, v, zb
	real(dl), dimension(:,:),allocatable :: cm, inv_cm
	integer :: n
        real :: xi, expo
    contains

    procedure :: LogLikeDataParams => w_LnLike
    end type wLikelihood

    public wLikelihood, wLikelihood_Add
    contains

!...................................................................

    subroutine wLikelihood_Add(LikeList, nb, zbins, l, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(wLikelihood), pointer :: this
    integer, intent(in) :: nb
    real, intent(in) :: l
    real, intent(in), dimension(nb) :: zbins
    integer :: i, j

    character(LEN=:), allocatable :: data_file
    integer file_unit
    
	if (Ini%Read_Logical('use_priorwde',.false.)) then	
		allocate(this)
  	
		this%n=nb
		if (debugging) write(*,*) 'numero di bin', this%n
                this%zb=zbins
		if (debugging) write(*,*) 'array redshift binnati', this%zb
		this%xi=l
		if (debugging) write(*,*) 'valore di xi', this%xi     !att. assunto essere uguale a correlation lenght in ingresso
		this%expo=3
		if (debugging) write(*,*) 'valore esponenziale della matrice', this%expo      !capire che esponente dare
                
        	if (allocated(this%v_wf) .eqv. .false.) allocate (this%v_w(this%n), this%v_wf(this%n), this%cm(this%n,this%n) , this%inv_cm(this%n,this%n))

		do i=1,this%n
			do j=1,this%n
				this%cm(i,j) = 1._dl/(1+((sign(this%zb(i)-this%zb(j),0.9_dl)/this%xi)**this%expo))
			end do
		end do

		if (debugging) write(*,*) this%cm

        	this%needs_background_functions = .true.
        	call LikeList%Add(this)   !added to the list of likelihoods
    	end if

    end subroutine wLikelihood_Add

!-------------------------------------------------------------------

    real(mcp) function w_LnLike(this, CMB, DataParams)
    Class(wLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) DataParams(:)                              !capire cos'è

    integer :: i,j
    real :: chi2, mean

        mean=0
	do i=1,this%n
		mean=mean+CMB%binw(i)
	end do
	mean=mean/this%n
	if (debugging) write(*,*) 'valore fiduciale', mean
	
	do i=1,this%n
		this%v_w(i)=CMB%binw(i)
		this%v_wf(i)=mean
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
!	call Matrix_InverseAsymm(this%inv_cm)

	do i=1,CMB%numbins
		do j=1,CMB%numbins
			chi2 = chi2 + this%v(i)*this%inv_cm(i,j)*this%v(j)
		end do
	end do

	if (debugging) then
		open(78, file='chi2_priorwde.dat', status='unknown', position='append')
		write(78,*) this%v_w, this%v, chi2
        	close(78)
	end if

	w_Lnlike = chi2/2._dl 

    end function  w_LnLike


    end module priorwde
