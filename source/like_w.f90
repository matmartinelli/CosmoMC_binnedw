    ! covariance matrix from paper 1703.05297v1

    module priorwde
    use CosmologyTypes    !da cui pesco CAMBParams
    use Likelihood_Cosmology	!da cui pesco LogLike..
    use MatrixUtils		!per avere l'inversione della matrice di covarianza
    use precision
    implicit none
    private

    logical :: debugging=.false.

    !likelihood variables
    type, extends(TCosmoCalcLikelihood) :: wLikelihood
        integer  :: prior_shape                          !shape of theoretical prior
        integer  :: modelclass                           !assumed model
        real(dl) :: prior_n, prior_xi                    !correlation parameters
        real(dl) :: prior_alpha, prior_beta, prior_gamma !auto-correlation parameters
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
    integer, intent(in) :: nb
    real, intent(in) :: l
    real, intent(in), dimension(nb) :: zbins
    integer :: i, j

    character(LEN=:), allocatable :: data_file
    integer file_unit
    
    if (Ini%Read_Logical('use_priorwde',.false.)) then
       allocate(this)

       !READING PRIOR SETTINGS
       this%prior_shape = Ini%Read_Int('prior_shape')
       this%modelclass  = Ini%Read_Int('theory_model')

       !READING CORRELATION PARAMETERS
       this%prior_n     = Ini%Read_Double('prior_n')
       this%prior_xi    = Ini%Read_Double('prior_xi')

       !READING AUTO CORRELATION PARAMETERS
       this%prior_alpha = Ini%Read_Double('prior_alpha')
       this%prior_beta  = Ini%Read_Double('prior_beta')
       this%prior_gamma = Ini%Read_Double('prior_gamma')

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
    real(dl), dimension(CMB%numbins)             :: diff_vec
    real(dl), dimension(CMB%numbins,CMB%numbins) :: covmat, inv_covmat
    real(dl), dimension(CMB%numbins)             :: autocorr
    real(dl)                                     :: distance, autodist
    integer, parameter                           :: exp_prior=1, CPZ_prior=2
    integer, parameter                           :: quintessence=1, GBD=2, horndeski=3

    mean=0

    !COMPUTING MEAN OF W_I VALUES-----------
    do i=1,CMB%numbins
       mean=mean+CMB%binw(i)
    end do
    mean=mean/CMB%numbins
    if (debugging) write(*,*) 'valore fiduciale', mean


    if (debugging) write(*,*) 'la dimensione dei vettori e', CMB%numbins
    if (debugging) write(*,*) 'il vettore v_w', CMB%binw(i)
    if (debugging) write(*,*) 'il vettore v_wfid', mean

    if (debugging) write(*,*) 'la matrice e', covmat

    !COMPUTING ARRAY OF w_i-mean
    do i=1,CMB%numbins
       diff_vec(i) = CMB%binw(i)-mean
    end do

    !COMPUTING AUTOCORRELATION
    do i=1,CMB%numbins
       if (this%modelclass.eq.quintessence) then
          autodist = -1._dl+1._dl/(1+CMB%binz(i))
       else if ((this%modelclass.eq.GBD).or.(this%modelclass.eq.horndeski)) then
          autodist = log(1._dl/(1+CMB%binz(i)))
       else
          write(*,*) 'MODEL CHOICE 1-3'
          write(*,*) 'YOUR CHOICE DOES NOT EXIST'
          stop
       end if
       autocorr(i) = this%prior_alpha+this%prior_beta*exp(this%prior_gamma*autodist)   
    end do

    !COMPUTING COV MAT AND ITS INVERSE
    if (this%prior_shape.eq.exp_prior) then
       write(*,*) 'NOTHING HERE YET, COME BACK LATER!'
       stop
    else if (this%prior_shape.eq.CPZ_prior) then
       do i=1,CMB%numbins
          do j=1,CMB%numbins
             if (this%modelclass.eq.quintessence) then
                distance = abs((1./(1+CMB%binz(i)))-(1./(1+CMB%binz(j))))
             else if ((this%modelclass.eq.GBD).or.(this%modelclass.eq.horndeski)) then
                distance = abs(log(1./(1+CMB%binz(i)))-log(1./(1+CMB%binz(j))))
             else
                write(*,*) 'MODEL CHOICE 1-3'
                write(*,*) 'YOUR CHOICE DOES NOT EXIST'
                stop
             end if
             covmat(i,j) = sqrt(autocorr(i)*autocorr(j))*1._dl/(1+((distance/this%prior_xi)**this%prior_n))
          end do
       end do
    else
       write(*,*) 'BAD CHOICE OF PRIOR'
       write(*,*) 'CHOOSE AN EXISTING ONE!!!'
       stop
    end if

    inv_covmat(:,:) = covmat(:,:)
    call Matrix_Inverse(inv_covmat)

    !COMPUTING CHI2
    chi2 = 0._dl

    chi2 = dot_product( diff_vec, MatMul(inv_covmat,diff_vec))

    if (debugging) then
       open(78, file='chi2_priorwde.dat', status='unknown', position='append')
       write(78,*) CMB%binw, diff_vec, chi2
       close(78)
    end if

    w_Lnlike = chi2/2._dl 

    if (feedback.gt.0) write(*,*) 'Prior Like =', w_Lnlike

    end function  w_LnLike


    end module priorwde
