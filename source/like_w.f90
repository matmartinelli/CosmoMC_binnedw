    ! covariance matrix from paper 1703.05297v1

    module priorwde
    use CosmologyTypes    !da cui pesco CAMBParams
    use Likelihood_Cosmology	!da cui pesco LogLike..
    use MatrixUtils		!per avere l'inversione della matrice di covarianza
    use precision
    implicit none
    private

    logical :: debugging=.false.
    logical :: debugging_paramsprior=.false.

    !likelihood variables
    type, extends(TCosmoCalcLikelihood) :: wLikelihood
        integer  :: prior_shape                          !shape of theoretical prior
        integer  :: modelclass                           !assumed model
        real :: prior_n, prior_xi                    !correlation parameters
        real :: prior_alpha, prior_beta, prior_gamma !auto-correlation parameters
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
    integer, parameter                           :: exp_prior=1, CPZ_prior=2
    integer, parameter                           :: quintessence=1, GBD=2, horndeski=3

    
    if (Ini%Read_Logical('use_priorwde',.false.)) then
       allocate(this)

       !READING PRIOR SETTINGS
       this%prior_shape = Ini%Read_Int('prior_shape')
       this%modelclass  = Ini%Read_Int('theory_model')


       !CORRELATION PARAMETERS
       if ((this%prior_shape.eq.CPZ_prior).or.(this%prior_shape.eq.exp_prior)) then
          if (this%prior_shape.eq.exp_prior) then
              if (this%modelclass.eq.quintessence) then
	   	   this%prior_n     = 1.8
	           this%prior_xi    = 0.7
              else if (this%modelclass.eq.GBD) then
		   this%prior_n     = 1.3
	           this%prior_xi    = 0.3
	      else if (this%modelclass.eq.horndeski) then
		   this%prior_n     = 1.2
	           this%prior_xi    = 0.3
	      else 
                   write(*,*) 'MODEL CHOICE 1-3'
	           write(*,*) 'YOUR CHOICE DOES NOT EXIST'
	           stop
              end if
          else
	      this%prior_n=2
              if (this%modelclass.eq.quintessence) then
	 	   this%prior_xi    = 0.6
              else if ((this%modelclass.eq.GBD).or.(this%modelclass.eq.horndeski)) then
                   this%prior_xi    = 0.2
              else 
                   write(*,*) 'MODEL CHOICE 1-3'
	           write(*,*) 'YOUR CHOICE DOES NOT EXIST'
	           stop
              end if
          end if
       else
          write(*,*) 'BAD CHOICE OF PRIORSHAPE'
          write(*,*) 'CHOOSE AN EXISTING ONE!!!'
          stop
       end if


       !AUTOCORRELATION PARAMETERS
       if (this%modelclass.eq.quintessence) then
      	   this%prior_alpha = 0.03
           this%prior_beta  = 0.3
           this%prior_gamma = 6.5
       else if (this%modelclass.eq.GBD) then
	   this%prior_alpha = 0.05
           this%prior_beta  = 0.8
           this%prior_gamma = 1.8
       else if (this%modelclass.eq.horndeski) then
	   this%prior_alpha = 0.05
           this%prior_beta  = 0.8
           this%prior_gamma = 2
       else 
           write(*,*) 'MODEL CHOICE 1-3'
           write(*,*) 'YOUR CHOICE DOES NOT EXIST'
           stop
       end if

       !PRINTING MODELS AND PARAMETERS
       if (debugging_paramsprior) then 
	   write(*,*) 'prior', this%prior_shape
	   write(*,*) 'model', this%modelclass
	   write(*,*) 'prior_n', this%prior_n
	   write(*,*) 'prior_xi', this%prior_xi
	   write(*,*) 'prior_alpha', this%prior_alpha
	   write(*,*) 'prior_beta', this%prior_beta
	   write(*,*) 'prior_gamma', this%prior_gamma
       end if



       this%needs_background_functions = .true.
       call LikeList%Add(this)   !added to the list of likelihoods
    end if

    end subroutine wLikelihood_Add

!-------------------------------------------------------------------

    real(mcp) function w_LnLike(this, CMB, DataParams)
    Class(wLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) DataParams(:)                             

    integer :: i,j,d
    real :: chi2
    real(dl), dimension(:), allocatable          :: diff_vec, wfid, wi, gpreds
    real(dl), dimension(:,:), allocatable        :: covmat, inv_covmat
    real(dl), dimension(:), allocatable        :: autocorr
    real(dl)                                     :: distance, autodist
    integer, parameter                           :: exp_prior=1, CPZ_prior=2
    integer, parameter                           :: quintessence=1, GBD=2, horndeski=3


    if(CMB%mode.ne.3) then
        d=CMB%numbins
        allocate (diff_vec(d), wfid(d), wi(d), gpreds(d))
        allocate (covmat(d,d),inv_covmat(d,d),autocorr(d))
        do i=1,CMB%numbins
           wi(i)=CMB%binw(i)
        end do
        gpreds(1) = (CMB%binz(1))/2.
        do i=2, CMB%numbins
           gpreds(i) = (CMB%binz(i)+CMB%binz(i-1))/2.
        end do
    else
        d=CMB%numbins+1
        allocate (diff_vec(d), wfid(d), wi(d), gpreds(d))
        allocate (covmat(d,d),inv_covmat(d,d),autocorr(d))
        wi(1)=CMB%binw0
        do i=1,CMB%numbins
           wi(i+1)=CMB%binw(i)
        end do
        gpreds(1) = 0
        do i=1, CMB%numbins
           gpreds(i+1) = (CMB%binz(i)+CMB%binz(i-1))/2.
        end do
    end if


    if (debugging) write(*,*) 'la dimensione dei vettori e', d
    if (debugging) write(*,*) 'il vettore v_w', wi
    
    
    !COMPUTING MEAN OF W_I VALUES-----------
    wfid(1) = (wi(1)+wi(2))/2
    do i=2,d
       wfid(i)=(wi(i-1)+wi(i)+wi(i+1))/3
    end do
    wfid(d)= (wi(d-1)+wi(d))/2
    if (debugging) write(*,*) 'il vettore v_wfid', wfid


    !COMPUTING ARRAY OF w_i-w_fid
    do i=1,d
       diff_vec(i) = wi(i)-wfid(i)
    end do

    if (debugging) write(*,*) 'il vettore differenza e', diff_vec

    !COMPUTING AUTOCORRELATION
    do i=1,d
       if (this%modelclass.eq.quintessence) then
          autodist = -1._dl+(1._dl/(1+gpreds(i)))
          if (debugging) write(*,*) 'quintessence autodistance', autodist
       else if ((this%modelclass.eq.GBD).or.(this%modelclass.eq.horndeski)) then
          autodist = log(1._dl/(1+gpreds(i)))
          if (debugging) write(*,*) 'horndeski autodistance', autodist
       else
          write(*,*) 'MODEL CHOICE 1-3'
          write(*,*) 'YOUR CHOICE DOES NOT EXIST'
          stop
       end if
       autocorr(i) = this%prior_alpha+(this%prior_beta*(exp(this%prior_gamma*autodist)))   
       if (debugging) write(*,*) 'autocorr',i, autocorr(i)
    end do


    !COMPUTING COV MAT AND ITS INVERSE
    if ((this%prior_shape.eq.CPZ_prior).or.(this%prior_shape.eq.exp_prior)) then
       do i=1,d
          do j=1,d
             if (this%modelclass.eq.quintessence) then
                distance = abs((1./(1+gpreds(i)))-(1./(1+gpreds(j))))
                if (debugging) write(*,*) 'quintessence distance', distance
             else if ((this%modelclass.eq.GBD).or.(this%modelclass.eq.horndeski)) then
                distance = abs(log(1./(1+gpreds(i)))-log(1./(1+gpreds(j))))
                if (debugging) write(*,*) 'horndeski distance', distance
             else
                write(*,*) 'MODEL CHOICE 1-3'
                write(*,*) 'YOUR CHOICE DOES NOT EXIST'
                stop
             end if
	     if (this%prior_shape.eq.CPZ_prior) then
            	covmat(i,j) = sqrt(autocorr(i)*autocorr(j))/(1+((distance/this%prior_xi)**this%prior_n))
                if (debugging) write(*,*) 'autocorr',i, autocorr(i)
                if (debugging) write(*,*) 'autocorr',j, autocorr(j)
                if (debugging) write(*,*) 'CPZ Cor',i,j, 1._dl/(1+((distance/(this%prior_xi))**this%prior_n))
                if (debugging) write(*,*) 'CPZ C',i,j, covmat(i,j)
	     else
		covmat(i,j) = sqrt(autocorr(i)*autocorr(j))*exp(-((distance/this%prior_xi)**this%prior_n))
                if (debugging) write(*,*) 'exp Cor',i,j, exp(-((distance/this%prior_xi)**this%prior_n))
                if (debugging) write(*,*) 'exp C',i,j, covmat(i,j)
	     end if
          end do
       end do
    else
       write(*,*) 'BAD CHOICE OF PRIOR'
       write(*,*) 'CHOOSE AN EXISTING ONE!!!'
       stop
    end if

    if (debugging) write(*,*) 'la matrice e'
    if (debugging) write(*,*) covmat

    inv_covmat(:,:) = covmat(:,:)
    call Matrix_Inverse(inv_covmat)

    if (debugging) write(*,*) 'la matrice inversa dopo e'
    if (debugging) write(*,*) inv_covmat


    
    !COMPUTING CHI2
    chi2 = 0._dl

!    chi2 = dot_product( diff_vec, MatMul(inv_covmat,diff_vec))

	do i=1,d
		do j=1,d
                        if (debugging) write(*,*) 'elemento', i, j, 'del chi2', diff_vec(i)*inv_covmat(i,j)*diff_vec(j)
                        if (j.eq.i) then
			    chi2 = chi2 + diff_vec(i)*inv_covmat(i,j)*diff_vec(j)/2._dl
                        else
                            chi2 = chi2 + diff_vec(i)*inv_covmat(i,j)*diff_vec(j)
                        end if
		end do
	end do

!    if (debugging) then
!       open(78, file='chi2_priorwde.dat', status='unknown', position='append')
!       write(78,*) CMB%binw, diff_vec, chi2
!       close(78)
!    end if

    w_Lnlike = chi2/2._dl 

    if (feedback.gt.0) write(*,*) 'Prior Like =', w_Lnlike

    end function  w_LnLike


    end module priorwde
