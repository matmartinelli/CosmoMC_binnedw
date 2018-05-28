    module DataLikelihoodList
    use likelihood
    use settings
    use CosmologyTypes
    implicit none

    contains

    subroutine SetDataLikelihoods(Ini)
    use HST
    use snovae
    use CMBLikelihoods
    use bao
    use mpk
    use wigglez
    use szcounts !Anna
    use wl
    use ElementAbundances
!FGmod-------
    use priorwde
!------------    
    
    class(TSettingIni), intent(in) :: Ini
!FGmod-------
    integer :: init_nbin, i
    real :: corr_l
    real, dimension(:), allocatable :: zbin_read
    character(LEN=100) :: binnum
    character(LEN=100) :: num, par

    call Ini%Read('numbins',init_nbin)

    allocate(zbin_read(init_nbin))
    do i=1,init_nbin
       write(binnum, *) i
       zbin_read(i) = Ini%Read_Double('param[binz'//trim(adjustl(binnum))//']')
    end do
    corr_l= Ini%Read_Double('param[corr_l]')
!------------
    CosmoSettings%get_sigma8 = Ini%Read_Logical('get_sigma8',.false.)

    call CMBLikelihood_Add(DataLikelihoods, Ini)

    call AbundanceLikelihood_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    call MPKLikelihood_Add(DataLikelihoods, Ini)

    if (use_mpk) call WiggleZLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

    call SZLikelihood_Add(DataLikelihoods, Ini) !Anna

    call WLLikelihood_Add(DataLikelihoods, Ini)


!FGmod-------
    call wLikelihood_Add(DataLikelihoods, init_nbin, zbin_read, corr_l, Ini)
!------------

    end subroutine SetDataLikelihoods


    end module DataLikelihoodList
