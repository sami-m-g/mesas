! -*- f90 -*-
subroutine solve(J_ts, Q_ts, SAS_lookup, P_list, weights_ts, sT_init_ts, dt, &
        verbose, debug, warning, jacobian,&
        mT_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
        n_substeps, numcomponent_list, numbreakpt_list, numflux, numsol, max_age, &
        timeseries_length, numcomponent_total, numbreakpt_total, &
        sT_ts, pQ_ts, WaterBalance_ts, &
        mT_ts, mQ_ts, mR_ts, C_Q_ts, ds_ts, dm_ts, dC_ts, SoluteBalance_ts)
    implicit none

    ! Start by declaring and initializing all the variables we will be using
    integer, intent(in) :: n_substeps, numflux, numsol, max_age, &
            timeseries_length, numcomponent_total, numbreakpt_total
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning, jacobian
    real(8), intent(in), dimension(0:timeseries_length - 1) :: J_ts
    real(8), intent(in), dimension(0:numflux - 1, 0:timeseries_length - 1) :: Q_ts
    real(8), intent(in), dimension(0:numcomponent_total - 1, 0:timeseries_length - 1) :: weights_ts
    real(8), intent(in), dimension(0:numbreakpt_total - 1, 0:timeseries_length - 1) :: SAS_lookup
    real(8), intent(in), dimension(0:numbreakpt_total - 1, 0:timeseries_length - 1) :: P_list
    real(8), intent(in), dimension(0:numsol - 1, 0:timeseries_length - 1) :: C_J_ts
    real(8), intent(in), dimension(0:numflux - 1, 0:numsol - 1, 0:timeseries_length - 1) :: alpha_ts
    real(8), intent(in), dimension(0:numsol - 1, 0:timeseries_length - 1) :: k1_ts
    real(8), intent(in), dimension(0:numsol - 1, 0:timeseries_length - 1) :: C_eq_ts
    real(8), intent(in), dimension(0:numsol - 1) :: C_old
    real(8), intent(in), dimension(0:max_age - 1) :: sT_init_ts
    real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mT_init_ts
    integer, intent(in), dimension(0:numflux - 1) :: numcomponent_list
    integer, intent(in), dimension(0:numcomponent_total - 1) :: numbreakpt_list
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: C_Q_ts
    real(8), intent(out), dimension(0:numbreakpt_total-1, 0:numflux - 1, 0:numsol - 1, 0:timeseries_length - 1) :: dC_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length) :: sT_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numsol - 1) :: mT_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numbreakpt_total-1) :: ds_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1) :: pQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: mR_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1) :: WaterBalance_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: SoluteBalance_ts
    real(8), dimension(0:numbreakpt_total-1, 0:numflux - 1, 0:timeseries_length-1) :: dW_ts
    real(8), dimension(0:numflux - 1, 0:timeseries_length - 1) :: P_old
    integer, dimension(0:numcomponent_total) :: breakpt_index_list
    integer, dimension(0:numflux) :: component_index_list
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_top_start
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_bot_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_top
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: PQcum_bot
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: PQcum_top
    integer, dimension(0:numcomponent_total-1, 0:timeseries_length * n_substeps - 1) :: leftbreakpt_bot
    integer, dimension(0:numcomponent_total-1, 0:timeseries_length * n_substeps - 1) :: leftbreakpt_top
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: pQ_temp
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: pQ_aver
    real(8), dimension(0:numflux - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mQ_temp
    real(8), dimension(0:numflux - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mQ_aver
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mR_temp
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mR_aver
    real(8), dimension(0:numbreakpt_total - 1, 0:timeseries_length * n_substeps - 1) :: fs_temp
    real(8), dimension(0:numbreakpt_total - 1, 0:timeseries_length * n_substeps - 1) :: fs_aver
    real(8), dimension(0:numbreakpt_total - 1, 0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: fsQ_temp
    real(8), dimension(0:numbreakpt_total - 1, 0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: fsQ_aver
    real(8), dimension(0:numbreakpt_total - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fm_temp
    real(8), dimension(0:numbreakpt_total - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fm_aver
    real(8), dimension(0:numbreakpt_total - 1, 0:numflux-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fmQ_temp
    real(8), dimension(0:numbreakpt_total - 1, 0:numflux-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fmQ_aver
    real(8), dimension(0:numbreakpt_total - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fmR_temp
    real(8), dimension(0:numbreakpt_total - 1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: fmR_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mT_start
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mT_temp
    real(8), dimension(0:numbreakpt_total-1, 0:timeseries_length * n_substeps - 1) :: ds_start
    real(8), dimension(0:numbreakpt_total-1, 0:timeseries_length * n_substeps - 1) :: ds_temp
    real(8), dimension(0:numbreakpt_total-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: dm_start
    real(8), dimension(0:numbreakpt_total-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: dm_temp
    real(8) :: one8, norm
    real(8) :: dS, dP, dSe, dPe, dSs, dPs
    real(8) :: h, hr
    integer, dimension(4) :: rk_coeff
    real(8), dimension(5) :: rk_time
    character(len = 128) :: tempdebugstring
    integer :: iT_substep, iT, jt, iT_s, jt_s, iT_prev, jt_substep
    integer :: iq, s, M, N, ip, ic, c, rk
    integer :: carry
    integer :: leftbreakpt
    real(8) :: PQcum_component
    real(8) :: STcum_in
    integer :: jt_this, topbot
    integer :: na
    integer :: ia
    integer :: i
    real(8) :: dif, grad
    logical :: foundit
    real :: start, finish, finish2

    !call f_verbose('...Initializing arrays...')
    one8 = 1.0
    rk_time = (/0.0D0, 0.5D0, 0.5D0, 1.0D0, 1.0D0/)
    rk_coeff = (/1, 2, 2, 1/)
    norm = 1.0 / n_substeps / n_substeps


    C_Q_ts = 0.
    sT_ts = 0.
    mT_ts = 0.
    ds_ts = 0.
    dm_ts = 0.
    dC_ts = 0.
    dW_ts = 0.
    pQ_ts = 0.
    mQ_ts = 0.
    mR_ts = 0.
    WaterBalance_ts = 0.
    SoluteBalance_ts = 0.
    P_old = 0.
    breakpt_index_list = 0
    component_index_list = 0
    STcum_top_start = 0.
    STcum_bot_start = 0.
    STcum_bot = 0.
    STcum_top = 0.
    PQcum_bot = 0.
    PQcum_top = 0.
    leftbreakpt_bot = 0
    leftbreakpt_top = 0
    pQ_temp = 0.
    pQ_aver = 0.
    mQ_temp = 0.
    mQ_aver = 0.
    mR_temp = 0.
    mR_aver = 0.
    fs_temp = 0.
    fs_aver = 0.
    fsQ_temp = 0.
    fsQ_aver = 0.
    fm_temp = 0.
    fm_aver = 0.
    fmQ_temp = 0.
    fmQ_aver = 0.
    fmR_temp = 0.
    fmR_aver = 0.
    sT_start = 0.
    sT_temp = 0.
    mT_start = 0.
    mT_temp = 0.
    ds_start = 0.
    ds_temp = 0.
    dm_start = 0.
    dm_temp = 0.
    iT_prev = -1

    ! The list of probabilities in each sas function is a 1-D array.
    ! breakpt_index_list gives the starting index of the probabilities (P) associated
    ! with each flux
    breakpt_index_list(0) = 0
    component_index_list(0) = 0
    do iq = 0, numflux - 1
        component_index_list(iq + 1) = component_index_list(iq) + numcomponent_list(iq)
        do ic = component_index_list(iq), component_index_list(iq+1) - 1
            breakpt_index_list(ic + 1) = breakpt_index_list(ic) + numbreakpt_list(ic)
        enddo
    enddo
    !call f_debug('breakpt_index_list', one8 * breakpt_index_list(:))

    ! modify the number of ages and the timestep by a facotr of n_substeps
    M = max_age * n_substeps
    N = timeseries_length * n_substeps
    h = dt / n_substeps

    !call f_verbose('...Setting initial conditions...')
    sT_ts(:, 0) = sT_init_ts
    do s = 0, numsol - 1
        mT_ts(:, 0, s) = mT_init_ts(:, s)
    end do


    !$acc data &
    !$acc copyin(component_index_list(:numflux)) &
    !$acc copyin(rk_coeff(:5)) &
    !$acc copyin(SAS_lookup(:,:)) &
    !$acc copyin(J_ts(:)) &
    !$acc copyin(breakpt_index_list(:)) &
    !$acc copyin(k1_ts(:numsol-1,:)) &
    !$acc copyin(C_J_ts(:numsol-1,:)) &
    !$acc copyin(weights_ts(:,:)) &
    !$acc copyin(alpha_ts(:numflux-1,:numsol-1,:)) &
    !$acc copyin(rk_time(:)) &
    !$acc copyin(Q_ts(:numflux-1,:)) &
    !$acc copyin(P_list(:,:)) &
    !$acc copyin(C_eq_ts(:numsol-1,:)) &
    !$acc copyin(numbreakpt_list(:)) &
    !$acc copyin(jacobian) &
    !$acc copyin(sT_init_ts(:),mT_init_ts(:,:numsol-1)) &
    !$acc copyin(dm_start(:numbreakpt_total-1,:numsol-1,:n-1)) &
    !$acc copyin(ds_start(:numbreakpt_total-1,:n-1)) &
    !$acc copyin(sT_start(:n-1)) &
    !$acc copyin(mT_start(:numsol-1,:n-1)) &
    !$acc create(fmR_aver(:numbreakpt_total-1,:numsol-1,:n-1)) &
    !$acc create(fm_aver(:,:numsol-1,:n-1)) &
    !$acc create(fmQ_aver(:numbreakpt_total-1,:numflux-1,:numsol-1,:n-1)) &
    !$acc create(fsQ_aver(:numbreakpt_total-1,:numflux-1,:n-1)) &
    !$acc create(mR_aver(:numsol-1,:n-1)) &
    !$acc create(pQ_aver(:numflux-1,:n-1)) &
    !$acc create(mQ_aver(:numflux-1,:numsol-1,:n-1)) &
    !$acc create(fs_aver(:,:n-1)) &
    !$acc create(fmQ_temp(:numbreakpt_total-1,:numflux-1,:numsol-1,:n-1)) &
    !$acc create(mQ_temp(:numflux-1,:numsol-1,:n-1)) &
    !$acc create(fmR_temp(:numbreakpt_total-1,:numsol-1,:n-1)) &
    !$acc create(fs_temp(:,:n-1)) &
    !$acc create(fm_temp(:,:numsol-1,:n-1)) &
    !$acc create(fmQ_temp(:numbreakpt_total-1,:,:,:n-1)) &
    !$acc create(fsQ_temp(:numbreakpt_total-1,:numflux-1,:n-1)) &
    !$acc create(pQ_temp(:numflux-1,:n-1)) &
    !$acc create(mR_temp(:numsol-1,:n-1)) &
    !$acc create(dm_temp(:numbreakpt_total-1,:numsol-1,:n-1)) &
    !$acc create(mT_temp(:numsol-1,:n-1)) &
    !$acc create(sT_temp(:n-1)) &
    !$acc create(ds_temp(:numbreakpt_total-1,:n-1)) &
    !$acc create(leftbreakpt_top(:,:n-1)) &
    !$acc create(leftbreakpt_bot(:,:n-1)) &
    !$acc create(STcum_bot(:n-1)) &
    !$acc create(PQcum_bot(:numflux-1,:n-1)) &
    !$acc create(STcum_top(:n-1)) &
    !$acc create(PQcum_top(:numflux-1,:n-1)) &
    !$acc copyin(STcum_bot_start(:)) &
    !$acc copyin(STcum_top_start(:)) &
    !$acc copyin(dW_ts(:,:numflux-1,:)) &
    !$acc copyin(dC_ts(:,:numflux-1,:numsol-1,:)) &
    !$acc copyin(sT_ts) &
    !$acc copyin(mT_ts) &
    !$acc copyin(ds_ts) &
    !$acc copyin(dm_ts) &
    !$acc copyin(dC_ts) &
    !$acc copyin(dW_ts) &
    !$acc copyin(pQ_ts) &
    !$acc copyin(mQ_ts) &
    !$acc copyin(mR_ts)

    !call f_verbose('...Starting main loop...')
    do iT = 0, max_age - 1

        ! Start the substep loop
        do iT_substep = 0, n_substeps - 1

            iT_s = iT * n_substeps +  iT_substep
            call cpu_time(start)
            !$acc kernels loop &
            !$acc present(component_index_list(:numflux)) &
            !$acc present(rk_coeff(:5)) &
            !$acc present(SAS_lookup(:,:)) &
            !$acc present(J_ts(:)) &
            !$acc present(breakpt_index_list(:)) &
            !$acc present(k1_ts(:numsol-1,:)) &
            !$acc present(C_J_ts(:numsol-1,:)) &
            !$acc present(weights_ts(:,:)) &
            !$acc present(alpha_ts(:numflux-1,:numsol-1,:)) &
            !$acc present(rk_time(:)) &
            !$acc present(Q_ts(:numflux-1,:)) &
            !$acc present(P_list(:,:)) &
            !$acc present(C_eq_ts(:numsol-1,:)) &
            !$acc present(numbreakpt_list(:)) &
            !$acc present(jacobian) &
            !$acc present(dm_start(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(ds_start(:numbreakpt_total-1,:n-1)) &
            !$acc present(sT_start(:n-1)) &
            !$acc present(mT_start(:numsol-1,:n-1)) &
            !$acc present(fmR_aver(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(fm_aver(:,:numsol-1,:n-1)) &
            !$acc present(fmQ_aver(:numbreakpt_total-1,:numflux-1,:numsol-1,:n-1)) &
            !$acc present(fsQ_aver(:numbreakpt_total-1,:numflux-1,:n-1)) &
            !$acc present(mR_aver(:numsol-1,:n-1)) &
            !$acc present(pQ_aver(:numflux-1,:n-1)) &
            !$acc present(mQ_aver(:numflux-1,:numsol-1,:n-1)) &
            !$acc present(fs_aver(:,:n-1)) &
            !$acc present(fmQ_temp(:numbreakpt_total-1,:numflux-1,:numsol-1,:n-1)) &
            !$acc present(mQ_temp(:numflux-1,:numsol-1,:n-1)) &
            !$acc present(fmR_temp(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(fs_temp(:,:n-1)) &
            !$acc present(fm_temp(:,:numsol-1,:n-1)) &
            !$acc present(fmQ_temp(:numbreakpt_total-1,:,:,:n-1)) &
            !$acc present(fsQ_temp(:numbreakpt_total-1,:numflux-1,:n-1)) &
            !$acc present(pQ_temp(:numflux-1,:n-1)) &
            !$acc present(mR_temp(:numsol-1,:n-1)) &
            !$acc present(dm_temp(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(mT_temp(:numsol-1,:n-1)) &
            !$acc present(sT_temp(:n-1)) &
            !$acc present(ds_temp(:numbreakpt_total-1,:n-1)) &
            !$acc present(leftbreakpt_top(:,:n-1)) &
            !$acc present(leftbreakpt_bot(:,:n-1)) &
            !$acc present(STcum_bot(:n-1)) &
            !$acc present(PQcum_bot(:numflux-1,:n-1)) &
            !$acc present(STcum_top(:n-1)) &
            !$acc present(PQcum_top(:numflux-1,:n-1)) &
            !$acc present(STcum_bot_start(:)) &
            !$acc present(STcum_top_start(:))
            do c = 0, N - 1
                jt_s = mod(c + iT_s, N)
                jt_substep = mod(jt_s, n_substeps)
                jt = (jt_s-jt_substep) / n_substeps

                pQ_aver(:, c)  = 0
                mQ_aver(:, :, c)  = 0
                mR_aver(:, c)  = 0
                if (jacobian) then
                    fs_aver(:, c)  = 0
                    fsQ_aver(:, :, c) = 0
                    fm_aver(:, :, c)  = 0
                    fmQ_aver(:, :, :, c) = 0
                    fmR_aver(:, :, c) = 0
                end if


                if ((iT_s>0).and.(c==N - iT_s)) then
                    sT_start(c) = sT_init_ts(iT_prev)
                    mT_start(:, c) = mT_init_ts(iT_prev, :)
                    if (jacobian) then
                        ds_start(:, c) = 0.
                        dm_start(:, :, c) = 0.
                    end if
                end if

                sT_temp(c) = sT_start(c)
                mT_temp(:, c) = mT_start(:, c)
                if (jacobian) then
                    ds_temp(:, c) = ds_start(:, c)
                    dm_temp(:, :, c) = dm_start(:, :, c)
                end if


                ! This is the Runge-Kutta 4th order algorithm

                do rk = 1, 5
                    hr = h * rk_time(rk)
                    if (rk>1) then
                        ! ########################## vv NEW STATE vv ##########################
                        ! Calculate the new age-ranked storage
                        sT_temp(c) = sT_start(c) ! Initial value
                        mT_temp(:, c) = mT_start(:, c) + mR_temp(:, c) * hr ! Initial value + reaction
                        ! Fluxes in & out
                        if (iT_s == 0) then
                            sT_temp(c) = sT_temp(c) + J_ts(jt) * hr / h
                            mT_temp(:, c) = mT_temp(:, c) + J_ts(jt) * C_J_ts(:, jt) * (hr/h)
                        end if
                        do iq = 0, numflux - 1
                            sT_temp(c) = sT_temp(c) - Q_ts(iq, jt) * pQ_temp(iq, c) * hr
                            mT_temp(:, c) = mT_temp(:, c) - mQ_temp(iq, :, c) * hr
                        end do
                        if (sT_temp(c)<0) then
                            !call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
                        end if
                        if (jacobian) then
                            ! Calculate new parameter sensitivity
                            do iq = 0, numflux - 1
                                ds_temp(:, c) = ds_start(:, c) - (fsQ_temp(:, iq, c)) * hr
                                dm_temp(:, :, c) = dm_start(:, :, c) - (fmQ_temp(:, :, iq, c)) * hr
                            end do
                            ds_temp(:, c) = ds_temp(:, c) - fs_temp(:, c) * hr
                            dm_temp(:, :, c) = dm_temp(:, :, c) - fm_temp(:, :, c) * hr + fmR_temp(:, :, c) * hr
                        end if
                        ! ########################## ^^ NEW STATE ^^ ##########################
                    end if
                    if (rk<5) then

                        ! ########################## vv GET FLUX vv ##########################
                        ! First get the cumulative age-ranked storage
                        if ((iT_s==0).and.(hr==0)) then
                            STcum_top(c) = 0
                            PQcum_top(:, c) = 0
                            leftbreakpt_top(:, c) = -1
                            STcum_bot(c) = 0
                            PQcum_bot(:, c) = 0
                            leftbreakpt_bot(:, c) = -1
                            pQ_temp(:, c) = 0
                        else
                            if (iT_s==0) then
                                STcum_top(c) = 0
                                STcum_bot(c) = STcum_top(c) + sT_temp(c) * hr
                            else
                                STcum_top(c) = STcum_top_start(jt_s) * (1-hr/h) + STcum_bot_start(jt_s+1) * (hr/h)
                                STcum_bot(c) = STcum_top(c) + sT_temp(c) * h
                            end if

                            if (jt_s==N) then
                                jt_this = timeseries_length - 1
                            else
                                jt_this = jt
                            end if
                            PQcum_top(:, c) = 0
                            PQcum_bot(:, c) = 0
                            do topbot = 1, 2
                                ! Main lookup loop
                                if (topbot==1) then
                                    STcum_in = STcum_top(c)
                                else
                                    STcum_in = STcum_bot(c)
                                end if
                                do iq = 0, numflux - 1
                                    do ic = component_index_list(iq), component_index_list(iq+1) - 1

                                        na = numbreakpt_list(ic)
                                        ia = 0

                                        if (STcum_in.le.SAS_lookup(breakpt_index_list(ic), jt_this)) then
                                            PQcum_component = P_list(breakpt_index_list(ic), jt_this)
                                            ia = -1
                                        else if (STcum_in.ge.SAS_lookup(breakpt_index_list(ic + 1) - 1, jt_this)) then
                                            PQcum_component = P_list(breakpt_index_list(ic + 1) - 1, jt_this)
                                            ia = na - 1
                                        else
                                            foundit = .FALSE.
                                            do i = 0, na - 1
                                                if (STcum_in.lt.SAS_lookup(breakpt_index_list(ic) + i, jt_this)) then
                                                    ia = i - 1
                                                    foundit = .TRUE.
                                                    exit
                                                endif
                                            enddo
                                            if (.not. foundit) then
                                                !call f_warning('I could not find the ST value. This should never happen!!!')
                                                PQcum_component = P_list(breakpt_index_list(ic + 1) - 1, jt_this)
                                                ia = na - 1
                                            else
                                                dif = STcum_in - SAS_lookup(breakpt_index_list(ic) + ia, jt_this)
                                                grad = (P_list(breakpt_index_list(ic) + ia + 1, jt_this) &
                                                        - P_list(breakpt_index_list(ic) + ia, jt_this)) &
                                                        / (SAS_lookup(breakpt_index_list(ic) + ia + 1, jt_this) &
                                                                - SAS_lookup(breakpt_index_list(ic) + ia, jt_this))
                                                PQcum_component = P_list(breakpt_index_list(ic) + ia, jt_this) + dif * grad
                                            endif
                                        endif

                                        if (topbot==1) then
                                            PQcum_top(iq, c) = PQcum_top(iq, c) + weights_ts(ic, jt_this) * PQcum_component
                                            leftbreakpt_top(ic, c) = ia
                                        else
                                            PQcum_bot(iq, c) = PQcum_bot(iq, c) + weights_ts(ic, jt_this) * PQcum_component
                                            leftbreakpt_bot(ic, c) = ia
                                        end if
                                    enddo
                                enddo
                            end do

                            if (iT_s==0) then
                                pQ_temp(:, c) = (PQcum_bot(:, c) - PQcum_top(:, c)) / hr
                            else
                                pQ_temp(:, c) = (PQcum_bot(:, c) - PQcum_top(:, c)) / h
                            end if
                        end if

                        do iq = 0, numflux - 1
                            if (sT_temp(c)==0) then
                                pQ_temp(iq, c) = 0
                            end if
                        end do

                        ! Solute mass flux accounting

                        do iq = 0, numflux - 1
                            do s = 0, numsol - 1

                                ! Get the mass flux out
                                if (sT_temp(c)>0) then
                                    mQ_temp(iq, s, c) = mT_temp( s, c) * alpha_ts(iq, s, jt) * Q_ts(iq, jt) &
                                            * pQ_temp(iq, c) / sT_temp(c)

                                    ! unless there is nothing in storage
                                else
                                    mQ_temp(iq, s, c) = 0.

                                end if
                            enddo
                        enddo

                        ! Reaction mass accounting
                        ! If there are first-order reactions, get the total mass rate
                        mR_temp(:, c) = k1_ts(:, jt) * (C_eq_ts(:, jt) * sT_temp(c) - mT_temp(:, c))

                        if (jacobian) then
                            fs_temp(:, c) = 0.
                            fsQ_temp(:, :, c) = 0.
                            if (sT_temp( c)>0) then
                                do iq = 0, numflux - 1
                                    fsQ_temp(:, iq, c) = fsQ_temp(:, iq, c) &
                                            + ds_temp(:, c) * pQ_temp(iq, c) * Q_ts(iq, jt) / sT_temp(c)
                                end do
                            end if
                            do iq = 0, numflux - 1
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    ! sensitivity to point before the start
                                    if ((leftbreakpt_top(ic, c)>=0).and.(leftbreakpt_top(ic, c)<numbreakpt_list(ic)-1)) then
                                        ip = breakpt_index_list(ic) + leftbreakpt_top(ic, c)
                                        !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                        dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                                        dP = P_list(ip+1, jt) - P_list(ip, jt)
                                        !call f_debug('dP/dS start    ', (/dP/dS/))
                                        fs_temp(ip, c) = fs_temp(ip, c) &
                                                + dP / (dS*dS) * sT_temp(c) * weights_ts(ic, jt) * Q_ts(iq, jt)
                                    end if
                                    ! sensitivity to point after the end
                                    if ((leftbreakpt_bot(ic, c)+1>0).and.(leftbreakpt_bot(ic, c)+1<=numbreakpt_list(ic)-1)) then
                                        ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, c) + 1
                                        !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                        dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                                        dP = P_list(ip, jt) - P_list(ip-1, jt)
                                        !call f_debug('dP/dS end      ', (/dP/dS/))
                                        fs_temp(ip, c) = fs_temp(ip, c) &
                                                - dP / (dS*dS) * sT_temp(c) * weights_ts(ic, jt) * Q_ts(iq, jt)
                                    end if
                                    ! sensitivity to point within
                                    if (leftbreakpt_bot(ic, c)>leftbreakpt_top(ic, c)) then
                                        !call f_debug('leftbreakpt_bot, _start', &
                                        !(/leftbreakpt_bot(ic, c)*one8, leftbreakpt_top(ic, c)*one8/))
                                        do leftbreakpt=leftbreakpt_top(ic, c)+1, leftbreakpt_bot(ic, c)
                                            ip = breakpt_index_list(ic) + leftbreakpt
                                            !call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                            if (leftbreakpt>0) then
                                                dSs = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                                                dPs = P_list(ip, jt) - P_list(ip-1, jt)
                                            else
                                                dSs = 1.
                                                dPs = 0.
                                            end if
                                            if (leftbreakpt<numbreakpt_list(ic)-1) then
                                                dSe = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                                                dPe = P_list(ip+1, jt) - P_list(ip, jt)
                                            else
                                                dSe = 1.
                                                dPe = 0.
                                            end if
                                            !call f_debug('dP/dS middle   ', (/dPe/dSe , dPs/dSs/))
                                            fs_temp(ip, c) = fs_temp(ip, c) &
                                                    - (dPe/dSe - dPs/dSs) / h * weights_ts(ic, jt) * Q_ts(iq, jt)
                                        end do
                                    end if
                                end do
                            end do
                            fm_temp(:, :, c) = 0
                            fmQ_temp(:, :, :, c) = 0
                            do iq = 0, numflux - 1
                                do ip = 0, numbreakpt_total - 1
                                    if (sT_temp(c)>0) then
                                        fmQ_temp(ip, iq, :, c) = fmQ_temp(ip, iq, :, c) &
                                                + dm_temp(ip, :, c) * alpha_ts(iq,:, jt) * Q_ts(iq, jt) &
                                                        * pQ_temp(iq, c) / sT_temp(c)
                                    end if
                                    fmR_temp(ip, :, c) = fmR_temp(ip, :, c) &
                                            + k1_ts(:, jt) * (C_eq_ts(:, jt) * ds_temp(ip, c) - dm_temp(ip, :, c))
                                end do
                            end do
                            do iq = 0, numflux - 1
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    ! sensitivity to point before the start
                                    if ((leftbreakpt_top(ic, c)>=0).and.(leftbreakpt_top(ic, c)<numbreakpt_list(ic)-1)) then
                                        ip = breakpt_index_list(ic) + leftbreakpt_top(ic, c)
                                        dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                                        dP = P_list(ip+1, jt) - P_list(ip, jt)
                                        fm_temp(ip, :, c) = fm_temp(ip, :, c) &
                                                + dP / (dS*dS) * mT_temp(:, c)&
                                                        * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                                    end if
                                    ! sensitivity to point after the end
                                    if ((leftbreakpt_bot(ic, c) + 1>0).and.(leftbreakpt_bot(ic, c) + 1<=numbreakpt_list(ic)-1)) then
                                        ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, c) + 1
                                        dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                                        dP = P_list(ip, jt) - P_list(ip-1, jt)
                                        fm_temp(ip, :, c) = fm_temp(ip, :, c) &
                                                - dP / (dS*dS) * mT_temp(:, c)&
                                                        * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                                    end if
                                    ! sensitivity to point within
                                    if (leftbreakpt_bot(ic, c)>leftbreakpt_top(ic, c)) then
                                        do leftbreakpt = leftbreakpt_top(ic, c)+1, leftbreakpt_bot(ic, c)
                                            ip = breakpt_index_list(ic) + leftbreakpt
                                            if (leftbreakpt>0) then
                                                dSs = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                                                dPs = P_list(ip, jt) - P_list(ip-1, jt)
                                            else
                                                dSs = 1.
                                                dPs = 0.
                                            end if
                                            if (leftbreakpt<numbreakpt_list(ic)-1) then
                                                dSe = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                                                dPe = P_list(ip+1, jt) - P_list(ip, jt)
                                            else
                                                dSe = 1.
                                                dPe = 0.
                                            end if
                                            fm_temp(ip, :, c) = fm_temp(ip, :, c) &
                                                    - (dPe/dSe - dPs/dSs) * mT_temp(:, c) / sT_temp(c) / h &
                                                            * weights_ts(ic, jt) * Q_ts(iq, jt)
                                        end do
                                    end if
                                end do
                            end do
                        end if
                        ! ########################## ^^ GET FLUX ^^ ##########################

                        pQ_aver(:, c)  = pQ_aver(:, c)  + rk_coeff(rk) * pQ_temp(:, c) / 6.
                        mQ_aver(:, :, c)  = mQ_aver(:, :, c)  + rk_coeff(rk) * mQ_temp(:, :, c) / 6.
                        mR_aver(:, c)  = mR_aver(:, c)  + rk_coeff(rk) * mR_temp(:, c) / 6.
                        if (jacobian) then
                            fs_aver(:, c)  = fs_aver(:, c)  + rk_coeff(rk) * fs_temp(:, c) / 6.
                            fsQ_aver(:, :, c) = fsQ_aver(:, :, c) + rk_coeff(rk) * fsQ_temp(:, :, c) / 6.
                            fm_aver(:, :, c)  = fm_aver(:, :, c)  + rk_coeff(rk) * fm_temp(:, :, c) / 6.
                            fmQ_aver(:, :, :, c) = fmQ_aver(:, :, :, c) + rk_coeff(rk) * fmQ_temp(:, :, :, c) / 6.
                            fmR_aver(:, :, c) = fmR_aver(:, :, c) + rk_coeff(rk) * fmR_temp(:, :, c) / 6.
                        end if

                    end if
                    if (rk==4) then
                        ! zero out the probabilities if there is no outflux this timestep
                        do iq = 0, numflux - 1
                            if (Q_ts(iq, jt)==0) then
                                pQ_aver(iq, c) = 0.
                                mQ_aver(iq, :, c) = 0.
                            end if
                        end do
                        pQ_temp(:, c)  = pQ_aver(:, c)
                        mQ_temp(:, :, c)  = mQ_aver(:, :, c)
                        mR_temp(:, c)  = mR_aver(:, c)
                        if (jacobian) then
                            fs_temp(:, c)  = fs_aver(:, c)
                            fsQ_temp(:, :, c) = fsQ_aver(:, :, c)
                            fm_temp(:, :, c)  = fm_aver(:, :, c)
                            fmQ_temp(:, :, :, c) = fmQ_aver(:, :, :, c)
                            fmR_temp(:, :, c) = fmR_aver(:, :, c)
                        end if
                    end if
                end do

                ! Update the state with the new estimates
                sT_start(c) = sT_temp(c)
                mT_start(:, c) = mT_temp(:, c)
                if (jacobian) then
                    ds_start(:, c) = ds_temp(:, c)
                    dm_start(:, :, c) = dm_temp(:, :, c)
                end if
            end do
            call cpu_time(finish)
            !print '("Calc time for ",i6," steps = ",f6.3," milliseconds.")',N,1000*(finish-start)
            ! Aggregate data from substep to timestep
            !$acc kernels &
            !$acc present(fmR_aver(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(fm_aver(:,:numsol-1,:n-1)) &
            !$acc present(fmQ_aver(:numbreakpt_total-1,:numflux-1,:numsol-1,:n-1)) &
            !$acc present(fsQ_aver(:numbreakpt_total-1,:numflux-1,:n-1)) &
            !$acc present(mR_aver(:numsol-1,:n-1)) &
            !$acc present(pQ_aver(:numflux-1,:n-1)) &
            !$acc present(mQ_aver(:numflux-1,:numsol-1,:n-1)) &
            !$acc present(fs_aver(:,:n-1)) &
            !$acc present(dm_start(:numbreakpt_total-1,:numsol-1,:n-1)) &
            !$acc present(ds_start(:numbreakpt_total-1,:n-1)) &
            !$acc present(sT_start(:n-1)) &
            !$acc present(mT_start(:numsol-1,:n-1)) &
            !$acc present(sT_ts) &
            !$acc present(mT_ts) &
            !$acc present(ds_ts) &
            !$acc present(dm_ts) &
            !$acc present(dC_ts) &
            !$acc present(dW_ts) &
            !$acc present(pQ_ts) &
            !$acc present(mQ_ts) &
            !$acc present(mR_ts) &
            !$acc present(STcum_bot_start(:)) &
            !$acc present(STcum_top_start(:))
            STcum_top_start(0) = STcum_bot_start(0)
            STcum_bot_start(0) = STcum_bot_start(0) + sT_init_ts(iT) * h
            do jt_substep = 0, n_substeps - 1
                !$acc loop independent
                do jt = 0, timeseries_length - 1
                    jt_s = jt * n_substeps +  jt_substep
                    c = mod(N + jt_s - iT_s, N)

                    STcum_top_start(jt_s+1) = STcum_bot_start(jt_s+1)
                    STcum_bot_start(jt_s+1) = STcum_bot_start(jt_s+1) + sT_start(c) * h


                    ! Get the timestep-averaged transit time distribution
                    !carry = ((iT<max_age-1).and.(jt_substep<iT_substep)) ! This gives incorrect result with PGI compiler
                    if ((iT<max_age-1).and.(jt_substep<iT_substep)) then
                        carry = 1
                    else
                        carry = 0
                    end if

                    pQ_ts(iT+carry, jt, :) =    pQ_ts(iT+carry, jt, :) + pQ_aver(:, c) * norm
                    mQ_ts(iT+carry, jt, :, :) = mQ_ts(iT+carry, jt, :, :) + mQ_aver(:, :, c) * norm
                    mR_ts(iT+carry, jt, :) =    mR_ts(iT+carry, jt, :) + mR_aver(:, c) * norm

                    if (jacobian) then
                        do iq = 0, numflux - 1
                            if (Q_ts(iq, jt)>0) then
                                dW_ts(:, iq, jt) = dW_ts(:, iq, jt) + fsQ_aver(:, iq, c) / Q_ts(iq, jt) * norm * dt
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                        dW_ts(ip, iq, jt) = dW_ts(ip, iq, jt) + fs_aver(ip, c) / Q_ts(iq, jt) * norm * dt
                                    enddo
                                enddo
                                dC_ts(:, iq, :, jt) = dC_ts(:, iq, :, jt) &
                                        + fmQ_aver(:, iq, :, c) / Q_ts(iq, jt) * norm * dt
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                        dC_ts(ip, iq, :, jt) = dC_ts(ip, iq, :, jt) &
                                                + fm_aver(ip, :, c) / Q_ts(iq, jt) * norm * dt
                                    enddo
                                enddo
                            end if
                        enddo
                    end if


                    ! Extract substep state at timesteps
                    ! age-ranked storage at the end of the timestep
                    if (jt_substep==n_substeps-1) then
                        sT_ts(iT,jt+1) = sT_ts(iT,jt+1) + sT_start(c) / n_substeps
                        ! parameter sensitivity
                        if (jacobian) then
                            do ip = 0, numbreakpt_total - 1
                                ds_ts(iT,jt+1, ip) = ds_ts(iT,jt+1, ip) + ds_start(ip, c) / n_substeps
                            enddo
                        end if
                        ! Age-ranked solute mass
                        do s = 0, numsol - 1
                            mT_ts(iT,jt+1, s) = mT_ts(iT,jt+1, s) + mT_start( s, c) / n_substeps
                            ! parameter sensitivity
                            if (jacobian) then
                                do ip = 0, numbreakpt_total - 1
                                    dm_ts(iT,jt+1, ip, s) = dm_ts(iT,jt+1, ip, s) + dm_start(ip, s, c) / n_substeps
                                enddo
                            end if
                        enddo
                    end if
                enddo
            enddo
            !$acc end kernels

            !call f_debug('sT_ts(iT, :)     ', sT_ts(iT, :))
            !call f_debug('pQ_aver0         ', pQ_aver(0,:))
            !call f_debug('pQ_aver1         ', pQ_aver(1,:))
            !call f_debug('pQ_ts(iT, :, 0)'  , pQ_ts(iT, :, 0))
            !call f_debug('pQ_ts(iT, :, 0)'  , pQ_ts(iT, :, 1))

            iT_prev = iT

            call cpu_time(finish2)
            !print '("Stor time for ",i6," steps = ",f6.3," milliseconds.")',N,1000*(finish2-finish)

        enddo
    enddo
    !$acc update self(sT_ts)
    !$acc update self(mT_ts)
    !$acc update self(ds_ts)
    !$acc update self(dm_ts)
    !$acc update self(dC_ts)
    !$acc update self(dW_ts)
    !$acc update self(pQ_ts)
    !$acc update self(mQ_ts)
    !$acc update self(mR_ts)

    ! Calculate a water balance
    ! Difference of starting and ending age-ranked storage
    do iT = 0, max_age - 1
        do jt = 0, timeseries_length - 1
            if (iT==0) then
                WaterBalance_ts(iT, jt) = J_ts(jt) - sT_ts(iT, jt+1)
            else
                WaterBalance_ts(iT, jt) = sT_ts(iT-1, jt) - sT_ts(iT, jt+1)
            end if
            ! subtract time-averaged water fluxes
            do iq = 0, numflux - 1
                WaterBalance_ts(iT, jt) = WaterBalance_ts(iT, jt) - (Q_ts(iq, jt) * pQ_ts(iT, jt, iq)) * dt
            end do

            ! Calculate a solute balance
            ! Difference of starting and ending age-ranked mass
            if (iT==0) then
                do s = 0, numsol - 1
                    SoluteBalance_ts(iT, jt, s) = C_J_ts( s, jt) * J_ts(jt) - mT_ts(iT, jt+1, s)
                end do
            else
                SoluteBalance_ts(iT, jt, :) = mT_ts(iT-1, jt, :) - mT_ts(iT, jt+1, :)
            end if
            ! Subtract timestep-averaged mass fluxes
            do iq = 0, numflux - 1
                SoluteBalance_ts(iT, jt, :) = SoluteBalance_ts(iT, jt, :) - (mQ_ts(iT, jt, iq, :)) * dt
            end do
            ! Reacted mass
            SoluteBalance_ts(iT, jt, :) = SoluteBalance_ts(iT, jt, :) + mR_ts(iT, jt, :) * dt
        enddo

        ! Print some updates
        if (mod(iT, 10).eq.0) then
            write (tempdebugstring, *) '...Done ', (iT), &
                    'of', (timeseries_length)
            !call f_verbose(tempdebugstring)
        endif

    enddo ! End of main loop
    !$acc end data region

    !call f_verbose('...Finalizing...')

    ! get the old water fraction
    P_old = 1 - transpose(sum(pQ_ts, DIM=1)) * dt

    do s = 0, numsol - 1
        do iq = 0, numflux - 1

            where (Q_ts(iq, :)>0)
                ! From the age-ranked mass
                C_Q_ts(:, iq, s) = sum(mQ_ts(:, :, iq, s), DIM=1) / Q_ts(iq, :) * dt

                ! From the old water concentration
                C_Q_ts(:, iq, s) = C_Q_ts(:, iq, s) + alpha_ts( iq, s, :) * C_old(s) * P_old(iq, :)

            end where

            if (jacobian) then
                do ip = 0, numbreakpt_total - 1
                    where (Q_ts(iq, :)>0)
                        dC_ts(ip, iq, s, :) = dC_ts(ip, iq, s, :) - C_old(s) * dW_ts(ip, iq, :)
                    end where
                end do
            end if


        enddo
    enddo


    !call f_verbose('...Finished...')

contains


    !subroutine f_debug_blank()
        !! Prints a blank line
        !!$acc routine seq
        !if (debug) then
            !print *, ''
        !endif
    !end subroutine f_debug_blank


    !subroutine f_debug(debugstring, debugdblepr)
        !! Prints debugging information
        !implicit none
        !!$acc routine seq
        !character(len = *), intent(in) :: debugstring
        !real(8), dimension(:), intent(in) :: diebugdblepr
        !if (debug) then
            !print 1, debugstring, diebugdblepr
            !1 format (A26, *(f16.10))
        !endif
    !end subroutine f_debug


    !subroutine f_warning(debugstring)
        !! Prints informative information
        !implicit none
        !!$acc routine seq
        !character(len = *), intent(in) :: debugstring
        !if (warning) then
            !print *, debugstring
        !endif
    !end subroutine f_warning


    !subroutine f_verbose(debugstring)
        !! Prints informative information
        !implicit none
        !!$acc routine seq
        !character(len = *), intent(in) :: debugstring
        !if (verbose) then
            !print *, debugstring
        !endif
    !end subroutine f_verbose

end subroutine solve
