! -*- f90 -*-
subroutine solveSAS(J_ts, Q_ts, SAS_args, SAS_lookup, P_list, weights_ts, sT_init_ts, dt, &
        verbose, debug, warning, jacobian,&
        mT_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
        n_substeps, numcomponent_list, numargs_list, numbreakpt_list, numflux, numsol, max_age, &
        timeseries_length, numcomponent_total, numargs_total, numbreakpt_total, &
        sT_ts, pQ_ts, WaterBalance_ts, &
        mT_ts, mQ_ts, mR_ts, C_Q_ts, ds_ts, dm_ts, dC_ts, SoluteBalance_ts)
    use cdf_gamma_mod
    !use cdf_beta_mod
    !use cdf_normal_mod
    implicit none

    ! Start by declaring and initializing all the variables we will be using
    integer, intent(in) :: n_substeps, numflux, numsol, max_age, &
            timeseries_length, numcomponent_total, numargs_total, numbreakpt_total
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning, jacobian
    real(8), intent(in), dimension(0:timeseries_length - 1) :: J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1) :: Q_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numcomponent_total - 1) :: weights_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numargs_total - 1) :: SAS_args
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numbreakpt_total - 1) :: SAS_lookup
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numbreakpt_total - 1) :: P_list
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: k1_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_eq_ts
    real(8), intent(in), dimension(0:numsol - 1) :: C_old
    real(8), intent(in), dimension(0:max_age - 1) :: sT_init_ts
    real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mT_init_ts
    integer, intent(in), dimension(0:numflux - 1) :: numcomponent_list
    integer, intent(in), dimension(0:numcomponent_total - 1) :: numbreakpt_list
    integer, intent(in), dimension(0:numcomponent_total - 1) :: numargs_list
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: C_Q_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numbreakpt_total-1, 0:numflux - 1, 0:numsol - 1) :: dC_ts
    real(8), intent(out), dimension(0:timeseries_length, 0:max_age - 1) :: sT_ts
    real(8), intent(out), dimension(0:timeseries_length, 0:numsol - 1, 0:max_age - 1) :: mT_ts
    real(8), intent(out), dimension(0:timeseries_length, 0:numbreakpt_total-1, 0:max_age - 1) :: ds_ts
    real(8), intent(out), dimension(0:timeseries_length, 0:numbreakpt_total-1, 0:numsol - 1, 0:max_age - 1) :: dm_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:max_age - 1) :: pQ_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1, 0:max_age - 1) :: mQ_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numsol - 1, 0:max_age - 1) :: mR_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:max_age - 1) :: WaterBalance_ts
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numsol - 1, 0:max_age - 1) :: SoluteBalance_ts
    real(8), dimension(0:timeseries_length-1, 0:numbreakpt_total-1, 0:numflux - 1) :: dW_ts
    real(8), dimension(0:timeseries_length - 1, 0:numflux - 1) :: P_old
    integer, dimension(0:numcomponent_total) :: breakpt_index_list
    integer, dimension(0:numcomponent_total) :: args_index_list
    integer, dimension(0:numcomponent_total) :: component_type
    integer, dimension(0:numflux) :: component_index_list
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_top_start
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_bot_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_top
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_top
    integer, dimension(0:timeseries_length * n_substeps - 1, 0:numcomponent_total-1) :: leftbreakpt_bot
    integer, dimension(0:timeseries_length * n_substeps - 1, 0:numcomponent_total-1) :: leftbreakpt_top
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1) :: fs_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1) :: fs_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numflux - 1) :: fsQ_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numflux - 1) :: fsQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numsol - 1) :: fm_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numsol - 1) :: fm_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numflux-1, 0:numsol - 1) :: fmQ_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numflux-1, 0:numsol - 1) :: fmQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numsol - 1) :: fmR_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total - 1, 0:numsol - 1) :: fmR_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mT_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1) :: ds_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1) :: ds_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_in
    integer, dimension(0:timeseries_length * n_substeps - 1) :: jt
    integer, dimension(0:timeseries_length * n_substeps - 1) :: jt_s
    integer :: iT_substep, iT, iT_s, iT_prev, jt_substep, jt_i
    real(8) :: one8, norm
    real(8) :: dS, dP, dSe, dPe, dSs, dPs
    real(8) :: h, hr
    real(8), dimension(4) :: rk_coeff
    real(8), dimension(5) :: rk_time
    character(len = 128) :: tempdebugstring
    integer :: iq, s, M, N, ip, ic, c, rk
    integer :: carry
    integer :: leftbreakpt
    real(8) :: PQcum_component, X, scale, loc, a_arg, b_arg
    integer :: CDFLIB90_STATUS
    logical :: CDFLIB90_CHECKINPUT
    integer :: jt_this, topbot
    integer :: na
    integer :: ia
    integer :: i
    real(8) :: dif, grad
    logical :: foundit
    real(8) :: start, finish
    real(8), dimension(37) :: runtime
    CDFLIB90_CHECKINPUT = .true.
    runtime = 0.

    !$acc data &
    !$acc copyin(component_index_list) &
    !$acc copyin(rk_coeff) &
    !$acc copyin(SAS_args) &
    !$acc copyin(SAS_lookup) &
    !$acc copyin(J_ts) &
    !$acc copyin(breakpt_index_list) &
    !$acc copyin(args_index_list) &
    !$acc copyin(k1_ts) &
    !$acc copyin(C_J_ts) &
    !$acc copyin(weights_ts) &
    !$acc copyin(alpha_ts) &
    !$acc copyin(rk_time) &
    !$acc copyin(Q_ts) &
    !$acc copyin(P_list) &
    !$acc copyin(C_eq_ts) &
    !$acc copyin(numbreakpt_list) &
    !$acc copyin(numargs_list) &
    !$acc copyin(jacobian) &
    !$acc copyin(sT_init_ts,mT_init_ts) &
    !$acc copyin(dm_start) &
    !$acc copyin(ds_start) &
    !$acc copyin(sT_start) &
    !$acc copyin(mT_start) &
    !$acc create(fmR_aver) &
    !$acc create(fm_aver) &
    !$acc create(fmQ_aver) &
    !$acc create(fsQ_aver) &
    !$acc create(mR_aver) &
    !$acc create(pQ_aver) &
    !$acc create(mQ_aver) &
    !$acc create(fs_aver) &
    !$acc create(fmQ_temp) &
    !$acc create(mQ_temp) &
    !$acc create(fmR_temp) &
    !$acc create(fs_temp) &
    !$acc create(fm_temp) &
    !$acc create(fsQ_temp) &
    !$acc create(pQ_temp) &
    !$acc create(mR_temp) &
    !$acc create(dm_temp) &
    !$acc create(mT_temp) &
    !$acc create(sT_temp) &
    !$acc create(ds_temp) &
    !$acc create(leftbreakpt_top) &
    !$acc create(leftbreakpt_bot) &
    !$acc create(STcum_bot) &
    !$acc create(PQcum_bot) &
    !$acc create(STcum_top) &
    !$acc create(PQcum_top) &
    !$acc copyin(STcum_bot_start) &
    !$acc copyin(STcum_top_start) &
    !$acc copyin(dW_ts) &
    !$acc copyin(dC_ts) &
    !$acc copyin(sT_ts) &
    !$acc copyin(mT_ts) &
    !$acc copyin(ds_ts) &
    !$acc copyin(dm_ts) &
    !$acc copyin(dC_ts) &
    !$acc copyin(dW_ts) &
    !$acc copyin(pQ_ts) &
    !$acc copyin(mQ_ts) &
    !$acc copyin(mR_ts) &
    !$acc create(STcum_in) &
    !$acc create(jt) &
    !$acc create(jt_s) &
    !$acc copyin(n_substeps, N) &
    !$acc create(iT, iT_substep, iT_s, hr, rk)
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
    args_index_list = 0
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


    call f_verbose('...Initializing arrays...')
    one8 = 1.0
    rk_time = (/0.0D0, 0.5D0, 0.5D0, 1.0D0, 1.0D0/)
    rk_coeff = (/1./6, 2./6, 2./6, 1./6/)
    norm = 1.0 / n_substeps / n_substeps


    ! The list of probabilities in each sas function is a 1-D array.
    ! breakpt_index_list gives the starting index of the probabilities (P) associated
    ! with each flux
    breakpt_index_list(0) = 0
    args_index_list(0) = 0
    component_index_list(0) = 0
    do iq = 0, numflux - 1
        component_index_list(iq + 1) = component_index_list(iq) + numcomponent_list(iq)
        do ic = component_index_list(iq), component_index_list(iq+1) - 1
            breakpt_index_list(ic + 1) = breakpt_index_list(ic) + numbreakpt_list(ic)
            args_index_list(ic + 1) = args_index_list(ic) + numargs_list(ic)
            component_type(ic) = int(SAS_args(0, args_index_list(ic)))
        enddo
    enddo
    call f_debug('breakpt_index_list', one8 * breakpt_index_list(:))

    ! modify the number of ages and the timestep by a facotr of n_substeps
    M = max_age * n_substeps
    N = timeseries_length * n_substeps
    h = dt / n_substeps

    call f_verbose('...Setting initial conditions...')
    sT_ts(0, :) = sT_init_ts
    do s = 0, numsol - 1
        mT_ts(0, s, :) = mT_init_ts(:, s)
    end do

    !$acc update device(component_index_list)
    !$acc update device(rk_coeff)
    !$acc update device(SAS_args)
    !$acc update device(SAS_lookup)
    !$acc update device(J_ts)
    !$acc update device(breakpt_index_list)
    !$acc update device(k1_ts)
    !$acc update device(C_J_ts)
    !$acc update device(weights_ts)
    !$acc update device(alpha_ts)
    !$acc update device(rk_time)
    !$acc update device(Q_ts)
    !$acc update device(P_list)
    !$acc update device(C_eq_ts)
    !$acc update device(numbreakpt_list)
    !$acc update device(jacobian)
    !$acc update device(sT_init_ts,mT_init_ts)
    !$acc update device(dm_start)
    !$acc update device(ds_start)
    !$acc update device(sT_start)
    !$acc update device(mT_start)
    !$acc update device(fmR_aver)
    !$acc update device(fm_aver)
    !$acc update device(fmQ_aver)
    !$acc update device(fsQ_aver)
    !$acc update device(mR_aver)
    !$acc update device(pQ_aver)
    !$acc update device(mQ_aver)
    !$acc update device(fs_aver)
    !$acc update device(fmQ_temp)
    !$acc update device(mQ_temp)
    !$acc update device(fmR_temp)
    !$acc update device(fs_temp)
    !$acc update device(fm_temp)
    !$acc update device(fsQ_temp)
    !$acc update device(pQ_temp)
    !$acc update device(mR_temp)
    !$acc update device(dm_temp)
    !$acc update device(mT_temp)
    !$acc update device(sT_temp)
    !$acc update device(ds_temp)
    !$acc update device(leftbreakpt_top)
    !$acc update device(leftbreakpt_bot)
    !$acc update device(STcum_bot)
    !$acc update device(PQcum_bot)
    !$acc update device(STcum_top)
    !$acc update device(PQcum_top)
    !$acc update device(STcum_bot_start)
    !$acc update device(STcum_top_start)
    !$acc update device(dW_ts)
    !$acc update device(dC_ts)
    !$acc update device(sT_ts)
    !$acc update device(mT_ts)
    !$acc update device(ds_ts)
    !$acc update device(dm_ts)
    !$acc update device(dC_ts)
    !$acc update device(dW_ts)
    !$acc update device(pQ_ts)
    !$acc update device(mQ_ts)
    !$acc update device(mR_ts)
    !$acc update device(STcum_in)
    !$acc update device(jt)
    !$acc update device(jt_s)
    !$acc update device(n_substeps, N)
    !$acc update device(iT, iT_substep, iT_s, hr, rk)

    call f_verbose('...Starting main loop...')
    do iT = 0, max_age - 1

        ! Start the substep loop
        do iT_substep = 0, n_substeps - 1

            iT_s = iT * n_substeps +  iT_substep
            !$acc update device(iT, iT_substep, iT_s)

            call cpu_time(start)
            !$acc kernels
            !$acc loop independent
            do c = 0, N - 1
                jt_s(c) = mod(c + iT_s, N)
                jt_substep = mod(jt_s(c), n_substeps)
                jt(c) = (jt_s(c)-jt_substep) / n_substeps
            end do
            !$acc end kernels
            !$acc update self(jt)
            call cpu_time(finish)
            runtime(1) = runtime(1) + 1000*(finish-start)

            call cpu_time(start)
            !$acc kernels

            pQ_aver  = 0
            mQ_aver  = 0
            mR_aver  = 0
            if (iT_s>0) then
                sT_start(N - iT_s) = sT_init_ts(iT_prev)
                mT_start(N - iT_s, :) = mT_init_ts(iT_prev, :)
            end if

            !$acc loop independent
            do c = 0, N - 1
                sT_temp(c) = sT_start(c)
            end do
            !$acc loop independent
            do s = 0, numsol - 1
                !$acc loop independent
                do c = 0, N - 1
                    mT_temp(c, s) = mT_start(c, s)
                end do
            end do

            !$acc end kernels
            call cpu_time(finish)
            runtime(2) = runtime(2) + 1000*(finish-start)

            if (jacobian) then
                call cpu_time(start)
                !$acc kernels
                fs_aver  = 0
                fsQ_aver = 0
                fm_aver  = 0
                fmQ_aver = 0
                fmR_aver = 0
                !$acc end kernels
                call cpu_time(finish)
                runtime(3) = runtime(3) + 1000*(finish-start)
                if (iT_s>0) then
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do ip = 0, numbreakpt_total - 1
                        ds_start( N - iT_s, ip) = 0.
                    end do
                    !$acc loop independent
                    do ip = 0, numbreakpt_total - 1
                        !$acc loop independent
                        do s = 0, numsol - 1
                            dm_start( N - iT_s, ip, s) = 0.
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(4) = runtime(4) + 1000*(finish-start)
                end if
                call cpu_time(start)
                !$acc kernels
                !$acc loop independent
                do ip = 0, numbreakpt_total - 1
                    !$acc loop independent
                    do c = 0, N - 1
                        ds_temp(c, ip) = ds_start(c, ip)
                    end do
                end do
                !$acc loop independent
                do s = 0, numsol - 1
                    !$acc loop independent
                    do ip = 0, numbreakpt_total - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            dm_temp(c, ip, s) = dm_start(c, ip, s)
                        end do
                    end do
                end do
                !$acc end kernels
                call cpu_time(finish)
                runtime(5) = runtime(5) + 1000*(finish-start)
            end if


            ! This is the Runge-Kutta 4th order algorithm

            do rk = 1, 5
                hr = h * rk_time(rk)
                !$acc update device(hr, rk)
                if (rk>1) then


                ! ########################## vv NEW STATE vv ##########################
                ! Calculate the new age-ranked storage
                !$acc update self(pQ_temp, sT_temp)
                call f_debug('NEW STATE rk           ', (/rk*one8, iT_s*one8/))
                call f_debug('pQ_temp                ', pQ_temp(:, 0))
                call f_debug('sT_temp 0              ', sT_temp(:))
                call cpu_time(start)
                !$acc kernels
                !$acc loop independent
                do c = 0, N - 1
                    sT_temp(c) = sT_start(c) ! Initial value
                end do
                !$acc loop independent
                do s = 0, numsol - 1
                    !$acc loop independent
                    do c = 0, N - 1
                        mT_temp(c, s) = mT_start(c, s) + mR_temp(c, s) * hr ! Initial value + reaction
                    end do
                end do
                !$acc end kernels
                call cpu_time(finish)
                !$acc update self(pQ_temp, sT_temp)
                call f_debug('sT_temp 1              ', sT_temp(:))
                runtime(6) = runtime(6) + 1000*(finish-start)
                call cpu_time(start)
                ! Fluxes in & out
                if (iT_s == 0) then
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do c = 0, N - 1
                        sT_temp(c) = sT_temp(c) + J_ts(jt(c)) * hr / h
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            mT_temp(c, s) = mT_temp(c, s) + J_ts(jt(c)) * C_J_ts(jt(c), s) * (hr/h)
                        end do
                    end do
                    !$acc end kernels
                    !$acc update self(pQ_temp, sT_temp)
                    call f_debug('sT_temp 2              ', sT_temp(:))
                    call cpu_time(finish)
                    runtime(7) = runtime(7) + 1000*(finish-start)
                end if
                call cpu_time(start)
                !$acc kernels
                !$acc loop independent
                do iq = 0, numflux - 1
                    !$acc loop independent
                    do c = 0, N - 1
                        sT_temp(c) = sT_temp(c) - Q_ts(jt(c), iq) * pQ_temp(c, iq) * hr
                        if (sT_temp(c)<0) then
                            call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
                        end if
                    end do
                end do
                !$acc loop independent
                do s = 0, numsol - 1
                    !$acc loop independent
                    do iq = 0, numflux - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            mT_temp(c, s) = mT_temp(c, s) - mQ_temp( c, iq, s) * hr
                        end do
                    end do
                end do
                !$acc end kernels
                call cpu_time(finish)
                call f_debug('sT_temp 3              ', sT_temp(:))
                !$acc update self(pQ_temp, sT_temp)
                runtime(8) = runtime(8) + 1000*(finish-start)
                if (jacobian) then
                    ! Calculate new parameter sensitivity
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do iq = 0, numflux - 1
                        !$acc loop independent
                        do ip = 0, numbreakpt_total - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                ds_temp(c, ip) = ds_start(c, ip) - fsQ_temp(c, ip, iq) * hr
                            end do
                        end do
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do iq = 0, numflux - 1
                            !$acc loop independent
                            do ip = 0, numbreakpt_total - 1
                                !$acc loop independent
                                do c = 0, N - 1
                                    dm_temp(c, ip, s) = dm_start(c, ip, s) - fmQ_temp(c, ip, iq, s) * hr
                                end do
                            end do
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(9) = runtime(9) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do ip = 0, numbreakpt_total - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            ds_temp(c, ip) = ds_temp(c, ip) - fs_temp(c, ip) * hr
                        end do
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do ip = 0, numbreakpt_total - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                dm_temp(c, ip, s) = dm_temp(c, ip, s) - fm_temp(c, ip, s) * hr - fmR_temp(c, ip, s) * hr
                            end do
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(10) = runtime(10) + 1000*(finish-start)
                end if
                ! ########################## ^^ NEW STATE ^^ ##########################
            end if
                if (rk<5) then
                !$acc update self(pQ_temp, sT_temp)
                call f_debug('GET FLUX  rk           ', (/rk*one8, iT_s*one8/))
                call f_debug('sT_temp                ', sT_temp(:))
                call f_debug('pQ_temp start          ', pQ_temp(:, 0))

                ! ########################## vv GET FLUX vv ##########################
                ! First get the cumulative age-ranked storage
                if ((iT_s==0).and.(hr==0)) then
                    call cpu_time(start)
                    !$acc kernels
                    STcum_top = 0
                    PQcum_top = 0
                    leftbreakpt_top = -1
                    STcum_bot = 0
                    PQcum_bot = 0
                    leftbreakpt_bot = -1
                    pQ_temp = 0
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(11) = runtime(11) + 1000*(finish-start)
                else
                    if (iT_s==0) then
                        call cpu_time(start)
                        !$acc kernels
                        STcum_top = 0
                        !$acc loop independent
                        do c = 0, N - 1
                            STcum_bot(c) = STcum_top(c) + sT_temp(c) * hr
                        end do
                        !$acc end kernels
                        call cpu_time(finish)
                        runtime(12) = runtime(12) + 1000*(finish-start)
                    else
                        call cpu_time(start)
                        !$acc kernels
                        !$acc loop independent
                        do c = 0, N - 1
                            STcum_top(c) = STcum_top_start(jt_s(c)) * (1-hr/h) + STcum_bot_start(jt_s(c)+1) * (hr/h)
                            STcum_bot(c) = STcum_top(c) + sT_temp(c) * h
                        end do
                        !$acc end kernels
                        call cpu_time(finish)
                        runtime(13) = runtime(13) + 1000*(finish-start)
                    end if

                    call cpu_time(start)
                    !$acc kernels
                    PQcum_top = 0
                    PQcum_bot = 0
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(14) = runtime(14) + 1000*(finish-start)
                    do topbot = 0, 1
                        ! Main lookup loop
                        if (topbot==0) then
                            call cpu_time(start)
                            !$acc kernels
                            !$acc loop independent
                            do c = 0, N - 1
                                STcum_in(c) = STcum_top(c)
                            end do
                            !$acc end kernels
                            call cpu_time(finish)
                            runtime(15) = runtime(15) + 1000*(finish-start)
                        else
                            call cpu_time(start)
                            !$acc kernels
                            !$acc loop independent
                            do c = 0, N - 1
                            STcum_in(c) = STcum_bot(c)
                            end do
                            !$acc end kernels
                            call cpu_time(finish)
                            runtime(16) = runtime(16) + 1000*(finish-start)
                        end if
                        call cpu_time(start)
                        do iq = 0, numflux - 1
                            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                if (component_type(ic)==-1) then
                                    !$acc kernels
                                    !$acc loop independent
                                    do c = 0, N - 1
                                        if (STcum_in(c).le.SAS_lookup( jt(c), breakpt_index_list(ic))) then
                                            if (topbot==0) then
                                                PQcum_component = P_list( jt(c), breakpt_index_list(ic))
                                                PQcum_top(c, iq) = PQcum_top(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                                leftbreakpt_top(c, ic) = -1
                                            else
                                                PQcum_component = P_list( jt(c), breakpt_index_list(ic))
                                                PQcum_bot(c, iq) = PQcum_bot(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                                leftbreakpt_bot(c, ic) = -1
                                            end if
                                        else
                                            na = numbreakpt_list(ic)
                                            ia = 0
                                            foundit = .FALSE.
                                            do i = 0, na - 1
                                                if (STcum_in(c).lt.SAS_lookup( jt(c), breakpt_index_list(ic) + i)) then
                                                    ia = i - 1
                                                    foundit = .TRUE.
                                                    exit
                                                endif
                                            enddo
                                            if (.not. foundit) then
                                                call f_warning('I could not find the ST value. This should never happen!!!')
                                                PQcum_component = P_list( jt(c), breakpt_index_list(ic + 1) - 1)
                                                ia = na - 1
                                            else
                                                dif = STcum_in(c) - SAS_lookup( jt(c), breakpt_index_list(ic) + ia)
                                                grad = (P_list( jt(c), breakpt_index_list(ic) + ia + 1) &
                                                - P_list( jt(c), breakpt_index_list(ic) + ia)) &
                                                / (SAS_lookup( jt(c), breakpt_index_list(ic) + ia + 1) &
                                                - SAS_lookup( jt(c), breakpt_index_list(ic) + ia))
                                                PQcum_component = P_list( jt(c), breakpt_index_list(ic) + ia) + dif * grad
                                            endif
                                            if (topbot==0) then
                                                PQcum_top(c, iq) = PQcum_top(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                                leftbreakpt_top(c, ic) = ia
                                            else
                                                PQcum_bot(c, iq) = PQcum_bot(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                                leftbreakpt_bot(c, ic) = ia
                                            end if
                                        endif
                                    enddo
                                    !$acc end kernels
                                elseif (component_type(ic)==1) then
                                    !print *, iT_s, hr
                                    ! Gamma distribution
                                    !$acc kernels
                                    !$acc loop independent
                                    do c = 0, N - 1
                                        loc = SAS_args(jt(c), args_index_list(ic) + 1)
                                        scale = SAS_args(jt(c), args_index_list(ic) + 2)
                                        a_arg = SAS_args(jt(c), args_index_list(ic) + 3)
                                        X = (STcum_in(c) - loc) / scale
                                        PQcum_component = 0
                                        !print *, c, jt(c), loc, scale, a_arg
                                        if (X.gt.0) then
                                            PQcum_component = CUM_GAMMA(X, a_arg, one8, CDFLIB90_STATUS, CDFLIB90_CHECKINPUT)
                                        end if
                                        !print *, PQcum_component, CDFLIB90_STATUS
                                        !if (CDFLIB90_STATUS.ne.0) then
                                        !    if (warning) then
                                        !        print *, "CDFLIB90_STATUS=", CDFLIB90_STATUS
                                        !    endif
                                        !end if
                                        if (topbot==0) then
                                            PQcum_top(c, iq) = PQcum_top(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                        else
                                            PQcum_bot(c, iq) = PQcum_bot(c, iq) + weights_ts( jt(c), ic) * PQcum_component
                                        end if
                                    end do
                                    !$acc end kernels
                                endif
                            enddo
                        end do
                        call cpu_time(finish)
                        runtime(17) = runtime(17) + 1000*(finish-start)
                    end do

                    if (iT_s==0) then
                        call cpu_time(start)
                        !$acc kernels
                        !$acc loop independent
                        do iq = 0, numflux - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                pQ_temp(c, iq) = (PQcum_bot(c, iq) - PQcum_top(c, iq)) / hr
                            end do
                        end do
                        !$acc end kernels
                        call cpu_time(finish)
                        runtime(18) = runtime(18) + 1000*(finish-start)
                    else
                        call cpu_time(start)
                        !$acc kernels
                        !$acc loop independent
                        do iq = 0, numflux - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                pQ_temp(c, iq) = (PQcum_bot(c, iq) - PQcum_top(c, iq)) / h
                            end do
                        end do
                        !$acc end kernels
                        call cpu_time(finish)
                        runtime(19) = runtime(19) + 1000*(finish-start)
                    end if
                end if

                call cpu_time(start)
                !$acc kernels
                do iq = 0, numflux - 1
                    where (sT_temp==0)
                        pQ_temp( :, iq) = 0
                    end where
                end do
                !$acc end kernels
                call cpu_time(finish)
                runtime(20) = runtime(20) + 1000*(finish-start)

                ! Solute mass flux accounting

                ! Get the mass flux out
                call cpu_time(start)
                !$acc kernels
                do iq = 0, numflux - 1
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            if (sT_temp(c)>0) then
                                mQ_temp(c, iq, s) = mT_temp(c,  s) * alpha_ts(jt(c), iq, s) * Q_ts(jt(c), iq) &
                                                    * pQ_temp(c, iq) / sT_temp(c)

                                ! unless there is nothing in storage
                            else
                                mQ_temp(c, iq, s) = 0.
                            end if
                        end do
                    enddo
                enddo
                !$acc end kernels
                call cpu_time(finish)
                runtime(21) = runtime(21) + 1000*(finish-start)

                ! Reaction mass accounting
                ! If there are first-order reactions, get the total mass rate
                call cpu_time(start)
                !$acc kernels
                do c = 0, N - 1
                    mR_temp(c, :) = k1_ts(jt(c), :) * (C_eq_ts(jt(c), :) * sT_temp(c) - mT_temp(c, :))
                end do
                !$acc end kernels
                call cpu_time(finish)
                runtime(22) = runtime(22) + 1000*(finish-start)

                if (jacobian) then
                    call cpu_time(start)
                    !$acc kernels
                    fs_temp = 0.
                    fsQ_temp = 0.
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(23) = runtime(23) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do c = 0, N - 1
                        if (sT_temp(c)>0) then
                            do iq = 0, numflux - 1
                                fsQ_temp(c, :, iq) = fsQ_temp(c, :, iq) &
                                        + ds_temp(c, :) * pQ_temp(c, iq) * Q_ts(jt(c), iq) / sT_temp(c)
                            end do
                        end if
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(24) = runtime(24) + 1000*(finish-start)
                    call cpu_time(start)
                    do iq = 0, numflux - 1
                        do ic = component_index_list(iq), component_index_list(iq+1) - 1
                            !$acc kernels
                            !$acc loop independent
                            do c = 0, N - 1
                                ! sensitivity to point before the start
                                if ((leftbreakpt_top(c, ic)>=0).and.(leftbreakpt_top(c, ic)<numbreakpt_list(ic)-1)) then
                                    ip = breakpt_index_list(ic) + leftbreakpt_top( c, ic)
                                    call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                    dS = SAS_lookup(jt(c), ip+1) - SAS_lookup(jt(c), ip)
                                    dP = P_list(jt(c), ip+1) - P_list(jt(c), ip)
                                    call f_debug('dP/dS start    ', (/dP/dS/))
                                    fs_temp(c, ip) = fs_temp(c, ip) &
                                        + dP / (dS*dS) * sT_temp(c) * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                end if
                                ! sensitivity to point after the end
                                if ((leftbreakpt_bot(c, ic)+1>0).and.(leftbreakpt_bot(c, ic)+1<=numbreakpt_list(ic)-1)) then
                                    ip = breakpt_index_list(ic) + leftbreakpt_bot(c, ic) + 1
                                    call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                    dS = SAS_lookup(jt(c), ip) - SAS_lookup(jt(c), ip-1)
                                    dP = P_list(jt(c), ip) - P_list(jt(c), ip-1)
                                    call f_debug('dP/dS end      ', (/dP/dS/))
                                    fs_temp(c, ip) = fs_temp(c, ip) &
                                        - dP / (dS*dS) * sT_temp(c) * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                end if
                                ! sensitivity to point within
                                if (leftbreakpt_bot(c, ic)>leftbreakpt_top(c, ic)) then
                                    call f_debug('leftbreakpt_bot, _start', &
                                    (/leftbreakpt_bot(ic, c)*one8, leftbreakpt_top(ic, c)*one8/))
                                    do leftbreakpt=leftbreakpt_top(c, ic)+1, leftbreakpt_bot(c, ic)
                                        ip = breakpt_index_list(ic) + leftbreakpt
                                        call f_debug('iq, ic, ip, c ', (/iq*one8, ic*one8, ip*one8, c*one8/))
                                        if (leftbreakpt>0) then
                                            dSs = SAS_lookup(jt(c), ip) - SAS_lookup(jt(c), ip-1)
                                            dPs = P_list(jt(c), ip) - P_list(jt(c), ip-1)
                                        else
                                            dSs = 1.
                                            dPs = 0.
                                        end if
                                        if (leftbreakpt<numbreakpt_list(ic)-1) then
                                            dSe = SAS_lookup(jt(c), ip+1) - SAS_lookup(jt(c), ip)
                                            dPe = P_list(jt(c), ip+1) - P_list(jt(c), ip)
                                        else
                                            dSe = 1.
                                            dPe = 0.
                                        end if
                                        call f_debug('dP/dS middle   ', (/dPe/dSe , dPs/dSs/))
                                        fs_temp(c, ip) = fs_temp(c, ip) &
                                            - (dPe/dSe - dPs/dSs) / h * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                    end do
                                end if
                            end do
                            !$acc end kernels
                        end do
                    end do
                    call cpu_time(finish)
                    runtime(25) = runtime(25) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc kernels
                    fm_temp = 0
                    fmQ_temp = 0
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(26) = runtime(26) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do iq = 0, numflux - 1
                            !$acc loop independent
                            do ip = 0, numbreakpt_total - 1
                                !$acc loop independent
                                do c = 0, N - 1
                                    if (sT_temp(c)>0) then
                                        fmQ_temp(c, ip, iq, s) = fmQ_temp(c, ip, iq, s) &
                                            + dm_temp(c, ip, s) * alpha_ts(jt(c), iq,s) * Q_ts(jt(c), iq) &
                                            * pQ_temp(c, iq) / sT_temp(c)
                                    end if
                                end do
                            end do
                        end do
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do ip = 0, numbreakpt_total - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                fmR_temp(c, ip, s) = fmR_temp(c, ip, s) &
                                    + k1_ts(jt(c), s) * (C_eq_ts(jt(c), s) * ds_temp(c, ip) - dm_temp(c, ip, s))
                            end do
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(27) = runtime(27) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc loop independent
                    do iq = 0, numflux - 1
                        do ic = component_index_list(iq), component_index_list(iq+1) - 1
                            !$acc kernels
                            !$acc loop independent
                            do c = 0, N - 1
                            ! sensitivity to point before the start
                                if ((leftbreakpt_top(c, ic)>=0).and.(leftbreakpt_top(c, ic)<numbreakpt_list(ic)-1)) then
                                    ip = breakpt_index_list(ic) + leftbreakpt_top(c, ic)
                                    dS = SAS_lookup(jt(c), ip+1) - SAS_lookup(jt(c), ip)
                                    dP = P_list(jt(c), ip+1) - P_list(jt(c), ip)
                                    fm_temp(c, ip, :) = fm_temp(c, ip, :) &
                                        + dP / (dS*dS) * mT_temp(c, :)&
                                        * alpha_ts(jt(c), iq, :) * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                end if
                                ! sensitivity to point after the end
                                if ((leftbreakpt_bot(c, ic) + 1>0).and.(leftbreakpt_bot(c, ic) + 1<=numbreakpt_list(ic)-1)) then
                                    ip = breakpt_index_list(ic) + leftbreakpt_bot(c, ic) + 1
                                    dS = SAS_lookup(jt(c), ip) - SAS_lookup(jt(c), ip-1)
                                    dP = P_list(jt(c), ip) - P_list(jt(c), ip-1)
                                    fm_temp(c, ip, :) = fm_temp(c, ip, :) &
                                        - dP / (dS*dS) * mT_temp(c, :)&
                                        * alpha_ts(jt(c), iq, :) * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                end if
                                ! sensitivity to point within
                                if (leftbreakpt_bot(c, ic)>leftbreakpt_top(c, ic)) then
                                    do leftbreakpt = leftbreakpt_top(c, ic)+1, leftbreakpt_bot(c, ic)
                                        ip = breakpt_index_list(ic) + leftbreakpt
                                        if (leftbreakpt>0) then
                                            dSs = SAS_lookup(jt(c), ip) - SAS_lookup(jt(c), ip-1)
                                            dPs = P_list(jt(c), ip) - P_list(jt(c), ip-1)
                                        else
                                            dSs = 1.
                                            dPs = 0.
                                        end if
                                        if (leftbreakpt<numbreakpt_list(ic)-1) then
                                            dSe = SAS_lookup(jt(c), ip+1) - SAS_lookup(jt(c), ip)
                                            dPe = P_list(jt(c), ip+1) - P_list(jt(c), ip)
                                        else
                                            dSe = 1.
                                            dPe = 0.
                                        end if
                                        fm_temp(c, ip, :) = fm_temp(c, ip, :) &
                                            - (dPe/dSe - dPs/dSs) * mT_temp(c, :) / sT_temp(c) / h &
                                            * weights_ts(jt(c), ic) * Q_ts(jt(c), iq)
                                    end do
                                end if
                            end do
                            !$acc end kernels
                            call cpu_time(finish)
                            runtime(28) = runtime(28) + 1000*(finish-start)
                        end do
                    end do
                end if
                ! ########################## ^^ GET FLUX ^^ ##########################

                call cpu_time(start)
                !$acc kernels
                !$acc loop independent
                do iq = 0, numflux - 1
                    !$acc loop independent
                    do c = 0, N - 1
                        pQ_aver(c, iq)  = pQ_aver(c, iq)  + rk_coeff(rk) * pQ_temp(c, iq)
                    end do
                end do
                !$acc loop independent
                do s = 0, numsol - 1
                    !$acc loop independent
                    do iq = 0, numflux - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            mQ_aver(c, iq, s)  = mQ_aver(c, iq, s)  + rk_coeff(rk) * mQ_temp(c, iq, s)
                        end do
                    end do
                end do
                !$acc loop independent
                do s = 0, numsol - 1
                    !$acc loop independent
                    do c = 0, N - 1
                        mR_aver(c, s)  = mR_aver(c, s)  + rk_coeff(rk) * mR_temp(c, s)
                    end do
                end do
                !$acc end kernels
                call cpu_time(finish)
                runtime(29) = runtime(29) + 1000*(finish-start)
                if (jacobian) then
                    call cpu_time(start)
                    !$acc kernels
                    !$acc loop independent
                    do ip = 0, numbreakpt_total - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            fs_aver  (c, ip)= fs_aver(c, ip)+ rk_coeff(rk) * fs_temp(c, ip)
                        end do
                    end do
                    !$acc loop independent
                    do iq = 0, numflux - 1
                        !$acc loop independent
                        do ip = 0, numbreakpt_total - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                fsQ_aver (c, ip, iq)= fsQ_aver(c, ip, iq)+ rk_coeff(rk) * fsQ_temp(c, ip, iq)
                            end do
                        end do
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do ip = 0, numbreakpt_total - 1
                            !$acc loop independent
                            do c = 0, N - 1
                                fm_aver  (c, ip, s)= fm_aver(c, ip, s)+ rk_coeff(rk) * fm_temp(c, ip, s)
                                fmR_aver (c, ip, s)= fmR_aver(c, ip, s)+ rk_coeff(rk) * fmR_temp(c, ip, s)
                            end do
                        end do
                    end do
                    !$acc loop independent
                    do s = 0, numsol - 1
                        !$acc loop independent
                        do iq = 0, numflux - 1
                            !$acc loop independent
                            do ip = 0, numbreakpt_total - 1
                                !$acc loop independent
                                do c = 0, N - 1
                                    fmQ_aver (c, ip, iq, s)= fmQ_aver(c, ip, iq, s)+ rk_coeff(rk) * fmQ_temp(c, ip, iq, s)
                                end do
                            end do
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(30) = runtime(30) + 1000*(finish-start)
                end if

                !$acc update self(pQ_temp, sT_temp)
                call f_debug('pQ_temp end            ', pQ_temp(:, 0))

                end if
                if (rk==4) then
                    call f_debug('FINALIZE  rk           ', (/rk*one8, iT_s*one8/))
                    ! zero out the probabilities if there is no outflux this timestep
                    call cpu_time(start)
                    !$acc kernels
                    do iq = 0, numflux - 1
                        !$acc loop independent
                        do c = 0, N - 1
                            if (Q_ts( jt(c), iq)==0) then
                                pQ_aver(c, iq) = 0.
                                mQ_aver(c, iq, :) = 0.
                            end if
                        end do
                    end do
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(31) = runtime(31) + 1000*(finish-start)
                    call cpu_time(start)
                    !$acc kernels
                    pQ_temp  = pQ_aver
                    mQ_temp  = mQ_aver
                    mR_temp  = mR_aver
                    !$acc end kernels
                    call cpu_time(finish)
                    runtime(32) = runtime(32) + 1000*(finish-start)
                    if (jacobian) then
                        call cpu_time(start)
                        !$acc kernels
                        !$acc loop independent
                        do c = 0, N - 1
                            !$acc loop independent
                            do ip = 0, numbreakpt_total - 1
                                fs_temp(c, ip)  = fs_aver(c, ip)
                            end do
                            !$acc loop independent
                            do iq = 0, numflux - 1
                                !$acc loop independent
                                do ip = 0, numbreakpt_total - 1
                                    fsQ_temp(c, ip, iq) = fsQ_aver(c, ip, iq)
                                end do
                            end do
                            !$acc loop independent
                            do s = 0, numsol - 1
                                !$acc loop independent
                                do ip = 0, numbreakpt_total - 1
                                    fm_temp(c, ip, s)  = fm_aver(c, ip, s)
                                    fmR_temp(c, ip, s) = fmR_aver(c, ip, s)
                                end do
                                !$acc loop independent
                                do iq = 0, numflux - 1
                                    !$acc loop independent
                                    do ip = 0, numbreakpt_total - 1
                                        fmQ_temp(c, ip, iq, s) = fmQ_aver(c, ip, iq, s)
                                    end do
                                end do
                            end do
                        end do
                        !$acc end kernels
                        call cpu_time(finish)
                        runtime(33) = runtime(33) + 1000*(finish-start)
                    end if
                !$acc update self(pQ_aver, sT_temp)
                call f_debug('pQ_aver                ', pQ_aver( :, 0))
                end if
            end do
            call f_debug_blank()

            ! Update the state with the new estimates
            call cpu_time(start)
            !$acc kernels
            sT_start = sT_temp
            mT_start = mT_temp
            !$acc end kernels
            call cpu_time(finish)
            runtime(34) = runtime(34) + 1000*(finish-start)
            if (jacobian) then
                call cpu_time(start)
                !$acc kernels
                ds_start = ds_temp
                dm_start = dm_temp
                !$acc end kernels
                call cpu_time(finish)
                runtime(35) = runtime(35) + 1000*(finish-start)
            end if
            ! Aggregate data from substep to timestep
            call cpu_time(start)
            !$acc kernels
            STcum_top_start(0) = STcum_bot_start(0)
            STcum_bot_start(0) = STcum_bot_start(0) + sT_init_ts(iT) * h
            !$acc loop independent
            do c = 0, N - 1
                STcum_top_start(jt_s(c)+1) = STcum_bot_start(jt_s(c)+1)
                STcum_bot_start(jt_s(c)+1) = STcum_bot_start(jt_s(c)+1) + sT_start(c) * h
            end do
            ! Get the timestep-averaged transit time distribution
            !!$acc loop independent
            do jt_i = 0, timeseries_length - 1
                do jt_substep = 0, n_substeps - 1
                    c = mod(N + jt_i * n_substeps + jt_substep - iT_s, N)
                    !print *, c, jt_i, jt(c), jt_i * n_substeps + jt_substep, jt_s(c)
                    if (jt_substep<iT_substep) then
                        if (iT<max_age-1) then
                            pQ_ts(jt_i, :   , iT+1) = pQ_ts(jt_i, :   , iT+1) + pQ_aver(c, :   ) * norm
                            mQ_ts(jt_i, :, :, iT+1) = mQ_ts(jt_i, :, :, iT+1) + mQ_aver(c, :, :) * norm
                            mR_ts(jt_i, :   , iT+1) = mR_ts(jt_i, :   , iT+1) + mR_aver(c, :   ) * norm
                        end if
                    else
                        pQ_ts(jt_i, :   , iT) = pQ_ts(jt_i, :   , iT) + pQ_aver(c, :   ) * norm
                        mQ_ts(jt_i, :, :, iT) = mQ_ts(jt_i, :, :, iT) + mQ_aver(c, :, :) * norm
                        mR_ts(jt_i, :   , iT) = mR_ts(jt_i, :   , iT) + mR_aver(c, :   ) * norm
                    end if
                end do
            end do

            do jt_i = 0, timeseries_length - 1
                do jt_substep = 0, n_substeps - 1
                    c = mod(N + jt_i * n_substeps + jt_substep - iT_s, N)
                    if (jacobian) then
                        do iq = 0, numflux - 1
                            if (Q_ts(jt_i, iq)>0) then
                                dW_ts(jt_i, :, iq) = dW_ts(jt_i, :, iq) + fsQ_aver(c, :, iq) / Q_ts(jt_i, iq) * norm * dt
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                        dW_ts(jt_i, ip, iq) = dW_ts(jt_i, ip, iq) + fs_aver(c, ip) / Q_ts(jt_i, iq) * norm * dt
                                    enddo
                                enddo
                                dC_ts(jt_i, :, iq, :) = dC_ts(jt_i, :, iq, :) &
                                + fmQ_aver(c, :, iq, :) / Q_ts(jt_i, iq) * norm * dt
                                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                    do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                        dC_ts(jt_i, ip, iq, :) = dC_ts(jt_i, ip, iq, :) &
                                        + fm_aver(c, ip, :) / Q_ts(jt_i, iq) * norm * dt
                                    enddo
                                enddo
                            end if
                        enddo
                    end if

                    ! Extract substep state at timesteps
                    ! age-ranked storage at the end of the timestep
                    if (jt_substep==n_substeps-1) then
                        sT_ts(jt_i+1, iT) = sT_ts(jt_i+1, iT) + sT_start(c) / n_substeps
                        ! parameter sensitivity
                        if (jacobian) then
                            do ip = 0, numbreakpt_total - 1
                                ds_ts(jt_i+1, ip, iT) = ds_ts(jt_i+1, ip, iT) + ds_start(c, ip) / n_substeps
                            enddo
                        end if
                        ! Age-ranked solute mass
                        do s = 0, numsol - 1
                            mT_ts(jt_i+1, s, iT) = mT_ts(jt_i+1, s, iT) + mT_start(c,  s) / n_substeps
                            ! parameter sensitivity
                            if (jacobian) then
                                do ip = 0, numbreakpt_total - 1
                                    dm_ts(jt_i+1, ip, s, iT) = dm_ts(jt_i+1, ip, s, iT) + dm_start(c, ip, s) / n_substeps
                                enddo
                            end if
                        enddo
                    end if
                enddo
            enddo
            !$acc end kernels
            call cpu_time(finish)
            runtime(36) = runtime(36) + 1000*(finish-start)

            call f_debug('sT_ts(iT, :)     ', sT_ts(iT, :))
            call f_debug('pQ_aver0         ', pQ_aver(0,:))
            call f_debug('pQ_aver1         ', pQ_aver(1,:))
            call f_debug('pQ_ts(iT, :, 0)'  , pQ_ts(iT, :, 0))
            call f_debug('pQ_ts(iT, :, 0)'  , pQ_ts(iT, :, 1))

            iT_prev = iT

        enddo
        ! Print some updates
        if (mod(iT, 10).eq.0) then
            write (tempdebugstring, *) '...Done ', (iT), &
                    'of', (max_age)
            call f_verbose(tempdebugstring)
        endif
    enddo
    call cpu_time(start)
    !$acc update self(sT_ts)
    !$acc update self(mT_ts)
    !$acc update self(ds_ts)
    !$acc update self(dm_ts)
    !$acc update self(dC_ts)
    !$acc update self(dW_ts)
    !$acc update self(pQ_ts)
    !$acc update self(mQ_ts)
    !$acc update self(mR_ts)
    call cpu_time(finish)
    runtime(37) = runtime(37) + 1000*(finish-start)
    !do i = 1, 37
        !print '("Time for section ",i2," = :",f10.0,": milliseconds")', i, runtime(i)
    !end do

    ! Calculate a water balance
    ! Difference of starting and ending age-ranked storage
    !$acc kernels
    do iT = 0, max_age - 1
        do jt_i = 0, timeseries_length - 1
            if (iT==0) then
                WaterBalance_ts( jt_i, iT) = J_ts(jt_i) - sT_ts( jt_i+1, iT)
            else
                WaterBalance_ts( jt_i, iT) = sT_ts( jt_i, iT-1) - sT_ts( jt_i+1, iT)
            end if
            ! subtract time-averaged water fluxes
            do iq = 0, numflux - 1
                WaterBalance_ts( jt_i, iT) = WaterBalance_ts( jt_i, iT) - (Q_ts( jt_i, iq) * pQ_ts( jt_i, iq, iT)) * dt
            end do

            ! Calculate a solute balance
            ! Difference of starting and ending age-ranked mass
            if (iT==0) then
                do s = 0, numsol - 1
                    SoluteBalance_ts( jt_i, s, iT) = C_J_ts(jt_i, s) * J_ts(jt_i) - mT_ts( jt_i+1, s, iT)
                end do
            else
                SoluteBalance_ts( jt_i, :, iT) = mT_ts( jt_i, :, iT-1) - mT_ts( jt_i+1, :, iT)
            end if
            ! Subtract timestep-averaged mass fluxes
            do iq = 0, numflux - 1
                SoluteBalance_ts( jt_i, :, iT) = SoluteBalance_ts( jt_i, :, iT) - (mQ_ts( jt_i, iq, :, iT)) * dt
            end do
            ! Reacted mass
            SoluteBalance_ts( jt_i, :, iT) = SoluteBalance_ts( jt_i, :, iT) + mR_ts( jt_i, :, iT) * dt
        enddo


    enddo ! End of main loop
    !$acc end kernels
    !$acc end data region

    call f_verbose('...Finalizing...')

    ! get the old water fraction
    P_old = 1 - sum(pQ_ts, DIM=3) * dt

    do s = 0, numsol - 1
        do iq = 0, numflux - 1

            where (Q_ts(:, iq)>0)
                ! From the age-ranked mass
                C_Q_ts(:, iq, s) = sum(mQ_ts(:, iq, s, :), DIM=2) / Q_ts(:, iq) * dt

                ! From the old water concentration
                C_Q_ts(:, iq, s) = C_Q_ts(:, iq, s) + alpha_ts( :, iq, s) * C_old(s) * P_old(:, iq)

            end where

            if (jacobian) then
                do ip = 0, numbreakpt_total - 1
                    where (Q_ts(:, iq)>0)
                        dC_ts(:, ip, iq, s) = dC_ts(:, ip, iq, s) - C_old(s) * dW_ts(:, ip, iq)
                    end where
                end do
            end if


        enddo
    enddo


    call f_verbose('...Finished...')

contains


    subroutine f_debug_blank()
        ! Prints a blank line
        if (debug) then
            print *, ''
        endif
    end subroutine f_debug_blank

    subroutine f_debug(debugstring, debugdblepr)
        ! Prints debugging information
        implicit none
        character(len = *), intent(in) :: debugstring
        real(8), dimension(:), intent(in) :: debugdblepr
        if (debug) then
            print 1, debugstring, debugdblepr
            1 format (A26, *(f16.10))
        endif
    end subroutine f_debug


    subroutine f_warning(debugstring)
        ! Prints informative information
        implicit none
        !$acc routine seq
        character(len = *), intent(in) :: debugstring
        if (warning) then
            print *, debugstring
        endif
    end subroutine f_warning


    subroutine f_verbose(debugstring)
        ! Prints informative information
        implicit none
        !$acc routine seq
        character(len = *), intent(in) :: debugstring
        if (verbose) then
            print *, debugstring
        endif
    end subroutine f_verbose

end subroutine solveSAS
