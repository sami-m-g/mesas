! -*- f90 -*-
subroutine solve(J_ts, Q_ts, SAS_lookup, P_list, weights_ts, sT_init_ts, dt, &
        verbose, debug, warning, &
        mT_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
        n_substeps, numcomponent_list, numbreakpt_list, numflux, numsol, max_age, &
        timeseries_length, numcomponent_total, numbreakpt_total, &
        sT_ts, pQ_ts, WaterBalance_ts, &
        mT_ts, mQ_ts, mR_ts, C_Q_ts, ds_ts, dm_ts, dC_ts, SoluteBalance_ts)
    implicit none

    ! Start by declaring and initializing all the variables we will be using
    integer, intent(in) :: n_substeps, numflux, numsol, max_age, &
            timeseries_length, numcomponent_total, numbreakpt_total
    real(8), intent(in), dimension(0:timeseries_length - 1) :: J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1) :: Q_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numcomponent_total - 1) :: weights_ts
    real(8), intent(in), dimension(0:numbreakpt_total - 1, 0:timeseries_length - 1) :: SAS_lookup
    real(8), intent(in), dimension(0:numbreakpt_total - 1, 0:timeseries_length - 1) :: P_list
    real(8), intent(in), dimension(0:max_age - 1) :: sT_init_ts
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning
    real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mT_init_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: k1_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_eq_ts
    real(8), intent(in), dimension(0:numsol - 1) :: C_old
    integer, intent(in), dimension(0:numflux - 1) :: numcomponent_list
    integer, intent(in), dimension(0:numcomponent_total - 1) :: numbreakpt_list
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: C_Q_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length) :: sT_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numsol - 1) :: mT_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numbreakpt_total-1) :: ds_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_ts
    real(8), intent(out), dimension(0:timeseries_length-1, 0:numbreakpt_total-1, 0:numflux - 1, 0:numsol - 1) :: dC_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1) :: pQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: mR_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1) :: WaterBalance_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: SoluteBalance_ts
    integer :: k, iT, jt, ik, jk, i_prev, js
    real(8) :: h
    real(8), dimension(0:timeseries_length-1, 0:numbreakpt_total-1, 0:numflux - 1) :: dW_ts
    real(8), dimension(0:timeseries_length - 1, 0:numflux - 1) :: P_old
    integer, dimension(0:numcomponent_total) :: breakpt_index_list
    integer, dimension(0:numflux) :: component_index_list
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: J_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: C_J_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: Q_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: C_eq_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: k1_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numcomponent_total - 1) :: weights_ss
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_top
    real(8), dimension(0:timeseries_length * n_substeps, 0:numflux - 1) :: PQcum_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_top
    integer, dimension(0:timeseries_length * n_substeps, 0:numcomponent_total-1) :: leftbreakpt_prev
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
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mT_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mT_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1) :: ds_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1) :: ds_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1) :: ds_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numbreakpt_total-1, 0:numsol - 1) :: dm_end
    real(8) :: one8, norm
    real(8) :: dS, dP, dSe, dPe, dSs, dPs
    character(len = 128) :: tempdebugstring
    integer :: iq, s, M, N, ip, ic
    logical :: carryover

    call f_verbose('...Initializing arrays...')
    one8 = 1.0

    C_Q_ts(:, :, :) = 0.
    sT_ts(:, :) = 0.
    mT_ts(:, :, :) = 0.
    ds_ts(:, :, :) = 0.
    dm_ts(:, :, :, :) = 0.
    dC_ts(:, :, :, :) = 0.
    dW_ts(:, :, :) = 0.
    pQ_ts(:, :, :) = 0.
    mQ_ts(:, :, :, :) = 0.
    mR_ts(:, :, :) = 0.
    WaterBalance_ts(:, :) = 0.
    SoluteBalance_ts(:, :, :) = 0.
    P_old(:, :) = 0.
    breakpt_index_list(:) = 0
    component_index_list(:) = 0
    STcum_prev(:) = 0.
    STcum_bot(:) = 0.
    STcum_top(:) = 0.
    PQcum_prev(:, :) = 0.
    PQcum_bot(:, :) = 0.
    PQcum_top(:, :) = 0.
    leftbreakpt_prev(:, :) = 0
    leftbreakpt_bot(:, :) = 0
    leftbreakpt_top(:, :) = 0
    pQ_temp(:, :) = 0.
    pQ_aver(:, :) = 0.
    mQ_temp(:, :, :) = 0.
    mQ_aver(:, :, :) = 0.
    mR_temp(:, :) = 0.
    mR_aver(:, :) = 0.
    fs_temp(:, :) = 0.
    fs_aver(:, :) = 0.
    fsQ_temp(:, :, :) = 0.
    fsQ_aver(:, :, :) = 0.
    fm_temp(:, :, :) = 0.
    fm_aver(:, :, :) = 0.
    fmQ_temp(:, :, :, :) = 0.
    fmQ_aver(:, :, :, :) = 0.
    fmR_temp(:, :, :) = 0.
    fmR_aver(:, :, :) = 0.
    sT_start(:) = 0.
    sT_temp(:) = 0.
    sT_end(:) = 0.
    mT_start(:, :) = 0.
    mT_temp(:, :) = 0.
    mT_end(:, :) = 0.
    ds_start(:, :) = 0.
    ds_temp(:, :) = 0.
    ds_end(:, :) = 0.
    dm_start(:, :, :) = 0.
    dm_temp(:, :, :) = 0.
    dm_end(:, :, :) = 0.
    i_prev = -1

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

    J_ss = reshape(spread(J_ts, 1, n_substeps), (/N/))
    C_J_ss = reshape(spread(C_J_ts, 1, n_substeps), (/N, numsol/))
    Q_ss = reshape(spread(Q_ts, 1, n_substeps), (/N, numflux/))
    alpha_ss = reshape(spread(alpha_ts, 1, n_substeps), (/N, numflux, numsol/))
    C_eq_ss = reshape(spread(C_eq_ts, 1, n_substeps), (/N, numsol/))
    k1_ss = reshape(spread(k1_ts, 1, n_substeps), (/N, numsol/))
    weights_ss = reshape(spread(weights_ts, 1, n_substeps), (/N, numcomponent_total/))
    !call f_debug('J_ss           ', J_ss)
    do s = 0, numsol - 1
        !call f_debug('C_J_ss           ', C_J_ss(:, s))
    end do
    do iq = 0, numflux - 1
        !call f_debug('Q_ss           ', Q_ss(:,iq))
    end do

    call f_verbose('...Setting initial conditions...')
    sT_ts(:, 0) = sT_init_ts
    do s = 0, numsol - 1
        mT_ts(:, 0, s) = mT_init_ts(:, s)
    end do

    call f_verbose('...Starting main loop...')
    do iT = 0, max_age - 1

        ! Start the substep loop
        do k = 0, n_substeps - 1

            ik = iT * n_substeps + k

            call f_debug_blank()
            call f_debug_blank()
            call f_debug('Agestep, Substep', (/ iT * one8, k * one8/))
            call f_debug_blank()

            ! Copy the state variables from the end of the previous substep as the start of this one
            ! but shifted by one substep
            sT_start(1:N - 1) = sT_end(0:N - 2)
            mT_start(1:N - 1, :) = mT_end(0:N - 2, :)
            ds_start(1:N - 1, :) = ds_end(0:N - 2, :)
            dm_start(1:N - 1, :, :) = dm_end(0:N - 2, :, :)
            ! Initialize the value at t=0
            if (ik>0) then
                sT_start(0) = sT_init_ts(i_prev)
                mT_start(0, :) = mT_init_ts(i_prev, :)
                ds_start(0, :) = 0.
                dm_start(0, :, :) = 0.
            end if

            ! These will hold the evolving state variables
            ! They are global variables modified by the new_state function

            ! This is the Runge-Kutta 4th order algorithm

            call f_debug('RK', (/1._8/))
            sT_temp = sT_start
            mT_temp = mT_start
            ds_temp = ds_start
            dm_temp = dm_start
            call update_aver(0)

            call get_flux(0.0D0)
            call update_aver(1)
            !call f_debug_blank()
            call f_debug('RK', (/2._8/))
            call new_state(h / 2)

            call get_flux(h / 2)
            call update_aver(2)
            !call f_debug_blank()
            call f_debug('RK', (/3._8/))
            call new_state(h / 2)

            call get_flux(h / 2)
            call update_aver(2)
            !call f_debug_blank()
            call f_debug('RK', (/4._8/))
            call new_state(h)

            call get_flux(h)
            call update_aver(1)

            ! zero out the probabilities if there is no outflux this timestep
            where (Q_ss==0)
                pQ_aver = 0.
            end where
            do s = 0, numsol - 1
                where (Q_ss==0)
                    mQ_aver(:, :, s) = 0.
                end where
            end do

            call f_debug('RK final', (/4._8/))
            call f_debug('pQ_aver        ', pQ_aver(:,0))
            ! Update the state with the new estimates
            call update_aver(-1)
            call new_state(h)
            sT_end = sT_temp
            mT_end = mT_temp
            ds_end = ds_temp
            dm_end = dm_temp

            ! output some debugging info if desired
            do iq = 0, numflux - 1
                !call f_debug('pQ_aver        ', pQ_aver(:, iq))
            enddo
            do s = 0, numsol - 1
                !call f_debug('mT_end         ', mT_end(:, s))
                !call f_debug('mR_aver        ', mR_aver(:, s))
                do iq = 0, numflux - 1
                    !call f_debug('mQ_aver        ', mQ_aver(:, iq, s))
                enddo
            enddo

            ! Aggregate flux data from substep to timestep

            ! Get the timestep-averaged transit time distribution
            norm = 1.0 / n_substeps / n_substeps
            carryover = ((n_substeps>1) .and. (k>0) .and. (iT<max_age-1))
            do iq = 0, numflux - 1
                do jt = 0, timeseries_length - 1
                    js = jt * n_substeps
                    pQ_ts(iT, jt, iq) = pQ_ts(iT, jt, iq) + sum(pQ_aver((js+k):(js+n_substeps-1), iq)) * norm
                    if (carryover) then
                        pQ_ts(iT+1, jt, iq) = pQ_ts(iT+1, jt, iq) + sum(pQ_aver((js):(js+k-1), iq)) * norm
                    endif
                enddo
                do ip = 0, numbreakpt_total - 1
                    do jt = 0, timeseries_length - 1
                        if (Q_ts(jt, iq)>0) then
                            js = jt * n_substeps
                            dW_ts(jt, ip, iq) = dW_ts(jt, ip, iq) &
                                    + sum(fsQ_aver((js+k):(js+n_substeps-1), ip, iq))/Q_ts(jt, iq) * norm * dt
                            if (carryover) then
                                dW_ts(jt, ip, iq) = dW_ts(jt, ip, iq) &
                                        + sum(fsQ_aver((js):(js+k-1), ip, iq))/Q_ts(jt, iq) * norm * dt
                            endif
                        endif
                    enddo
                enddo
                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                    do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                        do jt = 0, timeseries_length - 1
                            if (Q_ts(jt, iq)>0) then
                                js = jt * n_substeps
                                dW_ts(jt, ip, iq) = dW_ts(jt, ip, iq) &
                                        + sum(fs_aver((js+k):(js+n_substeps-1), ip))/Q_ts(jt, iq) * norm * dt
                                if (carryover) then
                                    dW_ts(jt, ip, iq) = dW_ts(jt, ip, iq) &
                                            + sum(fs_aver((js):(js+k-1), ip))/Q_ts(jt, iq) * norm * dt
                                endif
                            endif
                        enddo
                    enddo
                enddo
                do s = 0, numsol - 1
                    do jt = 0, timeseries_length - 1
                        js = jt * n_substeps
                        mQ_ts(iT, jt, iq, s) = mQ_ts(iT, jt, iq, s) + sum(mQ_aver((js+k):(js+n_substeps-1), iq, s)) * norm
                        if (carryover) then
                            mQ_ts(iT+1, jt, iq, s) = mQ_ts(iT+1, jt, iq, s) + sum(mQ_aver((js):(js+k-1), iq, s)) * norm
                        endif
                    enddo
                    do ip = 0, numbreakpt_total - 1
                        do jt = 0, timeseries_length - 1
                            if (Q_ts(jt, iq)>0) then
                                js = jt * n_substeps
                                dC_ts(jt, ip, iq, s) = dC_ts(jt, ip, iq, s) &
                                        + sum(fmQ_aver((js+k):(js+n_substeps-1), ip, iq, s))/Q_ts(jt, iq) * norm * dt
                                if (carryover) then
                                    dC_ts(jt, ip, iq, s) = dC_ts(jt, ip, iq, s) &
                                            + sum(fmQ_aver((js):(js+k-1), ip, iq, s))/Q_ts(jt, iq) * norm * dt
                                endif
                            endif
                        enddo
                    enddo
                    do ic = component_index_list(iq), component_index_list(iq+1) - 1
                        do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                            do jt = 0, timeseries_length - 1
                                if (Q_ts(jt, iq)>0) then
                                    js = jt * n_substeps
                                    dC_ts(jt, ip, iq, s) = dC_ts(jt, ip, iq, s) &
                                            + sum(fm_aver((js+k):(js+n_substeps-1), ip, s))/Q_ts(jt, iq) * norm * dt
                                    if (carryover) then
                                        dC_ts(jt, ip, iq, s) = dC_ts(jt, ip, iq, s) &
                                                + sum(fm_aver((js):(js+k-1), ip, s))/Q_ts(jt, iq) * norm * dt
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            do s = 0, numsol - 1
                do jt = 0, timeseries_length - 1
                    js = jt * n_substeps
                    mR_ts(iT, jt, s) = mR_ts(iT, jt, s) + sum(mR_aver((js+k):(js+n_substeps-1), s)) * norm
                    if (carryover) then
                        mR_ts(iT+1, jt, s) = mR_ts(iT+1, jt, s) + sum(mR_aver((js):(js+k-1), s)) * norm
                    endif
                enddo
            enddo

            ! Extract substep state at timesteps
            ! age-ranked storage at the end of the timestep
            sT_ts(iT, 1:) = sT_ts(iT, 1:) + sT_end(n_substeps-1:N-1:n_substeps) / n_substeps
            ! parameter sensitivity
            do ip = 0, numbreakpt_total - 1
                ds_ts(iT, 1:, ip) = ds_ts(iT, 1:, ip) + ds_end(n_substeps-1:N-1:n_substeps, ip) / n_substeps
            enddo
            ! Age-ranked solute mass
            do s = 0, numsol - 1
                mT_ts(iT, 1:, s) = mT_ts(iT, 1:, s) + mT_end(n_substeps-1:N-1:n_substeps, s) / n_substeps
                ! parameter sensitivity
                do ip = 0, numbreakpt_total - 1
                    dm_ts(iT, 1:, ip, s) = dm_ts(iT, 1:, ip, s) + dm_end(n_substeps-1:N-1:n_substeps, ip, s) / n_substeps
                enddo
            enddo

            ! Update the cumulative instantaneous trackers
            if (ik>0) then
                STcum_prev(1:N) = STcum_prev(1:N) + sT_end * h
                STcum_prev(0) = STcum_prev(0) + sT_init_ts(i_prev) * h
                call get_SAS(STcum_prev, PQcum_prev, leftbreakpt_prev, N+1)
            end if
            i_prev = iT

        enddo

        ! Calculate a water balance
        ! Difference of starting and ending age-ranked storage
        if (iT==0) then
            WaterBalance_ts(iT, :) = J_ts - sT_ts(iT, 1:)
        else
            WaterBalance_ts(iT, :) = sT_ts(iT-1, 0:timeseries_length-1) - sT_ts(iT, 1:timeseries_length)
        end if
        ! subtract time-averaged water fluxes
        WaterBalance_ts(iT, :) = WaterBalance_ts(iT, :) - sum(Q_ts * pQ_ts(iT, :, :), DIM=2) * dt

        ! Calculate a solute balance
        ! Difference of starting and ending age-ranked mass
        if (iT==0) then
            do s = 0, numsol - 1
                SoluteBalance_ts(iT, :, s) = C_J_ts(:, s) * J_ts - mT_ts(iT, 1:, s)
            end do
        else
            SoluteBalance_ts(iT, :, :) = mT_ts(iT-1, 0:timeseries_length-1, :) &
                    - mT_ts(iT, 1:timeseries_length, :)
        end if
        ! Subtract timestep-averaged mass fluxes
        SoluteBalance_ts(iT, :, :) = SoluteBalance_ts(iT, :, :) - sum(mQ_ts(iT, :, :, :), DIM=2) * dt
        ! Reacted mass
        SoluteBalance_ts(iT, :, :) = SoluteBalance_ts(iT, :, :) + mR_ts(iT, :, :) * dt

        ! Print some updates
        if (mod(iT, 10).eq.0) then
            write (tempdebugstring, *) '...Done ', (iT), &
                    'of', (timeseries_length)
            call f_verbose(tempdebugstring)
        endif

    enddo ! End of main loop

    call f_verbose('...Finalizing...')

    ! get the old water fraction
    P_old = 1 - sum(pQ_ts, DIM=1) * dt
    !call f_debug('P_old', (/P_old/))

    ! Estimate the outflow concentration
    do iq = 0, numflux - 1
        do s = 0, numsol - 1
            do iT = 0, max_age - 1
                !call f_debug('mQ_ts          ', mQ_ts(iT, :, iq, s))
            enddo
        enddo
    enddo
    do s = 0, numsol - 1
        do iq = 0, numflux - 1

            where (Q_ts(:,iq)>0)
                ! From the age-ranked mass
                C_Q_ts(:, iq, s) = sum(mQ_ts(:, :, iq, s), DIM=1) / Q_ts(:,iQ) * dt

                ! From the old water concentration
                C_Q_ts(:, iq, s) = C_Q_ts(:, iq, s) + alpha_ts(:, iq, s) * C_old(s) * P_old(:,iq)


            end where

            do ip = 0, numbreakpt_total - 1
                where (Q_ts(:,iq)>0)
                    dC_ts(:, ip, iq, s) = dC_ts(:, ip, iq, s) - C_old(s) * dW_ts(:, ip, iq)
                end where
            end do

            !call f_debug('C_Q_ts         ', C_Q_ts(:, iq, s))

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
            1 format (A16, *(f16.10))
        endif
    end subroutine f_debug


    subroutine f_warning(debugstring)
        ! Prints informative information
        implicit none
        character(len = *), intent(in) :: debugstring
        if (warning) then
            print *, debugstring
        endif
    end subroutine f_warning


    subroutine f_verbose(debugstring)
        ! Prints informative information
        implicit none
        character(len = *), intent(in) :: debugstring
        if (verbose) then
            print *, debugstring
        endif
    end subroutine f_verbose


    subroutine get_SAS(STcum_in, PQcum_out, leftbreakpt_out, n_array)
        ! Call the sas function and get the transit time distribution
        integer, intent(in) :: n_array
        real(8), intent(in), dimension(0:n_array - 1) :: STcum_in
        real(8), intent(out), dimension(0:n_array - 1, 0:numflux - 1) :: PQcum_out
        integer, intent(out), dimension(0:n_array - 1, 0:numcomponent_total - 1) :: leftbreakpt_out
        real(8), dimension(0:n_array - 1) :: PQcum_component
        ! Main lookup loop
        PQcum_out(:, :) = 0.
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                !call f_debug('breakpt_index_list        ', (/breakpt_index_list(ic)*one8, breakpt_index_list(ic + 1)-1*one8/))
                PQcum_component(:) = 0
                do jt = 0, timeseries_length - 1
                    jk = jt * n_substeps
                    !call f_debug('getting entries', (/jk*one8, (jk+n_substeps-1)*one8, ic*one8, jt*one8/))
                    !call f_debug('SAS_lookup     ', SAS_lookup(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt))
                    !call f_debug('P_list         ', P_list(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt))
                    call lookup(&
                            SAS_lookup(breakpt_index_list(ic):breakpt_index_list(ic + 1) - 1, jt), &
                            P_list(breakpt_index_list(ic):breakpt_index_list(ic + 1) - 1, jt), &
                            STcum_in(jk:jk+n_substeps-1), &
                            PQcum_component(jk:jk+n_substeps-1), &
                            leftbreakpt_out(jk:jk+n_substeps-1, ic), &
                            leftbreakpt_prev(jk:jk+n_substeps-1, ic), &
                            numbreakpt_list(ic), n_substeps)
                end do
                PQcum_out(:N-1, iq) = PQcum_out(:N-1, iq) + weights_ss(:N-1, ic) * PQcum_component(:N-1)
                if (n_array>N) then
                    jt = timeseries_length - 1
                    !call f_debug('getting last entry', (/n_array*one8, N*one8, ic*one8, jt*one8/))
                    !call f_debug('SAS_lookup     ', SAS_lookup(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt))
                    !call f_debug('P_list         ', P_list(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt))
                    call lookup(&
                            SAS_lookup(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt), &
                            P_list(breakpt_index_list(ic):breakpt_index_list(ic + 1)-1, jt),  &
                            STcum_in(N:N), &
                            PQcum_component(N:N), &
                            leftbreakpt_out(N:N, ic), &
                            leftbreakpt_prev(N:N, ic), &
                            numbreakpt_list(ic), 1)
                    PQcum_out(N, iq) = PQcum_out(N, iq) + weights_ss(N-1, ic) * PQcum_component(N)
                end if
            enddo
            !call f_debug('STcum_in       ', STcum_in(:))
            !call f_debug('PQcum_out      ', PQcum_out(:, iq))
        enddo
        end subroutine get_SAS

    subroutine lookup(xa, ya, x, y, ia, i0, na, n)
        ! A simple lookup table
        implicit none
        integer, intent(in) :: na, n
        real(8), intent(in), dimension(0:na - 1) :: xa
        real(8), intent(in), dimension(0:na - 1) :: ya
        real(8), intent(in), dimension(0:n - 1) :: x
        real(8), intent(inout), dimension(0:n - 1) :: y
        integer, intent(inout), dimension(0:n - 1) :: ia
        integer, intent(inout), dimension(0:n - 1) :: i0
        integer :: i, j
        real(8) :: dif, grad
        logical :: foundit
        do j = 0, n - 1
            if (x(j).le.xa(0)) then
                y(j) = ya(0)
                ia(j) = -1
            else if (x(j).ge.xa(na - 1)) then
                y(j) = ya(na - 1)
                ia(j) = na - 1
            else
                foundit = .FALSE.
                do i = 0, na - 1
                    if (x(j).lt.xa(i)) then
                        ia(j) = i - 1
                        foundit = .TRUE.
                        exit
                    endif
                enddo
                if (.not. foundit) then
                    call f_warning('I could not find the ST value. This should never happen!!!')
                    y(j) = ya(na - 1)
                    ia(j) = na - 1
                else
                    i = ia(j)
                    dif = x(j) - xa(i)
                    grad = (ya(i + 1) - ya(i)) / (xa(i + 1) - xa(i))
                    y(j) = ya(i) + dif * grad
                endif
            endif
        enddo
    end subroutine

    subroutine update_aver(coeff)
        ! Calculates the fluxes in the given the curent state
        implicit none
        integer, intent(in) :: coeff
        ! Average RK4 estimated change in the state variables
        if (coeff>0) then
            pQ_aver  = pQ_aver  + coeff * pQ_temp / 6.
            mQ_aver  = mQ_aver  + coeff * mQ_temp / 6.
            mR_aver  = mR_aver  + coeff * mR_temp / 6.
            fs_aver  = fs_aver  + coeff * fs_temp / 6.
            fsQ_aver = fsQ_aver + coeff * fsQ_temp / 6.
            fm_aver  = fm_aver  + coeff * fm_temp / 6.
            fmQ_aver = fmQ_aver + coeff * fmQ_temp / 6.
            fmR_aver = fmR_aver + coeff * fmR_temp / 6.
        else if (coeff==-1) then
            pQ_temp  = pQ_aver
            mQ_temp  = mQ_aver
            mR_temp  = mR_aver
            fs_temp  = fs_aver
            fsQ_temp = fsQ_aver
            fm_temp  = fm_aver
            fmQ_temp = fmQ_aver
            fmR_temp = fmR_aver
        else
            pQ_aver  = 0.
            mQ_aver  = 0.
            mR_aver  = 0.
            fs_aver  = 0.
            fsQ_aver = 0.
            fm_aver  = 0.
            fmQ_aver = 0.
            fmR_aver = 0.
        end if
    end subroutine

    subroutine get_flux(hr)
        ! Calculates the fluxes in the given the curent state
        implicit none
        real(8), intent(in) :: hr
        integer iq, s, ip, kk, leftbreakpt
        !call f_debug('get_flux', (/hr/))

        ! Use the SAS function lookup table to convert age-rank storage to the fraction of discharge of age T at each t

        ! First get the cumulative age-ranked storage
        if (ik==0) then
            if (hr==0) then
                STcum_top = 0
                PQcum_top = 0
                leftbreakpt_top = -1
                STcum_bot = 0
                PQcum_bot = 0
                leftbreakpt_bot = -1
                pQ_temp = 0
            else
                STcum_top = 0
                PQcum_top = 0
                leftbreakpt_top = -1
                STcum_bot = 0 + sT_temp * hr
                call get_SAS(STcum_bot, PQcum_bot, leftbreakpt_bot, N)
                pQ_temp = (PQcum_bot - PQcum_top) / hr
            end if
        else
            STcum_top = STcum_prev(0:N-1) * (1-hr/h) + (STcum_prev(1:N) + sT_end * h) * (hr/h)
            call get_SAS(STcum_top, PQcum_top, leftbreakpt_top, N)
            STcum_bot = STcum_top + sT_temp * h
            call get_SAS(STcum_bot, PQcum_bot, leftbreakpt_bot, N)
            pQ_temp = (PQcum_bot - PQcum_top) / h
        end if

        do iq = 0, numflux - 1
            where (sT_temp(:)==0)
                pQ_temp(:, iq) = 0
            end where
        end do

        ! Solute mass flux accounting
        !call f_debug('getting mQ_temp', (/0._8/))

        do iq = 0, numflux - 1
            do s = 0, numsol - 1

                ! Get the mass flux out
                where (sT_temp(:)>0)
                    mQ_temp(:, iq, s) = mT_temp(:, s) * alpha_ss(:, iq, s) * Q_ss(:, iq) * pQ_temp(:, iq) / sT_temp(:)

                    ! unless there is nothing in storage
                elsewhere
                    mQ_temp(:, iq, s) = 0.

                end where
            enddo
        enddo

        do s = 0, numsol - 1
            ! Reaction mass accounting

            ! If there are first-order reactions, get the total mass rate
            mR_temp(:, s) = k1_ss(:, s) * (C_eq_ss(:, s) * sT_temp(:) - mT_temp(:, s))

            !call f_debug('mQ_temp         ', mQ_temp(:, 0, s))
            !call f_debug('mE_temp         ', mQ_temp(:, 1, s))
            !call f_debug('mR_temp         ', mR_temp(:, s))

        enddo

        fs_temp(:, :) = 0.
        fsQ_temp(:, :, :) = 0.
        do iq = 0, numflux - 1
            do ip = 0, numbreakpt_total - 1
                where (sT_temp(:)>0)
                    fsQ_temp(:, ip, iq) = fsQ_temp(:, ip, iq) + ds_temp(:, ip) * pQ_temp(:, iq) * Q_ss(:, iq) / sT_temp
                end where
            end do
        end do
        do iq = 0, numflux - 1
            call f_debug_blank()
            !call f_debug_blank()
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                do jt = 0, timeseries_length - 1
                    do kk = 0, n_substeps - 1
                        jk = jt * n_substeps + kk
                        ! sensitivity to point before the start
                        if ((leftbreakpt_top(jk, ic)>=0).and.(leftbreakpt_top(jk, ic)<numbreakpt_list(ic)-1)) then
                            ip = breakpt_index_list(ic) + leftbreakpt_top(jk, ic)
                            !call f_debug('iq, ic, ip, jk ', (/iq*one8, ic*one8, ip*one8, jk*one8/))
                            dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                            dP = P_list(ip+1, jt) - P_list(ip, jt)
                            !call f_debug('dP/dS start    ', (/dP/dS/))
                            fs_temp(jk, ip) = fs_temp(jk, ip) &
                                    + dP / (dS*dS) * sT_temp(jk) * weights_ss(jk, ic) * Q_ss(jk, iq)
                        end if
                        ! sensitivity to point after the end
                        if ((leftbreakpt_bot(jk, ic)+1>0).and.(leftbreakpt_bot(jk, ic)+1<=numbreakpt_list(ic)-1)) then
                            ip = breakpt_index_list(ic) + leftbreakpt_bot(jk, ic) + 1
                            !call f_debug('iq, ic, ip, jk ', (/iq*one8, ic*one8, ip*one8, jk*one8/))
                            dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                            dP = P_list(ip, jt) - P_list(ip-1, jt)
                            !call f_debug('dP/dS end      ', (/dP/dS/))
                            fs_temp(jk, ip) = fs_temp(jk, ip) &
                                    - dP / (dS*dS) * sT_temp(jk) * weights_ss(jk, ic) * Q_ss(jk, iq)
                        end if
                        ! sensitivity to point within
                        if (leftbreakpt_bot(jk, ic)>leftbreakpt_top(jk, ic)) then
                            call f_debug('leftbreakpt_bot, _start', &
                                    (/leftbreakpt_bot(jk, ic)*one8, leftbreakpt_top(jk, ic)*one8/))
                            do leftbreakpt=leftbreakpt_top(jk, ic)+1, leftbreakpt_bot(jk, ic)
                                ip = breakpt_index_list(ic) + leftbreakpt
                                call f_debug('iq, ic, ip, jk ', (/iq*one8, ic*one8, ip*one8, jk*one8/))
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
                                call f_debug('dP/dS middle   ', (/dPe/dSe , dPs/dSs/))
                                fs_temp(jk, ip) = fs_temp(jk, ip) &
                                        - (dPe/dSe - dPs/dSs) / h * weights_ss(jk, ic) * Q_ss(jk, iq)
                            end do
                        end if
                    end do
                end do
            end do
        end do
        fm_temp(:, :, :) = 0
        fmQ_temp(:, :, :, :) = 0
        do s = 0, numsol - 1
            do iq = 0, numflux - 1
                do ip = 0, numbreakpt_total - 1
                    where (sT_temp(:)>0)
                        fmQ_temp(:, ip, iq, s) = fmQ_temp(:, ip, iq, s) &
                                + dm_temp(:, ip, s) * alpha_ss(:,iq,s) * Q_ss(:, iq) * pQ_temp(:, iq) / sT_temp
                    end where
                    fmR_temp(:, ip, s) = fmR_temp(:, ip, s) &
                            + k1_ss(:, s) * (C_eq_ss(:, s) * ds_temp(:, ip) - dm_temp(:, ip, s))
                end do
            end do
            do iq = 0, numflux - 1
                do ic = component_index_list(iq), component_index_list(iq+1) - 1
                    do jt = 0, timeseries_length - 1
                        do kk = 0, n_substeps - 1
                            jk = jt * n_substeps + kk
                            ! sensitivity to point before the start
                            if ((leftbreakpt_top(jk, ic)>=0).and.(leftbreakpt_top(jk, ic)<numbreakpt_list(ic)-1)) then
                                ip = breakpt_index_list(ic) + leftbreakpt_top(jk, ic)
                                dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                                dP = P_list(ip+1, jt) - P_list(ip, jt)
                                fm_temp(jk, ip, s) = fm_temp(jk, ip, s) &
                                        + dP / (dS*dS) * mT_temp(jk, s)&
                                                * alpha_ss(jk, iq, s) * weights_ss(jk, ic) * Q_ss(jk, iq)
                            end if
                            ! sensitivity to point after the end
                            if ((leftbreakpt_bot(jk, ic) + 1>0).and.(leftbreakpt_bot(jk, ic) + 1<=numbreakpt_list(ic)-1)) then
                                ip = breakpt_index_list(ic) + leftbreakpt_bot(jk, ic) + 1
                                dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                                dP = P_list(ip, jt) - P_list(ip-1, jt)
                                fm_temp(jk, ip, s) = fm_temp(jk, ip, s) &
                                        - dP / (dS*dS) * mT_temp(jk, s)&
                                                * alpha_ss(jk, iq, s) * weights_ss(jk, ic) * Q_ss(jk, iq)
                            end if
                            ! sensitivity to point within
                            if (leftbreakpt_bot(jk, ic)>leftbreakpt_top(jk, ic)) then
                                do leftbreakpt = leftbreakpt_top(jk, ic)+1, leftbreakpt_bot(jk, ic)
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
                                    fm_temp(jk, ip, s) = fm_temp(jk, ip, s) &
                                            - (dPe/dSe - dPs/dSs) * mT_temp(jk, s) / sT_temp(jk) / h &
                                                    * weights_ss(jk, ic) * Q_ss(jk, iq)
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do

        !call f_debug('get_flux finished', (/0._8/))
        !call f_debug_blank()

    end subroutine get_flux


    subroutine new_state(hr)
        ! Calculates the state given the fluxes

        real(8), intent(in) :: hr

        !call f_debug('new_state', (/0._8/))

        ! Calculate the new age-ranked storage
        sT_temp = sT_start ! Initial value
        ! Fluxes in & out
        if (ik == 0) then
            sT_temp = sT_temp + J_ss * hr / h
        end if
        do iq = 0, numflux - 1
            sT_temp = sT_temp - Q_ss(:, iq) * hr * pQ_temp(:, iq)
        enddo
        !call f_debug('sT_start ns    ', sT_start(:))
        !call f_debug('sT_temp ns      ', sT_temp(:))
        if (ANY(sT_temp<0)) then
                call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
        end if

        ! Print some debugging info
        do s = 0, numsol - 1
            !call f_debug('mT_start NS    ', mT_start(:, s))
            !call f_debug('mR_temp NS       ', mR_temp(:, s))
            do iq = 0, numflux - 1
                !call f_debug('mQ_temp NS       ', mQ_temp(:, iq, s))
            enddo
        enddo

        ! Calculate the new age-ranked mass
        do s = 0, numsol - 1
            ! Initial value + reaction
            mT_temp(:, s) = mT_start(:, s) + mR_temp(:, s) * hr
            ! Flux in
            if (ik==0) then
                mT_temp(:, s) = mT_temp(:, s) + J_ss(:) * C_J_ss(:, s) * (hr/h)
            end if
            ! Fluxes out
            do iq = 0, numflux - 1
                mT_temp(:, s) = mT_temp(:, s) - mQ_temp(:, iq, s) * hr
            enddo
        enddo

        do s = 0, numsol - 1
            !call f_debug('mT_temp NS      ', mT_temp(:, s))
            !call f_debug('C_J_ss NS      ', (/C_J_ss(:, s)/))
        enddo
        !call f_debug('J_ss NS        ', (/J_ss(:)/))
        !call f_debug('new_state finished', (/0._8/))
        !call f_debug_blank()

        ! Calculate new parameter sensitivity
        ds_temp = ds_start - fs_temp * hr - sum(fsQ_temp, dim=3) * hr
        dm_temp = dm_start - fm_temp * hr - sum(fmQ_temp, dim=3) * hr + fmR_temp * hr

    end subroutine new_state

end subroutine solve
