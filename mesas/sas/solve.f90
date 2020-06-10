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
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning
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
    integer :: kT, iT, jt, ik, jk, i_prev, kc
    real(8) :: h
    real(8), dimension(0:numbreakpt_total-1, 0:numflux - 1, 0:timeseries_length-1) :: dW_ts
    real(8), dimension(0:numflux - 1, 0:timeseries_length - 1) :: P_old
    integer, dimension(0:numcomponent_total) :: breakpt_index_list
    integer, dimension(0:numflux) :: component_index_list
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_bot
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_top
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps) :: PQcum_prev
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: PQcum_bot
    real(8), dimension(0:numflux - 1, 0:timeseries_length * n_substeps - 1) :: PQcum_top
    integer, dimension(0:numcomponent_total-1, 0:timeseries_length * n_substeps) :: leftbreakpt_prev
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
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_end
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mT_start
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mT_temp
    real(8), dimension(0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: mT_end
    real(8), dimension(0:numbreakpt_total-1, 0:timeseries_length * n_substeps - 1) :: ds_start
    real(8), dimension(0:numbreakpt_total-1, 0:timeseries_length * n_substeps - 1) :: ds_temp
    real(8), dimension(0:numbreakpt_total-1, 0:timeseries_length * n_substeps - 1) :: ds_end
    real(8), dimension(0:numbreakpt_total-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: dm_start
    real(8), dimension(0:numbreakpt_total-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: dm_temp
    real(8), dimension(0:numbreakpt_total-1, 0:numsol - 1, 0:timeseries_length * n_substeps - 1) :: dm_end
    real(8) :: one8, norm
    real(8) :: dS, dP, dSe, dPe, dSs, dPs
    character(len = 128) :: tempdebugstring
    integer :: iq, s, M, N, ip, ic
    integer :: carry

    call f_verbose('...Initializing arrays...')
    one8 = 1.0

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
    STcum_prev = 0.
    STcum_bot = 0.
    STcum_top = 0.
    PQcum_prev = 0.
    PQcum_bot = 0.
    PQcum_top = 0.
    leftbreakpt_prev = 0
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
    sT_end = 0.
    mT_start = 0.
    mT_temp = 0.
    mT_end = 0.
    ds_start = 0.
    ds_temp = 0.
    ds_end = 0.
    dm_start = 0.
    dm_temp = 0.
    dm_end = 0.
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

    call f_verbose('...Setting initial conditions...')
    sT_ts(:, 0) = sT_init_ts
    do s = 0, numsol - 1
        mT_ts(:, 0, s) = mT_init_ts(:, s)
    end do

    call f_verbose('...Starting main loop...')
    do iT = 0, max_age - 1

        ! Start the substep loop
        do kT = 0, n_substeps - 1

            ik = iT * n_substeps + kT

            call f_debug_blank()
            call f_debug_blank()
            call f_debug('Agestep, Substep', (/ iT * one8, kT * one8/))
            call f_debug_blank()

            ! Initialize the value at t=0
            if (ik>0) then
                sT_start(0) = sT_init_ts(i_prev)
                mT_start(:, 0) = mT_init_ts(i_prev, :)
                ds_start(:, 0) = 0.
                dm_start(:, :, 0) = 0.
            end if

            ! Copy the state variables from the end of the previous substep as the start of this one
            ! but shifted by one substep
            sT_start(1:N - 1) = sT_end(0:N - 2)
            mT_start(:, 1:N - 1) = mT_end(:, 0:N - 2)
            ds_start(:, 1:N - 1) = ds_end(:, 0:N - 2)
            dm_start(:, :, 1:N - 1) = dm_end(:, :, 0:N - 2)

            ! These will hold the evolving state variables
            ! They are global variables modified by the new_state function

            do jt = 0, timeseries_length - 1
                do kc = 0, n_substeps - 1
                    jk = jt * n_substeps + kc

                    ! This is the Runge-Kutta 4th order algorithm

                    call f_debug('RK', (/1._8/))
                    sT_temp(jk) = sT_start(jk)
                    mT_temp(:, jk) = mT_start(:, jk)
                    ds_temp(:, jk) = ds_start(:, jk)
                    dm_temp(:, :, jk) = dm_start(:, :, jk)
                    call update_aver(0)

                    call get_flux(0.0D0)
                    call update_aver(1)
                    call new_state(h / 2)

                    call get_flux(h / 2)
                    call update_aver(2)
                    call new_state(h / 2)

                    call get_flux(h / 2)
                    call update_aver(2)
                    call new_state(h)

                    call get_flux(h)
                    call update_aver(1)

                    ! zero out the probabilities if there is no outflux this timestep
                    do iq = 0, numflux - 1
                        if (Q_ts(iq, jt)==0) then
                            pQ_aver(iq, jk) = 0.
                            mQ_aver(iq, :, jk) = 0.
                        end if
                    end do

                    ! Update the state with the new estimates
                    call update_aver(-1)
                    call new_state(h)
                    sT_end(jk) = sT_temp(jk)
                    mT_end(:, jk) = mT_temp(:, jk)
                    ds_end(:, jk) = ds_temp(:, jk)
                    dm_end(:, :, jk) = dm_temp(:, :, jk)

                    ! Aggregate flux data from substep to timestep

                    ! Get the timestep-averaged transit time distribution
                    norm = 1.0 / n_substeps / n_substeps
                    if ((iT<max_age-1).and.(kc<kT)) then
                        carry = 1
                    else
                        carry = 0
                    end if
                    pQ_ts(iT+carry, jt, :) = pQ_ts(iT+carry, jt, :) + pQ_aver(:, jk) * norm
                    mQ_ts(iT+carry, jt, :, :) = mQ_ts(iT+carry, jt, :, :) + mQ_aver(:, :, jk) * norm
                    mR_ts(iT+carry, jt, :) = mR_ts(iT+carry, jt, :) + mR_aver(:, jk) * norm
                    do iq = 0, numflux - 1
                        if (Q_ts(iq, jt)>0) then
                            dW_ts(:, iq, jt) = dW_ts(:, iq, jt) + fsQ_aver(:, iq, jk) / Q_ts(iq, jt) * norm * dt
                            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                    dW_ts(ip, iq, jt) = dW_ts(ip, iq, jt) + fs_aver(ip, jk) / Q_ts(iq, jt) * norm * dt
                                enddo
                            enddo
                            dC_ts(:, iq, :, jt) = dC_ts(:, iq, :, jt) &
                                    + fmQ_aver(:, iq, :, jk) / Q_ts(iq, jt) * norm * dt
                            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                    dC_ts(ip, iq, :, jt) = dC_ts(ip, iq, :, jt) &
                                            + fm_aver(ip, :, jk) / Q_ts(iq, jt) * norm * dt
                                enddo
                            enddo
                        end if
                    enddo

                    ! Extract substep state at timesteps
                    ! age-ranked storage at the end of the timestep
                    if (kc==n_substeps-1) then
                        sT_ts(iT, jt+1) = sT_ts(iT, jt+1) + sT_end(jk) / n_substeps
                        ! parameter sensitivity
                        do ip = 0, numbreakpt_total - 1
                            ds_ts(iT, jt+1, ip) = ds_ts(iT, jt+1, ip) + ds_end(ip, jk) / n_substeps
                        enddo
                        ! Age-ranked solute mass
                        do s = 0, numsol - 1
                            mT_ts(iT, jt+1, s) = mT_ts(iT, jt+1, s) + mT_end( s, jk) / n_substeps
                            ! parameter sensitivity
                            do ip = 0, numbreakpt_total - 1
                                dm_ts(iT, jt+1, ip, s) = dm_ts(iT, jt+1, ip, s) + dm_end(ip, s, jk) / n_substeps
                            enddo
                        enddo
                    end if
                enddo
            enddo

            ! Update the cumulative instantaneous trackers
            if (ik>0) then
                STcum_prev(0) = STcum_prev(0) + sT_init_ts(i_prev) * h
                do jk = 0, N - 1
                    STcum_prev(jk+1) = STcum_prev(jk+1) + sT_end(jk) * h
                end do
                do jk = 0, N
                    call get_SAS(STcum_prev(jk), PQcum_prev(:, jk), leftbreakpt_prev(:, jk), jk)
                end do
            end if
            i_prev = iT

        enddo

        ! Calculate a water balance
        ! Difference of starting and ending age-ranked storage
        do jt = 0, timeseries_length - 1
            if (iT==0) then
                WaterBalance_ts(iT, jt) = J_ts(jt) - sT_ts(iT, jt+1)
            else
                WaterBalance_ts(iT, jt) = sT_ts(iT-1, jt) - sT_ts(iT, jt+1)
            end if
            ! subtract time-averaged water fluxes
            WaterBalance_ts(iT, jt) = WaterBalance_ts(iT, jt) - sum(Q_ts(:, jt) * pQ_ts(iT, jt, :)) * dt

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
            SoluteBalance_ts(iT, jt, :) = SoluteBalance_ts(iT, jt, :) - sum(mQ_ts(iT, jt, :, :), DIM=1) * dt
            ! Reacted mass
            SoluteBalance_ts(iT, jt, :) = SoluteBalance_ts(iT, jt, :) + mR_ts(iT, jt, :) * dt
        enddo

        ! Print some updates
        if (mod(iT, 10).eq.0) then
            write (tempdebugstring, *) '...Done ', (iT), &
                    'of', (timeseries_length)
            call f_verbose(tempdebugstring)
        endif

    enddo ! End of main loop

    call f_verbose('...Finalizing...')

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

            do ip = 0, numbreakpt_total - 1
                where (Q_ts(iq, :)>0)
                    dC_ts(ip, iq, s, :) = dC_ts(ip, iq, s, :) - C_old(s) * dW_ts(ip, iq, :)
                end where
            end do


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


    subroutine get_SAS(STcum_in, PQcum_out, leftbreakpt_out, jks)
        ! Call the sas function and get the transit time distribution
        integer, intent(in) :: jks
        real(8), intent(in) :: STcum_in
        real(8), intent(out), dimension(0:numflux - 1) :: PQcum_out
        integer, intent(out), dimension(0:numcomponent_total - 1) :: leftbreakpt_out
        real(8) :: PQcum_component
        integer :: jts
        if (jks==N) then
            jts = timeseries_length - 1
        else
            jts = jks/n_substeps
        end if
        ! Main lookup loop
        PQcum_out(:) = 0.
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                PQcum_component = 0
                call lookup(&
                        SAS_lookup(breakpt_index_list(ic):breakpt_index_list(ic + 1) - 1, jts), &
                        P_list(breakpt_index_list(ic):breakpt_index_list(ic + 1) - 1, jts), &
                        STcum_in, &
                        PQcum_component, &
                        leftbreakpt_out(ic), &
                        leftbreakpt_prev(ic, jks), &
                        numbreakpt_list(ic), 1)
                PQcum_out(iq) = PQcum_out(iq) + weights_ts(ic, jts) * PQcum_component
            enddo
        enddo
        end subroutine get_SAS

    subroutine lookup(xa, ya, x, y, ia, i0, na, n)
        ! A simple lookup table
        implicit none
        integer, intent(in) :: na, n
        real(8), intent(in), dimension(0:na - 1) :: xa
        real(8), intent(in), dimension(0:na - 1) :: ya
        real(8), intent(in) :: x
        real(8), intent(inout) :: y
        integer, intent(inout) :: ia
        integer, intent(inout) :: i0
        integer :: i
        real(8) :: dif, grad
        logical :: foundit
        if (x.le.xa(0)) then
            y = ya(0)
            ia = -1
        else if (x.ge.xa(na - 1)) then
            y = ya(na - 1)
            ia = na - 1
        else
            foundit = .FALSE.
            do i = 0, na - 1
                if (x.lt.xa(i)) then
                    ia = i - 1
                    foundit = .TRUE.
                    exit
                endif
            enddo
            if (.not. foundit) then
                call f_warning('I could not find the ST value. This should never happen!!!')
                y = ya(na - 1)
                ia = na - 1
            else
                dif = x - xa(ia)
                grad = (ya(ia + 1) - ya(ia)) / (xa(ia + 1) - xa(ia))
                y = ya(ia) + dif * grad
            endif
        endif
    end subroutine

    subroutine update_aver(coeff)
        ! Calculates the fluxes in the given the curent state
        implicit none
        integer, intent(in) :: coeff
        ! Average RK4 estimated change in the state variables
        if (coeff>0) then
            pQ_aver(:, jk)  = pQ_aver(:, jk)  + coeff * pQ_temp(:, jk) / 6.
            mQ_aver(:, :, jk)  = mQ_aver(:, :, jk)  + coeff * mQ_temp(:, :, jk) / 6.
            mR_aver(:, jk)  = mR_aver(:, jk)  + coeff * mR_temp(:, jk) / 6.
            fs_aver(:, jk)  = fs_aver(:, jk)  + coeff * fs_temp(:, jk) / 6.
            fsQ_aver(:, :, jk) = fsQ_aver(:, :, jk) + coeff * fsQ_temp(:, :, jk) / 6.
            fm_aver(:, :, jk)  = fm_aver(:, :, jk)  + coeff * fm_temp(:, :, jk) / 6.
            fmQ_aver(:, :, :, jk) = fmQ_aver(:, :, :, jk) + coeff * fmQ_temp(:, :, :, jk) / 6.
            fmR_aver(:, :, jk) = fmR_aver(:, :, jk) + coeff * fmR_temp(:, :, jk) / 6.
        else if (coeff==-1) then
            pQ_temp(:, jk)  = pQ_aver(:, jk)
            mQ_temp(:, :, jk)  = mQ_aver(:, :, jk)
            mR_temp(:, jk)  = mR_aver(:, jk)
            fs_temp(:, jk)  = fs_aver(:, jk)
            fsQ_temp(:, :, jk) = fsQ_aver(:, :, jk)
            fm_temp(:, :, jk)  = fm_aver(:, :, jk)
            fmQ_temp(:, :, :, jk) = fmQ_aver(:, :, :, jk)
            fmR_temp(:, :, jk) = fmR_aver(:, :, jk)
        else
            pQ_aver(:, jk)  = 0.
            mQ_aver(:, :, jk)  = 0.
            mR_aver(:, jk)  = 0.
            fs_aver(:, jk)  = 0.
            fsQ_aver(:, :, jk) = 0.
            fm_aver(:, :, jk)  = 0.
            fmQ_aver(:, :, :, jk) = 0.
            fmR_aver(:, :, jk) = 0.
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
                STcum_top(jk) = jk
                PQcum_top(:, jk) = 0
                leftbreakpt_top(:, jk) = -1
                STcum_bot(jk) = 0
                PQcum_bot(:, jk) = 0
                leftbreakpt_bot(:, jk) = -1
                pQ_temp(:, jk) = 0
            else
                STcum_top(jk) = 0
                PQcum_top(:, jk) = 0
                leftbreakpt_top(:, jk) = -1
                STcum_bot(jk) = 0 + sT_temp(jk) * hr
                call get_SAS(STcum_bot(jk), PQcum_bot(:, jk), leftbreakpt_bot(:, jk), jk)
                pQ_temp(:, jk) = (PQcum_bot(:, jk) - PQcum_top(:, jk)) / hr
            end if
        else
            STcum_top(jk) = STcum_prev(jk) * (1-hr/h) + (STcum_prev(jk+1) + sT_end(jk) * h) * (hr/h)
            call get_SAS(STcum_top(jk), PQcum_top(:, jk), leftbreakpt_top(:, jk), jk)
            STcum_bot(jk) = STcum_top(jk) + sT_temp(jk) * h
            call get_SAS(STcum_bot(jk), PQcum_bot(:, jk), leftbreakpt_bot(:, jk), jk)
            pQ_temp(:, jk) = (PQcum_bot(:, jk) - PQcum_top(:, jk)) / h
        end if

        do iq = 0, numflux - 1
            if (sT_temp(jk)==0) then
                pQ_temp(iq, jk) = 0
            end if
        end do

        ! Solute mass flux accounting

        do iq = 0, numflux - 1
            do s = 0, numsol - 1

                ! Get the mass flux out
                if (sT_temp(jk)>0) then
                    mQ_temp(iq, s, jk) = mT_temp( s, jk) * alpha_ts(iq, s, jt) * Q_ts(iq, jt) &
                                                * pQ_temp(iq, jk) / sT_temp(jk)

                    ! unless there is nothing in storage
                else
                    mQ_temp(iq, s, jk) = 0.

                end if
            enddo
        enddo

        ! Reaction mass accounting
        ! If there are first-order reactions, get the total mass rate
        mR_temp(:, jk) = k1_ts(:, jt) * (C_eq_ts(:, jt) * sT_temp(jk) - mT_temp(:, jk))

        fs_temp(:, jk) = 0.
        fsQ_temp(:, :, jk) = 0.
        if (sT_temp( jk)>0) then
            do iq = 0, numflux - 1
                fsQ_temp(:, iq, jk) = fsQ_temp(:, iq, jk) + ds_temp(:, jk) * pQ_temp(iq, jk) * Q_ts(iq, jt) / sT_temp(jk)
            end do
        end if
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                ! sensitivity to point before the start
                if ((leftbreakpt_top(ic, jk)>=0).and.(leftbreakpt_top(ic, jk)<numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_top(ic, jk)
                    !call f_debug('iq, ic, ip, jk ', (/iq*one8, ic*one8, ip*one8, jk*one8/))
                    dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                    dP = P_list(ip+1, jt) - P_list(ip, jt)
                    !call f_debug('dP/dS start    ', (/dP/dS/))
                    fs_temp(ip, jk) = fs_temp(ip, jk) &
                            + dP / (dS*dS) * sT_temp(jk) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point after the end
                if ((leftbreakpt_bot(ic, jk)+1>0).and.(leftbreakpt_bot(ic, jk)+1<=numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, jk) + 1
                    !call f_debug('iq, ic, ip, jk ', (/iq*one8, ic*one8, ip*one8, jk*one8/))
                    dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                    dP = P_list(ip, jt) - P_list(ip-1, jt)
                    !call f_debug('dP/dS end      ', (/dP/dS/))
                    fs_temp(ip, jk) = fs_temp(ip, jk) &
                            - dP / (dS*dS) * sT_temp(jk) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point within
                if (leftbreakpt_bot(ic, jk)>leftbreakpt_top(ic, jk)) then
                    call f_debug('leftbreakpt_bot, _start', &
                            (/leftbreakpt_bot(ic, jk)*one8, leftbreakpt_top(ic, jk)*one8/))
                    do leftbreakpt=leftbreakpt_top(ic, jk)+1, leftbreakpt_bot(ic, jk)
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
                        fs_temp(ip, jk) = fs_temp(ip, jk) &
                                - (dPe/dSe - dPs/dSs) / h * weights_ts(ic, jt) * Q_ts(iq, jt)
                    end do
                end if
            end do
        end do
        fm_temp(:, :, jk) = 0
        fmQ_temp(:, :, :, jk) = 0
        do iq = 0, numflux - 1
            do ip = 0, numbreakpt_total - 1
                if (sT_temp(jk)>0) then
                    fmQ_temp(ip, iq, :, jk) = fmQ_temp(ip, iq, :, jk) &
                            + dm_temp(ip, :, jk) * alpha_ts(iq,:, jt) * Q_ts(iq, jt) * pQ_temp(iq, jk) / sT_temp(jk)
                end if
                fmR_temp(ip, :, jk) = fmR_temp(ip, :, jk) &
                        + k1_ts(:, jt) * (C_eq_ts(:, jt) * ds_temp(ip, jk) - dm_temp(ip, :, jk))
            end do
        end do
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                ! sensitivity to point before the start
                if ((leftbreakpt_top(ic, jk)>=0).and.(leftbreakpt_top(ic, jk)<numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_top(ic, jk)
                    dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                    dP = P_list(ip+1, jt) - P_list(ip, jt)
                    fm_temp(ip, :, jk) = fm_temp(ip, :, jk) &
                            + dP / (dS*dS) * mT_temp(:, jk)&
                                    * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point after the end
                if ((leftbreakpt_bot(ic, jk) + 1>0).and.(leftbreakpt_bot(ic, jk) + 1<=numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, jk) + 1
                    dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                    dP = P_list(ip, jt) - P_list(ip-1, jt)
                    fm_temp(ip, :, jk) = fm_temp(ip, :, jk) &
                            - dP / (dS*dS) * mT_temp(:, jk)&
                                    * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point within
                if (leftbreakpt_bot(ic, jk)>leftbreakpt_top(ic, jk)) then
                    do leftbreakpt = leftbreakpt_top(ic, jk)+1, leftbreakpt_bot(ic, jk)
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
                        fm_temp(ip, :, jk) = fm_temp(ip, :, jk) &
                                - (dPe/dSe - dPs/dSs) * mT_temp(:, jk) / sT_temp(jk) / h &
                                        * weights_ts(ic, jt) * Q_ts(iq, jt)
                    end do
                end if
            end do
        end do

    end subroutine get_flux


    subroutine new_state(hr)
        ! Calculates the state given the fluxes

        real(8), intent(in) :: hr

        !call f_debug('new_state', (/0._8/))

        ! Calculate the new age-ranked storage
        sT_temp(jk) = sT_start(jk) ! Initial value
        mT_temp(:, jk) = mT_start(:, jk) + mR_temp(:, jk) * hr ! Initial value + reaction
        ! Fluxes in & out
        if (ik == 0) then
            sT_temp(jk) = sT_temp(jk) + J_ts(jt) * hr / h
            mT_temp(:, jk) = mT_temp(:, jk) + J_ts(jt) * C_J_ts(:, jt) * (hr/h)
        end if
        sT_temp(jk) = sT_temp(jk) - sum(Q_ts(:, jt) * pQ_temp(:, jk)) * hr
        mT_temp(:, jk) = mT_temp(:, jk) - sum(mQ_temp(:, :, jk), dim=1) * hr
        if (sT_temp(jk)<0) then
            call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
        end if

        ! Calculate new parameter sensitivity
        ds_temp(:, jk) = ds_start(:, jk) - fs_temp(:, jk) * hr - sum(fsQ_temp(:, :, jk), dim=2) * hr
        dm_temp(:, :, jk) = dm_start(:, :, jk) - fm_temp(:, :, jk) * hr - sum(fmQ_temp(:, :, :, jk), dim=3) * hr &
                + fmR_temp(:, :, jk) * hr

    end subroutine new_state

end subroutine solve
