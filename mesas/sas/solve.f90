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
    integer :: iT_substep, iT, jt, iT_k, jt_k, i_prev, jt_substep
    real(8) :: h
    real(8), dimension(0:numbreakpt_total-1, 0:numflux - 1, 0:timeseries_length-1) :: dW_ts
    real(8), dimension(0:numflux - 1, 0:timeseries_length - 1) :: P_old
    integer, dimension(0:numcomponent_total) :: breakpt_index_list
    integer, dimension(0:numflux) :: component_index_list
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_prev
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
        do iT_substep = 0, n_substeps - 1

            iT_k = iT * n_substeps + iT_substep

            call f_debug_blank()
            call f_debug_blank()
            call f_debug('Agestep, Substep', (/ iT * one8, iT_substep * one8/))
            call f_debug_blank()

            ! Initialize the value at t=0
            if (iT_k>0) then
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
                do jt_substep = 0, n_substeps - 1
                    jt_k = jt * n_substeps + jt_substep

                    ! This is the Runge-Kutta 4th order algorithm

                    call f_debug('RK', (/1._8/))
                    sT_temp(jt_k) = sT_start(jt_k)
                    mT_temp(:, jt_k) = mT_start(:, jt_k)
                    ds_temp(:, jt_k) = ds_start(:, jt_k)
                    dm_temp(:, :, jt_k) = dm_start(:, :, jt_k)
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
                            pQ_aver(iq, jt_k) = 0.
                            mQ_aver(iq, :, jt_k) = 0.
                        end if
                    end do

                    ! Update the state with the new estimates
                    call update_aver(-1)
                    call new_state(h)
                    sT_end(jt_k) = sT_temp(jt_k)
                    mT_end(:, jt_k) = mT_temp(:, jt_k)
                    ds_end(:, jt_k) = ds_temp(:, jt_k)
                    dm_end(:, :, jt_k) = dm_temp(:, :, jt_k)

                    ! Aggregate flux data from substep to timestep

                    ! Get the timestep-averaged transit time distribution
                    norm = 1.0 / n_substeps / n_substeps
                    if ((iT<max_age-1).and.(jt_substep<iT_substep)) then
                        carry = 1
                    else
                        carry = 0
                    end if
                    pQ_ts(iT+carry, jt, :) = pQ_ts(iT+carry, jt, :) + pQ_aver(:, jt_k) * norm
                    mQ_ts(iT+carry, jt, :, :) = mQ_ts(iT+carry, jt, :, :) + mQ_aver(:, :, jt_k) * norm
                    mR_ts(iT+carry, jt, :) = mR_ts(iT+carry, jt, :) + mR_aver(:, jt_k) * norm
                    do iq = 0, numflux - 1
                        if (Q_ts(iq, jt)>0) then
                            dW_ts(:, iq, jt) = dW_ts(:, iq, jt) + fsQ_aver(:, iq, jt_k) / Q_ts(iq, jt) * norm * dt
                            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                    dW_ts(ip, iq, jt) = dW_ts(ip, iq, jt) + fs_aver(ip, jt_k) / Q_ts(iq, jt) * norm * dt
                                enddo
                            enddo
                            dC_ts(:, iq, :, jt) = dC_ts(:, iq, :, jt) &
                                    + fmQ_aver(:, iq, :, jt_k) / Q_ts(iq, jt) * norm * dt
                            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                                do ip = breakpt_index_list(ic), breakpt_index_list(ic+1) - 1
                                    dC_ts(ip, iq, :, jt) = dC_ts(ip, iq, :, jt) &
                                            + fm_aver(ip, :, jt_k) / Q_ts(iq, jt) * norm * dt
                                enddo
                            enddo
                        end if
                    enddo

                    ! Extract substep state at timesteps
                    ! age-ranked storage at the end of the timestep
                    if (jt_substep==n_substeps-1) then
                        sT_ts(iT, jt+1) = sT_ts(iT, jt+1) + sT_end(jt_k) / n_substeps
                        ! parameter sensitivity
                        do ip = 0, numbreakpt_total - 1
                            ds_ts(iT, jt+1, ip) = ds_ts(iT, jt+1, ip) + ds_end(ip, jt_k) / n_substeps
                        enddo
                        ! Age-ranked solute mass
                        do s = 0, numsol - 1
                            mT_ts(iT, jt+1, s) = mT_ts(iT, jt+1, s) + mT_end( s, jt_k) / n_substeps
                            ! parameter sensitivity
                            do ip = 0, numbreakpt_total - 1
                                dm_ts(iT, jt+1, ip, s) = dm_ts(iT, jt+1, ip, s) + dm_end(ip, s, jt_k) / n_substeps
                            enddo
                        enddo
                    end if
                enddo
            enddo

            ! Update the cumulative instantaneous trackers
            if (iT_k>0) then
                STcum_prev(0) = STcum_prev(0) + sT_init_ts(i_prev) * h
                do jt_k = 0, N - 1
                    STcum_prev(jt_k+1) = STcum_prev(jt_k+1) + sT_end(jt_k) * h
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
                        numbreakpt_list(ic))
                PQcum_out(iq) = PQcum_out(iq) + weights_ts(ic, jts) * PQcum_component
            enddo
        enddo
    end subroutine get_SAS

    subroutine lookup(xa, ya, x, y, ia, na)
        ! A simple lookup table
        implicit none
        integer, intent(in) :: na
        real(8), intent(in), dimension(0:na - 1) :: xa
        real(8), intent(in), dimension(0:na - 1) :: ya
        real(8), intent(in) :: x
        real(8), intent(inout) :: y
        integer, intent(inout) :: ia
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
            pQ_aver(:, jt_k)  = pQ_aver(:, jt_k)  + coeff * pQ_temp(:, jt_k) / 6.
            mQ_aver(:, :, jt_k)  = mQ_aver(:, :, jt_k)  + coeff * mQ_temp(:, :, jt_k) / 6.
            mR_aver(:, jt_k)  = mR_aver(:, jt_k)  + coeff * mR_temp(:, jt_k) / 6.
            fs_aver(:, jt_k)  = fs_aver(:, jt_k)  + coeff * fs_temp(:, jt_k) / 6.
            fsQ_aver(:, :, jt_k) = fsQ_aver(:, :, jt_k) + coeff * fsQ_temp(:, :, jt_k) / 6.
            fm_aver(:, :, jt_k)  = fm_aver(:, :, jt_k)  + coeff * fm_temp(:, :, jt_k) / 6.
            fmQ_aver(:, :, :, jt_k) = fmQ_aver(:, :, :, jt_k) + coeff * fmQ_temp(:, :, :, jt_k) / 6.
            fmR_aver(:, :, jt_k) = fmR_aver(:, :, jt_k) + coeff * fmR_temp(:, :, jt_k) / 6.
        else if (coeff==-1) then
            pQ_temp(:, jt_k)  = pQ_aver(:, jt_k)
            mQ_temp(:, :, jt_k)  = mQ_aver(:, :, jt_k)
            mR_temp(:, jt_k)  = mR_aver(:, jt_k)
            fs_temp(:, jt_k)  = fs_aver(:, jt_k)
            fsQ_temp(:, :, jt_k) = fsQ_aver(:, :, jt_k)
            fm_temp(:, :, jt_k)  = fm_aver(:, :, jt_k)
            fmQ_temp(:, :, :, jt_k) = fmQ_aver(:, :, :, jt_k)
            fmR_temp(:, :, jt_k) = fmR_aver(:, :, jt_k)
        else
            pQ_aver(:, jt_k)  = 0.
            mQ_aver(:, :, jt_k)  = 0.
            mR_aver(:, jt_k)  = 0.
            fs_aver(:, jt_k)  = 0.
            fsQ_aver(:, :, jt_k) = 0.
            fm_aver(:, :, jt_k)  = 0.
            fmQ_aver(:, :, :, jt_k) = 0.
            fmR_aver(:, :, jt_k) = 0.
        end if
    end subroutine

    subroutine get_flux(hr)
        ! Calculates the fluxes in the given the curent state
        implicit none
        real(8), intent(in) :: hr
        integer iq, s, ip, leftbreakpt
        !call f_debug('get_flux', (/hr/))

        ! Use the SAS function lookup table to convert age-rank storage to the fraction of discharge of age T at each t

        ! First get the cumulative age-ranked storage
        if (iT_k==0) then
            if (hr==0) then
                STcum_top(jt_k) = jt_k
                PQcum_top(:, jt_k) = 0
                leftbreakpt_top(:, jt_k) = -1
                STcum_bot(jt_k) = 0
                PQcum_bot(:, jt_k) = 0
                leftbreakpt_bot(:, jt_k) = -1
                pQ_temp(:, jt_k) = 0
            else
                STcum_top(jt_k) = 0
                PQcum_top(:, jt_k) = 0
                leftbreakpt_top(:, jt_k) = -1
                STcum_bot(jt_k) = 0 + sT_temp(jt_k) * hr
                call get_SAS(STcum_bot(jt_k), PQcum_bot(:, jt_k), leftbreakpt_bot(:, jt_k), jt_k)
                pQ_temp(:, jt_k) = (PQcum_bot(:, jt_k) - PQcum_top(:, jt_k)) / hr
            end if
        else
            STcum_top(jt_k) = STcum_prev(jt_k) * (1-hr/h) + (STcum_prev(jt_k+1) + sT_end(jt_k) * h) * (hr/h)
            call get_SAS(STcum_top(jt_k), PQcum_top(:, jt_k), leftbreakpt_top(:, jt_k), jt_k)
            STcum_bot(jt_k) = STcum_top(jt_k) + sT_temp(jt_k) * h
            call get_SAS(STcum_bot(jt_k), PQcum_bot(:, jt_k), leftbreakpt_bot(:, jt_k), jt_k)
            pQ_temp(:, jt_k) = (PQcum_bot(:, jt_k) - PQcum_top(:, jt_k)) / h
        end if

        do iq = 0, numflux - 1
            if (sT_temp(jt_k)==0) then
                pQ_temp(iq, jt_k) = 0
            end if
        end do

        ! Solute mass flux accounting

        do iq = 0, numflux - 1
            do s = 0, numsol - 1

                ! Get the mass flux out
                if (sT_temp(jt_k)>0) then
                    mQ_temp(iq, s, jt_k) = mT_temp( s, jt_k) * alpha_ts(iq, s, jt) * Q_ts(iq, jt) &
                            * pQ_temp(iq, jt_k) / sT_temp(jt_k)

                    ! unless there is nothing in storage
                else
                    mQ_temp(iq, s, jt_k) = 0.

                end if
            enddo
        enddo

        ! Reaction mass accounting
        ! If there are first-order reactions, get the total mass rate
        mR_temp(:, jt_k) = k1_ts(:, jt) * (C_eq_ts(:, jt) * sT_temp(jt_k) - mT_temp(:, jt_k))

        fs_temp(:, jt_k) = 0.
        fsQ_temp(:, :, jt_k) = 0.
        if (sT_temp( jt_k)>0) then
            do iq = 0, numflux - 1
                fsQ_temp(:, iq, jt_k) = fsQ_temp(:, iq, jt_k) + ds_temp(:, jt_k) * pQ_temp(iq, jt_k) * Q_ts(iq, jt) / sT_temp(jt_k)
            end do
        end if
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                ! sensitivity to point before the start
                if ((leftbreakpt_top(ic, jt_k)>=0).and.(leftbreakpt_top(ic, jt_k)<numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_top(ic, jt_k)
                    !call f_debug('iq, ic, ip, jt_k ', (/iq*one8, ic*one8, ip*one8, jt_k*one8/))
                    dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                    dP = P_list(ip+1, jt) - P_list(ip, jt)
                    !call f_debug('dP/dS start    ', (/dP/dS/))
                    fs_temp(ip, jt_k) = fs_temp(ip, jt_k) &
                            + dP / (dS*dS) * sT_temp(jt_k) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point after the end
                if ((leftbreakpt_bot(ic, jt_k)+1>0).and.(leftbreakpt_bot(ic, jt_k)+1<=numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, jt_k) + 1
                    !call f_debug('iq, ic, ip, jt_k ', (/iq*one8, ic*one8, ip*one8, jt_k*one8/))
                    dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                    dP = P_list(ip, jt) - P_list(ip-1, jt)
                    !call f_debug('dP/dS end      ', (/dP/dS/))
                    fs_temp(ip, jt_k) = fs_temp(ip, jt_k) &
                            - dP / (dS*dS) * sT_temp(jt_k) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point within
                if (leftbreakpt_bot(ic, jt_k)>leftbreakpt_top(ic, jt_k)) then
                    call f_debug('leftbreakpt_bot, _start', &
                            (/leftbreakpt_bot(ic, jt_k)*one8, leftbreakpt_top(ic, jt_k)*one8/))
                    do leftbreakpt=leftbreakpt_top(ic, jt_k)+1, leftbreakpt_bot(ic, jt_k)
                        ip = breakpt_index_list(ic) + leftbreakpt
                        call f_debug('iq, ic, ip, jt_k ', (/iq*one8, ic*one8, ip*one8, jt_k*one8/))
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
                        fs_temp(ip, jt_k) = fs_temp(ip, jt_k) &
                                - (dPe/dSe - dPs/dSs) / h * weights_ts(ic, jt) * Q_ts(iq, jt)
                    end do
                end if
            end do
        end do
        fm_temp(:, :, jt_k) = 0
        fmQ_temp(:, :, :, jt_k) = 0
        do iq = 0, numflux - 1
            do ip = 0, numbreakpt_total - 1
                if (sT_temp(jt_k)>0) then
                    fmQ_temp(ip, iq, :, jt_k) = fmQ_temp(ip, iq, :, jt_k) &
                            + dm_temp(ip, :, jt_k) * alpha_ts(iq,:, jt) * Q_ts(iq, jt) * pQ_temp(iq, jt_k) / sT_temp(jt_k)
                end if
                fmR_temp(ip, :, jt_k) = fmR_temp(ip, :, jt_k) &
                        + k1_ts(:, jt) * (C_eq_ts(:, jt) * ds_temp(ip, jt_k) - dm_temp(ip, :, jt_k))
            end do
        end do
        do iq = 0, numflux - 1
            do ic = component_index_list(iq), component_index_list(iq+1) - 1
                ! sensitivity to point before the start
                if ((leftbreakpt_top(ic, jt_k)>=0).and.(leftbreakpt_top(ic, jt_k)<numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_top(ic, jt_k)
                    dS = SAS_lookup(ip+1, jt) - SAS_lookup(ip, jt)
                    dP = P_list(ip+1, jt) - P_list(ip, jt)
                    fm_temp(ip, :, jt_k) = fm_temp(ip, :, jt_k) &
                            + dP / (dS*dS) * mT_temp(:, jt_k)&
                                    * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point after the end
                if ((leftbreakpt_bot(ic, jt_k) + 1>0).and.(leftbreakpt_bot(ic, jt_k) + 1<=numbreakpt_list(ic)-1)) then
                    ip = breakpt_index_list(ic) + leftbreakpt_bot(ic, jt_k) + 1
                    dS = SAS_lookup(ip, jt) - SAS_lookup(ip-1, jt)
                    dP = P_list(ip, jt) - P_list(ip-1, jt)
                    fm_temp(ip, :, jt_k) = fm_temp(ip, :, jt_k) &
                            - dP / (dS*dS) * mT_temp(:, jt_k)&
                                    * alpha_ts(iq, :, jt) * weights_ts(ic, jt) * Q_ts(iq, jt)
                end if
                ! sensitivity to point within
                if (leftbreakpt_bot(ic, jt_k)>leftbreakpt_top(ic, jt_k)) then
                    do leftbreakpt = leftbreakpt_top(ic, jt_k)+1, leftbreakpt_bot(ic, jt_k)
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
                        fm_temp(ip, :, jt_k) = fm_temp(ip, :, jt_k) &
                                - (dPe/dSe - dPs/dSs) * mT_temp(:, jt_k) / sT_temp(jt_k) / h &
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
        sT_temp(jt_k) = sT_start(jt_k) ! Initial value
        mT_temp(:, jt_k) = mT_start(:, jt_k) + mR_temp(:, jt_k) * hr ! Initial value + reaction
        ! Fluxes in & out
        if (iT_k == 0) then
            sT_temp(jt_k) = sT_temp(jt_k) + J_ts(jt) * hr / h
            mT_temp(:, jt_k) = mT_temp(:, jt_k) + J_ts(jt) * C_J_ts(:, jt) * (hr/h)
        end if
        sT_temp(jt_k) = sT_temp(jt_k) - sum(Q_ts(:, jt) * pQ_temp(:, jt_k)) * hr
        mT_temp(:, jt_k) = mT_temp(:, jt_k) - sum(mQ_temp(:, :, jt_k), dim=1) * hr
        if (sT_temp(jt_k)<0) then
            call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
        end if

        ! Calculate new parameter sensitivity
        ds_temp(:, jt_k) = ds_start(:, jt_k) - fs_temp(:, jt_k) * hr - sum(fsQ_temp(:, :, jt_k), dim=2) * hr
        dm_temp(:, :, jt_k) = dm_start(:, :, jt_k) - fm_temp(:, :, jt_k) * hr - sum(fmQ_temp(:, :, :, jt_k), dim=3) * hr &
                + fmR_temp(:, :, jt_k) * hr

    end subroutine new_state

end subroutine solve
