! -*- f90 -*-
subroutine solve(J_ts, Q_ts, SAS_lookup, P_list, sT_init_ts, dt, &
        verbose, debug, warning, &
        mS_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
        n_substeps, numflux, numsol, max_age, &
        timeseries_length, nP_list, nP_total, &
        sT_ts, pQ_ts, WaterBalance_ts, &
        mS_ts, mQ_ts, mR_ts, C_Q_ts, SoluteBalance_ts)
    implicit none

    ! Start by declaring and initializing all the variables we will be using
    integer, intent(in) :: n_substeps, numflux, numsol, max_age, &
            timeseries_length, nP_total
    real(8), intent(in), dimension(0:timeseries_length - 1) :: J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1) :: Q_ts
    real(8), intent(in), dimension(0:nP_total - 1, 0:timeseries_length - 1) :: SAS_lookup
    real(8), intent(in), dimension(0:nP_total - 1, 0:timeseries_length - 1) :: P_list
    real(8), intent(in), dimension(0:max_age - 1) :: sT_init_ts
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning
    real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mS_init_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_J_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: k1_ts
    real(8), intent(in), dimension(0:timeseries_length - 1, 0:numsol - 1) :: C_eq_ts
    real(8), intent(in), dimension(0:numsol - 1) :: C_old
    integer, intent(in), dimension(0:numflux - 1) :: nP_list
    real(8), intent(out), dimension(0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: C_Q_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length) :: sT_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length, 0:numsol - 1) :: mS_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1) :: pQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: mR_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1) :: WaterBalance_ts
    real(8), intent(out), dimension(0:max_age - 1, 0:timeseries_length - 1, 0:numsol - 1) :: SoluteBalance_ts
    integer :: k, iT, jt, ik, jk, i_prev, js
    real(8) :: h
    real(8), dimension(0:timeseries_length - 1, 0:numflux - 1) :: P_old
    integer, dimension(0:numflux) :: iP_list
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: J_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: C_J_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: Q_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: alpha_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: C_eq_ss
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: k1_ss
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_temp
    real(8), dimension(0:timeseries_length * n_substeps, 0:numflux - 1) :: PQcum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_end
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_prev
    real(8), dimension(0:timeseries_length * n_substeps, 0:numsol - 1) :: MScum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_prev
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_prev
    real(8) :: one8, weight
    character(len = 128) :: tempdebugstring
    integer iq, s, M, N
    logical :: carryover

    call f_verbose('...Initializing arrays...')
    one8 = 1.0

    C_Q_ts(:, :, :) = 0.
    sT_ts(:, :) = 0.
    mS_ts(:, :, :) = 0.
    pQ_ts(:, :, :) = 0.
    mQ_ts(:, :, :, :) = 0.
    mR_ts(:, :, :) = 0.
    WaterBalance_ts(:, :) = 0.
    SoluteBalance_ts(:, :, :) = 0.
    P_old(:, :) = 0.
    iP_list(:) = 0
    STcum_start(:) = 0.
    STcum_temp(:) = 0.
    PQcum_start(:, :) = 0.
    PQcum_temp(:, :) = 0.
    sT_start(:) = 0.
    sT_temp(:) = 0.
    sT_end(:) = 0.
    sT_prev(:) = 0.
    MScum_start(:, :) = 0.
    pQ1(:, :) = 0.
    pQ2(:, :) = 0.
    pQ3(:, :) = 0.
    pQ4(:, :) = 0.
    pQ_aver(:, :) = 0.
    pQ_end(:, :) = 0.
    pQ_prev(:, :) = 0.
    mQ1(:, :, :) = 0.
    mQ2(:, :, :) = 0.
    mQ3(:, :, :) = 0.
    mQ4(:, :, :) = 0.
    mQ_aver(:, :, :) = 0.
    mQ_end(:, :, :) = 0.
    mQ_prev(:, :, :) = 0.
    mR1(:, :) = 0.
    mR2(:, :) = 0.
    mR3(:, :) = 0.
    mR4(:, :) = 0.
    mR_aver(:, :) = 0.
    mR_end(:, :) = 0.
    mR_prev(:, :) = 0.
    mS_start(:, :) = 0.
    mS_temp(:, :) = 0.
    mS_end (:, :) = 0.
    mS_prev (:, :) = 0.
    i_prev = -1

    ! The list of probabilities in each sas function is a 1-D array.
    ! iP_list gives the starting index of the probabilities (P) associated
    ! with each flux
    iP_list(0) = 0
    do iq = 0, numflux - 1
        iP_list(iq + 1) = iP_list(iq) + nP_list(iq)
    enddo
    !call f_debug('iP_list', one8 * iP_list(:))

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
    call f_debug('J_ss           ', J_ss)
    do s = 0, numsol - 1
        call f_debug('C_J_ss           ', C_J_ss(:, s))
    end do
    do iq = 0, numflux - 1
        call f_debug('Q_ss           ', Q_ss(:,iq))
    end do

    call f_verbose('...Setting initial conditions...')
    sT_ts(:, 0) = sT_init_ts
    do s = 0, numsol - 1
        mS_ts(:, 0, s) = mS_init_ts(:, s)
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
            mS_start(1:N - 1, :) = mS_end(0:N - 2, :)
            ! Initialize the value at t=0
            if (ik>0) then
                sT_start(0) = sT_init_ts(i_prev)
                mS_start(0, :) = mS_init_ts(i_prev, :)
            end if

            ! These will hold the evolving state variables
            ! They are global variables modified by the new_state function

            ! This is the Runge-Kutta 4th order algorithm

            call f_debug('RK', (/1._8/))
            sT_temp = sT_start
            mS_temp = mS_start
            call get_flux(0.0D0, sT_temp, mS_temp, pQ1, mQ1, mR1)
            call f_debug_blank()
            call f_debug('RK', (/2._8/))
            call new_state(h / 2, sT_temp, mS_temp, pQ1, mQ1, mR1)
            call get_flux(h / 2, sT_temp, mS_temp, pQ2, mQ2, mR2)
            call f_debug_blank()
            call f_debug('RK', (/3._8/))
            call new_state(h / 2, sT_temp, mS_temp, pQ2, mQ2, mR2)
            call get_flux(h / 2, sT_temp, mS_temp, pQ3, mQ3, mR3)
            call f_debug_blank()
            call f_debug('RK', (/4._8/))
            call new_state(h, sT_temp, mS_temp, pQ3, mQ3, mR3)
            call get_flux(h, sT_temp, mS_temp, pQ4, mQ4, mR4)

            ! Average RK4 estimated change in the state variables
            pQ_aver = (pQ1 + 2 * pQ2 + 2 * pQ3 + pQ4) / 6.
            mQ_aver = (mQ1 + 2 * mQ2 + 2 * mQ3 + mQ4) / 6.
            mR_aver = (mR1 + 2 * mR2 + 2 * mR3 + mR4) / 6.

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
            ! Update the state with the new estimates
            call new_state(h, sT_end, mS_end, pQ_aver, mQ_aver, mR_aver)

            ! output some debugging info if desired
            do iq = 0, numflux - 1
                call f_debug('pQ_aver        ', pQ_aver(:, iq))
            enddo
            do s = 0, numsol - 1
                !call f_debug('mS_end         ', mS_end(:, s))
                !call f_debug('mR_aver        ', mR_aver(:, s))
                do iq = 0, numflux - 1
                    !call f_debug('mQ_aver        ', mQ_aver(:, iq, s))
                enddo
            enddo

            ! Aggregate flux data from substep to timestep

            ! Get the timestep-averaged transit time distribution
            call f_debug_blank()
            call f_debug('## Aggregating ', (/iT*one8/))
            weight = 1.0 / n_substeps / n_substeps
            do jt = 0, timeseries_length - 1
                js = jt * n_substeps
                carryover = ((n_substeps>1) .and. (k>0) .and. (iT<max_age-1))
                call f_debug('## timestep    ', (/js*one8, (js+k)*one8, (js+n_substeps-1)*one8/))
                call f_debug('pQ_aver        ', pQ_aver((js+k):(js+n_substeps-1), 0))
                do iq = 0, numflux - 1
                    pQ_ts(iT, jt, iq) = pQ_ts(iT, jt, iq) + sum(pQ_aver((js+k):(js+n_substeps-1), iq)) * weight
                    if (carryover) then
                        pQ_ts(iT+1, jt, iq) = pQ_ts(iT+1, jt, iq) + sum(pQ_aver((js):(js+k-1), iq)) * weight
                    endif
                    do s = 0, numsol - 1
                        mQ_ts(iT, jt, iq, s) = mQ_ts(iT, jt, iq, s) + sum(mQ_aver((js+k):(js+n_substeps-1), iq, s)) * weight
                        if (carryover) then
                            mQ_ts(iT+1, jt, iq, s) = mQ_ts(iT+1, jt, iq, s) + sum(mQ_aver((js):(js+k-1), iq, s)) * weight
                        endif
                    enddo
                enddo
                do s = 0, numsol - 1
                    mR_ts(iT, jt, s) = mR_ts(iT, jt, s) + sum(mR_aver((js+k):(js+n_substeps-1), s)) * weight
                    if (carryover) then
                        mR_ts(iT+1, jt, s) = mR_ts(iT+1, jt, s) + sum(mR_aver((js):(js+k-1), s)) * weight
                    endif
                enddo
            enddo
            call f_debug('pQ_ts          ', pQ_ts(iT, :, 0))
            if (iT<max_age-1) then
                call f_debug('pQ_ts iT+1     ', pQ_ts(iT+1, :, 0))
            end if

            ! Extract substep state at timesteps
            ! age-ranked storage at the end of the timestep
            sT_ts(iT, 1:) = sT_ts(iT, 1:) + sT_end(n_substeps-1:N-1:n_substeps) / n_substeps
            ! Age-ranked solute mass
            do s = 0, numsol - 1
                mS_ts(iT, 1:, s) = mS_ts(iT, 1:, s) + mS_end(n_substeps-1:N-1:n_substeps, s) / n_substeps
                call f_debug('mS_ts          ', mS_ts(iT, :, s))
            enddo

            ! Update the cumulative instantaneous trackers
            call get_flux(h, sT_end, mS_end, pQ_end, mQ_end, mR_end)
            if (ik>0) then
                STcum_start(1:N)    = STcum_start(1:N)    + sT_end * h
                MScum_start(1:N, :) = MScum_start(1:N, :) + mS_end * h
                PQcum_start(1:N, :) = PQcum_start(1:N, :) + pQ_end * h
                STcum_start(0) = STcum_start(0)    + sT_init_ts(i_prev) * h
                MScum_start(0, :) = MScum_start(0, :) + mS_init_ts(i_prev, :) * h
                PQcum_start(0, :) = PQcum_start(0, :) + pQ1(0, :) * h
                call f_debug('PQcum_start 1  ', PQcum_start(:, 0))
                PQcum_start(:, :) = 0
                call get_SAS(STcum_start, PQcum_start, N+1)
                call f_debug('PQcum_start 2  ', PQcum_start(:, 0))
            end if

            sT_prev = sT_end
            mS_prev = mS_end
            pQ_prev = pQ_end
            mR_prev = mR_end
            mQ_prev = mQ_end

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
                SoluteBalance_ts(iT, :, s) = C_J_ts(:, s) * J_ts - mS_ts(iT, 1:, s)
            end do
        else
            SoluteBalance_ts(iT, :, :) = mS_ts(iT-1, 0:timeseries_length-1, :) &
                    - mS_ts(iT, 1:timeseries_length, :)
        end if
        ! Subtract timestep-averaged mass fluxes
        SoluteBalance_ts(iT, :, :) = SoluteBalance_ts(iT, :, :) - sum(mQ_ts(iT, :, :, :), DIM=2) * dt
        ! Reacted mass
        SoluteBalance_ts(iT, :, :) = SoluteBalance_ts(iT, :, :) + mR_ts(iT, :, :) * dt

        ! Print some updates
        if (mod(jt, 1000).eq.1000) then
            write (tempdebugstring, *) '...Done ', char(jt), &
                    'of', char(timeseries_length)
            call f_verbose(tempdebugstring)
        endif

    enddo ! End of main loop

    call f_verbose('...Finalizing...')

    ! get the old water fraction
    P_old = 1 - sum(pQ_ts, DIM=1) * dt
    call f_debug('P_old', (/P_old/))

    ! Estimate the outflow concentration
    do iq = 0, numflux - 1
        do s = 0, numsol - 1
            do iT = 0, max_age - 1
                call f_debug('mQ_ts          ', mQ_ts(iT, :, iq, s))
            enddo
        enddo
    enddo
    do iq = 0, numflux - 1
        do s = 0, numsol - 1

            where (Q_ts(:,iq)>0)
                ! From the age-ranked mass
                C_Q_ts(:, iq, s) = sum(mQ_ts(:, :, iq, s), DIM=1) / Q_ts(:,iQ)

                ! From the old water concentration
                C_Q_ts(:, iq, s) = C_Q_ts(:, iq, s) + alpha_ts(:, iq, s) * C_old(s) * P_old(:,iq)
            end where

        call f_debug('C_Q_ts         ', C_Q_ts(:, iq, s))

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
            1 format (A16, *(f10.5))
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


    subroutine get_SAS(STcum_toget, PQcum_toget, n_array)
        ! Call the sas function and get the transit time distribution
        integer, intent(in) :: n_array
        real(8), intent(in), dimension(0:n_array - 1) :: STcum_toget
        real(8), intent(out), dimension(0:n_array - 1, 0:numflux - 1) :: PQcum_toget
        ! Main lookup loop
        do iq = 0, numflux - 1
            do jt = 0, timeseries_length - 1
                jk = jt * n_substeps
                call f_debug('getting entries', (/jk*one8, (jk+n_substeps-1)*one8/))
                call lookup(&
                        SAS_lookup(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                        P_list(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                        STcum_toget(jk:jk+n_substeps-1), PQcum_toget(jk:jk+n_substeps-1, iq), &
                        nP_list(iq), n_substeps)
            end do
            if (n_array>N) then
                call f_debug('getting last entry', (/n_array*one8, N*one8/))
                jt = timeseries_length - 1
                call lookup(&
                        SAS_lookup(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                        P_list(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                        STcum_toget(N:N), PQcum_toget(N:N, iq), &
                        nP_list(iq), 1)
            end if
        enddo
        end subroutine get_SAS

    subroutine lookup(xa, ya, x, y, na, n)
        ! A simple lookup table
        implicit none
        integer, intent(in) :: na, n
        real(8), intent(in), dimension(0:na - 1) :: xa
        real(8), intent(in), dimension(0:na - 1) :: ya
        real(8), intent(in), dimension(0:n - 1) :: x
        real(8), intent(inout), dimension(0:n - 1) :: y
        integer :: i, j, i0
        real(8) :: dif, grad
        logical :: foundit
        i0 = 0
        do j = 0, n - 1
            if (x(j).le.xa(0)) then
                y(j) = ya(0)
            else if (x(j).ge.xa(na - 1)) then
                y(j) = ya(na - 1)
            else
                foundit = .FALSE.
                do i = i0, na - 1
                    if (x(j).lt.xa(i)) then
                        i0 = i - 1
                        foundit = .TRUE.
                        exit
                    endif
                enddo
                if (.not. foundit) then
                    y(j) = ya(na - 1)
                else
                    dif = x(j) - xa(i0)
                    grad = (ya(i0 + 1) - ya(i0)) / (xa(i0 + 1) - xa(i0))
                    y(j) = ya(i0) + dif * grad
                endif
            endif
        enddo
    end subroutine


    subroutine get_flux(hr, sT_in, mS_in, pQ_out, mQ_out, mR_out)
        ! Calculates the fluxes in the given the curent state
        implicit none
        real(8), intent(in) :: hr
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1) :: sT_in
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_in
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_out
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_out
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_out
        integer iq, s
        call f_debug('get_flux', (/hr/))

        ! Use the SAS function lookup table to convert age-rank storage to the fraction of discharge of age T at each t

        ! First get the cumulative age-ranked storage
        if (ik==0) then
            if (hr==0) then
                STcum_temp = 0
            else
                STcum_temp = 0 + sT_in * hr
            end if
        else
            STcum_temp = STcum_start(0:N-1) * (1-hr/h) + (STcum_start(1:N) + sT_prev * h) * (hr/h) + sT_in * h
        end if

        ! Print some debuggin info
        !call f_debug('get_flux', (/0._8/))
        call f_debug('STcum_start    ', STcum_start(:))
        call f_debug('sT_prev        ', sT_prev(:))
        call f_debug('sT_in          ', sT_in(:))
        call f_debug('STcum_temp     ', STcum_temp(:))
        do s = 0, numsol - 1
            call f_debug('mS_in          ', mS_in(:, s))
        enddo
        call f_debug_blank()

        call get_SAS(STcum_temp, PQcum_temp, N)

        ! Get the pdf form
        call f_debug('getting pQ_out ', (/0._8/))
        if (ik==0) then
            if (hr==0) then
                pQ_out = 0
            else
                pQ_out = (PQcum_temp - 0) / hr
            end if
        else
            pQ_out = (PQcum_temp - (PQcum_start(0:N-1, :) * (1-hr/h) + (PQcum_start(1:N, :) + pQ_prev * h ) * (hr/h))) /h
        end if
        do iq = 0, numflux - 1
            where (sT_in(:)==0)
                pQ_out(:, iq) = 0
            end where
        end do

        do iq = 0, numflux - 1
            call f_debug('PQcum_start    ', PQcum_start(:, iq))
            call f_debug('pQ_prev        ', pQ_prev(:, iq))
            call f_debug('PQcum_temp     ', PQcum_temp(:, iq))
            call f_debug('pQ_out         ', pQ_out(:, iq))
        end do
        call f_debug_blank()

        ! Solute mass flux accounting
        call f_debug('getting mQ_out', (/0._8/))
        do iq = 0, numflux - 1
            do s = 0, numsol - 1

                ! Get the mass flux out
                where (sT_in(:)>0)
                    mQ_out(:, iq, s) = mS_in(:, s) / sT_in(:) * alpha_ss(:, iq, s) &
                            * Q_ss(:, iq) * pQ_out(:, iq)

                    ! unless there is nothing in storage
                elsewhere
                    mQ_out(:, iq, s) = 0.

                end where
            enddo

        enddo

        ! Reaction mass accounting
        !call f_debug('getting mR_out', (/0._8/))
        do s = 0, numsol - 1

            ! If there are first-order reactions, get the total mass rate
            mR_out(:, s) = k1_ss(:, s) * (C_eq_ss(:, s) * sT_in(:) - mS_in(:, s))

            call f_debug('mQ_out         ', mQ_out(:, 0, s))
            !call f_debug('mE_out         ', mQ_out(:, 1, s))
            !call f_debug('mR_out         ', mR_out(:, s))

        enddo
        call f_debug('get_flux finished', (/0._8/))
        call f_debug_blank()

    end subroutine get_flux


    subroutine new_state(hr, sT_out, mS_out, pQ_in, mQ_in, mR_in)
        ! Calculates the state given the fluxes

        real(8), intent(in) :: hr
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1) :: sT_out
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_out
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_in
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_in
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_in

        call f_debug('new_state', (/0._8/))

        ! Calculate the new age-ranked storage
        sT_out = sT_start ! Initial value
        ! Fluxes in & out
        if (ik == 0) then
            sT_out = sT_out + J_ss * hr / h
        end if
        do iq = 0, numflux - 1
            sT_out = sT_out - Q_ss(:, iq) * hr * pQ_in(:, iq)
        enddo
        call f_debug('sT_start ns    ', sT_start(:))
        call f_debug('sT_out ns      ', sT_out(:))
        if (ANY(sT_out<0)) then
                call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
        end if

        ! Print some debugging info
        do s = 0, numsol - 1
            call f_debug('mS_start NS    ', mS_start(:, s))
            !call f_debug('mR_in NS       ', mR_in(:, s))
            do iq = 0, numflux - 1
                call f_debug('mQ_in NS       ', mQ_in(:, iq, s))
            enddo
        enddo

        ! Calculate the new age-ranked mass
        do s = 0, numsol - 1
            ! Initial value + reaction
            mS_out(:, s) = mS_start(:, s) + mR_in(:, s) * hr
            ! Flux in
            if (ik==0) then
                mS_out(:, s) = mS_out(:, s) + J_ss(:) * C_J_ss(:, s)
            end if
            ! Fluxes out
            do iq = 0, numflux - 1
                mS_out(:, s) = mS_out(:, s) - mQ_in(:, iq, s) * hr
            enddo
        enddo

        do s = 0, numsol - 1
            call f_debug('mS_out NS      ', mS_out(:, s))
            !call f_debug('C_J_ss NS      ', (/C_J_ss(:, s)/))
        enddo
        !call f_debug('J_ss NS        ', (/J_ss(:)/))
        call f_debug('new_state finished', (/0._8/))
        call f_debug_blank()

    end subroutine new_state

end subroutine solve
