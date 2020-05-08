! -*- f90 -*-
subroutine solve(J_ts, Q_ts, SAS_lookup, P_list, sT_init, dt, &
        verbose, debug, warning, &
        mS_init, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
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
    real(8), intent(in), dimension(0:max_age - 1) :: sT_init
    real(8), intent(in) :: dt
    logical, intent(in) :: verbose, debug, warning
    real(8), intent(in), dimension(0:max_age - 1, 0:numsol - 1) :: mS_init
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
    integer :: k, iT, jt, ik, jk
    real(8) :: h
    real(8), dimension(0:timeseries_length - 1, 0:numflux - 1) :: P_old
    integer, dimension(0:numflux) :: iP_list
    real(8), dimension(0:timeseries_length * n_substeps) :: STcum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: STcum_temp
    real(8), dimension(0:timeseries_length * n_substeps, 0:numflux - 1) :: PQcum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: PQcum_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_start
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1) :: sT_end
    real(8), dimension(0:timeseries_length * n_substeps, 0:numsol - 1) :: MScum_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR1
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR2
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR3
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR4
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_aver
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_end
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_start
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_temp
    real(8), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_end
    real(8) :: one8
    character(len = 128) :: tempdebugstring
    integer iq, s, M, N

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
    MScum_start(:, :) = 0.
    pQ1(:, :) = 0.
    pQ2(:, :) = 0.
    pQ3(:, :) = 0.
    pQ4(:, :) = 0.
    pQ_aver(:, :) = 0.
    pQ_end(:, :) = 0.
    mQ1(:, :, :) = 0.
    mQ2(:, :, :) = 0.
    mQ3(:, :, :) = 0.
    mQ4(:, :, :) = 0.
    mQ_aver(:, :, :) = 0.
    mQ_end(:, :, :) = 0.
    mR1(:, :) = 0.
    mR2(:, :) = 0.
    mR3(:, :) = 0.
    mR4(:, :) = 0.
    mR_aver(:, :) = 0.
    mR_end(:, :) = 0.
    mS_start(:, :) = 0.
    mS_temp(:, :) = 0.
    mS_end (:, :) = 0.

    ! The list of probabilities in each sas function is a 1-D array.
    ! iP_list gives the starting index of the probabilities (P) associated
    ! with each flux
    iP_list(0) = 0
    do iq = 0, numflux - 1
        iP_list(iq + 1) = iP_list(iq) + nP_list(iq)
    enddo
    call f_debug('iP_list', one8 * iP_list(:))

    ! modify the number of ages and the timestep by a facotr of n_substeps
    M = max_age * n_substeps
    N = timeseries_length * n_substeps
    h = dt / n_substeps

    ! The time boundary conditions are set with (in python) ST_init, which becomes STcum_init
    ! and mS_init, which becomes mS_init in fortran
    call f_verbose('...Setting initial conditions...')
    sT_ts(:, 0) = sT_init
    mS_ts(:, 0, s) = mS_init(:, s)

    call f_verbose('...Starting main loop...')
    do iT = 0, max_age - 1

        ! Start the substep loop
        do k = 0, n_substeps - 1

            ik = iT * n_substeps + k

            call f_debug('Agestep, Substep', (/ iT * one8, k * one8/))

            ! Copy the state variables from the end of the previous substep as the start of this one
            ! but shifted by one substep
            sT_start(1:N - 1) = sT_end(0:N - 2)
            mS_start(1:N - 1, :) = mS_end(0:N - 2, :)
            ! Initialize the value at t=0
            if (ik>0) then
                sT_start(0) = sT_init(ik - 1)
                mS_start(0, :) = mS_init(ik - 1, :)
            end if

            ! These will hold the evolving state variables
            ! They are global variables modified by the new_state function

            ! This is the Runge-Kutta 4th order algorithm

            call f_debug('RK', (/1._8/))
            sT_temp = sT_start
            mS_temp = mS_start
            call get_flux(0.0D0, sT_temp, mS_temp, pQ1, mQ1, mR1)
            call f_debug('RK', (/2._8/))
            call new_state(h / 2, sT_temp, mS_temp, pQ1, mQ1, mR1)
            call get_flux(h / 2, sT_temp, mS_temp, pQ2, mQ2, mR2)
            call f_debug('RK', (/3._8/))
            call new_state(h / 2, sT_temp, mS_temp, pQ2, mQ2, mR2)
            call get_flux(h / 2, sT_temp, mS_temp, pQ3, mQ3, mR3)
            call f_debug('RK', (/4._8/))
            call new_state(h, sT_temp, mS_temp, pQ3, mQ3, mR3)
            call get_flux(h, sT_temp, mS_temp, pQ4, mQ4, mR4)

            ! Average RK4 estimated change in the state variables
            pQ_aver = (pQ1 + 2 * pQ2 + 2 * pQ3 + pQ4) / 6.
            mQ_aver = (mQ1 + 2 * mQ2 + 2 * mQ3 + mQ4) / 6.
            mR_aver = (mR1 + 2 * mR2 + 2 * mR3 + mR4) / 6.

            ! zero out the probabilities if there is no outflux this timestep
            where (Q_ts==0)
                pQ_aver = 0.
            end where
            do s = 0, numsol - 1
                where (Q_ts==0)
                    mQ_aver(:, :, s) = 0.
                end where
            end do

            ! Update the cumulative instantaneous trackers before we overwrite the '_end' variables
            if (ik>0) then
                call get_flux(h, sT_end, mS_end, pQ_end, mQ_end, mR_end)
                STcum_start(1:N)    = STcum_start(1:N)    + sT_end * h
                MScum_start(1:N, :) = MScum_start(1:N, :) + mS_end * h
                PQcum_start(1:N, :) = PQcum_start(1:N, :) + pQ_end * h
                STcum_start(0) = STcum_start(0)    + sT_start(0) * h
                MScum_start(0, :) = MScum_start(0, :) + mS_start(0, :) * h
                PQcum_start(0, :) = PQcum_start(0, :) + pQ1(0, :) * h
            end if

            ! Update the state with the new estimates
            call new_state(h, sT_end, mS_end, pQ_aver, mQ_aver, mR_aver)

            ! output some debugging info if desired
            do iq = 0, numflux - 1
                call f_debug('pQ_aver', pQ_aver(:, iq))
            enddo
            do s = 0, numsol - 1
                call f_debug('mS_end', mS_end(:, s))
                call f_debug('mR_aver', mR_aver(:, s))
                do iq = 0, numflux - 1
                    call f_debug('mQ_aver', mQ_aver(:, iq, s))
                enddo
            enddo

            ! Aggregate flux data from substep to timestep

            ! Get the timestep-averaged transit time distribution
            do iq = 0, numflux - 1
                pQ_ts(iT, 0, iq) = pQ_ts(iT, 0, iq) &
                        + sum(pQ_aver(0:k, iq)) / n_substeps
                pQ_ts(iT, 1:, iq) = pQ_ts(iT, 1:, iq) &
                        + sum(reshape(pQ_aver(k + 1:N - n_substeps + k, iq), &
                                [n_substeps, timeseries_length - 1]), 1) / n_substeps
            enddo

            ! Get the timestep-averaged mass flux
            do iq = 0, numflux - 1
                do s = 0, numsol - 1
                    mQ_ts(iT, 0, iq, s) = mQ_ts(iT, 0, iq, s) &
                            + sum(mQ_aver(0:k, iq, s)) / n_substeps
                    mQ_ts(iT, 1:, iq, s) = mQ_ts(iT, 1:, iq, s)&
                            + sum(reshape(mQ_aver(k + 1:N - n_substeps + k, iq, s), &
                                    [n_substeps, timeseries_length - 1]), 1) / n_substeps
                enddo
            enddo

            ! Get the timestep-averaged reaction rates
            do s = 0, numsol - 1
                mR_ts(iT, 0, s) = mR_ts(iT, 0, s) &
                        + sum(mR_aver(0:k, s)) / n_substeps
                mR_ts(iT, 1:, s) = mR_ts(iT, 1:, s) &
                        + sum(reshape(mR_aver(k + 1:N - n_substeps + k, s), &
                                [n_substeps, timeseries_length - 1]), 1) / n_substeps
            enddo

        enddo

        ! Extract substep state at timesteps
        ! age-ranked storage at the end of the timestep
        sT_ts(iT, 1:) = sT_end(n_substeps-1:N-1:n_substeps)
        ! Age-ranked solute mass
        do s = 0, numsol - 1
            mS_ts(iT, 1:, s) = mS_end(n_substeps-1:N-1:n_substeps, s)
        enddo

        ! Print some updates
        if (mod(jt, 1000).eq.1000) then
            write (tempdebugstring, *) '...Done ', char(jt), &
                    'of', char(timeseries_length)
            call f_verbose(tempdebugstring)
        endif

    enddo ! End of main loop

    ! get the old water fraction
    P_old = 1 - sum(pQ_ts, DIM=1) * dt
    call f_debug('P_old', (/P_old/))

    ! Estimate the outflow concentration
    do s = 0, numsol - 1
        where (Q_ts>0)

            ! From the age-ranked mass
            C_Q_ts(:, :, s) = C_Q_ts(:, :, s) + sum(mQ_ts(:, :, :, s), DIM=1) * dt / Q_ts

            ! From the old water concentration
            C_Q_ts(:, :, s) = C_Q_ts(:, :, s) + alpha_ts(:, :, s) * C_old(s) * P_old

        end where

        call f_debug('C_Q_ts', C_Q_ts(:, iq, s))

    enddo

    ! Calculate a water balance based on these
    ! Difference of starting and ending age-ranked storage
    WaterBalance_ts(0, :) = J_ts * dt - sT_ts(0, 1:)
    WaterBalance_ts(1:, :) = sT_ts(:max_age-2, 0:timeseries_length-1) - sT_ts(1:max_age-1, 1:timeseries_length)
    ! subtract time-averaged water fluxes
    WaterBalance_ts = WaterBalance_ts - sum(pQ_ts(:max_age-1, :, :), DIM=3) * dt

    ! Calculate a solute balance based on these
    do s = 0, numsol - 1
        ! Difference of starting and ending age-ranked mass
        SoluteBalance_ts(0, :, s) = C_J_ts(:, s) * J_ts * dt - mS_ts(0, 1:, s)
        SoluteBalance_ts(1:, :, s) = mS_ts(0:max_age-2, 0:timeseries_length-1, s) &
                                   - mS_ts(1:max_age-1, 1:timeseries_length, s)
    enddo
    ! Subtract timestep-averaged mass fluxes
    SoluteBalance_ts = SoluteBalance_ts - sum(mQ_ts, DIM=3) * dt
    ! Reacted mass
    SoluteBalance_ts = SoluteBalance_ts + mR_ts * dt

    call f_verbose('...Finished...')

contains


    subroutine f_debug(debugstring, debugdblepr)
        ! Prints debugging information
        implicit none
        character(len = *), intent(in) :: debugstring
        real(8), dimension(:), intent(in) :: debugdblepr
        if (debug) then
            print *, debugstring, debugdblepr
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


    subroutine get_flux(hr, sT_temp, mS_temp, pQ_temp, mQ_temp, mR_temp)
        ! Calculates the fluxes in the given the curent state
        implicit none
        real(8), intent(in) :: hr
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_temp
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_temp
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_temp
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_temp
        integer iq, s

        ! Use the SAS function lookup table to convert age-rank storage to the fraction of discharge of age T at each t

        ! First get the cumulative age-ranked storage
        if (ik==0) then
            if (hr==0) then
                STcum_temp = 0
            else
                STcum_temp = 0 + sT_temp * hr
            end if
        else
            STcum_temp = STcum_start(0:N-1) * (1-hr/h) + (STcum_start(1:N) + sT_end * h) * (hr/h) + sT_temp * h
        end if

        ! Print some debuggin info
        call f_debug('get_flux', (/0._8/))
        call f_debug('sT_temp', sT_temp)
        call f_debug('STcum_temp', STcum_temp)
        do s = 0, numsol - 1
            call f_debug('mS_temp', mS_temp(:, s))
        enddo

        ! Main lookup loop
        do iq = 0, numflux - 1
            do jt = 0, timeseries_length - 1
                do k = 0, n_substeps - 1
                    jk = jt * n_substeps + k
                    ! Call the sas function and get the transit time distribution
                    call lookup(&
                            SAS_lookup(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                            P_list(iP_list(iq):iP_list(iq + 1) - 1, jt), &
                            STcum_temp(jk), PQcum_temp(jk, iq), &
                            nP_list(iq), N)
                end do
            end do
            call f_debug('PQcum_temp', PQcum_temp(:, iq))
        enddo

        ! Get the pdf form
        if (ik==0) then
            if (hr==0) then
                pQ_temp = 0
            else
                pQ_temp = (PQcum_temp - 0) / hr
            end if
        else
            pQ_temp = (PQcum_temp - (PQcum_start(0:N-1, :) * (1-hr/h) + (PQcum_start(1:N, :) + pQ_end * h ) * (hr/h))) /h
        end if
        call f_debug('pQ_temp', pQ_temp(:, iq))

        ! Solute mass flux accounting
        do iq = 0, numflux - 1
            do s = 0, numsol - 1

                ! Get the mass flux out
                where (sT_temp(:)>0)
                    mQ_temp(:, iq, s) = mS_temp(:, s) / sT_temp(:) * alpha_ts(:, iq, s) &
                            * Q_ts(:, iq) * pQ_temp(:, iq)

                    ! unless there is nothing in storage
                elsewhere
                    mQ_temp(:, iq, s) = 0.

                end where
            enddo

        enddo

        ! Reaction mass accounting
        do s = 0, numsol - 1

            ! If there are first-order reactions, get the total mass rate
            mR_temp(:, s) = k1_ts(:, s) * (C_eq_ts(:, s) * sT_temp(:) - mS_temp(:, s))

            call f_debug('mQ_temp', mQ_temp(:, 0, s))
            call f_debug('mE_temp', mQ_temp(:, 1, s))
            call f_debug('mR_temp', mR_temp(:, s))

        enddo

    end subroutine get_flux


    subroutine new_state(hr, sT_temp, mS_temp, pQ_temp, mQ_temp, mR_temp)
        ! Calculates the state given the fluxes

        real(8), intent(in) :: hr
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1) :: sT_temp
        real(8), intent(out), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mS_temp
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1) :: pQ_temp
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numflux - 1, 0:numsol - 1) :: mQ_temp
        real(8), intent(in), dimension(0:timeseries_length * n_substeps - 1, 0:numsol - 1) :: mR_temp

        call f_debug('new_state', (/0._8/))

        ! Calculate the new age-ranked storage
        sT_temp = sT_start ! Initial value
        ! Fluxes in & out
        if (ik == 0) then
            sT_temp = sT_temp + J_ts * hr / h
            do iq = 0, numflux - 1
                sT_temp = sT_temp - Q_ts(:, iq) * hr * pQ_temp(:, iq) * hr/h
            enddo
        else
            do iq = 0, numflux - 1
                sT_temp = sT_temp - Q_ts(:, iq) * hr * pQ_temp(:, iq)
            enddo
        end if
        if (ANY(sT_temp<0)) then
                call f_warning('WARNING: A value of sT is negative. Try increasing the number of substeps')
        end if

        ! Print some debugging info
        do s = 0, numsol - 1
            call f_debug('mS_start NS', mS_start(:, s))
            call f_debug('mR_temp NS', mR_temp(:, s))
            do iq = 0, numflux - 1
                call f_debug('mQ_temp NS', mQ_temp(:, iq, s))
            enddo
        enddo

        ! Calculate the new age-ranked mass
        do s = 0, numsol - 1
            ! Initial value + reaction
            mS_temp(:, s) = mS_start(:, s) + mR_temp(:, s) * hr
            ! Fluxes out
            do iq = 0, numflux - 1
                mS_temp(:, s) = mS_temp(:, s) - mQ_temp(:, iq, s) * hr
            enddo
            ! Flux in
            mS_temp(:, s) = mS_temp(:, s) + J_ts(:) * C_J_ts(:, s) * hr
        enddo

        do s = 0, numsol - 1
            call f_debug('mS_temp NS', mS_temp(:, s))
            call f_debug('C_J_ts NS', (/C_J_ts(:, s)/))
        enddo
        call f_debug('J_ts NS', (/J_ts(:)/))

    end subroutine new_state

end subroutine solve
