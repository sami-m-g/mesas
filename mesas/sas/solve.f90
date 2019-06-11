! -*- f90 -*-
      subroutine f_solve(J_ts, Q_ts, SAS_lookup, P_list, STcum_init_ts, dt,  &
                 verbose, debug, full_outputs,&
                 CS_init_ts, C_J_ts, alpha_ts, k1_ts, C_eq_ts, C_old, &
                 n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_list, nP_total, &
                 STcum_ts, PQcum_ts, WaterBalance_ts, &
                 MScum_ts, MQcum_ts, MRcum_ts, C_Q_ts, SoluteBalance_ts)
      implicit none
      real(8), intent(in), dimension(0:timeseries_length-1) :: J_ts
      real(8), intent(in), dimension(0:timeseries_length-1, 0:numflux-1) :: Q_ts
      real(8), intent(in), dimension(0:nP_total-1, 0:timeseries_length-1) :: SAS_lookup
      real(8), intent(in), dimension(0:nP_total-1) :: P_list
      real(8), intent(in), dimension(0:max_age) :: STcum_init_ts
      real(8), intent(in) :: dt
      logical, intent(in) :: verbose, debug, full_outputs
      real(8), intent(in), dimension(0:max_age-1,0:numsol-1) :: CS_init_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: C_J_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numflux-1,0:numsol-1) &
                 :: alpha_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: k1_ts
      real(8), intent(in), dimension(0:timeseries_length-1,0:numsol-1) :: C_eq_ts
      real(8), intent(in), dimension(0:numsol-1) :: C_old
      integer, intent(in) ::  n_substeps, numflux, numsol, max_age, &
                 timeseries_length, nP_total
      integer, intent(in), dimension(0:numflux-1) :: nP_list
      real(8), intent(out), dimension(0:timeseries_length-1, 0:numflux-1, 0:numsol-1) :: C_Q_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length) :: STcum_ts
      real(8), intent(out), dimension(0:max_age-1, 0:timeseries_length-1) :: &
                 WaterBalance_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numflux-1) :: PQcum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numsol-1) :: MScum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numflux-1, 0:numsol-1) :: MQcum_ts
      real(8), intent(out), dimension(0:max_age, 0:timeseries_length, 0:numsol-1) :: MRcum_ts
      real(8), intent(out), dimension(0:max_age-1, 0:timeseries_length-1, 0:numsol-1) :: SoluteBalance_ts
      integer :: k, i
      real(8) :: h, P_old
      integer, dimension(0:numflux) :: iP_list
      real(8), dimension(0:max_age*n_substeps) :: STcum_init
      real(8), dimension(0:max_age*n_substeps, 0:numflux-1) :: PQcum_init
      real(8), dimension(0:max_age*n_substeps, 0:numflux-1) :: PQcum_temp
      real(8), dimension(0:max_age*n_substeps) :: STcum_temp
      real(8), dimension(0:max_age*n_substeps-1) :: sT_start
      real(8), dimension(0:max_age*n_substeps-1) :: sT_temp
      real(8), dimension(0:max_age*n_substeps-1) :: sT_end
      real(8), dimension(0:max_age*n_substeps) :: STcum_end
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ1
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ2
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ3
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ4
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ_substep
      real(8), dimension(0:max_age-1, 0:numflux-1) :: pQ_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ1
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ2
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ3
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ4
      real(8), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_substep
      real(8), dimension(0:max_age-1, 0:numflux-1, 0:numsol-1) :: mQ_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR1
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR2
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR3
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR4
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR_substep
      real(8), dimension(0:max_age-1, 0:numsol-1) :: mR_ts
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_start
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_end
      real(8), dimension(0:max_age*n_substeps, 0:numsol-1) :: MScum_end
      real(8), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mS_temp
      real(8) :: one8
      character(len=128) :: tempdebugstring
      integer is, ie, iq, s, M
      one8 = 1.0
      call f_verbose('...Initializing arrays...')
      C_Q_ts(:, :, :) = 0.
      STcum_ts(:, :) = 0.
      WaterBalance_ts(:, :) = 0.
      PQcum_ts(:, :, :) = 0.
      MScum_ts(:, :, :) = 0.
      MQcum_ts(:, :, :, :) = 0.
      MRcum_ts(:, :, :) = 0.
      SoluteBalance_ts(:, :, :) = 0.
      STcum_init(:) = 0.
      PQcum_init(:, :) = 0.
      PQcum_temp(:,:) = 0.
      STcum_temp(:) = 0.
      sT_start(:) = 0.
      sT_temp(:) = 0.
      sT_end(:) = 0.
      STcum_end(:) = 0.
      pQ1(:, :) = 0.
      pQ2(:, :) = 0.
      pQ3(:, :) = 0.
      pQ4(:, :) = 0.
      pQ_substep(:, :) = 0.
      pQ_ts(:, :) = 0.
      mQ1(:, :, :) = 0.
      mQ2(:, :, :) = 0.
      mQ3(:, :, :) = 0.
      mQ4(:, :, :) = 0.
      mQ_substep(:, :, :) = 0.
      mQ_ts(:, :, :) = 0.
      mR1(:, :) = 0.
      mR2(:, :) = 0.
      mR3(:, :) = 0.
      mR4(:, :) = 0.
      mR_substep(:, :) = 0.
      mR_ts(:, :) = 0.
      mS_start(:, :) = 0.
      mS_end(:, :) = 0.
      MScum_end(:, :) = 0.
      mS_temp(:, :) = 0.
      iP_list = cumsum_int(nP_list, numflux)
      call f_debug('iP_list', one8*iP_list(:))
      M = max_age * n_substeps
      h = dt / n_substeps
      call f_verbose('...Setting initial conditions...')
      do i=0,max_age-1
        is = i*n_substeps
        ie = (i+1)*n_substeps-1
        sT_end(is:ie) = (STcum_init_ts(i+1) - STcum_init_ts(i))/n_substeps
        do s=0,numsol-1
          mS_end(is:ie,s) = (CS_init_ts(i,s) * (STcum_init_ts(i+1) - STcum_init_ts(i)))&
                          /n_substeps
        enddo
      enddo
      STcum_init = cumsum(sT_end, M)
      do iq=0,numflux-1
        call lookup(SAS_lookup(iP_list(iq):iP_list(iq+1)-1,0), P_list(iP_list(iq):iP_list(iq+1)-1), &
                    STcum_init, PQcum_init(:,iq), nP_list(iq), M+1)
      enddo
      STcum_ts(1:max_age,0) = STcum_init(n_substeps:M:n_substeps)
      do iq=0,numflux-1
        PQcum_ts(:,0,iq) = PQcum_init(0:M:n_substeps, iq)
      enddo
      do s=0,numsol-1
        MScum_end(:,s) = cumsum(mS_end(:,s), M)
        MScum_ts(:,0,s) = MScum_end(0:M:n_substeps, s)
      enddo
      call f_verbose('...Starting main loop...')
      do i=0,timeseries_length-1
        pQ_ts(:, :) = 0.
        mQ_ts(:, :, :) = 0.
        mR_ts(:, :) = 0.
        do k=0,n_substeps-1
          call f_debug('Timestep, Substep', (/ i*one8, k*one8/))
          sT_start(1:M-1) = sT_end(0:M-2)
          sT_start(0) = 0.
          mS_start(1:M-1,:) = mS_end(0:M-2,:)
          mS_start(0,:) = 0.
          sT_temp = sT_start
          mS_temp = mS_start
          call f_debug('RK', (/1._8/))
          call get_flux(pQ1, mQ1, mR1)
          call new_state(h/2, pQ1, mQ1, mR1)
          call f_debug('RK', (/2._8/))
          call get_flux(pQ2, mQ2, mR2)
          call new_state(h/2, pQ2, mQ2, mR2)
          call f_debug('RK', (/3._8/))
          call get_flux(pQ3, mQ3, mR3)
          call new_state(h, pQ3, mQ3, mR3)
          call f_debug('RK', (/4._8/))
          call get_flux(pQ4, mQ4, mR4)
          pQ_substep = (pQ1 + 2*pQ2 + 2*pQ3 + pQ4) / 6.
          mQ_substep = (mQ1 + 2*mQ2 + 2*mQ3 + mQ4) / 6.
          mR_substep = (mR1 + 2*mR2 + 2*mR3 + mR4) / 6.
          do iq=0,numflux-1
            if (Q_ts(i,iq)==0) then
              pQ_substep(:,iq) = 0.
            endif
          enddo
          call new_state(h, pQ_substep, mQ_substep, mR_substep)
          sT_end = sT_temp
          mS_end = mS_temp
          do iq=0,numflux-1
            call f_debug('pQ_substep', pQ_substep(:,iq))
          enddo
          do s=0,numsol-1
            call f_debug('mS_end', mS_end(:,s))
            call f_debug('mR_substep', mR_substep(:,s))
            do iq=0,numflux-1
              call f_debug('mQ_substep', mQ_substep(:,iq,s))
            enddo
          enddo
          do iq=0,numflux-1
            P_old = 1 - sum(pQ_substep(:,iq))
            call f_debug('P_old', (/P_old/))
            do s=0,numsol-1
              if (Q_ts(i,iq)>0) then
                C_Q_ts(i,iq,s) = C_Q_ts(i,iq,s) + sum(mQ_substep(:,iq,s)) / &
                 (dt * Q_ts(i,iq)) / n_substeps
                C_Q_ts(i,iq,s) = C_Q_ts(i,iq,s) + alpha_ts(i,iq,s) * C_old(s) * P_old/n_substeps
                call f_debug('C_Q_ts', C_Q_ts(i:i,iq,s))
              endif
            enddo
          enddo
          if (full_outputs) then
            do iq=0,numflux-1
              pQ_ts(0, iq) = pQ_ts(0, iq) &
                           + sum(pQ_substep(0:k, iq))/n_substeps
              pQ_ts(1:max_age-1, iq) = pQ_ts(1:max_age-1, iq) &
                   + sum(reshape(pQ_substep(k+1:M-n_substeps+k, iq),&
                                 [n_substeps,max_age-1]),1)/n_substeps
            enddo
            do iq=0,numflux-1
              do s=0,numsol-1
                mQ_ts(0, iq, s) = mQ_ts(0, iq, s) &
                           + sum(mQ_substep(0:k, iq, s))/n_substeps
                mQ_ts(1:max_age-1, iq, s) = mQ_ts(1:max_age-1, iq, s)&
                           + sum(reshape(mQ_substep(k+1:M-n_substeps+k, iq, s), &
                                 [n_substeps,max_age-1]),1)/n_substeps
              enddo
            enddo
            do s=0,numsol-1
              mR_ts(0, s) = mR_ts(0, s) &
                           + sum(mR_substep(0:k, s))/n_substeps
              mR_ts(1:max_age-1, s) = mR_ts(1:max_age-1, s) &
                           + sum(reshape(mR_substep(k+1:M-n_substeps+k, s), &
                                 [n_substeps,max_age-1]),1)/n_substeps
            enddo
          endif
        enddo
        if (full_outputs) then
          STcum_end = cumsum(sT_end, M)
          STcum_ts(1:max_age, i+1) = STcum_end(n_substeps:M:n_substeps)
          do iq=0,numflux-1
            PQcum_ts(0:max_age,i+1,iq) = cumsum(pQ_ts(0:max_age-1,iq), max_age)
          enddo
          WaterBalance_ts(1:max_age-1, i) = diff(STcum_ts(0:max_age-1, i), max_age) - &
                   diff(STcum_ts(1:max_age, i+1), max_age)
          WaterBalance_ts(0, i) = J_ts(i) * dt - STcum_ts(1, i+1)
          do iq=0,numflux-1
            WaterBalance_ts(1:max_age-1, i) = WaterBalance_ts(1:max_age-1, i) - dt &
                   * (Q_ts(i,iq) * diff(PQcum_ts(1:max_age,i+1,iq), max_age))
            WaterBalance_ts(0, i) = WaterBalance_ts(0, i) - dt * Q_ts(i,iq) * &
                   (PQcum_ts(1,i+1,iq) - PQcum_ts(0,i+1,iq))
          enddo
          do s=0,numsol-1
            MScum_end(:,s) = cumsum(mS_end(:,s), M)
            MScum_ts(:,i+1,s) = MScum_end(0:M:n_substeps,s)
            MRcum_ts(:,i+1,s) = cumsum(mR_ts(:,s), max_age)
            do iq=0,numflux-1
              MQcum_ts(:,i+1,iq,s) = cumsum(mQ_ts(:,iq,s),&
                   max_age)
            enddo
          enddo
          do s=0,numsol-1
            SoluteBalance_ts(1:max_age-1,i,s)=(&
                   diff(MScum_ts(0:max_age-1,i,s), max_age) - &
                   diff(MScum_ts(1:max_age,i+1,s), max_age) &
                   + dt * diff(MRcum_ts(1:max_age-1,i+1,s), max_age))
            SoluteBalance_ts(0,i,s) = (C_J_ts(i,s) * J_ts(i) * dt - MScum_ts(0,i+1,s) + &
                   dt * (MRcum_ts(0,i+1,s) - MRcum_ts(0,i+1,s)))
            do iq=0,numflux-1
              SoluteBalance_ts(0:max_age-1,i,s) = (SoluteBalance_ts(0:max_age-1,i,s&
                   ) - dt * diff(MQcum_ts(0:max_age-1,i+1,iq,s), max_age+1))
              SoluteBalance_ts(0,i,s) = (SoluteBalance_ts(0,i,s) - dt * (MQcum_ts(0,&
                   i +1,iq,s) - MQcum_ts(0,i+1,iq,s)))
            enddo
          enddo
        endif
        if (mod(i,1000).eq.1000) then
            write (tempdebugstring,*) '...Done ', char(i), &
                                  'of', char(timeseries_length)
            call f_verbose(tempdebugstring)
        endif
      enddo
      call f_verbose('...Finished...')

      contains

      subroutine f_debug(debugstring, debugdblepr)
      implicit none
      character(len=*), intent(in) :: debugstring
      real(8), dimension(:), intent(in) :: debugdblepr
      if (debug) then
          print *, debugstring, debugdblepr
          !if (size(debugdblepr)==0) then
          !    call dblepr(debugstring, -1, 0, 0)
          !else
          !    call dblepr(debugstring, -1, debugdblepr, size(debugdblepr))
          !endif
      endif
      end subroutine f_debug

      subroutine f_verbose(debugstring)
      implicit none
      character(len=*), intent(in) :: debugstring
      if (verbose) then
          print *, debugstring
          !call dblepr(debugstring, -1, 0, 0)
      endif
      end subroutine f_verbose

      subroutine lookup(xa, ya, x, y, na, n)
      implicit none
      integer, intent(in) ::  na, n
      real(8), intent(in), dimension(0:na-1) :: xa
      real(8), intent(in), dimension(0:na-1) :: ya
      real(8), intent(in), dimension(0:n-1) :: x
      real(8), intent(out), dimension(0:n-1) :: y
      integer :: i, j, i0
      real(8) :: dif, grad
      logical :: foundit
      i0 = 0
      do j=0,n-1
        if (x(j).le.xa(0)) then
            y(j) = ya(0)
        else
            foundit = .FALSE.
            do i=i0,na-1
                if (x(j).lt.xa(i)) then
                    i0 = i-1
                    foundit = .TRUE.
                    exit
                endif
            enddo
            if (.not. foundit) then
                y(j) = ya(na-1)
            else
                dif = x(j) - xa(i0)
                grad = (ya(i0+1)-ya(i0))/(xa(i0+1)-xa(i0))
                y(j) = ya(i0) + dif * grad
            endif
        endif
      enddo
      end subroutine

      function diff(arr, n)
      implicit none
      integer, intent(in) ::  n
      real(8), intent(in), dimension(0:n-1) :: arr
      real(8), dimension(0:n-1-1) :: diff
      integer :: i
      do i=0,n-1-1
        diff(i) = arr(i+1) - arr(i)
      enddo
      end function diff

      function cumsum_int(arr, n)
      implicit none
      integer, intent(in) ::  n
      integer, intent(in), dimension(0:n-1) :: arr
      integer, dimension(0:n) :: cumsum_int
      integer :: i
      cumsum_int(0) = 0
      do i=0,n-1
        cumsum_int(i+1) = cumsum_int(i) + arr(i)
      enddo
    end function cumsum_int

      function cumsum(arr, n)
      implicit none
      integer, intent(in) ::  n
      real(8), intent(in), dimension(0:n-1) :: arr
      real(8), dimension(0:n) :: cumsum
      integer :: i
      cumsum(0) = 0
      do i=0,n-1
        cumsum(i+1) = cumsum(i) + arr(i)
      enddo
      end function cumsum

      subroutine get_flux(pQ_temp, mQ_temp, mR_temp)
      implicit none
      real(8), intent(out), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ_temp
      real(8), intent(out), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_temp
      real(8), intent(out), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR_temp
      integer M, iq, s
      M = max_age * n_substeps
      STcum_temp = cumsum(sT_temp, M)
      call f_debug('get_flux', (/0._8/))
      call f_debug('sT_temp', sT_temp)
      call f_debug('STcum_temp', STcum_temp)
      do s=0,numsol-1
        call f_debug('mS_temp', mS_temp(:,s))
      enddo
      do iq=0,numflux-1
        call lookup(SAS_lookup(iP_list(iq):iP_list(iq+1)-1,i), P_list(iP_list(iq):iP_list(iq+1)-1), STcum_temp, &
                 PQcum_temp(:,iq), nP_list(iq), M+1)
        pQ_temp(:,iq) = diff(PQcum_temp(:,iq), M+1)
        call f_debug('PQcum_temp', PQcum_temp(:,iq))
        call f_debug('pQ_temp', pQ_temp(:,iq))
        do s=0,numsol-1
          where (sT_temp(:)>0)
            mQ_temp(:,iq,s) = mS_temp(:,s) / sT_temp(:) * alpha_ts(i,iq,s) &
                 * Q_ts(i,iq) * pQ_temp(:,iq)
          elsewhere
            mQ_temp(:,iq,s) = 0.
          end where
        enddo
      enddo
      do s=0,numsol-1
         mR_temp(:,s) = k1_ts(i,s) * (C_eq_ts(i,s) * sT_temp - mS_temp(:,s))
         call f_debug('mQ_temp', mQ_temp(:,0,s))
         call f_debug('mE_temp', mQ_temp(:,1,s))
         call f_debug('mR_temp', mR_temp(:,s))
      enddo
      end subroutine get_flux

      subroutine new_state(hr, pQ_temp, mQ_temp, mR_temp)
      real(8), intent(in) :: hr
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numflux-1) :: pQ_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numflux-1, 0:numsol-1) :: mQ_temp
      real(8), intent(inout), dimension(0:max_age*n_substeps-1, 0:numsol-1) :: mR_temp
      call f_debug('new_state', (/0._8/))
      M = max_age * n_substeps
      sT_temp = sT_start
      do iq=0,numflux-1
        sT_temp = sT_temp - Q_ts(i,iq) * pQ_temp(:,iq) * hr
      enddo
      sT_temp(0) = sT_temp(0) + J_ts(i) * hr
      do s=0,numsol-1
        call f_debug('mS_start NS', mS_start(:,s))
        call f_debug('mR_temp NS', mR_temp(:,s))
        do iq=0,numflux-1
          call f_debug('mQ_temp NS', mQ_temp(:,iq,s))
        enddo
      enddo
      do s=0,numsol-1
        mS_temp(:,s) = mS_start(:,s) + mR_temp(:,s) * hr
        do iq=0,numflux-1
          mS_temp(:,s) = mS_temp(:,s) - mQ_temp(:,iq,s) * hr
        enddo
        mS_temp(0,s) = mS_temp(0,s) + J_ts(i) * C_J_ts(i,s) * hr
      enddo
      do s=0,numsol-1
        call f_debug('mS_temp NS', mS_temp(:,s))
        call f_debug('C_J_ts NS', (/C_J_ts(i,s)/))
      enddo
      call f_debug('J_ts NS', (/J_ts(i)/))
      end subroutine new_state

      end subroutine f_solve
