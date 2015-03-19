!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-10-13 15:24:24 pbrowne>
!!!
!!!    module to hold all the information to control the the main program
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> module pf_control holds all the information to control the main program
module pf_control
    implicit none
    type, public :: pf_control_type
        integer :: nens !<the total number of ensemble members
        real(kind = kind(1.0D0)), allocatable, dimension(:) :: weight !< the negative log of the weights of the particles
        integer :: time_obs !< the number of observations we will assimilate
        integer :: time_bwn_obs !< the number of model timesteps between observations
        real(kind = kind(1.0D0)) :: nudgefac !< the nudging factor
        logical :: gen_data !< true generates synthetic obs for a twin experiment
        logical :: gen_Q !< true attempts to build up \f$Q\f$ from
        !<long model run
        logical :: human_readable !< unused
        integer :: timestep = 0 !< the current timestep as the model progresses
        real(kind = kind(1.0D0)), allocatable, dimension(:,:) :: psi !< state vector of ensemble members on this mpi process
        real(kind = kind(1.0D0)), allocatable, dimension(:) :: mean !< mean state vector
        real(kind = kind(1.0D0)) :: nfac !< standard deviation of normal distribution in mixture density
        real(kind = kind(1.0D0)) :: ufac !< half width of the uniform distribution in mixture density
        real(kind = kind(1.0D0)) :: efac
        real(kind = kind(1.0D0)) :: keep !< proportion of particles to keep in EWPF EW step
        real(kind = kind(1.0D0)) :: time !< dunno
        real(kind = kind(1.0D0)) :: Qscale !< scalar to
        !!multiply Q by
        real(kind = kind(1.0d0)) :: rho !< enkf inflation factor
        !! so that \f$P_f =
        !! (1+\rho)P_f\f$
        real(kind = kind(1.0d0)) :: len !< R localisation
        !! length scale
        integer :: couple_root !< empire master processor
        logical :: use_talagrand !< switch if true outputs rank histograms
        logical :: use_weak !< switch unused
        logical :: use_mean !< switch if true outputs ensemble mean
        logical :: use_var !< switch if true outputs ensemble variance
        logical :: use_traj !< switch if true outputs trajectories
        logical :: use_rmse !< switch if true outputs Root Mean Square Errors
        integer, dimension(:,:), allocatable :: talagrand !< storage for rank histograms
        integer :: count !< number of ensemble members associated with this MPI process
        integer, allocatable, dimension(:) :: particles !< particles associates with this MPI process
        character(2) :: filter!type !< which filter to use
        !< currently this has a number of
        !<options:
        !< - SE -- a stochastic ensemble
        !< - SI -- the SIR filter
        !< - ET -- the L-ETKF
        !< - EW -- the Equivalent Weights
        !< particle filter
        character(1) :: init !< which method to initialise ensemble
        !< currently this has a number of
        !< options:
        !< - N -- perturb around the model
        !< initial conditions with random noise
        !< distributed \f$\mathcal{N}(0,I)\f$
        !< - P -- perturb around the model
        !< initial conditions with random noise
        !< distributed \f$\mathcal{N}(0,Q)\f$
        !< - R -- read model states from
        !< rstrt folder where each ensemble member
        !< is stored in the file rstrt/##.state
        !< - S -- read model states from
        !< start folder where each ensemble member
        !< is stored in the file start/##.state

    end type pf_control_type
    type(pf_control_type), save :: pf !< the derived data type holding all controlling data

    TYPE, PUBLIC :: model_params_type
        INTEGER :: obs_dim !< size of the observation space
        INTEGER :: state_dim !< dimension of the model
    END TYPE model_params_type
    TYPE(model_params_type), SAVE :: model_params

contains
    !> subroutine to ensure pf_control data is ok
    subroutine set_pf_controls
        integer :: ios
        write(6, '(A)') 'Opening pf_parameters.dat'

        call parse_pf_parameters

        pf % efac = 0.001/pf % nens
        write(6, '(A)') 'pf_parameters.dat successfully read to control pf code.'
        call flush(6)
        if (pf % human_readable .and. pf % gen_data) then
            open(64, file = 'pf_data', iostat = ios, action = 'read', status = 'replace')
            if (ios .ne. 0) stop 'Error checking pf_data'
            close(64)
        end if




    end subroutine set_pf_controls


    !>subroutine to read the namelist file and save it to pf datatype
    !!Here we read pf_parameters.dat
    !!
    !! pf_parameters.dat is a fortran namelist file. As such, within
    !! it there must be a line beginning
    !!
    !! &pf_params
    !!
    !! To make it (probably) work, ensure there is a forward slash on
    !! the penultimate line and a blank line to end the file
    !!
    !! This is just the fortran standard for namelists though.
    !!
    !!
    !! On to the content...in any order, the pf_parameters.dat may
    !! contain the following things:
    !! 
    !! Integers:
    !! - \link pf_control::pf_control_type::time_obs time_obs \endlink
    !! - \link pf_control::pf_control_type::time_bwn_obs
    !! time_bwn_obs\endlink
    !!
    !! Reals, double precision:
    !! - \link pf_control::pf_control_type::nudgefac nudgefac\endlink
    !! - \link pf_control::pf_control_type::nfac nfac\endlink
    !! - \link pf_control::pf_control_type::ufac ufac\endlink
    !! - \link pf_control::pf_control_type::Qscale Qscale \endlink
    !! - \link pf_control::pf_control_type::keep keep  \endlink
    !! - \link pf_control::pf_control_type::rho rho  \endlink
    !! - \link pf_control::pf_control_type::len len  \endlink
    !!
    !! 2 Characters:
    !! - \link pf_control::pf_control_type::filter filter\endlink
    !!
    !! 1 Character:
    !! - \link pf_control::pf_control_type::init init\endlink
    !!
    !! Logicals:
    !! - \link pf_control::pf_control_type::gen_Q gen_Q\endlink
    !! - \link pf_control::pf_control_type::gen_data gen_data\endlink
    !! - \link pf_control::pf_control_type::use_talagrand use_talagrand\endlink
    !! - \link pf_control::pf_control_type::use_weak use_weak\endlink
    !! - \link pf_control::pf_control_type::use_var use_var\endlink
    !! - \link pf_control::pf_control_type::use_traj use_traj\endlink
    !! - \link pf_control::pf_control_type::use_rmse use_rmse\endlink
    !! - \link pf_control::pf_control_type::human_readable human_readable\endlink
    subroutine parse_pf_parameters
        implicit none
        character(*), parameter :: filename = 'pf_parameters.dat'
        integer :: ios

        integer :: time_obs = -1
        integer :: time_bwn_obs = -1
        real(kind = kind(1.0D0)) :: nudgefac = -1.0d0
        logical :: gen_data, gen_Q, human_readable
        real(kind = kind(1.0D0)) :: nfac = -1.0d0
        real(kind = kind(1.0D0)) :: ufac = -1.0d0
        real(kind = kind(1.0D0)) :: Qscale = -1.0d0
        real(kind = kind(1.0D0)) :: rho = 0.0d0
        real(kind = kind(1.0d0)) :: len = -1.0d0
        real(kind = kind(1.0D0)) :: keep
        logical :: use_talagrand, use_weak, use_mean, use_var, use_traj&
        &, use_rmse
        character(2) :: filter = '++'
        character(1) :: init = '+'

        namelist/pf_params/time_obs, time_bwn_obs, &
        &nudgefac, &
        &gen_data, gen_Q, &
        &human_readable, &
        &nfac, &
        &keep, &
        &ufac, &
        &Qscale, &
        &rho, &
        &len, &
        &use_talagrand, use_weak, use_mean, use_var, use_traj, use_rmse, &
        &filter, &
        &init


        gen_data = .false.
        gen_Q = .false.
        human_readable = .false.
        use_talagrand = .false.
        use_weak = .false.
        use_mean = .false.
        use_var = .false.
        use_traj = .false.
        use_rmse = .false.

        open(32, file = filename, iostat = ios, action = 'read'&
                &, status = 'old')
        if (ios .ne. 0) stop 'Cannot open pf_parameters.dat'
        read(32, nml = pf_params)
        !      print*,time_obs,time_bwn_obs,nudgefac,gen_data
        close(32)

        if (time_obs .gt. - 1) then
            print*, 'read time_obs = ', time_obs
            pf % time_obs = time_obs
        end if
        if (time_bwn_obs .gt. - 1) then
            print*, 'read time_bwn_obs = ', time_bwn_obs
            pf % time_bwn_obs = time_bwn_obs
        end if
        if (nudgefac .gt. - 1.0d0) then
            print*, 'read nudgefac = ', nudgefac
            pf % nudgefac = nudgefac
        end if
        !logical :: gen_data,gen_Q,human_readable
        if (nfac .gt. - 1.0d0) then
            print*, 'read nfac = ', nfac
            pf % nfac = nfac
        end if

        if (keep .gt. - 1.0d0) then
            print*, 'read keep = ', keep
            pf % keep = keep
        end if

        if (ufac .gt. - 1.0d0) then
            print*, 'read ufac = ', ufac
            pf % ufac = ufac
        end if

        !real(kind=kind(1.0D0)) :: ufac=-1.0d0
        if (Qscale .gt. - 1.0d0) then
            print*, 'read Qscale = ', Qscale
            pf % Qscale = Qscale
        end if


        if (rho .gt. 0.0d0) then
            print*, 'read rho = ', rho
            pf % rho = rho
        elseif (rho .lt. 0.0d0) then
            print*, 'read rho = ', rho, ' WARNING ABOUT THAT ONE! rho is normally positive'
            pf % rho = rho
        else
            pf % rho = rho
        end if


        if (len .ge. 0.0d0) then
            print*, 'read len = ', len
            pf % len = len
        else
            pf % len = len
        end if


        !logical ::
        !use_talagrand,use_weak,use_mean,use_var,use_traj,use_rmse


        !ensure that if we are generating the data then the EWPF is selected
        if (gen_data) then
            filter = 'EW'
            end if


            if (filter .ne. '++') then
                print*, 'read type = ', filter
                pf % filter = filter
            end if





            !let us verify pf%filter
            if (pf % filter .eq. 'EW') then
                print*, 'Running the equivalent weights particle filter'
            elseif (pf % filter .eq. 'SE') then
                print*, 'Running a stochastic ensemble'
            elseif (pf % filter .eq. 'SI') then
                print*, 'Running the SIR particle filter'
            elseif (pf % filter .eq. 'ET') then
                print*, 'Running the Ensemble Transform Kalman Filter'
                !print*,'Error: The ETKF is not implemented here'
                !stop
            elseif (pf % filter .eq. 'EA') then
                print*, 'Running the Ensemble Adjustment Kalman Filter'
                print*, 'Error: The EAKF is not implemented here yet'
                stop
            else
                print*, 'Error: Incorrect filter type selected'
                print*, 'Please ensure that pf%filter in pf_parameters.dat is either:'
                print*, 'EW                  the equivalent weights particle filter'
                print*, 'SE                  a stochastic ensemble'
                print*, 'SI                  the SIR particle filter'
                print*, 'ET                  the Ensemble Transform Kalman Filter'
                print*, 'EA                  the Ensemble Adjustment Kalman Filter'
                stop
            end if


            if (init .ne. '+') then
                print*, 'read init = ', init
                pf % init = init
            end if

            pf % gen_data = gen_data
            pf % gen_Q = gen_Q
            pf % human_readable = human_readable
            pf % use_talagrand = use_talagrand
            pf % use_weak = use_weak
            pf % use_mean = use_mean
            pf % use_var = use_var
            pf % use_traj = use_traj
            pf % use_rmse = use_rmse





        end subroutine parse_pf_parameters


        !> subroutine to allocate space for the filtering code
        subroutine allocate_pf
            use sizes
            use histogram_data
            integer :: st
            allocate(pf % weight(pf % nens), stat = st)
            if (st .ne. 0) stop 'Error in allocating pf%weight'
            pf % weight = -log(1.0D0/pf % nens)
            allocate(pf % psi(state_dim, pf % count), stat = st)
            if (st .ne. 0) stop 'Error in allocating pf%psi'

            if (pf % use_talagrand) then
                allocate(pf % talagrand(rhn_n, pf % nens + 1), stat = st)
                if (st .ne. 0) stop 'Error in allocating pf%talagrand'
                pf % talagrand = 0
            end if

            !      allocate(pf%particles(pf%count),stat=st)
            !      if(st .ne. 0) stop 'Error in allocating pf%particles'

        end subroutine allocate_pf

        !> subroutine to deallocate space for the filtering code
        subroutine deallocate_pf
            deallocate(pf % weight)
            deallocate(pf % psi)
            if (allocated(pf % talagrand)) deallocate(pf % talagrand)
            deallocate(pf % particles)
        end subroutine deallocate_pf


    end module pf_control
