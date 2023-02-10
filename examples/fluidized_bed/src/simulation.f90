!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use lpt_class,         only: lpt
  use lowmach_class,     only: lowmach
  use vdscalar_class,    only: vdscalar
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get an LPT solver, a lowmach solver, and corresponding time tracker
  type(vdscalar), dimension(2), public :: sc
  type(lowmach),     public :: fs
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile,lptfile,tfile,adsfile,scfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rhof,dRHOdt
  real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
  real(WP), dimension(:,:,:,:), allocatable :: srcSClp
  real(WP) :: visc,rho,diffusivity,inlet_velocity

  !> Scalar indices
  integer, parameter :: ind_T  =1
  integer, parameter :: ind_CO2=2
  integer :: nscalar=2
  real(WP), dimension(2) :: SC_inlet

  !> Thermo-chemistry parameters
  real(WP), parameter :: Rcst=8.314_WP         !< J/(mol.K)
  real(WP), parameter :: W_N2  = 28.0135e-3_WP !< kg/mol
  real(WP), parameter :: W_CO2 = 44.0095e-3_WP !< kg/mol
  real(WP), parameter :: W_H2O = 18.0153e-3_WP !< kg/mol
  real(WP) :: Pthermo
  real(WP) :: mean

  !> Wallclock time for monitoring
  type :: timer
     real(WP) :: time_in
     real(WP) :: time
     real(WP) :: percent
  end type timer
  type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_sc,wt_rest

  ! !> Event for post-processing
  ! type(event) :: ppevt

contains

  !> Define here our equation of state - rho(P,T,Yk)
  !> This just updates sc%rho
  subroutine get_sc_rho
    real(WP) :: T,Y_CO2,Y_N2,Wmix
    integer :: i,j,k
    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
          do i=fs%cfg%imino_,fs%cfg%imaxo_
             if (fs%cfg%VF(i,j,k).gt.0.0_WP) then
                T=sc(ind_T)%SC(i,j,k)
                Y_CO2=min(max(sc(ind_CO2)%SC(i,j,k),0.0_WP),1.0_WP)
                Y_N2=1.0_WP-Y_CO2
                Wmix=1.0_WP/(Y_N2/W_N2+Y_CO2/W_CO2)
                sc(1)%rho(i,j,k)=Pthermo*Wmix/(Rcst*T)*(1.0_WP-lp%VF(i,j,k))
             else
                sc(1)%rho(i,j,k)=1.0_WP
             end if
          end do
       end do
    end do
    ! Update all other scalars rhos
    do i=2,nscalar
       sc(i)%rho=sc(1)%rho
    end do
  end subroutine get_sc_rho

  !> Compute viscosity based on Sutherland's law
  subroutine get_viscosity
    real(WP) :: T
    real(WP), parameter :: Tref=273.11_WP
    real(WP), parameter :: Sref=110.56_WP
    integer :: i,j,k
    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
          do i=fs%cfg%imino_,fs%cfg%imaxo_
             if (fs%cfg%VF(i,j,k).gt.0.0_WP) then
                T=sc(ind_T)%SC(i,j,k)
                fs%visc(i,j,k)=visc*(T/Tref)**1.5_WP*(Tref+Sref)/(T+Sref)
             else
                fs%visc(i,j,k)=1.0_WP
             end if
          end do
       end do
    end do
  end subroutine get_viscosity


  !> Compute thermal diffusivity based on Sutherland's law
  subroutine get_thermal_diffusivity
    real(WP) :: T
    real(WP), parameter :: Tref=273.11_WP
    real(WP), parameter :: Sref=110.56_WP
    integer :: i,j,k
    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
          do i=fs%cfg%imino_,fs%cfg%imaxo_
             if (fs%cfg%VF(i,j,k).gt.0.0_WP) then
                T=sc(ind_T)%SC(i,j,k)
                sc(ind_T)%diff(i,j,k)=diffusivity*(T/Tref)**1.5_WP*(Tref+Sref)/(T+Sref)
             else
                sc(ind_T)%diff(i,j,k)=0.0_WP
             end if
          end do
       end do
    end do
  end subroutine get_thermal_diffusivity


   ! !> Specialized subroutine that outputs useful statistics
   ! subroutine postproc_stat()
   !    use string,    only: str_medium
   !    use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
   !    use parallel,  only: MPI_REAL_WP
   !    implicit none
   !    integer :: iunit,ierr,i,j,k
   !    real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_
   !    character(len=str_medium) :: filename,timestamp
   !    ! Allocate vertical line storage
   !    allocate(Uavg (fs%cfg%imin:fs%cfg%imax)); Uavg =0.0_WP
   !    allocate(Uavg_(fs%cfg%imin:fs%cfg%imax)); Uavg_=0.0_WP
   !    allocate(vol_ (fs%cfg%imin:fs%cfg%imax)); vol_ =0.0_WP
   !    allocate(vol  (fs%cfg%imin:fs%cfg%imax)); vol  =0.0_WP
   !    ! Integrate all data over x and z
   !    do k=fs%cfg%kmin_,fs%cfg%kmax_
   !       do j=fs%cfg%jmin_,fs%cfg%jmax_
   !          do i=fs%cfg%imin_,fs%cfg%imax_
   !             vol_(i) = vol_(i)+fs%cfg%vol(i,j,k)
   !             Uavg_(i)=Uavg_(i)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)
   !          end do
   !       end do
   !    end do
   !    ! All-reduce the data
   !    call MPI_ALLREDUCE( vol_, vol,fs%cfg%nx,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%nx,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
   !    do i=fs%cfg%imin,fs%cfg%imax
   !       if (vol(i).gt.0.0_WP) then
   !          Uavg(i)=Uavg(i)/vol(i)
   !       else
   !          Uavg(i)=0.0_WP
   !       end if
   !    end do
   !    ! If root, print it out
   !    if (fs%cfg%amRoot) then
   !       filename='Uavg_'
   !       write(timestamp,'(es12.5)') time%t
   !       open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
   !       write(iunit,'(a12,3x,a12)') 'Height','Uavg'
   !       do i=fs%cfg%imin,fs%cfg%imax
   !          write(iunit,'(es12.5,3x,es12.5)') fs%cfg%xm(i),Uavg(i)
   !       end do
   !       close(iunit)
   !    end if
   !    ! Deallocate work arrays
   !    deallocate(Uavg,Uavg_,vol,vol_)
   ! end subroutine postproc_stat

  !> Function that localizes the left (x-) of the domain
  function left_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imin) isIn=.true.
  end function left_of_domain

  !> Function that localizes the right (x+) of the domain
  function right_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imax+1) isIn=.true.
  end function right_of_domain


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none

    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max time',time%tmax)
      call param_read('Max cfl number',time%cflmax)
      time%dt=time%dtmax
      time%itmax=2
    end block initialize_timetracker


    ! Initialize timers
    initialize_timers: block
      wt_total%time=0.0_WP; wt_total%percent=0.0_WP
      wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
      wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
      wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
      wt_sc%time=0.0_WP;    wt_sc%percent=0.0_WP
      wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
    end block initialize_timers


    ! Create a low Mach flow solver with bconds
    create_flow_solver: block
      use ils_class,     only: pcg_pfmg,pcg_amg
      use lowmach_class, only: dirichlet,clipped_neumann
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
      ! Define boundary conditions
      call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
      call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )

      ! ! Assign constant density
      ! call param_read('Density',rho); fs%rho=rho
      ! Assign constant viscosity
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Assign acceleration of gravity
      call param_read('Gravity',fs%gravity)
      ! Configure pressure solver
      call param_read('Pressure iteration',fs%psolv%maxit)
      call param_read('Pressure tolerance',fs%psolv%rcvg)
      ! Configure implicit velocity solver
      call param_read('Implicit iteration',fs%implicit%maxit)
      call param_read('Implicit tolerance',fs%implicit%rcvg)
      ! Setup the solver
      call fs%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)
    end block create_flow_solver

    ! Create scalar solvers for T and CO2
    create_scalar: block
      use ils_class,      only: gmres
      use vdscalar_class, only: bcond,quick,dirichlet,neumann
      integer :: ii
      ! Create scalar solvers
      sc(ind_T)  =vdscalar(cfg=cfg,scheme=quick,name='T')
      sc(ind_CO2)=vdscalar(cfg=cfg,scheme=quick,name='CO2')
      ! Assign constant diffusivity
      call param_read('Thermal diffusivity',diffusivity)
      sc(ind_T)%diff=diffusivity
      call param_read('CO2 diffusivity',diffusivity)
      sc(ind_CO2)%diff=diffusivity      
      do ii=1,nscalar
         ! Configure implicit scalar solver
         sc(ii)%implicit%maxit=fs%implicit%maxit; sc(ii)%implicit%rcvg=fs%implicit%rcvg
         ! Define boundary conditions
         call sc(ii)%add_bcond(name='inflow',type=dirichlet,locator=left_of_domain,dir='-x')
         call sc(ii)%add_bcond(name='outflow',type=neumann,locator=right_of_domain,dir='+x')
         ! Setup the solver
         call sc(ii)%setup(implicit_ils=gmres)
      end do
    end block create_scalar

    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(dRHOdt  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resSC   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcSClp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,nscalar))
      allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(rhof    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays

    ! Initialize our LPT solver
    initialize_lpt: block
      use random, only: random_uniform
      use mathtools, only: Pi
      real(WP) :: dp,Hbed,VFavg,Tp,Lpart,Lparty,Lpartz,Volp
      integer :: i,ix,iy,iz,np,npx,npy,npz
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Set scalar information
      call lp%scalar_init(sc=sc)
      ! Get drag model from the input
      call param_read('Drag model',lp%drag_model,default='Tenneti')
      ! Get adsorption model from the input
      call param_read('Adsorption model',lp%ads_model,default='NONE')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from the input
      call param_read('Particle diameter',dp)
      ! Set filter scale to 3.5*dx
      lp%filter_width=3.5_WP*cfg%min_meshsize

      ! Root process initializes particles uniformly
      call param_read('Bed height',Hbed)
      call param_read('Particle volume fraction',VFavg)
      call param_read('Particle temperature',Tp,default=298.15_WP)
      if (lp%cfg%amRoot) then
         ! Particle volume
         Volp = Pi/6.0_WP*dp**3
         ! Get number of particles
         Lpart = (Volp/VFavg)**(1.0_WP/3.0_WP)
         npx=int(Hbed/Lpart)
         npy = int(cfg%yL/Lpart)
         Lparty = cfg%yL/real(npy,WP)
         npz = int(cfg%zL/Lpart)
         Lpartz = cfg%zL/real(npz,WP)
         np = npx*npy*npz
         call lp%resize(np)
         ! Distribute particles
         do i=1,np
            ! Give position
            ix = (i-1)/(npy*npz)
            iy = (i-1-npy*npz*ix)/npz
            iz = i-1-npy*npz*ix-npz*iy
            lp%p(i)%pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart
            lp%p(i)%pos(2) = lp%cfg%y(lp%cfg%jmin)+(real(iy,WP)+0.5_WP)*Lparty
            lp%p(i)%pos(3) = lp%cfg%z(lp%cfg%kmin)+(real(iz,WP)+0.5_WP)*Lpartz

            ! Give id
            !lp%p(i)%id=-1
            lp%p(i)%id=int(i,8)
            ! Set the diameter
            lp%p(i)%d=dp
            ! Set the temperature
            lp%p(i)%T=Tp
            ! Give zero mass of carbamate
            lp%p(i)%Mc=0.0_WP
            ! Give zero velocity
            lp%p(i)%vel=0.0_WP
            ! Give zero collision force
            lp%p(i)%Acol=0.0_WP
            lp%p(i)%Tcol=0.0_WP
            ! Give zero dt
            lp%p(i)%dt=0.0_WP
            ! Locate the particle on the mesh
            lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
            ! Activate the particle
            lp%p(i)%flag=0
         end do
      end if
      call lp%sync()

      ! Get initial particle volume fraction
      call lp%update_VF()
      ! Set collision timescale
      call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',lp%e_n)
      call param_read('Wall restitution',lp%e_w,default=lp%e_n)
      call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
      ! Set gravity
      call param_read('Gravity',lp%gravity)
      if (lp%cfg%amRoot) then
         print*,"===== Particle Setup Description ====="
         print*,'Number of particles', np
         print*,'Mean volume fraction',VFavg
      end if
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=2,nvec=3,name='lpt')
      pmesh%varname(1)='diameter'
      pmesh%varname(2)='Mc'
      pmesh%vecname(1)='velocity'
      pmesh%vecname(2)='ang_vel'
      pmesh%vecname(3)='Fcol'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=lp%p(i)%d
         pmesh%var(2,i)=lp%p(i)%Mc
         pmesh%vec(:,1,i)=lp%p(i)%vel
         pmesh%vec(:,2,i)=lp%p(i)%angVel
         pmesh%vec(:,3,i)=lp%p(i)%Acol
      end do
    end block create_pmesh

    ! Initialize our scalar fields
    initialize_scalar: block
      use vdscalar_class, only: bcond
      type(bcond), pointer :: mybc
      real(WP) :: CO2i,Ti
      integer :: i,j,k,n,ii
      ! Read in the intial values
      call param_read('Initial T',Ti)
      call param_read('Initial CO2',CO2i)
      call param_read('Pressure',Pthermo)
      call param_read('Inlet CO2',SC_inlet(ind_CO2))
      call param_read('Inlet T',SC_inlet(ind_T))
      ! Assign values
      sc(ind_T)%SC=Ti
      sc(ind_CO2)%SC=CO2i
      do ii=1,nscalar
         ! Initialize the scalars at the inlet
         call sc(ii)%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc(ii)%SC(i,j,k)=SC_inlet(ii)
         end do
         ! Apply all other boundary conditions
         call sc(ii)%apply_bcond(time%t,time%dt)
         ! Get density from pressure, temperature, and mass fraction
         call get_sc_rho()
         call sc(ii)%rho_multiply()
      end do
      ! Update transport variables
      call get_viscosity
      call get_thermal_diffusivity
      ! Update fluid density
      rhof=sc(1)%rho/(1.0_WP-lp%VF)
    end block initialize_scalar

    ! Initialize our velocity field
    initialize_velocity: block
      use lowmach_class, only: bcond
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      ! Zero initial field
      fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
      ! Set inflow velocity/momentum
      call param_read('Inlet velocity',inlet_velocity)
      call fs%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         fs%rhoU(i,j,k)=sc(1)%rho(i,j,k)*inlet_velocity
         fs%U(i,j,k)   =inlet_velocity
      end do
      ! Set density
      fs%rho=sc(1)%rho
      ! Form momentum
      call fs%rho_multiply
      ! Apply all other boundary conditions
      call fs%apply_bcond(time%t,time%dt)
      call fs%interp_vel(Ui,Vi,Wi)
      call fs%get_div(drhodt=dRHOdt)
      ! Compute MFR through all boundary conditions
      call fs%get_mfr()
    end block initialize_velocity


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='fluidized_bed')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('epsp',lp%VF)
      call ens_out%add_scalar('density',rhof)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('temperature',sc(ind_T)%SC)
      call ens_out%add_scalar('CO2',sc(ind_CO2)%SC)
      call ens_out%add_scalar('viscosity',fs%visc)
      call ens_out%add_scalar('diffusivity',sc(ind_T)%diff)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create monitor file
    create_monitor: block
      use string, only: str_medium
      integer :: ii
      character(len=str_medium) :: str
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call lp%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      call lp%get_max()
      do ii=1,nscalar
         call sc(ii)%get_max()
         call sc(ii)%get_int()
      end do
      ! Create simulation monitor
      mfile=monitor(fs%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(fs%Umax,'Umax')
      call mfile%add_column(fs%Vmax,'Vmax')
      call mfile%add_column(fs%Wmax,'Wmax')
      call mfile%add_column(fs%Pmax,'Pmax')
      call mfile%add_column(fs%divmax,'Maximum divergence')
      call mfile%add_column(fs%psolv%it,'Pressure iteration')
      call mfile%add_column(fs%psolv%rerr,'Pressure error')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(fs%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
      call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
      call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
      call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
      call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
      call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
      call cflfile%add_column(lp%CFL_col,'Collision CFL')
      call cflfile%write()
      ! Create scalar monitor
      scfile=monitor(fs%cfg%amRoot,'scalar')
      call scfile%add_column(time%n,'Timestep number')
      call scfile%add_column(time%t,'Time')
      call scfile%add_column(sc(1)%rhomin,'RHOmin')
      call scfile%add_column(sc(1)%rhomax,'RHOmax')
      call scfile%add_column(sc(1)%rhoint,'RHOint')
      do ii=1,nscalar
         str=trim(sc(ii)%name)//'min'
         call scfile%add_column(sc(ii)%SCmin,trim(str))
         str=trim(sc(ii)%name)//'max'
         call scfile%add_column(sc(ii)%SCmax,trim(str))
         str=trim(sc(ii)%name)//'int'
         call scfile%add_column(sc(ii)%SCint,trim(str))
      end do
      call adsfile%write()
      ! Create adsorption monitor
      adsfile=monitor(fs%cfg%amRoot,'adsorption')
      call adsfile%add_column(time%n,'Timestep number')
      call adsfile%add_column(time%t,'Time')
      call adsfile%add_column(Pthermo,'Pthermo')
      call adsfile%add_column(sc(ind_CO2)%SCmax,'YCO2max')
      call adsfile%add_column(sc(ind_CO2)%SCmin,'YCO2min')
      call adsfile%add_column(lp%Mcmin,'Particle CO2min')
      call adsfile%add_column(lp%Mcmax,'Particle CO2max')
      call adsfile%add_column(lp%Mcmean,'Particle CO2mean')
      call adsfile%write()
      ! Create LPT monitor
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%VFmean,'VFp mean')
      call lptfile%add_column(lp%VFmax,'VFp max')
      call lptfile%add_column(lp%np,'Particle number')
      call lptfile%add_column(lp%Umin,'Particle Umin')
      call lptfile%add_column(lp%Umax,'Particle Umax')
      call lptfile%add_column(lp%Vmin,'Particle Vmin')
      call lptfile%add_column(lp%Vmax,'Particle Vmax')
      call lptfile%add_column(lp%Wmin,'Particle Wmin')
      call lptfile%add_column(lp%Wmax,'Particle Wmax')
      call lptfile%add_column(lp%dmin,'Particle dmin')
      call lptfile%add_column(lp%dmax,'Particle dmax')
      call lptfile%write()
      ! Create timing monitor
      tfile=monitor(amroot=fs%cfg%amRoot,name='timing')
      call tfile%add_column(time%n,'Timestep number')
      call tfile%add_column(time%t,'Time')
      call tfile%add_column(wt_total%time,'Total [s]')
      call tfile%add_column(wt_vel%time,'Velocity [s]')
      call tfile%add_column(wt_vel%percent,'Velocity [%]')
      call tfile%add_column(wt_pres%time,'Pressure [s]')
      call tfile%add_column(wt_pres%percent,'Pressure [%]')
      call tfile%add_column(wt_lpt%time,'LPT [s]')
      call tfile%add_column(wt_lpt%percent,'LPT [%]')
      call tfile%add_column(wt_sc%time,'Scalar [s]')
      call tfile%add_column(wt_sc%percent,'Scalar [%]')
      call tfile%add_column(wt_rest%time,'Rest [s]')
      call tfile%add_column(wt_rest%percent,'Rest [%]')
      call tfile%write()
    end block create_monitor


    ! ! Create a specialized post-processing file
    ! create_postproc: block
    !   ! Create event for data postprocessing
    !   ppevt=event(time=time,name='Postproc output')
    !   call param_read('Postproc output period',ppevt%tper)
    !   ! Perform the output
    !   if (ppevt%occurs()) call postproc_stat()
    ! end block create_postproc

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use mathtools, only: twoPi
    use parallel, only: parallel_time
    implicit none
    integer :: ii

    ! Perform time integration
    do while (.not.time%done())
       
       ! Initial wallclock time
       wt_total%time_in=parallel_time()

       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call lp%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old scalar
       do ii=1,nscalar
          sc(ii)%rhoold=sc(ii)%rho
          sc(ii)%SCold =sc(ii)%SC
       end do

       ! Remember old velocity, and momentum
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Get fluid stress
       wt_lpt%time_in=parallel_time()
       call fs%get_div_stress(resU,resV,resW)

       ! Predict new density by linear extrapolation
       ! Is this needed?
       ! sc(1)%rho=sc(1)%rho+time%dt*dRHOdt
       ! do ii=2,nscalar
       !    sc(i)%rho=sc(1)%rho
       ! end do

       ! Collide and advance particles
       call lp%collide(dt=time%dtmid)
       rhof=sc(1)%rho/(1.0_WP-lp%VF)
       call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=rhof,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
            T=sc(ind_T)%SC,YCO2=sc(ind_CO2)%SC,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcSC=srcSClp)
       
       wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)
 
          ! ============= SCALAR SOLVER =======================
          wt_sc%time_in=parallel_time()
          do ii=1,2
             ! Build mid-time scalar
             sc(ii)%SC=0.5_WP*(sc(ii)%SC+sc(ii)%SCold)

             ! Explicit calculation of drhoSC/dt from scalar equation
             call sc(ii)%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
             ! call fs%cfg%integrate(resSC,integral=mean); mean=mean/fs%cfg%vol_total
             ! print *,ii,'resSC',mean

             ! Assemble explicit residual
             resSC=time%dt*resSC-(2.0_WP*sc(ii)%rho*sc(ii)%SC-(sc(ii)%rho+sc(ii)%rhoold)*sc(ii)%SCold)

             ! Add in mass source from particle
             resSC = resSC + srcSClp(:,:,:,ii)
             ! call fs%cfg%integrate(resSC,integral=mean); mean=mean/fs%cfg%vol_total
             ! print *,ii,'resSC',mean

             ! Form implicit residual
             !call sc(ii)%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)

             ! Apply this residual
             sc(ii)%SC=2.0_WP*sc(ii)%SC-sc(ii)%SCold+resSC!/sc(ii)%rho

             ! Apply other boundary conditions on the resulting field
             call sc(ii)%apply_bcond(time%t,time%dt)

             ! Apply scalar boundary conditions
             scalar_bcond: block
               use vdscalar_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ni
               call sc(ii)%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc(ii)%SC(i,j,k)=SC_inlet(ii)
               end do
             end block scalar_bcond

             call sc(ii)%rho_multiply()
          end do
          
          wt_sc%time=wt_sc%time+parallel_time()-wt_sc%time_in
          ! ===================================================

          ! Update dependent variables
          resSC = sc(1)%rho
          call get_sc_rho
          ! do ii=1,nscalar
          !    sc(ii)%sc=sc(ii)%sc*resSC/sc(ii)%rho
          ! end do

          call get_viscosity
          call get_thermal_diffusivity
          
          ! Build mid-density
          fs%rho=0.5_WP*(sc(1)%rho+sc(1)%rhoold)
          wt_vel%time_in=parallel_time()

          ! Compute rate-of-change of density accounting for particles
          ! CHECK THIS dRHOdt
          dRHOdt=(sc(1)%RHO-sc(1)%RHOold)/time%dtmid-srcSClp(:,:,:,ind_CO2)/time%dtmid

          ! Build mid-time velocity and momentum
          fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
          fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
          fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

          ! Add momentum source terms
          call fs%addsrc_gravity(resU,resV,resW)

          ! Assemble explicit residual
          resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
          resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
          resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

          ! Add momentum source term from lpt
          add_lpt_src: block
            integer :: i,j,k
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                     resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                     resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                  end do
               end do
            end do
          end block add_lpt_src
    
          ! Form implicit residuals
          call fs%solve_implicit(time%dtmid,resU,resV,resW)
    
          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW
    
          ! Apply other boundary conditions and update momentum
          call fs%apply_bcond(time%tmid,time%dtmid)
          call fs%rho_multiply()
          call fs%apply_bcond(time%tmid,time%dtmid)

          ! Reset Dirichlet BCs
          dirichlet_velocity: block
            use lowmach_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: n,i,j,k
            call fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%rhoU(i,j,k)=fs%rho(i,j,k)*inlet_velocity
               fs%U(i,j,k)   =inlet_velocity
            end do
          end block dirichlet_velocity

          wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

          ! Solve Poisson equation
          wt_pres%time_in=parallel_time()
          call fs%correct_mfr(drhodt=dRHOdt)
          call fs%get_div(drhodt=dRHOdt)
          fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
          fs%psolv%sol=0.0_WP
          call fs%psolv%solve()
          call fs%shift_p(fs%psolv%sol)

          ! Correct momentum and rebuild velocity
          call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
          fs%P=fs%P+fs%psolv%sol
          fs%rhoU=fs%rhoU-time%dtmid*resU
          fs%rhoV=fs%rhoV-time%dtmid*resV
          fs%rhoW=fs%rhoW-time%dtmid*resW
          call fs%rho_divide
          wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       wt_vel%time_in=parallel_time()
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=dRHOdt)
       wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=lp%p(i)%d
               pmesh%var(2,i)=lp%p(i)%Mc
               pmesh%vec(:,1,i)=lp%p(i)%vel
               pmesh%vec(:,2,i)=lp%p(i)%angVel
               pmesh%vec(:,3,i)=lp%p(i)%Acol
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call fs%get_max()
       call lp%get_max()
       do ii=1,nscalar
          call sc(ii)%get_max()
          call sc(ii)%get_int()
       end do
       call mfile%write()
       call cflfile%write()
       call lptfile%write()
       call adsfile%write()
       call scfile%write()

       ! ! Specialized post-processing
       ! if (ppevt%occurs()) call postproc_stat()

       ! Monitor timing
       wt_total%time=parallel_time()-wt_total%time_in
       wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
       wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
       wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
       wt_sc%percent=wt_sc%time/wt_total%time*100.0_WP
       wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_lpt%time-wt_sc%time
       wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
       call tfile%write()
       wt_total%time=0.0_WP; wt_total%percent=0.0_WP
       wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
       wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
       wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
       wt_sc%time=0.0_WP;    wt_sc%percent=0.0_WP
       wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP

    end do

  end subroutine simulation_run


  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    ! monitor
    ! ensight
    ! bcond
    ! timetracker

    ! Deallocate work arrays
    deallocate(resU,resV,resW,srcUlp,srcVlp,srcWlp,Ui,Vi,Wi,dRHOdt,rhof)

  end subroutine simulation_final

end module simulation



