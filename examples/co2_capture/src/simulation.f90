!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg,D
  use lpt_class,         only: lpt
  use fft2d_class,       only: fft2d
  use ddadi_class,       only: ddadi
  use lowmach_class,     only: lowmach
  use vdscalar_class,    only: vdscalar
  use gp_class,          only: gpibm
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get a scalar solver, an LPT solver, a lowmach solver, a ghost point IBM solver
  !> and corresponding time tracker, plus a couple of linear solvers
  type(vdscalar), dimension(2), public :: sc
  type(fft2d),       public :: ps
  type(ddadi),       public :: vs
  type(lowmach),     public :: fs
  type(ddadi), dimension(:), allocatable, public :: ss
  type(lpt),         public :: lp
  type(gpibm),       public :: gp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile,scfile,lptfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rhof
  real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp,ghost
  real(WP), dimension(:,:,:,:), allocatable :: srcSClp
  real(WP) :: visc,diffusivity,rhoUin

  !> Scalar indices
  integer, parameter :: ind_T  =1
  integer, parameter :: ind_CO2=2
  integer :: nscalar=2
  real(WP), dimension(2) :: SC0,SCin,SCmin,SCmax

  !> Thermo-chemistry parameters
  real(WP), parameter :: Rcst=8.314_WP         !< J/(mol.K)
  real(WP), parameter :: W_N2  = 28.0135e-3_WP !< kg/mol
  real(WP), parameter :: W_CO2 = 44.0095e-3_WP !< kg/mol
  real(WP), parameter :: W_H2O = 18.0153e-3_WP !< kg/mol
  real(WP) :: Pthermo,fCp

  ! Ghost point IBM
  integer :: gp_no
  integer, dimension(:,:),allocatable :: gp_map

contains


  !> Define here our equation of state - rho(P,T,Yk)
  !> This just updates sc%rho
  subroutine get_sc_rho()
    real(WP) :: T,Y_CO2,Y_N2,Wmix
    integer :: i,j,k
    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
          do i=fs%cfg%imino_,fs%cfg%imaxo_
             T=sc(ind_T)%SC(i,j,k)
             Y_CO2=min(max(sc(ind_CO2)%SC(i,j,k),0.0_WP),1.0_WP)
             Y_N2=1.0_WP-Y_CO2
             Wmix=1.0_WP/(Y_N2/W_N2+Y_CO2/W_CO2)
             sc(1)%rho(i,j,k)=Pthermo*Wmix/(Rcst*T)*(1.0_WP-lp%VF(i,j,k))
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
             T=sc(ind_T)%SC(i,j,k)
             fs%visc(i,j,k)=visc*(T/Tref)**1.5_WP*(Tref+Sref)/(T+Sref)
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
             T=sc(ind_T)%SC(i,j,k)
             sc(ind_T)%diff(i,j,k)=diffusivity*(T/Tref)**1.5_WP*(Tref+Sref)/(T+Sref)
          end do
       end do
    end do
  end subroutine get_thermal_diffusivity


  !> Function that localizes the left (x-) of the domain
  function left_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imin.and.j.ge.pg%jmin.and.j.le.pg%jmax.and.k.ge.pg%kmin.and.k.le.pg%kmax) isIn=.true.
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


  !> Function that localizes the left (x-) of the domain for scalars
  function left_of_domain_sc(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    real(WP) :: radius
    logical :: isIn
    isIn=.false.
    !if (i.eq.pg%imin-1) isIn=.true.
    radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
    if (radius.le.0.5_WP*D.and.i.eq.pg%imin) isIn=.true.
  end function left_of_domain_sc


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

    ! Create a low Mach flow solver with bconds
    create_flow_solver: block
      use lowmach_class, only: dirichlet,clipped_neumann
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
      ! Define boundary conditions
      call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
      call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
      ! Assign constant viscosity
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Assign heat capacity
      call param_read('Heat capacity',fCp)
      ! Configure pressure solver
      ps=fft2d(cfg=cfg,name='Pressure',nst=7)
      ! Configure implicit velocity solver
      vs=ddadi(cfg=cfg,name='Velocity',nst=7)
      ! Setup the solver
      call fs%setup(pressure_solver=ps,implicit_solver=vs)
    end block create_flow_solver


    ! Create scalar solvers for T and CO2
    create_scalar: block
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
      ! Configure implicit scalar solver
      allocate(ss(nscalar))
      do ii=1,nscalar
         ! Define boundary conditions
         call sc(ii)%add_bcond(name='inflow',type=dirichlet,locator=left_of_domain_sc)
         call sc(ii)%add_bcond(name='outflow',type=neumann,locator=right_of_domain,dir='+x')
         ! Setup the solver
         ss(ii)=ddadi(cfg=cfg,name='Scalar',nst=13)
         call sc(ii)%setup(implicit_solver=ss(ii))
      end do
    end block create_scalar


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resSC   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(srcSClp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:nscalar))
      allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(rhof    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(ghost   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize our LPT solver
    initialize_lpt: block
      use random, only: random_uniform
      use mathtools, only: Pi
      real(WP) :: dp,VFavg,Tp,Lpart,Lparty,Lpartz,Volp
      real(WP), dimension(3) :: pos
      integer :: i,ix,iy,iz,np_,np,npx,npy,npz
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from the input
      call param_read('Particle diameter',dp)
      ! Set particle temperature
      call param_read('Particle temperature',Tp,default=298.15_WP)
      ! Get particle heat capacity from the input
      !call param_read('Particle heat capacity',lp%pCp)
      ! Ger filter width
      call param_read('Filter width',lp%filter_width)
      ! Set volume fraction
      call param_read('Particle volume fraction',VFavg)
      ! Root process initializes particles uniformly
      if (lp%cfg%amRoot) then
         ! Particle volume
         Volp = Pi/6.0_WP*dp**3
         ! Get number of particles
         Lpart = (Volp/VFavg)**(1.0_WP/3.0_WP)
         npx=int(cfg%xL/Lpart)
         npy = int(cfg%yL/Lpart)
         Lparty = cfg%yL/real(npy,WP)
         npz = int(cfg%zL/Lpart)
         Lpartz = cfg%zL/real(npz,WP)
         np_ = npx*npy*npz
         np=0
         ! Remove particles outside of the cylinder
         do i=1,np_
            ! Give position
            ix = (i-1)/(npy*npz)
            iy = (i-1-npy*npz*ix)/npz
            iz = i-1-npy*npz*ix-npz*iy
            pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart
            pos(2) = lp%cfg%y(lp%cfg%jmin)+(real(iy,WP)+0.5_WP)*Lparty
            pos(3) = lp%cfg%z(lp%cfg%kmin)+(real(iz,WP)+0.5_WP)*Lpartz
            if (sqrt(sum(pos(2:3)**2)).lt.0.5_WP*D-0.5_WP*dp) np=np+1
         end do
         call lp%resize(np)
         ! Distribute particles
         np=0
         do i=1,np_
            ! Give position
            ix = (i-1)/(npy*npz)
            iy = (i-1-npy*npz*ix)/npz
            iz = i-1-npy*npz*ix-npz*iy
            pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart
            pos(2) = lp%cfg%y(lp%cfg%jmin)+(real(iy,WP)+0.5_WP)*Lparty
            pos(3) = lp%cfg%z(lp%cfg%kmin)+(real(iz,WP)+0.5_WP)*Lpartz
            if (sqrt(sum(pos(2:3)**2)).ge.0.5_WP*D-0.5_WP*dp) cycle
            np=np+1
            lp%p(np)%pos(1)=pos(1)
            lp%p(np)%pos(2)=pos(2)
            lp%p(np)%pos(3)=pos(3)
            
            ! Give id (fixed)
            lp%p(np)%id=-1
            ! Set the diameter
            lp%p(np)%d=dp
            ! Set the temperature
            !lp%p(np)%T=Tp
            ! Give zero velocity
            lp%p(np)%vel=0.0_WP
            ! Give zero collision force
            lp%p(np)%Acol=0.0_WP
            lp%p(np)%Tcol=0.0_WP
            ! Give zero dt
            lp%p(np)%dt=0.0_WP
            ! Locate the particle on the mesh
            lp%p(np)%ind=lp%cfg%get_ijk_global(lp%p(np)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
            ! Activate the particle
            lp%p(np)%flag=0
         end do
      end if
      call lp%sync()

      ! Get initial particle volume fraction
      call lp%update_VF()
      if (lp%cfg%amRoot) then
         print*,"===== Particle Setup Description ====="
         print*,'Number of particles', np
         print*,'Mean volume fraction',VFavg
      end if
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=1,nvec=0,name='lpt')
      pmesh%varname(1)='diameter'
      !pmesh%varname(2)='temperature'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=lp%p(i)%d
         !pmesh%var(2,i)=lp%p(i)%Tp
      end do
    end block create_pmesh


    ! Initialize our scalar fields
    initialize_scalar: block
      use vdscalar_class, only: bcond
      type(bcond), pointer :: mybc
      integer :: i,j,k,n,ii
      ! Read in the intial values
      call param_read('Pressure',Pthermo)
      call param_read('Initial T',SC0(ind_T))
      call param_read('Initial CO2',SC0(ind_CO2))
      call param_read('Inlet T',SCin(ind_T))
      call param_read('Inlet CO2',SCin(ind_CO2))
      call param_read('Min T',SCmin(ind_T),default=-huge(1.0_WP))
      call param_read('Max T',SCmax(ind_T),default=-huge(1.0_WP))
      call param_read('Min CO2',SCmin(ind_CO2),default=-huge(1.0_WP))
      call param_read('Max CO2',SCmax(ind_CO2),default=huge(1.0_WP))
      ! Assign values
      sc(ind_T)%SC=SC0(ind_T)
      sc(ind_CO2)%SC=SC0(ind_CO2)
      do ii=1,nscalar
         ! Initialize the scalars at the inlet
         call sc(ii)%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc(ii)%SC(i,j,k)=SCin(ii)
         end do
         ! Apply all other boundary conditions
         call sc(ii)%apply_bcond(time%t,time%dt)
      end do
      ! Get density from pressure, temperature, and mass fraction
      call get_sc_rho()
      do ii=1,nscalar
         call sc(ii)%rho_multiply()
      end do
      ! Update transport variables
      call get_viscosity
      call get_thermal_diffusivity
    end block initialize_scalar


    ! Initialize our velocity field
    initialize_velocity: block
      integer :: n,i,j,k
      real(WP) :: T,Y_CO2,Y_N2,Wmix
      ! Set density
      fs%rho=sc(1)%rho
      rhof=fs%rho/(1.0_WP-lp%VF)
      ! Read inlet velocity
      call param_read('Inlet velocity',rhoUin)
      ! Set uniform momentum and velocity
      T=SCin(ind_T); Y_CO2=SCin(ind_CO2); Y_N2=1.0_WP-Y_CO2
      Wmix=1.0_WP/(Y_N2/W_N2+Y_CO2/W_CO2)
      rhoUin=Pthermo*Wmix/(Rcst*T)*rhoUin
      fs%rhoU=rhoUin; fs%rhoV=0.0_WP; fs%rhoW=0.0_WP
      call fs%rho_divide()
      call fs%interp_vel(Ui,Vi,Wi)
      resSC=0.0_WP
      call fs%get_div(drhodt=resSC)
      ! Compute MFR through all boundary conditions
      call fs%get_mfr()
    end block initialize_velocity


    ! Initialize the ghost points
    create_gp: block
      gp=gpibm(cfg=cfg,no=2)
      call gp%update()
      ghost=real(gp%label,WP)
    end block create_gp


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='pipe')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('levelset',cfg%Gib)
      call ens_out%add_scalar('density',rhof)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('temperature',sc(ind_T)%SC)
      call ens_out%add_scalar('CO2',sc(ind_CO2)%SC)
      call ens_out%add_scalar('viscosity',fs%visc)
      call ens_out%add_scalar('diffusivity',sc(ind_T)%diff)
      call ens_out%add_scalar('epsp',lp%VF)
      call ens_out%add_scalar('ghostpoints',ghost)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create a monitor file
    create_monitor: block
      use string, only: str_medium
      integer :: ii
      character(len=str_medium) :: str
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
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
      call scfile%write()
      ! Create LPT monitor
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%VFmean,'VFp mean')
      call lptfile%add_column(lp%VFmax,'VFp max')
      call lptfile%add_column(lp%np,'Particle number')
      !call lptfile%add_column(lp%drag_time,'Drag time')
      call lptfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    implicit none
    integer :: ii

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Particle solver
       lpt: block
         ! Get fluid stress and compute source terms
         call fs%get_div_stress(resU,resV,resW)
         call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=rhof,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
              srcU=srcUlp,srcV=srcVlp,srcW=srcWlp)
       end block lpt

       ! Remember old scalar
       do ii=1,nscalar
          sc(ii)%rhoold=sc(ii)%rho
          sc(ii)%SCold =sc(ii)%SC
       end do

       ! Remember old velocity and momentum
       fs%RHOold=fs%RHO
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          ! ============= SCALAR SOLVER =======================
          do ii=1,nscalar

             ! Build mid-time scalar
             sc(ii)%SC=0.5_WP*(sc(ii)%SC+sc(ii)%SCold)

             ! Explicit calculation of drhoSC/dt from scalar equation
             call sc(ii)%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

             ! Assemble explicit residual
             resSC=time%dt*resSC-(2.0_WP*sc(ii)%rho*sc(ii)%SC-(sc(ii)%rho+sc(ii)%rhoold)*sc(ii)%SCold)

             ! Form implicit residual
             call sc(ii)%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)

             ! Apply this residual
             sc(ii)%SC=2.0_WP*sc(ii)%SC-sc(ii)%SCold+resSC
             
             ! Clip
             where (sc(ii)%SC.lt.SCmin(ii)) sc(ii)%SC=SCmin(ii)
             where (sc(ii)%SC.gt.SCmax(ii)) sc(ii)%SC=SCmax(ii)

             ! Apply IB forcing to enforce Neumann at the pipe walls
             ibforcing_sc: block
               use gp_class, only: neumann
               where (cfg%Gib.lt.0.0_WP) sc(ii)%sc=SC0(ii)
               call gp%apply_bcond(type=neumann,BP=0.0_WP,A=sc(ii)%sc)
             end block ibforcing_sc

             ! Apply other boundary conditions on the resulting field
             call sc(ii)%apply_bcond(time%t,time%dt)

             ! Apply scalar boundary conditions
             scalar_bcond: block
               use vdscalar_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call sc(ii)%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc(ii)%SC(i,j,k)=SCin(ii)
               end do
             end block scalar_bcond
             call sc(ii)%rho_multiply()
          end do
          ! ===================================================

          ! ============ UPDATE PROPERTIES ====================
          ! Update dependent variables
          call get_sc_rho()
          call get_viscosity
          call get_thermal_diffusivity
          ! ===================================================


          ! ============ VELOCITY SOLVER ======================

          ! Build n+1 density
          fs%rho=0.5_WP*(sc(1)%rho+sc(1)%rhoold)

          ! Build mid-time velocity and momentum
          fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
          fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
          fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

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
                     !resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                     !resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcVlp(i,j-1:j,k))
                     !resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcWlp(i,j,k-1:k))
                  end do
               end do
            end do
          end block add_lpt_src

          ! Form implicit residuals
          call fs%solve_implicit(time%dt,resU,resV,resW)

          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW

          ! Apply IB forcing to enforce no-slip at the pipe walls
          ibforcing: block
            integer :: i,j,k
            real(WP) :: VFx,VFy,VFz
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                     VFy=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                     VFz=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                     fs%U(i,j,k)=VFx*fs%U(i,j,k)
                     fs%V(i,j,k)=VFy*fs%V(i,j,k)
                     fs%W(i,j,k)=VFz*fs%W(i,j,k)
                  end do
               end do
            end do
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
            call fs%cfg%sync(fs%W)
          end block ibforcing

          ! Apply other boundary conditions and update momentum
          call fs%apply_bcond(time%tmid,time%dtmid)
          call fs%rho_multiply()
          call fs%apply_bcond(time%tmid,time%dtmid)

          ! Reset Dirichlet BCs
          dirichlet_velocity: block
            use lowmach_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: n,i,j,k
            real(WP) :: VFx
            call fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
               fs%rhoU(i,j,k)=rhoUin*VFx
               fs%rhoV(i,j,k)=0.0_WP
               fs%rhoW(i,j,k)=0.0_WP
            end do
            call fs%rho_divide()
          end block dirichlet_velocity

          ! Solve Poisson equation
          resSC=(sc(1)%rho-sc(1)%rhoold)/time%dt
          call fs%correct_mfr(drhodt=resSC)
          call fs%get_div(drhodt=resSC)
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
          ! ===================================================

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute molecular density, interpolated velocity and divergence
       rhof=sc(1)%rho/(1.0_WP-lp%VF)
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=resSC)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=lp%p(i)%d
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
       call scfile%write()
       call lptfile%write()

    end do

  end subroutine simulation_run


  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    ! monitor
    ! ensight
    ! timetracker

    ! Deallocate work arrays
    deallocate(resU,resV,resW,resSC,srcUlp,srcVlp,srcWlp,srcSClp,Ui,Vi,Wi,rhof,gp_map)

  end subroutine simulation_final

end module simulation
