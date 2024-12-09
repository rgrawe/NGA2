!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use lpt_class,         only: lpt
  use fft2d_class,       only: fft2d
  use ddadi_class,       only: ddadi
  use lowmach_class,     only: lowmach
  use vdscalar_class,    only: vdscalar
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get a scalar solver, an LPT solver, a lowmach solver, and corresponding time tracker, plus a couple of linear solvers
  type(vdscalar),    public :: sc
  type(fft2d),       public :: ps
  type(ddadi),       public :: vs
  type(ddadi),       public :: ss
  type(lowmach),     public :: fs
  type(lpt),         public :: lp
  type(timetracker), public :: time
  
  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt
  
  !> Simulation monitor file
  type(monitor) :: mfile,cflfile,lptfile,scfile

  !> Position to start
  real(WP) :: x0
  logical  :: lagrange

  public :: simulation_init,simulation_run,simulation_final
  
  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rhof,Tf
  real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp,srcSClp
  real(WP), dimension(:,:,:), allocatable :: tmp
  real(WP) :: rhoUin,rho,dp,VFavg,Tp,SCin,fCp,diff_mult

contains

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

  !> Function that localizes the left (x-) of the domain for scalars
  function left_of_domain_sc(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    !if (i.eq.pg%imin-1) isIn=.true.
    if (pg%xm(i).le.x0) isIn=.true.
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
      real(WP) :: visc
      call param_read('Starting position',x0,default=0.0_WP)
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
      ! Define boundary conditions
      call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
      call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
      ! Assign constant viscosity
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Configure pressure solver
      ps=fft2d(cfg=cfg,name='Pressure',nst=7)
      ! Configure implicit velocity solver
      vs=ddadi(cfg=cfg,name='Velocity',nst=7)
      ! Setup the solver
      call fs%setup(pressure_solver=ps,implicit_solver=vs)
    end block create_flow_solver
    

    ! Create scalar solvers for T
    create_scalar: block
      use vdscalar_class, only: bcond,quick,dirichlet,neumann
      real(WP) :: diff
      ! Create scalar solver
      sc=vdscalar(cfg=cfg,scheme=quick,name='T')
      ! Assign constant diffusivity
      call param_read('Thermal diffusivity',diff)
      sc%diff=diff
      ! Assign heat capacity
      call param_read('Heat capacity',fCp)
      ! Define boundary conditions
      call sc%add_bcond(name='inflow',type=dirichlet,locator=left_of_domain_sc)
      call sc%add_bcond(name='outflow',type=neumann,locator=right_of_domain,dir='+x')
      ! Setup the solver
      ss=ddadi(cfg=cfg,name='Scalar',nst=13)
      call sc%setup(implicit_solver=ss)
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
      allocate(tmp     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp=0.0_WP
      allocate(srcSClp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Tf      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(rhof    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays

    ! Initialize our LPT solver
    initialize_lpt: block
      use random, only: random_uniform
      use mathtools, only: Pi
      real(WP) :: Lpart,Lparty,Lpartz,Volp,Vbox,Hbed
      integer :: i,j,k,ix,iy,iz,np,npx,npy,npz,nn,ii,jj,kk,ip,jp,kp
      logical :: lattice
      integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
      integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
      logical :: overlap
      
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from the input
      call param_read('Particle diameter',dp)
      ! Set particle temperature
      call param_read('Particle temperature',Tp,default=298.15_WP)
      ! Get particle heat capacity from the input
      call param_read('Particle heat capacity',lp%pCp)
      ! Set filter width to 3.5*dx
      call param_read('Filter width',lp%filter_width)
      ! Set height and volume fraction of the particle bed
      call param_read('Box height',Hbed)
      call param_read('Particle volume fraction',VFavg)
      !call param_read('Number of particles',np)
      ! Euler-Euler or Euler Lagrange
      call param_read('EL simulation', lagrange,default=.true.)
      ! Select particle to turn monitor
      call param_read('Tp index',lp%Tp_ind)
      ! Select whether to turn off particle
      call param_read('Tp off',lp%Tp_off,default=.false.)
      ! Select lattice or random configuration
      call param_read('Lattice configuration',lattice,default=.true.)

      if (lagrange) then
         ! Root process initializes particles
         if (lattice) then
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
               VFavg=np*Volp/(Hbed*cfg%yL*cfg%zL)
               ! Distribute particles
               do i=1,np
                  ! Give position
                  ix = (i-1)/(npy*npz)
                  iy = (i-1-npy*npz*ix)/npz
                  iz = i-1-npy*npz*ix-npz*iy
                  lp%p(i)%pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart+0.5_WP
                  lp%p(i)%pos(2) = lp%cfg%y(lp%cfg%jmin)+(real(iy,WP)+0.5_WP)*Lparty+0.5_WP
                  lp%p(i)%pos(3) = lp%cfg%z(lp%cfg%kmin)+(real(iz,WP)+0.5_WP)*Lpartz+0.5_WP
                  ! Locate the particle on the mesh
                  lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  ! Set the diameter
                  lp%p(i)%d=dp
                  ! Activate the particle
                  lp%p(i)%flag=0
                  ! Set the temperature
                  lp%p(i)%T=Tp
                  ! Give zero velocity
                  lp%p(i)%vel=0.0_WP
                  ! Give zero dt
                  lp%p(i)%dt=0.0_WP
                  ! Set ID
                  if (lp%p(i)%pos(1).le.x0) then
                     lp%p(i)%id=0
                  else
                     lp%p(i)%id=-1
                  end if              
               end do
            end if
         else
            Volp = Pi/6.0_WP*dp**3
            !Vbox=Hbed*cfg%yL*cfg%zL
            Vbox=0.0_WP
            do k=lp%cfg%kmin_,lp%cfg%kmax_
               do j=lp%cfg%jmin_,lp%cfg%jmax_
                  do i=lp%cfg%imin_,fs%cfg%imax_
                     Vbox=Vbox+lp%cfg%vol(i,j,k)*lp%cfg%VF(i,j,k)
                  end do
               end do
            end do
            np = int(VFavg*Vbox/Volp)
            call lp%resize(np)
            ! Allocate particle in cell arrays
            allocate(npic(     lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); npic=0
            allocate(ipic(1:40,lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); ipic=0
            VFavg=real(np,WP)*Volp/Vbox
            do i=1,np
               print *, i
               ! Set the diameter
               lp%p(i)%d=dp
               ! Give position (avoid overlap)
               overlap=.true.
               do while (overlap)
                  lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin_),lp%cfg%x(lp%cfg%imax_+1)),&
                       &       random_uniform(lp%cfg%y(lp%cfg%jmin_),lp%cfg%y(lp%cfg%jmax_+1)-dp),&
                       &       random_uniform(lp%cfg%z(lp%cfg%kmin_),lp%cfg%z(lp%cfg%kmax_+1)-dp)]
                  if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=0.0_WP
                  lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  overlap=.false.
                  do kk=lp%p(i)%ind(3)-1,lp%p(i)%ind(3)+1
                     do jj=lp%p(i)%ind(2)-1,lp%p(i)%ind(2)+1
                        do ii=lp%p(i)%ind(1)-1,lp%p(i)%ind(1)+1
                           do nn=1,npic(ii,jj,kk)
                              j=ipic(nn,ii,jj,kk)
                              if (sqrt(sum((lp%p(i)%pos-lp%p(j)%pos)**2)).lt.0.5_WP*(lp%p(i)%d+lp%p(j)%d)) overlap=.true.
                           end do
                        end do
                     end do
                  end do
               end do
               ! Activate the particle
               lp%p(i)%flag=0
               ip=lp%p(i)%ind(1); jp=lp%p(i)%ind(2); kp=lp%p(i)%ind(3)
               npic(ip,jp,kp)=npic(ip,jp,kp)+1
               ipic(npic(ip,jp,kp),ip,jp,kp)=i
               ! Set the temperature
               lp%p(i)%T=Tp
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Set ID
               if (lp%p(i)%pos(1).le.x0) then
                  lp%p(i)%id=0
               else
                  lp%p(i)%id=-1
               end if
            end do
         end if
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()       
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Euler-Lagrange'
            print*,'Number of particles', lp%np
            print*,'Mean volume fraction',VFavg
            print *, 'Location of particle being tracked',lp%p(lp%Tp_ind)%pos
         end if
      else
         np=0
         lp%VF=VFavg
         if (lp%cfg%amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Euler-Euler'
            print*,'Number of particles', lp%np
            print*,'Mean volume fraction',VFavg
         end if
      end if
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=2,nvec=0,name='lpt')
      pmesh%varname(1)='diameter'
      pmesh%varname(2)='Tp'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=lp%p(i)%d
         pmesh%var(2,i)=lp%p(i)%T
      end do
    end block create_pmesh

    
    ! Initialize our scalar fields
    initialize_scalar: block
      use vdscalar_class, only: bcond
      type(bcond), pointer :: mybc
      integer :: i,j,k,n
      real(WP) :: Ti
      ! Read in the intial values
      call param_read('Density',rho)
      call param_read('Initial T',Ti)
      call param_read('Inlet T',SCin)
      call param_read('Diffusivity multiplier',diff_mult)
      ! Assign values
      sc%SC=Ti
      ! Initialize the scalars at the inlet
      call sc%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         sc%SC(i,j,k)=SCin
      end do
      ! Apply all other boundary conditions
      call sc%apply_bcond(time%t,time%dt)
      ! Store density
      sc%rho=rho*(1.0_WP-lp%VF)
      call sc%rho_multiply()
    end block initialize_scalar

    ! Initialize our velocity field
    initialize_velocity: block
      use lowmach_class, only: bcond
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      ! Set density
      fs%rho=sc%rho
      rhof=fs%rho/(1.0_WP-lp%VF)
      ! Read inlet velocity
      call param_read('Inlet velocity',rhoUin)
      ! Set uniform momentum and velocity
      rhoUin=rho*rhoUin
      fs%rhoU=rhoUin; fs%rhoV=0.0_WP; fs%rhoW=0.0_WP
      call fs%rho_divide()
      call fs%apply_bcond(time%t,time%dt)
      call fs%interp_vel(Ui,Vi,Wi)
      resSC=0.0_WP
      call fs%get_div(drhodt=resSC)
      ! Compute MFR through all boundary conditions
      call fs%get_mfr()
    end block initialize_velocity

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='PTHF_verification_EL')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('epsp',lp%VF)
      call ens_out%add_scalar('density',rhof)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('temperature',sc%SC)
      call ens_out%add_scalar('viscosity',fs%visc)
      call ens_out%add_scalar('diffusivity',sc%diff)
      call ens_out%add_scalar('ptke',lp%ptke)
      call ens_out%add_scalar('diff_pt',lp%diff_pt)
      call ens_out%add_scalar('srcT',srcSClp)
      call ens_out%add_scalar('PTHF_source',tmp)
      call ens_out%add_scalar('T_filtered',Tf)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create monitor file
    create_monitor: block
      use string, only: str_medium
      integer :: ii
      character(len=str_medium) :: str
      ! Prepare some info about fields
      real(WP) :: cfl
      call lp%get_cfl(time%dt,cflc=time%cfl)
      call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
      call fs%get_max()
      call lp%get_max()
      call sc%get_max()
      call sc%get_int()
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
      call cflfile%add_column(lp%CFLpt_x,'PT xCFL')
      call cflfile%add_column(lp%CFLpt_y,'PT yCFL')
      call cflfile%add_column(lp%CFLpt_z,'PT zCFL')
      call cflfile%write()
      ! Create scalar monitor
      scfile=monitor(fs%cfg%amRoot,'scalar')
      call scfile%add_column(time%n,'Timestep number')
      call scfile%add_column(time%t,'Time')
      call scfile%add_column(sc%rhomin,'RHOmin')
      call scfile%add_column(sc%rhomax,'RHOmax')
      call scfile%add_column(sc%rhoint,'RHOint')
      str=trim(sc%name)//'min'
      call scfile%add_column(sc%SCmin,trim(str))
      str=trim(sc%name)//'max'
      call scfile%add_column(sc%SCmax,trim(str))
      str=trim(sc%name)//'int'
      call scfile%add_column(sc%SCint,trim(str))
      ! Create LPT monitor
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%VFmean,'VFp mean')
      call lptfile%add_column(lp%VFmax,'VFp max')
      call lptfile%add_column(lp%np,'Particle number')
      call lptfile%add_column(lp%Tp,'Particle Tp')
      call lptfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use parallel, only: parallel_time
    implicit none
    integer :: ii,i,j,k
    real(WP) :: cfl, vel_mag

    ! Perform time integration
    do while (.not.time%done())
       ! Increment time
       call lp%get_cfl(time%dt,cflc=time%cfl)
       call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
       call time%adjust_dt()
       call time%increment()

       ! Particle solver
       lpt: block
         use mathtools, only: Pi
         real(WP) :: fVF,pVF,Rep,frho,fvisc,fdiff,dt,Pr,Nu,theta
         real(WP), dimension(3) :: vel

         ! Get fluid stress
         call fs%get_div_stress(resU,resV,resW)
         ! Filter fluid quantities
         fs%Uold=fs%U; !call lp%filter(fs%Uold); call lp%filter(Ui)
         fs%Vold=fs%V; !call lp%filter(fs%Vold); call lp%filter(Vi)
         fs%Wold=fs%W; !call lp%filter(fs%Wold); call lp%filter(Wi)
         sc%SCold=sc%SC; !call lp%filter(sc%SCold)

         ! ! Re-Apply scalar boundary conditions if pre-filtering
         ! scalar_inlet: block
         !   use vdscalar_class, only: bcond
         !   type(bcond), pointer :: mybc
         !   integer :: n,i,j,k
         !   do ii=1,nscalar
         !      call sc(ii)%get_bcond('inflow',mybc)
         !      do n=1,mybc%itr%no_
         !         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !         sc(ii)%SCold(i,j,k)=SCin(ii)
         !      end do
         !   end do
         ! end block scalar_inlet

         ! Track filtered tempreature for monitoring
         Tf=sc%SCold

         ! Caculate source terms from particles
         if(lagrange) then            
            ! EL Source Term
            ! Advance particles
            call lp%advance(dt=time%dtmid,U=fs%Uold,V=fs%Vold,W=fs%Wold,rho=rhof,visc=fs%visc,diff=sc%diff,&
                 &          stress_x=resU,stress_y=resV,stress_z=resW,T=sc%SCold,&
                 &          srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcSC=srcSClp,fCp=fCp)
         else
            ! EE Source Term
            pVF=VFavg
            fVF=1.0_WP-pVF
            dt=time%dtmid
            ! Calculate source terms for each cell
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     ! Cell fluid properties
                     vel(1)=fs%U(i,j,k); vel(2)=fs%V(i,j,k); vel(3)=fs%W(i,j,k)
                     frho=rhof(i,j,k)
                     fvisc=fs%visc(i,j,k)
                     fdiff=sc%diff(i,j,k)
                     Rep=frho*fVF*norm2(vel)*dp/fvisc
                     
                     ! Compute Nu
                     Pr=fvisc/fdiff
                     Nu=(-0.46_WP+1.77_WP*fVF+0.69_WP*fVF**2)/fVf**3+(1.37_WP-2.4_WP*fVf+1.2_WP*fVf**2)*Rep**(0.7_WP)*Pr**(1.0_WP/3.0_WP)
                     theta=1-1.6_WP*pVF*fVF-3*pVF*fVF**4*exp(-Rep**0.4_WP*pVF)
                     srcSClp(i,j,k)=-dt*3.0_WP*Pi*pVF*Nu*fdiff*(sc%SC(i,j,k)-Tp)/(dp**2*2.0_WP*theta)
                  end do
               end do
            end do
         end if

         ! Compute PTKE and store source terms
         call lp%get_ptke(dt=time%dtmid,Ui=Ui,Vi=Vi,Wi=Wi,visc=fs%visc,rho=rhof,T=SC%SCold,fCp=fCp,&
             &           diff=sc%diff,srcU=resU,srcV=resV,srcW=resW,srcT=tmp,lagrange=lagrange)
         srcUlp=srcUlp+resU; srcVlp=srcVlp+resV; srcWlp=srcWlp+resW; srcSClp=srcSClp+tmp
         
       end block lpt

       ! Remember old scalar
       sc%rhoold=sc%rho
       sc%SCold=sc%SC
       
       ! Remember old velocity and momentum
       fs%RHOold=fs%RHO
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          ! ============= SCALAR SOLVER =======================
          sc%diff=sc%diff*diff_mult
          !lp%diff_pt=sc%diff
          ! Build mid-time scalar
          sc%SC=0.5_WP*(sc%SC+sc%SCold)
             
          ! Explicit calculation of drhoSC/dt from scalar equation
          call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

          ! Assemble explicit residual
          resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
             
          ! Heat & mass transfer from particles
          do k=sc%cfg%kmin_,sc%cfg%kmax_
             do j=sc%cfg%jmin_,sc%cfg%jmax_
                do i=sc%cfg%imin_,sc%cfg%imax_
                   resSC(i,j,k)=resSC(i,j,k)+srcSClp(i,j,k)
                end do
             end do
          end do

          ! Form implicit residual
          call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
             
          ! Apply this residual
          sc%SC=2.0_WP*sc%SC-sc%SCold+resSC

          ! Apply other boundary conditions on the resulting field
          call sc%apply_bcond(time%t,time%dt)

          ! Apply scalar boundary conditions
          scalar_bcond: block
            use vdscalar_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: n,i,j,k
            call sc%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               sc%SC(i,j,k)=SCin
            end do
          end block scalar_bcond
          call sc%rho_multiply()
          sc%diff=sc%diff/diff_mult

          ! ===================================================

          ! ============ VELOCITY SOLVER ======================
          
          ! Build n+1 density
          fs%rho=0.5_WP*(sc%rho+sc%rhoold)

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
               fs%rhoU(i,j,k)=rhoUin
               fs%rhoV(i,j,k)=0.0_WP
               fs%rhoW(i,j,k)=0.0_WP
            end do
            call fs%rho_divide()
          end block dirichlet_velocity

          ! Solve Poisson equation
          ! Compute rate-of-change of density
          resSC=0.0_WP
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

       ! Recompute molecular density, interpolated velocity, and divergence
       rhof=sc%rho/(1.0_WP-lp%VF)
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=resSC)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=lp%p(i)%d
               pmesh%var(2,i)=lp%p(i)%T
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call fs%get_max()
       call lp%get_max()
       call sc%get_max()
       call sc%get_int()
       call mfile%write()
       call cflfile%write()
       call lptfile%write()
       call scfile%write()

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
    deallocate(resU,resV,resW,resSC,srcUlp,srcVlp,srcWlp,srcSClp,Ui,Vi,Wi,rhof,tmp)

  end subroutine simulation_final

end module simulation
