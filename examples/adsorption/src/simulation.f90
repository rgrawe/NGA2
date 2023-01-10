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
  !use thermochem,        only: Rcst,lookup_molarmass
  implicit none
  private

  !> Get a scalar solver, an LPT solver, a lowmach solver, and corresponding time tracker
  type(vdscalar), dimension(2), public :: sc
  type(lowmach),     public :: fs
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile,adsfile
  type(monitor) :: lptfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,dRHOdt
  real(WP), dimension(:,:,:), allocatable :: resSC,N2_rho
  real(WP) :: visc,rho,diffusivity,CO2mass_init,CO2mass,N2mass,delta_CO2mass,total_mass

  !> Scalar indices
  integer, parameter :: ind_T  =1
  integer, parameter :: ind_CO2=2

  !> Thermo-chemistry parameters
  real(WP), parameter :: Rcst=8.314_WP         !< J/(mol.K)
  real(WP), parameter :: W_N2  = 28.0135e-3_WP !< kg/mol
  real(WP), parameter :: W_CO2 = 44.0095e-3_WP !< kg/mol
  real(WP), parameter :: W_H2O = 18.0153e-3_WP !< kg/mol
  real(WP) :: Pthermo,Pthermo_old,Wmix_old,mean_rho,Ti
  
contains

  !> Compute pressure from composition
  subroutine get_Pthermo
    real(WP) :: Wmix_avg,rho_avg,PCO2,PN2
    resSC=1.0_WP/(sc(ind_CO2)%SC/W_CO2+(1.0_WP-sc(ind_CO2)%SC)/W_N2)
    call fs%cfg%integrate(sc(ind_CO2)%rho*sc(ind_T)%sc/resSC,integral=Pthermo)
    Pthermo=Pthermo*Rcst/fs%cfg%vol_total
    ! call fs%cfg%integrate(sc(ind_CO2)%rho*sc(ind_T)%sc*sc(ind_CO2)%sc/W_CO2,integral=PCO2)
    ! PCO2=PCO2*Rcst/fs%cfg%vol_total
    ! print *,'PCO2',PCO2
    ! call fs%cfg%integrate(sc(ind_CO2)%rho*sc(ind_T)%sc*(1-sc(ind_CO2)%sc)/W_N2,integral=PN2)
    ! PN2=PN2*Rcst/fs%cfg%vol_total
    ! print *,'PN2',PN2    
  end subroutine get_Pthermo
  
  !> Define here our equation of state - rho(P,T,Yk)
   subroutine get_rho
     real(WP) :: T,Y_CO2,Y_N2,Wmix
     integer :: i,j,k
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           do i=fs%cfg%imino_,fs%cfg%imaxo_
              if (fs%cfg%VF(i,j,k).gt.0.0_WP) then
                 T=sc(ind_T)%SC(i,j,k)
                 Y_CO2=sc(ind_CO2)%SC(i,j,k)
                 Y_N2=1.0_WP-Y_CO2
                 Wmix=1.0_WP/(Y_N2/W_N2+Y_CO2/W_CO2)
                 sc(ind_CO2)%rho(i,j,k)=Pthermo*Wmix/(Rcst*T)*(1.0_WP-lp%VF(i,j,k))
                 sc(ind_T)%rho(i,j,k)=sc(ind_CO2)%rho(i,j,k)
              else
                 sc(ind_CO2)%rho(i,j,k)=1.0_WP
                 sc(ind_T)%rho(i,j,k)=1.0_WP
              end if
           end do
        end do
     end do
   end subroutine get_rho

   
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
                 fs%visc(i,j,k)=0.0_WP
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

   !> Calculate the total mass of CO2 and N2
   subroutine get_mass
     call fs%cfg%integrate(sc(ind_CO2)%rhoSC,integral=CO2mass)
     call fs%cfg%integrate((1.0_WP-sc(ind_CO2)%SC)*sc(ind_CO2)%rho,integral=N2mass)
     delta_CO2mass=CO2mass_init-CO2mass
     call fs%cfg%integrate(sc(ind_CO2)%rho,integral=total_mass)
   end subroutine get_mass
   
   !> Initialization of problem solver
   subroutine simulation_init
     use param, only: param_read
     implicit none

    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max time',time%tmax)
      time%dt=time%dtmax
      time%itmax=2
    end block initialize_timetracker
    
    ! Create a low Mach flow solver without bconds
    ! Density will be computed later
    create_flow_solver: block
      use ils_class,     only: pcg_pfmg,pcg_amg
      ! use lowmach_class, only: dirichlet,clipped_neumann
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
      ! Assign initial viscosity
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


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(dRHOdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays

    ! Create scalar solvers for T and CO2
    create_scalar: block
      use ils_class,      only: gmres
      use vdscalar_class, only: quick
      ! Create scalar solvers
      sc(ind_T)  =vdscalar(cfg=cfg,scheme=quick,name='Temperature')
      sc(ind_CO2)=vdscalar(cfg=cfg,scheme=quick,name='CO2')
      ! Assign constant diffusivity
      call param_read('Thermal diffusivity',diffusivity)
      sc(ind_T)%diff=diffusivity
      call param_read('CO2 diffusivity',diffusivity)
      sc(ind_CO2)%diff=diffusivity
      ! Configure implicit scalar solver
      sc(ind_T)%implicit%maxit=fs%implicit%maxit; sc(ind_T)%implicit%rcvg=fs%implicit%rcvg
      sc(ind_CO2)%implicit%maxit=fs%implicit%maxit; sc(ind_CO2)%implicit%rcvg=fs%implicit%rcvg
      ! Setup the solver
      call sc(ind_T)%setup(implicit_ils=gmres)
      call sc(ind_CO2)%setup(implicit_ils=gmres)
    end block create_scalar

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
      ! Get adsorption model from the input
      call param_read('Adsorption model',lp%ads_model,default='NONE')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from the input
      call param_read('Particle diameter',dp)
      ! Set filter scale to 3.5*dx
      lp%filter_width=3.5_WP*cfg%min_meshsize

      ! Root process initializes particle at center of box
      call param_read('Particle count',np)
      call param_read('Particle temperature',Tp,default=298.15_WP)
      if (lp%cfg%amRoot) then
         call lp%resize(np)
         ! Distribute particles
         do i=1,np
            ! Give position
            lp%p(i)%pos(1) = 0.5_WP*lp%cfg%xL
            lp%p(i)%pos(2) = 0.0_WP
            lp%p(i)%pos(3) = 0.0_WP
            ! Give id
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
            lp%p(i)%col=0.0_WP
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
      call param_read('Collision timescale',lp%Tcol,default=5.0_WP*time%dt)
      lp%Tcol=5.0_WP*time%dt
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',lp%e_n)
      ! Set gravity
      call param_read('Gravity',lp%gravity)
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=2,name='lpt')
      pmesh%varname(1)='diameter'
      pmesh%varname(2)='Mc'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=lp%p(i)%d
         pmesh%var(2,i)=lp%p(i)%Mc
      end do
    end block create_pmesh

    
    ! Initialize our scalar fields
    initialize_scalar: block
      real(WP) :: CO2i
      ! Read in the intial values
      call param_read('Initial T',Ti)
      call param_read('Initial CO2',CO2i)
      call param_read('Pressure',Pthermo)
      ! Assign values
      sc(ind_T)%SC=Ti
      sc(ind_CO2)%SC=CO2i
      ! Get density from pressure, temperature, and mass fraction
      call get_rho()
      fs%rho=sc(ind_CO2)%rho
      call fs%cfg%integrate(sc(ind_CO2)%rho,integral=mean_rho); mean_rho=mean_rho/fs%cfg%vol_total
      call sc(ind_T)%rho_multiply()
      call sc(ind_CO2)%rho_multiply()
      ! Get initial mass of CO2 in box
      call fs%cfg%integrate(sc(ind_CO2)%rhoSC,integral=CO2mass_init)
      ! Update transport variables
      call get_viscosity
      call get_thermal_diffusivity
    end block initialize_scalar


    ! Initialize our velocity field
    initialize_velocity: block
      fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
      ! Form momentum
      call fs%rho_multiply
      call fs%interp_vel(Ui,Vi,Wi)
      dRHOdt=0.0_WP
      call fs%get_div(drhodt=dRHOdt)
    end block initialize_velocity

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='adsorption')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('epsp',lp%VF)
      call ens_out%add_scalar('Tf',lp%Tf)
      call ens_out%add_scalar('Yf',lp%Yf)
      call ens_out%add_scalar('density',fs%rho)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('temperature',sc(ind_T)%SC)
      call ens_out%add_scalar('CO2',sc(ind_CO2)%SC)
      call ens_out%add_scalar('viscosity',fs%visc)
      call ens_out%add_scalar('diffusivity',sc(ind_T)%diff)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      call lp%get_max()
      call sc(ind_T)%get_max()
      call sc(ind_CO2)%get_max()
      call get_mass()
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
      call mfile%add_column(sc(ind_T)%SCmax,'Tmax')
      call mfile%add_column(sc(ind_T)%SCmin,'Tmin')
      call mfile%write()
      ! Create adsorption monitor
      adsfile=monitor(fs%cfg%amRoot,'adsorption')
      call adsfile%add_column(time%n,'Timestep number')
      call adsfile%add_column(time%t,'Time')
      call adsfile%add_column(Pthermo,'Pthermo')
      call adsfile%add_column(mean_rho,'mean rho')
      call adsfile%add_column(sc(ind_CO2)%SCmax,'YCO2max')
      call adsfile%add_column(sc(ind_CO2)%SCmin,'YCO2min')
      call adsfile%add_column(total_mass,'Total FLuid Mass')
      call adsfile%add_column(N2mass,'Total N2 Mass')
      call adsfile%add_column(CO2mass,'Total CO2 Mass')
      call adsfile%add_column(delta_CO2mass,'Total CO2 Loss')
      !call adsfile%add_column(lp%Mcmin,'Particle CO2min')
      !call adsfile%add_column(lp%Mcmax,'Particle CO2max')
      call adsfile%add_column(lp%Mcmean,'Particle CO2mean')
      call adsfile%write()
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
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use mathtools, only: twoPi
    implicit none
    integer :: ii
    real(WP) :: total_CO2_mass

    ! Perform time integration
    do while (.not.time%done())
       
       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old scalar
       sc(ind_T)%rhoold=sc(ind_T)%rho
       sc(ind_T)%SCold =sc(ind_T)%SC
       sc(ind_CO2)%rhoold=sc(ind_CO2)%rho
       sc(ind_CO2)%SCold =sc(ind_CO2)%SC

       ! Remember old density, velocity, and momentum
       fs%rhoold=fs%rho
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Collide and advance particles
       call lp%collide(dt=time%dtmid)
       resU=fs%rho/(1.0_WP-lp%VF)
       call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=fs%visc,T=sc(ind_T)%SC,YCO2=sc(ind_CO2)%SC)
       
       ! Perform sub-iterations
       do while (time%it.le.time%itmax)
          ! ============= SCALAR SOLVER =======================
          do ii=1,2
             ! Build mid-time scalar
             sc(ii)%SC=0.5_WP*(sc(ii)%SC+sc(ii)%SCold)

             ! Explicit calculation of drhoSC/dt from scalar equation
             call sc(ii)%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

             ! Assemble explicit residual
             resSC=time%dt*resSC-(2.0_WP*sc(ii)%rho*sc(ii)%SC-(sc(ii)%rho+sc(ii)%rhoold)*sc(ii)%SCold)
             
             ! Add in mass source from particle
             resSC = resSC + lp%srcSC(:,:,:,ii)
             print *,'resSC1',resSC, ii
             
             ! Form implicit residual
             call sc(ii)%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
             print *,'resSC2',resSC, ii

             ! Apply this residual
             sc(ii)%SC=2.0_WP*sc(ii)%SC-sc(ii)%SCold+resSC

             ! Apply other boundary conditions on the resulting field
             call sc(ii)%apply_bcond(time%t,time%dt)
          end do
          ! ===================================================
          
          ! Rescale density and scalars
          rescale_scalars: block
            real(WP) :: drho,rho0
            real(WP) :: mass
            call sc(ind_CO2)%rho_multiply()
            call fs%cfg%integrate(sc(ind_CO2)%rhoold,integral=mean_rho); mean_rho=mean_rho/fs%cfg%vol_total
            rho0=mean_rho
            call fs%cfg%integrate(lp%srcSC(:,:,:,ind_CO2),integral=drho); drho=drho/fs%cfg%vol_total
            mean_rho=mean_rho+drho
            
            ! Rescale density
            sc(ind_CO2)%rho=sc(ind_CO2)%rhoold*mean_rho/rho0
            sc(ind_T)%rho=sc(ind_CO2)%rho
            
            ! Rescale scalars
            sc(ind_CO2)%SC=sc(ind_CO2)%rhoSC/sc(ind_CO2)%rho
            sc(ind_T)%SC=Ti

          end block rescale_scalars
          
          ! ============ UPDATE PROPERTIES ====================
          ! Update pressure
          call get_Pthermo
          ! Update density
          !call get_rho
          ! Compute rate-of-change of density accounting for particles
          dRHOdt=(sc(ind_CO2)%RHO-sc(ind_CO2)%RHOold)/time%dtmid-lp%srcSC(:,:,:,ind_CO2)/time%dtmid
          ! Update transport variables
          call get_viscosity
          call get_thermal_diffusivity
          ! Backup rhoSC
          !do ii=1,2
          !   resSC=sc(ii)%rho*sc(ii)%SC
          !   ! Rescale scalar for conservation
          !   sc(ii)%SC=resSC/sc(ii)%rho
          !end do
          ! UPDATE THE VISCOSITY
          ! call get_viscosity
          ! UPDATE THE DIFFUSIVITY
          ! call get_diffusivity
          ! ===================================================

          ! Build n+1 density
          fs%rho=0.5_WP*(fs%rho+fs%rhoold)

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

          ! Add momentum source terms from lpt
          add_lpt_src: block
            integer :: i,j,k
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                    resU(i,j,k) =resU(i,j,k) +sum(fs%itpr_x(:,i,j,k)*lp%srcU(i-1:i,j,k))
                    resV(i,j,k) =resV(i,j,k) +sum(fs%itpr_y(:,i,j,k)*lp%srcV(i,j-1:j,k))
                    resW(i,j,k) =resW(i,j,k) +sum(fs%itpr_z(:,i,j,k)*lp%srcW(i,j,k-1:k))
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
          call fs%rho_multiply()
          !call sc(ind_CO2)%rho_multiply()

          ! Solve Poisson equation
          !dRhodt=0.0_WP
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

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=dRHOdt)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=lp%p(i)%d
               pmesh%var(2,i)=lp%p(i)%Mc
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call fs%get_max()
       call sc(ind_T)%get_max()
       call sc(ind_CO2)%get_max()
       call get_mass()
       call lp%get_max()
       call mfile%write()
       call cflfile%write()
       call adsfile%write()
       call lptfile%write()

       stop
    end do

    ! Ouput particle bed
    !call lp%write(filename='part.file')

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
    deallocate(resU,resV,resW,resSC,Ui,Vi,Wi,dRHOdt)

  end subroutine simulation_final

end module simulation
