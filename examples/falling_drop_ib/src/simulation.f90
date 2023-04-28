!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use hypre_str_class,   only: hypre_str
  use tpns_class,        only: tpns
  use vfs_class,         only: vfs
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
  type(hypre_str),   public :: ps
  type(hypre_str),   public :: vs
  type(tpns),        public :: fs
  type(vfs),         public :: vf
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight) :: ens_out
  type(event)   :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile

  public :: simulation_init,simulation_run,simulation_final

  !> Private work arrays
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
  real(WP), dimension(:,:,:), allocatable :: Uslip,Vslip,Wslip

  !> Droplet definition
  real(WP), dimension(3) :: center
  real(WP) :: radius

contains


  !> Function that defines a level set function for a falling drop problem
  function levelset_falling_drop(xyz,t) result(G)
    implicit none
    real(WP), dimension(3),intent(in) :: xyz
    real(WP), intent(in) :: t
    real(WP) :: G
    ! Create the droplet
    G=radius-sqrt(sum((xyz-center)**2))
  end function levelset_falling_drop


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize time tracker with 2 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max cfl number',time%cflmax)
      call param_read('Max time',time%tmax)
      time%dt=time%dtmax
      time%itmax=2
    end block initialize_timetracker


    ! Initialize our VOF solver and field
    create_and_initialize_vof: block
      use mms_geom,  only: cube_refine_vol
      use vfs_class, only: lvira,VFhi,VFlo
      integer :: i,j,k,n,si,sj,sk
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP) :: vol,area
      integer, parameter :: amr_ref_lvl=4
      ! Create a VOF solver
      vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
      ! Initialize to a droplet and a pool
      center=[0.0_WP,0.075_WP,0.0_WP]; radius=0.01_WP
      do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         do j=vf%cfg%jmino_,vf%cfg%jmaxo_
            do i=vf%cfg%imino_,vf%cfg%imaxo_
               ! Set cube vertices
               n=0
               do sk=0,1
                  do sj=0,1
                     do si=0,1
                        n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                     end do
                  end do
               end do
               ! Call adaptive refinement code to get volume and barycenters recursively
               vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
               call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_falling_drop,0.0_WP,amr_ref_lvl)
               vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
               if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                  vf%Lbary(:,i,j,k)=v_cent
                  vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
               else
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end if
            end do
         end do
      end do
      ! Update the band
      call vf%update_band()
      ! Perform interface reconstruction from VOF field
      call vf%build_interface()
      ! Set interface planes at the boundaries
      call vf%set_full_bcond()
      ! Create discontinuous polygon mesh from IRL interface
      call vf%polygonalize_interface()
      ! Calculate distance from polygons
      call vf%distance_from_polygon()
      ! Calculate subcell phasic volumes
      call vf%subcell_vol()
      ! Calculate curvature
      call vf%get_curvature()
      ! Reset moments to guarantee compatibility with interface reconstruction
      call vf%reset_volume_moments()
    end block create_and_initialize_vof


    ! Create a two-phase flow solver without bconds
    create_and_initialize_flow_solver: block
      use hypre_str_class, only: pcg_pfmg
      ! Create flow solver
      fs=tpns(cfg=cfg,name='Two-phase NS')
      ! Assign constant viscosity to each phase
      call param_read('Liquid dynamic viscosity',fs%visc_l)
      call param_read('Gas dynamic viscosity',fs%visc_g)
      ! Assign constant density to each phase
      call param_read('Liquid density',fs%rho_l)
      call param_read('Gas density',fs%rho_g)
      ! Read in surface tension coefficient
      call param_read('Surface tension coefficient',fs%sigma)
      ! Assign acceleration of gravity
      call param_read('Gravity',fs%gravity)
      ! Configure pressure solver
      ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
      call param_read('Pressure iteration',ps%maxit)
      call param_read('Pressure tolerance',ps%rcvg)
      ! Configure implicit velocity solver
      vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
      call param_read('Implicit iteration',vs%maxit)
      call param_read('Implicit tolerance',vs%rcvg)
      ! Setup the solver
      call fs%setup(pressure_solver=ps,implicit_solver=vs)
      ! Zero initial field
      fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
      ! Calculate cell-centered velocities and divergence
      call fs%interp_vel(Ui,Vi,Wi)
      call fs%get_div()
    end block create_and_initialize_flow_solver


    ! Slip velocity model for contact line
    init_slip_velocity: block
      use mathtools, only: Pi
      allocate(Uslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Uslip=0.0_WP
      allocate(Vslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vslip=0.0_WP
      allocate(Wslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wslip=0.0_WP
      call param_read('Static contact angle',fs%contact_angle); fs%contact_angle=fs%contact_angle*Pi/180.0_WP
    end block init_slip_velocity


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='FallingDrop')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('VOF',vf%VF)
      call ens_out%add_scalar('pressure',fs%P)
      call ens_out%add_scalar('curvature',vf%curv)
      call ens_out%add_scalar('Gib',cfg%Gib)
      call ens_out%add_vector('slip',Uslip,Vslip,Wslip)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      call vf%get_max()
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
      call mfile%add_column(vf%VFmax,'VOF maximum')
      call mfile%add_column(vf%VFmin,'VOF minimum')
      call mfile%add_column(vf%VFint,'VOF integral')
      call mfile%add_column(fs%divmax,'Maximum divergence')
      call mfile%add_column(fs%psolv%it,'Pressure iteration')
      call mfile%add_column(fs%psolv%rerr,'Pressure error')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(fs%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(fs%CFLst,'STension CFL')
      call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
      call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
      call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
      call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
      call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
      call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
      call cflfile%write()
    end block create_monitor


  end subroutine simulation_init


  !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
  subroutine simulation_run
    use tpns_class, only: static_contact
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old VOF
       vf%VFold=vf%VF

       ! Remember old velocity
       fs%Uold=fs%U
       fs%Vold=fs%V
       fs%Wold=fs%W

       ! Apply time-varying Dirichlet conditions
       ! This is where time-dpt Dirichlet would be enforced

       ! Prepare old staggered density (at n)
       call fs%get_olddensity(vf=vf)

       ! VOF solver step
       call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

       ! Prepare new staggered viscosity (at n+1)
       call fs%get_viscosity(vf=vf)

       ! Evaluate IB slip velocity due to contact lien
       calc_slip_velocity: block
         use irl_fortran_interface
         use mathtools, only: normalize
         integer :: i,j,k
         real(WP), dimension(3) :: nw,ni,nt
         real(WP), parameter :: beta=10.0_WP
         real(WP) :: slip
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  ! Zero out Uslip if there's no interface or wall in the cell
                  Uslip(i,j,k)=0.0_WP; Vslip(i,j,k)=0.0_WP; Wslip(i,j,k)=0.0_WP
                  !if (getNumberOfVertices(vf%interface_polygon(1,i,j,k)).eq.0.or.cfg%VF(i,j,k).eq.1.0_WP) cycle
                  ! There is an interface and a wall, find the slip velocity
                  !nw=cfg%Nib(:,i,j,k); ni=calculateNormal(vf%interface_polygon(1,i,j,k))
                  !slip=beta*fs%sigma*(cos(fs%contact_angle)-dot_product(ni,nw))
                  ! Mess around with simple test
                  slip=0.1_WP
                  nw=cfg%Nib(:,i,j,k); ni=[0.0_WP,0.0_WP,1.0_WP]
                  ! Find slip direction
                  nt=normalize(ni-dot_product(ni,nw)*nw)
                  ! Set slip velocity
                  Uslip(i,j,k)=slip*nt(1); Vslip(i,j,k)=slip*nt(2); Wslip(i,j,k)=slip*nt(3)
               end do
            end do
         end do
       end block calc_slip_velocity

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          ! Build mid-time velocity
          fs%U=0.5_WP*(fs%U+fs%Uold)
          fs%V=0.5_WP*(fs%V+fs%Vold)
          fs%W=0.5_WP*(fs%W+fs%Wold)

          ! Preliminary mass and momentum transport step at the interface
          call fs%prepare_advection_upwind(dt=time%dt)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

          ! Add momentum source terms
          call fs%addsrc_gravity(resU,resV,resW)

          ! Assemble explicit residual
          resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
          resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
          resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

          ! Apply IB forcing to enforce BC at the pipe walls
          !ibforcing: block
          !   integer :: i,j,k
          !   do k=fs%cfg%kmin_,fs%cfg%kmax_
          !      do j=fs%cfg%jmin_,fs%cfg%jmax_
          !         do i=fs%cfg%imin_,fs%cfg%imax_
          !            resU(i,j,k)=resU(i,j,k)-(1.0_WP-sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k)))*fs%rho_U(i,j,k)*fs%U(i,j,k)
          !            resV(i,j,k)=resV(i,j,k)-(1.0_WP-sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k)))*fs%rho_V(i,j,k)*fs%V(i,j,k)
          !            resW(i,j,k)=resW(i,j,k)-(1.0_WP-sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k)))*fs%rho_W(i,j,k)*fs%W(i,j,k)
          !         end do
          !      end do
          !   end do
          !end block ibforcing

          ! Form implicit residuals
          call fs%solve_implicit(time%dt,resU,resV,resW)

          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW

          ! Apply direct IB forcing
          ibforcing: block
            integer :: i,j,k
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     ! Set IB velocity to zero
                     fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                     fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                     fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                     ! Also add slip velocity for contact line
                     fs%U(i,j,k)=fs%U(i,j,k)+(1.0_WP-sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k)))*sum(fs%itpr_x(:,i,j,k)*Uslip(i-1:i,j,k))
                     fs%V(i,j,k)=fs%V(i,j,k)+(1.0_WP-sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k)))*sum(fs%itpr_y(:,i,j,k)*Vslip(i,j-1:j,k))
                     fs%W(i,j,k)=fs%W(i,j,k)+(1.0_WP-sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k)))*sum(fs%itpr_z(:,i,j,k)*Wslip(i,j,k-1:k))
                  end do
               end do
            end do
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
            call fs%cfg%sync(fs%W)
          end block ibforcing

          ! Apply other boundary conditions
          call fs%apply_bcond(time%t,time%dt)

          ! Solve Poisson equation
          call fs%update_laplacian()
          call fs%correct_mfr()
          call fs%get_div()
          call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
          fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
          test: block
            real(WP) :: int
            call cfg%integrate(A=fs%psolv%rhs,integral=int)
            if (cfg%amroot) print*,'>>>   initial integral of pressure RHS=',int
            fs%psolv%rhs=fs%psolv%rhs-int/cfg%fluid_vol
            call cfg%integrate(A=fs%psolv%rhs,integral=int)
            if (cfg%amroot) print*,'>>> corrected integral of pressure RHS=',int
          end block test
          fs%psolv%sol=0.0_WP
          call fs%psolv%solve()
          call fs%shift_p(fs%psolv%sol)

          ! Correct velocity
          call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
          fs%P=fs%P+fs%psolv%sol
          fs%U=fs%U-time%dt*resU/fs%rho_U
          fs%V=fs%V-time%dt*resV/fs%rho_V
          fs%W=fs%W-time%dt*resW/fs%rho_W

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div()

       ! Output to ensight
       if (ens_evt%occurs()) call ens_out%write_data(time%t)

       ! Perform and output monitoring
       call fs%get_max()
       call vf%get_max()
       call mfile%write()
       call cflfile%write()

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
    deallocate(resU,resV,resW,Ui,Vi,Wi)

  end subroutine simulation_final


end module simulation
