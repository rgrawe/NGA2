!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Dout
   use hypre_str_class,   only: hypre_str
   use fft3d_class,       only: fft3d
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp),      public :: fs
   type(fft3d),       public :: ps
   type(hypre_str),   public :: vs
   type(sgsmodel),    public :: sgs
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile

   public :: simulation_init,simulation_run,simulation_final

   !> Work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: visc,bforce,rhoUaxial_tgt,rhoUaxial_avg,rhoUtheta_tgt,rhoUtheta_avg
   real(WP) :: swirl_number,swirl_number_tgt


contains


   !> Function that computes swirl number
   function get_swirl_number(U,V,W,R) result(SN)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: U,V,W !< Cell-centered velocity field
      real(WP), intent(in) :: R     !< Typically the outer radius
      real(WP) :: SN
      integer :: i,j,k,ierr
      real(WP) :: theta,radius
      real(WP) :: myaxialflux,mythetaflux,axialflux,thetaflux
      myaxialflux=0.0_WP; mythetaflux=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               radius=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
               theta=atan2(cfg%ym(j),cfg%zm(k))
               myaxialflux=myaxialflux+cfg%vol(i,j,k)*cfg%VF(i,j,k)*fs%rho*Ui(i,j,k)**2
               mythetaflux=mythetaflux+cfg%vol(i,j,k)*cfg%VF(i,j,k)*fs%rho*Ui(i,j,k)*(Vi(i,j,k)*cos(theta)-Wi(i,j,k)*sin(theta))
            end do
         end do
      end do
      call MPI_ALLREDUCE(myaxialflux,axialflux,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(mythetaflux,thetaflux,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      SN=thetaflux/(R*max(axialflux,epsilon(1.0_WP)))
   end function get_swirl_number


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


      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         !ps%maxlevel=14
         !call param_read('Pressure iteration',ps%maxit)
         !call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Uaxial,Utheta,amp,theta,radius
         ! Zero out velocity
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Read velocity field parameters
         call param_read('Bulk axial velocity',Uaxial)
         call param_read('Bulk theta velocity',Utheta)
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         ! Initialize velocity
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  ! U velocity
                  fs%U(i,j,k)=+Uaxial+Uaxial*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Uaxial*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)*cos(8.0_WP*twoPi*fs%cfg%ym(j)/fs%cfg%yL)
                  ! V velocity
                  radius=sqrt(cfg%y(j)**2+cfg%zm(k)**2)
                  theta=atan2(cfg%y(j),cfg%zm(k))
                  fs%V(i,j,k)=+Utheta*radius*cos(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  ! W velocity
                  radius=sqrt(cfg%ym(j)**2+cfg%z(k)**2)
                  theta=atan2(cfg%ym(j),cfg%z(k))
                  fs%W(i,j,k)=-Utheta*radius*sin(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
         ! Get target rhoUaxial
         call cfg%integrate(A=fs%rho*Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/cfg%fluid_vol
         rhoUaxial_tgt=rhoUaxial_avg
         ! Get target rhoUtheta
         resU=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  radius=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
                  theta=atan2(cfg%ym(j),cfg%zm(k))
                  if (radius.gt.0.0_WP) resU(i,j,k)=(Vi(i,j,k)*cos(theta)-Wi(i,j,k)*sin(theta))/radius
               end do
            end do
         end do
         call cfg%integrate(A=resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/cfg%fluid_vol
         rhoUtheta_tgt=rhoUtheta_avg
         ! Compute swirl number and coeff
         swirl_number=get_swirl_number(U=Ui,V=Vi,W=Wi,R=0.5_WP*Dout)
         swirl_number_tgt=swirl_number
      end block initialize_velocity


      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pipe')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('levelset',cfg%Gib)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('visc_sgs',sgs%visc)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(rhoUaxial_avg,'Average rhoUaxial')
         call mfile%add_column(rhoUtheta_avg,'Average rhoUtheta')
         call mfile%add_column(swirl_number,'Swirl number')
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
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Turbulence modeling
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            resU=fs%rho
            call fs%get_gradu(gradU)
            call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
            fs%visc=visc+sgs%visc
         end block sgs_modeling

         ! Calculate body forcing
         calc_bodyforcing: block
            integer :: i,j,k
            real(WP) :: theta,radius
            call cfg%integrate(A=fs%rho*Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/cfg%fluid_vol
            resU=0.0_WP
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     radius=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
                     theta=atan2(cfg%ym(j),cfg%zm(k))
                     if (radius.gt.0.0_WP) resU(i,j,k)=(Vi(i,j,k)*cos(theta)-Wi(i,j,k)*sin(theta))/radius
                  end do
               end do
            end do
            call cfg%integrate(A=resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/cfg%fluid_vol
         end block calc_bodyforcing

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW

            ! Add body forcing (do we need to time it by the volume?)
            add_bodyforcing: block
               integer :: i,j,k
               real(WP) :: theta,radius
               resU=resU+(rhoUaxial_tgt-rhoUaxial_avg)
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        radius=sqrt(cfg%y(j)**2+cfg%zm(k)**2)
                        theta=atan2(cfg%y(j),cfg%zm(k))
                        resV(i,j,k)=resV(i,j,k)+(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*cos(theta)
                        radius=sqrt(cfg%ym(j)**2+cfg%z(k)**2)
                        theta=atan2(cfg%ym(j),cfg%z(k))
                        resW(i,j,k)=resW(i,j,k)-(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*sin(theta)
                     end do
                  end do
               end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
            end block add_bodyforcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                        fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                        fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibforcing

            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)

            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Compute swirl number
         swirl_number=get_swirl_number(U=Ui,V=Vi,W=Wi,R=0.5_WP*Dout)

         ! Perform and output monitoring
         call fs%get_max()
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
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)

   end subroutine simulation_final

end module simulation
