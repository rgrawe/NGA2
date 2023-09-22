!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,             only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use lowmach_class,     only: lowmach
   use sgsmodel_class,    only: sgsmodel
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use datafile_class,     only: datafile
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single-phase incompressible flow solver and corresponding time tracker
   type(lowmach),     public :: fs
   type(timetracker), public :: time
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile,lptfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi,rho0,dRHOdt
   real(WP), dimension(:,:,:),   allocatable :: srcUlp,srcVlp,srcWlp
   real(WP), dimension(:,:,:),   allocatable :: tmp1,tmp2,tmp3
   real(WP), dimension(:,:,:,:), allocatable :: SR

   !> Max timestep size for LPT
   real(WP) :: lp_dt,lp_dt_max

   !> Fluid viscosity and density
   real(WP) :: visc,rho

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW

   !> Event for post-processing
   type(event) :: ppevt
   real(WP), dimension(:), allocatable :: Uavg,Uavg_,Vavg,Vavg_,Wavg_,Wavg,U2_,U2,V2_,V2,W2_,W2
   real(WP), dimension(:), allocatable :: vol,vol_,rhoVol,rhoVol_,VF,ptke_,ptke,sgsvisc_,sgsvisc
   
 contains

   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_vel()
     use string,    only: str_medium
     use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
     use parallel,  only: MPI_REAL_WP
     implicit none
      integer :: iunit,ierr,i,j,k
      character(len=str_medium) :: filename,timestamp
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j) = vol_(j)+fs%cfg%vol(i,j,k)
               rhoVol_(j) = rhoVol_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*time%dt
               Uavg_(j)=Uavg_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%U(i,j,k)*time%dt
               Vavg_(j)=Vavg_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%V(i,j,k)*time%dt
               Wavg_(j)=Wavg_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%W(i,j,k)*time%dt
               U2_(j)=U2_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%U(i,j,k)**2*time%dt
               V2_(j)=V2_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%V(i,j,k)**2*time%dt
               W2_(j)=W2_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*fs%W(i,j,k)**2*time%dt
               ptke_(j)=ptke_(j)+fs%cfg%vol(i,j,k)*fs%RHO(i,j,k)*lp%ptke(i,j,k)*time%dt
               sgsvisc_(j)=sgsvisc_(j)+fs%cfg%vol(i,j,k)*sgs%visc(i,j,k)*time%dt
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(vol_, vol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(rhoVol_, rhoVol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Vavg_,Vavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Wavg_,Wavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(U2_,U2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(V2_,V2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(W2_,W2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(ptke_,ptke,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(sgsvisc_,sgsvisc,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      do j=fs%cfg%jmin,fs%cfg%jmax
         if (vol(j).gt.0.0_WP) then
            VF(j)=rhoVol(j)/vol(j)/rho
            Uavg(j)=Uavg(j)/rhoVol(j)
            Vavg(j)=Vavg(j)/rhoVol(j)
            Wavg(j)=Wavg(j)/rhoVol(j)
            U2(j)=U2(j)/rhoVol(j)
            V2(j)=V2(j)/rhoVol(j)
            W2(j)=W2(j)/rhoVol(j)
            ptke(j)=ptke(j)/rhoVol(j)
            sgsvisc(j)=sgsvisc(j)/vol(j)
         else
            VF(j)=0.0_WP
            Uavg(j)=0.0_WP
            Vavg(j)=0.0_WP
            Wavg(j)=0.0_WP
            U2(j)=0.0_WP
            V2(j)=0.0_WP
            W2(j)=0.0_WP
            ptke(j)=0.0_WP
            sgsvisc(j)=0.0_WP
         end if
      end do
      ! If root, print it out
      if (fs%cfg%amRoot) then
         filename='Uavg_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'Height','VFavg','Uavg','Uvar','Vvar','Wvar','ptke','mu_t'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') fs%cfg%ym(j),VF(j),Uavg(j),&
                 U2(j)-Uavg(j)**2,V2(j)-Vavg(j)**2,W2(j)-Wavg(j)**2,ptke(j),sgsvisc(j)
         end do
         close(iunit)
      end if
    end subroutine postproc_vel

     !> Function that localizes the bottom (y-) of the domain
   function bottom_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmin) isIn=.true.
   end function bottom_of_domain

   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain


   !> Initialization of problem solver
   subroutine simulation_init
     use param, only: param_read
     implicit none

     ! Initialize time tracker with 2 subiterations
     initialize_timetracker: block
       time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
       end block initialize_timetracker
       
     ! Handle restart/saves here
     restart_and_save: block
       character(len=str_medium) :: timestamp
       ! Create event for saving restart files
       save_evt=event(time,'Restart output')
       call param_read('Restart output period',save_evt%tper)
       ! Check if we are restarting
       call param_read(tag='Restart from',val=timestamp,short='r',default='')
       restarted=.false.; if (len_trim(timestamp).gt.0) restarted=.true.
       if (restarted) then
          ! If we are, read the name of the directory
          call param_read('Restart from',timestamp,'r')
          ! Read the datafile
          df=datafile(pg=cfg,fdata='restart/data_'//trim(adjustl(timestamp)))
       else
          ! Prepare a new directory for storing files for restart
          call execute_command_line('mkdir -p restart')
          ! If we are not restarting, we will still need a datafile for saving restart files
          df=datafile(pg=cfg,filename=trim(cfg%name),nval=4,nvar=3)
          df%valname(1)='t'
          df%valname(2)='dt'
          df%valname(3)='meanU'
          df%valname(4)='meanW'
          df%varname(1)='U'
          df%varname(2)='V'
          df%varname(3)='W'
       end if
     end block restart_and_save

     ! Revisit timetracker to adjust time and time step values if this is a restart
     update_timetracker: block
       if (restarted) then
          call df%pullval(name='t' ,val=time%t )
          call df%pullval(name='dt',val=time%dt)
          time%told=time%t-time%dt
       end if
     end block update_timetracker

       ! Create a low Mach flow solver with bconds
       create_flow_solver: block
         use lowmach_class,   only: dirichlet
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define boundary conditions
         call fs%add_bcond(name='bottom',type=dirichlet,locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top',type=dirichlet,locator=top_of_domain,face='y',dir=+1,canCorrect=.true. )
         ! Assign constant density
         call param_read('Density',rho); fs%rho=rho
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
       end block create_flow_solver


       ! Allocate work arrays
       allocate_work_arrays: block
         allocate(dRHOdt  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcUlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcVlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcWlp  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho0    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR    (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp1    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp2    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp3    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
       end block allocate_work_arrays

       ! Initialize our LPT solver
       initialize_lpt: block
         use random, only: random_uniform
         use mathtools, only: Pi
         real(WP) :: dp
         integer :: i,j,np,ierr
         logical :: overlap
         character(len=str_medium) :: timestamp
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the input
         call param_read('Drag model',lp%drag_model,default='Tenneti')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get number of particles from the input
         call param_read('Number of particles',np)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Maximum timestep size used for particles
         call param_read('Particle timestep size',lp_dt_max,default=huge(1.0_WP))
         lp_dt=lp_dt_max
         if (restarted) then
            call param_read('Restart from',timestamp,'r')
            ! Read the part file
            call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
            do i=1,lp%np_
               if (lp%p(i)%pos(2).lt.(lp%cfg%y(lp%cfg%jmin)+0.51_WP*lp%p(i)%d).and.(lp%p(i)%id.ne.-1)) then
                  lp%p(i)%id=-1
                  lp%p(i)%vel=0.0_WP
                  lp%p(i)%angVel=0.0_WP
               end if
            end do
         else
            if (lp%cfg%amRoot) then
               call lp%resize(np)
               ! Distribute particles
               do i=1,np
                  ! Set the diameter
                  lp%p(i)%d=dp
                  ! Give position (avoid overlap)
                  overlap=.true.
                  do while(overlap)
                     lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin)+0.5_WP*dp,lp%cfg%x(lp%cfg%imax+1)-0.5_WP*dp),&
                          &       random_uniform(lp%cfg%y(lp%cfg%jmin)+0.5_WP*dp,lp%cfg%y(lp%cfg%jmin)+0.75_WP*lp%cfg%yL),&
                          &       random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
                     if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin_)
                     overlap=.false.
                     check: do j=1,i-1
                        if (sqrt(sum((lp%p(i)%pos-lp%p(j)%pos)**2)).lt.0.5_WP*(lp%p(i)%d+lp%p(j)%d)) then
                           overlap=.true.
                           exit check
                        end if
                     end do check
                  end do
                  !print *, real(i,WP)/real(np,WP)*100.0_WP,'%'
                  ! Give id
                  lp%p(i)%id=int(i,8)
                  ! Set the temperature
                  lp%p(i)%T=298.15_WP
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

            if (lp%cfg%amRoot) then
               print*,"===== Particle Setup Description ====="
               print*,'Number of particles', np
            end if
         end if
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         ! Get initial particle volume fraction
         call lp%update_VF()
       end block initialize_lpt


       ! Create partmesh object for Lagrangian particle output
       create_pmesh: block
         integer :: i
        pmesh=partmesh(nvar=1,nvec=1,name='lpt')
        pmesh%varname(1)='diameter'
        pmesh%vecname(1)='velocity'
        call lp%update_partmesh(pmesh)
        do i=1,lp%np_
           pmesh%var(1,i)=lp%p(i)%d
           pmesh%vec(:,1,i)=lp%p(i)%vel
        end do
      end block create_pmesh

      ! Initialize our velocity field
      initialize_velocity: block
        use param, only: param_read
        use lowmach_class, only: bcond
        integer :: n,i,j,k
        type(bcond), pointer :: mybc
        call param_read('Ubulk',Ubulk)
        call param_read('Wbulk',Wbulk,default=0.0_WP)
        if (restarted) then
           call df%pullvar(name='U',var=fs%U)
           call df%pullvar(name='V',var=fs%V)
           call df%pullvar(name='W',var=fs%W)
           call df%pullval(name='meanU',val=meanU)
           call df%pullval(name='meanW',val=meanW)
        else
           ! Zero initial field
           fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
           ! Initialize velocity based on specified bulk
           where (fs%umask.eq.0) fs%U=Ubulk
           where (fs%wmask.eq.0) fs%W=Wbulk
           meanU=Ubulk
           meanW=Wbulk
        end if
        ! Set no-slip walls
        call fs%get_bcond('bottom',mybc)
        do n=1,mybc%itr%no_
           i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
           fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
        end do
        call fs%get_bcond('top',mybc)
        do n=1,mybc%itr%no_
           i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
           fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
        end do
        ! Set density from particle volume fraction and store initial density
        fs%rho=rho*(1.0_WP-lp%VF)
        rho0=rho
        ! Form momentum
        call fs%rho_multiply
        ! Apply all other boundary conditions
        call fs%interp_vel(Ui,Vi,Wi)
        dRHOdt=0.0_WP
        call fs%get_div(drhodt=dRHOdt)
        ! Compute MFR through all boundary conditions
        call fs%get_mfr()
      end block initialize_velocity


      ! Create an LES model
      create_sgs: block
        sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('div',fs%div)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('sgs_visc',sgs%visc)
         call ens_out%add_scalar('VF',lp%VF)
         call ens_out%add_scalar('PTKE',lp%ptke)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
       end block create_ensight

       
      ! Create a monitor file
      create_monitor: block
        ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call lp%get_max()
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
         ! Create LPT monitor
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp_dt,'Particle dt')
         call lptfile%add_column(lp%VFmean,'VFp mean')
         call lptfile%add_column(lp%VFmax,'VFp max')
         call lptfile%add_column(lp%VFvar,'VFp var')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%ncol,'Collision number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%write()
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
         ! Create forcing monitor
         forcefile=monitor(fs%cfg%amRoot,'forcing')
         call forcefile%add_column(time%n,'Timestep number')
         call forcefile%add_column(time%t,'Time')
         call forcefile%add_column(meanU,'Bulk U')
         call forcefile%add_column(meanW,'Bulk W')
         call forcefile%write()
      end block create_monitor
      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Allocate vertical line storage
         allocate(vol_ (fs%cfg%jmin:fs%cfg%jmax)); vol_ =0.0_WP
         allocate(vol  (fs%cfg%jmin:fs%cfg%jmax)); vol  =0.0_WP
         allocate(rhoVol_(fs%cfg%jmin:fs%cfg%jmax)); rhoVol_ =0.0_WP
         allocate(rhoVol (fs%cfg%jmin:fs%cfg%jmax)); rhoVol  =0.0_WP
         allocate(VF(fs%cfg%jmin:fs%cfg%jmax)); VF =0.0_WP
         allocate(Uavg (fs%cfg%jmin:fs%cfg%jmax)); Uavg =0.0_WP
         allocate(Uavg_(fs%cfg%jmin:fs%cfg%jmax)); Uavg_=0.0_WP
         allocate(Vavg (fs%cfg%jmin:fs%cfg%jmax)); Vavg =0.0_WP
         allocate(Vavg_(fs%cfg%jmin:fs%cfg%jmax)); Vavg_=0.0_WP
         allocate(Wavg (fs%cfg%jmin:fs%cfg%jmax)); Wavg =0.0_WP
         allocate(Wavg_(fs%cfg%jmin:fs%cfg%jmax)); Wavg_=0.0_WP
         allocate(U2 (fs%cfg%jmin:fs%cfg%jmax)); U2 =0.0_WP
         allocate(U2_(fs%cfg%jmin:fs%cfg%jmax)); U2_=0.0_WP
         allocate(V2 (fs%cfg%jmin:fs%cfg%jmax)); V2 =0.0_WP
         allocate(V2_(fs%cfg%jmin:fs%cfg%jmax)); V2_=0.0_WP
         allocate(W2 (fs%cfg%jmin:fs%cfg%jmax)); W2 =0.0_WP
         allocate(W2_(fs%cfg%jmin:fs%cfg%jmax)); W2_=0.0_WP
         allocate(ptke_(fs%cfg%jmin:fs%cfg%jmax)); ptke_=0.0_WP
         allocate(ptke (fs%cfg%jmin:fs%cfg%jmax)); ptke =0.0_WP
         allocate(sgsvisc (fs%cfg%jmin:fs%cfg%jmax)); sgsvisc =0.0_WP
         allocate(sgsvisc_(fs%cfg%jmin:fs%cfg%jmax)); sgsvisc_=0.0_WP
         ! Perform the output
         if (ppevt%occurs()) call postproc_vel()
       end block create_postproc


     end subroutine simulation_init

     !To Do - Particle settling


   !> Time integrate our problem
   subroutine simulation_run
     implicit none
     real(WP) :: cfl
     
      ! Perform time integration
     do while (.not.time%done())

         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old density, velocity, and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Particle update
         lpt: block
           integer :: i
           real(WP) :: dt_done,mydt
           ! ‘Glue’ particles to bottom wall
           ! do i=1,lp%np_
           !    if (lp%p(i)%pos(2).lt.(lp%cfg%y(lp%cfg%jmin)+0.5_WP*lp%p(i)%d).and.(lp%p(i)%id.ne.-1)) then
           !       lp%p(i)%id=-1
           !       lp%p(i)%vel=0.0_WP
           !       lp%p(i)%angVel=0.0_WP
           !    end if
           ! end do
           ! Get fluid stress
           call fs%get_div_stress(resU,resV,resW)
           ! Zero-out LPT source terms
           srcUlp=0.0_WP; srcVlp=0.0_WP; srcWlp=0.0_WP
           ! Sub-iteratore
           call lp%get_cfl(lp_dt,cflc=cfl,cfl=cfl)
           if (cfl.gt.0.0_WP) lp_dt=min(lp_dt*time%cflmax/cfl,lp_dt_max)
           dt_done=0.0_WP
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
              mydt=min(lp_dt,time%dtmid-dt_done)
              ! Collide and advance particles
              call lp%collide(dt=mydt)
              call lp%advance(dt=mydt,U=fs%U,V=fs%V,W=fs%W,rho=rho0,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
                   srcU=tmp1,srcV=tmp2,srcW=tmp3)
              srcUlp=srcUlp+tmp1
              srcVlp=srcVlp+tmp2
              srcWlp=srcWlp+tmp3
              ! Increment
              dt_done=dt_done+mydt
           end do
           ! Compute PTKE and store source terms
           call lp%get_ptke(dt=time%dtmid,Ui=Ui,Vi=Vi,Wi=Wi,visc=fs%visc,rho=rho0,srcU=tmp1,srcV=tmp2,srcW=tmp3)
           srcUlp=srcUlp+tmp1
           srcVlp=srcVlp+tmp2
           srcWlp=srcWlp+tmp3
           ! Update density based on particle volume fraction
           fs%rho=rho*(1.0_WP-lp%VF)
           dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid
         end block lpt

         ! Turbulence modeling
         sgs_modeling: block
           use sgsmodel_class, only: dynamic_smag
           call fs%get_strainrate(SR)
           call sgs%get_visc(type=dynamic_smag,dt=time%dtold,rho=rho0,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
           fs%visc=visc+sgs%visc
         end block sgs_modeling

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

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

          ! Add body forcing
            forcing: block
               use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
               use parallel, only: MPI_REAL_WP
               integer :: i,j,k,ierr
               real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
               myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) then
                           myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                           myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                        end if
                        if (fs%wmask(i,j,k).eq.0) then
                           myW   =myW   +fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))
                           myWvol=myWvol+fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)
                        end if
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               where (fs%umask.eq.0) resU=resU+fs%rho*(Ubulk-meanU)
               call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanW=meanW/Wvol
               where (fs%wmask.eq.0) resW=resW+fs%rho*(Wbulk-meanW)
            end block forcing

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
            call fs%get_bcond('bottom',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%rhoU(i,j,k)=0.0_WP; fs%rhoV(i,j,k)=0.0_WP; fs%rhoW(i,j,k)=0.0_WP
               fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
            end do
             call fs%get_bcond('top',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%rhoU(i,j,k)=0.0_WP; fs%rhoV(i,j,k)=0.0_WP; fs%rhoW(i,j,k)=0.0_WP
               fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
            end do
          end block dirichlet_velocity


          ! Solve Poisson equation
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
                 pmesh%vec(:,1,i)=lp%p(i)%vel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call forcefile%write()
         call lptfile%write()

         ! Specialized post-processing
         if (ppevt%occurs()) call postproc_vel()

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              character(len=str_medium) :: timestamp
              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Populate df and write it
              call df%pushval(name='t' ,val=time%t     )
              call df%pushval(name='dt',val=time%dt    )
              call df%pushvar(name='U' ,var=fs%U       )
              call df%pushvar(name='V' ,var=fs%V       )
              call df%pushvar(name='W' ,var=fs%W       )
              call df%pushval(name='meanU', val=meanU  )
              call df%pushval(name='meanW', val=meanW  )
              call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
              ! Write particle file
              call lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
            end block save_restart
         end if

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
      deallocate(resU,resV,resW,srcUlp,srcVlp,srcWlp,Ui,Vi,Wi,dRHOdt,SR,tmp1,tmp2,tmp3)
      deallocate(vol,vol_,rhoVol,rhoVol_,VF,Uavg,Uavg_,Vavg,Vavg_,Wavg,Wavg_,U2_,U2,V2_,V2,W2_,W2,ptke_,ptke,sgsvisc_,sgsvisc)

   end subroutine simulation_final





end module simulation
