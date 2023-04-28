!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use df_class,          only: dfibm
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a direct forcing solver, a vfs and tpns solver, and corresponding time tracker
   type(hypre_str),   public :: ps,vs
   type(vfs),         public :: vf
   type(tpns),        public :: fs
   type(dfibm),       public :: df
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,ibmfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Problem definition
	real(WP), dimension(3) :: center
	real(WP) :: radius

contains


   !> Function that defines a level set function for a falling drop problem
	function levelset_falling_drop(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G
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
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      

      ! Initialize our VOF solver and field
	   create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize to a droplet and a pool
         center=[0.0_WP,0.06_WP,0.0_WP]
         radius=0.01_WP
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
      
      
      ! Initialize our direct forcing solver
      initialize_df: block
         use mathtools, only: twoPi,arctan
         integer :: i,j,k,np
         real(WP) :: Dcyl,Ycyl,amp,freq,theta,r

         ! Create solver
         df=dfibm(cfg=cfg,name='IBM')

         ! Read cylinder properties
         call param_read('Number of markers',np)
         call param_read('Cylinder diameter',Dcyl)
         call param_read('Cylinder position',Ycyl)
         call param_read('Perturbation amp',amp,default=0.05_WP*Dcyl)
         call param_read('Perturbation freq',freq,default=0.0_WP*Dcyl)
         ! Root process initializes marker particles
         if (df%cfg%amRoot) then
            call df%resize(np)
            ! Distribute marker particles
            do i=1,np
               ! Set various parameters for the marker
               df%p(i)%id =1
               df%p(i)%vel=0.0_WP
               ! Set position
               theta=real(i-1,WP)*twoPi/real(np,WP)
               r=0.5_WP*Dcyl+amp*sin(freq*theta)
               df%p(i)%pos(1)=r*cos(theta)
               df%p(i)%pos(2)=r*sin(theta)
               df%p(i)%pos(3)=0.0_WP
               ! Assign element area
               df%p(i)%dA=twoPi*r/real(np,WP)*df%cfg%zL
               ! Assign outward normal vector
               df%p(i)%norm=df%p(i)%pos/r
               ! Shift cylinder
               df%p(i)%pos(2)=df%p(i)%pos(2)+Ycyl
               ! Locate the particle on the mesh
               df%p(i)%ind=df%cfg%get_ijk_global(df%p(i)%pos,[df%cfg%imin,df%cfg%jmin,df%cfg%kmin])
               ! Activate the particle
               df%p(i)%flag=0
            end do
         end if
         call df%sync()
         
         ! Get initial volume fraction
         call df%update_VF()
         
         ! All processes initialize IBM objects
         call df%setup_obj()
         
         ! Define levelset (only used for visualization)
         do k=df%cfg%kmin_,df%cfg%kmax_
            do j=df%cfg%jmin_,df%cfg%jmax_
               do i=df%cfg%imin_,df%cfg%imax_
                  theta=arctan(df%cfg%xm(i),df%cfg%ym(j)-Ycyl)
                  r=0.5_WP*Dcyl+amp*sin(freq*theta)
                  df%G(i,j,k)=sqrt((df%cfg%xm(i))**2+(df%cfg%ym(j)-Ycyl)**2)-r
               end do
            end do
         end do
         call df%cfg%sync(df%G)
         
         if (df%cfg%amRoot) then
            print*,"===== Direct Forcing Setup Description ====="
            print*,'Number of marker particles', df%np
            print*,'Number of IBM objects', df%nobj
            print*,'Particle spacing / dx', twoPi*r/real(np,WP)/df%cfg%min_meshsize
         end if

      end block initialize_df
      
      
      ! Create partmesh object for marker particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='ibm')
         pmesh%varname(1)='area'
         pmesh%vecname(1)='velocity'
         call df%update_partmesh(pmesh)
         do i=1,df%np_
            pmesh%var(1,i)=df%p(i)%dA
            pmesh%vec(:,1,i)=df%p(i)%vel
         end do
      end block create_pmesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ib_drop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('markers',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('GIB',df%G)
         call ens_out%add_scalar('VOF',vf%VF)
		   call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('pressure',fs%P)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call df%get_max()
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
         ! Create IBM monitor
         ibmfile=monitor(amroot=df%cfg%amRoot,name='ibm')
         call ibmfile%add_column(time%n,'Timestep number')
         call ibmfile%add_column(time%t,'Time')
         call ibmfile%add_column(df%VFmin,'VF min')
         call ibmfile%add_column(df%VFmax,'VF max')
         call ibmfile%add_column(df%Fx,'Fx')
         call ibmfile%add_column(df%Fy,'Fy')
         call ibmfile%add_column(df%Fz,'Fz')
         call ibmfile%add_column(df%np,'Marker number')
         call ibmfile%add_column(df%Umin,'Marker Umin')
         call ibmfile%add_column(df%Umax,'Marker Umax')
         call ibmfile%add_column(df%Vmin,'Marker Vmin')
         call ibmfile%add_column(df%Vmax,'Marker Vmax')
         call ibmfile%add_column(df%Wmin,'Marker Wmin')
         call ibmfile%add_column(df%Wmax,'Marker Wmax')
         call ibmfile%write()
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
         
         ! Remember old VOF
		   vf%VFold=vf%VF

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
		   call fs%get_olddensity(vf=vf)
			
			! VOF solver step
		   call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
			
			! Prepare new staggered viscosity (at n+1)
		   call fs%get_viscosity(vf=vf)
         
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
				
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Add momentum source term from direct forcing
            ibm_correction: block
               integer :: i,j,k
               ! Get IB source terms
               resU=vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
               call df%get_source(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU)
               ! Interpolate to the staggered cells and synchronize
               resU=0.0_WP; resV=0.0_WP; resW=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=sum(fs%itpr_x(:,i,j,k)*df%srcU(i-1:i,j,k))
                        resV(i,j,k)=sum(fs%itpr_y(:,i,j,k)*df%srcV(i,j-1:j,k))
                        resW(i,j,k)=sum(fs%itpr_z(:,i,j,k)*df%srcW(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
               ! Increment velocity field
               fs%U=fs%U+resU
               fs%V=fs%V+resV
               fs%W=fs%W+resW
            end block ibm_correction
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
				call fs%correct_mfr()
				call fs%get_div()
				call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
				fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
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
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call df%update_partmesh(pmesh)
               do i=1,df%np_
                  pmesh%var(1,i)=df%p(i)%dA
                  pmesh%vec(:,1,i)=df%p(i)%vel
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call df%get_max()
         call mfile%write()
         call cflfile%write()
         call ibmfile%write()
         
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
