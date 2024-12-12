!> Basic Lagrangian particle solver class:
!> Provides support for Lagrangian-transported objects
module lpt_class
  use precision,      only: WP
  use string,         only: str_medium
  use config_class,   only: config
  use diag_class,     only: diag
  use ddadi_class,    only: ddadi
  use vdscalar_class, only: vdscalar
  use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
  implicit none
  private


  ! Expose type/constructor/methods
  public :: lpt


  !> Memory adaptation parameter
  real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
  real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor

  !> I/O chunk size to read at a time
  integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing

  !> Basic particle object definition
  type :: part
     !> MPI_INTEGER8 data
     integer(kind=8) :: id                !< Particle ID
     !> MPI_DOUBLE_PRECISION data
     real(WP) :: d                        !< Particle diameter
     real(WP), dimension(3) :: pos        !< Particle center coordinates
     real(WP), dimension(3) :: vel        !< Velocity of particle
     real(WP), dimension(3) :: angVel     !< Angular velocity of particle
     real(WP), dimension(3) :: Acol       !< Collision acceleration
     real(WP), dimension(3) :: Tcol       !< Collision torque
     real(WP) :: T                        !< Temperature
     real(WP) :: Mc                       !< Mass of adsorbed CO2
     real(WP) :: MacroCO2                 !< Mass of CO2 in the macropores
     real(WP) :: dt                       !< Time step size for the particle
     !> MPI_INTEGER data
     integer , dimension(3) :: ind        !< Index of cell containing particle center
     integer  :: flag                     !< Control parameter (0=normal, 1=done->will be removed)
  end type part
  !> Number of blocks, block length, and block types in a particle
  integer, parameter                         :: part_nblock=3
  integer           , dimension(part_nblock) :: part_lblock=[1,20,4]
  type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_INTEGER8,MPI_DOUBLE_PRECISION,MPI_INTEGER]
  !> MPI_PART derived datatype and size
  type(MPI_Datatype) :: MPI_PART
  integer :: MPI_PART_SIZE

  !> Lagrangian particle tracking solver object definition
  type :: lpt

     ! This is our underlying config
     class(config), pointer :: cfg                       !< This is the config the solver is build for

     type(diag) :: tridiag                               !< Tridiagonal solver for implicit filter
     type(ddadi) :: implicit                             !< Implicit solver for filtering

     ! This is the name of the solver
     character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)

     ! Particle data
     integer :: np                                       !< Global number of particles
     integer :: np_                                      !< Local number of particles
     integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
     type(part), dimension(:), allocatable :: p          !< Array of particles of type part

     ! Overlap particle (i.e., ghost) data
     integer :: ng_                                      !< Local number of ghosts
     type(part), dimension(:), allocatable :: g          !< Array of ghosts of type part

     ! CFL numbers
     real(WP) :: CFLp_x,CFLp_y,CFLp_z,CFL_col            !< CFL numbers
     real(WP) :: CFLpt_x,CFLpt_y,CFLpt_z                 !< CFL numbers

     ! Particle density
     real(WP) :: rho                                     !< Density of particle
     ! Particle heat capacity
     real(WP) :: pCp                                     !< Heat capacity of particle

     ! Gravitational acceleration
     real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity

     ! Solver parameters
     real(WP) :: nstep=1                                 !< Number of substeps (default=1)
     logical :: use_lift=.false.                         !< Compute lift force on particles
     character(len=str_medium), public :: ads_model      !< Adsorption model
     character(len=str_medium), public :: inert_gas      !< Inert gas

     ! Collisional parameters
     logical :: use_col=.true.                           !< Flag for collisions
     real(WP) :: tau_col                                 !< Characteristic collision time scale
     real(WP) :: e_n                                     !< Normal restitution coefficient
     real(WP) :: e_w                                     !< Wall restitution coefficient
     real(WP) :: mu_f                                    !< Friction coefficient
     real(WP) :: clip_col=0.2_WP                         !< Maximum allowable overlap
     real(WP), dimension(:,:,:),   allocatable :: Wdist  !< Signed wall distance - naive for now (could be redone with FMM)
     real(WP), dimension(:,:,:,:), allocatable :: Wnorm  !< Wall normal function - naive for now (could be redone with FMM)

     ! Injection parameters
     real(WP) :: mfr                                     !< Mass flow rate for particle injection
     real(WP), dimension(3) :: inj_pos                   !< Center location to inject particles
     real(WP), dimension(3) :: inj_vel                   !< Celocity assigned during injection
     real(WP) :: inj_T                                   !< Temperature assigned during injection
     real(WP) :: inj_dmean                               !< Mean diameter assigned during injection
     real(WP) :: inj_dsd                                 !< STD diameter assigned during injection
     real(WP) :: inj_dmin                                !< Min diameter assigned during injection
     real(WP) :: inj_dmax                                !< Max diameter assigned during injection
     real(WP) :: inj_dshift                              !< Diameter shift assigned during injection
     real(WP) :: inj_d                                   !< Diameter to inject particles within

     ! Monitoring info
     real(WP) :: VFmin,VFmax,VFmean,VFvar                !< Volume fraction info
     real(WP) :: Mcmin,Mcmax,Mcmean                      !< Particle CO2 info
     real(WP) :: Tmin,Tmax,Tmean,Tvar                    !< Temperature info
     real(WP) :: Mmin,Mmax,Mmean,Mvar                    !< Mass info
     real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
     real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
     real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
     integer  :: np_new,np_out                           !< Number of new and removed particles
     integer  :: ncol=0                                  !< Number of collisions

     ! Particle volume fraction
     real(WP), dimension(:,:,:), allocatable :: VF       !< Particle volume fraction, cell-centered

     ! Filtering operation
     logical :: implicit_filter                          !< Solve implicitly
     real(WP) :: filter_width                            !< Characteristic filter width
     real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z    !< Divergence operator
     real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z    !< Gradient operator

     ! Pseudo turbulence
     real(WP), dimension(:,:,:), allocatable :: ptke     !< Pseudo-turbulent kinetic energy, cell-centered
     real(WP), dimension(:,:,:), allocatable :: diff_pt  !< Pseudo-turbulent diffusivity, cell-centered
     real(WP), dimension(:,:,:,:), allocatable :: itpr_x,itpr_y,itpr_z !< Interpolation from cell face to cell center

     ! Scalar indices
     integer :: nscalar,ind_CO2,ind_T

   contains
     procedure :: update_partmesh                        !< Update a partmesh object using current particles
     procedure :: collide                                !< Evaluate interparticle collision force
     procedure :: advance                                !< Step forward the particle ODEs
     procedure :: get_rhs                                !< Compute rhs of particle odes
     procedure :: resize                                 !< Resize particle array to given size
     procedure :: resize_ghost                           !< Resize ghost array to given size
     procedure :: recycle                                !< Recycle particle array by removing flagged particles
     procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
     procedure :: share                                  !< Share particles across interprocessor boundaries
     procedure :: read                                   !< Parallel read particles from file
     procedure :: write                                  !< Parallel write particles to file
     procedure :: get_max                                !< Extract various monitoring data
     procedure :: get_cfl                                !< Calculate maximum CFL
     procedure :: update_VF                              !< Compute particle volume fraction
     procedure :: get_ptke                               !< Compute pseudo-turbulent kinetic energy (and Reynolds stress)
     procedure :: filter                                 !< Apply volume filtering to field
     procedure :: inject                                 !< Inject particles at a prescribed boundary
     procedure :: scalar_init                            !< Give lpt access to scalar information
  end type lpt


  !> Declare lpt solver constructor
  interface lpt
     procedure constructor
  end interface lpt

contains


  !> Default constructor for lpt solver
  function constructor(cfg,name) result(self)
    implicit none
    type(lpt) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), optional :: name
    integer :: i,j,k,l

    ! Set the name for the solver
    if (present(name)) self%name=trim(adjustl(name))

    ! Point to pgrid object
    self%cfg=>cfg

    ! Allocate variables
    allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
    self%np_=0; self%np=0
    call self%resize(0)
    self%np_new=0; self%np_out=0

    ! Initialize MPI derived datatype for a particle
    call prepare_mpi_part()

    ! Allocate VF and PTKE on cfg mesh
    allocate(self%VF  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF=0.0_WP
    allocate(self%ptke(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%ptke=0.0_WP
    allocate(self%diff_pt(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%diff_pt=0.0_WP

    ! Zero friction by default
    self%mu_f=0.0_WP

    ! Allocate finite volume divergence operators
    allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    ! Create divergence operator to cell center [xm,ym,zm]
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
             self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
             self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
          end do
       end do
    end do

    ! Allocate finite difference velocity gradient operators
    allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
    allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
    allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
    ! Create gradient coefficients to cell faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
             self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
             self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do

    ! Allocate finite difference density interpolation coefficients
    allocate(self%itpr_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
    allocate(self%itpr_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
    allocate(self%itpr_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
    ! Create density interpolation coefficients to cell face
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             self%itpr_x(:,i,j,k)=self%cfg%dxmi(i)*[self%cfg%xm(i)-self%cfg%x(i),self%cfg%x(i)-self%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
             self%itpr_y(:,i,j,k)=self%cfg%dymi(j)*[self%cfg%ym(j)-self%cfg%y(j),self%cfg%y(j)-self%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
             self%itpr_z(:,i,j,k)=self%cfg%dzmi(k)*[self%cfg%zm(k)-self%cfg%z(k),self%cfg%z(k)-self%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do
    
    ! Loop over the domain and zero divergence in walls
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                self%div_x(:,i,j,k)=0.0_WP
                self%div_y(:,i,j,k)=0.0_WP
                self%div_z(:,i,j,k)=0.0_WP
             end if
          end do
       end do
    end do

    ! Zero out gradient to wall faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
             ! if (self%cfg%VF(i,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
             ! if (self%cfg%VF(i,j,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
             ! if (self%cfg%VF(i,j,k).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
          end do
       end do
    end do

    ! Adjust interpolation coefficients to cell faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             ! Linear interpolation in x
             if (self%cfg%VF(i,j,k).eq.1.and.self%cfg%VF(i-1,j,k).lt.1) self%itpr_x(:,i,j,k)=[0.0_WP,1.0_WP]
             if (self%cfg%VF(i,j,k).lt.1.and.self%cfg%VF(i-1,j,k).eq.1) self%itpr_x(:,i,j,k)=[1.0_WP,0.0_WP]
             ! Linear interpolation in y
             if (self%cfg%VF(i,j,k).eq.1.and.self%cfg%VF(i,j-1,k).lt.1) self%itpr_y(:,i,j,k)=[0.0_WP,1.0_WP]
             if (self%cfg%VF(i,j,k).lt.1.and.self%cfg%VF(i,j-1,k).eq.1) self%itpr_y(:,i,j,k)=[1.0_WP,0.0_WP]
             ! Linear interpolation in z
             if (self%cfg%VF(i,j,k).eq.1.and.self%cfg%VF(i,j,k-1).lt.1) self%itpr_z(:,i,j,k)=[0.0_WP,1.0_WP]
             if (self%cfg%VF(i,j,k).lt.1.and.self%cfg%VF(i,j,k-1).eq.1) self%itpr_z(:,i,j,k)=[1.0_WP,0.0_WP]
          end do
       end do
    end do

    ! Adjust metrics to account for lower dimensionality
    if (self%cfg%nx.eq.1) then
       self%div_x=0.0_WP
       self%grd_x=0.0_WP
    end if
    if (self%cfg%ny.eq.1) then
       self%div_y=0.0_WP
       self%grd_y=0.0_WP
    end if
    if (self%cfg%nz.eq.1) then
       self%div_z=0.0_WP
       self%grd_z=0.0_WP
    end if

    ! Create implicit solver object for filtering
    self%implicit=ddadi(cfg=self%cfg,name='Filter',nst=7)
    self%implicit%stc(1,:)=[ 0, 0, 0]
    self%implicit%stc(2,:)=[+1, 0, 0]
    self%implicit%stc(3,:)=[-1, 0, 0]
    self%implicit%stc(4,:)=[ 0,+1, 0]
    self%implicit%stc(5,:)=[ 0,-1, 0]
    self%implicit%stc(6,:)=[ 0, 0,+1]
    self%implicit%stc(7,:)=[ 0, 0,-1]
    call self%implicit%init()

    ! Set filter width to zero by default
    self%filter_width=0.0_WP

    ! Solve implicitly by default
    self%implicit_filter=.true.

    ! Generate a wall distance/norm function
    allocate(self%Wdist(  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
    allocate(self%Wnorm(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
    ! First pass to set correct sign
    do k=self%cfg%kmino_,self%cfg%kmaxo_
       do j=self%cfg%jmino_,self%cfg%jmaxo_
          do i=self%cfg%imino_,self%cfg%imaxo_
             if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                self%Wdist(i,j,k)=-sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
             else
                self%Wdist(i,j,k)=+sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
             end if
             self%Wnorm(:,i,j,k)=0.0_WP
          end do
       end do
    end do
    ! Second pass to compute local distance (get 2 closest cells)
    do l=1,2
       do k=self%cfg%kmino_,self%cfg%kmaxo_
          do j=self%cfg%jmino_,self%cfg%jmaxo_
             do i=self%cfg%imino_+l,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i-l,j,k).lt.0.0_WP) then
                   ! There is a wall at x(i)
                   if (abs(self%cfg%xm(i  )-self%cfg%x(i-l+1)).lt.abs(self%Wdist(i  ,j,k))) then
                      self%Wdist(i  ,j,k)=sign(self%cfg%xm(i  )-self%cfg%x(i-l+1),self%Wdist(i  ,j,k))
                      self%Wnorm(:,i  ,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
                   end if
                   if (abs(self%cfg%xm(i-l)-self%cfg%x(i)).lt.abs(self%Wdist(i-l,j,k))) then
                      self%Wdist(i-l,j,k)=sign(self%cfg%xm(i-l)-self%cfg%x(i),self%Wdist(i-l,j,k))
                      self%Wnorm(:,i-l,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
                   end if
                end if
             end do
          end do
       end do
       do k=self%cfg%kmino_,self%cfg%kmaxo_
          do j=self%cfg%jmino_+l,self%cfg%jmaxo_
             do i=self%cfg%imino_,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i,j-l,k).lt.0.0_WP) then
                   ! There is a wall at y(j)
                   if (abs(self%cfg%ym(j  )-self%cfg%y(j-l+1)).lt.abs(self%Wdist(i,j  ,k))) then
                      self%Wdist(i,j  ,k)=sign(self%cfg%ym(j  )-self%cfg%y(j-l+1),self%Wdist(i,j  ,k))
                      self%Wnorm(:,i,j  ,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
                   end if
                   if (abs(self%cfg%ym(j-l)-self%cfg%y(j)).lt.abs(self%Wdist(i,j-l,k))) then
                      self%Wdist(i,j-l,k)=sign(self%cfg%ym(j-l)-self%cfg%y(j),self%Wdist(i,j-l,k))
                      self%Wnorm(:,i,j-l,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
                   end if
                end if
             end do
          end do
       end do
       do k=self%cfg%kmino_+l,self%cfg%kmaxo_
          do j=self%cfg%jmino_,self%cfg%jmaxo_
             do i=self%cfg%imino_,self%cfg%imaxo_
                if (self%Wdist(i,j,k)*self%Wdist(i,j,k-l).lt.0.0_WP) then
                   ! There is a wall at z(k)
                   if (abs(self%cfg%zm(k  )-self%cfg%z(k-l+1)).lt.abs(self%Wdist(i,j,k  ))) then
                      self%Wdist(i,j,k  )=sign(self%cfg%zm(k  )-self%cfg%z(k-l+1),self%Wdist(i,j,k  ))
                      self%Wnorm(:,i,j,k  )=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
                   end if
                   if (abs(self%cfg%zm(k-l)-self%cfg%z(k)).lt.abs(self%Wdist(i,j,k-l))) then
                      self%Wdist(i,j,k-l)=sign(self%cfg%zm(k-l)-self%cfg%z(k),self%Wdist(i,j,k-l))
                      self%Wnorm(:,i,j,k-l)=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
                   end if
                end if
             end do
          end do
       end do
    end do
    call self%cfg%sync(self%Wdist)
    call self%cfg%sync(self%Wnorm)
    
    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (self%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end function constructor
  

  !> Resolve collisional interaction between particles
  !> Requires tau_col, e_n, e_w and mu_f to be set beforehand
  subroutine collide(this,dt)
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
    integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell

    ! Start by zeroing out the collision force
    zero_force: block
      integer :: i
      do i=1,this%np_
         this%p(i)%Acol=0.0_WP
         this%p(i)%Tcol=0.0_WP
      end do
    end block zero_force

    ! Return if not used
    if (.not.this%use_col) return

    ! Then share particles across overlap
    call this%share()

    ! We can now assemble particle-in-cell information
    pic_prep: block
      use mpi_f08
      integer :: i,ip,jp,kp,ierr
      integer :: mymax_npic,max_npic

      ! Allocate number of particle in cell
      allocate(npic(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); npic=0

      ! Count particles and ghosts per cell
      do i=1,this%np_
         ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
         npic(ip,jp,kp)=npic(ip,jp,kp)+1
      end do
      do i=1,this%ng_
         ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
         npic(ip,jp,kp)=npic(ip,jp,kp)+1
      end do

      ! Get maximum number of particle in cell
      mymax_npic=maxval(npic); call MPI_ALLREDUCE(mymax_npic,max_npic,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)

      ! Allocate pic map
      allocate(ipic(1:max_npic,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); ipic=0

      ! Assemble pic map
      npic=0
      do i=1,this%np_
         ip=this%p(i)%ind(1); jp=this%p(i)%ind(2); kp=this%p(i)%ind(3)
         npic(ip,jp,kp)=npic(ip,jp,kp)+1
         ipic(npic(ip,jp,kp),ip,jp,kp)=i
      end do
      do i=1,this%ng_
         ip=this%g(i)%ind(1); jp=this%g(i)%ind(2); kp=this%g(i)%ind(3)
         npic(ip,jp,kp)=npic(ip,jp,kp)+1
         ipic(npic(ip,jp,kp),ip,jp,kp)=-i
      end do

    end block pic_prep

    ! Finally, calculate collision force
    collision_force: block
      use mpi_f08
      use mathtools, only: Pi,normalize,cross_product
      integer :: i1,i2,ii,jj,kk,nn,ierr
      real(WP) :: d1,m1,d2,m2,d12,m12
      real(WP), dimension(3) :: r1,v1,w1,r2,v2,w2,v12,n12,f_n,t12,f_t
      real(WP) :: k_n,eta_n,k_coeff,eta_coeff,k_coeff_w,eta_coeff_w,rnv,r_influ,delta_n,rtv
      real(WP), parameter :: aclipnorm=1.0e-6_WP
      real(WP), parameter :: acliptan=1.0e-9_WP
      real(WP), parameter :: rcliptan=0.05_WP

      ! Reset collision counter
      this%ncol=0

      ! Precompute coefficients for k and eta
      k_coeff=(Pi**2+log(this%e_n)**2)/this%tau_col**2
      eta_coeff=-2.0_WP*log(this%e_n)/this%tau_col
      k_coeff_w=(Pi**2+log(this%e_w)**2)/this%tau_col**2
      eta_coeff_w=-2.0_WP*log(this%e_w)/this%tau_col

      ! Loop over all local particles
      collision: do i1=1,this%np_

         ! Cycle if id<=0
         if (this%p(i1)%id.le.0) cycle collision
         
         ! Store particle data
         r1=this%p(i1)%pos
         v1=this%p(i1)%vel
         w1=this%p(i1)%angVel
         d1=this%p(i1)%d
         m1=this%rho*Pi/6.0_WP*d1**3+this%p(i1)%mc

         ! First collide with walls
         d12=this%cfg%get_scalar(pos=this%p(i1)%pos,i0=this%p(i1)%ind(1),j0=this%p(i1)%ind(2),k0=this%p(i1)%ind(3),S=this%Wdist,bc='d')
         n12=this%Wnorm(:,this%p(i1)%ind(1),this%p(i1)%ind(2),this%p(i1)%ind(3))
         n12=-normalize(n12+[epsilon(1.0_WP),epsilon(1.0_WP),epsilon(1.0_WP)])
         rnv=dot_product(v1,n12)
         r_influ=min(2.0_WP*abs(rnv)*dt,0.2_WP*d1)
         delta_n=min(0.5_WP*d1+r_influ-d12,this%clip_col*0.5_WP*d1)

         ! Assess if there is collision
         if (delta_n.gt.0.0_WP) then
            ! Normal collision
            k_n=m1*k_coeff_w
            eta_n=m1*eta_coeff_w
            f_n=-k_n*delta_n*n12-eta_n*rnv*n12
            ! Tangential collision
            f_t=0.0_WP
            if (this%mu_f.gt.0.0_WP) then
               t12 = v1-rnv*n12+cross_product(0.5_WP*d1*w1,n12)
               rtv = sqrt(sum(t12*t12))
               if (rnv*dt/d1.gt.aclipnorm) then
                  if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
               else
                  if (rtv*dt/d1.lt.acliptan) rtv=0.0_WP
               end if
               if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
            end if
            ! Calculate collision force
            f_n=f_n/m1; f_t=f_t/m1
            this%p(i1)%Acol=this%p(i1)%Acol+f_n+f_t
            ! Calculate collision torque
            this%p(i1)%Tcol=this%p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t)
         end if

         ! Loop over nearest cells
         do kk=this%p(i1)%ind(3)-1,this%p(i1)%ind(3)+1
            do jj=this%p(i1)%ind(2)-1,this%p(i1)%ind(2)+1
               do ii=this%p(i1)%ind(1)-1,this%p(i1)%ind(1)+1

                  ! Loop over particles in that cell
                  do nn=1,npic(ii,jj,kk)

                     ! Get index of neighbor particle
                     i2=ipic(nn,ii,jj,kk)

                     ! Get relevant data from correct storage
                     if (i2.gt.0) then
                        r2=this%p(i2)%pos
                        v2=this%p(i2)%vel
                        w2=this%p(i2)%angVel
                        d2=this%p(i2)%d
                        m2=this%rho*Pi/6.0_WP*d2**3+this%p(i2)%mc
                     else if (i2.lt.0) then
                        i2=-i2
                        r2=this%g(i2)%pos
                        v2=this%g(i2)%vel
                        w2=this%g(i2)%angVel
                        d2=this%g(i2)%d
                        m2=this%rho*Pi/6.0_WP*d2**3+this%g(i2)%mc
                     end if

                     ! Compute relative information
                     d12=norm2(r1-r2)
                     if (d12.lt.10.0_WP*epsilon(d12)) cycle !< this should skip auto-collision
                     n12=(r2-r1)/d12
                     v12=v1-v2
                     rnv=dot_product(v12,n12)
                     r_influ=min(abs(rnv)*dt,0.1_WP*(d1+d2))
                     delta_n=min(0.5_WP*(d1+d2)+r_influ-d12,this%clip_col*0.5_WP*(d1+d2))

                     ! Assess if there is collision
                     if (delta_n.gt.0.0_WP) then
                        ! Normal collision
                        m12=m1*m2/(m1+m2)
                        k_n=m12*k_coeff
                        eta_n=m12*eta_coeff
                        f_n=-k_n*delta_n*n12-eta_n*rnv*n12
                        ! Tangential collision
                        f_t=0.0_WP
                        if (this%mu_f.gt.0.0_WP) then
                           t12 = v12-rnv*n12+cross_product(0.5_WP*(d1*w1+d2*w2),n12)
                           rtv = sqrt(sum(t12*t12))
                           if (rnv*dt*2.0_WP/(d1+d2).gt.aclipnorm) then
                              if (rtv/rnv.lt.rcliptan) rtv=0.0_WP
                           else
                              if (rtv*dt*2.0_WP/(d1+d2).lt.acliptan) rtv=0.0_WP
                           end if
                           if (rtv.gt.0.0_WP) f_t=-this%mu_f*sqrt(sum(f_n*f_n))*t12/rtv
                        end if
                        ! Calculate collision force
                        f_n=f_n/m1; f_t=f_t/m1
                        this%p(i1)%Acol=this%p(i1)%Acol+f_n+f_t
                        ! Calculate collision torque
                        this%p(i1)%Tcol=this%p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t)
                        ! Add up the collisions
                        this%ncol=this%ncol+1
                     end if
                     
                  end do

               end do
            end do
         end do

         ! Deal with dimensionality
         if (this%cfg%nx.eq.1) then
            this%p(i1)%Acol(1)=0.0_WP
            this%p(i1)%Tcol(2)=0.0_WP
            this%p(i1)%Tcol(3)=0.0_WP
         end if
         if (this%cfg%ny.eq.1) then
            this%p(i1)%Tcol(1)=0.0_WP
            this%p(i1)%Acol(2)=0.0_WP
            this%p(i1)%Tcol(3)=0.0_WP
         end if
         if (this%cfg%nz.eq.1) then
            this%p(i1)%Tcol(1)=0.0_WP
            this%p(i1)%Tcol(2)=0.0_WP
            this%p(i1)%Acol(3)=0.0_WP
         end if

      end do collision

      ! Determine total number of collisions
      call MPI_ALLREDUCE(this%ncol,nn,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%ncol=nn/2

    end block collision_force

  end subroutine collide


  !> Advance the particle equations by a specified time step dt
  !> p%id=0 => no coll, no solve
  !> p%id=-1=> no coll, no move
  subroutine advance(this,dt,U,V,W,rho,visc,diff,stress_x,stress_y,stress_z,T,YCO2,srcU,srcV,srcW,srcSC)
    use mpi_f08, only : MPI_SUM,MPI_INTEGER
    use mathtools, only: Pi
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: diff      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: T         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: YCO2      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcU   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcV   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcW   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,:), intent(inout), optional :: srcSC
    integer :: i,j,k,ierr
    real(WP) :: mydt,dt_done,Ip,dmdt,dcdt,dm,dTdt,dTdt_r,dTemp,mass,fT,fCp
    real(WP), dimension(3) :: acc,dmom
    type(part) :: myp,pold

    ! Zero out source term arrays
    if (present(srcU)) srcU=0.0_WP
    if (present(srcV)) srcV=0.0_WP
    if (present(srcW)) srcW=0.0_WP
    if (present(srcSC)) srcSC=0.0_WP
    
    ! Zero out number of particles removed
    this%np_out=0

    ! Advance the equations
    do i=1,this%np_
       ! Avoid particles with id=0
       if (this%p(i)%id.eq.0) cycle
       ! Create local copy of particle
       myp=this%p(i)
       ! Time-integrate until dt_done=dt
       dt_done=0.0_WP
       do while (dt_done.lt.dt)
          ! Decide the timestep size
          mydt=min(myp%dt,dt-dt_done)
          ! Remember the particle
          pold=myp
          ! Particle moment of inertia per unit mass
          Ip = 0.1_WP*myp%d**2
          ! Advance with Euler prediction
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,diff=diff,stress_x=stress_x,stress_y=stress_y,stress_z=stress_z,T=T,YCO2=YCO2,&
               p=myp,acc=acc,opt_dt=myp%dt,dmdt=dmdt,dcdt=dcdt,dTdt=dTdt,dTdt_r=dTdt_r,fCp=fCp)
          myp%pos=pold%pos+0.5_WP*mydt*myp%vel
          myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity+myp%Acol)
          myp%angVel=pold%angVel+0.5_WP*mydt*myp%Tcol/Ip
          myp%Mc=pold%Mc+0.5_WP*mydt*dmdt
          myp%MacroCO2=pold%MacroCO2+0.5_WP*mydt*dcdt
          myp%T=pold%T+0.5_WP*mydt*dTdt
          ! Correct with midpoint rule
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,diff=diff,stress_x=stress_x,stress_y=stress_y,stress_z=stress_z,T=T,YCO2=YCO2,&
               p=myp,acc=acc,opt_dt=myp%dt,dmdt=dmdt,dcdt=dcdt,dTdt=dTdt,dTdt_r=dTdt_r,fCp=fCp)
          myp%pos=pold%pos+mydt*myp%vel
          myp%vel=pold%vel+mydt*(acc+this%gravity+myp%Acol)
          myp%angVel=pold%angVel+mydt*myp%Tcol/Ip
          myp%Mc=pold%Mc+mydt*dmdt
          myp%MacroCO2=pold%MacroCO2+mydt*dcdt
          myp%T=pold%T+mydt*dTdt
          ! Relocalize
          myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
          ! Send source term back to the mesh
          mass=this%rho*Pi/6.0_WP*myp%d**3+myp%Mc
          dm=mydt*dmdt
          dTemp=(mass*this%pCp*(dTdt-dTdt_r)*mydt)/fCp+(dm*myp%T*this%pCp)/fCp
          dmom=mydt*acc*mass+myp%vel*dm
          if (this%cfg%nx.gt.1.and.present(srcU)) call this%cfg%set_scalar(Sp=-dmom(1),pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=srcU,bc='n')
          if (this%cfg%ny.gt.1.and.present(srcV)) call this%cfg%set_scalar(Sp=-dmom(2),pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=srcV,bc='n')
          if (this%cfg%nz.gt.1.and.present(srcW)) call this%cfg%set_scalar(Sp=-dmom(3),pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=srcW,bc='n')
          if (present(srcSC)) then
             call this%cfg%set_scalar(Sp=-dTemp,pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=srcSC(:,:,:,this%ind_T  ),bc='n')
             call this%cfg%set_scalar(Sp=-dm   ,pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=srcSC(:,:,:,this%ind_CO2),bc='n')
          end if
          
          ! Increment
          dt_done=dt_done+mydt
       end do
       ! Correct the position to take into account periodicity
       if (this%cfg%xper) myp%pos(1)=this%cfg%x(this%cfg%imin)+modulo(myp%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
       if (this%cfg%yper) myp%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(myp%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
       if (this%cfg%zper) myp%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(myp%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
       ! Handle particles that have left the domain
       if (myp%pos(1).lt.this%cfg%x(this%cfg%imin).or.myp%pos(1).gt.this%cfg%x(this%cfg%imax+1)) myp%flag=1
       if (myp%pos(2).lt.this%cfg%y(this%cfg%jmin).or.myp%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) myp%flag=1
       if (myp%pos(3).lt.this%cfg%z(this%cfg%kmin).or.myp%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) myp%flag=1
       ! Relocalize the particle
       myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
       ! Count number of particles removed
       if (myp%flag.eq.1) this%np_out=this%np_out+1
       ! Copy back to particle
       if (myp%id.eq.-1) then
          this%p(i)%T=myp%T
          this%p(i)%Mc=myp%Mc
       else
          this%p(i)=myp
       end if
    end do
    ! Communicate particles
    call this%sync()

    ! Sum up particles removed
    call MPI_ALLREDUCE(this%np_out,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=i

    ! Divide source arrays by volume, sum at boundaries, and volume filter if present
    if (present(srcU)) then
       srcU=srcU/this%cfg%vol; call this%cfg%syncsum(srcU); call this%filter(srcU)
    end if
    if (present(srcV)) then
       srcV=srcV/this%cfg%vol; call this%cfg%syncsum(srcV); call this%filter(srcV)
    end if
    if (present(srcW)) then
       srcW=srcW/this%cfg%vol; call this%cfg%syncsum(srcW); call this%filter(srcW)
    end if
    if (present(srcSC)) then
       do i=1,this%nscalar
          srcSC(:,:,:,i)=srcSC(:,:,:,i)/this%cfg%vol; call this%cfg%syncsum(srcSC(:,:,:,i)); call this%filter(srcSC(:,:,:,i))
       end do
    end if

    ! Recompute volume fraction
    call this%update_VF()

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine advance


  !> Calculate RHS of the particle ODEs
  subroutine get_rhs(this,U,V,W,rho,visc,diff,stress_x,stress_y,stress_z,T,YCO2,p,acc,opt_dt,dmdt,dcdt,dTdt,dTdt_r,fCp)
    use mathtools, only: Pi
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: diff      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: T         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: YCO2      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    type(part), intent(in) :: p
    real(WP), dimension(3), intent(out) :: acc
    real(WP), intent(out) :: opt_dt,dmdt,dcdt,dTdt,fCp,dTdt_r
    real(WP) :: fvisc,fdiff,frho,pVF,fVF,fT,fYCO2,rhop,tau,Re,H_a,fCp_inert,W_inert
    real(WP), dimension(3) :: fvel,fstress
    real(WP), parameter :: W_He = 4.002602e-3_WP
    real(WP), parameter :: W_CO2=44.0095e-3_WP
    real(WP), parameter :: W_H2O=18.0153e-3_WP
    real(WP), parameter :: W_N2=28.0135e-3_WP
    real(WP), parameter :: fCp_CO2=844.0_WP
    real(WP), parameter :: fCp_He=5190.0_WP
    real(WP), parameter :: fCp_N2=1040.0_WP

    ! Select inert gas
    select case(trim(this%inert_gas))
    case('Helium','helium','HELIUM')
       W_inert=W_He
       fCp_inert=fCp_He
    case('Nitrogen','nitrogen','NITROGEN')
       W_inert=W_N2
       fCp_inert=fCp_N2
    end select
    

    ! Interpolate fluid quantities to particle location
    interpolate: block
      ! Interpolate the fluid phase velocity to the particle location
      fvel=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=U,V=V,W=W)
      ! Interpolate the fluid phase stress to the particle location
      fstress=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=stress_x,V=stress_y,W=stress_z)
      ! Interpolate the fluid phase viscosity to the particle location
      fvisc=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=visc,bc='n')
      fvisc=fvisc+epsilon(1.0_WP)
      ! Interpolate the fluid phase thermal diffusivity to the particle location
      fdiff=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=diff,bc='n')
      ! Interpolate the fluid phase density to the particle location
      frho=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=rho,bc='n')
      ! Interpolate the particle volume fraction to the particle location
      pVF=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=this%VF,bc='n')
      fVF=1.0_WP-pVF
      ! Interpolate the fluid temperature to the particle location
      fT=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=T,bc='n')
      ! Interpolate the CO2 mass fraction to the particle location
      fYCO2=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=YCO2,bc='n')
    end block interpolate

    ! Particle Reynolds number
    Re=fVF*frho*norm2(p%vel-fvel)*p%d/fvisc+epsilon(1.0_WP)

    ! Effective particle density
    rhop=this%rho+6.0_WP*p%mc/(pi*p%d**3)

    ! Particle response time
    tau=rhop*p%d**2/(18.0_WP*fvisc)

    ! Compute acceleration due to drag
    compute_drag: block
      real(WP) :: corr,b1,b2
      ! Tenneti and Subramaniam (2011)
      b1=5.81_WP*pVF/fVF**3+0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
      b2=pVF**3*Re*(0.95_WP+0.61_WP*pVF**3/fVF**2)
      corr=fVF*((1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2)           
      ! Return acceleration and optimal timestep size
      acc=(fvel-p%vel)*corr/tau+fstress/rhop
      opt_dt=tau/corr/real(this%nstep,WP)
    end block compute_drag

    ! Compute mass transfer
    mass_transfer: block
      use messager, only: die
      real(WP) :: Cs,A_d,E_d,dS_d,dH_d,PCO2,Pg,k_d,Keq_d,Cc,dcsdt !Lee parameters
      real(WP) :: qs,b_t,dH_t,b0_t,T0_t,th_t,th0_t,alpha_t,mc_s,k_t,k0_t,Ea_t,mc_clip !Bos parameters
      real(WP) :: T_ref,D_ref,D_m,D_k,D_e,D_c,d_p,eps_p,tau_p,k_ldf,PCO2_clip,a1,a2,b01,b02,b1,b2,E1R,E2R,kf,Bi,c_CO2
      real(WP), parameter :: Rcst=8.314_WP     ! J/(mol.K)

      select case(trim(this%ads_model))
      case('LEE','Lee')
         ! NETL 32-D sorbent
         ! Original form taken from Lee et al. 28th International Pittsburgh Coal Conference (2011)
         ! Updated parameters taken from Suh et al. International Journal of Greenhouse Gas Control (2013)

         ! Available sites
         Cs=2531.0_WP

         ! Rate law constants for dry reaction
         A_d=10.0_WP**(2.0_WP)
         E_d=57700.0_WP
         dS_d=-174.6_WP
         dH_d=-64700.0_WP

         ! Temperature dependent parameters
         PCO2=fYCO2/W_CO2*frho*Rcst*fT
         Pg=((1.0_WP-fYCO2)/W_inert+fYCO2/W_CO2)*frho*Rcst*fT
         k_d=A_d*p%T/Cs*exp(-E_d/(Rcst*p%T))
         Keq_d=exp(dS_d/Rcst)*exp(-dH_d/(Rcst*p%T))/Pg
         
         ! Concentration in mol/m^3
         Cc=p%Mc/(Pi/6.0_WP*p%d**3*W_CO2)
         
         ! Reaction Rates in concentrations
         dcsdt=k_d*((Cs-2.0_WP*Cc)**2*PCO2-1.0_WP/Keq_d*Cc*(Cc))
         
         ! Mass transfer
         dmdt=dcsdt*Pi/6.0_WP*p%d**3*W_CO2

         ! Heat of adsorption
         H_a=64.7e3_WP ! J/mol

      case('BOS','Bos')
         ! Lewatit VP OC 1065 sorbent
         ! From Bos et al. Chemical Engineering Journal (2019)

         ! Kinetic parameters
         qs=3.4_WP ! mol CO2/kg sorbent
         mc_s=qs*this%rho*W_CO2*(Pi/6.0_WP*p%d**3) ! kg CO2/sorbent particle
         b0_t=93.0e-5_WP ! 1/Pa
         dH_t=95.3e3_WP ! J/mol
         T0_t=353.15_WP ! K
         th0_t=0.37_WP ! dimensionless
         alpha_t=0.33_WP ! dimensionless
         k0_t=3.5e-2_WP ! mol/(kg*Pa*s)
         Ea_t=15.2e3_WP ! J/mol

         ! Temperature dependence
         PCO2=fYCO2/W_CO2*frho*Rcst*fT
         b_t=b0_t*exp(dH_t/(Rcst*T0_t)*(T0_t/p%T-1))
         th_t=th0_t+alpha_t*(1-T0_t/p%T)
         k_t=k0_t*exp(-Ea_t/(Rcst*p%T))
         mc_clip=max(p%mc,0.0_WP)

         ! Mass transfer
         dmdt=this%rho*W_CO2*(Pi/6.0_WP*p%d**3)*k_t*((1-(mc_clip/mc_s)**th_t)**(1/th_t)*PCO2-mc_clip/(mc_s*b_t))

         ! Heat of adsorption
         H_a=75.0e3_WP ! J/mol

      case('Z13X','z13x')
         ! Zeolite 13X sorbent

         ! Equilibrium expression
         ! Son et al. Journal of Chemical & Engineering Data (2018)
         PCO2=fYCO2/W_CO2*frho*Rcst*fT
         PCO2_clip=max(PCO2,0.0_WP)
         a1=2.219_WP !mol/kg
         a2=3.602_WP !mol/kg
         b01=5.866e-11_WP  !1/Pa
         b02=4.436e-11_WP  !1/Pa
         E1R=5404.0_WP !K
         E2R=4197.0_WP !K

         b1=b01*exp(E1R/p%T)
         b2=b02*exp(E2R/p%T)

         mc_s=W_CO2*(Pi/6.0_WP*p%d**3)*this%rho*(a1*b1*PCO2_clip/(1+b1*PCO2_clip)+a2*b2*PCO2_clip/(1+b2*PCO2_clip))

         ! Molecular diffusivity
         T_ref=298.0_WP !K
         D_ref=1.67e-5_WP !m^2/s
         D_m=D_ref*(p%T/T_ref)**1.5_WP

         ! Knudsen diffusivity
         d_p=1e-9_WP !m, Assumed
         D_k=d_p/3.0_WP*(8.0_WP*Rcst*p%T/(Pi*W_CO2))**(1.0_WP/2.0_WP)
         
         ! Effective diffusivity
         eps_p=0.2_WP !Assumed
         tau_p=2.5_WP !Assumed
         D_e=1.0_WP/(tau_p/eps_p*(1/D_m+1/D_k))
         mc_clip=max(p%mc,0.0_WP)

         !Estimate k_ldf
         k_ldf=15.0_WP*D_e/((p%d/2)**2.0_WP)
         dmdt=k_ldf*(mc_s-mc_clip)

         ! Heat of adsorption
         H_a=37.9e3_WP !J/mol Son et al. Journal of Chemical & Engineering Data (2018)

      case('Z13XAPG','z13xapg')
         !A modfied Zeolite 13X (APG) equilibrium and bi-LDF expression
         !Equilibrium data from Wang et al. Chemical Engineering Journal (2012)
         !See model from Wilkins and Rajendran, Adsorption (2019) and Purdue, Journal of Physical Chemistry C (2018)

         c_CO2=max(fYCO2/W_CO2*frho,0.0_WP) !mol/m^3
         a1= 3.684_WP !mol/kg
         a2= 1.074_WP !mol/kg
         b01= 2.975e-6_WP !m^3/mol
         b02= 6.516e-6_WP !m^3/mol
         E1R= 27.4e3_WP !J/mol
         E2R= 34.04e3_WP !J/mol

         b1=b01*exp(E1R/(Rcst*p%T)) !m^3/mol
         b2=b02*exp(E2R/(Rcst*p%T)) !m^3/mol

         qs=(a1*b1*p%MacroCO2/(1+b1*p%MacroCO2)+a2*b2*p%MacroCO2/(1+b2*p%MacroCO2)) !mol CO2/kg sorbent
         mc_s=W_CO2*(Pi/6.0_WP*p%d**3)*this%rho*qs !kg CO2

         ! Molecular diffusivity
         T_ref=298.0_WP !K
         D_ref=1.67e-5_WP !m^2/s
         D_m=D_ref*(p%T/T_ref)**1.5_WP

         ! Knudsen diffusivity
         d_p=1.5e-6_WP !m
         D_k=d_p/3.0_WP*(8.0_WP*Rcst*p%T/(Pi*W_CO2))**(1.0_WP/2.0_WP)
         
         ! Effective diffusivity in the pores
         eps_p=0.3_WP !Assumed
         tau_p=2.0_WP !Assumed
         D_e=1.0_WP/(tau_p*(1/D_m+1/D_k))

         ! Crystal diffusivity
         D_c=6.7e-13_WP*exp(11.7e3_WP/(Rcst*p%T))

         ! Calculate Biot number
         kf=3.17e-2_WP
         Bi=p%d/2*kf/(5.0_WP*eps_p*D_e)

         ! Clip adsorbed CO2
         mc_clip=max(p%mc,0.0_WP)
         
         ! Concentration in macropores
         dcdt=15.0_WP*D_e/(p%d/2.0_WP)**2*Bi/(Bi+1.0_WP)*(c_CO2-p%MacroCO2)-this%rho/eps_p*(15.0_WP*D_c/(d_p/2.0_WP)**2*(qs-mc_clip/(W_CO2*(Pi/6.0_WP*p%d**3)*this%rho)))

         ! CO2 adsorbed in micropores
         dmdt=15.0_WP*D_c/(d_p/2.0_WP)**2*(mc_s-mc_clip)

         ! Heat of adsorption
         H_a=31.904e3_WP

         case('Z13XAPGLDF','z13xapgldf')
         !A modfied Zeolite 13X (APG) equilibrium and LDF expression
         !Equilibrium data from Wang et al. Chemical Engineering Journal (2012)
         !See model from Wilkins and Rajendran, Adsorption (2019) and Purdue, Journal of Physical Chemistry C (2018)

         PCO2=fYCO2/W_CO2*frho*Rcst*fT
         PCO2_clip=max(PCO2,0.0_WP)
         a1= 3.684 !mol/kg
         a2= 1.074 !mol/kg
         b01= 2.975e-6 !m^3/mol
         b02= 6.516e-6 !m^3/mol
         E1R= 27.4e3 !J/mol
         E2R= 34.04e3 !J/mol

         b1=b01/(Rcst*p%T)*exp(E1R/(Rcst*p%T)) !1/Pa
         b2=b02/(Rcst*p%T)*exp(E2R/(Rcst*p%T)) !1/Pa

         qs=(a1*b1*PCO2_clip/(1+b1*PCO2_clip)+a2*b2*PCO2_clip/(1+b2*PCO2_clip)) !mol CO2/kg sorbent
         mc_s=W_CO2*(Pi/6.0_WP*p%d**3)*this%rho*qs !kg CO2
         qs=qs+epsilon(1.0_WP)

         ! Molecular diffusivity
         T_ref=298.0_WP !K
         D_ref=1.67e-5_WP !m^2/s
         D_m=D_ref*(p%T/T_ref)**1.5_WP

         ! Knudsen diffusivity
         d_p=1.5e-6_WP !m
         D_k=d_p/3.0_WP*(8.0_WP*Rcst*p%T/(Pi*W_CO2))**(1.0_WP/2.0_WP)
         
         ! Effective diffusivity
         eps_p=0.3_WP 
         tau_p=2.0_WP 
         D_e=1.0_WP/(tau_p/eps_p*(1/D_m+1/D_k))
         mc_clip=max(p%mc,0.0_WP)

         !Estimate k_ldf
         k_ldf=15.0_WP*D_e/((p%d/2)**2.0_WP)*PCO2_clip/(Rcst*fT)/(qs*this%rho)
         dmdt=k_ldf*(mc_s-mc_clip)

         ! Heat of adsorption
         H_a=31.904e3_WP
         
      case('NONE','None')
         dmdt=0.0_WP
      case default
         call die('lpt_class get_rhs: Unknown ads model')
      end select
    end block mass_transfer

    ! TO DO Compute intraparticle effects

    ! Compute heat transfer
    compute_heat_transfer: block
      real(WP) :: Pr,Nu,theta
      
      !Particle dTdt from reaction      
      dTdt_r=H_a/(rhop*W_CO2*(Pi/6.0_WP*p%d**3)*this%pCp)*dmdt

      !Find mixture heat capacity
      fCp=fYCO2*fCp_CO2+(1-fYCO2)*fCp_inert
      
      Pr=fvisc/fdiff
      !Nu=(7.0_WP-10.0_WP*fVF+5.0_WP*fVF**2)*(1.0_WP+0.7_WP*Re**(0.2_WP)*Pr**(1.0_WP/3.0_WP))& ! Gunn (1978)
      !     + (1.33_WP-2.4_WP*fVF+1.2_WP*fVF**2)*Re**(0.7_WP)*Pr**(1.0_WP/3.0_WP)
      Nu=(-0.46_WP+1.77_WP*fVF+0.69_WP*fVF**2)/fVf**3+(1.37_WP-2.4_WP*fVf+1.2_WP*fVf**2)*Re**(0.7_WP)*Pr**(1.0_WP/3.0_WP) ! Sun (2015)
      theta=1.0_WP-1.6_WP*pVF*fVF-3*pVF*fVF**4*exp(-Re**0.4_WP*pVF)
      dTdt=Nu*fCp*(fT-p%T)/(3.0_WP*Pr*this%pCp*tau)*pi/(4.0_WP*theta)+dTdt_r !<- Includes Sun's bulk Nu correction factor

      
    end block compute_heat_transfer

  end subroutine get_rhs


  !> Update particle volume fraction using our current particles
  subroutine update_VF(this)
    use mathtools, only: Pi
    implicit none
    class(lpt), intent(inout) :: this
    integer :: i
    real(WP) :: Vp
    ! Reset volume fraction
    this%VF=0.0_WP
    ! Transfer particle volume
    do i=1,this%np_
       ! Skip inactive particle
       if (this%p(i)%flag.eq.1) cycle
       ! Transfer particle volume
       Vp=Pi/6.0_WP*this%p(i)%d**3
       call this%cfg%set_scalar(Sp=Vp,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%VF,bc='n')
    end do
    this%VF=this%VF/this%cfg%vol
    ! Sum at boundaries
    call this%cfg%syncsum(this%VF)
    ! Apply volume filter
    call this%filter(this%VF)
  end subroutine update_VF

  !> Compute pseudo-turbulent kinetic energy and Reynolds stresses
  !> Mehrabadi et al. (2015) "Pseudo-turbulent gas-phase velocity 
  !> fluctuations in homogeneous gas–solid flow:
  !> fixed particle assemblies and freely evolving suspensions"
  subroutine get_ptke(this,dt,Ui,Vi,Wi,visc,rho,T,diff,Y,srcU,srcV,srcW,srcT,srcY)
    use mathtools, only: Pi
    use messager, only: die
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Ui     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Vi     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: Wi     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho    !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: T      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: diff   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: Y      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcU   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcV   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcW   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcT   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout), optional :: srcY   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(:,:,:), allocatable :: slip_x,slip_y,slip_z,Re
    real(WP), dimension(:,:,:,:), allocatable :: PTRS,PTHF,PTMF,alpha
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
    integer :: i,j,k
    real(WP) :: Vp,frho,fvisc,pVF,Rep,b_par,b_perp,Nu,Pr,fVF,alpha_par,alpha_perp
    real(WP), dimension(3) :: fvel,slip
    real(WP), parameter :: a=0.523_WP,b=0.305_WP,c=0.144_WP,d=3.511_WP,e=1.801_WP,f=0.005_WP

    ! Rotation variables
    integer :: ind
    real(WP), parameter :: one_third=1.0_WP/3.0_WP
    real(WP) :: U1_dot,U2_dot,U3_dot,buf
    real(WP), dimension(3) :: U1,U2,U3
    real(WP), dimension(3,3) :: Q,temp,bij,alphaij

    ! Check for consistency
    if (present(srcY).and..not.present(srcT)) call die('[lpt get_ptke] srcY requires srcT to be present')
    if (present(srcT).and.(.not.present(T).or..not.present(diff))) call die('[lpt get_ptke] srcT requires T and diff')

    ! Allocate arrays
    allocate(slip_x(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); slip_x=0.0_WP
    allocate(slip_y(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); slip_y=0.0_WP
    allocate(slip_z(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); slip_z=0.0_WP
    allocate(Re(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); Re=0.0_WP
    allocate(PTRS(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    if (present(srcU).or.present(srcV).or.present(srcW).or.present(srcT).or.present(srcY)) then
       allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    end if
    if (present(srcT)) allocate(alpha(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    if (present(srcT)) allocate(PTHF(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    if (present(srcY)) allocate(PTMF(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

    ! Loop over all particles
    do i=1,this%np_
       ! Skip inactive particle
       if (this%p(i)%flag.eq.1.or.this%p(i)%id.eq.0) cycle

       ! Interpolate fluid quantities to particle location
       fvel(1)=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=Ui,bc='n')
       fvel(2)=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=Vi,bc='n')
       fvel(3)=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=Wi,bc='n')
       fvisc=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=visc,bc='n')
       fvisc=fvisc+epsilon(1.0_WP)
       frho=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=rho,bc='n')
       pVF=this%cfg%get_scalar(pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=this%VF,bc='n')
       pVF=pVF+epsilon(1.0_WP)

       ! Compute particle slip and Reynolds number
       slip=this%p(i)%vel-fvel
       Rep=(1.0_WP-pVF)*frho*norm2(slip)*this%p(i)%d/fvisc

       ! Transfer to the grid
       Vp=Pi/6.0_WP*this%p(i)%d**3
       call this%cfg%set_scalar(Sp=Vp*slip(1)/pVF,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=slip_x,bc='n')
       call this%cfg%set_scalar(Sp=Vp*slip(2)/pVF,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=slip_y,bc='n')
       call this%cfg%set_scalar(Sp=Vp*slip(3)/pVF,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=slip_z,bc='n')
       call this%cfg%set_scalar(Sp=Vp*Rep/pVF    ,pos=this%p(i)%pos,i0=this%p(i)%ind(1),j0=this%p(i)%ind(2),k0=this%p(i)%ind(3),S=Re    ,bc='n')
    end do
    slip_x=slip_x/this%cfg%vol
    slip_y=slip_y/this%cfg%vol
    slip_z=slip_z/this%cfg%vol
    Re=Re/this%cfg%vol
    ! Sum at boundaries
    call this%cfg%syncsum(slip_x)
    call this%cfg%syncsum(slip_y)
    call this%cfg%syncsum(slip_z)
    call this%cfg%syncsum(Re)
    ! Apply volume filter
    call this%filter(slip_x)
    call this%filter(slip_y)
    call this%filter(slip_z)
    call this%filter(Re)

    PTRS=0.0_WP
    do k=this%cfg%kmino_,this%cfg%kmaxo_
       do j=this%cfg%jmino_,this%cfg%jmaxo_
          do i=this%cfg%imino_,this%cfg%imaxo_
             ! Get Eulerian particle data (remove volume fraction)
             slip(1)=slip_x(i,j,k)
             slip(2)=slip_y(i,j,k)
             slip(3)=slip_z(i,j,k)
             Rep=Re(i,j,k)
             pVF=this%VF(i,j,k)+epsilon(1.0_WP)
             fVF=1.0_WP-pVF
             
             ! Compute PTKE
             this%ptke(i,j,k) = 0.5_WP*(sum(slip**2))*(2.0_WP*pVF + 2.5_WP*pVF*(1.0_WP-pVF)**3*exp(-pVF*sqrt(Rep)))

             !  Assume isotropic in 2D
             if (this%cfg%nx.eq.1) then
                bij = 0.0_WP
                bij(2,2) = 0.5_WP
                bij(3,3) = 0.5_WP
             else if (this%cfg%ny.eq.1) then
                bij = 0.0_WP
                bij(1,1) = 0.5_WP
                bij(3,3) = 0.5_WP
             else if (this%cfg%nz.eq.1) then
                bij = 0.0_WP
                bij(1,1) = 0.5_WP
                bij(2,2) = 0.5_WP
             else
                ! Apply rotation in 3D
                b_par  = a/(1.0_WP+b*exp(-c*Rep))*exp(-d*pVF/(1.0_WP+e*exp(-f*Rep)))
                b_perp = -0.5_WP*b_par                

                ! Add trace
                b_par  = b_par  + one_third
                b_perp = b_perp + one_third

                ! Compute the reference axes
                U1 = slip
                U1_dot = dot_product(U1,U1)+epsilon(1.0_WP)

                ! Generate an orthonormal set based on max slip component
                ind = maxloc(abs(U1),dim=1)
                select case (ind)
                case(1)
                   ! Max slip in x-direction
                   U2     = -(U1(2)/U1_dot)*U1
                   U2(2)  = 1.0_WP + U2(2)
                   U2_dot = dot_product(U2,U2)+epsilon(1.0_WP)
                case(2)
                   ! Max slip in y-direction
                   U2     = -(U1(1)/U1_dot)*U1
                   U2(1)  = 1.0_WP + U2(1)
                   U2_dot = dot_product(U2,U2)+epsilon(1.0_WP)
                case(3)
                   ! Max slip in z-direction
                   U2     = -(U1(1)/U1_dot)*U1
                   U2(1)  = 1.0_WP + U2(1)
                   U2_dot = dot_product(U2,U2)+epsilon(1.0_WP)                
                end select

                ! U3 right-hand coordinate system (cross product)
                U3(1) = U1(2)*U2(3) - U1(3)*U2(2)
                U3(2) = U1(3)*U2(1) - U1(1)*U2(3)
                U3(3) = U1(1)*U2(2) - U1(2)*U2(1)
                U3_dot = dot_product(U3,U3)+epsilon(1.0_WP)

                ! Normalize basis vectors
                U1 = U1/sqrt(U1_dot)
                U2 = U2/sqrt(U2_dot)
                U3 = U3/sqrt(U3_dot)

                ! Construct rotation matrices
                Q(:,1) = U1
                Q(:,2) = U2
                Q(:,3) = U3

                ! Multiply Q by b^dag
                temp(1,1) = Q(1,1)*b_par
                temp(2,1) = Q(2,1)*b_par
                temp(3,1) = Q(3,1)*b_par

                temp(1,2) = Q(1,2)*b_perp
                temp(2,2) = Q(2,2)*b_perp
                temp(3,2) = Q(3,2)*b_perp

                temp(1,3) = Q(1,3)*b_perp
                temp(2,3) = Q(2,3)*b_perp
                temp(3,3) = Q(3,3)*b_perp

                ! Multiply Q*b^dag (temp) by Q^T  to get bij tensor
                bij(1,1) = temp(1,1)*Q(1,1) + temp(1,2)*Q(1,2) + temp(1,3)*Q(1,3)
                bij(1,2) = temp(1,1)*Q(2,1) + temp(1,2)*Q(2,2) + temp(1,3)*Q(2,3)
                bij(1,3) = temp(1,1)*Q(3,1) + temp(1,2)*Q(3,2) + temp(1,3)*Q(3,3)

                bij(2,1) = temp(2,1)*Q(1,1) + temp(2,2)*Q(1,2) + temp(2,3)*Q(1,3)
                bij(2,2) = temp(2,1)*Q(2,1) + temp(2,2)*Q(2,2) + temp(2,3)*Q(2,3)
                bij(2,3) = temp(2,1)*Q(3,1) + temp(2,2)*Q(3,2) + temp(2,3)*Q(3,3)

                bij(3,1) = temp(3,1)*Q(1,1) + temp(3,2)*Q(1,2) + temp(3,3)*Q(1,3)
                bij(3,2) = temp(3,1)*Q(2,1) + temp(3,2)*Q(2,2) + temp(3,3)*Q(2,3)
                bij(3,3) = temp(3,1)*Q(3,1) + temp(3,2)*Q(3,2) + temp(3,3)*Q(3,3)
             end if
             
             ! Store the Reynolds stress
             buf=fVF*rho(i,j,k)*this%ptke(i,j,k)
             PTRS(1,i,j,k) = 2.0_WP*buf*bij(1,1)
             PTRS(2,i,j,k) = 2.0_WP*buf*bij(2,2)
             PTRS(3,i,j,k) = 2.0_WP*buf*bij(3,3)
             PTRS(4,i,j,k) = 2.0_WP*buf*bij(1,2)
             PTRS(5,i,j,k) = 2.0_WP*buf*bij(2,3)
             PTRS(6,i,j,k) = 2.0_WP*buf*bij(1,3)

             ! Store the pseudo-turbulent diffusivity
             if (present(T)) then
                ! Compute Nusselt Number
                Pr=visc(i,j,k)/diff(i,j,k)
                Nu=(-0.46_WP+1.77_WP*fVF+0.69_WP*fVF**2)/fVf**3+(1.37_WP-2.4_WP*fVf+1.2_WP*fVf**2)*Rep**(0.7_WP)*Pr**(1.0_WP/3.0_WP)
                ! Pseudo-turbulent diffusivity
                alpha_par=diff(i,j,k)*(2.0_WP*Rep*(Rep+1.4_WP)*Pr**2*exp(-0.002089_WP*Rep)/(3.0_WP*Pi*Nu)*&
                        (fVF*(-5.11_WP*pVF+10.1_WP*pVF**2-10.85_WP*pVF**3)+1-exp(-10.96_WP*pVF))/&
                        ((1.17_WP*pVF-0.2021_WP*pVF**(1.0_WP/2.0_WP)+0.08568_WP*pVF**(1.0_WP/4.0_WP))*fVF**2*(1.0_WP-1.6_WP*pVF*fVF-3.0_WP*pVF*fVF**4*exp(-Rep**0.4_WP*pVF))))
                this%diff_pt(i,j,k)=alpha_par/rho(i,j,k)
                !  Assume isotropic in 2D
                if (this%cfg%nx.eq.1) then
                   alphaij = 0.0_WP
                   alphaij(2,2) = alpha_par
                   alphaij(3,3) = alpha_par
                else if (this%cfg%ny.eq.1) then
                   alphaij = 0.0_WP
                   alphaij(1,1) = alpha_par
                   alphaij(3,3) = alpha_par
                else if (this%cfg%nz.eq.1) then
                   alphaij = 0.0_WP
                   alphaij(1,1) = alpha_par
                   alphaij(2,2) = alpha_par
                else
                   ! Apply rotation in 3D
                   alpha_perp=b_perp/b_par*alpha_par

                   ! Multiply Q by alpha^dag
                   temp(1,1) = Q(1,1)*alpha_par
                   temp(2,1) = Q(2,1)*alpha_par
                   temp(3,1) = Q(3,1)*alpha_par

                   temp(1,2) = Q(1,2)*alpha_perp
                   temp(2,2) = Q(2,2)*alpha_perp
                   temp(3,2) = Q(3,2)*alpha_perp

                   temp(1,3) = Q(1,3)*alpha_perp
                   temp(2,3) = Q(2,3)*alpha_perp
                   temp(3,3) = Q(3,3)*alpha_perp

                   ! Multiply Q*alpha^dag (temp) by Q^T  to get alphaij tensor
                   alphaij(1,1) = temp(1,1)*Q(1,1) + temp(1,2)*Q(1,2) + temp(1,3)*Q(1,3)
                   alphaij(1,2) = temp(1,1)*Q(2,1) + temp(1,2)*Q(2,2) + temp(1,3)*Q(2,3)
                   alphaij(1,3) = temp(1,1)*Q(3,1) + temp(1,2)*Q(3,2) + temp(1,3)*Q(3,3)

                   alphaij(2,1) = temp(2,1)*Q(1,1) + temp(2,2)*Q(1,2) + temp(2,3)*Q(1,3)
                   alphaij(2,2) = temp(2,1)*Q(2,1) + temp(2,2)*Q(2,2) + temp(2,3)*Q(2,3)
                   alphaij(2,3) = temp(2,1)*Q(3,1) + temp(2,2)*Q(3,2) + temp(2,3)*Q(3,3)

                   alphaij(3,1) = temp(3,1)*Q(1,1) + temp(3,2)*Q(1,2) + temp(3,3)*Q(1,3)
                   alphaij(3,2) = temp(3,1)*Q(2,1) + temp(3,2)*Q(2,2) + temp(3,3)*Q(2,3)
                   alphaij(3,3) = temp(3,1)*Q(3,1) + temp(3,2)*Q(3,2) + temp(3,3)*Q(3,3)
                end if
                ! Store alpha
                alpha(1,i,j,k)=alphaij(1,1)
                alpha(2,i,j,k)=alphaij(2,2)
                alpha(3,i,j,k)=alphaij(3,3)
                alpha(4,i,j,k)=alphaij(1,2)
                alpha(5,i,j,k)=alphaij(3,2)
                alpha(6,i,j,k)=alphaij(1,3)
             end if
          end do
       end do
    end do
    ! Compute PTHF
    if (present(T)) then
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                PTHF(1,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(1,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*T(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(4,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*T(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(6,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*T(i,j,k-1:k)))
                PTHF(2,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(4,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*T(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(2,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*T(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(5,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*T(i,j,k-1:k)))
                PTHF(3,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(6,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*T(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(5,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*T(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(3,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*T(i,j,k-1:k)))
             end do
          end do
       end do
       call this%cfg%sync(PTHF)
    end if
    
    if (present(Y)) then
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                PTMF(1,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(1,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*Y(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(4,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*Y(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(6,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*Y(i,j,k-1:k)))
                PTMF(2,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(4,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*Y(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(2,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*Y(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(5,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*Y(i,j,k-1:k)))
                PTMF(3,i,j,k)=fVf*(sum(this%itpr_x(:,i,j,k)*alpha(6,i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*Y(i-1:i,j,k))+&
                                   sum(this%itpr_y(:,i,j,k)*alpha(5,i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*Y(i,j-1:j,k))+&
                                   sum(this%itpr_z(:,i,j,k)*alpha(3,i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*Y(i,j,k-1:k)))
             end do
          end do
       end do
       call this%cfg%sync(PTMF)
    end if

    ! Return source terms
    if (present(srcU)) then
       ! Interpolate PTRS to cell face
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                FX(i,j,k)=sum(this%itpr_x(:,i,j,k)*PTRS(1,i-1:i,j,k))
                FY(i,j,k)=sum(this%itpr_y(:,i,j,k)*PTRS(4,i,j-1:j,k))
                FZ(i,j,k)=sum(this%itpr_z(:,i,j,k)*PTRS(6,i,j,k-1:k))
             end do
          end do
       end do
       ! Take divergence
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                srcU(i,j,k)=-dt*(sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1)))
             end do
          end do
       end do
       call this%cfg%sync(srcU)
    end if
    if (present(srcV)) then
       ! Interpolate PTRS to cell face
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                FX(i,j,k)=sum(this%itpr_x(:,i,j,k)*PTRS(4,i-1:i,j,k))
                FY(i,j,k)=sum(this%itpr_y(:,i,j,k)*PTRS(2,i,j-1:j,k))
                FZ(i,j,k)=sum(this%itpr_z(:,i,j,k)*PTRS(5,i,j,k-1:k))
             end do
          end do
       end do
       ! Take divergence
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                srcV(i,j,k)=-dt*(sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1)))
             end do
          end do
       end do
       call this%cfg%sync(srcV)
    end if
    if (present(srcW)) then
       ! Interpolate PTRS to cell face
       do k=this%cfg%kmin_,this%cfg%kmax_+1
          do j=this%cfg%jmin_,this%cfg%jmax_+1
             do i=this%cfg%imin_,this%cfg%imax_+1
                FX(i,j,k)=sum(this%itpr_x(:,i,j,k)*PTRS(6,i-1:i,j,k))
                FY(i,j,k)=sum(this%itpr_y(:,i,j,k)*PTRS(5,i,j-1:j,k))
                FZ(i,j,k)=sum(this%itpr_z(:,i,j,k)*PTRS(3,i,j,k-1:k))
             end do
          end do
       end do
       ! Take divergence
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                srcW(i,j,k)=-dt*(sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1)))
             end do
          end do
       end do
       call this%cfg%sync(srcW)
    end if
    if (present(srcT)) then
       ! Take divergence
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                srcT(i,j,k)=dt*(sum(this%div_x(:,i,j,k)*PTHF(1,i:i+1,j,k))+sum(this%div_y(:,i,j,k)*PTHF(2,i,j:j+1,k))+sum(this%div_z(:,i,j,k)*PTHF(3,i,j,k:k+1)))
             end do
          end do
       end do
       call this%cfg%sync(srcT)
    end if
    if (present(srcY)) then
       ! Take divergence
       do k=this%cfg%kmin_,this%cfg%kmax_
          do j=this%cfg%jmin_,this%cfg%jmax_
             do i=this%cfg%imin_,this%cfg%imax_
                srcY(i,j,k)=dt*(sum(this%div_x(:,i,j,k)*PTMF(1,i:i+1,j,k))+sum(this%div_y(:,i,j,k)*PTMF(2,i,j:j+1,k))+sum(this%div_z(:,i,j,k)*PTMF(3,i,j,k:k+1)))
             end do
          end do
       end do
       call this%cfg%sync(srcY)
    end if

    ! Clean up
    deallocate(slip_x,slip_y,slip_z,Re,PTRS)
    if (allocated(FX)) deallocate(FX)
    if (allocated(FY)) deallocate(FY)
    if (allocated(FZ)) deallocate(FZ)
    if (allocated(alpha)) deallocate(alpha)
    if (allocated(PTHF)) deallocate(PTHF)
    if (allocated(PTMF)) deallocate(PTMF)
  end subroutine get_ptke

  !> Stores scalar information
  subroutine scalar_init(this,sc)
    use messager, only: die
    implicit none
    class(lpt), intent(inout) :: this
    class(vdscalar), dimension(:), target, intent(in) :: sc
    character(len=str_medium) :: name
    integer :: i
    this%nscalar=size(sc)
    this%ind_T = 0
    this%ind_CO2 = 0
    do i=1,this%nscalar
       select case(trim(sc(i)%name))
       case('Temperature','T')
          this%ind_T=i
       case('CO2')
          this%ind_CO2=i
       end select
    end do
    if (this%ind_T.eq.0) call die('Scalar_init: No temperature found')
    if (this%ind_CO2.eq.0) call die('Scalar_init: No CO2 found')
  end subroutine scalar_init


  !> Laplacian filtering operation
  subroutine filter(this,A)
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP) :: filter_coeff
    integer :: i,j,k,n,nstep
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

    ! Return without filtering if filter width is zero
    if (this%filter_width.le.0.0_WP) return

    ! Recompute filter coeff and number of explicit steps needed
    filter_coeff=max(this%filter_width**2-this%cfg%min_meshsize**2,0.0_WP)/(16.0_WP*log(2.0_WP))

    if (this%implicit_filter) then  !< Apply filter implicitly
       if (.not.this%implicit%setup_done) then
          ! Prepare diffusive operator (only need to do this once)
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   this%implicit%opr(1,i,j,k)=1.0_WP-(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x(-1,i+1,j,k)+&
                   &                                  this%div_x( 0,i,j,k)*filter_coeff*this%grd_x( 0,i  ,j,k)+&
                   &                                  this%div_y(+1,i,j,k)*filter_coeff*this%grd_y(-1,i,j+1,k)+&
                   &                                  this%div_y( 0,i,j,k)*filter_coeff*this%grd_y( 0,i,j  ,k)+&
                   &                                  this%div_z(+1,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k+1)+&
                   &                                  this%div_z( 0,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k  ))
                   this%implicit%opr(2,i,j,k)=      -(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x( 0,i+1,j,k))
                   this%implicit%opr(3,i,j,k)=      -(this%div_x( 0,i,j,k)*filter_coeff*this%grd_x(-1,i  ,j,k))
                   this%implicit%opr(4,i,j,k)=      -(this%div_y(+1,i,j,k)*filter_coeff*this%grd_y( 0,i,j+1,k))
                   this%implicit%opr(5,i,j,k)=      -(this%div_y( 0,i,j,k)*filter_coeff*this%grd_y(-1,i,j  ,k))
                   this%implicit%opr(6,i,j,k)=      -(this%div_z(+1,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k+1))
                   this%implicit%opr(7,i,j,k)=      -(this%div_z( 0,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k  ))
                end do
             end do
          end do
       end if
       ! Solve the linear system
       call this%implicit%setup()
       this%implicit%rhs=A
       this%implicit%sol=0.0_WP
       call this%implicit%solve()
       A=this%implicit%sol
       call this%cfg%sync(A)
       
    else  !< Apply filter explicitly
       ! Allocate flux arrays
       allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       nstep=ceiling(6.0_WP*filter_coeff/this%cfg%min_meshsize**2)
       filter_coeff=filter_coeff/real(nstep,WP)
       do n=1,nstep
          ! Diffusive flux of A
          do k=this%cfg%kmin_,this%cfg%kmax_+1
             do j=this%cfg%jmin_,this%cfg%jmax_+1
                do i=this%cfg%imin_,this%cfg%imax_+1
                   FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
                   FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
                   FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
                end do
             end do
          end do
          ! Divergence of fluxes
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   A(i,j,k)=A(i,j,k)+sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
                end do
             end do
          end do
       end do
       ! Deallocate flux arrays
       deallocate(FX,FY,FZ)
    end if

  end subroutine filter


  !> Inject particles from a prescribed location with given mass flowrate
  !> Requires injection parameters to be set beforehand
  subroutine inject(this,dt)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    use mathtools, only: Pi
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt              !< Timestep size over which to advance
    real(WP) :: inj_min(3),inj_max(3)          !< Min/max extents of injection
    real(WP) :: Mgoal,Madded,Mtmp,buf          !< Mass flow rate parameters
    real(WP), save :: previous_error=0.0_WP    !< Store mass left over from previous timestep
    integer(kind=8) :: maxid_,maxid            !< Keep track of maximum particle id
    integer :: i,j,np0_,np2,np_tmp,count,ierr
    integer, dimension(:), allocatable :: nrecv
    type(part), dimension(:), allocatable :: p2
    type(MPI_Status) :: status
    logical :: overlap
    
    ! Initial number of particles
    np0_=this%np_
    this%np_new=0

    ! Get the particle mass that should be added to the system
    Mgoal  = this%mfr*dt+previous_error
    Madded = 0.0_WP

    ! Determine id to assign to particle
    maxid_=0
    do i=1,this%np_
       maxid_=max(maxid_,this%p(i)%id)
    end do
    call MPI_ALLREDUCE(maxid_,maxid,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr)

    ! Communicate nearby particles to check for overlap
    if (this%use_col) then
       allocate(nrecv(this%cfg%nproc))
       count=0
       inj_min(1)=this%cfg%x(this%cfg%imino)
       inj_max(1)=this%inj_pos(1)+this%inj_dmax
       inj_min(2)=this%inj_pos(2)-0.5_WP*this%inj_d-this%inj_dmax
       inj_max(2)=this%inj_pos(2)+0.5_WP*this%inj_d+this%inj_dmax
       inj_min(3)=this%inj_pos(3)-0.5_WP*this%inj_d-this%inj_dmax
       inj_max(3)=this%inj_pos(3)+0.5_WP*this%inj_d+this%inj_dmax
       do i=1,this%np_
          if ( this%p(i)%pos(1).gt.inj_min(1).and.this%p(i)%pos(1).lt.inj_max(1) .and.&
               this%p(i)%pos(2).gt.inj_min(2).and.this%p(i)%pos(2).lt.inj_max(2) .and.&
               this%p(i)%pos(3).gt.inj_min(3).and.this%p(i)%pos(3).lt.inj_max(3)) count=count+1
       end do
       call MPI_GATHER(count,1,MPI_INTEGER,nrecv,1,MPI_INTEGER,0,this%cfg%comm,ierr)
       if (this%cfg%amRoot) then
          np2=sum(nrecv)
          allocate(p2(np2))
       else
          allocate(p2(count))
       end if
       count=0
       do i=1,this%np_
          if ( this%p(i)%pos(1).gt.inj_min(1).and.this%p(i)%pos(1).lt.inj_max(1) .and.&
               this%p(i)%pos(2).gt.inj_min(2).and.this%p(i)%pos(2).lt.inj_max(2) .and.&
               this%p(i)%pos(3).gt.inj_min(3).and.this%p(i)%pos(3).lt.inj_max(3)) then
             count=count+1
             p2(count)=this%p(i)
          end if
       end do
       if (this%cfg%amRoot) then
          do i=2,this%cfg%nproc
             if (nrecv(i).gt.0) then
                call MPI_recv(p2(sum(nrecv(1:i-1))+1:sum(nrecv(1:i))),nrecv(i),MPI_PART,i-1,0,this%cfg%comm,status,ierr)
             end if
          end do
       else
          if (count.gt.0) call MPI_send(p2,count,MPI_PART,0,0,this%cfg%comm,ierr)
       end if
       deallocate(nrecv)
    end if

    ! Add new particles until desired mass is achieved
    do while (Madded.lt.Mgoal)

       if (this%cfg%amRoot) then
          ! Initialize parameters
          Mtmp = 0.0_WP
          np_tmp = 0
          ! Loop while the added volume is not sufficient
          do while (Mtmp.lt.Mgoal-Madded)
             ! Increment counter
             np_tmp=np_tmp+1
             count = np0_+np_tmp
             ! Create space for new particle
             call this%resize(count)
             ! Generate a diameter
             this%p(count)%d=get_diameter()
             ! Set various parameters for the particle
             this%p(count)%id    =maxid+int(np_tmp,8)
             this%p(count)%dt    =0.0_WP
             this%p(count)%Acol  =0.0_WP
             this%p(count)%Tcol  =0.0_WP
             this%p(count)%T     =this%inj_T
             this%p(count)%angVel=0.0_WP
             ! Give a position at the injector to the particle
             this%p(count)%pos=get_position()
             overlap=.false.
             ! Check overlap with particles recently injected
             if (this%use_col) then
                do j=1,np_tmp-1
                   if (norm2(this%p(count)%pos-this%p(j)%pos).lt.0.5_WP*(this%p(count)%d+this%p(j)%d)) overlap=.true.
                end do
                ! Check overlap with all other particles
                if (.not.overlap) then
                   do j=1,np2
                      if (norm2(this%p(count)%pos-p2(j)%pos).lt.0.5_WP*(this%p(count)%d+p2(j)%d)) overlap=.true.
                   end do
                end if
             end if

             if (overlap) then
                ! Try again
                np_tmp=np_tmp-1
             else
                ! Localize the particle
                this%p(count)%ind(1)=this%cfg%imin; this%p(count)%ind(2)=this%cfg%jmin; this%p(count)%ind(3)=this%cfg%kmin
                this%p(count)%ind=this%cfg%get_ijk_global(this%p(count)%pos,this%p(count)%ind)
                ! Give it a velocity
                this%p(count)%vel=this%inj_vel
                ! Make it an "official" particle
                this%p(count)%flag=0
                ! Update the added mass for the timestep
                Mtmp = Mtmp + this%rho*Pi/6.0_WP*this%p(count)%d**3
             end if
          end do
       end if
       ! Communicate particles
       call this%sync()
       ! Loop through newly created particles
       buf=0.0_WP
       do i=np0_+1,this%np_
          ! Remove if out of bounds
          if (this%cfg%VF(this%p(i)%ind(1),this%p(i)%ind(2),this%p(i)%ind(3)).le.0.0_WP) this%p(i)%flag=1
          if (this%p(i)%flag.eq.0) then
             ! Update the added mass for the timestep
             buf = buf + this%rho*Pi/6.0_WP*this%p(i)%d**3
             ! Update the max particle id
             maxid = max(maxid,this%p(i)%id)
             ! Increment counter
             this%np_new=this%np_new+1
          end if
       end do
       ! Total mass added
       call MPI_ALLREDUCE(buf,Mtmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); Madded=Madded+Mtmp
       ! Clean up particles
       call this%recycle()
       ! Update initial npart
       np0_=this%np_
       ! Maximum particle id
       call MPI_ALLREDUCE(maxid,maxid_,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr); maxid=maxid_
    end do

    ! Remember the error
    previous_error = Mgoal-Madded

    ! Sum up injected particles
    call MPI_ALLREDUCE(this%np_new,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_new=i

    if (allocated(p2)) deallocate(p2)

  contains

    ! Compute particle diameter
    function get_diameter() result(dp)
      use random, only: random_lognormal
      implicit none
      real(WP) :: dp
      dp=random_lognormal(m=this%inj_dmean-this%inj_dshift,sd=this%inj_dsd)+this%inj_dshift
      do while (dp.gt.this%inj_dmax+epsilon(1.0_WP).or.dp.lt.this%inj_dmin-epsilon(1.0_WP))
         dp=random_lognormal(m=this%inj_dmean-this%inj_dshift,sd=this%inj_dsd)+this%inj_dshift
      end do
    end function get_diameter

    ! Position for bulk injection of particles
    function get_position() result(pos)
      use random, only: random_uniform
      use mathtools, only: twoPi
      implicit none
      real(WP), dimension(3) :: pos
      real(WP) :: rand,r,theta
      ! Set x position
      pos(1) = this%inj_pos(1)
      ! Random y & z position within a circular region
      if (this%cfg%nz.eq.1) then
         pos(2)=random_uniform(lo=this%inj_pos(2)-0.5_WP*this%inj_d,hi=this%inj_pos(3)+0.5_WP*this%inj_d)
         pos(3) = this%cfg%zm(this%cfg%kmin_)
      else
         rand=random_uniform(lo=0.0_WP,hi=1.0_WP)
         r=0.5_WP*this%inj_d*sqrt(rand) !< sqrt(rand) avoids accumulation near the center
         call random_number(rand)
         theta=random_uniform(lo=0.0_WP,hi=twoPi)
         pos(2) = this%inj_pos(2)+r*sin(theta)
         pos(3) = this%inj_pos(3)+r*cos(theta)
      end if
    end function get_position

  end subroutine inject
  
  
    !> Calculate the CFL
  subroutine get_cfl(this,dt,cflc,cfl)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(in)  :: dt
    real(WP), intent(out) :: cflc
    real(WP), optional :: cfl
    integer :: i,j,k,ierr
    real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,my_CFL_col,my_CFLpt_x,my_CFLpt_y,my_CFLpt_z

    ! Set the CFLs to zero
    my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP; my_CFL_col=0.0_WP
    do i=1,this%np_
       my_CFLp_x=max(my_CFLp_x,abs(this%p(i)%vel(1))*this%cfg%dxi(this%p(i)%ind(1)))
       my_CFLp_y=max(my_CFLp_y,abs(this%p(i)%vel(2))*this%cfg%dyi(this%p(i)%ind(2)))
       my_CFLp_z=max(my_CFLp_z,abs(this%p(i)%vel(3))*this%cfg%dzi(this%p(i)%ind(3)))
       my_CFL_col=max(my_CFL_col,sqrt(sum(this%p(i)%vel**2))/this%p(i)%d)
    end do
    my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt

    my_CFLpt_x=0.0_WP; my_CFLpt_y=0.0_WP; my_CFLpt_z=0.0_WP
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             my_CFLpt_x=max(my_CFLpt_x,4.0_WP*this%diff_pt(i,j,k)*this%cfg%dxi(i)**2)
             my_CFLpt_y=max(my_CFLpt_y,4.0_WP*this%diff_pt(i,j,k)*this%cfg%dyi(j)**2)
             my_CFLpt_z=max(my_CFLpt_z,4.0_WP*this%diff_pt(i,j,k)*this%cfg%dzi(k)**2)
          end do
       end do
    end do
    my_CFLpt_x=my_CFLpt_x*dt; my_CFLpt_y=my_CFLpt_y*dt; my_CFLpt_z=my_CFLpt_z*dt

    ! Get the parallel max
    call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLpt_x,this%CFLpt_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLpt_y,this%CFLpt_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLpt_z,this%CFLpt_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

    ! Return the maximum convective CFL
    cflc=max(this%CFLp_x,this%CFLp_y,this%CFLp_z,this%CFLpt_x,this%CFLpt_y,this%CFLpt_z)

    ! Compute collision CFL
    my_CFL_col=10.0_WP*my_CFL_col*dt
    call MPI_ALLREDUCE(my_CFL_col,this%CFL_col,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

    ! If asked for, also return the maximum overall CFL
    if (present(CFL)) cfl=max(cflc,this%CFL_col)
    
  end subroutine get_cfl


  !> Extract various monitoring data from particle field
  subroutine get_max(this)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    class(lpt), intent(inout) :: this
    real(WP) :: buf,safe_np
    integer :: i,j,k,ierr

    ! Create safe np
    safe_np=real(max(this%np,1),WP)

    ! Particle min/max/mean
    this%Tmin=huge(1.0_WP); this%Tmax=-huge(1.0_WP); this%Tmean=0.0_WP
    this%Mmin=huge(1.0_WP); this%Mmax=-huge(1.0_WP); this%Mmean=0.0_WP
    this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
    this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
    this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
    do i=1,this%np_
       this%Tmin=min(this%Tmin,this%p(i)%T     ); this%Tmax=max(this%Tmax,this%p(i)%T     ); this%Tmean=this%Tmean+this%p(i)%T
       this%Mmin=min(this%Mmin,this%p(i)%Mc    ); this%Mmax=max(this%Mmax,this%p(i)%Mc    ); this%Mmean=this%Mmean+this%p(i)%Mc
       this%Umin=min(this%Umin,this%p(i)%vel(1)); this%Umax=max(this%Umax,this%p(i)%vel(1)); this%Umean=this%Umean+this%p(i)%vel(1)
       this%Vmin=min(this%Vmin,this%p(i)%vel(2)); this%Vmax=max(this%Vmax,this%p(i)%vel(2)); this%Vmean=this%Vmean+this%p(i)%vel(2)
       this%Wmin=min(this%Wmin,this%p(i)%vel(3)); this%Wmax=max(this%Wmax,this%p(i)%vel(3)); this%Wmean=this%Wmean+this%p(i)%vel(3)
    end do
    call MPI_ALLREDUCE(this%Tmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Tmin =buf
    call MPI_ALLREDUCE(this%Tmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Tmax =buf
    call MPI_ALLREDUCE(this%Tmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Tmean=buf/safe_np
    call MPI_ALLREDUCE(this%Mmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Mmin =buf
    call MPI_ALLREDUCE(this%Mmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Mmax =buf
    call MPI_ALLREDUCE(this%Mmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Mmean=buf/safe_np
    call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
    call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
    call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/safe_np
    call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
    call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
    call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/safe_np
    call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
    call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
    call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/safe_np

    ! Particle variance
    this%Tvar=0.0_WP
    this%Uvar=0.0_WP
    this%Vvar=0.0_WP
    this%Wvar=0.0_WP
    do i=1,this%np_
       this%Tvar=this%Tvar+(this%p(i)%T     -this%Tmean)**2.0_WP
       this%Mvar=this%Mvar+(this%p(i)%Mc    -this%Mmean)**2.0_WP
       this%Uvar=this%Uvar+(this%p(i)%vel(1)-this%Umean)**2.0_WP
       this%Vvar=this%Vvar+(this%p(i)%vel(2)-this%Vmean)**2.0_WP
       this%Wvar=this%Wvar+(this%p(i)%vel(3)-this%Wmean)**2.0_WP
    end do
    call MPI_ALLREDUCE(this%Tvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Tvar=buf/safe_np
    call MPI_ALLREDUCE(this%Mvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Mvar=buf/safe_np
    call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/safe_np
    call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/safe_np
    call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/safe_np

    ! Get mean, max, and min volume fraction
    this%VFmean=0.0_WP
    this%VFmax =-huge(1.0_WP)
    this%VFmin =+huge(1.0_WP)
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             this%VFmean=this%VFmean+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*this%VF(i,j,k)
             this%VFmax =max(this%VFmax,this%VF(i,j,k))
             this%VFmin =min(this%VFmin,this%VF(i,j,k))
          end do
       end do
    end do
    call MPI_ALLREDUCE(this%VFmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFmean=buf/this%cfg%vol_total
    call MPI_ALLREDUCE(this%VFmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%VFmax =buf
    call MPI_ALLREDUCE(this%VFmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%VFmin =buf

    ! Get volume fraction variance
    this%VFvar=0.0_WP
    do k=this%cfg%kmin_,this%cfg%kmax_
       do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
             this%VFvar=this%VFvar+this%cfg%VF(i,j,k)*this%cfg%vol(i,j,k)*(this%VF(i,j,k)-this%VFmean)**2.0_WP
          end do
       end do
    end do
    call MPI_ALLREDUCE(this%VFvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%VFvar=buf/this%cfg%vol_total

    ! Get CO2 mass mean, max and min
    this%Mcmin=huge(1.0_WP); this%Mcmax=-huge(1.0_WP); this%Mcmean=0.0_WP
    do i=1,this%np_
       this%Mcmin=min(this%Mcmin,this%p(i)%Mc); this%Mcmax=max(this%Mcmax,this%p(i)%Mc     ); this%Mcmean=this%Mcmean+this%p(i)%Mc
    end do
    call MPI_ALLREDUCE(this%Mcmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Mcmin =buf
    call MPI_ALLREDUCE(this%Mcmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Mcmax =buf
    call MPI_ALLREDUCE(this%Mcmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Mcmean=buf/safe_np
  end subroutine get_max


  !> Update particle mesh using our current particles
  subroutine update_partmesh(this,pmesh)
    use partmesh_class, only: partmesh
    implicit none
    class(lpt), intent(inout) :: this
    class(partmesh), intent(inout) :: pmesh
    integer :: i
    ! Reset particle mesh storage
    call pmesh%reset()
    ! Nothing else to do if no particle is present
    if (this%np_.eq.0) return
    ! Copy particle info
    call pmesh%set_size(this%np_)
    do i=1,this%np_
       pmesh%pos(:,i)=this%p(i)%pos
    end do
  end subroutine update_partmesh


  !> Creation of the MPI datatype for particle
  subroutine prepare_mpi_part()
    use mpi_f08
    use messager, only: die
    implicit none
    integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
    integer(MPI_ADDRESS_KIND) :: lb,extent
    type(MPI_Datatype) :: MPI_PART_TMP
    integer :: i,mysize,ierr
    ! Prepare the displacement array
    disp(1)=0
    do i=2,part_nblock
       call MPI_Type_size(part_tblock(i-1),mysize,ierr)
       disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(part_lblock(i-1),MPI_ADDRESS_KIND)
    end do
    ! Create and commit the new type
    call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
    call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
    call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
    call MPI_Type_commit(MPI_PART,ierr)
    ! If a problem was encountered, say it
    if (ierr.ne.0) call die('[lpt prepare_mpi_part] MPI Particle type creation failed')
    ! Get the size of this type
    call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
  end subroutine prepare_mpi_part


  !> Synchronize particle arrays across processors
  subroutine sync(this)
    use mpi_f08
    implicit none
    class(lpt), intent(inout) :: this
    integer, dimension(0:this%cfg%nproc-1) :: nsend_proc,nrecv_proc
    integer, dimension(0:this%cfg%nproc-1) :: nsend_disp,nrecv_disp
    integer :: n,prank,ierr
    type(part), dimension(:), allocatable :: buf_send
    ! Recycle first to minimize communication load
    call this%recycle()
    ! Prepare information about what to send
    nsend_proc=0
    do n=1,this%np_
       prank=this%cfg%get_rank(this%p(n)%ind)
       nsend_proc(prank)=nsend_proc(prank)+1
    end do
    nsend_proc(this%cfg%rank)=0
    ! Inform processors of what they will receive
    call MPI_ALLtoALL(nsend_proc,1,MPI_INTEGER,nrecv_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    ! Prepare displacements for all-to-all
    nsend_disp(0)=0
    nrecv_disp(0)=this%np_   !< Directly add particles at the end of main array
    do n=1,this%cfg%nproc-1
       nsend_disp(n)=nsend_disp(n-1)+nsend_proc(n-1)
       nrecv_disp(n)=nrecv_disp(n-1)+nrecv_proc(n-1)
    end do
    ! Allocate buffer to send particles
    allocate(buf_send(sum(nsend_proc)))
    ! Pack the particles in the send buffer
    nsend_proc=0
    do n=1,this%np_
       ! Get the rank
       prank=this%cfg%get_rank(this%p(n)%ind)
       ! Skip particles still inside
       if (prank.eq.this%cfg%rank) cycle
       ! Pack up for sending
       nsend_proc(prank)=nsend_proc(prank)+1
       buf_send(nsend_disp(prank)+nsend_proc(prank))=this%p(n)
       ! Flag particle for removal
       this%p(n)%flag=1
    end do
    ! Allocate buffer for receiving particles
    call this%resize(this%np_+sum(nrecv_proc))
    ! Perform communication
    call MPI_ALLtoALLv(buf_send,nsend_proc,nsend_disp,MPI_PART,this%p,nrecv_proc,nrecv_disp,MPI_PART,this%cfg%comm,ierr)
    ! Deallocate buffer
    deallocate(buf_send)
    ! Recycle to remove duplicate particles
    call this%recycle()
  end subroutine sync


  !> Share particles across processor boundaries
  subroutine share(this,nover)
    use mpi_f08
    use messager, only: warn,die
    implicit none
    class(lpt), intent(inout) :: this
    integer, optional :: nover
    type(part), dimension(:), allocatable :: tosend
    type(part), dimension(:), allocatable :: torecv
    integer :: no,nsend,nrecv
    type(MPI_Status) :: status
    integer :: icnt,isrc,idst,ierr
    integer :: i,n

    ! Check overlap size
    if (present(nover)) then
       no=nover
       if (no.gt.this%cfg%no) then
          call warn('[lpt_class share] Specified overlap is larger than that of cfg - reducing no')
          no=this%cfg%no
       else if (no.le.0) then
          call die('[lpt_class share] Specified overlap cannot be less or equal to zero')
       end if
    else
       no=1
    end if

    ! Clean up ghost array
    call this%resize_ghost(n=0); this%ng_=0

    ! Share ghost particles to the left in x
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(1).lt.this%cfg%imin_+no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(1).lt.this%cfg%imin_+no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%xper.and.tosend(nsend)%ind(1).lt.this%cfg%imin+no) then
             tosend(nsend)%pos(1)=tosend(nsend)%pos(1)+this%cfg%xL
             tosend(nsend)%ind(1)=tosend(nsend)%ind(1)+this%cfg%nx
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,0,-1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

    ! Share ghost particles to the right in x
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(1).gt.this%cfg%imax_-no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(1).gt.this%cfg%imax_-no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%xper.and.tosend(nsend)%ind(1).gt.this%cfg%imax-no) then
             tosend(nsend)%pos(1)=tosend(nsend)%pos(1)-this%cfg%xL
             tosend(nsend)%ind(1)=tosend(nsend)%ind(1)-this%cfg%nx
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,0,+1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

    ! Share ghost particles to the left in y
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
             tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
             tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,1,-1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

    ! Share ghost particles to the right in y
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
             tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
             tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,1,+1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

    ! Share ghost particles to the left in z
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
             tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
             tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,2,-1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

    ! Share ghost particles to the right in z
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
    end do
    allocate(tosend(nsend))
    nsend=0
    do n=1,this%np_
       if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) then
          nsend=nsend+1
          tosend(nsend)=this%p(n)
          if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
             tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
             tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
          end if
       end if
    end do
    nrecv=0
    call MPI_CART_SHIFT(this%cfg%comm,2,+1,isrc,idst,ierr)
    call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
    allocate(torecv(nrecv))
    call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
    call this%resize_ghost(this%ng_+nrecv)
    this%g(this%ng_+1:this%ng_+nrecv)=torecv
    this%ng_=this%ng_+nrecv
    if (allocated(tosend)) deallocate(tosend)
    if (allocated(torecv)) deallocate(torecv)

  end subroutine share


  !> Adaptation of particle array size
  subroutine resize(this,n)
    implicit none
    class(lpt), intent(inout) :: this
    integer, intent(in) :: n
    type(part), dimension(:), allocatable :: tmp
    integer :: size_now,size_new
    ! Resize particle array to size n
    if (.not.allocated(this%p)) then
       ! Allocate directly to size n
       allocate(this%p(n))
       this%p(1:n)%flag=1
    else
       ! Update from a non-zero size to another non-zero size
       size_now=size(this%p,dim=1)
       if (n.gt.size_now) then
          size_new=max(n,int(real(size_now,WP)*coeff_up))
          allocate(tmp(size_new))
          tmp(1:size_now)=this%p
          tmp(size_now+1:)%flag=1
          call move_alloc(tmp,this%p)
       else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
          allocate(tmp(n))
          tmp(1:n)=this%p(1:n)
          call move_alloc(tmp,this%p)
       end if
    end if
  end subroutine resize


  !> Adaptation of ghost array size
  subroutine resize_ghost(this,n)
    implicit none
    class(lpt), intent(inout) :: this
    integer, intent(in) :: n
    type(part), dimension(:), allocatable :: tmp
    integer :: size_now,size_new
    ! Resize ghost array to size n
    if (.not.allocated(this%g)) then
       ! Allocate directly to size n
       allocate(this%g(n))
       this%g(1:n)%flag=1
    else
       ! Update from a non-zero size to another non-zero size
       size_now=size(this%g,dim=1)
       if (n.gt.size_now) then
          size_new=max(n,int(real(size_now,WP)*coeff_up))
          allocate(tmp(size_new))
          tmp(1:size_now)=this%g
          tmp(size_now+1:)%flag=1
          call move_alloc(tmp,this%g)
       else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
          allocate(tmp(n))
          tmp(1:n)=this%g(1:n)
          call move_alloc(tmp,this%g)
       end if
    end if
  end subroutine resize_ghost


  !> Clean-up of particle array by removing flag=1 particles
  subroutine recycle(this)
    implicit none
    class(lpt), intent(inout) :: this
    integer :: new_size,i,ierr
    ! Compact all active particles at the beginning of the array
    new_size=0
    if (allocated(this%p)) then
       do i=1,size(this%p,dim=1)
          if (this%p(i)%flag.ne.1) then
             new_size=new_size+1
             if (i.ne.new_size) then
                this%p(new_size)=this%p(i)
                this%p(i)%flag=1
             end if
          end if
       end do
    end if
    ! Resize to new size
    call this%resize(new_size)
    ! Update number of particles
    this%np_=new_size
    call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    this%np=sum(this%np_proc)
  end subroutine recycle


  !> Parallel write particles to file
  subroutine write(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(lpt), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: i,ierr,iunit

    ! Root serial-writes the file header
    if (this%cfg%amRoot) then
       ! Open the file
       open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
       if (ierr.ne.0) call die('[lpt write] Problem encountered while serial-opening data file: '//trim(filename))
       ! Number of particles and particle object size
       write(iunit) this%np,MPI_PART_SIZE
       ! Done with the header
       close(iunit)
    end if

    ! The rest is done in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[lpt write] Problem encountered while parallel-opening data file: '//trim(filename))

    ! Get current position
    call MPI_FILE_GET_POSITION(ifile,offset,ierr)

    ! Compute the offset and write
    do i=1,this%cfg%rank
       offset=offset+int(this%np_proc(i),MPI_OFFSET_KIND)*int(MPI_PART_SIZE,MPI_OFFSET_KIND)
    end do
    if (this%np_.gt.0) call MPI_FILE_WRITE_AT(ifile,offset,this%p,this%np_,MPI_PART,status,ierr)

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Wrote ",i0," particles to file [",a,"] on partitioned grid [",a,"]")') this%np,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine write


  !> Parallel read particles to file
  subroutine read(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(lpt), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
    integer :: i,j,ierr,npadd,psize,nchunk,cnt
    integer, dimension(:,:), allocatable :: ppp

    ! First open the file in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[lpt read] Problem encountered while reading data file: '//trim(filename))

    ! Read file header first
    call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
    call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)

    ! Remember current position
    call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)

    ! Check compatibility of particle type
    if (psize.ne.MPI_PART_SIZE) call die('[lpt read] Particle type unreadable')

    ! Naively share reading task among all processors
    nchunk=int(npadd/(this%cfg%nproc*part_chunk_size))+1
    allocate(ppp(this%cfg%nproc,nchunk))
    ppp=int(npadd/(this%cfg%nproc*nchunk))
    cnt=0
    out:do j=1,nchunk
       do i=1,this%cfg%nproc
          cnt=cnt+1
          if (cnt.gt.mod(npadd,this%cfg%nproc*nchunk)) exit out
          ppp(i,j)=ppp(i,j)+1
       end do
    end do out

    ! Read by chunk
    do j=1,nchunk
       ! Find offset
       offset=header_offset+int(MPI_PART_SIZE,MPI_OFFSET_KIND)*int(sum(ppp(1:this%cfg%rank,:))+sum(ppp(this%cfg%rank+1,1:j-1)),MPI_OFFSET_KIND)
       ! Resize particle array
       call this%resize(this%np_+ppp(this%cfg%rank+1,j))
       ! Read this file
       call MPI_FILE_READ_AT(ifile,offset,this%p(this%np_+1:this%np_+ppp(this%cfg%rank+1,j)),ppp(this%cfg%rank+1,j),MPI_PART,status,ierr)
       ! Most general case: relocate every droplet
       do i=this%np_+1,this%np_+ppp(this%cfg%rank+1,j)
          this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
       end do
       ! Exchange all that
       call this%sync()
    end do

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Read ",i0," particles from file [",a,"] on partitioned grid [",a,"]")') npadd,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine read


end module lpt_class
