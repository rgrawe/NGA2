!> Incompressible scalar solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes constant diffusivity and density.
module scalar_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: scalar,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available advection schemes for scalar transport
   integer, parameter, public :: upwind=0            !< First order upwind scheme
   integer, parameter, public :: quick=1             !< Quick scheme
   integer, parameter, public :: bquick=2            !< BQuick scheme
   
   
   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Constant density scalar solver object definition
   type :: scalar
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_SCALAR'  !< Solver name (default=UNNAMED_SCALAR)
      
      ! Constant property fluid, but diffusivity is still a field due to LES modeling
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP), dimension(:,:,:), allocatable :: diff     !< These is our constant+SGS dynamic diffusivity for the scalar
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Scalar variable
      real(WP), dimension(:,:,:), allocatable :: SC       !< SC array
      
      ! Old scalar variable
      real(WP), dimension(:,:,:), allocatable :: SCold    !< SCold array
      
      ! Implicit scalar solver
      class(linsol), pointer :: implicit                  !< Iterative linear solver object for an implicit prediction of the scalar residual
      integer, dimension(:,:,:), allocatable :: stmap     !< Inverse map from stencil shift to index location
      
      ! Metrics
      integer :: scheme                                   !< Advection scheme for scalar
      integer :: nst                                      !< Scheme order (and elemental stencil size)
      integer :: stp1,stp2                                !< Plus interpolation stencil extent for scalar advection
      integer :: stm1,stm2                                !< Minus interpolation stencil extent for scalar advection
      real(WP), dimension(:,:,:,:), allocatable :: itp_xp,itp_yp,itp_zp  !< Plus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: itp_xm,itp_ym,itp_zm  !< Minus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: div_x ,div_y ,div_z   !< Divergence for SC
      real(WP), dimension(:,:,:,:), allocatable :: grd_x ,grd_y ,grd_z   !< Scalar gradient for SC
      real(WP), dimension(:,:,:,:), allocatable :: itp_x ,itp_y ,itp_z   !< Second order interpolation for SC diffusivity
      
      ! Bquick requires additional storage
	   real(WP), dimension(:,:,:,:), allocatable :: bitp_xp,bitp_yp,bitp_zp  !< Plus interpolation for SC  - backup
      real(WP), dimension(:,:,:,:), allocatable :: bitp_xm,bitp_ym,bitp_zm  !< Minus interpolation for SC - backup

      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for modifying SC metrics
      
      ! Monitoring quantities
      real(WP) :: SCmax,SCmin,SCint                       !< Maximum and minimum, integral scalar
      
   contains
      procedure :: print=>scalar_print                    !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the scalar solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: metric_reset                           !< Reset adaptive metrics like bquick
      procedure :: metric_adjust                          !< Adjust adaptive metrics like bquick
      procedure :: get_drhoSCdt                           !< Calculate drhoSC/dt
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: get_int                                !< Calculate integral field values
      procedure :: solve_implicit                         !< Solve for the scalar residuals implicitly
   end type scalar
   
   
   !> Declare scalar solver constructor
   interface scalar
      procedure constructor
   end interface scalar
   
contains
   
   
   !> Default constructor for scalar solver
   function constructor(cfg,scheme,name) result(self)
      use messager, only: die
      implicit none
      type(scalar) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: scheme
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%SC   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SC   =0.0_WP
      allocate(self%SCold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SCold=0.0_WP
      allocate(self%diff (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%diff =0.0_WP
      
      ! Prepare advection scheme
      self%scheme=scheme
      select case (self%scheme)
      case (upwind)
         ! Check current overlap
		   if (self%cfg%no.lt.1) call die('[scalar constructor] Scalar transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
		   self%nst=1
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case (quick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[scalar constructor] Scalar transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case (bquick)
         ! Check current overlap
	      if (self%cfg%no.lt.2) call die('[scalar constructor] Scalar transport scheme requires larger overlap')
	      ! Set interpolation stencil sizes
	      self%nst=3
	      self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
	      self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case default
         call die('[scalar constructor] Unknown scalar transport scheme selected')
      end select
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for SC
      allocate(self%mask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(self%cfg%imax+1:,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(:,self%cfg%jmax+1:,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(:,:,self%cfg%kmax+1:)=2
      end if
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask)
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(scalar), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference diffusivity interpolation coefficients
      allocate(this%itp_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create diffusivity interpolation coefficients to cell face
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%itp_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate finite difference scalar interpolation coefficients
      allocate(this%itp_xp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_xm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_yp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_ym(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_zp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%itp_zm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create scalar interpolation coefficients to cell faces
      select case (this%scheme)
      case (upwind)
         this%itp_xp=1.0_WP; this%itp_xm=1.0_WP
         this%itp_yp=1.0_WP; this%itp_ym=1.0_WP
         this%itp_zp=1.0_WP; this%itp_zm=1.0_WP
      case (quick,bquick)
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itp_xp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itp_xm(:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itp_yp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itp_ym(:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itp_zp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itp_zm(:,i,j,k))
               end do
            end do
         end do
      end select
      
      ! Allocate finite volume divergence operators
      allocate(this%div_x(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_y(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_z(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%div_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%div_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grd_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%grd_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%grd_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%grd_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in x from [xm,ym,zm] to [x,ym,zm]
               this%grd_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in y from [xm,ym,zm] to [xm,y,zm]
               this%grd_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(scalar), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask)
      
      ! Adjust interpolation coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Linear interpolation in x
               if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).gt.0) this%itp_x(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i-1,j,k).eq.0) this%itp_x(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in y
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).gt.0) this%itp_y(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j-1,k).eq.0) this%itp_y(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in z
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).gt.0) this%itp_z(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j,k-1).eq.0) this%itp_z(:,i,j,k)=[1.0_WP,0.0_WP]
            end do
         end do
      end do
      
      ! Adjust scalar interpolation to reflect Dirichlet boundaries
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face
               if (this%mask(i-1,j,k).eq.2) then
                  this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm(-1,i,j,k)=1.0_WP
                  this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i  ,j,k).eq.2) then
                  this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm( 0,i,j,k)=1.0_WP
                  this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp( 0,i,j,k)=1.0_WP
               end if
               ! Y face
               if (this%mask(i,j-1,k).eq.2) then
                  this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym(-1,i,j,k)=1.0_WP
                  this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j  ,k).eq.2) then
                  this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym( 0,i,j,k)=1.0_WP
                  this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp( 0,i,j,k)=1.0_WP
               end if
               ! Z face
               if (this%mask(i,j,k-1).eq.2) then
                  this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm(-1,i,j,k)=1.0_WP
                  this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j,k  ).eq.2) then
                  this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm( 0,i,j,k)=1.0_WP
                  this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp( 0,i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to SC divergence
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).gt.0) then
                  this%div_x(:,i,j,k)=0.0_WP
                  this%div_y(:,i,j,k)=0.0_WP
                  this%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) this%grd_x(:,i,j,k)=0.0_WP     !< FD gradient in x of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) this%grd_y(:,i,j,k)=0.0_WP     !< FD gradient in y of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) this%grd_z(:,i,j,k)=0.0_WP     !< FD gradient in z of SC
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%div_x=0.0_WP
         this%grd_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%div_y=0.0_WP
         this%grd_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%div_z=0.0_WP
         this%grd_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the scalar solver now that bconds have been defined
   subroutine setup(this,implicit_solver)
      implicit none
      class(scalar), intent(inout) :: this
      class(linsol), target, intent(in), optional :: implicit_solver
      integer :: count,st
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
      ! Bquick needs to remember the quick coefficients
	   select case (this%scheme)
      case (bquick)
         allocate(this%bitp_xp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_xp=this%itp_xp
         allocate(this%bitp_xm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_xm=this%itp_xm
         allocate(this%bitp_yp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_yp=this%itp_yp
         allocate(this%bitp_ym(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_ym=this%itp_ym
         allocate(this%bitp_zp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_zp=this%itp_zp
         allocate(this%bitp_zm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)); this%bitp_zm=this%itp_zm
      end select
      
      ! Prepare implicit solver if it had been provided
      if (present(implicit_solver)) then
         
         ! Point to implicit solver linsol object
         this%implicit=>implicit_solver
         
         ! Set dynamic stencil map for the scalar solver
         count=1; this%implicit%stc(count,:)=[0,0,0]
         do st=1,abs(this%stp1)
            count=count+1; this%implicit%stc(count,:)=[+st,0,0]
            count=count+1; this%implicit%stc(count,:)=[-st,0,0]
            count=count+1; this%implicit%stc(count,:)=[0,+st,0]
            count=count+1; this%implicit%stc(count,:)=[0,-st,0]
            count=count+1; this%implicit%stc(count,:)=[0,0,+st]
            count=count+1; this%implicit%stc(count,:)=[0,0,-st]
         end do
         
         ! Set the diagonal to 1 to make sure all cells participate in solver
         this%implicit%opr(1,:,:,:)=1.0_WP
         
         ! Initialize the implicit velocity solver
         call this%implicit%init()
         
      end if
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(scalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      external :: locator
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      character(len=2), optional :: dir
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      if (present(dir)) then
         select case (lowercase(dir))
         case ('+x','x+','xp','px'); new_bc%dir=1
         case ('-x','x-','xm','mx'); new_bc%dir=2
         case ('+y','y+','yp','py'); new_bc%dir=3
         case ('-y','y-','ym','my'); new_bc%dir=4
         case ('+z','z+','zp','pz'); new_bc%dir=5
         case ('-z','z-','zm','mz'); new_bc%dir=6
         case default; call die('[scalar add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[scalar apply_bcond] Neumann requires a direction')
         new_bc%dir=0
      end if
      new_bc%itr=iterator(this%cfg,new_bc%name,locator,'c')
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         do n=1,new_bc%itr%n_
            i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
            this%mask(i,j,k)=2
         end do
      case (neumann)
         ! No modification - this assumes Neumann is only applied at walls or domain boundaries
      case default
         call die('[scalar apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(scalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[scalar get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)           ! Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann)             ! Apply Neumann condition
               
               ! Implement based on bcond direction
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%SC(i,j,k)=this%SC(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[scalar apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%SC)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit rhoSC time derivative based on rhoU/rhoV/rhoW
   subroutine get_drhoSCdt(this,drhoSCdt,rhoU,rhoV,rhoW)
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoSCdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoU     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoV     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoW     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Flux of rhoSC
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itp_xp(:,i,j,k)*this%SC(i+this%stp1:i+this%stp2,j,k)) &
               &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itp_xm(:,i,j,k)*this%SC(i+this%stm1:i+this%stm2,j,k)) &
               &         +sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*this%SC(i-1:i,j,k))
               ! Fluxes on y-face
               FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itp_yp(:,i,j,k)*this%SC(i,j+this%stp1:j+this%stp2,k)) &
               &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itp_ym(:,i,j,k)*this%SC(i,j+this%stm1:j+this%stm2,k)) &
               &         +sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*this%SC(i,j-1:j,k))
               ! Fluxes on z-face
               FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itp_zp(:,i,j,k)*this%SC(i,j,k+this%stp1:k+this%stp2)) &
               &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itp_zm(:,i,j,k)*this%SC(i,j,k+this%stm1:k+this%stm2)) &
               &         +sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*this%SC(i,j,k-1:k))
            end do
         end do
      end do
      ! Time derivative of rhoSC
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoSCdt(i,j,k)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &               sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &               sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync residual
      call this%cfg%sync(drhoSCdt)
   end subroutine get_drhoSCdt
   
   
   !> Calculate the min and max of our SC field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(scalar), intent(inout) :: this
      integer :: ierr
      real(WP) :: my_SCmax,my_SCmin
      my_SCmax=maxval(this%SC); call MPI_ALLREDUCE(my_SCmax,this%SCmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_SCmin=minval(this%SC); call MPI_ALLREDUCE(my_SCmin,this%SCmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max
   
   
   !> Calculate the integral of our SC field
   subroutine get_int(this)
      implicit none
      class(scalar), intent(inout) :: this
      call this%cfg%integrate(this%SC,integral=this%SCint)
   end subroutine get_int
   
   
   !> Solve for implicit scalar residual
   subroutine solve_implicit(this,dt,resSC,rhoU,rhoV,rhoW)
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resSC !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoU  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoV  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoW  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,sti,std
      
      ! Prepare convective operator
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)+abs(rhoU(i+std,j,k)))*this%itp_xp(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)+abs(rhoV(i,j+std,k)))*this%itp_yp(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)+abs(rhoW(i,j,k+std)))*this%itp_zp(sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)-abs(rhoU(i+std,j,k)))*this%itp_xm(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)-abs(rhoV(i,j+std,k)))*this%itp_ym(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)-abs(rhoW(i,j,k+std)))*this%itp_zm(sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
      
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x(-1,i+1,j,k)+&
               &                                                                this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x( 0,i  ,j,k)+&
               &                                                                this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y(-1,i,j+1,k)+&
               &                                                                this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y( 0,i,j  ,k)+&
               &                                                                this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z(-1,i,j,k+1)+&
               &                                                                this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z(-1,i,j,k  ))
            end do
         end do
      end do
      
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resSC
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resSC=this%implicit%sol
      
      ! Sync up residual
      call this%cfg%sync(resSC)
      
   end subroutine solve_implicit
   
   
   !> Metric resetting for adaptive discretization like bquick
   subroutine metric_reset(this)
	  implicit none
	  class(scalar), intent(inout) :: this
	  select case (this%scheme)
	  case (bquick)
	     this%itp_xp=this%bitp_xp
	     this%itp_xm=this%bitp_xm
	     this%itp_yp=this%bitp_yp
	     this%itp_ym=this%bitp_ym
	     this%itp_zp=this%bitp_zp
	     this%itp_zm=this%bitp_zm
	  end select
   end subroutine metric_reset
   
   
   !> Adjust adaptive metrics like bquick
   subroutine metric_adjust(this,SC,SCmin,SCmax)
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SC !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), optional :: SCmin
      real(WP), optional :: SCmax
      integer :: i,j,k
      select case (this%scheme)
      case (bquick)
         if (present(SCmin)) then
            do k=this%cfg%kmin_,this%cfg%kmax_+1
               do j=this%cfg%jmin_,this%cfg%jmax_+1
                  do i=this%cfg%imin_,this%cfg%imax_+1
                     if (SC(i,j,k).lt.SCmin) then
                        !this%itp_xp(i:i+1,j,k)=[0.0_WP,1.0_WP,0.0_WP]
                        !this%itp_xm(i:i+1,j,k)=[0.0_WP,1.0_WP,0.0_WP]
                        !this%itp_yp(i,j:j+1,k)=[0.0_WP,1.0_WP,0.0_WP]
                        !this%itp_ym(i,j:j+1,k)=[0.0_WP,1.0_WP,0.0_WP]
                        !this%itp_zp(i,j,k:k+1)=[0.0_WP,1.0_WP,0.0_WP]
                        !this%itp_zm(i,j,k:k+1)=[0.0_WP,1.0_WP,0.0_WP]
					      end if
                  end do
               end do
            end do
         end if
         if (present(SCmax)) then
            do k=this%cfg%kmin_,this%cfg%kmax_+1
               do j=this%cfg%jmin_,this%cfg%jmax_+1
                  do i=this%cfg%imin_,this%cfg%imax_+1
                     if (SC(i,j,k).gt.SCmax) then
                        !this%itp_xp(i:i+1,j,k)=[0.0_WP,1.0_WP,0.0_WP]
					         !this%itp_xm(i:i+1,j,k)=[0.0_WP,1.0_WP,0.0_WP]
					         !this%itp_yp(i,j:j+1,k)=[0.0_WP,1.0_WP,0.0_WP]
					         !this%itp_ym(i,j:j+1,k)=[0.0_WP,1.0_WP,0.0_WP]
					         !this%itp_zp(i,j,k:k+1)=[0.0_WP,1.0_WP,0.0_WP]
					         !this%itp_zm(i,j,k:k+1)=[0.0_WP,1.0_WP,0.0_WP]
				         end if
                  end do
               end do
            end do
         end if
      end select
   end subroutine metric_adjust
   
   
   !> Print out info for scalar solver
   subroutine scalar_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(scalar), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Constant density scalar solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine scalar_print
   
   
end module scalar_class
