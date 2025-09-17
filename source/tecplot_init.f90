module tecplot_wrappers
!-----------------------------------
! tecplot_init.f90
!-----------------------------------
! Author: Juan Diego Colmenares
!-----------------------------------
!
!	This module calls on the functions of the tecio64.a library
! 	to create the TECPLOT binary file. See the Tecplot 360 Data 
! 	Format Guide for more information.
!
!---------------------------------------------------------------------
	
	USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_CHAR

    	implicit none
    	public

    	real(8) 	  :: SolTime
    	integer, private :: VIsDouble,  &
    	        	     FileType,   &
    	    		     ZoneType,   &
    	    		     ParentZn,   &
    	    		     IsBlock,    &
    	    		     Debug,      &
    	      		     ICellMax,   &
    	    		     JCellMax,   &
    	    		     KCellMax,   &
    	    		     NFConns,    &
    	    		     FNMode,     &
    	    		     ShrConn,    &
    	    		     III

 	integer :: StrandID
    	integer, dimension(:), pointer :: NullPtr
    	integer, private               :: I, k

    	include 'tecio.f90'

    	contains

    	subroutine iniTec360GridFile(name, isDouble)
	! 
	! Create a file contaning only information about the
	! grid. This may be used for visualization of the initial geometry.
	!
    	    character(len=*), intent(in) :: name
    	    integer, intent(in), optional :: isDouble
    	    integer :: FileType

    	    NullPtr => null()
    	 	Debug   = 0
    	  	FileType = 1

    	  	if (present(isDouble)) then
    	    	VIsDouble = isDouble
    	    else
    	    	VIsDouble = 0
    	    end if

    	    I = TecIni112( name//C_NULL_CHAR,        &
    	                  'X Y Z'//C_NULL_CHAR,      &
    	                  name//'.plt'//C_NULL_CHAR, &
    	                  '.'//C_NULL_CHAR, 	     &
    	    	      FileType, 	     &
    	    	      Debug, 		     &
    	    	      VIsDouble )

    	end subroutine iniTec360GridFile

    	subroutine iniTec360File(name, variables, isDouble)
	! 
	! Create a file contaning only information about the grid and cell
	! cell centered solution values.
	!
    	  	character(len=*), intent(in) :: name, variables
    	  	integer, intent(in), optional :: isDouble

    	  	NullPtr => null()
    	  	Debug   = 0
    	  	FileType = 0
    	  	if (present(isDouble)) then
    	    	VIsDouble = isDouble
    	  	else
    	    	VIsDouble = 0
    	  	end if

    	  	I = TecIni112(name//C_NULL_CHAR, &
    	                variables//C_NULL_CHAR, &
    	                name//'.plt'//C_NULL_CHAR, &
    	                '.'//C_NULL_CHAR, &
    	                FileType, &
    	                Debug, &
    	                VIsDouble)

    	end subroutine iniTec360File

    	subroutine newTec360GridZone(name, nodes, panels, isDouble)
	! 
	! Create a zone for a Tecplot file contaning only information about the
	! grid. This may be used for visualization of the geometry.
	!
    	    real(8), dimension(:,:), intent(in) :: nodes
    	    integer, dimension(:,:),  intent(in) :: panels
    	    character(len=*), intent(in) :: name
    	    integer, intent(in), optional, target :: isDouble
    	    integer :: num_nodes, num_panels
    	    integer, pointer :: ptr

    	    num_nodes = size(nodes,2)
    	    num_panels = size(panels,2)

    	    if (present(isDouble)) then
    	        ptr => isDouble
    	    else
    	        ptr = 0
    	    end if

    	    ZoneType = 3
    	    ParentZn = 0
    	    IsBlock = 1
    	    ICellMax = 0
    	    JCellMax = 0
    	    KCellMax = 0
    	    NFConns = 0
    	    FNMode = 0
    	    ShrConn = 0

    	    I = TecZne112(name//C_NULL_CHAR, &
    	                ZoneType, &
    	                num_nodes, &
    	                num_panels, &
    	                0, &
    	                ICellMax, &
    	                JCellMax, &
    	                KCellMax, &
    	                SolTime, &
    	                StrandID, &
    	                ParentZn, &
    	                IsBlock, &
    	                NFConns, &
    	                FNMode, &
    	                0, &
    	                0, &
    	                0, &
    	                nullptr, &
    	                nullptr, &
    	                nullptr, &
    	                ShrConn)
    	  	III = num_nodes
    	  	I   = TecDat112(III, nodes(1,:), ptr)
    	  	I   = TecDat112(III, nodes(2,:), ptr)
    	  	I   = TecDat112(III, nodes(3,:), ptr)
    	  	I   = TecNod112(panels)

    	end subroutine newTec360GridZone

    	subroutine newTec360Zone(name, nodes, panels, variables,&
    	    	             ValueLocation, isDouble)
    	    ! 	     
    	    ! create zone for tecplot binary file. This may be used for general
	    ! purposes and is currently not used by GUAVLAM.
    	    !
    	    real(8), dimension(:,:), intent(in) :: nodes, variables
    	    integer, dimension(:,:),  intent(in) :: panels
    	    character(len=*), intent(in) :: name
    	    integer :: num_nodes, num_panels
    	    integer, dimension(:), intent(in) :: ValueLocation
    	    integer, intent(in), optional:: isDouble

    	    num_nodes = size(nodes,2)
    	    num_panels = size(panels,2)

    	    Debug   = 0
    	    if (present(isDouble)) then
    	        VIsDouble = isDouble
    	    else
    	        VIsDouble = 0
    	    end if
    	    ZoneType = 3
    	    ParentZn = 0
    	    IsBlock = 1
    	    ICellMax = 0
    	    JCellMax = 0
    	    KCellMax = 0
    	    NFConns = 0
    	    FNMode = 0
    	    ShrConn = 0

    	    I = TecZne112(name//C_NULL_CHAR, &
    	                ZoneType, &
    	                num_nodes, &
    	                num_panels, &
    	                0, &
    	                ICellMax, &
    	                JCellMax, &
    	                KCellMax, &
    	                SolTime, &
    	                StrandID, &
    	                ParentZn, &
    	                IsBlock, &
    	                NFConns, &
    	                FNMode, &
    	                0, &
    	                0, &
    	                0, &
    	                nullptr, &
    	                ValueLocation, &
    	                nullptr, &
    	                ShrConn)
    	  	III = num_nodes
    	  	I   = TecDat112(III,nodes(1,:),1)
    	  	I   = TecDat112(III,nodes(2,:),1)
    	  	I   = TecDat112(III,nodes(3,:),1)

    	    do k = 1, size(variables,2)
    	        select case (ValueLocation(k+3))
    	        case (0)
    	            III = num_panels
    	        case (1)
    	            III = num_nodes
    	        end select
    	        I   = TecDat112(III,variables(:,k),1)
    	    end do
    	   	I   = TecNod112(panels)

    	end subroutine newTec360Zone


    	subroutine newSurfaceZone(name, nodes, normals, loads, panels,&
    	    	             ValueLocation, isDouble)
    	    ! 	     
    	    ! create zone for tecplot binary file containing the nodes of
	    ! surface, as well as the normal vector components and the loads the
	    ! bound vortex sheets.
    	    !
    	    real(8), dimension(:,:), intent(in) :: nodes, normals, loads
    	    integer, dimension(:,:),  intent(in) :: panels
    	    character(len=*), intent(in) :: name
    	    integer :: num_nodes, num_panels
    	    integer, dimension(:), intent(in) :: ValueLocation
    	    integer, intent(in), optional:: isDouble

    	    num_nodes = size(nodes,2)
    	    num_panels = size(panels,2)

    	    Debug   = 0
    	    if (present(isDouble)) then
    	        VIsDouble = isDouble
    	    else
    	        VIsDouble = 0
    	    end if
    	    ZoneType = 3
    	    ParentZn = 0
    	    IsBlock = 1
    	    ICellMax = 0
    	    JCellMax = 0
    	    KCellMax = 0
    	    NFConns = 0
    	    FNMode = 0
    	    ShrConn = 0

    	    I = TecZne112(name//C_NULL_CHAR, &
    	                ZoneType, &
    	                num_nodes, &
    	                num_panels, &
    	                0, &
    	                ICellMax, &
    	                JCellMax, &
    	                KCellMax, &
    	                SolTime, &
    	                StrandID, &
    	                ParentZn, &
    	                IsBlock, &
    	                NFConns, &
    	                FNMode, &
    	                0, &
    	                0, &
    	                0, &
    	                nullptr, &
    	                ValueLocation, &
    	                nullptr, &
    	                ShrConn)
    	  	III = num_nodes
    	  	I   = TecDat112(III,nodes(1,:),1)
    	  	I   = TecDat112(III,nodes(2,:),1)
    	  	I   = TecDat112(III,nodes(3,:),1)

    	        III = num_panels
    	        I   = TecDat112(III,normals(1,:),1)
    	        I   = TecDat112(III,normals(2,:),1)
    	        I   = TecDat112(III,normals(3,:),1)
    	        I   = TecDat112(III,loads(:,1),1)

    	   	I   = TecNod112(panels)

    	end subroutine newSurfaceZone

    	subroutine newWakeZone(name, nodes, panels, ValueLocation, isDouble)
    	    ! 	     
    	    ! create zone for tecplot binary file containing the nodes of a wake
	    ! structure.
    	    !
    	    real(8), dimension(:,:), intent(in) :: nodes
    	    integer, dimension(:,:),  intent(in) :: panels
    	    character(len=*), intent(in) :: name
    	    integer :: num_nodes, num_panels
    	    integer, dimension(:), intent(in) :: ValueLocation
    	    integer, intent(in), optional:: isDouble

    	    num_nodes = size(nodes,2)
    	    num_panels = size(panels,2)

    	    Debug   = 0
    	    if (present(isDouble)) then
    	        VIsDouble = isDouble
    	    else
    	        VIsDouble = 0
    	    end if
    	    ZoneType = 3
    	    ParentZn = 0
    	    IsBlock = 1
    	    ICellMax = 0
    	    JCellMax = 0
    	    KCellMax = 0
    	    NFConns = 0
    	    FNMode = 0
    	    ShrConn = 0

    	    I = TecZne112(name//C_NULL_CHAR, &
    	                ZoneType, &
    	                num_nodes, &
    	                num_panels, &
    	                0, &
    	                ICellMax, &
    	                JCellMax, &
    	                KCellMax, &
    	                SolTime, &
    	                StrandID, &
    	                ParentZn, &
    	                IsBlock, &
    	                NFConns, &
    	                FNMode, &
    	                0, &
    	                0, &
    	                0, &
    	                [0,0,0,1,1,1,1], &
    	                ValueLocation, &
    	                nullptr, &
    	                ShrConn)
    	  	III = num_nodes
    	  	I   = TecDat112(III,nodes(1,:),1)
    	  	I   = TecDat112(III,nodes(2,:),1)
    	  	I   = TecDat112(III,nodes(3,:),1)

    	   	I   = TecNod112(panels)

    	end subroutine newWakeZone

end module tecplot_wrappers
