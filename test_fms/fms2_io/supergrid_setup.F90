!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!>@author Jessica Liptak
!>@description This module contains routines required for test_io_supergrid
module supergrid_setup
use mpp_domains_mod, only : mpp_domains_init, mpp_define_domains, mpp_define_io_domain, mpp_define_layout
use mpp_domains_mod, only : mpp_get_domain_extents, domain2D
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain,mpp_get_global_domain
use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE, NORTH, EAST, CORNER, CENTER
use mpp_domains_mod, only : mpp_get_io_domain
use mpp_mod, only : mpp_npes
use fms2_io_mod, only : file_exists

implicit none

!> The MOM_domain_type contains information about the domain decompositoin.
type, public :: MOM_domain_type
  type(domain2D), pointer :: mpp_domain => NULL() !< The FMS domain with halos
                                !! on this processor, centered at h points.
  integer :: niglobal           !< The total horizontal i-domain size.
  integer :: njglobal           !< The total horizontal j-domain size.
  integer :: nihalo             !< The i-halo size in memory.
  integer :: njhalo             !< The j-halo size in memory.
  logical :: symmetric          !< True if symmetric memory is used with
                                !! this domain.
  integer :: layout(2)          !< This domain's processor layout.  This is
                                !! saved to enable the construction of related
                                !! new domains with different resolutions or
                                !! other properties.
  integer :: io_layout(2)       !< The IO-layout used with this domain.
  integer :: X_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the i-direction in a define_domain call.
  integer :: Y_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the j-direction in a define_domain call.
  logical :: thin_halo_updates  !< If true, optional arguments may be used to
                                !! specify the width of the halos that are
                                !! updated with each call.
  logical, pointer :: maskmap(:,:) => NULL() !< A pointer to an array indicating
                                !! which logical processors are actually used for
                                !! the ocean code. The other logical processors
                                !! would be contain only land points and are not
                                !! assigned to actual processors. This need not be
                                !! assigned if all logical processors are used.
end type MOM_domain_type

!> Container for horizontal index ranges for data, computational and global domains
type, public :: hor_index_type
  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  logical :: symmetric  !< True if symmetric memory is used.

  integer :: niglobal !< The global number of h-cells in the i-direction
  integer :: njglobal !< The global number of h-cells in the j-direction

  integer :: turns      !< Number of quarter-turn rotations from input to model
end type hor_index_type

public :: domains_init, hor_index_init, create_supergrid_domain

contains

!> Define a domain for the supergrid (SGdom)
! For a model grid size (ni, nj), the supergrid size is (2*ni+, 2*nj+1)
! For example, a 3x2 model grid defined by 6 centroids (c) has a
 ! 7x5 supergrid defined by vertices (*) and face midpoints (.)
!  *--.--*--.--*--.--*
!  |     |     |     |
!  .  c  .  c  .  c  .
!  |     |     |     |
!  *--.--*--.--*--.--*
!  |     |     |     |
!  .  c  .  c  .  c  .
!  |     |     |     |
!  *--.--*--.--*--.--*
subroutine create_supergrid_domain(G, SGdom)
  type(MOM_domain_type), intent(in) :: G !< Supergrid domain
  type(MOM_domain_type), intent(inout) :: SGdom !< supergrid domain
  ! local
  integer :: npei, npej ! number of pes in the i, j directions
  integer, dimension(:), allocatable :: exni, exnj ! domain extents
  integer :: global_indices(4), i
  integer, dimension(2) :: io_layout ! I/O layout.

  npei = G%layout(1) ; npej = G%layout(2)
  allocate(exni(npei)) ; allocate(exnj(npej))
  call mpp_get_domain_extents(G%mpp_domain, exni, exnj)
  if (.not.(associated(SGdom%mpp_domain))) allocate(SGdom%mpp_domain)
  SGdom%nihalo = 2*G%nihalo+1
  SGdom%njhalo = 2*G%njhalo+1
  SGdom%niglobal = 2*G%niglobal
  SGdom%njglobal = 2*G%njglobal
  SGdom%layout(:) = G%layout(:)
  SGdom%io_layout(:) = G%io_layout(:)
  ! define extent of the global domain
  global_indices(1) = 1+SGdom%nihalo
  global_indices(2) = SGdom%niglobal+SGdom%nihalo
  global_indices(3) = 1+SGdom%njhalo
  global_indices(4) = SGdom%njglobal+SGdom%njhalo
  exni(:) = 2*exni(:) ; exnj(:) = 2*exnj(:)
  call mpp_define_domains(global_indices, SGdom%layout, SGdom%mpp_domain, &
                          xflags=G%X_FLAGS, yflags=G%Y_FLAGS, &
                          xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
                          xextent=exni,yextent=exnj, symmetry=.true., &
                          name="MOM_supergrid")

  !do i=1,npei
 !   write(*, '(A,I2)'), "exni = ", exni(i)
 ! enddo

  !do i=1,npej
  !  write(*, '(A,I2)'), "exnj = ", exnj(i)
  !enddo

  !write(*, '(A,I2,A,I2)'), "SGdom%nihalo, SGdom%njhalo= ", 2*G%nihalo+1," ", 2*G%njhalo+1
  !write(*, '(A,I2,A,I2)'), "SGdom%niglobal, SGdom%njglobal= ", 2*G%niglobal," ", 2*G%njglobal
  !write(*, '(A,I2,A,I2,A,I2,A,I2)'), "global_indices ", global_indices(1)," ", global_indices(2), " ", global_indices(3), " ", global_indices(4)
  ! define the IO domain
  call mpp_define_io_domain(SGdom%mpp_domain, SGdom%io_layout)
  if (allocated(exni)) deallocate(exni)
  if (allocated(exnj)) deallocate(exnj)
end subroutine create_supergrid_domain

!> Returns various data that has been stored in a MOM_domain_type
subroutine get_domain_extent(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
  isg, ieg, jsg, jeg, idg_offset, jdg_offset, &
  symmetric, local_indexing, index_offset)
  type(MOM_domain_type), &
    intent(in)  :: Domain !< The MOM domain from which to extract information
    integer, intent(out) :: isc    !< The start i-index of the computational domain
    integer, intent(out) :: iec    !< The end i-index of the computational domain
    integer, intent(out) :: jsc    !< The start j-index of the computational domain
    integer, intent(out) :: jec    !< The end j-index of the computational domain
    integer, intent(out) :: isd    !< The start i-index of the data domain
    integer, intent(out) :: ied    !< The end i-index of the data domain
    integer, intent(out) :: jsd    !< The start j-index of the data domain
    integer, intent(out) :: jed    !< The end j-index of the data domain
    integer, intent(out) :: isg    !< The start i-index of the global domain
    integer, intent(out) :: ieg    !< The end i-index of the global domain
    integer, intent(out) :: jsg    !< The start j-index of the global domain
    integer, intent(out) :: jeg    !< The end j-index of the global domain
    integer, intent(out) :: idg_offset !< The offset between the corresponding global and data i-index spaces.
    integer, intent(out) :: jdg_offset !< The offset between the corresponding global and data j-index spaces.
    logical, intent(out) :: symmetric  !< True if symmetric memory is used.
    logical, optional, intent(in)  :: local_indexing !< If true, local tracer array indices start at 1 as in most MOM6 code.
    integer, optional, intent(in)  :: index_offset   !< A fixed additional offset to all indices. This
                !! can be useful for some types of debugging with
                !! dynamic memory allocation.
    ! Local variables
    integer :: ind_off
    logical :: local

  local = .true. ; if (present(local_indexing)) local = local_indexing
  ind_off = 0 ; if (present(index_offset)) ind_off = index_offset

  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(Domain%mpp_domain, isg, ieg, jsg, jeg)

  ! This code institutes the MOM convention that local array indices start at 1.
  if (local) then
    idg_offset = isd-1 ; jdg_offset = jsd-1
    isc = isc-isd+1 ; iec = iec-isd+1 ; jsc = jsc-jsd+1 ; jec = jec-jsd+1
    ied = ied-isd+1 ; jed = jed-jsd+1
    isd = 1 ; jsd = 1
  else
    idg_offset = 0 ; jdg_offset = 0
  endif
  if (ind_off /= 0) then
    idg_offset = idg_offset + ind_off ; jdg_offset = jdg_offset + ind_off
    isc = isc + ind_off ; iec = iec + ind_off
    jsc = jsc + ind_off ; jec = jec + ind_off
    isd = isd + ind_off ; ied = ied + ind_off
    jsd = jsd + ind_off ; jed = jed + ind_off
  endif
  symmetric = Domain%symmetric

end subroutine get_domain_extent

!> Sets various index values in a hor_index_type.
subroutine hor_index_init(Domain, HI, local_indexing, index_offset)
  type(MOM_domain_type),  intent(in)    :: Domain     !< The MOM domain from which to extract information.
  type(hor_index_type),   intent(inout) :: HI         !< A horizontal index type to populate with data
  logical, optional,      intent(in)    :: local_indexing !< If true, all tracer data domains start at 1
  integer, optional,      intent(in)    :: index_offset   !< A fixed additional offset to all indices

  ! get_domain_extent ensures that domains start at 1 for compatibility between
  ! static and dynamically allocated arrays.

  HI%symmetric=.true.
  call get_domain_extent(Domain, HI%isc, HI%iec, HI%jsc, HI%jec, &
                         HI%isd, HI%ied, HI%jsd, HI%jed, &
                         HI%isg, HI%ieg, HI%jsg, HI%jeg, &
                         HI%idg_offset, HI%jdg_offset, HI%symmetric, &
                         local_indexing=local_indexing)
  HI%niglobal = Domain%niglobal
  HI%njglobal = Domain%njglobal

  HI%IscB = HI%isc ; HI%JscB = HI%jsc
  HI%IsdB = HI%isd ; HI%JsdB = HI%jsd
  HI%IsgB = HI%isg ; HI%JsgB = HI%jsg
  if (HI%symmetric) then
    HI%IscB = HI%isc-1 ; HI%JscB = HI%jsc-1
    HI%IsdB = HI%isd-1 ; HI%JsdB = HI%jsd-1
    HI%IsgB = HI%isg-1 ; HI%JsgB = HI%jsg-1
  endif
  HI%IecB = HI%iec ; HI%JecB = HI%jec
  HI%IedB = HI%ied ; HI%JedB = HI%jed
  HI%IegB = HI%ieg ; HI%JegB = HI%jeg

  HI%turns = 0
end subroutine hor_index_init

!> Initalizes a MOM_domain_type variable and optionally returns data describing
!! properties of the domain type.
subroutine domains_init(MOM_dom, symmetric, static_memory, user_def_io_layout, &
  NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, X_FLAGS, Y_FLAGS, &
  min_halo)
  type(MOM_domain_type),           intent(inout) :: MOM_dom !< Structure with the domain information
  logical, optional,               intent(in)    :: symmetric    !< If present, this specifies
  !! whether this domain is symmetric
  logical, optional,               intent(in)    :: static_memory !< If present and true, this
  !! domain type is set up for static memory and error checking of
  !! various input values is performed against those in the input file.
  integer, optional, dimension(2), intent(in)    :: user_def_io_layout !< optional user-defined io_layout
  integer, optional,               intent(in)    :: NIHALO  !< Default halo sizes, required with static memory.
  integer, optional,               intent(in)    :: NJHALO  !< Default halo sizes, required with static memory.
  integer, optional,               intent(in)    :: NIGLOBAL!< Total domain sizes, required with static memory.
  integer, optional,               intent(in)    :: NJGLOBAL!< Total domain sizes, required with static memory.
  integer, optional,               intent(in)    :: NIPROC  !< Processor counts, required with static memory.
  integer, optional,               intent(in)    :: NJPROC  !< Processor counts, required with static memory.
  integer, optional,               intent(in)    :: X_FLAGS !< Parameters to use for x axis
  integer, optional,               intent(in)    :: Y_FLAGS !< Parameters to use for y axis
  integer, dimension(2), optional, intent(inout) :: min_halo!< If present, this sets the minimum halo size for this domain in the i- and j-
                                                            !! directions, and returns the actual halo size used.
  ! local
  integer, dimension(4) :: global_indices
  integer, dimension(2) :: layout = (/ 0, 0 /)
  integer, dimension(2) :: io_layout = (/ 0, 0 /)
  integer :: nihalo_dflt, njhalo_dflt
  integer :: pe, proc_used
  logical :: reentrant_x, reentrant_y, tripolar_N, is_static
  logical  :: mask_table_exists
  character(len=128) :: mask_table, inputdir
  character(len=64)  :: dom_name, inc_nm
  character(len=200) :: mesg
  integer :: xsiz, ysiz, nip_parsed, njp_parsed
  integer :: isc,iec,jsc,jec ! The bounding indices of the computational domain.
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=40) :: nihalo_nm, njhalo_nm, layout_nm, io_layout_nm, masktable_nm
  character(len=40) :: niproc_nm, njproc_nm
  integer :: xhalo_d2,yhalo_d2

  if (.not.associated(MOM_dom%mpp_domain)) then
    allocate(MOM_dom%mpp_domain)
  endif
  ! define the halo sizes
  nihalo_dflt = 4 ; njhalo_dflt = 4
  if (present(NIHALO)) nihalo_dflt = NIHALO
  if (present(NJHALO)) njhalo_dflt = NJHALO
  ! define the min_halo values
  if (present(min_halo)) then
    MOM_dom%nihalo = max(MOM_dom%nihalo, min_halo(1))
    min_halo(1) = MOM_dom%nihalo
    MOM_dom%njhalo = max(MOM_dom%njhalo, min_halo(2))
    min_halo(2) = MOM_dom%njhalo
  endif
  ! determine whether the domain is symmetric
  MOM_dom%symmetric = .true.
  if (present(symmetric)) MOM_dom%symmetric = symmetric
  ! define the global indices
  if (present(NIGLOBAL)) MOM_dom%niglobal = NIGLOBAL
  if (present(NJGLOBAL)) MOM_dom%njglobal = NJGLOBAL
  global_indices(1) = 1 ; global_indices(2) = MOM_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = MOM_dom%njglobal
  ! define the layout
  if (present(NIPROC)) layout(1) = NIPROC
  if (present(NJPROC)) layout(2) = NJPROC
  proc_used = mpp_npes()
  if ( layout(1)==0 .and. layout(2)==0 ) call mpp_define_layout(global_indices, proc_used, layout)
  if ( layout(1)/=0 .and. layout(2)==0 ) layout(2) = proc_used/layout(1)
  if ( layout(1)==0 .and. layout(2)/=0 ) layout(1) = proc_used/layout(2)
  MOM_dom%layout = layout
   ! check for a mask table file
  mask_table = "mask_table.nc"
  mask_table_exists = file_exists(mask_table)
  MOM_dom%X_FLAGS = 0
  MOM_dom%Y_FLAGS = 0
  if (present(X_FLAGS)) MOM_dom%X_FLAGS = X_FLAGS
  if (present(Y_FLAGS)) MOM_dom%Y_FLAGS = Y_FLAGS
  dom_name = "MOM"
  call mpp_domains_init()
  ! define the domain
  if (mask_table_exists) then
    call mpp_define_domains(global_indices, layout, MOM_dom%mpp_domain, &
                xflags=MOM_dom%X_FLAGS, yflags=MOM_dom%Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name, &
                maskmap=MOM_dom%maskmap )
  else
    call mpp_define_domains(global_indices, layout, MOM_dom%mpp_domain, &
                xflags=MOM_dom%X_FLAGS, yflags=MOM_dom%Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name)
  endif
  ! define the io layout
  if (present(user_def_io_layout)) io_layout = user_def_io_layout
  call mpp_define_io_domain(MOM_dom%mpp_domain, io_layout)
  MOM_dom%io_layout = io_layout

end subroutine domains_init

end module supergrid_setup
