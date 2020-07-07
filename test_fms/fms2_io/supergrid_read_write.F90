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
!>@description This module contains the read and write routines required for test_io_supergrid
module supergrid_read_write
use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64, error_unit, output_unit
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain,mpp_get_global_domain
use mpp_domains_mod, only : NORTH, EAST, CORNER, CENTER, domain2D, mpp_get_io_domain
use mpp_domains_mod, only : mpp_update_domains, mpp_get_compute_domains
use mpp_mod, only : mpp_chksum, mpp_pe, mpp_error, mpp_root_pe, mpp_get_current_pelist
use mpp_mod, only : mpp_npes, mpp_sync_self, mpp_send, mpp_recv
use mpp_parameter_mod, only : FATAL, WARNING, NOTE, EVENT_SEND
use mpp_parameter_mod, only : To_East => WUPDATE, To_West => EUPDATE, Omit_Corners => EDGEUPDATE
use mpp_parameter_mod, only : To_North => SUPDATE, To_South => NUPDATE
use fms2_io_mod, only : file_exists, read_data, write_data, check_if_open, get_dimension_size, get_global_io_domain_indices
use fms2_io_mod, only : open_file, close_file, register_axis, register_field, FmsNetcdfDomainFile_t
use fms_io_mod, only : old_read_data => read_data
use netcdf_io_mod, only : netcdf_io_init
use supergrid_setup, only : MOM_domain_type, hor_index_type

implicit none

public :: test_symmetric_read_write
contains
!> Create an array, write it to a netcdf file, read it back in on the supergrid domain using fms_io
!! and fms2_io read_data interfaces, and check that the post-processed array from the new IO
!! matches the post-processed array from the old IO
subroutine test_symmetric_read_write(HI, G, SGdom)
  type(hor_index_type), intent(in) :: HI !< structure with the horizontal index information
  type(MOM_domain_type), intent(in) :: G ! regular domain
  type(MOM_domain_type), intent(inout) :: SGdom ! supergrid domain
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netdf file object
  real, dimension(2*HI%isd-3:2*HI%ied+1,2*HI%jsd-3:2*HI%jed+1) :: tmpC, tmpB
  real, dimension(HI%isd:HI%isd,HI%jsd:HI%jed) :: geoLonT, geoLonT_new
  real, allocatable :: tmpA(:,:)
  integer :: ni, nj,i2, j2, i, j, isg, ieg, jsg,jeg, isc, iec, jsc, jec, is, ie, js, je
  integer :: c(2), e(2), isd, ied, jsd, jed, xd_size, yd_size
  real, allocatable :: nxp(:), nyp(:)
  character(len=10) :: nc_mode
  logical :: file_open_success
  integer(kind=int64) :: old_chksum, new_chksum
  type(domain2d), pointer :: io_domain

  !Get the sizes of the global I/O domain and allocate buffer.
  io_domain => mpp_get_io_domain(SGdom%mpp_domain)
  call mpp_get_global_domain(io_domain, xsize=ni, ysize=nj, position=CORNER)
  if (.not.(allocated(tmpA))) allocate(tmpA(ni, nj))
  tmpA = 999.0
  ! write the data
  if (.not.(file_exists("some_data.nc"))) call write_2darray_to_ncfile(tmpA,SGdom)
  !call mpp_get_data_domain(SGdom%mpp_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, position=CORNER)
  !call mpp_get_compute_domain(SGdom%mpp_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=CORNER)
  !write(*,'(A,I2,A,I2,A,I2,A, I2)'), "isg, ieg, jsg, jeg",isg," ", ieg, " ", jsg, " ", jeg
  !write(*,'(A,I2,A,I2,A,I2,A, I2)'), "isc, iec, jsc, jec",isc," ", iec, " ", jsc, " ", jec
  !ie = iec - isg + 1
  !je = jec - jsg + 1
  !is = isc - isg + 1
  !js = jsc - jsg + 1
  !out_chksum = mpp_chksum(tmpA, pelist=(/mpp_pe()/))
  ! allocate an array that is the same size as tmpC, but not indexed with domain scheme
  tmpB(:,:) = 999.0; tmpC(:,:)=999.0
  ! read in the data using fms_io read_data
  write(output_unit, '(A,I2,A)'), "PE", mpp_pe()," is reading in data using fms_io"
  call old_read_data("./some_data.nc", 'x', tmpB, SGdom%mpp_domain, position=CORNER)
  ! Read in the array using fms2_io read_data
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, "some_data.nc", "read", SGdom%mpp_domain, is_restart=.false.)
  !call get_dimension_size(fileobj, 'nxp', dimsize_x)
  !call get_dimension_size(fileobj, 'nyp', dimsize_y)
  call register_axis(fileobj, 'nxp', 'x', domain_position=EAST)
  call register_axis(fileobj, 'nyp', 'y', domain_position=NORTH)
  write(output_unit, '(A,I2,A)'), "PE", mpp_pe()," is reading in data using fms2_io"
  call read_data(fileobj,'x', tmpC)
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! compare checksums
  old_chksum = mpp_chksum(tmpB, pelist=(/mpp_pe()/))
  new_chksum = mpp_chksum(tmpC, pelist=(/mpp_pe()/))
  if (new_chksum .ne. old_chksum) call mpp_error(WARNING, &
    "test_symmetric_read_write: tmpB and tmpC checksums do not match.")
  ! post process arrays and transfer every other value geoLonT and geoLonT_new
  geoLonT(:,:) = 0.0; geoLonT_new(:,:) = 0.0
  call pass_var_2d(tmpB, SGdom, position=CORNER)
  call extrapolate_metric(tmpB, 2*(HI%jsc-HI%jsd)+2, missing=999.0)
  call pass_var_2d(tmpC, SGdom, position=CORNER)
  call extrapolate_metric(tmpC, 2*(HI%jsc-HI%jsd)+2, missing=999.0)

  do j=HI%jsd,HI%jed ; do i=HI%isd,HI%ied ; i2 = 2*i ; j2 = 2*j
    geoLonT(i,j) = tmpB(i2-1,j2-1)
    geoLonT_new(i,j) = tmpC(i2-1,j2-1)
    write(output_unit,"(A,F7.5,A,F7.5)"), "geolont, geolont_new ",geoLonT(i,j)," ",geoLonT_new(i,j)
  enddo ; enddo
  ! compare checksums
  old_chksum = mpp_chksum(geoLonT, pelist=(/mpp_pe()/))
  new_chksum = mpp_chksum(geoLonT_new, pelist=(/mpp_pe()/))
  if (new_chksum .ne. old_chksum) call mpp_error(WARNING, &
    "test_symmetric_read_write: geoLonT checksums do not match.")

    !write(output_unit, '(A)'), "geoLonT
  
    !do i=1,size(geoLonT,1)
    !  write(output_unit,'(I2)') (geoLonT(i,j), j=1,size(geoLonT,2))
    !enddo
  !endif
  if (allocated(tmpA)) deallocate(tmpA)
  if (allocated(nxp)) deallocate(nxp)
  if (allocated(nyp)) deallocate(nyp)
  if (associated(io_domain)) nullify(io_domain)
end subroutine test_symmetric_read_write

subroutine write_2darray_to_ncfile(array, MOM_domain)
  real, dimension(:,:), intent(inout) :: array ! 2d array of values to write to the file
  type(MOM_domain_type), intent(in) :: MOM_domain ! MOM domain structure
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netdf file object
  character(len=10) :: nc_mode, nc_file_type
  logical :: file_open_success
  integer :: i, k, ni, nj 
  integer, allocatable, dimension(:) :: nxp, nyp
  ! generate the data
  call random_number(array)
  ! create the axes
  allocate(nxp(ni))
  allocate(nyp(nj))
  k=0.0
  do i=1,ni
    k=k+1.0
    nxp(i)=k
 !   write(output_unit,'(A,I2,A,F5.2)'), "index is ",i," nxp is ", nxp(i)
  enddo
  k=0.0
  do i=1,nj
    k=k+2.0
    nyp(i)=k
  !  write(output_unit,'(A,I2,A,F5.2)'), "index is ",i," nyp is ", nyp(i)
  enddo
  ! write the array to a file
  nc_file_type = "64bit"
  nc_mode = "write"
  if (file_exists("some_data.nc")) nc_mode = "overwrite"
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, "./some_data.nc", trim(nc_mode), MOM_domain%mpp_domain, is_restart=.false.)
  call register_axis(fileobj, 'nxp', 'x', domain_position=EAST)
  call register_axis(fileobj, 'nyp', 'y', domain_position=NORTH)
  call register_field(fileobj, 'x', "double", dimensions=(/'nxp', 'nyp'/))
  call write_data(fileobj, 'x', array)
  if (check_if_open(fileobj)) call close_file(fileobj)

end subroutine write_2darray_to_ncfile

!> Extrapolates missing metric data into all the halo regions.
subroutine extrapolate_metric(var, jh, missing)
  real, dimension(:,:), intent(inout) :: var     !< The array in which to fill in halos
  integer,              intent(in)    :: jh      !< The size of the halos to be filled
  real,       optional, intent(in)    :: missing !< The missing data fill value, 0 by default.
  ! Local variables
  real :: badval
  integer :: i,j

  badval = 0.0 ; if (present(missing)) badval = missing

  ! Fill in southern halo by extrapolating from the computational domain
  do j=lbound(var,2)+jh,lbound(var,2),-1 ; do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i,j+1)-var(i,j+2)
  enddo ; enddo

  ! Fill in northern halo by extrapolating from the computational domain
  do j=ubound(var,2)-jh,ubound(var,2) ; do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i,j-1)-var(i,j-2)
  enddo ; enddo

  ! Fill in western halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2) ; do i=lbound(var,1)+jh,lbound(var,1),-1
    if (var(i,j)==badval) var(i,j) = 2.0*var(i+1,j)-var(i+2,j)
  enddo ; enddo

  ! Fill in eastern halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2) ; do i=ubound(var,1)-jh,ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i-1,j)-var(i-2,j)
  enddo ; enddo

end subroutine extrapolate_metric

!> pass_var_2d does a halo update for a two-dimensional array.
subroutine pass_var_2d(array, MOM_dom, sideflag, complete, position, halo, inner_halo, clock)
  real, dimension(:,:),  intent(inout) :: array    !< The array which is having its halos points
                                                   !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be sent.
  integer,     optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  logical,     optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                   !! halo updates should be completed before
                                                   !! progress resumes.  Omitting complete is the
                                                   !! same as setting complete to .true.
  integer,     optional, intent(in)    :: position !< An optional argument indicating the position.
                                                   !!  This is usally CORNER, but is CENTER
                                                   !! by default.
  integer,     optional, intent(in)    :: halo     !< The size of the halo to update - the full halo
                                                   !! by default.
  integer,     optional, intent(in)    :: inner_halo !< The size of an inner halo to avoid updating,
                                                   !! or 0 to avoid updating symmetric memory
                                                   !! computational domain points.  Setting this >=0
                                                   !! also enforces that complete=.true.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  ! Local variables
  real, allocatable, dimension(:,:) :: tmp
  integer :: pos, i_halo, j_halo
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  integer :: inner, i, j, isfw, iefw, isfe, iefe, jsfs, jefs, jsfn, jefn
  integer :: dirflag
  logical :: block_til_complete
  integer, parameter :: To_All = To_East + To_West + To_North + To_South !< A flag for passing in all directions


  !if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true. ; if (present(complete)) block_til_complete = complete
  pos = CENTER ; if (present(position)) pos = position

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    ! Store the original values.
    allocate(tmp(size(array,1), size(array,2)))
    tmp(:,:) = array(:,:)
    block_til_complete = .true.
  endif ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position, &
                        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position)
  endif

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    call mpp_get_compute_domain(MOM_dom%mpp_domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(MOM_dom%mpp_domain, isd, ied, jsd, jed)
    ! Convert to local indices for arrays starting at 1.
    isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
    jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
    i_halo = min(inner_halo, isc-1) ; j_halo = min(inner_halo, jsc-1)

    ! Figure out the array index extents of the eastern, western, northern and southern regions to copy.
    if (pos == CENTER) then
      if (size(array,1) == ied) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call mpp_error(FATAL, "pass_var_2d: wrong i-size for CENTER array.") ; endif
      if (size(array,2) == jed) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call mpp_error(FATAL, "pass_var_2d: wrong j-size for CENTER array.") ; endif
    elseif (pos == CORNER) then
      if (size(array,1) == ied) then
        isfw = max(isc - (i_halo+1), 1) ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      elseif (size(array,1) == ied+1) then
        isfw = isc - i_halo ; iefw = isc+1 ; isfe = iec+1 ; iefe = min(iec + 1 + i_halo, ied+1)
      else ; call mpp_error(FATAL, "pass_var_2d: wrong i-size for CORNER array.") ; endif
      if (size(array,2) == jed) then
        jsfs = max(jsc - (j_halo+1), 1) ; jefs = jsc ; jsfn = jec ; jefn = jec + j_halo
      elseif (size(array,2) == jed+1) then
        jsfs = jsc - j_halo ; jefs = jsc+1 ; jsfn = jec+1 ; jefn = min(jec + 1 + j_halo, jed+1)
      else ; call mpp_error(FATAL, "pass_var_2d: wrong j-size for CORNER array.") ; endif
    else
      call mpp_error(FATAL, "pass_var_2d: Unrecognized position")
    endif

    ! Copy back the stored inner halo points
    do j=jsfs,jefn ; do i=isfw,iefw ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefn ; do i=isfe,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefs ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfn,jefn ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo

    deallocate(tmp)
  endif ; endif

  !if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif
end subroutine pass_var_2d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_2d(buf, sizes)

  real(kind=real64), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2)))
end subroutine allocate_array_real64_2d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_2d(section, array, s, c)

  real(kind=real64), dimension(:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ) = section(:,:)
end subroutine put_array_section_real64_2d

!> @brief Get a section of larger array.
subroutine get_array_section_real64_2d(section, array, s, c)

  real(kind=real64), dimension(:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  section(:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 )
end subroutine get_array_section_real64_2d


end module supergrid_read_write
