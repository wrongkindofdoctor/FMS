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
!>@description This program tests the functionality of fms2-io to read a data array
!! into a supergrid domain.
program test_io_supergrid

use mpp_domains_mod, only : mpp_domains_init, mpp_define_domains, mpp_define_io_domain, mpp_define_layout
use mpp_domains_mod, only : mpp_domains_exit, mpp_get_domain_extents, mpp_deallocate_domain, domain2D
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain,mpp_get_global_domain
use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE, NORTH, EAST, CORNER, CENTER
use mpp_domains_mod, only : mpp_get_io_domain
use mpp_mod, only : mpp_init, mpp_exit, mpp_npes
use fms_io_mod, only : fms_io_init, fms_io_exit
use netcdf_io_mod, only : netcdf_io_init
use supergrid_setup, only : MOM_domain_type, hor_index_type, domains_init, hor_index_init, create_supergrid_domain
use supergrid_read_write, only : test_symmetric_read_write, test_symmetric_read_write_OG
implicit none

type(MOM_domain_type) :: G !< structure with the information for the base domain
type(MOM_domain_type) :: SGdom !< supergrid domain
type(hor_index_type) :: HI !< structure with the horizontal index information
character(len=10) :: nc_file_type
integer :: ncchksz = 1024*64
integer :: header_buffer_val = 16384

nc_file_type = "64bit"

call mpp_init()
call domains_init(G, symmetric=.false., user_def_io_layout=(/1,1/), &
                  NIHALO=2, NJHALO=2, NIGLOBAL=4, NJGLOBAL=5, NIPROC=4, NJPROC=1, &
                  X_FLAGS=CYCLIC_GLOBAL_DOMAIN, Y_FLAGS=FOLD_NORTH_EDGE)
call hor_index_init(G, HI)
call create_supergrid_domain(G, SGdom)
call fms_io_init()
call netcdf_io_init(ncchksz, header_buffer_val, nc_file_type)
call test_symmetric_read_write_OG(HI, G, SGdom)
call fms_io_exit()
call mpp_domains_exit()
call mpp_exit()

end program test_io_supergrid
