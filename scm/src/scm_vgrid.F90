!> \file scm_vgrid.f90
!!  Contains the vertical grid setup routines.

module scm_vgrid

use scm_kinds, only: sp, dp, qp, kind_scm_dp, kind_scm_sp
use scm_physical_constants, only : con_cp, con_rocp, con_fvirt, con_g, con_rd

implicit none

private
public get_FV3_vgrid, calc_pres_exner_geopotential, calc_geopotential

logical :: verbose = .true.

contains

!> \ingroup SCM
!! @{
!! \defgroup vgrid scm_vgrid
!! @{
!! Contains the vertical grid setup routines.

!GJF: most of this was obtained from FV3/atmos_cubed_sphere/tools/fv_eta.F90 from FV3
!GJF: current as of March 2022
!GJF: assuming proprocessor variable USE_VAR_ETA is NOT set
subroutine get_FV3_vgrid(scm_input, scm_state)
  use scm_type_defs, only: scm_input_type, scm_state_type

      type(scm_input_type), intent(in) :: scm_input
      type(scm_state_type), intent(inout) :: scm_state

      integer :: km           ! vertical dimension
      integer :: ks           ! number of pure p layers
      real(kind=dp) :: ptop         ! model top (Pa)

! local
      real(kind=dp)               :: pres_sfc_inv, p_ref, ak_tmp, bk_tmp

      real:: p0=1000.E2
      real:: pc=200.E2

      real :: pint = 100.E2
      real :: stretch_fac = 1.03
      integer :: auto_routine = 0

      integer  k, last_index, mid_index, ierr, dummy, n_levels_file

      character(len=80)     :: line
      character(len=16)     :: file_format
      integer :: nx, ny
      real(kind_scm_dp), allocatable :: pres_l_row(:), pres_i(:,:)
      real(kind_scm_dp), parameter :: zero_dp = 0.0
      ! added for forcing initialized pressure to be single precision for
      ! single and double precision runs
      real(kind_scm_sp), parameter :: zero_sp = 0.0
      real(kind_scm_sp), allocatable :: pres_l_row_sp(:)

#include "fv_eta.h"

      km = scm_state%n_levels
      ptop = 1.
! Definition: press(i,j,k) = ak(k) + bk(k) * ps(i,j)

      if (trim(scm_state%npz_type) == 'superC' .or. trim(scm_state%npz_type) == 'superK') then
        auto_routine = 1
        select case (km)
          case (20)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (24)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (30)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (40)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (50)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (60)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (80)
             ptop = 56.e2
             pint = ptop
             stretch_fac = 1.03
          case (90)          ! super-duper cell
             ptop = 40.e2
             stretch_fac = 1.025
             auto_routine = 2
        end select
      else if (trim(scm_state%npz_type) == 'input') then

        !> - Open the appropriate file.
        open(unit=1, file=scm_state%vert_coord_file, status='old', action='read', iostat=ierr)
        if(ierr /= 0) then
          write(*,*) 'There was an error opening the file ', scm_state%vert_coord_file, ' in the run directory. &
            Error code = ',ierr
          error stop
        endif

        !> The file being read in must have the following format:
        !!       include a single line description: number of coefficients, number of layers
        !!       ak/bk pairs, with each pair occupying a single line
        !!       the pairs must be ordered from surface to TOA
        !!       the pairs define the interfaces of the grid to create levels-1 layer

        !> - The first line contains the number of coefficients and number of levels
        read(1,*) dummy, n_levels_file
        if (n_levels_file /= scm_state%n_levels) then
          write(*,*) 'There is a mismatch in the number of levels expected and the number of coefficients supplied in the file ',scm_state%vert_coord_file
          error stop
        end if
        !> - Read in the coefficient data.
        do k=1, km+1
          read(1,*)  scm_state%a_k(k), scm_state%b_k(k)
        end do
        close(1)

        ! flip scm_state%a_k, scm_state%b_k in vertical (a_k and b_k are expected to be TOA-to-surface at this point)
        mid_index = (km+1)/2
        last_index = km+1
        do k = 1, mid_index
            ak_tmp = scm_state%a_k(k)
            bk_tmp = scm_state%b_k(k)
            scm_state%a_k(k) = scm_state%a_k(last_index)
            scm_state%b_k(k) = scm_state%b_k(last_index)
            scm_state%a_k(last_index) = ak_tmp
            scm_state%b_k(last_index) = bk_tmp
            last_index = last_index - 1
        end do

      else
        select case (km)
          case (5,10) ! does this work????
            ! Equivalent Shallow Water: for modon test
            ptop = 500.e2
            ks = 0
            do k=1,km+1
               scm_state%b_k(k) = real(k-1) / real (km)
               scm_state%a_k(k) = ptop*(1.-scm_state%b_k(k))
            enddo
          case (24)
            ks = 5
            do k=1,km+1
              scm_state%a_k(k) = a24(k)
              scm_state%b_k(k) = b24(k)
            enddo
          case (26)
            ks = 7
            do k=1,km+1
              scm_state%a_k(k) = a26(k)
              scm_state%b_k(k) = b26(k)
            enddo
          case (30)          ! For Baroclinic Instability Test
            ptop = 2.26e2
            pint = 250.E2
            stretch_fac = 1.03
            auto_routine = 1
          case (31)               ! N = 4, M=2
            if (trim(scm_state%npz_type) == 'lowtop') then
              ptop = 300.
              pint = 100.E2
              stretch_fac = 1.035
              auto_routine = 5
            else
              ptop = 100.
              stretch_fac = 1.035
              auto_routine = 1
            endif
          case (32)
            if (trim(scm_state%npz_type) == 'old32') then
              ks = 13              ! high-res trop_32 setup
              do k=1,km+1
                scm_state%a_k(k) = a32old(k)
                scm_state%b_k(k) = b32old(k)
              enddo
            elseif (trim(scm_state%npz_type) == 'lowtop') then
              ptop = 100.
              stretch_fac = 1.035
              auto_routine = 1
            else
              ks = 7
              do k=1,km+1
                scm_state%a_k(k) = a32(k)
                scm_state%b_k(k) = b32(k)
              enddo
            endif
        !miz
          case (33)
            ks = 7
            do k=1,km+1
              scm_state%a_k(k) = a33(k)
              scm_state%b_k(k) = b33(k)
            enddo
            !miz
          case (39)               ! N = 5
            ptop = 100.
            stretch_fac = 1.035
            auto_routine = 1
          case (40)
            ptop = 50.e2   ! For super cell test
            pint = 300.E2
            stretch_fac = 1.03
            auto_routine = 1
          case (41)
            ptop = 100.
            pint = 100.E2
            stretch_fac = 1.035
            auto_routine = 1
          case (47)
            if (trim(scm_state%npz_type) == 'lowtop') then
              ptop = 100.
              stretch_fac = 1.035
              auto_routine = 1
            else
              !         ks = 27       ! high-res trop-strat
              ks = 20       ! Oct 23, 2012
              do k=1,km+1
                scm_state%a_k(k) = a47(k)
                scm_state%b_k(k) = b47(k)
              enddo
            endif
          case (48)
            ks = 28
            do k=1,km+1
              scm_state%a_k(k) = a48(k)
              scm_state%b_k(k) = b48(k)
            enddo
          case (50)
            ! ! *Very-low top: for idealized super-cell simulation:
            ! ptop = 50.e2
            ! pint = 250.E2
            ! stretch_fac = 1.03
            ! auto_routine = 1
            ks = 19
            do k=1,km+1
              scm_state%a_k(k) = a50(k)
              scm_state%b_k(k) = b50(k)
            enddo
          case (51)
            if (trim(scm_state%npz_type) == 'lowtop') then
              ptop = 100.
              stretch_fac = 1.03
              auto_routine = 1
            elseif (trim(scm_state%npz_type) == 'meso') then
              ptop = 20.E2
              pint = 100.E2
              stretch_fac = 1.05
              auto_routine = 1
            elseif (trim(scm_state%npz_type) == 'meso2') then
              ptop = 1.E2
              pint = 100.E2
              stretch_fac = 1.05
              auto_routine = 6
            else
              ptop = 100.
              pint = 100.E2
              stretch_fac = 1.035
              auto_routine = 1
            endif
          case (52)
            if (trim(scm_state%npz_type) == 'rce') then
              ptop = 30.e2    ! for special DPM RCE experiments
              stretch_fac = 1.03
              auto_routine = 1
            else
              ks = 35         ! pint = 223
              do k=1,km+1
                scm_state%a_k(k) = a52(k)
                scm_state%b_k(k) = b52(k)
              enddo
            endif
          case (54)
            ks = 11         ! pint =  109.4
            do k=1,km+1
              scm_state%a_k(k) = a54(k)
              scm_state%b_k(k) = b54(k)
            enddo
          case (55)  ! Mid-top:             ! N = 7
            ptop = 10.
            pint = 100.E2
            stretch_fac = 1.035
            auto_routine = 1
          case (56)
            ks = 26
            do k=1,km+1
              scm_state%a_k(k) = a56(k)
              scm_state%b_k(k) = b56(k)
            enddo
          case (60)
            if (trim(scm_state%npz_type) == 'gfs') then
              ks = 20
              do k=1,km+1
                scm_state%a_k(k) = a60gfs(k)
                scm_state%b_k(k) = b60gfs(k)
              enddo
            else if (trim(scm_state%npz_type) == 'BCwave') then
              ptop = 3.e2
              !            pint = 250.E2
              pint = 300.E2    ! revised for Moist test
              stretch_fac = 1.03
              auto_routine = 1
            else if (trim(scm_state%npz_type) == 'meso') then
              ptop = 40.e2
              pint = 250.E2
              stretch_fac = 1.03
              auto_routine = 1
            else
              ks = 19
              do k=1,km+1
                scm_state%a_k(k) = a60(k)
                scm_state%b_k(k) = b60(k)
              enddo
            endif
          case (63)
            if (trim(scm_state%npz_type) == 'meso') then
              ks = 11
              do k=1,km+1
                scm_state%a_k(k) = a63meso(k)
                scm_state%b_k(k) = b63meso(k)
              enddo
            elseif (trim(scm_state%npz_type) == 'hitop') then
              ptop = 1.   ! high top
              pint = 100.E2
              stretch_fac = 1.035
              auto_routine = 1
            else!if (trim(scm_state%npz_type) == 'gfs') then
              !Used for SHiELD
              ! GFS L64 equivalent setting
              ks = 23
              do k=1,km+1
                scm_state%a_k(k) = a63(k)
                scm_state%b_k(k) = b63(k)
              enddo
            endif
          case (64)
            if (trim(scm_state%npz_type) == 'gfs') then
              ks = 23
              do k=1,km+1
                scm_state%a_k(k) = a64gfs(k)
                scm_state%b_k(k) = b64gfs(k)
              enddo
            else
              ks = 46
              do k=1,km+1
                scm_state%a_k(k) = a64(k)
                scm_state%b_k(k) = b64(k)
              enddo
            endif
         ! xi chen's l65
         case (65)
            ks = 29
            do k=1,km+1
               scm_state%a_k(k) = a65(k)
               scm_state%b_k(k) = b65(k)
            enddo

            !-->cjg
          case (68)
            ks = 27
            do k=1,km+1
              scm_state%a_k(k) = a68(k)
              scm_state%b_k(k) = b68(k)
            enddo
          case (71)               ! N = 9
            ptop = 1.
            stretch_fac = 1.03
            auto_routine = 1
         ! kgao: introduce EMC's L75 config
         case (75)
            if (trim(scm_state%npz_type) == 'emc') then
               ! EMC's L75 config
               ks = 12
               do k=1,km+1
                  scm_state%a_k(k) = a75(k)
                  scm_state%b_k(k) = b75(k)
               enddo
            else
               ! HS-SGO test configuration
               pint = 100.E2
               ptop = 10.E2
               stretch_fac = 1.035
               auto_routine = 6
            endif

          case (79)               ! N = 10, M=5
            if (trim(scm_state%npz_type) == 'gcrm') then
              pint = 100.E2
              ptop = 3.E2
              stretch_fac = 1.035
              auto_routine = 6
            else
              ptop = 1.
              stretch_fac = 1.03
              auto_routine = 1
            endif
         ! kgao L88
         case (88)
            ks = 20 !19 bug fix
            do k=1,km+1
               scm_state%a_k(k) = a88(k)
               scm_state%b_k(k) = b88(k)
            enddo
         case (90)          ! super-duper cell
            ptop = 40.e2
            stretch_fac = 1.025
            auto_routine = 2
          case (91) ! NGGPS_GFS
            pint = 100.E2
            ptop = 40.
            stretch_fac = 1.029
            auto_routine = 6
          case (95)
            ! Mid-top settings:
            pint = 100.E2
            ptop = 20.
            stretch_fac = 1.029
            auto_routine = 6
          case (96)
            ks = 27
            do k=1,km+1
              scm_state%a_k(k) = a96(k)
              scm_state%b_k(k) = b96(k)
            enddo
            !<--cjg
          case (100)
            ks = 38
            do k=1,km+1
              scm_state%a_k(k) = a100(k)
              scm_state%b_k(k) = b100(k)
            enddo
          case (104)
            ks = 73
            do k=1,km+1
              scm_state%a_k(k) = a104(k)
              scm_state%b_k(k) = b104(k)
            enddo
      ! IFS-like L125
          case (125)
            ks = 33
            ptop = a125(1)
            pint = a125(ks+1)
            do k=1,km+1
              scm_state%a_k(k) = a125(k)
              scm_state%b_k(k) = b125(k)
            enddo
          case (127)               ! N = 10, M=5
             if (trim(scm_state%npz_type) == 'hitop') then
                ptop = 1.
                stretch_fac = 1.03
                auto_routine = 2
             elseif (trim(scm_state%npz_type) == 'gfs') then
                ks = 39
                ptop = a127(1)
                pint = a127(ks+1)
                do k=1,km+1
                   scm_state%a_k(k) = a127(k)
                   scm_state%b_k(k) = b127(k)
                enddo
             else
                ptop = 1.
                pint = 75.E2
                stretch_fac = 1.028
                auto_routine = 6
             endif
          case (151)
            !LES applications
            ptop = 75.e2
            pint = 500.E2
            stretch_fac = 1.01
            auto_routine = 3
          case default
            if(trim(scm_state%npz_type) == 'hitop') then
              ptop = 1.
              pint = 100.E2
            elseif(trim(scm_state%npz_type) == 'midtop') then
              ptop = 10.
              pint = 100.E2
            elseif(trim(scm_state%npz_type) == 'lowtop') then
              ptop = 1.E2
              pint = 100.E2
            endif

            if (trim(scm_state%npz_type) == 'gfs') then
              auto_routine = 6
            elseif(trim(scm_state%npz_type) == 'les') then
              auto_routine = 3
            elseif(trim(scm_state%npz_type) == 'mountain_wave') then
              auto_routine = 4
            elseif (km > 79) then
              auto_routine = 2
            else
              auto_routine = 1
            endif
        end select
      endif ! superC/superK

      select case (auto_routine)
        case (1)
          call var_hi(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint, stretch_fac)
        case (2)
          call var_hi2(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint, stretch_fac)
        case (3)
          call var_les(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint, stretch_fac)
        case (4)
          call mount_waves(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint)
        case (5)
          call var_dz(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint, stretch_fac)
        case (6)
          call var_gfs(km, scm_state%a_k, scm_state%b_k, ptop, ks, pint, stretch_fac)
      end select

      call check_eta_levels (scm_state%a_k, scm_state%b_k)

      if (verbose) then
         write(*, '(A4, A13, A13, A11)') 'klev', 'ak', 'bk', 'p_ref'
         do k=1,km+1
            write(*,'(I4, F13.5, F13.5, F11.2)') k, scm_state%a_k(k), scm_state%b_k(k), 1000.E2*scm_state%b_k(k) + scm_state%a_k(k)
         enddo
      endif

      !> - Calculate interface pressures, sigma, and exner function.

      ! flip scm_state%a_k, scm_state%b_k in vertical
      mid_index = (km+1)/2
      last_index = km+1
      do k = 1, mid_index
          ak_tmp = scm_state%a_k(k)
          bk_tmp = scm_state%b_k(k)
          scm_state%a_k(k) = scm_state%a_k(last_index)
          scm_state%b_k(k) = scm_state%b_k(last_index)
          scm_state%a_k(last_index) = ak_tmp
          scm_state%b_k(last_index) = bk_tmp
          last_index = last_index - 1
      end do

      p_ref = scm_input%input_pres_surf(1)
      pres_sfc_inv = 1.0/p_ref

      nx = size(scm_state%pres_i, 1)
      ny = size(scm_state%pres_i, 2)
      allocate(pres_i(nx,ny), source=zero_dp)
      do k=1, km+1
        pres_i(:,k) = scm_state%a_k(k) + scm_state%b_k(k)*p_ref
        scm_state%si(:,k) = scm_state%a_k(k)*pres_sfc_inv + scm_state%b_k(k)
        scm_state%exner_i(:,k) = (scm_state%pres_i(:,k)/1.0E5)**con_rocp
      end do
      scm_state%pres_i = pres_i

      !> - Calculate layer center pressures, sigma, and exner function.
      allocate(pres_l_row(nx), source=zero_dp)
      allocate(pres_l_row_sp(nx), source=zero_sp)
      do k=1, km
        pres_l_row_sp = ((1.0/(con_rocp+1.0))*&
          (pres_i(:,k)**(con_rocp+1.0) - pres_i(:,k+1)**(con_rocp+1.0))/ &
          (pres_i(:,k) - pres_i(:,k+1)))**(1.0/con_rocp)
        scm_state%pres_l(:,k) = pres_l_row_sp
        scm_state%sl(:,k) = 0.5*(scm_state%si(:,k) + scm_state%si(:,k+1))

        scm_state%exner_l(:,k) = (scm_state%pres_l(:,k)/1.0E5)**con_rocp

      end do

end subroutine get_FV3_vgrid

subroutine var_hi(km, ak, bk, ptop, ks, pint, s_rate)
 integer, intent(in):: km
 real,    intent(in):: ptop
 real,    intent(in):: s_rate        !< between [1. 1.1]
 real,    intent(out):: ak(km+1), bk(km+1)
 real,    intent(inout):: pint
 integer, intent(out):: ks
! Local
 real, parameter:: p00 = 1.E5
 real, dimension(km+1):: ze, pe1, peln, eta
 real, dimension(km):: dz, s_fac, dlnp
 real ztop, t0, dz0, sum1, tmp1
 real ep, es, alpha, beta, gama
!---- Tunable parameters:
 integer:: k_inc = 15   !<number of layers from bottom up to near const dz region
 real:: s0 = 0.10 !< lowest layer stretch factor
!-----------------------
 real:: s_inc
 integer  k

    pe1(1) = ptop
    peln(1) = log(pe1(1))
    pe1(km+1) = p00
    peln(km+1) = log(pe1(km+1))

    t0 = 270.
    ztop = con_rd/con_g*t0*(peln(km+1) - peln(1))

     s_inc = (1.-s0) / real(k_inc)
     s_fac(km)  = s0

     do k=km-1, km-k_inc, -1
        s_fac(k)  = s_fac(k+1) + s_inc
     enddo

     s_fac(km-k_inc-1) = 0.5*(s_fac(km-k_inc) + s_rate)

#ifdef HIWPP
     do k=km-k_inc-2, 4, -1
        s_fac(k) = s_rate * s_fac(k+1)
     enddo
     s_fac(3) = 0.5*(1.15+s_rate)*s_fac(4)
     s_fac(2) = 1.15 *s_fac(3)
     s_fac(1) = 1.3 *s_fac(2)
#else
     do k=km-k_inc-2, 9, -1
        s_fac(k) = s_rate * s_fac(k+1)
     enddo

     s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
     s_fac(7) = 1.1 *s_fac(8)
     s_fac(6) = 1.15*s_fac(7)
     s_fac(5) = 1.2 *s_fac(6)
     s_fac(4) = 1.3 *s_fac(5)
     s_fac(3) = 1.4 *s_fac(4)
     s_fac(2) = 1.45 *s_fac(3)
     s_fac(1) = 1.5 *s_fac(2)
#endif

     sum1 = 0.
     do k=1,km
        sum1 = sum1 + s_fac(k)
     enddo

     dz0 = ztop / sum1

     do k=1,km
        dz(k) = s_fac(k) * dz0
     enddo

     ze(km+1) = 0.
     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo

! Re-scale dz with the stretched ztop
     do k=1,km
        dz(k) = dz(k) * (ztop/ze(1))
     enddo

     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo
!     ze(1) = ztop

     if ( verbose ) then
          write(*,*) 'var_hi: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
!          do k=1,km
!             write(*,*) k, s_fac(k)
!          enddo
     endif

     call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
     do k=1,km
         dz(k) = ze(k) - ze(k+1)
       dlnp(k) = con_g*dz(k) / (con_rd*t0)
     enddo
     do k=2,km
        peln(k) = peln(k-1) + dlnp(k-1)
         pe1(k) = exp(peln(k))
     enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
     ks = 0
     do k=2,km
        if ( pint < pe1(k)) then
             ks = k-1
             exit
        endif
     enddo
     if ( verbose ) then
        write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
        write(*,*) 'ptop =', ptop
     endif
     pint = pe1(ks+1)

#ifdef NO_UKMO_HB
     do k=1,ks+1
        ak(k) = pe1(k)
        bk(k) = 0.
     enddo

     do k=ks+2,km+1
        bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
        ak(k) =  pe1(k) - bk(k) * pe1(km+1)
     enddo
     bk(km+1) = 1.
     ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
     do k=1,km+1
        eta(k) = pe1(k) / pe1(km+1)
     enddo

     ep =  eta(ks+1)
     es =  eta(km)
!     es =  1.
     alpha = (ep**2-2.*ep*es) / (es-ep)**2
     beta  = 2.*ep*es**2 / (es-ep)**2
     gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
     do k=1,ks+1
        ak(k) = eta(k)*1.e5
        bk(k) = 0.
     enddo

     do k=ks+2, km
        ak(k) = alpha*eta(k) + beta + gama/eta(k)
        ak(k) = ak(k)*1.e5
     enddo
        ak(km+1) = 0.

     do k=ks+2, km
        bk(k) = (pe1(k) - ak(k))/pe1(km+1)
     enddo
        bk(km+1) = 1.
#endif

     if ( verbose ) then
         write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
         do k=1,km
            write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
         enddo
         tmp1 = ak(ks+1)
         do k=ks+1,km
            tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
         enddo
         write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
     endif
end subroutine var_hi

subroutine var_hi2(km, ak, bk, ptop, ks, pint, s_rate)
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
 ! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 270.
     ztop = con_rd/con_g*t0*(peln(km+1) - peln(1))

      s_fac(km  ) = 0.15
      s_fac(km-1) = 0.20
      s_fac(km-2) = 0.30
      s_fac(km-3) = 0.40
      s_fac(km-4) = 0.50
      s_fac(km-5) = 0.60
      s_fac(km-6) = 0.70
      s_fac(km-7) = 0.80
      s_fac(km-8) = 0.90
      s_fac(km-9) = 0.95
      s_fac(km-10) = 0.5*(s_fac(km-9) + s_rate)

      do k=km-11, 8, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(7) = 0.5*(1.1+s_rate)*s_fac(9)
      s_fac(6) = 1.05*s_fac(7)
      s_fac(5) = 1.1*s_fac(6)
      s_fac(4) = 1.15*s_fac(5)
      s_fac(3) = 1.2*s_fac(4)
      s_fac(2) = 1.3*s_fac(3)
      s_fac(1) = 1.4*s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

 ! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
 !     ze(1) = ztop

      if ( verbose ) write(*,*) 'var_hi2: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

 ! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = con_g*dz(k) / (con_rd*t0)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo

 ! Pe(k) = ak(k) + bk(k) * PS
 ! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( verbose ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
 ! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
 !     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

 ! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( verbose ) then
          write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
          do k=1,km
             write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
          enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif


end subroutine var_hi2

subroutine var_les(km, ak, bk, ptop, ks, pint, s_rate)
  implicit none
  integer, intent(in):: km
  real,    intent(in):: ptop
  real,    intent(in):: s_rate        !< between [1. 1.1]
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(inout):: pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, s_fac, dlnp, pm, dp, dk
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama
  real, parameter:: akap = 2./7.
!---- Tunable parameters:
  integer:: k_inc = 10   !< number of layers from bottom up to near const dz region
  real:: s0 = 0.8     !< lowest layer stretch factor
!-----------------------
  real:: s_inc
  integer  k

     pe1(1) = ptop
     peln(1) = log(pe1(1))
     pe1(km+1) = p00
     peln(km+1) = log(pe1(km+1))

     t0 = 273.
     ztop = con_rd/con_g*t0*(peln(km+1) - peln(1))

      s_inc = (1.-s0) / real(k_inc)
      s_fac(km)  = s0

      do k=km-1, km-k_inc, -1
         s_fac(k)  = s_fac(k+1) + s_inc
      enddo

      s_fac(km-k_inc-1) = 0.5*(s_fac(km-k_inc) + s_rate)

      do k=km-k_inc-2, 5, -1
         s_fac(k) = s_rate * s_fac(k+1)
      enddo

      s_fac(4) = 0.5*(1.1+s_rate)*s_fac(5)
      s_fac(3) = 1.1 *s_fac(4)
      s_fac(2) = 1.1 *s_fac(3)
      s_fac(1) = 1.1 *s_fac(2)

      sum1 = 0.
      do k=1,km
         sum1 = sum1 + s_fac(k)
      enddo

      dz0 = ztop / sum1

      do k=1,km
         dz(k) = s_fac(k) * dz0
      enddo

      ze(km+1) = 0.
      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo

! Re-scale dz with the stretched ztop
      do k=1,km
         dz(k) = dz(k) * (ztop/ze(1))
      enddo

      do k=km,1,-1
         ze(k) = ze(k+1) + dz(k)
      enddo
!     ze(1) = ztop

      if ( verbose ) then
           write(*,*) 'var_les: computed model top (m)=', ztop, ' bottom/top dz=', dz(km), dz(1)
!           do k=1,km
!              write(*,*) k, s_fac(k)
!           enddo
      endif

      call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 2)

! Given z --> p
      do k=1,km
          dz(k) = ze(k) - ze(k+1)
        dlnp(k) = con_g*dz(k) / (con_rd*t0)
        !write(*,*) k, dz(k)
      enddo
      do k=2,km
         peln(k) = peln(k-1) + dlnp(k-1)
          pe1(k) = exp(peln(k))
      enddo


! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo
      if ( verbose ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
      endif
      pint = pe1(ks+1)

      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo

      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.

      if ( verbose ) then
 !         write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
 !         do k=1,km
 !            pm(k) = 0.5*(pe1(k)+pe1(k+1))/100.
 !            write(*,*) k, pm(k), dz(k)
 !         enddo
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
          write(*,800) (pm(k), k=km,1,-1)
      endif

    do k=1,km
       dp(k) = (pe1(k+1) - pe1(k))/100.
       dk(k) =  pe1(k+1)**akap - pe1(k)**akap
    enddo

800 format(1x,5(1x,f9.4))

end subroutine var_les

subroutine mount_waves(km, ak, bk, ptop, ks, pint)
  integer, intent(in):: km
  real,    intent(out):: ak(km+1), bk(km+1)
  real,    intent(out):: ptop, pint
  integer, intent(out):: ks
! Local
  real, parameter:: p00 = 1.E5
  real, dimension(km+1):: ze, pe1, peln, eta
  real, dimension(km):: dz, dlnp
  real ztop, t0, dz0, sum1, tmp1
  real ep, es, alpha, beta, gama, s_fac
  integer  k, k500

  pint = 300.e2
!      s_fac = 1.05
!      dz0 = 500.
  if ( km <= 60 ) then
       s_fac = 1.0
       dz0 = 500.
  else
       s_fac = 1.
       dz0 = 250.
  endif

! Basic parameters for HIWPP mountain waves
   t0 = 300.
! ztop = 20.0e3; 500-m resolution in halft of the vertical domain
! ztop = real(km-1)*500.
!-----------------------
! Compute temp ptop based on isothermal atm
! ptop = p00*exp(-grav*ztop/(rdgas*t0))

! Lowest half has constant resolution
     ze(km+1) = 0.
     do k=km, km-19, -1
        ze(k) = ze(k+1) + dz0
     enddo

! Stretching from 10-km and up:
     do k=km-20, 3,  -1
        dz0 = s_fac * dz0
        ze(k) = ze(k+1) + dz0
     enddo
     ze(2) = ze(3) + sqrt(2.)*dz0
     ze(1) = ze(2) + 2.0*dz0

!    call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
     do k=1,km
         dz(k) = ze(k) - ze(k+1)
       dlnp(k) = con_g*dz(k) / (con_rd*t0)
     enddo

      pe1(km+1) = p00
     peln(km+1) = log(p00)
     do k=km,1,-1
        peln(k) = peln(k+1) - dlnp(k)
         pe1(k) = exp(peln(k))
     enddo

! Comnpute new ptop
     ptop = pe1(1)

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
      ks = 0
      do k=2,km
         if ( pint < pe1(k)) then
              ks = k-1
              exit
         endif
      enddo

      if ( verbose ) then
         write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
         write(*,*) 'Modified ptop =', ptop, ' ztop=', ze(1)/1000.
         do k=1,km
            write(*,*) k, 'ze =', ze(k)/1000.
         enddo
      endif
      pint = pe1(ks+1)

#ifdef NO_UKMO_HB
      do k=1,ks+1
         ak(k) = pe1(k)
         bk(k) = 0.
      enddo

      do k=ks+2,km+1
         bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
         ak(k) =  pe1(k) - bk(k) * pe1(km+1)
      enddo
      bk(km+1) = 1.
      ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
      do k=1,km+1
         eta(k) = pe1(k) / pe1(km+1)
      enddo
      ep =  eta(ks+1)
      es =  eta(km)
!     es =  1.
      alpha = (ep**2-2.*ep*es) / (es-ep)**2
      beta  = 2.*ep*es**2 / (es-ep)**2
      gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
      do k=1,ks+1
         ak(k) = eta(k)*1.e5
         bk(k) = 0.
      enddo

      do k=ks+2, km
         ak(k) = alpha*eta(k) + beta + gama/eta(k)
         ak(k) = ak(k)*1.e5
      enddo
         ak(km+1) = 0.

      do k=ks+2, km
         bk(k) = (pe1(k) - ak(k))/pe1(km+1)
      enddo
         bk(km+1) = 1.
#endif

      if ( verbose ) then
          tmp1 = ak(ks+1)
          do k=ks+1,km
             tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
          enddo
          write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
      endif

end subroutine mount_waves

subroutine var_dz(km, ak, bk, ptop, ks, pint, s_rate)
 integer, intent(in):: km
 real,    intent(in):: ptop
 real,    intent(in):: s_rate        !< between [1. 1.1]
 real,    intent(out):: ak(km+1), bk(km+1)
 real,    intent(inout):: pint
 integer, intent(out):: ks
! Local
 real, parameter:: p00 = 1.E5
 real, dimension(km+1):: ze, pe1, peln, eta
 real, dimension(km):: dz, s_fac, dlnp
 real ztop, t0, dz0, sum1, tmp1
 real ep, es, alpha, beta, gama
 integer  k

    pe1(1) = ptop
    peln(1) = log(pe1(1))
    pe1(km+1) = p00
    peln(km+1) = log(pe1(km+1))

    t0 = 270.
    ztop = con_rd/con_g*t0*(peln(km+1) - peln(1))

     s_fac(km  ) = 0.10
     s_fac(km-1) = 0.20
     s_fac(km-2) = 0.30
     s_fac(km-3) = 0.40
     s_fac(km-4) = 0.50
     s_fac(km-5) = 0.60
     s_fac(km-6) = 0.70
     s_fac(km-7) = 0.80
     s_fac(km-8) = 0.90
     s_fac(km-9) = 0.95
     s_fac(km-10) = 0.5*(s_fac(km-9) + s_rate)

     do k=km-11, 9, -1
        s_fac(k) = min(10.0, s_rate * s_fac(k+1) )
     enddo

     s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
     s_fac(7) = 1.1 *s_fac(8)
     s_fac(6) = 1.15*s_fac(7)
     s_fac(5) = 1.2 *s_fac(6)
     s_fac(4) = 1.3 *s_fac(5)
     s_fac(3) = 1.4 *s_fac(4)
     s_fac(2) = 1.5 *s_fac(3)
     s_fac(1) = 1.6 *s_fac(2)

     sum1 = 0.
     do k=1,km
        sum1 = sum1 + s_fac(k)
     enddo

     dz0 = ztop / sum1

     do k=1,km
        dz(k) = s_fac(k) * dz0
     enddo

     ze(km+1) = 0.
     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo

! Re-scale dz with the stretched ztop
     do k=1,km
        dz(k) = dz(k) * (ztop/ze(1))
     enddo

     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo
!     ze(1) = ztop

     if ( verbose ) write(*,*) 'var_dz: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
     call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 1)

! Given z --> p
     do k=1,km
         dz(k) = ze(k) - ze(k+1)
       dlnp(k) = con_g*dz(k) / (con_rd*t0)
     enddo
     do k=2,km
        peln(k) = peln(k-1) + dlnp(k-1)
         pe1(k) = exp(peln(k))
     enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
     ks = 0
     do k=2,km
        if ( pint < pe1(k)) then
             ks = k-1
             exit
        endif
     enddo
     if ( verbose ) then
        write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
        write(*,*) 'ptop =', ptop
     endif
     pint = pe1(ks+1)

#ifdef NO_UKMO_HB
     do k=1,ks+1
        ak(k) = pe1(k)
        bk(k) = 0.
     enddo

     do k=ks+2,km+1
        bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
        ak(k) =  pe1(k) - bk(k) * pe1(km+1)
     enddo
     bk(km+1) = 1.
     ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
     do k=1,km+1
        eta(k) = pe1(k) / pe1(km+1)
     enddo

     ep =  eta(ks+1)
     es =  eta(km)
!     es =  1.
     alpha = (ep**2-2.*ep*es) / (es-ep)**2
     beta  = 2.*ep*es**2 / (es-ep)**2
     gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
     do k=1,ks+1
        ak(k) = eta(k)*1.e5
        bk(k) = 0.
     enddo

     do k=ks+2, km
        ak(k) = alpha*eta(k) + beta + gama/eta(k)
        ak(k) = ak(k)*1.e5
     enddo
        ak(km+1) = 0.

     do k=ks+2, km
        bk(k) = (pe1(k) - ak(k))/pe1(km+1)
     enddo
        bk(km+1) = 1.
#endif

     if ( verbose ) then
         write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
         do k=1,km
            write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
         enddo
         tmp1 = ak(ks+1)
         do k=ks+1,km
            tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
         enddo
         write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
     endif
end subroutine var_dz

subroutine var_gfs(km, ak, bk, ptop, ks, pint, s_rate)
 integer, intent(in):: km
 real,    intent(in):: ptop
 real,    intent(in):: s_rate        !< between [1. 1.1]
 real,    intent(out):: ak(km+1), bk(km+1)
 real,    intent(inout):: pint
 integer, intent(out):: ks
! Local
 real, parameter:: p00 = 1.E5
 real, dimension(km+1):: ze, pe1, peln, eta
 real, dimension(km):: dz, s_fac, dlnp
 real ztop, t0, dz0, sum1, tmp1
 real ep, es, alpha, beta, gama
!---- Tunable parameters:
 integer:: k_inc = 25   !< number of layers from bottom up to near const dz region
 real:: s0 = 0.13 !< lowest layer stretch factor
!-----------------------
 real:: s_inc
 integer  k

    pe1(1) = ptop
    peln(1) = log(pe1(1))
    pe1(km+1) = p00
    peln(km+1) = log(pe1(km+1))

    t0 = 270.
    ztop = con_rd/con_g*t0*(peln(km+1) - peln(1))

     s_inc = (1.-s0) / real(k_inc)
     s_fac(km)  = s0

     do k=km-1, km-k_inc, -1
        s_fac(k)  = s_fac(k+1) + s_inc
     enddo

     do k=km-k_inc-1, 9, -1
        s_fac(k) = s_rate * s_fac(k+1)
     enddo
     s_fac(8) = 0.5*(1.1+s_rate)*s_fac(9)
     s_fac(7) = 1.10*s_fac(8)
     s_fac(6) = 1.15*s_fac(7)
     s_fac(5) = 1.20*s_fac(6)
     s_fac(4) = 1.26*s_fac(5)
     s_fac(3) = 1.33*s_fac(4)
     s_fac(2) = 1.41*s_fac(3)
     s_fac(1) = 1.60*s_fac(2)

     sum1 = 0.
     do k=1,km
        sum1 = sum1 + s_fac(k)
     enddo

     dz0 = ztop / sum1

     do k=1,km
        dz(k) = s_fac(k) * dz0
     enddo

     ze(km+1) = 0.
     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo

! Re-scale dz with the stretched ztop
     do k=1,km
        dz(k) = dz(k) * (ztop/ze(1))
     enddo

     do k=km,1,-1
        ze(k) = ze(k+1) + dz(k)
     enddo
!     ze(1) = ztop

     if ( verbose ) then
          write(*,*) 'var_gfs: computed model top (m)=', ztop*0.001, ' bottom/top dz=', dz(km), dz(1)
!          do k=1,km
!             write(*,*) k, s_fac(k)
!          enddo
     endif

!     call sm1_edge(1, 1, 1, 1, km, 1, 1, ze, 3)

! Given z --> p
     do k=1,km
         dz(k) = ze(k) - ze(k+1)
       dlnp(k) = con_g*dz(k) / (con_rd*t0)
     enddo
     do k=2,km
        peln(k) = peln(k-1) + dlnp(k-1)
         pe1(k) = exp(peln(k))
     enddo

! Pe(k) = ak(k) + bk(k) * PS
! Locate pint and KS
     ks = 0
     do k=2,km
        if ( pint < pe1(k)) then
             ks = k-1
             exit
        endif
     enddo
     if ( verbose ) then
        write(*,*) 'For (input) PINT=', 0.01*pint, ' KS=', ks, 'pint(computed)=', 0.01*pe1(ks+1)
        write(*,*) 'ptop =', ptop
     endif
     pint = pe1(ks+1)

#ifdef NO_UKMO_HB
     do k=1,ks+1
        ak(k) = pe1(k)
        bk(k) = 0.
     enddo

     do k=ks+2,km+1
        bk(k) = (pe1(k) - pint) / (pe1(km+1)-pint)  ! bk == sigma
        ak(k) =  pe1(k) - bk(k) * pe1(km+1)
     enddo
     bk(km+1) = 1.
     ak(km+1) = 0.
#else
! Problematic for non-hydrostatic
     do k=1,km+1
        eta(k) = pe1(k) / pe1(km+1)
     enddo

     ep =  eta(ks+1)
     es =  eta(km)
!     es =  1.
     alpha = (ep**2-2.*ep*es) / (es-ep)**2
     beta  = 2.*ep*es**2 / (es-ep)**2
     gama = -(ep*es)**2 / (es-ep)**2

! Pure pressure:
     do k=1,ks+1
        ak(k) = eta(k)*1.e5
        bk(k) = 0.
     enddo

     do k=ks+2, km
        ak(k) = alpha*eta(k) + beta + gama/eta(k)
        ak(k) = ak(k)*1.e5
     enddo
        ak(km+1) = 0.

     do k=ks+2, km
        bk(k) = (pe1(k) - ak(k))/pe1(km+1)
     enddo
        bk(km+1) = 1.
#endif

     if ( verbose ) then
         write(*,*) 'KS=', ks, 'PINT (mb)=', pint/100.
         do k=1,km
            write(*,*) k, 0.5*(pe1(k)+pe1(k+1))/100., dz(k)
         enddo
         tmp1 = ak(ks+1)
         do k=ks+1,km
            tmp1 = max(tmp1, (ak(k)-ak(k+1))/max(1.E-5, (bk(k+1)-bk(k))) )
         enddo
         write(*,*) 'Hybrid Sigma-P: minimum allowable surface pressure (hpa)=', tmp1/100.
     endif
end subroutine var_gfs

subroutine sm1_edge(is, ie, js, je, km, i, j, ze, ntimes)
  integer, intent(in):: is, ie, js, je, km
  integer, intent(in):: ntimes, i, j
  real, intent(inout):: ze(is:ie,js:je,km+1)
 ! local:
  real, parameter:: df = 0.25
  real dz(km)
  real flux(km+1)
  integer k, n, k1, k2

      k2 = km-1
      do k=1,km
         dz(k) = ze(i,j,k+1) - ze(i,j,k)
      enddo

   do n=1,ntimes
      k1 = 2 + (ntimes-n)

      flux(k1  ) = 0.
      flux(k2+1) = 0.
      do k=k1+1,k2
         flux(k) = df*(dz(k) - dz(k-1))
      enddo

      do k=k1,k2
         dz(k) = dz(k) - flux(k) + flux(k+1)
      enddo
   enddo

   do k=km,1,-1
      ze(i,j,k) = ze(i,j,k+1) - dz(k)
   enddo

end subroutine sm1_edge

subroutine check_eta_levels(ak, bk)
 real,    intent(in) :: ak(:)
 real,    intent(in) :: bk(:)
 !--- local variables
 real :: ph1, tmp
 integer :: nlev, k
 logical :: monotonic

 nlev = size(ak(:))

 monotonic = .true.
 ph1 = ak(1)
 do k=2,nlev
    tmp = ak(k) + bk(k)*1000.E2
    if (tmp <= ph1) then
      monotonic = .false.
      exit
    endif
    ph1 = tmp
 enddo

 if (.not. monotonic) then
   if (verbose) then
      write(*, '(A4, A13, A13, A11)') 'klev', 'ak', 'bk', 'p_ref'
      do k=1,nlev
         write(*,'(I4, F13.5, F13.5, F11.2)') k, ak(k), bk(k), ak(k) + bk(k)*1000.E2
      enddo
   endif
   write(*,*) 'FV3 check_eta_levels: ak/bk pairs do not provide a monotonic vertical coordinate'
   error stop
 endif

end subroutine check_eta_levels

!> This subroutine calculates the pressure and exner function at grid centers and interface levels given a surface pressure and interface-level GFS grid coefficients.
!! This subroutine should be called to update the pressures of the model levels as the surface pressure of the column changes.
subroutine calc_pres_exner_geopotential(time_level, scm_state)
  use scm_type_defs, only: scm_state_type

  integer, intent(in) :: time_level
  type(scm_state_type), intent(inout) :: scm_state

  real(kind=dp)               :: pres_sfc_inv, tem, dgeopotential_lower_half, dgeopotential_upper_half
  integer               :: i,k

  !> - Calculate interface pressures, sigma, and exner function.
  do i=1, scm_state%n_cols
    pres_sfc_inv = 1.0/scm_state%pres_surf(i)
    do k=1, scm_state%n_levels+1
      scm_state%pres_i(i,k) = scm_state%a_k(k) + scm_state%b_k(k)*scm_state%pres_surf(i)
      scm_state%si(i,k) = scm_state%a_k(k)*pres_sfc_inv + scm_state%b_k(k)
      scm_state%exner_i(i,k) = (scm_state%pres_i(i,k)*1.0E-5)**con_rocp
    end do
  end do

  !> - Calculate layer center pressures, sigma, and exner function.
  do i=1, scm_state%n_cols
    do k=1, scm_state%n_levels
      scm_state%pres_l(i,k) = ((1.0/(con_rocp+1.0))*&
        (scm_state%pres_i(i,k)**(con_rocp+1.0) - scm_state%pres_i(i,k+1)**(con_rocp+1.0))/ &
        (scm_state%pres_i(i,k) - scm_state%pres_i(i,k+1)))**(1.0/con_rocp)
      scm_state%sl(i,k) = 0.5*(scm_state%si(i,k) + scm_state%si(i,k+1))
      scm_state%exner_l(i,k) = (scm_state%pres_l(i,k)*1.0E-5)**con_rocp
    end do
  end do

  do i=1, scm_state%n_cols
    scm_state%geopotential_i(i,1) = 0.0
  end do

  do i=1, scm_state%n_cols
    do k=1, scm_state%n_levels
      tem = con_cp*scm_state%state_T(i,k,time_level)*&
        (1.0 + con_fvirt*max(scm_state%state_tracer(i,k,scm_state%water_vapor_index,time_level), 0.0))/scm_state%exner_l(i,k)
      dgeopotential_lower_half = (scm_state%exner_i(i,k) - scm_state%exner_l(i,k))*tem
      dgeopotential_upper_half = (scm_state%exner_l(i,k) - scm_state%exner_i(i,k+1))*tem
      scm_state%geopotential_l(i,k) = scm_state%geopotential_i(i,k) + dgeopotential_lower_half
      scm_state%geopotential_i(i,k+1) = scm_state%geopotential_l(i,k) + dgeopotential_upper_half
    end do
  end do

end subroutine calc_pres_exner_geopotential

subroutine calc_geopotential(time_level, scm_state)
  use scm_type_defs, only: scm_state_type

  integer, intent(in) :: time_level
  type(scm_state_type), intent(inout) :: scm_state

  integer i,k
  real(kind=dp) :: tem, dgeopotential_lower_half, dgeopotential_upper_half

  do i=1, scm_state%n_cols
    scm_state%geopotential_i(i,1) = 0.0
  end do

  do i=1, scm_state%n_cols
    do k=1, scm_state%n_levels
      tem = con_cp*scm_state%state_T(i,k,time_level)*(1.0 + &
        con_fvirt*max(scm_state%state_tracer(i,k,scm_state%water_vapor_index,time_level), 0.0))/scm_state%exner_l(i,k)
      dgeopotential_lower_half = (scm_state%exner_i(i,k) - scm_state%exner_l(i,k))*tem
      dgeopotential_upper_half = (scm_state%exner_l(i,k) - scm_state%exner_i(i,k+1))*tem
      scm_state%geopotential_l(i,k) = scm_state%geopotential_i(i,k) + dgeopotential_lower_half
      scm_state%geopotential_i(i,k+1) = scm_state%geopotential_l(i,k) + dgeopotential_upper_half
    end do
  end do

end subroutine calc_geopotential
!> @}
!> @}
end module scm_vgrid
