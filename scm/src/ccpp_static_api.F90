
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated API for the CCPP static build
!!
!
module ccpp_static_api

   use ccpp_HAFS_v0_hwrf_cap, only: HAFS_v0_hwrf_tsinit_cap
   use ccpp_HAFS_v0_hwrf_cap, only: HAFS_v0_hwrf_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_cap, only: HAFS_v0_hwrf_init_cap
   use ccpp_HAFS_v0_hwrf_cap, only: HAFS_v0_hwrf_run_cap
   use ccpp_HAFS_v0_hwrf_cap, only: HAFS_v0_hwrf_final_cap
   use ccpp_HAFS_v0_hwrf_time_vary_cap, only: HAFS_v0_hwrf_time_vary_tsinit_cap
   use ccpp_HAFS_v0_hwrf_time_vary_cap, only: HAFS_v0_hwrf_time_vary_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_time_vary_cap, only: HAFS_v0_hwrf_time_vary_init_cap
   use ccpp_HAFS_v0_hwrf_time_vary_cap, only: HAFS_v0_hwrf_time_vary_run_cap
   use ccpp_HAFS_v0_hwrf_time_vary_cap, only: HAFS_v0_hwrf_time_vary_final_cap
   use ccpp_HAFS_v0_hwrf_radiation_cap, only: HAFS_v0_hwrf_radiation_tsinit_cap
   use ccpp_HAFS_v0_hwrf_radiation_cap, only: HAFS_v0_hwrf_radiation_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_radiation_cap, only: HAFS_v0_hwrf_radiation_init_cap
   use ccpp_HAFS_v0_hwrf_radiation_cap, only: HAFS_v0_hwrf_radiation_run_cap
   use ccpp_HAFS_v0_hwrf_radiation_cap, only: HAFS_v0_hwrf_radiation_final_cap
   use ccpp_HAFS_v0_hwrf_physics_cap, only: HAFS_v0_hwrf_physics_tsinit_cap
   use ccpp_HAFS_v0_hwrf_physics_cap, only: HAFS_v0_hwrf_physics_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_physics_cap, only: HAFS_v0_hwrf_physics_init_cap
   use ccpp_HAFS_v0_hwrf_physics_cap, only: HAFS_v0_hwrf_physics_run_cap
   use ccpp_HAFS_v0_hwrf_physics_cap, only: HAFS_v0_hwrf_physics_final_cap
   use ccpp_HAFS_v0_hwrf_ps_cap, only: HAFS_v0_hwrf_ps_tsinit_cap
   use ccpp_HAFS_v0_hwrf_ps_cap, only: HAFS_v0_hwrf_ps_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_ps_cap, only: HAFS_v0_hwrf_ps_init_cap
   use ccpp_HAFS_v0_hwrf_ps_cap, only: HAFS_v0_hwrf_ps_run_cap
   use ccpp_HAFS_v0_hwrf_ps_cap, only: HAFS_v0_hwrf_ps_final_cap
   use ccpp_HAFS_v0_hwrf_ps_time_vary_cap, only: HAFS_v0_hwrf_ps_time_vary_tsinit_cap
   use ccpp_HAFS_v0_hwrf_ps_time_vary_cap, only: HAFS_v0_hwrf_ps_time_vary_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_ps_time_vary_cap, only: HAFS_v0_hwrf_ps_time_vary_init_cap
   use ccpp_HAFS_v0_hwrf_ps_time_vary_cap, only: HAFS_v0_hwrf_ps_time_vary_run_cap
   use ccpp_HAFS_v0_hwrf_ps_time_vary_cap, only: HAFS_v0_hwrf_ps_time_vary_final_cap
   use ccpp_HAFS_v0_hwrf_ps_radiation_cap, only: HAFS_v0_hwrf_ps_radiation_tsinit_cap
   use ccpp_HAFS_v0_hwrf_ps_radiation_cap, only: HAFS_v0_hwrf_ps_radiation_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_ps_radiation_cap, only: HAFS_v0_hwrf_ps_radiation_init_cap
   use ccpp_HAFS_v0_hwrf_ps_radiation_cap, only: HAFS_v0_hwrf_ps_radiation_run_cap
   use ccpp_HAFS_v0_hwrf_ps_radiation_cap, only: HAFS_v0_hwrf_ps_radiation_final_cap
   use ccpp_HAFS_v0_hwrf_ps_physics_cap, only: HAFS_v0_hwrf_ps_physics_tsinit_cap
   use ccpp_HAFS_v0_hwrf_ps_physics_cap, only: HAFS_v0_hwrf_ps_physics_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_ps_physics_cap, only: HAFS_v0_hwrf_ps_physics_init_cap
   use ccpp_HAFS_v0_hwrf_ps_physics_cap, only: HAFS_v0_hwrf_ps_physics_run_cap
   use ccpp_HAFS_v0_hwrf_ps_physics_cap, only: HAFS_v0_hwrf_ps_physics_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_cap, only: HAFS_v0_hwrf_thompson_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_cap, only: HAFS_v0_hwrf_thompson_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_cap, only: HAFS_v0_hwrf_thompson_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_cap, only: HAFS_v0_hwrf_thompson_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_cap, only: HAFS_v0_hwrf_thompson_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_time_vary_cap, only: HAFS_v0_hwrf_thompson_time_vary_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_time_vary_cap, only: HAFS_v0_hwrf_thompson_time_vary_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_time_vary_cap, only: HAFS_v0_hwrf_thompson_time_vary_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_time_vary_cap, only: HAFS_v0_hwrf_thompson_time_vary_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_time_vary_cap, only: HAFS_v0_hwrf_thompson_time_vary_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_radiation_cap, only: HAFS_v0_hwrf_thompson_radiation_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_radiation_cap, only: HAFS_v0_hwrf_thompson_radiation_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_radiation_cap, only: HAFS_v0_hwrf_thompson_radiation_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_radiation_cap, only: HAFS_v0_hwrf_thompson_radiation_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_radiation_cap, only: HAFS_v0_hwrf_thompson_radiation_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_physics_cap, only: HAFS_v0_hwrf_thompson_physics_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_physics_cap, only: HAFS_v0_hwrf_thompson_physics_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_physics_cap, only: HAFS_v0_hwrf_thompson_physics_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_physics_cap, only: HAFS_v0_hwrf_thompson_physics_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_physics_cap, only: HAFS_v0_hwrf_thompson_physics_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_cap, only: HAFS_v0_hwrf_thompson_ps_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_cap, only: HAFS_v0_hwrf_thompson_ps_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_cap, only: HAFS_v0_hwrf_thompson_ps_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_cap, only: HAFS_v0_hwrf_thompson_ps_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_cap, only: HAFS_v0_hwrf_thompson_ps_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_time_vary_cap, only: HAFS_v0_hwrf_thompson_ps_time_vary_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_time_vary_cap, only: HAFS_v0_hwrf_thompson_ps_time_vary_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_time_vary_cap, only: HAFS_v0_hwrf_thompson_ps_time_vary_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_time_vary_cap, only: HAFS_v0_hwrf_thompson_ps_time_vary_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_time_vary_cap, only: HAFS_v0_hwrf_thompson_ps_time_vary_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_radiation_cap, only: HAFS_v0_hwrf_thompson_ps_radiation_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_radiation_cap, only: HAFS_v0_hwrf_thompson_ps_radiation_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_radiation_cap, only: HAFS_v0_hwrf_thompson_ps_radiation_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_radiation_cap, only: HAFS_v0_hwrf_thompson_ps_radiation_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_radiation_cap, only: HAFS_v0_hwrf_thompson_ps_radiation_final_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_physics_cap, only: HAFS_v0_hwrf_thompson_ps_physics_tsinit_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_physics_cap, only: HAFS_v0_hwrf_thompson_ps_physics_tsfinal_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_physics_cap, only: HAFS_v0_hwrf_thompson_ps_physics_init_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_physics_cap, only: HAFS_v0_hwrf_thompson_ps_physics_run_cap
   use ccpp_HAFS_v0_hwrf_thompson_ps_physics_cap, only: HAFS_v0_hwrf_thompson_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_cap, only: SCM_GFS_v15p2_tsinit_cap
   use ccpp_SCM_GFS_v15p2_cap, only: SCM_GFS_v15p2_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_cap, only: SCM_GFS_v15p2_init_cap
   use ccpp_SCM_GFS_v15p2_cap, only: SCM_GFS_v15p2_run_cap
   use ccpp_SCM_GFS_v15p2_cap, only: SCM_GFS_v15p2_final_cap
   use ccpp_SCM_GFS_v15p2_time_vary_cap, only: SCM_GFS_v15p2_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_time_vary_cap, only: SCM_GFS_v15p2_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_time_vary_cap, only: SCM_GFS_v15p2_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_time_vary_cap, only: SCM_GFS_v15p2_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_time_vary_cap, only: SCM_GFS_v15p2_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_radiation_cap, only: SCM_GFS_v15p2_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_radiation_cap, only: SCM_GFS_v15p2_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_radiation_cap, only: SCM_GFS_v15p2_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_radiation_cap, only: SCM_GFS_v15p2_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_radiation_cap, only: SCM_GFS_v15p2_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_physics_cap, only: SCM_GFS_v15p2_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_physics_cap, only: SCM_GFS_v15p2_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_physics_cap, only: SCM_GFS_v15p2_physics_init_cap
   use ccpp_SCM_GFS_v15p2_physics_cap, only: SCM_GFS_v15p2_physics_run_cap
   use ccpp_SCM_GFS_v15p2_physics_cap, only: SCM_GFS_v15p2_physics_final_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_cap, only: SCM_GFS_v15p2_ACM_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_cap, only: SCM_GFS_v15p2_ACM_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_cap, only: SCM_GFS_v15p2_ACM_ps_init_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_cap, only: SCM_GFS_v15p2_ACM_ps_run_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_cap, only: SCM_GFS_v15p2_ACM_ps_final_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_time_vary_cap, only: SCM_GFS_v15p2_ACM_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_time_vary_cap, only: SCM_GFS_v15p2_ACM_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_time_vary_cap, only: SCM_GFS_v15p2_ACM_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_time_vary_cap, only: SCM_GFS_v15p2_ACM_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_time_vary_cap, only: SCM_GFS_v15p2_ACM_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_radiation_cap, only: SCM_GFS_v15p2_ACM_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_radiation_cap, only: SCM_GFS_v15p2_ACM_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_radiation_cap, only: SCM_GFS_v15p2_ACM_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_radiation_cap, only: SCM_GFS_v15p2_ACM_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_radiation_cap, only: SCM_GFS_v15p2_ACM_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_physics_cap, only: SCM_GFS_v15p2_ACM_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_physics_cap, only: SCM_GFS_v15p2_ACM_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_physics_cap, only: SCM_GFS_v15p2_ACM_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_physics_cap, only: SCM_GFS_v15p2_ACM_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_ACM_ps_physics_cap, only: SCM_GFS_v15p2_ACM_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_FA_cap, only: SCM_GFS_v15p2_FA_tsinit_cap
   use ccpp_SCM_GFS_v15p2_FA_cap, only: SCM_GFS_v15p2_FA_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_FA_cap, only: SCM_GFS_v15p2_FA_init_cap
   use ccpp_SCM_GFS_v15p2_FA_cap, only: SCM_GFS_v15p2_FA_run_cap
   use ccpp_SCM_GFS_v15p2_FA_cap, only: SCM_GFS_v15p2_FA_final_cap
   use ccpp_SCM_GFS_v15p2_FA_time_vary_cap, only: SCM_GFS_v15p2_FA_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_FA_time_vary_cap, only: SCM_GFS_v15p2_FA_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_FA_time_vary_cap, only: SCM_GFS_v15p2_FA_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_FA_time_vary_cap, only: SCM_GFS_v15p2_FA_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_FA_time_vary_cap, only: SCM_GFS_v15p2_FA_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_FA_radiation_cap, only: SCM_GFS_v15p2_FA_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_FA_radiation_cap, only: SCM_GFS_v15p2_FA_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_FA_radiation_cap, only: SCM_GFS_v15p2_FA_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_FA_radiation_cap, only: SCM_GFS_v15p2_FA_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_FA_radiation_cap, only: SCM_GFS_v15p2_FA_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_FA_physics_cap, only: SCM_GFS_v15p2_FA_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_FA_physics_cap, only: SCM_GFS_v15p2_FA_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_FA_physics_cap, only: SCM_GFS_v15p2_FA_physics_init_cap
   use ccpp_SCM_GFS_v15p2_FA_physics_cap, only: SCM_GFS_v15p2_FA_physics_run_cap
   use ccpp_SCM_GFS_v15p2_FA_physics_cap, only: SCM_GFS_v15p2_FA_physics_final_cap
   use ccpp_SCM_GFS_v15p2_MYJ_cap, only: SCM_GFS_v15p2_MYJ_tsinit_cap
   use ccpp_SCM_GFS_v15p2_MYJ_cap, only: SCM_GFS_v15p2_MYJ_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_MYJ_cap, only: SCM_GFS_v15p2_MYJ_init_cap
   use ccpp_SCM_GFS_v15p2_MYJ_cap, only: SCM_GFS_v15p2_MYJ_run_cap
   use ccpp_SCM_GFS_v15p2_MYJ_cap, only: SCM_GFS_v15p2_MYJ_final_cap
   use ccpp_SCM_GFS_v15p2_MYJ_time_vary_cap, only: SCM_GFS_v15p2_MYJ_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_MYJ_time_vary_cap, only: SCM_GFS_v15p2_MYJ_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_MYJ_time_vary_cap, only: SCM_GFS_v15p2_MYJ_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_MYJ_time_vary_cap, only: SCM_GFS_v15p2_MYJ_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_MYJ_time_vary_cap, only: SCM_GFS_v15p2_MYJ_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_MYJ_radiation_cap, only: SCM_GFS_v15p2_MYJ_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_MYJ_radiation_cap, only: SCM_GFS_v15p2_MYJ_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_MYJ_radiation_cap, only: SCM_GFS_v15p2_MYJ_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_MYJ_radiation_cap, only: SCM_GFS_v15p2_MYJ_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_MYJ_radiation_cap, only: SCM_GFS_v15p2_MYJ_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_MYJ_physics_cap, only: SCM_GFS_v15p2_MYJ_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_MYJ_physics_cap, only: SCM_GFS_v15p2_MYJ_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_MYJ_physics_cap, only: SCM_GFS_v15p2_MYJ_physics_init_cap
   use ccpp_SCM_GFS_v15p2_MYJ_physics_cap, only: SCM_GFS_v15p2_MYJ_physics_run_cap
   use ccpp_SCM_GFS_v15p2_MYJ_physics_cap, only: SCM_GFS_v15p2_MYJ_physics_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_cap, only: SCM_GFS_v15p2_RRTMGP_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_cap, only: SCM_GFS_v15p2_RRTMGP_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_cap, only: SCM_GFS_v15p2_RRTMGP_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_cap, only: SCM_GFS_v15p2_RRTMGP_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_cap, only: SCM_GFS_v15p2_RRTMGP_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_physics_cap, only: SCM_GFS_v15p2_RRTMGP_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_physics_cap, only: SCM_GFS_v15p2_RRTMGP_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_physics_cap, only: SCM_GFS_v15p2_RRTMGP_physics_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_physics_cap, only: SCM_GFS_v15p2_RRTMGP_physics_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_physics_cap, only: SCM_GFS_v15p2_RRTMGP_physics_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_cap, only: SCM_GFS_v15p2_RRTMGP_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_cap, only: SCM_GFS_v15p2_RRTMGP_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_cap, only: SCM_GFS_v15p2_RRTMGP_ps_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_cap, only: SCM_GFS_v15p2_RRTMGP_ps_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_cap, only: SCM_GFS_v15p2_RRTMGP_ps_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v15p2_RRTMGP_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_radiation_cap, only: SCM_GFS_v15p2_RRTMGP_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_physics_cap, only: SCM_GFS_v15p2_RRTMGP_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_physics_cap, only: SCM_GFS_v15p2_RRTMGP_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_physics_cap, only: SCM_GFS_v15p2_RRTMGP_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_physics_cap, only: SCM_GFS_v15p2_RRTMGP_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_RRTMGP_ps_physics_cap, only: SCM_GFS_v15p2_RRTMGP_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_cap, only: SCM_GFS_v15p2_YSU_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_cap, only: SCM_GFS_v15p2_YSU_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_cap, only: SCM_GFS_v15p2_YSU_ps_init_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_cap, only: SCM_GFS_v15p2_YSU_ps_run_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_cap, only: SCM_GFS_v15p2_YSU_ps_final_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_time_vary_cap, only: SCM_GFS_v15p2_YSU_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_time_vary_cap, only: SCM_GFS_v15p2_YSU_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_time_vary_cap, only: SCM_GFS_v15p2_YSU_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_time_vary_cap, only: SCM_GFS_v15p2_YSU_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_time_vary_cap, only: SCM_GFS_v15p2_YSU_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_radiation_cap, only: SCM_GFS_v15p2_YSU_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_radiation_cap, only: SCM_GFS_v15p2_YSU_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_radiation_cap, only: SCM_GFS_v15p2_YSU_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_radiation_cap, only: SCM_GFS_v15p2_YSU_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_radiation_cap, only: SCM_GFS_v15p2_YSU_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_physics_cap, only: SCM_GFS_v15p2_YSU_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_physics_cap, only: SCM_GFS_v15p2_YSU_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_physics_cap, only: SCM_GFS_v15p2_YSU_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_physics_cap, only: SCM_GFS_v15p2_YSU_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_YSU_ps_physics_cap, only: SCM_GFS_v15p2_YSU_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_cap, only: SCM_GFS_v15p2_no_nsst_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_cap, only: SCM_GFS_v15p2_no_nsst_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_cap, only: SCM_GFS_v15p2_no_nsst_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_cap, only: SCM_GFS_v15p2_no_nsst_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_cap, only: SCM_GFS_v15p2_no_nsst_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_radiation_cap, only: SCM_GFS_v15p2_no_nsst_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_radiation_cap, only: SCM_GFS_v15p2_no_nsst_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_radiation_cap, only: SCM_GFS_v15p2_no_nsst_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_radiation_cap, only: SCM_GFS_v15p2_no_nsst_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_radiation_cap, only: SCM_GFS_v15p2_no_nsst_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_physics_cap, only: SCM_GFS_v15p2_no_nsst_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_physics_cap, only: SCM_GFS_v15p2_no_nsst_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_physics_cap, only: SCM_GFS_v15p2_no_nsst_physics_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_physics_cap, only: SCM_GFS_v15p2_no_nsst_physics_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_physics_cap, only: SCM_GFS_v15p2_no_nsst_physics_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_cap, only: SCM_GFS_v15p2_no_nsst_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_cap, only: SCM_GFS_v15p2_no_nsst_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_cap, only: SCM_GFS_v15p2_no_nsst_ps_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_cap, only: SCM_GFS_v15p2_no_nsst_ps_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_cap, only: SCM_GFS_v15p2_no_nsst_ps_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_time_vary_cap, only: SCM_GFS_v15p2_no_nsst_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_radiation_cap, only: SCM_GFS_v15p2_no_nsst_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_radiation_cap, only: SCM_GFS_v15p2_no_nsst_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_radiation_cap, only: SCM_GFS_v15p2_no_nsst_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_radiation_cap, only: SCM_GFS_v15p2_no_nsst_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_radiation_cap, only: SCM_GFS_v15p2_no_nsst_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_physics_cap, only: SCM_GFS_v15p2_no_nsst_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_physics_cap, only: SCM_GFS_v15p2_no_nsst_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_physics_cap, only: SCM_GFS_v15p2_no_nsst_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_physics_cap, only: SCM_GFS_v15p2_no_nsst_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_no_nsst_ps_physics_cap, only: SCM_GFS_v15p2_no_nsst_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_noahmp_cap, only: SCM_GFS_v15p2_noahmp_tsinit_cap
   use ccpp_SCM_GFS_v15p2_noahmp_cap, only: SCM_GFS_v15p2_noahmp_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_noahmp_cap, only: SCM_GFS_v15p2_noahmp_init_cap
   use ccpp_SCM_GFS_v15p2_noahmp_cap, only: SCM_GFS_v15p2_noahmp_run_cap
   use ccpp_SCM_GFS_v15p2_noahmp_cap, only: SCM_GFS_v15p2_noahmp_final_cap
   use ccpp_SCM_GFS_v15p2_noahmp_time_vary_cap, only: SCM_GFS_v15p2_noahmp_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_noahmp_time_vary_cap, only: SCM_GFS_v15p2_noahmp_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_noahmp_time_vary_cap, only: SCM_GFS_v15p2_noahmp_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_noahmp_time_vary_cap, only: SCM_GFS_v15p2_noahmp_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_noahmp_time_vary_cap, only: SCM_GFS_v15p2_noahmp_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_noahmp_radiation_cap, only: SCM_GFS_v15p2_noahmp_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_noahmp_radiation_cap, only: SCM_GFS_v15p2_noahmp_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_noahmp_radiation_cap, only: SCM_GFS_v15p2_noahmp_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_noahmp_radiation_cap, only: SCM_GFS_v15p2_noahmp_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_noahmp_radiation_cap, only: SCM_GFS_v15p2_noahmp_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_noahmp_physics_cap, only: SCM_GFS_v15p2_noahmp_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_noahmp_physics_cap, only: SCM_GFS_v15p2_noahmp_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_noahmp_physics_cap, only: SCM_GFS_v15p2_noahmp_physics_init_cap
   use ccpp_SCM_GFS_v15p2_noahmp_physics_cap, only: SCM_GFS_v15p2_noahmp_physics_run_cap
   use ccpp_SCM_GFS_v15p2_noahmp_physics_cap, only: SCM_GFS_v15p2_noahmp_physics_final_cap
   use ccpp_SCM_GFS_v15p2_ps_cap, only: SCM_GFS_v15p2_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ps_cap, only: SCM_GFS_v15p2_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ps_cap, only: SCM_GFS_v15p2_ps_init_cap
   use ccpp_SCM_GFS_v15p2_ps_cap, only: SCM_GFS_v15p2_ps_run_cap
   use ccpp_SCM_GFS_v15p2_ps_cap, only: SCM_GFS_v15p2_ps_final_cap
   use ccpp_SCM_GFS_v15p2_ps_time_vary_cap, only: SCM_GFS_v15p2_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ps_time_vary_cap, only: SCM_GFS_v15p2_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ps_time_vary_cap, only: SCM_GFS_v15p2_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_ps_time_vary_cap, only: SCM_GFS_v15p2_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_ps_time_vary_cap, only: SCM_GFS_v15p2_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_ps_radiation_cap, only: SCM_GFS_v15p2_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ps_radiation_cap, only: SCM_GFS_v15p2_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ps_radiation_cap, only: SCM_GFS_v15p2_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_ps_radiation_cap, only: SCM_GFS_v15p2_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_ps_radiation_cap, only: SCM_GFS_v15p2_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_ps_physics_cap, only: SCM_GFS_v15p2_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_ps_physics_cap, only: SCM_GFS_v15p2_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_ps_physics_cap, only: SCM_GFS_v15p2_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_ps_physics_cap, only: SCM_GFS_v15p2_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_ps_physics_cap, only: SCM_GFS_v15p2_ps_physics_final_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_cap, only: SCM_GFS_v15p2_saYSU_ps_tsinit_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_cap, only: SCM_GFS_v15p2_saYSU_ps_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_cap, only: SCM_GFS_v15p2_saYSU_ps_init_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_cap, only: SCM_GFS_v15p2_saYSU_ps_run_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_cap, only: SCM_GFS_v15p2_saYSU_ps_final_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_time_vary_cap, only: SCM_GFS_v15p2_saYSU_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_time_vary_cap, only: SCM_GFS_v15p2_saYSU_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_time_vary_cap, only: SCM_GFS_v15p2_saYSU_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_time_vary_cap, only: SCM_GFS_v15p2_saYSU_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_time_vary_cap, only: SCM_GFS_v15p2_saYSU_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_radiation_cap, only: SCM_GFS_v15p2_saYSU_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_radiation_cap, only: SCM_GFS_v15p2_saYSU_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_radiation_cap, only: SCM_GFS_v15p2_saYSU_ps_radiation_init_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_radiation_cap, only: SCM_GFS_v15p2_saYSU_ps_radiation_run_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_radiation_cap, only: SCM_GFS_v15p2_saYSU_ps_radiation_final_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_physics_cap, only: SCM_GFS_v15p2_saYSU_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_physics_cap, only: SCM_GFS_v15p2_saYSU_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_physics_cap, only: SCM_GFS_v15p2_saYSU_ps_physics_init_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_physics_cap, only: SCM_GFS_v15p2_saYSU_ps_physics_run_cap
   use ccpp_SCM_GFS_v15p2_saYSU_ps_physics_cap, only: SCM_GFS_v15p2_saYSU_ps_physics_final_cap
   use ccpp_SCM_GFS_v16_cap, only: SCM_GFS_v16_tsinit_cap
   use ccpp_SCM_GFS_v16_cap, only: SCM_GFS_v16_tsfinal_cap
   use ccpp_SCM_GFS_v16_cap, only: SCM_GFS_v16_init_cap
   use ccpp_SCM_GFS_v16_cap, only: SCM_GFS_v16_run_cap
   use ccpp_SCM_GFS_v16_cap, only: SCM_GFS_v16_final_cap
   use ccpp_SCM_GFS_v16_time_vary_cap, only: SCM_GFS_v16_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_time_vary_cap, only: SCM_GFS_v16_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_time_vary_cap, only: SCM_GFS_v16_time_vary_init_cap
   use ccpp_SCM_GFS_v16_time_vary_cap, only: SCM_GFS_v16_time_vary_run_cap
   use ccpp_SCM_GFS_v16_time_vary_cap, only: SCM_GFS_v16_time_vary_final_cap
   use ccpp_SCM_GFS_v16_radiation_cap, only: SCM_GFS_v16_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_radiation_cap, only: SCM_GFS_v16_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_radiation_cap, only: SCM_GFS_v16_radiation_init_cap
   use ccpp_SCM_GFS_v16_radiation_cap, only: SCM_GFS_v16_radiation_run_cap
   use ccpp_SCM_GFS_v16_radiation_cap, only: SCM_GFS_v16_radiation_final_cap
   use ccpp_SCM_GFS_v16_physics_cap, only: SCM_GFS_v16_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_physics_cap, only: SCM_GFS_v16_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_physics_cap, only: SCM_GFS_v16_physics_init_cap
   use ccpp_SCM_GFS_v16_physics_cap, only: SCM_GFS_v16_physics_run_cap
   use ccpp_SCM_GFS_v16_physics_cap, only: SCM_GFS_v16_physics_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_cap, only: SCM_GFS_v16_RRTMGP_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_cap, only: SCM_GFS_v16_RRTMGP_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_cap, only: SCM_GFS_v16_RRTMGP_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_cap, only: SCM_GFS_v16_RRTMGP_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_cap, only: SCM_GFS_v16_RRTMGP_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_time_vary_cap, only: SCM_GFS_v16_RRTMGP_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_time_vary_cap, only: SCM_GFS_v16_RRTMGP_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_time_vary_cap, only: SCM_GFS_v16_RRTMGP_time_vary_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_time_vary_cap, only: SCM_GFS_v16_RRTMGP_time_vary_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_time_vary_cap, only: SCM_GFS_v16_RRTMGP_time_vary_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_radiation_cap, only: SCM_GFS_v16_RRTMGP_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_radiation_cap, only: SCM_GFS_v16_RRTMGP_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_radiation_cap, only: SCM_GFS_v16_RRTMGP_radiation_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_radiation_cap, only: SCM_GFS_v16_RRTMGP_radiation_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_radiation_cap, only: SCM_GFS_v16_RRTMGP_radiation_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_physics_cap, only: SCM_GFS_v16_RRTMGP_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_physics_cap, only: SCM_GFS_v16_RRTMGP_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_physics_cap, only: SCM_GFS_v16_RRTMGP_physics_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_physics_cap, only: SCM_GFS_v16_RRTMGP_physics_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_physics_cap, only: SCM_GFS_v16_RRTMGP_physics_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_cap, only: SCM_GFS_v16_RRTMGP_ps_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_cap, only: SCM_GFS_v16_RRTMGP_ps_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_cap, only: SCM_GFS_v16_RRTMGP_ps_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_cap, only: SCM_GFS_v16_RRTMGP_ps_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_cap, only: SCM_GFS_v16_RRTMGP_ps_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v16_RRTMGP_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v16_RRTMGP_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v16_RRTMGP_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v16_RRTMGP_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_time_vary_cap, only: SCM_GFS_v16_RRTMGP_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_radiation_cap, only: SCM_GFS_v16_RRTMGP_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_radiation_cap, only: SCM_GFS_v16_RRTMGP_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_radiation_cap, only: SCM_GFS_v16_RRTMGP_ps_radiation_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_radiation_cap, only: SCM_GFS_v16_RRTMGP_ps_radiation_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_radiation_cap, only: SCM_GFS_v16_RRTMGP_ps_radiation_final_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_physics_cap, only: SCM_GFS_v16_RRTMGP_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_physics_cap, only: SCM_GFS_v16_RRTMGP_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_physics_cap, only: SCM_GFS_v16_RRTMGP_ps_physics_init_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_physics_cap, only: SCM_GFS_v16_RRTMGP_ps_physics_run_cap
   use ccpp_SCM_GFS_v16_RRTMGP_ps_physics_cap, only: SCM_GFS_v16_RRTMGP_ps_physics_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_cap, only: SCM_GFS_v16_no_nsst_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_cap, only: SCM_GFS_v16_no_nsst_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_cap, only: SCM_GFS_v16_no_nsst_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_cap, only: SCM_GFS_v16_no_nsst_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_cap, only: SCM_GFS_v16_no_nsst_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_time_vary_cap, only: SCM_GFS_v16_no_nsst_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_time_vary_cap, only: SCM_GFS_v16_no_nsst_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_time_vary_cap, only: SCM_GFS_v16_no_nsst_time_vary_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_time_vary_cap, only: SCM_GFS_v16_no_nsst_time_vary_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_time_vary_cap, only: SCM_GFS_v16_no_nsst_time_vary_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_radiation_cap, only: SCM_GFS_v16_no_nsst_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_radiation_cap, only: SCM_GFS_v16_no_nsst_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_radiation_cap, only: SCM_GFS_v16_no_nsst_radiation_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_radiation_cap, only: SCM_GFS_v16_no_nsst_radiation_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_radiation_cap, only: SCM_GFS_v16_no_nsst_radiation_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_physics_cap, only: SCM_GFS_v16_no_nsst_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_physics_cap, only: SCM_GFS_v16_no_nsst_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_physics_cap, only: SCM_GFS_v16_no_nsst_physics_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_physics_cap, only: SCM_GFS_v16_no_nsst_physics_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_physics_cap, only: SCM_GFS_v16_no_nsst_physics_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_cap, only: SCM_GFS_v16_no_nsst_ps_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_cap, only: SCM_GFS_v16_no_nsst_ps_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_cap, only: SCM_GFS_v16_no_nsst_ps_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_cap, only: SCM_GFS_v16_no_nsst_ps_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_cap, only: SCM_GFS_v16_no_nsst_ps_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_time_vary_cap, only: SCM_GFS_v16_no_nsst_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_time_vary_cap, only: SCM_GFS_v16_no_nsst_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_time_vary_cap, only: SCM_GFS_v16_no_nsst_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_time_vary_cap, only: SCM_GFS_v16_no_nsst_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_time_vary_cap, only: SCM_GFS_v16_no_nsst_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_radiation_cap, only: SCM_GFS_v16_no_nsst_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_radiation_cap, only: SCM_GFS_v16_no_nsst_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_radiation_cap, only: SCM_GFS_v16_no_nsst_ps_radiation_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_radiation_cap, only: SCM_GFS_v16_no_nsst_ps_radiation_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_radiation_cap, only: SCM_GFS_v16_no_nsst_ps_radiation_final_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_physics_cap, only: SCM_GFS_v16_no_nsst_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_physics_cap, only: SCM_GFS_v16_no_nsst_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_physics_cap, only: SCM_GFS_v16_no_nsst_ps_physics_init_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_physics_cap, only: SCM_GFS_v16_no_nsst_ps_physics_run_cap
   use ccpp_SCM_GFS_v16_no_nsst_ps_physics_cap, only: SCM_GFS_v16_no_nsst_ps_physics_final_cap
   use ccpp_SCM_GFS_v16_ps_cap, only: SCM_GFS_v16_ps_tsinit_cap
   use ccpp_SCM_GFS_v16_ps_cap, only: SCM_GFS_v16_ps_tsfinal_cap
   use ccpp_SCM_GFS_v16_ps_cap, only: SCM_GFS_v16_ps_init_cap
   use ccpp_SCM_GFS_v16_ps_cap, only: SCM_GFS_v16_ps_run_cap
   use ccpp_SCM_GFS_v16_ps_cap, only: SCM_GFS_v16_ps_final_cap
   use ccpp_SCM_GFS_v16_ps_time_vary_cap, only: SCM_GFS_v16_ps_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_ps_time_vary_cap, only: SCM_GFS_v16_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_ps_time_vary_cap, only: SCM_GFS_v16_ps_time_vary_init_cap
   use ccpp_SCM_GFS_v16_ps_time_vary_cap, only: SCM_GFS_v16_ps_time_vary_run_cap
   use ccpp_SCM_GFS_v16_ps_time_vary_cap, only: SCM_GFS_v16_ps_time_vary_final_cap
   use ccpp_SCM_GFS_v16_ps_radiation_cap, only: SCM_GFS_v16_ps_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_ps_radiation_cap, only: SCM_GFS_v16_ps_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_ps_radiation_cap, only: SCM_GFS_v16_ps_radiation_init_cap
   use ccpp_SCM_GFS_v16_ps_radiation_cap, only: SCM_GFS_v16_ps_radiation_run_cap
   use ccpp_SCM_GFS_v16_ps_radiation_cap, only: SCM_GFS_v16_ps_radiation_final_cap
   use ccpp_SCM_GFS_v16_ps_physics_cap, only: SCM_GFS_v16_ps_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_ps_physics_cap, only: SCM_GFS_v16_ps_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_ps_physics_cap, only: SCM_GFS_v16_ps_physics_init_cap
   use ccpp_SCM_GFS_v16_ps_physics_cap, only: SCM_GFS_v16_ps_physics_run_cap
   use ccpp_SCM_GFS_v16_ps_physics_cap, only: SCM_GFS_v16_ps_physics_final_cap
   use ccpp_SCM_GFS_v16_ugwpv1_cap, only: SCM_GFS_v16_ugwpv1_tsinit_cap
   use ccpp_SCM_GFS_v16_ugwpv1_cap, only: SCM_GFS_v16_ugwpv1_tsfinal_cap
   use ccpp_SCM_GFS_v16_ugwpv1_cap, only: SCM_GFS_v16_ugwpv1_init_cap
   use ccpp_SCM_GFS_v16_ugwpv1_cap, only: SCM_GFS_v16_ugwpv1_run_cap
   use ccpp_SCM_GFS_v16_ugwpv1_cap, only: SCM_GFS_v16_ugwpv1_final_cap
   use ccpp_SCM_GFS_v16_ugwpv1_time_vary_cap, only: SCM_GFS_v16_ugwpv1_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v16_ugwpv1_time_vary_cap, only: SCM_GFS_v16_ugwpv1_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v16_ugwpv1_time_vary_cap, only: SCM_GFS_v16_ugwpv1_time_vary_init_cap
   use ccpp_SCM_GFS_v16_ugwpv1_time_vary_cap, only: SCM_GFS_v16_ugwpv1_time_vary_run_cap
   use ccpp_SCM_GFS_v16_ugwpv1_time_vary_cap, only: SCM_GFS_v16_ugwpv1_time_vary_final_cap
   use ccpp_SCM_GFS_v16_ugwpv1_radiation_cap, only: SCM_GFS_v16_ugwpv1_radiation_tsinit_cap
   use ccpp_SCM_GFS_v16_ugwpv1_radiation_cap, only: SCM_GFS_v16_ugwpv1_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v16_ugwpv1_radiation_cap, only: SCM_GFS_v16_ugwpv1_radiation_init_cap
   use ccpp_SCM_GFS_v16_ugwpv1_radiation_cap, only: SCM_GFS_v16_ugwpv1_radiation_run_cap
   use ccpp_SCM_GFS_v16_ugwpv1_radiation_cap, only: SCM_GFS_v16_ugwpv1_radiation_final_cap
   use ccpp_SCM_GFS_v16_ugwpv1_physics_cap, only: SCM_GFS_v16_ugwpv1_physics_tsinit_cap
   use ccpp_SCM_GFS_v16_ugwpv1_physics_cap, only: SCM_GFS_v16_ugwpv1_physics_tsfinal_cap
   use ccpp_SCM_GFS_v16_ugwpv1_physics_cap, only: SCM_GFS_v16_ugwpv1_physics_init_cap
   use ccpp_SCM_GFS_v16_ugwpv1_physics_cap, only: SCM_GFS_v16_ugwpv1_physics_run_cap
   use ccpp_SCM_GFS_v16_ugwpv1_physics_cap, only: SCM_GFS_v16_ugwpv1_physics_final_cap
   use ccpp_SCM_GSD_v1_cap, only: SCM_GSD_v1_tsinit_cap
   use ccpp_SCM_GSD_v1_cap, only: SCM_GSD_v1_tsfinal_cap
   use ccpp_SCM_GSD_v1_cap, only: SCM_GSD_v1_init_cap
   use ccpp_SCM_GSD_v1_cap, only: SCM_GSD_v1_run_cap
   use ccpp_SCM_GSD_v1_cap, only: SCM_GSD_v1_final_cap
   use ccpp_SCM_GSD_v1_time_vary_cap, only: SCM_GSD_v1_time_vary_tsinit_cap
   use ccpp_SCM_GSD_v1_time_vary_cap, only: SCM_GSD_v1_time_vary_tsfinal_cap
   use ccpp_SCM_GSD_v1_time_vary_cap, only: SCM_GSD_v1_time_vary_init_cap
   use ccpp_SCM_GSD_v1_time_vary_cap, only: SCM_GSD_v1_time_vary_run_cap
   use ccpp_SCM_GSD_v1_time_vary_cap, only: SCM_GSD_v1_time_vary_final_cap
   use ccpp_SCM_GSD_v1_radiation_cap, only: SCM_GSD_v1_radiation_tsinit_cap
   use ccpp_SCM_GSD_v1_radiation_cap, only: SCM_GSD_v1_radiation_tsfinal_cap
   use ccpp_SCM_GSD_v1_radiation_cap, only: SCM_GSD_v1_radiation_init_cap
   use ccpp_SCM_GSD_v1_radiation_cap, only: SCM_GSD_v1_radiation_run_cap
   use ccpp_SCM_GSD_v1_radiation_cap, only: SCM_GSD_v1_radiation_final_cap
   use ccpp_SCM_GSD_v1_physics_cap, only: SCM_GSD_v1_physics_tsinit_cap
   use ccpp_SCM_GSD_v1_physics_cap, only: SCM_GSD_v1_physics_tsfinal_cap
   use ccpp_SCM_GSD_v1_physics_cap, only: SCM_GSD_v1_physics_init_cap
   use ccpp_SCM_GSD_v1_physics_cap, only: SCM_GSD_v1_physics_run_cap
   use ccpp_SCM_GSD_v1_physics_cap, only: SCM_GSD_v1_physics_final_cap
   use ccpp_SCM_GSD_v1_ps_cap, only: SCM_GSD_v1_ps_tsinit_cap
   use ccpp_SCM_GSD_v1_ps_cap, only: SCM_GSD_v1_ps_tsfinal_cap
   use ccpp_SCM_GSD_v1_ps_cap, only: SCM_GSD_v1_ps_init_cap
   use ccpp_SCM_GSD_v1_ps_cap, only: SCM_GSD_v1_ps_run_cap
   use ccpp_SCM_GSD_v1_ps_cap, only: SCM_GSD_v1_ps_final_cap
   use ccpp_SCM_GSD_v1_ps_time_vary_cap, only: SCM_GSD_v1_ps_time_vary_tsinit_cap
   use ccpp_SCM_GSD_v1_ps_time_vary_cap, only: SCM_GSD_v1_ps_time_vary_tsfinal_cap
   use ccpp_SCM_GSD_v1_ps_time_vary_cap, only: SCM_GSD_v1_ps_time_vary_init_cap
   use ccpp_SCM_GSD_v1_ps_time_vary_cap, only: SCM_GSD_v1_ps_time_vary_run_cap
   use ccpp_SCM_GSD_v1_ps_time_vary_cap, only: SCM_GSD_v1_ps_time_vary_final_cap
   use ccpp_SCM_GSD_v1_ps_radiation_cap, only: SCM_GSD_v1_ps_radiation_tsinit_cap
   use ccpp_SCM_GSD_v1_ps_radiation_cap, only: SCM_GSD_v1_ps_radiation_tsfinal_cap
   use ccpp_SCM_GSD_v1_ps_radiation_cap, only: SCM_GSD_v1_ps_radiation_init_cap
   use ccpp_SCM_GSD_v1_ps_radiation_cap, only: SCM_GSD_v1_ps_radiation_run_cap
   use ccpp_SCM_GSD_v1_ps_radiation_cap, only: SCM_GSD_v1_ps_radiation_final_cap
   use ccpp_SCM_GSD_v1_ps_physics_cap, only: SCM_GSD_v1_ps_physics_tsinit_cap
   use ccpp_SCM_GSD_v1_ps_physics_cap, only: SCM_GSD_v1_ps_physics_tsfinal_cap
   use ccpp_SCM_GSD_v1_ps_physics_cap, only: SCM_GSD_v1_ps_physics_init_cap
   use ccpp_SCM_GSD_v1_ps_physics_cap, only: SCM_GSD_v1_ps_physics_run_cap
   use ccpp_SCM_GSD_v1_ps_physics_cap, only: SCM_GSD_v1_ps_physics_final_cap
   use ccpp_SCM_RRFS_v1beta_cap, only: SCM_RRFS_v1beta_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_cap, only: SCM_RRFS_v1beta_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_cap, only: SCM_RRFS_v1beta_init_cap
   use ccpp_SCM_RRFS_v1beta_cap, only: SCM_RRFS_v1beta_run_cap
   use ccpp_SCM_RRFS_v1beta_cap, only: SCM_RRFS_v1beta_final_cap
   use ccpp_SCM_RRFS_v1beta_time_vary_cap, only: SCM_RRFS_v1beta_time_vary_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_time_vary_cap, only: SCM_RRFS_v1beta_time_vary_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_time_vary_cap, only: SCM_RRFS_v1beta_time_vary_init_cap
   use ccpp_SCM_RRFS_v1beta_time_vary_cap, only: SCM_RRFS_v1beta_time_vary_run_cap
   use ccpp_SCM_RRFS_v1beta_time_vary_cap, only: SCM_RRFS_v1beta_time_vary_final_cap
   use ccpp_SCM_RRFS_v1beta_radiation_cap, only: SCM_RRFS_v1beta_radiation_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_radiation_cap, only: SCM_RRFS_v1beta_radiation_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_radiation_cap, only: SCM_RRFS_v1beta_radiation_init_cap
   use ccpp_SCM_RRFS_v1beta_radiation_cap, only: SCM_RRFS_v1beta_radiation_run_cap
   use ccpp_SCM_RRFS_v1beta_radiation_cap, only: SCM_RRFS_v1beta_radiation_final_cap
   use ccpp_SCM_RRFS_v1beta_physics_cap, only: SCM_RRFS_v1beta_physics_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_physics_cap, only: SCM_RRFS_v1beta_physics_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_physics_cap, only: SCM_RRFS_v1beta_physics_init_cap
   use ccpp_SCM_RRFS_v1beta_physics_cap, only: SCM_RRFS_v1beta_physics_run_cap
   use ccpp_SCM_RRFS_v1beta_physics_cap, only: SCM_RRFS_v1beta_physics_final_cap
   use ccpp_SCM_RRFS_v1beta_ps_cap, only: SCM_RRFS_v1beta_ps_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_ps_cap, only: SCM_RRFS_v1beta_ps_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_ps_cap, only: SCM_RRFS_v1beta_ps_init_cap
   use ccpp_SCM_RRFS_v1beta_ps_cap, only: SCM_RRFS_v1beta_ps_run_cap
   use ccpp_SCM_RRFS_v1beta_ps_cap, only: SCM_RRFS_v1beta_ps_final_cap
   use ccpp_SCM_RRFS_v1beta_ps_time_vary_cap, only: SCM_RRFS_v1beta_ps_time_vary_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_ps_time_vary_cap, only: SCM_RRFS_v1beta_ps_time_vary_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_ps_time_vary_cap, only: SCM_RRFS_v1beta_ps_time_vary_init_cap
   use ccpp_SCM_RRFS_v1beta_ps_time_vary_cap, only: SCM_RRFS_v1beta_ps_time_vary_run_cap
   use ccpp_SCM_RRFS_v1beta_ps_time_vary_cap, only: SCM_RRFS_v1beta_ps_time_vary_final_cap
   use ccpp_SCM_RRFS_v1beta_ps_radiation_cap, only: SCM_RRFS_v1beta_ps_radiation_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_ps_radiation_cap, only: SCM_RRFS_v1beta_ps_radiation_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_ps_radiation_cap, only: SCM_RRFS_v1beta_ps_radiation_init_cap
   use ccpp_SCM_RRFS_v1beta_ps_radiation_cap, only: SCM_RRFS_v1beta_ps_radiation_run_cap
   use ccpp_SCM_RRFS_v1beta_ps_radiation_cap, only: SCM_RRFS_v1beta_ps_radiation_final_cap
   use ccpp_SCM_RRFS_v1beta_ps_physics_cap, only: SCM_RRFS_v1beta_ps_physics_tsinit_cap
   use ccpp_SCM_RRFS_v1beta_ps_physics_cap, only: SCM_RRFS_v1beta_ps_physics_tsfinal_cap
   use ccpp_SCM_RRFS_v1beta_ps_physics_cap, only: SCM_RRFS_v1beta_ps_physics_init_cap
   use ccpp_SCM_RRFS_v1beta_ps_physics_cap, only: SCM_RRFS_v1beta_ps_physics_run_cap
   use ccpp_SCM_RRFS_v1beta_ps_physics_cap, only: SCM_RRFS_v1beta_ps_physics_final_cap
   use ccpp_SCM_csawmg_cap, only: SCM_csawmg_tsinit_cap
   use ccpp_SCM_csawmg_cap, only: SCM_csawmg_tsfinal_cap
   use ccpp_SCM_csawmg_cap, only: SCM_csawmg_init_cap
   use ccpp_SCM_csawmg_cap, only: SCM_csawmg_run_cap
   use ccpp_SCM_csawmg_cap, only: SCM_csawmg_final_cap
   use ccpp_SCM_csawmg_time_vary_cap, only: SCM_csawmg_time_vary_tsinit_cap
   use ccpp_SCM_csawmg_time_vary_cap, only: SCM_csawmg_time_vary_tsfinal_cap
   use ccpp_SCM_csawmg_time_vary_cap, only: SCM_csawmg_time_vary_init_cap
   use ccpp_SCM_csawmg_time_vary_cap, only: SCM_csawmg_time_vary_run_cap
   use ccpp_SCM_csawmg_time_vary_cap, only: SCM_csawmg_time_vary_final_cap
   use ccpp_SCM_csawmg_radiation_cap, only: SCM_csawmg_radiation_tsinit_cap
   use ccpp_SCM_csawmg_radiation_cap, only: SCM_csawmg_radiation_tsfinal_cap
   use ccpp_SCM_csawmg_radiation_cap, only: SCM_csawmg_radiation_init_cap
   use ccpp_SCM_csawmg_radiation_cap, only: SCM_csawmg_radiation_run_cap
   use ccpp_SCM_csawmg_radiation_cap, only: SCM_csawmg_radiation_final_cap
   use ccpp_SCM_csawmg_physics_cap, only: SCM_csawmg_physics_tsinit_cap
   use ccpp_SCM_csawmg_physics_cap, only: SCM_csawmg_physics_tsfinal_cap
   use ccpp_SCM_csawmg_physics_cap, only: SCM_csawmg_physics_init_cap
   use ccpp_SCM_csawmg_physics_cap, only: SCM_csawmg_physics_run_cap
   use ccpp_SCM_csawmg_physics_cap, only: SCM_csawmg_physics_final_cap
   use ccpp_SCM_csawmg_ps_cap, only: SCM_csawmg_ps_tsinit_cap
   use ccpp_SCM_csawmg_ps_cap, only: SCM_csawmg_ps_tsfinal_cap
   use ccpp_SCM_csawmg_ps_cap, only: SCM_csawmg_ps_init_cap
   use ccpp_SCM_csawmg_ps_cap, only: SCM_csawmg_ps_run_cap
   use ccpp_SCM_csawmg_ps_cap, only: SCM_csawmg_ps_final_cap
   use ccpp_SCM_csawmg_ps_time_vary_cap, only: SCM_csawmg_ps_time_vary_tsinit_cap
   use ccpp_SCM_csawmg_ps_time_vary_cap, only: SCM_csawmg_ps_time_vary_tsfinal_cap
   use ccpp_SCM_csawmg_ps_time_vary_cap, only: SCM_csawmg_ps_time_vary_init_cap
   use ccpp_SCM_csawmg_ps_time_vary_cap, only: SCM_csawmg_ps_time_vary_run_cap
   use ccpp_SCM_csawmg_ps_time_vary_cap, only: SCM_csawmg_ps_time_vary_final_cap
   use ccpp_SCM_csawmg_ps_radiation_cap, only: SCM_csawmg_ps_radiation_tsinit_cap
   use ccpp_SCM_csawmg_ps_radiation_cap, only: SCM_csawmg_ps_radiation_tsfinal_cap
   use ccpp_SCM_csawmg_ps_radiation_cap, only: SCM_csawmg_ps_radiation_init_cap
   use ccpp_SCM_csawmg_ps_radiation_cap, only: SCM_csawmg_ps_radiation_run_cap
   use ccpp_SCM_csawmg_ps_radiation_cap, only: SCM_csawmg_ps_radiation_final_cap
   use ccpp_SCM_csawmg_ps_physics_cap, only: SCM_csawmg_ps_physics_tsinit_cap
   use ccpp_SCM_csawmg_ps_physics_cap, only: SCM_csawmg_ps_physics_tsfinal_cap
   use ccpp_SCM_csawmg_ps_physics_cap, only: SCM_csawmg_ps_physics_init_cap
   use ccpp_SCM_csawmg_ps_physics_cap, only: SCM_csawmg_ps_physics_run_cap
   use ccpp_SCM_csawmg_ps_physics_cap, only: SCM_csawmg_ps_physics_final_cap
   use gmtb_scm_type_defs, only: physics
   use ozne_def, only: levozp
   use ozne_def, only: oz_coeff
   use h2o_def, only: levh2o
   use h2o_def, only: h2o_coeff
   use gmtb_scm_physical_constants, only: con_t0c
   use gmtb_scm_physical_constants, only: con_p0
   use gmtb_scm_physical_constants, only: con_g
   use gmtb_scm_physical_constants, only: con_rd
   use gmtb_scm_physical_constants, only: con_eps
   use gmtb_scm_physical_constants, only: con_pi
   use gmtb_scm_physical_constants, only: con_rerth
   use gmtb_scm_physical_constants, only: con_omega
   use gmtb_scm_physical_constants, only: con_cp
   use gmtb_scm_physical_constants, only: con_rv
   use gmtb_scm_physical_constants, only: con_fvirt
   use gmtb_scm_physical_constants, only: con_ttp
   use gmtb_scm_physical_constants, only: con_hvap
   use gmtb_scm_physical_constants, only: con_hfus
   use gmtb_scm_physical_constants, only: con_epsm1
   use gmtb_scm_physical_constants, only: con_rog
   use gmtb_scm_physical_constants, only: con_rocp
   use GFS_typedefs, only: LTP
   use gmtb_scm_physical_constants, only: con_tice
   use GFS_typedefs, only: huge
   use gmtb_scm_physical_constants, only: con_jcal
   use gmtb_scm_physical_constants, only: con_rhw0
   use gmtb_scm_physical_constants, only: con_sbc
   use gmtb_scm_physical_constants, only: rlapse
   use gmtb_scm_physical_constants, only: rhowater
   use gmtb_scm_physical_constants, only: con_cliq
   use gmtb_scm_physical_constants, only: con_csol
   use gmtb_scm_physical_constants, only: con_epsq
   use ozne_def, only: oz_pres
   use h2o_def, only: h2o_pres
   use gmtb_scm_physical_constants, only: con_cvap
   use gmtb_scm_physical_constants, only: con_vonKarman
   use gmtb_scm_physical_constants, only: con_epsqs

   implicit none

   private
   public :: ccpp_physics_timestep_init,ccpp_physics_timestep_finalize,ccpp_physics_init,ccpp_physics_run,ccpp_physics_finalize

   contains

   subroutine ccpp_physics_timestep_init(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="HAFS_v0_hwrf") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ACM_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ACM_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ACM_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ACM_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ACM_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_FA") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_FA_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_FA_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_FA_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_FA_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_MYJ") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_MYJ_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_MYJ_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_MYJ_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_MYJ_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_YSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_YSU_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_YSU_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_YSU_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_YSU_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_noahmp") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_noahmp_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_noahmp_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_noahmp_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_noahmp_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_saYSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_saYSU_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_saYSU_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_saYSU_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_saYSU_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ugwpv1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ugwpv1_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ugwpv1_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ugwpv1_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ugwpv1_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_csawmg") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_csawmg_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_ps_time_vary_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_ps_radiation_tsinit_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_ps_physics_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_ps_tsinit_cap(physics=physics,cdata=cdata,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_timestep_init

   subroutine ccpp_physics_timestep_finalize(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="HAFS_v0_hwrf") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ACM_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ACM_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ACM_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ACM_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ACM_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_FA") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_FA_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_FA_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_FA_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_FA_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_MYJ") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_MYJ_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_MYJ_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_MYJ_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_MYJ_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_YSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_YSU_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_YSU_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_YSU_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_YSU_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_noahmp") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_noahmp_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_noahmp_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_noahmp_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_noahmp_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_saYSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_saYSU_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_saYSU_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_saYSU_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_saYSU_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ugwpv1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ugwpv1_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ugwpv1_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ugwpv1_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ugwpv1_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_ps_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_csawmg") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_tsfinal_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_csawmg_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_ps_time_vary_tsfinal_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_ps_radiation_tsfinal_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_ps_physics_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_ps_tsfinal_cap(cdata=cdata)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_timestep_finalize

   subroutine ccpp_physics_init(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="HAFS_v0_hwrf") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ACM_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ACM_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ACM_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ACM_ps_physics_init_cap(physics=physics,con_p0=con_p0,cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ACM_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_FA") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_FA_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_FA_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_FA_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_FA_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_MYJ") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_MYJ_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_MYJ_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_MYJ_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_MYJ_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_radiation_init_cap(physics=physics,cdata=cdata)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_radiation_init_cap(physics=physics,cdata=cdata)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_YSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_YSU_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_YSU_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_YSU_ps_physics_init_cap(physics=physics,con_p0=con_p0,cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_YSU_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_noahmp") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_noahmp_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_noahmp_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_noahmp_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_noahmp_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_saYSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_saYSU_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_saYSU_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_saYSU_ps_physics_init_cap(physics=physics,con_p0=con_p0,cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_saYSU_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_radiation_init_cap(physics=physics,cdata=cdata)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_ps_radiation_init_cap(physics=physics,cdata=cdata)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ugwpv1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ugwpv1_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ugwpv1_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ugwpv1_physics_init_cap(physics=physics,cdata=cdata,con_pi=con_pi,con_rerth=con_rerth,con_p0=con_p0, &
                  con_g=con_g,con_omega=con_omega,con_cp=con_cp,con_rd=con_rd,con_rv=con_rv, &
                  con_fvirt=con_fvirt)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ugwpv1_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_pi=con_pi,con_rerth=con_rerth,con_p0=con_p0, &
                  con_g=con_g,con_omega=con_omega,con_cp=con_cp,con_rd=con_rd,con_rv=con_rv, &
                  con_fvirt=con_fvirt)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_ps_physics_init_cap(physics=physics,con_p0=con_p0,cdata=cdata,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_physics_init_cap(cdata=cdata,physics=physics,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_ps_physics_init_cap(physics=physics,con_p0=con_p0,cdata=cdata,con_g=con_g,con_rd=con_rd,con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_eps=con_eps)

         end if

      else if (trim(suite_name)=="SCM_csawmg") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_rv=con_rv, &
                  con_cp=con_cp,con_eps=con_eps,con_ttp=con_ttp,con_hvap=con_hvap,con_hfus=con_hfus)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_rv=con_rv,con_cp=con_cp,con_eps=con_eps,con_ttp=con_ttp,con_hvap=con_hvap, &
                  con_hfus=con_hfus)

         end if

      else if (trim(suite_name)=="SCM_csawmg_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_ps_time_vary_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_ps_radiation_init_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_ps_physics_init_cap(physics=physics,cdata=cdata,con_p0=con_p0,con_g=con_g,con_rd=con_rd,con_rv=con_rv, &
                  con_cp=con_cp,con_eps=con_eps,con_ttp=con_ttp,con_hvap=con_hvap,con_hfus=con_hfus)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_ps_init_cap(cdata=cdata,physics=physics,levozp=levozp,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_coeff=h2o_coeff,con_t0c=con_t0c,con_p0=con_p0,con_g=con_g,con_rd=con_rd, &
                  con_rv=con_rv,con_cp=con_cp,con_eps=con_eps,con_ttp=con_ttp,con_hvap=con_hvap, &
                  con_hfus=con_hfus)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_init

   subroutine ccpp_physics_run(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="HAFS_v0_hwrf") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_eps=con_eps, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_epsm1=con_epsm1,con_rhw0=con_rhw0, &
                  con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater,con_cliq=con_cliq, &
                  con_csol=con_csol,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater, &
                  con_cliq=con_cliq,con_csol=con_csol,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_eps=con_eps, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_epsm1=con_epsm1,con_rhw0=con_rhw0, &
                  con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater,con_cliq=con_cliq, &
                  con_csol=con_csol,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater, &
                  con_cliq=con_cliq,con_csol=con_csol,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ACM_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ACM_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ACM_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ACM_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ACM_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_FA") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_FA_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_FA_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_FA_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_FA_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_MYJ") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_MYJ_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_MYJ_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_MYJ_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_rd=con_rd,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt,con_rhw0=con_rhw0, &
                  con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_MYJ_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_YSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_YSU_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_YSU_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_YSU_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_rv=con_rv,con_eps=con_eps,con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_YSU_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_pi=con_pi, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_rd=con_rd,con_sbc=con_sbc,con_t0c=con_t0c,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_sbc=con_sbc,con_t0c=con_t0c, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_noahmp") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_noahmp_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_noahmp_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_noahmp_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,rhowater=rhowater,con_t0c=con_t0c,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_noahmp_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater, &
                  con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_saYSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_saYSU_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_saYSU_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_saYSU_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_rv=con_rv,con_eps=con_eps,con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_saYSU_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_pi=con_pi, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_rv=con_rv,con_hfus=con_hfus,con_eps=con_eps,con_epsm1=con_epsm1,con_pi=con_pi, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_epsqs=con_epsqs,con_g=con_g,con_rd=con_rd,con_pi=con_pi,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_hfus=con_hfus, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_rd=con_rd,con_sbc=con_sbc,con_t0c=con_t0c,con_rv=con_rv,con_hfus=con_hfus, &
                  con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_sbc=con_sbc,con_t0c=con_t0c, &
                  con_rv=con_rv,con_hfus=con_hfus,con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_rv=con_rv,con_hfus=con_hfus,con_eps=con_eps,con_epsm1=con_epsm1,con_pi=con_pi, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_hfus=con_hfus, &
                  con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap,con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_rv=con_rv,con_hfus=con_hfus,con_eps=con_eps,con_epsm1=con_epsm1,con_pi=con_pi, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_rv=con_rv,con_hfus=con_hfus, &
                  con_pi=con_pi,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap,con_t0c=con_t0c)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ugwpv1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ugwpv1_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ugwpv1_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ugwpv1_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ugwpv1_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq, &
                  levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres, &
                  h2o_coeff=h2o_coeff)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal,con_eps=con_eps, &
                  con_epsm1=con_epsm1,con_fvirt=con_fvirt,con_rd=con_rd,con_rhw0=con_rhw0, &
                  con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,rhowater=rhowater, &
                  con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_eps=con_eps)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff)

         end if

      else if (trim(suite_name)=="SCM_csawmg") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_tice=con_tice,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_fvirt=con_fvirt,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0,con_sbc=con_sbc, &
                  con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c,con_rv=con_rv,con_omega=con_omega, &
                  con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o, &
                  h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_hvap=con_hvap,huge=huge,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,con_sbc=con_sbc,con_pi=con_pi,rlapse=rlapse,con_t0c=con_t0c, &
                  con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres, &
                  oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq, &
                  con_cvap=con_cvap)

         end if

      else if (trim(suite_name)=="SCM_csawmg_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_ps_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_ps_radiation_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP)
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_ps_physics_run_cap(physics=physics,cdata=cdata,con_g=con_g,con_cp=con_cp,con_hvap=con_hvap, &
                  huge=huge,con_rd=con_rd,con_fvirt=con_fvirt,con_vonKarman=con_vonKarman, &
                  con_pi=con_pi,con_rv=con_rv,con_omega=con_omega,con_epsq=con_epsq,levozp=levozp, &
                  oz_pres=oz_pres,oz_coeff=oz_coeff,levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff, &
                  con_cliq=con_cliq,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_t0c=con_t0c)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_ps_run_cap(physics=physics,cdata=cdata,con_eps=con_eps,con_epsm1=con_epsm1,con_fvirt=con_fvirt, &
                  con_rog=con_rog,con_rocp=con_rocp,con_rd=con_rd,LTP=LTP,con_g=con_g,con_cp=con_cp, &
                  con_hvap=con_hvap,huge=huge,con_vonKarman=con_vonKarman,con_pi=con_pi,con_rv=con_rv, &
                  con_omega=con_omega,con_epsq=con_epsq,levozp=levozp,oz_pres=oz_pres,oz_coeff=oz_coeff, &
                  levh2o=levh2o,h2o_pres=h2o_pres,h2o_coeff=h2o_coeff,con_cliq=con_cliq,con_cvap=con_cvap, &
                  con_t0c=con_t0c)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_run

   subroutine ccpp_physics_finalize(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="HAFS_v0_hwrf") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="HAFS_v0_hwrf_thompson_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = HAFS_v0_hwrf_thompson_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = HAFS_v0_hwrf_thompson_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = HAFS_v0_hwrf_thompson_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = HAFS_v0_hwrf_thompson_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ACM_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ACM_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ACM_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ACM_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ACM_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_FA") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_FA_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_FA_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_FA_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_FA_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_MYJ") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_MYJ_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_MYJ_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_MYJ_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_MYJ_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_RRTMGP_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_RRTMGP_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_YSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_YSU_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_YSU_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_YSU_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_YSU_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_no_nsst_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_no_nsst_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_noahmp") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_noahmp_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_noahmp_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_noahmp_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_noahmp_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v15p2_saYSU_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v15p2_saYSU_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v15p2_saYSU_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v15p2_saYSU_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v15p2_saYSU_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_RRTMGP_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_RRTMGP_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_RRTMGP_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_RRTMGP_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_RRTMGP_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_no_nsst_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_no_nsst_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_no_nsst_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_no_nsst_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_no_nsst_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GFS_v16_ugwpv1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v16_ugwpv1_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v16_ugwpv1_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GFS_v16_ugwpv1_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v16_ugwpv1_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_GSD_v1_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GSD_v1_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GSD_v1_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_GSD_v1_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GSD_v1_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_RRFS_v1beta_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_RRFS_v1beta_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_RRFS_v1beta_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_RRFS_v1beta_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_RRFS_v1beta_ps_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_csawmg") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_final_cap(cdata=cdata)

         end if

      else if (trim(suite_name)=="SCM_csawmg_ps") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_csawmg_ps_time_vary_final_cap(cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_csawmg_ps_radiation_final_cap()
            else if (trim(group_name)=="physics") then
               ierr = SCM_csawmg_ps_physics_final_cap(cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_csawmg_ps_final_cap(cdata=cdata)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_finalize

end module ccpp_static_api
