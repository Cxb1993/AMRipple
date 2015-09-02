module params_module
  use bl_types

  real(dp_t) :: dens_ratio
  integer    :: pres_id,f_id,uvel_id,vvel_id,wvel_id,rho_id
  real(dp_t) :: max_deltaT,CFL
  real(dp_t) :: anum
  real(dp_t), parameter :: dp_tiny = tiny(anum)

end module params_module
