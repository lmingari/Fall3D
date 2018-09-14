subroutine write_DBS_data(fname,nx,ny,nz,itime, &
           topg,LDU,SOIL,vegfra,u,v,w,T,Tp,p,qv, &
           ro,pblh,ust,rmol,znt,spd10,L,smoi,prat)
  !********************************************************************
  !*
  !*   Writes DBS data in NetCDF format. This routine MUST be called
  !*   after write_WRF_grid and for each time step
  !*
  !********************************************************************
  use KindType, ONLY :  ip,rp
  use netcdf
  use DBS_nc
  implicit none
  !
  character(len=*) :: fname
  integer(ip)      :: itime,nx,ny,nz
  real(rp)         :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),T(nx,ny,nz)
  real(rp)         :: Tp(nx,ny,nz), p(nx,ny,nz),qv(nx,ny,nz),ro(nx,ny,nz)
  real(rp)         :: topg(nx,ny),LDU(nx,ny),SOIL(nx,ny),vegfra(nx,ny),pblh(nx,ny)
  real(rp)         :: ust(nx,ny),rmol(nx,ny),znt(nx,ny),spd10(nx,ny)
  real(rp)         :: L(nx,ny),smoi(nx,ny),prat(nx,ny)
  !
  !*** Open netCDF file and get ncID
  !
  if( nf90_open(TRIM(fname),NF90_WRITE, ncID) /= 0 ) call runend('write_dbs_data : Error in nf90_open')
  !
  !*** Inquire var ID and write variables
  !
  !
  !***   2D variables
  !
  if( nf90_inq_varid(ncID,topg_DBS_name,topg_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting topg_DBS_ID')
  if( nf90_put_var(ncID, topg_DBS_ID, topg, start=(/1,1/), count=(/nx,ny/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable topg')
  call writim('topog')
  !
  if( nf90_inq_varid(ncID,LDU_DBS_name,LDU_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting ldu_DBS_ID')
  if( nf90_put_var(ncID, LDU_DBS_ID, LDU, start=(/1,1/), count=(/nx,ny/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable LDU')
  call writim('LDU')
  !
  if( nf90_inq_varid(ncID,SOIL_DBS_name,SOIL_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting soil_DBS_ID')
  if( nf90_put_var(ncID, SOIL_DBS_ID, SOIL, start=(/1,1/), count=(/nx,ny/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable SOIL')
  call writim('SOIL')
  !
  if( nf90_inq_varid(ncID,vegfra_DBS_name,vegfra_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting vegfra_DBS_ID')
  if( nf90_put_var(ncID, vegfra_DBS_ID, vegfra, start=(/1,1/), count=(/nx,ny/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable VEGFRA')
  call writim('VEGFRA')
  !
  if( nf90_inq_varid(ncID,pblh_DBS_name,pblh_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting pblh_DBS_ID')
  if( nf90_put_var(ncID, pblh_DBS_ID, pblh, start=(/1,1,itime/), count=(/nx,ny,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable pblh')
  call writim('pblh')
  !
  if( nf90_inq_varid(ncID,ust_DBS_name,ust_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting ust_DBS_ID')
  if( nf90_put_var(ncID, ust_DBS_ID, ust, start=(/1,1,itime/), count=(/nx,ny,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable ust')
  call writim('ust')
  !
  if( nf90_inq_varid(ncID,smoi_DBS_name,smoi_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting smoi_DBS_ID')
  if( nf90_put_var(ncID, smoi_DBS_ID, smoi, start=(/1,1,itime/), count=(/nx,ny,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable smoi')
  call writim('smoi')
  !
  if( nf90_inq_varid(ncID,rmol_DBS_name,rmol_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting rmol_DBS_ID')
  if( nf90_put_var(ncID, rmol_DBS_ID, rmol, start=(/1,1,itime/), count=(/nx,ny,1/)) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable rmol')
  call writim('rmol')
  !
  if( nf90_inq_varid(ncID,znt_DBS_name,znt_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting znt_DBS_ID')
  if( nf90_put_var(ncID, znt_DBS_ID, znt, start=(/1,1,itime/), count=(/nx,ny,1/)) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable znt')
  call writim('znt')
  !
  if( nf90_inq_varid(ncID,spd10_DBS_name,spd10_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting spd10_DBS_ID')
  if( nf90_put_var(ncID, spd10_DBS_ID, spd10, start=(/1,1,itime/), count=(/nx,ny,1/)) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable spd10')
  call writim('spd10')
  !
  if( nf90_inq_varid(ncID,L_DBS_name,L_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting L_DBS_ID')
  if( nf90_put_var(ncID, L_DBS_ID, L, start=(/1,1,itime/), count=(/nx,ny,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable L')
  call writim('L')
  !
  if( nf90_inq_varid(ncID,prat_DBS_name,prat_DBS_ID) /= 0) &
       call runend('write_DBS_data : Error getting prat_DBS_ID')
  if( nf90_put_var(ncID, prat_DBS_ID, prat, start=(/1,1,itime/), count=(/nx,ny,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable prat')
  call writim('prat')
  !
  !***   3D variables
  !
  if( nf90_inq_varid(ncID,u_DBS_name,u_DBS_ID) /= 0) call runend('write_DBS_data : Error getting u_DBS_ID')
  if( nf90_put_var(ncID, u_DBS_ID, u, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable u')
  call writim('u')
  !
  if( nf90_inq_varid(ncID,v_DBS_name,v_DBS_ID) /= 0) call runend('write_DBS_data : Error getting v_DBS_ID')
  if( nf90_put_var(ncID, v_DBS_ID, v, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable v')
  call writim('v')
  !
  if( nf90_inq_varid(ncID,w_DBS_name,w_DBS_ID) /= 0) call runend('write_DBS_data : Error getting w_DBS_ID')
  if( nf90_put_var(ncID, w_DBS_ID, w, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable w')
  call writim('w')
  !
  if( nf90_inq_varid(ncID,T_DBS_name,T_DBS_ID) /= 0) call runend('write_DBS_data : Error getting T_DBS_ID')
  if( nf90_put_var(ncID, T_DBS_ID, T, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable T')
  call writim('T')
  !
  if( nf90_inq_varid(ncID,Tp_DBS_name,Tp_DBS_ID) /= 0) call runend('write_DBS_data : Error getting Tp_DBS_ID')
  if( nf90_put_var(ncID, Tp_DBS_ID, Tp, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable Tp')
  call writim('Tp')
  !
  if( nf90_inq_varid(ncID,p_DBS_name,p_DBS_ID) /= 0) call runend('write_DBS_data : Error getting p_DBS_ID')
  if( nf90_put_var(ncID, p_DBS_ID, p, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable P')
  call writim('p')
  !
  if( nf90_inq_varid(ncID,qv_DBS_name,qv_DBS_ID) /= 0) call runend('write_DBS_data : Error getting qv_DBS_ID')
  if( nf90_put_var(ncID, qv_DBS_ID, qv, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable qv')
  call writim('qv')
  !
  if( nf90_inq_varid(ncID,ro_DBS_name,ro_DBS_ID) /= 0) call runend('write_DBS_data : Error getting ro_DBS_ID')
  if( nf90_put_var(ncID, ro_DBS_ID, ro, start=(/1,1,1,itime/), count=(/nx,ny,nz,1/) ) /= 0 ) &
       call wriwar('write_dbs_data: error in nf90_put_var for variable qv')
  call writim('ro')
  !
  !*** Close the file
  !
  if( nf90_close(ncID) /= 0) call runend('write_dbs_data : Error in nf90_close')
  !
  return
end subroutine write_dbs_data
