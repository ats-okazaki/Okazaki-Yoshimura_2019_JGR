    program main
    implicit none

! param
    integer,parameter   :: ystr=1871
    integer,parameter   :: yend=2007
    integer,parameter   :: ijmax=128*64
    integer,parameter   :: yymax=yend-ystr+1
! input
    real*4  :: prcp(ijmax,12,yymax)
    real*4  :: prcp01(ijmax,12,yymax)
    real*4  :: evap(ijmax,12,yymax)
    real*4  :: evap01(ijmax,12,yymax)
    real*4  :: sst(ijmax,12,yymax)
    real*4  :: dow(ijmax),dowd(ijmax)
! output
    real*4  :: d_coral(ijmax,12,yymax) !! annual mean of isotope ration in speleothem (VPDB)
    real*4  :: d_coral2d(ijmax,yymax) !! annual mean of isotope ration in speleothem (VPDB)
! internal work
    real*4  :: pterm,eterm,dterm
    real*4,parameter    :: q=20000  !! depth of upper layer [mm]
    real*4,parameter    :: d=0.4*q  !! mixing from the lower layer [mm/mon]
    real*4,parameter    :: a=-0.20  !! dT/d(d_coral)
    real*4  :: d_ow,d_owd
    real*4  :: wgt
    real*4  :: tau          !! mean transit time [years]
    real*4  :: strmax
    real*4  :: k1,k2,init
    integer :: t,yy,mm,ij,dend
    integer,parameter   :: vmiss=-999.d0
    real*4,parameter    :: r_ini=0.985
    character*10    :: arg1,arg2,arg3
    character   :: ouf*128

    open(11,file='prcp.mon.grd',access='direct',recl=4*ijmax)   ! precipitation
    open(12,file='prcp01.mon.grd',access='direct',recl=4*ijmax) ! iso in precipitation
    open(13,file='evap.mon.grd',access='direct',recl=4*ijmax)   ! evaporatin
    open(14,file='evap01.mon.grd',access='direct',recl=4*ijmax) ! iso in evaporation
    open(15,file='hadisst_sst.grd',access='direct',recl=4*ijmax)! SST
    open(16,file='ow01sfc.grd',access='direct',recl=4*ijmax)    ! iso in SST (LeGrande & Schmidt, 2006)
    open(17,file='ow01deep.grd',access='direct',recl=4*ijmax)   ! iso in SST (LeGrande & Schmidt, 2006)
    do yy=1,yymax
      do mm=1,12
        read(11,rec=(yy-1)*12+mm) (prcp(ij,mm,yy),ij=1,ijmax)
        read(12,rec=(yy-1)*12+mm) (prcp01(ij,mm,yy),ij=1,ijmax)
        read(13,rec=(yy-1)*12+mm) (evap(ij,mm,yy),ij=1,ijmax)
        read(14,rec=(yy-1)*12+mm) (evap01(ij,mm,yy),ij=1,ijmax)
        read(15,rec=(yy-1)*12+mm) (sst(ij,mm,yy),ij=1,ijmax)
      enddo
    enddo
    read(16,rec=1) (dow(ij),ij=1,ijmax)
    read(17,rec=1) (dowd(ij),ij=1,ijmax)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    !! unit
    do yy=1,yymax
      do mm=1,12
        do ij=1,ijmax
          prcp(ij,mm,yy)=prcp(ij,mm,yy)*24*60*60    !! kg/m2/s --> mm/d
          prcp01(ij,mm,yy)=prcp01(ij,mm,yy)*24*60*60    !! kg/m2/s --> mm/d
          evap(ij,mm,yy)=evap(ij,mm,yy)/2500000*24*60*60    !! W/m2 --> mm/d
          evap01(ij,mm,yy)=evap01(ij,mm,yy)/2500000*24*60*60    !! W/m2 --> mm/d
        enddo
      enddo
    enddo

    do ij=1,ijmax
      d_ow=dow(ij)
      d_owd=dowd(ij)
      do yy=1,yymax
        do mm=1,12
          call getday(dend,mm,yy)
          pterm=((prcp01(ij,mm,yy)-prcp(ij,mm,yy))*1000.d0-d_ow*prcp(ij,mm,yy))/q*dend
          eterm=((evap01(ij,mm,yy)-evap(ij,mm,yy))*1000.d0-d_ow*evap(ij,mm,yy))/q*dend
          dterm=(d_owd-d_ow)*d/q
          d_ow=d_ow+pterm-eterm+dterm
          d_coral(ij,mm,yy)=d_ow+a*(sst(ij,mm,yy)-273.15)   ! ctrl
        enddo
        d_coral2d(ij,yy)=0.d0
        do mm=1,12
          d_coral2d(ij,yy)=d_coral2d(ij,yy)+d_coral(ij,mm,yy)
        enddo
        d_coral2d(ij,yy)=d_coral2d(ij,yy)/12
      enddo
    enddo
    do ij=1,ijmax
      if (dow(ij).eq.vmiss.or.dowd(ij).eq.vmiss) then
        do yy=1,yymax
          d_coral2d(ij,yy)=vmiss
        enddo
      endif
    enddo

    ouf='coral.ann.grd'
    open(21,file=ouf,access='direct',recl=4*ijmax)
    do yy=1,yymax
      write(21,rec=yy)(d_coral2d(ij,yy),ij=1,ijmax)
    enddo
    close(21)

    end
        
      
!**********************************************************************
    subroutine getday(dd,mm,yy)
    implicit none
!input
    integer :: mm,yy
!output
    integer :: dd

    if(mm.eq.1.or.mm.eq.3.or.mm.eq.5.or.mm.eq.7.or.mm.eq.8.or.mm.eq.10.or.mm.eq.12)then
      dd=31
    elseif(mm.eq.4.or.mm.eq.6.or.mm.eq.9.or.mm.eq.11)then
      dd=30
    else
      if(mod(yy,4).eq.0)then
        if(mod(yy,100).eq.0)then
          if(mod(yy,400).eq.0)then
            dd=29
          else
            dd=28
          endif
        else
          dd=29
        endif
      else
        dd=28
      endif
    endif

    return
    end subroutine
