*Roden Model for MIROC5
*Calculates the del18O value of cellulose based on a user-defined growing season
target='celld.ann.bin'
ystr=1871
yend=2007
stt=1
edt=(yend-ystr)*12+12
tt=stt

*Data
dir=''
ptop=100
nppfile=''

'!rm 'target

* This section should be modified by user
* ---------------------------------------
* net primary production (NPP)
'open 'nppfile
* isotope ratio in vapor (kg/kg)
'open 'dir'/q01'
* vapor (kg/kg)
'open 'dir'/qt'
* surface pressure (hPa)
'open 'dir'/Ps'
* 2-m temperature (K)
'open 'dir'/T2'
* 2-m relative humidity
'open 'dir'/rh2'
* isotope ratio in precipitation (kg/m**2/s)
'open 'dir'/prcp01'
* recipitation (W/m**2)
'open 'dir'/prcp'
'open 'dir'/evap01'
'open 'dir'/evap'


'set dfile 1'
'set lon 0 360'
'set lat -90 90'
'set t 1 12'
'set z 1'
'define npx=lterp(npp,t2.5(t=1))'
'modify npx seasonal'

'set dfile 2'
'set lon 0 360'
'set lat -90 90'
'set t 1'
'npxan=sum(npx,t=1,t=12)'
'set lev 1000'
'set t 'stt' 'edt

**define variables name corresponding to IsoGSM
'define pwat1clm=vint(ps.4,q01.2,'ptop')'
'define pwatclm=vint(ps.4,qt.3,'ptop')'
'define pressfc=ps.4*100'
'define tmp2m=t2.5'
'define rhx=rh2.6*100'
'define prate1sfc=prcp01.7'
'define pratesfc=prcp.8'

*'define dvx=lterp(tr2.1/spfh.1,tmp2m)'
'set z 1'

**define first month of growing season, enter month after "="
gsbn=6
'define gsb='gsbn
**define last month of growing season, enter month after "="
gsen=9
'define gse='gsen
**definition of "water year", which spans from the month after the end of the growing 
**season to the end of the next growing season
wybn=gsen-11
wyen=gsbn
'define wyb=gse-11'
'define wye=gsb'
**define stomatal conductance (mol*m-2*s-1)
'define sc=0.0001'
**define boundary layer conductance (mol*m-2*s-1)
'define bc=0.0004'
*define difference between air temperature and leaf temperature
*'define dt=6'
'define dt=1'

* This section should not be modifed by user!
* -------------------------------------------
'define ak=1.032'
'define kb=1.021'
'define ro=0.0020052'
'define ff=27'
'define fre=42/100'

*define growing season temperature
'define tmp1=tmp2m'
'define tmp=tmp1-273.1'

*define growing season temperature at leaf surface
'define atmp=tmp+dt'

*define growing season relative humidity
'define rh=rhx'

*define growing season barometric pressure
'define bp=pressfc'

*define average growing season atmospheric water del18O value
'define dl=(pwat1clm/pwatclm-1)'
'define dl1=dl*1000'
'define dv=dl1'

**define weighted average source water isotopic value, which is 
**the amount-weighted del18O value from the end of the previous
**year's growing season to the end of the current growing season
'define sw=(ave(prate1sfc,t-3,t+0)/ave(pratesfc,t-3,t+0)-1)*1000'

*Calculate Leaf vapor pressure (ei) in kPa
'define ei=(101325*EXP((((-0.1299*(1-(373.16/(273.16+atmp)))-0.6445)*(1-(373.16/(273.16+atmp)))-1.976)*(1-(373.16/(273.16+atmp)))+13.3185)*(1-(373.16/(273.16+atmp)))))/1000'
*Calculate saturation vapor pressure (esat) in kPa
'define esat=(101325*EXP((((-0.1299*(1-(373.16/(273.16+tmp)))-0.6445)*(1-(373.16/(273.16+tmp)))-1.976)*(1-(373.16/(273.16+tmp)))+13.3185)*(1-(373.16/(273.16+tmp)))))/1000'
*calculate ambient vapor pressure (ea) in kPa
'define ea=(rh/100)*esat'
*Calculate leaf conductance of water vapor (gi) in mol*m-2*s-1
'define gi=1/(1/sc+1/bc)'
*calculate leaf transpiration (et) in mol*m-2*s-1
'define et=((ei-ea)/bp)*gi'
*calculate leaf water vapor fraction (wi) as a molar fraction
'define wi=ei/bp'
*calculate ambient water vapor (wa) as a molar fraction
'define wa=ea/bp'
*calculate leaf surface water vapor (ws) as a molar fraction
'define ws=((sc*wi)-et*(1 -wi/2))/(sc - et/2)'
*calculate vapor pressure at leaf surface (es) in kPa
'define es=ws*bp'
*calculate temperature-dependent equilibrium fractionation factor between vapor and liquid water at leaf surface (ef)
'define ef=EXP((1.137*1000/(273.16+atmp)/(273.16+atmp))-(0.4156/(273.16+atmp))-0.0020667)'
*calculate isotopic composition of source water as an isotopic ratio, R-value (rsw)
'define rsw=ro*(1+sw/1000)'
*calculate isotopic composition of atmospheric water vapor as an isotopic ratio, R-value, (ratm)
'define ratm=ro*(1+dv/1000)'
*calculate isotopic composition of leaf water as an isotopic ratio, R-value (lho)
'define lho=ef*((ak*rsw*((ei-es)/ei))+(kb*rsw*(es-ea)/ei)+(ratm*ea/ei))'
*calculate delta value for leaf water (ldho)
'define ldho=((lho/ro)-1)*1000'
*calculate delta value for annual integrated cellulose (celld)
'define celld=fre*(sw+ff)+(1-fre)*(ldho+ff)'

*weight by npp
'define celldw=ave(celld*npx,t+0,t+11)/ave(npx,t+0,t+11)'


* Output
* ------
'set gxout fwrite'
'set fwrite 'target
'set x 1 128'
'set y 1 64'
'set z 1'

tt=stt
while tt <= edt
'set t 'tt
'd celldw'
tt=tt+12
endwhile

'disable fwrite'

'quit'
