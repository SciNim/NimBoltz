## COMMON blocks of magboltz-11.21.f, translated by f2nim; maintained by hand since.

import ../frt, ../tables
export tables

const nscale* {.intdefine.} = 40_000_000
  ## NSCALE of SETUP/SETUPT: collisions per unit of NMAX in the input deck.
  ## The regression suite builds with -d:nscale=4000000 for speed.

type
  CommonBfld* = object  ## COMMON /BFLD/
    eovb*: float64
    wb*: float64
    btheta*: float64
    bmag*: float64
  CommonInpt* = object  ## COMMON /INPT/
    ngas*: int
    nstep*: int
    naniso*: int
    efinal*: float64
    estep*: float64
    akt*: float64
    ary*: float64
    tempc*: float64
    torr*: float64
    ipen*: int
  CommonSetp* = object  ## COMMON /SETP/
    tmax*: float64
    small*: float64
    api*: float64
    estart*: float64
    theta*: float64
    phi*: float64
    rstart*: float64
    efield*: float64
    nmax*: int
  CommonCtowns* = object  ## COMMON /CTOWNS/
    alpha*: float64
    att*: float64
  CommonThrm* = object  ## COMMON /THRM/
    amgas*: array[1..6, float64]
    vtmb*: array[1..6, float64]
    tcfmx*: float64
    tcfmxg*: array[1..6, float64]
    ithrm*: int
  CommonNames* = object  ## COMMON /NAMES/
    nameg*: array[1..6, string]
  CommonScript* = object  ## COMMON /SCRIPT/
    dscrpt*: Arr2[string, 6, 300]
  CommonCnsts* = object  ## COMMON /CNSTS/
    echarg*: float64
    emass*: float64
    amu*: float64
    pir2*: float64
  CommonRatio* = object  ## COMMON /RATIO/
    an1*: float64
    an2*: float64
    an3*: float64
    an4*: float64
    an5*: float64
    an6*: float64
    an*: float64
    frac*: array[1..6, float64]
  CommonGasn* = object  ## COMMON /GASN/
    ngasn*: array[1..6, int]
  CommonMix2* = object  ## COMMON /MIX2/ (energy grid; cross sections dropped)
    eg*: array[1..4000, float64]
    eroot*: array[1..4000, float64]
  CommonCnsts1* = object  ## COMMON /CNSTS1/
    const1*: float64
    const2*: float64
    const3*: float64
    const4*: float64
    const5*: float64
  CommonLarget* = object  ## COMMON /LARGET/
    tcfmax*: array[1..6, float64]
  CommonAnist* = object  ## COMMON /ANIST/
    niso*: int
  CommonFred* = object  ## COMMON /FRED/
    fcion*: array[1..4000, float64]
    fcatt*: array[1..4000, float64]
  CommonMratio* = object  ## COMMON /MRATIO/
    van1*: float64
    van2*: float64
    van3*: float64
    van4*: float64
    van5*: float64
    van6*: float64
    van*: float64
  CommonVel* = object  ## COMMON /VEL/
    wx*: float64
    wy*: float64
    wz*: float64
  CommonVelerr* = object  ## COMMON /VELERR/
    dwx*: float64
    dwy*: float64
    dwz*: float64
  CommonCtwner* = object  ## COMMON /CTWNER/
    alper*: float64
    atter*: float64
  CommonDifvel* = object  ## COMMON /DIFVEL/
    difln*: float64
    diftr*: float64
  CommonDiferl* = object  ## COMMON /DIFERL/
    dfler*: float64
    dfter*: float64
  CommonDiflab* = object  ## COMMON /DIFLAB/
    difxx*: float64
    difyy*: float64
    difzz*: float64
    difyz*: float64
    difxy*: float64
    difxz*: float64
  CommonDiferb* = object  ## COMMON /DIFERB/
    dxxer*: float64
    dyyer*: float64
    dzzer*: float64
    dyzer*: float64
    dxyer*: float64
    dxzer*: float64
  CommonGer* = object  ## COMMON /GER/
    twopi*: float64
  CommonTtrm* = object  ## COMMON /TTRM/
    ztot*: float64
    ttot*: float64
    ztots*: float64
    ttots*: float64
  CommonOutptt* = object  ## COMMON /OUTPTT/
    time*: array[1..300, float64]
    icoll*: Arr2[int, 6, 5]
    spec*: array[1..4000, float64]
    tmax1*: float64
    ave*: float64
    den*: float64
    xid*: float64
    x*: float64
    y*: float64
    z*: float64
    st*: float64
    nnull*: int
    icoln*: Arr2[int, 6, 300]
  CommonDecor* = object  ## COMMON /DECOR/
    ncolm*: int
    ncorln*: int
    ncorst*: int
  CommonDens* = object  ## COMMON /DENS/
    den*: array[1..4000, float64]
  CommonRanm* = object  ## COMMON /RANM/
    rnmx*: array[1..6, float64]
  CommonSint* = object  ## COMMON /SINT/
    simf*: array[1..4000, float64]
  CommonSstout* = object  ## COMMON /SSTOUT/
    vdout*: float64
    vderr*: float64
    wsout*: float64
    wserr*: float64
    dlout*: float64
    dlerr*: float64
    dtout*: float64
    dterr*: float64
    alphsst*: float64
    alpherr*: float64
    attsst*: float64
    atterr*: float64
  CommonTofout* = object  ## COMMON /TOFOUT/
    ralpha*: float64
    ralper*: float64
    tofene*: float64
    tofener*: float64
    tofwv*: float64
    tofwver*: float64
    tofdl*: float64
    tofdler*: float64
    tofdt*: float64
    tofdter*: float64
    tofwr*: float64
    tofwrer*: float64
    rattof*: float64
    ratofer*: float64
  CommonCtcalc* = object  ## COMMON /CTCALC/
    zplane1*: float64
    zplane2*: float64
    zplane3*: float64
    zplane4*: float64
    zplane5*: float64
    zplane6*: float64
    zplane7*: float64
    zplane8*: float64
    izfinal*: int
  CommonCion* = object  ## COMMON /CION/
    alphast*: float64
    vdst*: float64
    tstep*: float64
    zstep*: float64
    tfinal*: float64
    zfinal*: float64
    itfinal*: int
    iprim*: int
  CommonTplout* = object  ## COMMON /TPLOUT/
    etpl*: array[1..8, float64]
    xtpl*: array[1..8, float64]
    ytpl*: array[1..8, float64]
    ztpl*: array[1..8, float64]
    ttpl*: array[1..8, float64]
    xxtpl*: array[1..8, float64]
    yytpl*: array[1..8, float64]
    zztpl*: array[1..8, float64]
    vztpl*: array[1..8, float64]
    netpl*: array[1..8, int]
    attoint*: float64
    attert*: float64
    aioert*: float64
  CommonPttof* = object  ## COMMON /PTTOF/
    ri*: array[1..8, float64]
    ept*: array[1..8, float64]
    vzpt*: array[1..8, float64]
    ttest*: array[1..8, float64]
  CommonSplout* = object  ## COMMON /SPLOUT/
    espl*: array[1..8, float64]
    xspl*: array[1..8, float64]
    yspl*: array[1..8, float64]
    zspl*: array[1..8, float64]
    tspl*: array[1..8, float64]
    xxspl*: array[1..8, float64]
    yyspl*: array[1..8, float64]
    zzspl*: array[1..8, float64]
    vzspl*: array[1..8, float64]
    tssum*: array[1..8, float64]
    tssum2*: array[1..8, float64]
    attoion*: float64
    attioer*: float64
    attater*: float64
    nesst*: array[1..9, int]
  CommonSpl1* = object  ## COMMON /SPL1/
    tmspl*: array[1..8, float64]
    ttmspl*: array[1..8, float64]
    rspl*: array[1..8, float64]
    rrspl*: array[1..8, float64]
    rrspm*: array[1..8, float64]
  CommonTofgout* = object  ## COMMON /TOFGOUT/
    ralpha*: float64
    ralper*: float64
    tofene*: float64
    tofener*: float64
    tofwvz*: float64
    tofwvzer*: float64
    tofwvy*: float64
    tofwvyer*: float64
    tofdzz*: float64
    tofdzzer*: float64
    tofdxx*: float64
    tofdxxer*: float64
    tofdyy*: float64
    tofdyyer*: float64
    tofdyz*: float64
    tofdyzer*: float64
    tofwrz*: float64
    tofwrzer*: float64
    tofwry*: float64
    tofwryer*: float64
    rattof*: float64
    ratofer*: float64
  CommonTploutg* = object  ## COMMON /TPLOUTG/
    etpl*: array[1..8, float64]
    xtpl*: array[1..8, float64]
    ytpl*: array[1..8, float64]
    ztpl*: array[1..8, float64]
    ttpl*: array[1..8, float64]
    xxtpl*: array[1..8, float64]
    yytpl*: array[1..8, float64]
    zztpl*: array[1..8, float64]
    yztpl*: array[1..8, float64]
    vztpl*: array[1..8, float64]
    vytpl*: array[1..8, float64]
    netpl*: array[1..8, int]
    attoint*: float64
    attert*: float64
    aioert*: float64
  CommonPttofg* = object  ## COMMON /PTTOFG/
    ri*: array[1..8, float64]
    ept*: array[1..8, float64]
    vzpt*: array[1..8, float64]
    vypt*: array[1..8, float64]
    ttest*: array[1..8, float64]
  CommonTofhout* = object  ## COMMON /TOFHOUT/
    ralpha*: float64
    ralper*: float64
    tofene*: float64
    tofener*: float64
    tofwvz*: float64
    tofwvzer*: float64
    tofwvy*: float64
    tofwvyer*: float64
    tofwvx*: float64
    tofwvxer*: float64
    tofdzz*: float64
    tofdzzer*: float64
    tofdxx*: float64
    tofdxxer*: float64
    tofdyy*: float64
    tofdyyer*: float64
    tofdyz*: float64
    tofdyzer*: float64
    tofdxy*: float64
    tofdxyer*: float64
    tofdxz*: float64
    tofdxzer*: float64
    tofwrz*: float64
    tofwrzer*: float64
    tofwry*: float64
    tofwryer*: float64
    tofwrx*: float64
    tofwrxer*: float64
    rattof*: float64
    ratofer*: float64
  CommonTplouth* = object  ## COMMON /TPLOUTH/
    etpl*: array[1..8, float64]
    xtpl*: array[1..8, float64]
    ytpl*: array[1..8, float64]
    ztpl*: array[1..8, float64]
    ttpl*: array[1..8, float64]
    xxtpl*: array[1..8, float64]
    yytpl*: array[1..8, float64]
    zztpl*: array[1..8, float64]
    yztpl*: array[1..8, float64]
    xztpl*: array[1..8, float64]
    xytpl*: array[1..8, float64]
    vztpl*: array[1..8, float64]
    vytpl*: array[1..8, float64]
    vxtpl*: array[1..8, float64]
    netpl*: array[1..8, int]
    attoint*: float64
    attert*: float64
    aioert*: float64
  CommonPttofh* = object  ## COMMON /PTTOFH/
    ri*: array[1..8, float64]
    ept*: array[1..8, float64]
    vzpt*: array[1..8, float64]
    vypt*: array[1..8, float64]
    vxpt*: array[1..8, float64]
    ttest*: array[1..8, float64]
  CommonScrip* = object  ## COMMON /SCRIP/
    dscrpt*: array[1..960, string]
  CommonLarge* = object  ## COMMON /LARGE/
    tcfmax*: array[1..8, float64]
  CommonAnis* = object  ## COMMON /ANIS/
    niso*: int
  CommonOutpt* = object  ## COMMON /OUTPT/
    time*: array[1..300, float64]
    icoll*: array[1..30, int]
    spec*: array[1..4000, float64]
    tmax1*: float64
    ave*: float64
    den*: float64
    xid*: float64
    x*: float64
    y*: float64
    z*: float64
    st*: float64
    nnull*: int
    icoln*: array[1..960, int]

var
  bfld*: CommonBfld
  inpt*: CommonInpt
  setp*: CommonSetp
  ctowns*: CommonCtowns
  thrm*: CommonThrm
  names*: CommonNames
  script*: CommonScript
  cnsts*: CommonCnsts
  ratio*: CommonRatio
  gasn*: CommonGasn
  mix2*: CommonMix2
  cnsts1*: CommonCnsts1
  larget*: CommonLarget  ## rest in gasTable
  gasTable*: array[1..maxGases, CollisionTable[maxGasProcesses]]  ## /LARGET/ + /ANIST/
  anist*: CommonAnist
  fred*: CommonFred
  mratio*: CommonMratio
  vel*: CommonVel
  velerr*: CommonVelerr
  ctwner*: CommonCtwner
  difvel*: CommonDifvel
  diferl*: CommonDiferl
  diflab*: CommonDiflab
  diferb*: CommonDiferb
  ger*: CommonGer
  ttrm*: CommonTtrm
  outptt*: CommonOutptt
  decor*: CommonDecor
  dens*: CommonDens
  ranm*: CommonRanm
  sint*: CommonSint
  sstout*: CommonSstout
  tofout*: CommonTofout
  ctcalc*: CommonCtcalc
  cion*: CommonCion
  tplout*: CommonTplout
  pttof*: CommonPttof
  splout*: CommonSplout
  spl1*: CommonSpl1
  tofgout*: CommonTofgout
  tploutg*: CommonTploutg
  pttofg*: CommonPttofg
  tofhout*: CommonTofhout
  tplouth*: CommonTplouth
  pttofh*: CommonPttofh
  scrip*: CommonScrip
  large*: CommonLarge  ## rest in mixTable
  mixTable*: CollisionTable[maxMixProcesses]  ## /LARGE/ + /ANIS/
  anis*: CommonAnis
  outpt*: CommonOutpt
