# -----------------------------------------------------------------------
# Program: ETC88.R  
#  Author: Hermine Maes
#    Date: 01 07 2010 
#
# Estimating Extended Twin-kinship (ET) ML Correlations
#
# Revision History
#   Hermine Maes -- 01 07 2010 updated & reformatted
#   Hermine Maes -- 08 31 2010 updated with changes from EC88ccT5.R
#   Hermine Maes -- 09 15 2010 updated with changes from EC88cc2TST5.R
# -----------------------------------------------------------------------

require(OpenMx)

# Prepare Data
# -----------------------------------------------------------------------!
EtVars <- c('famid','e1','e2','e3','e4','e5','e6','e7','e8','e9','e10',
	'e11','e12','e13','e14','e15','e16','e17','e18',
	'a1','a2','a3','a4','a5','a6','a7','a8','a9',
	'a10','a11','a12','a13','a14','a15','a16','a17','a18')
data(mzmData)
data(dzmData)
data(mzfData)
data(dzfData)
data(dzoData)

dataMZM <- mzmData[2:19]
dataDZM <- dzmData[2:19]
dataMZF <- mzfData[2:19]
dataDZF <- dzfData[2:19]
dataDZO <- dzoData[2:19]

# Prepare Labels
# -----------------------------------------------------------------------!
selVars <-c('e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18')
Mpsi <- c('FA','MO','MS','MS','FS','FS')
Mtwc <- c('MC1','MC1','FC1','FC1','MC2','MC2','FC2','FC2')
labMeansM <- c('MT1','MT2',Mpsi,'FS1','FS2',Mtwc) 
labMeansF <- c('FT1','FT2',Mpsi,'MS1','MS2',Mtwc) 
labMeansMF <- c('MT1','FT2',Mpsi,'FS1','MS2',Mtwc)
labMeans <- c('MT1','FT1','MT2','FT2','FA','MO','MS','FS','MS1','FS1','MS2','FS2','MC1','FC1','MC2','FC2')
labMeans10 <- c('MT','FT','FA','MO','MS','FS','MP','FP','MC','FC')

SDpsi <- c('SDm','SDf','SDm','SDm','SDf','SDf')
SDtwc <- c('SDm','SDm','SDf','SDf','SDm','SDm','SDf','SDf')
labSDevM <- c('SDm','SDm',SDpsi,'SDf','SDf',SDtwc)
labSDevF <- c('SDf','SDf',SDpsi,'SDm','SDm',SDtwc)
labSDevMF <- c('SDm','SDf',SDpsi,'SDf','SDm',SDtwc)
labSDevs <- c('SDm','SDf')

# Prepare Start Values
# -----------------------------------------------------------------------!
svMZ <- 0.8
svDZ <- 0.4
svSMZ <- 0.1
svSDZ <- 0.1
svSSMZ <- 0.1
svSSDZ <- 0.1
svSI <- 0.4
svSSI <- 0.1
svSP <- 0.1
svPO <- 0.4
svSPO <- 0.1
svGP <- 0.2
svAVMZ <- 0.2
svAVSI <- 0.1
svAVDZ <- 0.1
svSAVMZ <- 0.1
svSAVDZ <- 0.1
svCOMZ <- 0.2
svCODZ <- 0.1 
 
ns <- 2	# number of sibs/kids
nvr <- 18	# number of variables in expected means/covariances
nv <- 1	# number of variables analyzed

# Specify Model Matrices and Algebra
# -----------------------------------------------------------------------!
etc88Model <- mxModel("et",

    mxMatrix("Full", 1, 1, F, 9, name="x"),
    mxMatrix("Unit", ns,ns, name="xx_xx"),

    mxMatrix("Full", nv, nv, F, 1, name="mv"),
    mxMatrix("Full", nv, nv, F, 1, name="fv"),

    mxMatrix("Full", nv, nv, T, svMZ, "r1", -.2, .9, name="mzm"),
    mxMatrix("Full", nv, nv, T, svDZ, "r2", -.2, .9, name="dzm"),
    mxMatrix("Full", nv, nv, T, svMZ, "r3", -.2, .9, name="mzf"),
    mxMatrix("Full", nv, nv, T, svDZ, "r4", -.2, .9, name="dzf"),
    mxMatrix("Full", nv, nv, T, svDZ, "r5", -.2, .9, name="dzmf"),
    mxMatrix("Full", nv, nv, T, svSMZ, "r6", -.2, .9, name="mzmw"),
    mxMatrix("Full", nv, nv, T, svSDZ, "r7", -.2, .9, name="dzmw"),
    mxMatrix("Full", nv, nv, T, svSMZ, "r8", -.2, .9, name="mzfh"),
    mxMatrix("Full", nv, nv, T, svSDZ, "r9", -.2, .9, name="dzfh"),
    mxMatrix("Full", nv, nv, T, svSDZ, "r10", -.2, .9, name="dzmfh"),
    mxMatrix("Full", nv, nv, T, svSDZ, "r11", -.2, .9, name="dzfmw"),
    mxMatrix("Full", nv, nv, T, svSSMZ, "r12", -.2, .9, name="wmzmw"),
    mxMatrix("Full", nv, nv, T, svSSDZ, "r13", -.2, .9, name="wdzmw"),
    mxMatrix("Full", nv, nv, T, svSSMZ, "r14", -.2, .9, name="hmzfh"),
    mxMatrix("Full", nv, nv, T, svSSDZ, "r15", -.2, .9, name="hdzfh"),
    mxMatrix("Full", nv, nv, T, svSSDZ, "r16", -.2, .9, name="wdzmfh"),

    mxMatrix("Full", nv, nv, T, svSI, "r17", -.2, .9, name="sim"),
    mxMatrix("Full", nv, nv, T, svSI, "r18", -.2, .9, name="simf"),
    mxMatrix("Full", nv, nv, T, svSI, "r19", -.2, .9, name="sif"),
    mxMatrix("Full", nv, nv, T, svSSI, "r20", -.2, .9, name="simw"),
    mxMatrix("Full", nv, nv, T, svSSI, "r21", -.2, .9, name="sifmw"),
    mxMatrix("Full", nv, nv, T, svSSI, "r22", -.2, .9, name="simfh"),
    mxMatrix("Full", nv, nv, T, svSSI, "r23", -.2, .9, name="sifh"),
    
    mxAlgebra(rbind(cbind(mv,sim,simf,simf), cbind(t(sim),mv,simf,simf),
        cbind(t(simf),t(simf),fv,sif), cbind(t(simf),t(simf),t(sif),fv)), name="bs"),  #sibs
    mxAlgebra(cbind(sim,sim,simf,simf), name="mtsi"),  #male twin-sibs
    mxAlgebra(cbind(simf,simf,sif,sif), name="ftsi"),  #female twin-sibs
    mxAlgebra(rbind(simw,simw,sifmw,sifmw), name="wmtsi"),  #male twin wives-sibs
    mxAlgebra(rbind(simfh,simfh,sifh,sifh), name="hftsi"),  #female twin husbs-sibs

    mxMatrix("Full", nv, nv, T, svSP, "r24", -.2, .9, name="mf"),
    mxMatrix("Full", nv, nv, T, svPO, "r25", -.2, .9, name="fs"),
    mxMatrix("Full", nv, nv, T, svPO, "r26", -.2, .9, name="fd"),
    mxMatrix("Full", nv, nv, T, svPO, "r27", -.2, .9, name="ms"),
    mxMatrix("Full", nv, nv, T, svPO, "r28", -.2, .9, name="md"),
    mxMatrix("Full", nv, nv, T, svSPO, "r29", -.2, .9, name="fsw"),
    mxMatrix("Full", nv, nv, T, svSPO, "r30", -.2, .9, name="msw"),
    mxMatrix("Full", nv, nv, T, svSPO, "r31", -.2, .9, name="fdh"),
    mxMatrix("Full", nv, nv, T, svSPO, "r32", -.2, .9, name="mdh"),

    mxAlgebra(rbind(cbind(mv,mf),cbind(t(mf),fv)), name="sp"),  #spouse
    mxAlgebra(cbind(fs,ms), name="ps"),  #parent-son
    mxAlgebra(cbind(fd,md), name="pd"),  #parent-daughter
    mxAlgebra(cbind(fs,fs,fd,fd), name="fo"),  #father-offspring
    mxAlgebra(cbind(ms,ms,md,md), name="mo"),  #mother-offspring
    mxAlgebra(rbind(fsw,msw), name="psil"),  #parent-son in law
    mxAlgebra(rbind(fdh,mdh), name="pdil"),  #parent-daughter in law
    mxAlgebra(rbind(fo,mo), name="fmo"),   #father_mother-offspring

    mxMatrix("Full", nv, nv, T, svGP, "r33", -.2, .9, name="gfps"),
    mxMatrix("Full", nv, nv, T, svGP, "r34", -.2, .9, name="gmps"),
    mxMatrix("Full", nv, nv, T, svGP, "r35", -.2, .9, name="gfpd"),
    mxMatrix("Full", nv, nv, T, svGP, "r36", -.2, .9, name="gmpd"),
    mxMatrix("Full", nv, nv, T, svGP, "r37", -.2, .9, name="gfms"),
    mxMatrix("Full", nv, nv, T, svGP, "r38", -.2, .9, name="gmms"),
    mxMatrix("Full", nv, nv, T, svGP, "r39", -.2, .9, name="gfmd"),
    mxMatrix("Full", nv, nv, T, svGP, "r40", -.2, .9, name="gmmd"),
    mxMatrix("Full", nv, nv, T, svAVMZ, "r41", -.2, .9, name="unemzm"),
    mxMatrix("Full", nv, nv, T, svAVMZ, "r42", -.2, .9, name="unimzm"),
    mxMatrix("Full", nv, nv, T, svAVMZ, "r43", -.2, .9, name="anemzf"),
    mxMatrix("Full", nv, nv, T, svAVMZ, "r44", -.2, .9, name="animzf"),

    mxAlgebra(rbind(cbind(gfps,gfps,gfpd,gfpd),
        cbind(gmps,gmps,gmpd,gmpd)), name="gpcp"),  #grandparent-child paternal
    mxAlgebra(rbind(cbind(gfms,gfms,gfmd,gfmd),
        cbind(gmms,gmms,gmmd,gmmd)), name="gpcm"),  #grandparent-child maternal
    mxAlgebra(cbind(unemzm,unemzm,unimzm,unimzm), name="umzm"),  #uncles via mzm twins
    mxAlgebra(cbind(anemzf,anemzf,animzf,animzf), name="amzf"),  #uncles via mzf twins

    mxMatrix("Full", nv, nv, T, svAVSI, "r45", -.2, .9, name="unesim"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r46", -.2, .9, name="unisim"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r47", -.2, .9, name="anesifm"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r48", -.2, .9, name="anisifm"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r49", -.2, .9, name="unesimf"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r50", -.2, .9, name="unisimf"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r51", -.2, .9, name="anesif"),
    mxMatrix("Full", nv, nv, T, svAVSI, "r52", -.2, .9, name="anisif"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r53", -.2, .9, name="unedzm"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r54", -.2, .9, name="unidzm"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r55", -.2, .9, name="anedzf"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r56", -.2, .9, name="anidzf"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r57", -.2, .9, name="unedzmf"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r58", -.2, .9, name="unidzmf"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r59", -.2, .9, name="anedzmf"),
    mxMatrix("Full", nv, nv, T, svAVDZ, "r60", -.2, .9, name="anidzmf"),

    mxAlgebra(rbind(cbind(unesim,unisim),cbind(anesifm,anisifm)) %x% xx_xx, name="avsim"),  #avuncular via sim
    mxAlgebra(rbind(cbind(unesimf,unisimf),cbind(anesif,anisif)) %x% xx_xx, name="avsif"),  #avuncular via sif
    mxAlgebra(cbind(unedzm,unedzm,unidzm,unidzm), name="udzm"),  #uncles via dzm twins
    mxAlgebra(cbind(anedzf,anedzf,anidzf,anidzf), name="adzf"),  #aunts via dzf twins
    mxAlgebra(cbind(unedzmf,unedzmf,unidzmf,unidzmf), name="udzmf"),  #uncles via dzmf twins
    mxAlgebra(cbind(anedzmf,anedzmf,anidzmf,anidzmf), name="adzmf"),  #aunts via dzmf twins

    mxMatrix("Full", nv, nv, T, svSAVMZ, "r61", -.2, .9, name="wunemzm"),
    mxMatrix("Full", nv, nv, T, svSAVMZ, "r62", -.2, .9, name="wunimzm"),
    mxMatrix("Full", nv, nv, T, svSAVMZ, "r63", -.2, .9, name="hanemzf"),
    mxMatrix("Full", nv, nv, T, svSAVMZ, "r64", -.2, .9, name="hanimzf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r65", -.2, .9, name="wunedzm"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r66", -.2, .9, name="wunidzm"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r67", -.2, .9, name="hanedzf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r68", -.2, .9, name="hanidzf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r69", -.2, .9, name="wunedzmf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r70", -.2, .9, name="wunidzmf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r71", -.2, .9, name="hanedzmf"),
    mxMatrix("Full", nv, nv, T, svSAVDZ, "r72", -.2, .9, name="hanidzmf"),

    mxAlgebra(cbind(wunemzm,wunemzm,wunimzm,wunimzm), name="amzm"),  #aunts via mzm twins
    mxAlgebra(cbind(hanemzf,hanemzf,hanimzf,hanimzf), name="umzf"),  #uncles via mzf twins
    mxAlgebra(cbind(wunedzm,wunedzm,wunidzm,wunidzm), name="adzm"),  #aunts via dzm twins
    mxAlgebra(cbind(hanedzf,hanedzf,hanidzf,hanidzf), name="udzf"),  #uncles via dzf twins
    mxAlgebra(cbind(wunedzmf,wunedzmf,wunidzmf,wunidzmf), name="aidzmf"),  #aunts via dzmf twins
    mxAlgebra(cbind(hanedzmf,hanedzmf,hanidzmf,hanidzmf), name="uidzmf"),  #uncles via dzmf twins

    mxMatrix("Full", nv, nv, T, svCOMZ, "r73", -.2, .9, name="commzm"),
    mxMatrix("Full", nv, nv, T, svCOMZ, "r74", -.2, .9, name="comfmzm"),
    mxMatrix("Full", nv, nv, T, svCOMZ, "r75", -.2, .9, name="cofmzm"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r76", -.2, .9, name="comdzm"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r77", -.2, .9, name="comfdzm"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r78", -.2, .9, name="cofdzm"),
    mxMatrix("Full", nv, nv, T, svCOMZ, "r79", -.2, .9, name="commzf"),
    mxMatrix("Full", nv, nv, T, svCOMZ, "r80", -.2, .9, name="comfmzf"),
    mxMatrix("Full", nv, nv, T, svCOMZ, "r81", -.2, .9, name="cofmzf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r82", -.2, .9, name="comdzf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r83", -.2, .9, name="comfdzf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r84", -.2, .9, name="cofdzf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r85", -.2, .9, name="comdzmf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r86", -.2, .9, name="comfdzmf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r87", -.2, .9, name="cofmdzmf"),
    mxMatrix("Full", nv, nv, T, svCODZ, "r88", -.2, .9, name="cofdzmf"),

    mxAlgebra(rbind(cbind(commzm,comfmzm),cbind(comfmzm,cofmzm)) %x% xx_xx, name="comzm"),  #cousins mzm
    mxAlgebra(rbind(cbind(comdzm,comfdzm),cbind(comfdzm,cofdzm)) %x% xx_xx, name="codzm"),  #cousins dzm
    mxAlgebra(rbind(cbind(commzf,comfmzf),cbind(comfmzf,cofmzf)) %x% xx_xx, name="comzf"),  #cousins mzf
    mxAlgebra(rbind(cbind(comdzf,comfdzf),cbind(comfdzf,cofdzf)) %x% xx_xx, name="codzf"),  #cousins dzf
    mxAlgebra(rbind(cbind(comdzmf,comfdzmf),cbind(cofmdzmf,cofdzmf)) %x% xx_xx, name="codzmf"),  #cousins dzmf

# t1       t2       pa(fm)   si        s1        s2       k1        k2
# 1        1        2        4         1         1        4         4
    mxAlgebra(rbind(
        cbind( mv     , mzm    , ps     , mtsi    , mf      , mzmw   , fo      , umzm  ),
        cbind( t(mzm) , mv     , ps     , mtsi    , mzmw    , mf     , umzm    , fo    ),
        cbind( t(ps)  , t(ps)  , sp     , fmo     , psil    , psil   , gpcp    , gpcp  ),
        cbind( t(mtsi), t(mtsi), t(fmo) , bs      , wmtsi   , wmtsi  , avsim   , avsim ),
        cbind( t(mf)  , t(mzmw), t(psil), t(wmtsi), fv      , wmzmw  , mo      , amzm  ), 
        cbind( t(mzmw), t(mf)  , t(psil), t(wmtsi), t(wmzmw), fv     , amzm    , mo    ),
        cbind( t(fo)  , t(umzm), t(gpcp), t(avsim), t(mo)   , t(amzm), bs      , comzm ),
        cbind( t(umzm), t(fo)  , t(gpcp), t(avsim), t(amzm) , t(mo)  , t(comzm), bs    )), name="corMZM"),
    mxAlgebra(rbind(
        cbind( mv     , dzm    , ps     , mtsi    , mf      , dzmw   , fo      , udzm  ),
        cbind( t(dzm) , mv     , ps     , mtsi    , dzmw    , mf     , udzm    , fo    ),
        cbind( t(ps)  , t(ps)  , sp     , fmo     , psil    , psil   , gpcp    , gpcp  ),
        cbind( t(mtsi), t(mtsi), t(fmo) , bs      , wmtsi   , wmtsi  , avsim   , avsim ),
        cbind( t(mf)  , t(dzmw), t(psil), t(wmtsi), fv      , wdzmw  , mo      , adzm  ), 
        cbind( t(dzmw), t(mf)  , t(psil), t(wmtsi), t(wdzmw), fv     , adzm    , mo    ),
        cbind( t(fo)  , t(udzm), t(gpcp), t(avsim), t(mo)   , t(adzm), bs      , codzm ),
        cbind( t(udzm), t(fo)  , t(gpcp), t(avsim), t(adzm) , t(mo)  , t(codzm), bs    )), name="corDZM"),
    mxAlgebra(rbind(
        cbind( fv     , mzf    , pd     , ftsi    , mf      , mzfh   , mo      , amzf  ),
        cbind( t(mzf) , fv     , pd     , ftsi    , mzfh    , mf     , amzf    , mo    ),
        cbind( t(pd)  , t(pd)  , sp     , fmo     , pdil    , pdil   , gpcm    , gpcm  ),
        cbind( t(ftsi), t(ftsi), t(fmo) , bs      , hftsi   , hftsi  , avsif   , avsif ),
        cbind( t(mf)  , t(mzfh), t(pdil), t(hftsi), mv      , hmzfh  , fo      , umzf  ), 
        cbind( t(mzfh), t(mf)  , t(pdil), t(hftsi), t(hmzfh), mv     , umzf    , fo    ),
        cbind( t(mo)  , t(amzf), t(gpcm), t(avsif), t(fo)   , t(umzf), bs      , comzf ),
        cbind( t(amzf), t(mo)  , t(gpcm), t(avsif), t(umzf) , t(fo)  , t(comzf), bs    )), name="corMZF"),
    mxAlgebra(rbind(
        cbind( fv     , dzf    , pd     , ftsi    , mf      , dzfh   , mo      , adzf  ),
        cbind( t(dzf) , fv     , pd     , ftsi    , dzfh    , mf     , adzf    , mo    ),
        cbind( t(pd)  , t(pd)  , sp     , fmo     , pdil    , pdil   , gpcm    , gpcm  ),
        cbind( t(ftsi), t(ftsi), t(fmo) , bs      , hftsi   , hftsi  , avsif   , avsif ),
        cbind( t(mf)  , t(dzfh), t(pdil), t(hftsi), mv      , hdzfh  , fo      , udzf  ), 
        cbind( t(dzfh), t(mf)  , t(pdil), t(hftsi), t(hdzfh), mv     , udzf    , fo    ),
        cbind( t(mo)  , t(adzf), t(gpcm), t(avsif), t(fo)   , t(udzf), bs      , codzf ),
        cbind( t(adzf), t(mo)  , t(gpcm), t(avsif), t(udzf) , t(fo)  , t(codzf), bs    )), name="corDZF"),
    mxAlgebra(rbind(
        cbind( mv      , dzmf    , ps     , mtsi    , mf        , dzmfh    , fo       , udzmf  ),
        cbind( t(dzmf) , fv      , pd     , ftsi    , dzfmw     , mf       , adzmf    , mo     ),
        cbind( t(ps)   , t(pd)   , sp     , fmo     , psil      , pdil     , gpcp     , gpcm   ),
        cbind( t(mtsi) , t(ftsi) , t(fmo) , bs      , wmtsi     , hftsi    , avsim    , avsif  ),
        cbind( t(mf)   , t(dzfmw), t(psil), t(wmtsi), fv        , wdzmfh   , mo       , aidzmf ), 
        cbind( t(dzmfh), t(mf)   , t(pdil), t(hftsi), t(wdzmfh) , mv       , uidzmf   , fo     ),
        cbind( t(fo)   , t(adzmf), t(gpcp), t(avsim), t(mo)     , t(uidzmf), bs       , codzmf ),
        cbind( t(udzmf), t(mo)   , t(gpcm), t(avsif), t(aidzmf) , t(fo)    , t(codzmf), bs     )), name="corDZO"),

    mxMatrix("Full", 1, 16, T, 0, labMeans, -1, 1, name="means16"),
    mxMatrix("Full", 1, 2, T, 1.5, labSDevs, 0.0001, 10, name="sdev2"),

    mxMatrix("Full", 1, nvr, T, 0, labMeansM, -1, 1, name="meanMZM"),
    mxMatrix("Full", 1, nvr, T, 0, labMeansM, -1, 1, name="meanDZM"), 
    mxMatrix("Full", 1, nvr, T, 0, labMeansF, -1, 1, name="meanMZF"), 
    mxMatrix("Full", 1, nvr, T, 0, labMeansF, -1, 1, name="meanDZF"), 
    mxMatrix("Full", 1, nvr, T, 0, labMeansMF, -1, 1, name="meanDZO"), 
    mxMatrix("Diag", nvr, nvr, T, 1.5, labSDevM, 0.0001, 10, name="sdM"),
    mxMatrix("Diag", nvr, nvr, T, 1.5, labSDevF, 0.0001, 10, name="sdF"),
    mxMatrix("Diag", nvr, nvr, T, 1.5, labSDevMF, 0.0001, 10, name="sdMF"),
    
    mxAlgebra( sqrt(sdM)%&%corMZM, name="covMZM"),
    mxAlgebra( sqrt(sdM)%&%corDZM, name="covDZM"),
    mxAlgebra( sqrt(sdF)%&%corMZF, name="covMZF"),
    mxAlgebra( sqrt(sdF)%&%corDZF, name="covDZF"),
    mxAlgebra( sqrt(sdMF)%&%corDZO, name="covDZO"),
    
    mxAlgebra(rbind(
        cbind(mzm,mzf,x,x),                         #mztwins
        cbind(dzm,dzf,dzmf,x),                      #dztwins
        cbind(sim,sif,simf,x),                      #siblings
        cbind(fs,md,fd,ms),                         #parents 
        cbind(gfps,gmpd,gfpd,gmps),                 #grandparents
        cbind(gfms,gmmd,gfmd,gmms),                 #grandparents
        cbind(unemzm,x,unimzm,x),                   #avuncular mz twin
        cbind(x,animzf,x,anemzf),                   #avuncular mz twin
        cbind(unedzm,anidzmf,unidzm,anedzmf),       #avuncular dz twin
        cbind(unedzmf,anidzf,unidzmf,anedzf),       #avuncular dz twin
        cbind(unesim,anisifm,unisim,anesifm),       #avuncular sibs
        cbind(unesimf,anisif,unisimf,anesif),       #avuncular sibs
        cbind(commzm,cofmzm,comfmzm,x),             #cousins mzm
        cbind(commzf,cofmzf,comfmzf,x),             #cousins mzf
        cbind(comdzm,cofdzm,comfdzm,x),             #cousins dzm
        cbind(comdzf,cofdzf,comfdzf,x),             #cousins dzf
        cbind(comdzmf,cofdzmf,comfdzmf,cofmdzmf),   #cousins dzmf
        cbind(x,x,mf,x),                            #spouses
        cbind(x,x,mzmw,mzfh),                       #spouse w mz cotwin
        cbind(dzmfh,dzfmw,dzmw,dzfh),               #spouse w dz cotwin
        cbind(simfh,sifmw,simw,sifh),               #spouses w siblings
        cbind(fdh,msw,fsw,mdh),                     #spouses w parents
        cbind(hmzfh,wmzmw,x,x),                     #spouse w spouse of mz cotwin
        cbind(hdzfh,wdzmw,wdzmfh,x),                #spouse w spouse of dz cotwin
        cbind(x,wunimzm,x,wunemzm),                 #spouse w avuncular mz twin
        cbind(hanemzf,x,hanimzf,x),                 #spouse w avuncular mz twin
        cbind(hanedzmf,wunidzm,hanidzmf,wunedzm),   #spouse w avuncular dz twin
        cbind(hanedzf,wunidzmf,hanidzf,wunedzmf)),  #spouse w avuncular dz twin
        name="expCor"),
#    c('mztw','dztw','si','pc','gpp','gpm','amzp','amzm','adzp','adzm','asip','asim',
#        'comzm','comzf','codzm','codzf','codzmf','sp','smztw','sdztw','ssi','spa',
#        'sspmz','sspdz','samzp','samzm','sadzp','sadzm'),
#    c('Mm','Ff','Mf','Fm')'

#    mxMatrix("Full", 1, 5, F, .0001, name="cnstrPos"), 
#    mxAlgebra(cbind(min(Re(eigen(mzmef, only.values=T)$values)),
# 		   			min(Re(eigen(dzmef, only.values=T)$values)),
# 		   			min(Re(eigen(mzfef, only.values=T)$values)),
# 		   			min(Re(eigen(dzfef, only.values=T)$values)),
# 		   			min(Re(eigen(dzmfef, only.values=T)$values))), name="minCor"),
#    mxConstraint('minCor',">",'cnstrPos', name="constr"),
    mxModel("MZM",
	    mxData( observed=dataMZM, type="raw" ),
	    mxFitFunctionML(),mxExpectationNormal( covariance="et.covMZM", means="et.meanMZM", dimnames=selVars ) ),
    mxModel("DZM",
	    mxData( observed=dataDZM, type="raw" ),
	    mxFitFunctionML(),mxExpectationNormal( covariance="et.covDZM", means="et.meanDZM", dimnames=selVars ) ),
    mxModel("MZF",
	    mxData( observed=dataMZF, type="raw" ),
	    mxFitFunctionML(),mxExpectationNormal( covariance="et.covMZF", means="et.meanMZF", dimnames=selVars ) ),
    mxModel("DZF",
	    mxData( observed=dataDZF, type="raw" ),
	    mxFitFunctionML(),mxExpectationNormal( covariance="et.covDZF", means="et.meanDZF", dimnames=selVars ) ),
    mxModel("DZO",
	    mxData( observed=dataDZO, type="raw" ),
	    mxFitFunctionML(),mxExpectationNormal( covariance="et.covDZO", means="et.meanDZO", dimnames=selVars ) ),
	mxAlgebra( expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective, 		name="m2LL" ),
	mxFitFunctionAlgebra("m2LL")
    )
    
# Run Model and Generate Output
# -----------------------------------------------------------------------
etc88Model <- mxOption(etc88Model,"Standard Errors", "No")
etc88Model <- mxOption(etc88Model,"Calculate Hessian", "No")
etc88Fit <-  mxRun(etc88Model)
summary(etc88Fit)


#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(etc88Fit$output$Minus,161045.203,0.05)

omxCheckCloseEnough(etc88Fit$output$estimate[["SDm"]],1.3359,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["SDf"]],1.3381,0.001)

omxCheckCloseEnough(etc88Fit$output$estimate[["MT1"]],0.0112,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FT1"]],-0.0098,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MT2"]],-0.0029,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FT2"]],-0.0070,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FA"]],0.0216,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MO"]],-0.0051,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MS"]],-0.0070,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FS"]],-0.0087,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MS1"]],-0.0311,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FS1"]],0.0147,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MS2"]],0.0389,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FS2"]],0.0301,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MC1"]],-0.0441,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FC1"]],0.0064,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["MC2"]],-0.0062,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["FC2"]],-0.0194,0.001)

omxCheckCloseEnough(etc88Fit$output$estimate[["r1"]],0.7690,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r2"]],0.5700,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r3"]],0.7802,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r4"]],0.5711,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r5"]],0.5656,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r6"]],-0.0116,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r7"]],-0.0079,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r8"]],-0.0240,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r9"]],-0.0166,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r10"]],0.0395,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r11"]],-0.0185,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r12"]],-0.0765,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r13"]],0.0802,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r14"]],-0.0636,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r15"]],-0.0106,0.02)
omxCheckCloseEnough(etc88Fit$output$estimate[["r16"]],-0.0116,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r17"]],0.5479,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r18"]],0.5426,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r19"]],0.5506,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r20"]],0.0100,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r21"]],-0.0335,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r22"]],0.0046,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r23"]],0.0123,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r24"]],0.0002,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r25"]],0.4808,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r26"]],0.4749,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r27"]],0.4772,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r28"]],0.4772,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r29"]],-0.0171,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r30"]],-0.0379,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r31"]],0.0170,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r32"]],0.0056,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r33"]],0.2290,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r34"]],0.2380,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r35"]],0.2249,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r36"]],0.2152,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r37"]],0.2906,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r38"]],0.2304,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r39"]],0.2927,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r40"]],0.1820,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r41"]],0.4355,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r42"]],0.4039,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r43"]],0.3830,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r44"]],0.3901,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r45"]],0.2425,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r46"]],0.3052,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r47"]],0.2942,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r48"]],0.2703,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r49"]],0.3330,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r50"]],0.3094,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r51"]],0.3039,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r52"]],0.2393,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r53"]],0.3343,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r54"]],0.3363,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r55"]],0.3186,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r56"]],0.2931,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r57"]],0.2704,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r58"]],0.2629,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r59"]],0.2592,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r60"]],0.3215,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r61"]],-0.0454,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r62"]],-0.0147,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r63"]],0.0213,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r64"]],-0.0828,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r65"]],0.0350,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r66"]],-0.0466,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r67"]],-0.0637,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r68"]],0.0594,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r69"]],0.0177,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r70"]],-0.1113,0.02)
omxCheckCloseEnough(etc88Fit$output$estimate[["r71"]],-0.0262,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r72"]],0.0509,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r73"]],0.3495,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r74"]],0.1425,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r75"]],0.1506,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r76"]],0.2102,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r77"]],0.1577,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r78"]],0.2218,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r79"]],0.2716,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r80"]],0.1927,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r81"]],0.3341,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r82"]],0.2338,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r83"]],0.3007,0.01)
omxCheckCloseEnough(etc88Fit$output$estimate[["r84"]],0.4257,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r85"]],0.0313,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r86"]],0.1701,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r87"]],0.0906,0.001)
omxCheckCloseEnough(etc88Fit$output$estimate[["r88"]],0.0521,0.01)
