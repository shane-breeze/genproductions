def param_card(mPhi, mChi):
    return """
######################################################################
## PARAM_CARD AUTOMATICALY GENERATED BY THE UFO  #####################
######################################################################

###################################
## INFORMATION FOR SMINPUTS
###################################
Block SMINPUTS 
    1 1.279000e+02 # aEWM1 
    2 1.166370e-05 # Gf 
    3 1.184000e-01 # aS 

###################################
## INFORMATION FOR MASS
###################################
Block MASS 
    6 1.720000e+02 # MT 
   15 1.777000e+00 # MTA 
   23 9.118760e+01 # MZ 
   25 1.250000e+02 # MH 
   51 1.000000e+01 # MXc 
   55 {mPhi:.6e}
  5000001 1.000000e+01 # MXr 
  9100012 {mChi:.6e} # MXd 
##  Not dependent paramater.
## Those values should be edited following analytical the 
## analytical expression. Some generator could simply ignore 
## those values and use the analytical expression
  22 0.000000 # a : 0.0 
  24 79.824360 # W+ : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  21 0.000000 # g : 0.0 
  9000001 0.000000 # ghA : 0.0 
  9000003 79.824360 # ghWp : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  9000004 79.824360 # ghWm : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  82 0.000000 # ghG : 0.0 
  12 0.000000 # ve : 0.0 
  14 0.000000 # vm : 0.0 
  16 0.000000 # vt : 0.0 
  11 0.000000 # e- : 0.0 
  13 0.000000 # mu- : 0.0 
  2 0.000000 # u : 0.0 
  4 0.000000 # c : 0.0 
  1 0.000000 # d : 0.0 
  3 0.000000 # s : 0.0 
  5 0.000000 # b : 0.0 
  251 79.824360 # G+ : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  9000002 91.187600 # ghZ : MZ 
  250 91.187600 # G0 : MZ 
  18 {mChi:.6f} # Xd : MXd 

###################################
## INFORMATION FOR DECAY
###################################
DECAY   6 1.508336e+00 
DECAY  23 2.495200e+00 
DECAY  24 2.085000e+00 
DECAY  25 4.070000e-03 
DECAY  55 1.000000e+01 
##  Not dependent paramater.
## Those values should be edited following analytical the 
## analytical expression. Some generator could simply ignore 
## those values and use the analytical expression
DECAY  22 0.000000 # a : 0.0 
DECAY  21 0.000000 # g : 0.0 
DECAY  9000001 0.000000 # ghA : 0.0 
DECAY  82 0.000000 # ghG : 0.0 
DECAY  12 0.000000 # ve : 0.0 
DECAY  14 0.000000 # vm : 0.0 
DECAY  16 0.000000 # vt : 0.0 
DECAY  11 0.000000 # e- : 0.0 
DECAY  13 0.000000 # mu- : 0.0 
DECAY  15 0.000000 # ta- : 0.0 
DECAY  2 0.000000 # u : 0.0 
DECAY  4 0.000000 # c : 0.0 
DECAY  1 0.000000 # d : 0.0 
DECAY  3 0.000000 # s : 0.0 
DECAY  5 0.000000 # b : 0.0 
DECAY  5000001 0.000000 # Xr : 0.0 
DECAY  51 0.000000 # Xc : 0.0 
DECAY  18 0.000000 # Xd : 0.0 
DECAY  9000002 2.495200 # ghZ : WZ 
DECAY  9000003 2.085000 # ghWp : WW 
DECAY  9000004 2.085000 # ghWm : WW 
DECAY  250 2.495200 # G0 : WZ 
DECAY  251 2.085000 # G+ : WW 

###################################
## INFORMATION FOR CKMBLOCK
###################################
Block CKMBLOCK 
    1 2.277360e-01 # cabi 

###################################
## INFORMATION FOR DMINPUTS
###################################
Block DMINPUTS 
    1 0.000000e+00 # gVXc 
    2 1.000000e-99 # gVXd 
    3 1.001000e+00 # gAXd 
    4 1.000000e-99 # gVd11 
    5 1.000000e-99 # gVu11 
    6 1.000000e-99 # gVd22 
    7 1.000000e-99 # gVu22 
    8 1.000000e-99 # gVd33 
    9 1.000000e-99 # gVu33 
   10 2.500000e-01 # gAd11 
   11 2.500000e-01 # gAu11 
   12 2.500000e-01 # gAd22 
   13 2.500000e-01 # gAu22 
   14 2.500000e-01 # gAd33 
   15 2.500000e-01 # gAu33 
   16 0.000000e+00 # gVh 

###################################
## INFORMATION FOR LOOP
###################################
Block LOOP 
    1 9.118800e+01 # MU_R 

###################################
## INFORMATION FOR YUKAWA
###################################
Block YUKAWA 
    6 1.720000e+02 # ymt 
   15 1.777000e+00 # ymtau 
#===========================================================
# QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE)
#===========================================================

Block QNUMBERS 9000001  # ghA 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000002  # ghZ 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000003  # ghWp 
        1 3  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000004  # ghWm 
        1 -3  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 82  # ghG 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 8  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 250  # G0 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 251  # G+ 
        1 3  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 5000001  # Xr 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 51  # Xc 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 18  # Xd 
        1 0  # 3 times electric charge
        2 2  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 55  # Y1 
        1 0  # 3 times electric charge
        2 3  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
""".format(mPhi = mPhi, mChi = mChi)

def parameters(mPhi, mChi):
    return """
# This file was automatically created by FeynRules 2.3.7
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Mon 24 Aug 2015 13:37:17



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# This is a default parameter object representing the renormalization scale (MU_R).
MU_R = Parameter(name = 'MU_R',
                 nature = 'external',
                 type = 'real',
                 value = 91.188,
                 texname = '\\text{\\mu_r}',
                 lhablock = 'LOOP',
                 lhacode = [1])

# User-defined parameters.
cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.227736,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

gVXc = Parameter(name = 'gVXc',
                 nature = 'external',
                 type = 'real',
                 value = 0.,
                 texname = 'g_{\\text{VXc}}',
                 lhablock = 'DMINPUTS',
                 lhacode = [ 1 ])

gVXd = Parameter(name = 'gVXd',
                 nature = 'external',
                 type = 'real',
                 value = 1e-99,
                 texname = 'g_{\\text{VXd}}',
                 lhablock = 'DMINPUTS',
                 lhacode = [ 2 ])

gAXd = Parameter(name = 'gAXd',
                 nature = 'external',
                 type = 'real',
                 value = 1.001,
                 texname = 'g_{\\text{AXd}}',
                 lhablock = 'DMINPUTS',
                 lhacode = [ 3 ])

gVd11 = Parameter(name = 'gVd11',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vd11}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 4 ])

gVu11 = Parameter(name = 'gVu11',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vu11}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 5 ])

gVd22 = Parameter(name = 'gVd22',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vd22}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 6 ])

gVu22 = Parameter(name = 'gVu22',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vu22}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 7 ])

gVd33 = Parameter(name = 'gVd33',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vd33}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 8 ])

gVu33 = Parameter(name = 'gVu33',
                  nature = 'external',
                  type = 'real',
                  value = 1e-99,
                  texname = 'g_{\\text{Vu33}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 9 ])

gAd11 = Parameter(name = 'gAd11',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Ad11}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 10 ])

gAu11 = Parameter(name = 'gAu11',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Au11}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 11 ])

gAd22 = Parameter(name = 'gAd22',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Ad22}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 12 ])

gAu22 = Parameter(name = 'gAu22',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Au22}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 13 ])

gAd33 = Parameter(name = 'gAd33',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Ad33}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 14 ])

gAu33 = Parameter(name = 'gAu33',
                  nature = 'external',
                  type = 'real',
                  value = 0.25,
                  texname = 'g_{\\text{Au33}}',
                  lhablock = 'DMINPUTS',
                  lhacode = [ 15 ])

gVh = Parameter(name = 'gVh',
                nature = 'external',
                type = 'real',
                value = 0.,
                texname = 'g_{\\text{Vh}}',
                lhablock = 'DMINPUTS',
                lhacode = [ 16 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

MXr = Parameter(name = 'MXr',
                nature = 'external',
                type = 'real',
                value = 10.,
                texname = '\\text{MXr}',
                lhablock = 'MASS',
                lhacode = [ 5000001 ])

MXc = Parameter(name = 'MXc',
                nature = 'external',
                type = 'real',
                value = 10.,
                texname = '\\text{MXc}',
                lhablock = 'MASS',
                lhacode = [ 51 ])

MXd = Parameter(name = 'MXd',
                nature = 'external',
                type = 'real',
                value = """ + """{mChi:.1f},""".format(mChi=mChi) + """
                texname = '\\text{MXd}',
                lhablock = 'MASS',
                lhacode = [ 9100012 ])

MY1 = Parameter(name = 'MY1',
                nature = 'external',
                type = 'real',
                value = """ + """{mPhi:.1f},""".format(mPhi=mPhi) + """
                texname = '\\text{MY1}',
                lhablock = 'MASS',
                lhacode = [ 55 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00407,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

WY1 = Parameter(name = 'WY1',
                nature = 'external',
                type = 'real',
                value = 10.,
                texname = '\\text{WY1}',
                lhablock = 'DECAY',
                lhacode = [ 55 ])

CKM1x1 = Parameter(name = 'CKM1x1',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM1x1}')

CKM1x2 = Parameter(name = 'CKM1x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.sin(cabi)',
                   texname = '\\text{CKM1x2}')

CKM2x1 = Parameter(name = 'CKM2x1',
                   nature = 'internal',
                   type = 'complex',
                   value = '-cmath.sin(cabi)',
                   texname = '\\text{CKM2x1}')

CKM2x2 = Parameter(name = 'CKM2x2',
                   nature = 'internal',
                   type = 'complex',
                   value = 'cmath.cos(cabi)',
                   texname = '\\text{CKM2x2}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2a33}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3a33}')

MFU = Parameter(name = 'MFU',
               nature = 'internal',
               type = 'real',
               value = '0.002550',
               texname = '\\text{MFU}')

MFC = Parameter(name = 'MFC',
               nature = 'internal',
               type = 'real',
               value = '1.27',
               texname = '\\text{MFC}')

MFD = Parameter(name = 'MFD',
               nature = 'internal',
               type = 'real',
               value = '0.00504',
               texname = '\\text{MFD}')

MFS = Parameter(name = 'MFS',
               nature = 'internal',
               type = 'real',
               value = '0.101',
               texname = '\\text{MFS}')

MFB = Parameter(name = 'MFB',
               nature = 'internal',
               type = 'real',
               value = '4.7',
               texname = '\\text{MFB}')
"""

def restrict_test(mPhi, mChi):
    return """
######################################################################
## PARAM_CARD AUTOMATICALY GENERATED BY THE UFO  #####################
######################################################################

###################################
## INFORMATION FOR SMINPUTS
###################################
Block SMINPUTS 
    1 1.279000e+02 # aEWM1 
    2 1.166370e-05 # Gf 
    3 1.184000e-01 # aS 

###################################
## INFORMATION FOR MASS
###################################
Block MASS 
    6 1.720000e+02 # MT 
   15 1.777000e+00 # MTA 
   23 9.118760e+01 # MZ 
   25 1.250000e+02 # MH 
   51 1.000000e+01 # MXc 
   55 {mPhi:.6e} # MY1 
  5000001 1.000000e+01 # MXr 
  9100012 {mChi:.6e} # MXd 
##  Not dependent paramater.
## Those values should be edited following analytical the 
## analytical expression. Some generator could simply ignore 
## those values and use the analytical expression
  22 0.000000 # a : 0.0 
  24 79.824360 # W+ : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  21 0.000000 # g : 0.0 
  9000001 0.000000 # ghA : 0.0 
  9000003 79.824360 # ghWp : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  9000004 79.824360 # ghWm : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  82 0.000000 # ghG : 0.0 
  12 0.000000 # ve : 0.0 
  14 0.000000 # vm : 0.0 
  16 0.000000 # vt : 0.0 
  11 0.000000 # e- : 0.0 
  13 0.000000 # mu- : 0.0 
  2 0.000000 # u : 0.0 
  4 0.000000 # c : 0.0 
  1 0.000000 # d : 0.0 
  3 0.000000 # s : 0.0 
  5 0.000000 # b : 0.0 
  251 79.824360 # G+ : cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2)))) 
  9000002 91.187600 # ghZ : MZ 
  250 91.187600 # G0 : MZ 
  18 {mChi:.6f} # Xd : MXd 

###################################
## INFORMATION FOR DECAY
###################################
DECAY   6 1.508336e+00 
DECAY  23 2.495200e+00 
DECAY  24 2.085000e+00 
DECAY  25 4.070000e-03 
DECAY  55 1.000000e+01 
##  Not dependent paramater.
## Those values should be edited following analytical the 
## analytical expression. Some generator could simply ignore 
## those values and use the analytical expression
DECAY  22 0.000000 # a : 0.0 
DECAY  21 0.000000 # g : 0.0 
DECAY  9000001 0.000000 # ghA : 0.0 
DECAY  82 0.000000 # ghG : 0.0 
DECAY  12 0.000000 # ve : 0.0 
DECAY  14 0.000000 # vm : 0.0 
DECAY  16 0.000000 # vt : 0.0 
DECAY  11 0.000000 # e- : 0.0 
DECAY  13 0.000000 # mu- : 0.0 
DECAY  15 0.000000 # ta- : 0.0 
DECAY  2 0.000000 # u : 0.0 
DECAY  4 0.000000 # c : 0.0 
DECAY  1 0.000000 # d : 0.0 
DECAY  3 0.000000 # s : 0.0 
DECAY  5 0.000000 # b : 0.0 
DECAY  5000001 0.000000 # Xr : 0.0 
DECAY  51 0.000000 # Xc : 0.0 
DECAY  18 0.000000 # Xd : 0.0 
DECAY  9000002 2.495200 # ghZ : WZ 
DECAY  9000003 2.085000 # ghWp : WW 
DECAY  9000004 2.085000 # ghWm : WW 
DECAY  250 2.495200 # G0 : WZ 
DECAY  251 2.085000 # G+ : WW 

###################################
## INFORMATION FOR CKMBLOCK
###################################
Block CKMBLOCK 
    1 2.277360e-01 # cabi 

###################################
## INFORMATION FOR DMINPUTS
###################################
Block DMINPUTS 
    1 0.000000e+00 # gVXc 
    2 1.000000e-99 # gVXd 
    3 1.001000e+00 # gAXd 
    4 1.000000e-99 # gVd11 
    5 1.000000e-99 # gVu11 
    6 1.000000e-99 # gVd22 
    7 1.000000e-99 # gVu22 
    8 1.000000e-99 # gVd33 
    9 1.000000e-99 # gVu33 
   10 2.500000e-01 # gAd11 
   11 2.500000e-01 # gAu11 
   12 2.500000e-01 # gAd22 
   13 2.500000e-01 # gAu22 
   14 2.500000e-01 # gAd33 
   15 2.500000e-01 # gAu33 
   16 0.000000e+00 # gVh 

###################################
## INFORMATION FOR LOOP
###################################
Block LOOP 
    1 9.118800e+01 # MU_R 

###################################
## INFORMATION FOR YUKAWA
###################################
Block YUKAWA 
    6 1.720000e+02 # ymt 
   15 1.777000e+00 # ymtau 
#===========================================================
# QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE)
#===========================================================

Block QNUMBERS 9000001  # ghA 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000002  # ghZ 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000003  # ghWp 
        1 3  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 9000004  # ghWm 
        1 -3  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 82  # ghG 
        1 0  # 3 times electric charge
        2 -1  # number of spin states (2S+1)
        3 8  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 250  # G0 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 251  # G+ 
        1 3  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 5000001  # Xr 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 51  # Xc 
        1 0  # 3 times electric charge
        2 1  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 18  # Xd 
        1 0  # 3 times electric charge
        2 2  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 1  # Particle/Antiparticle distinction (0=own anti)
Block QNUMBERS 55  # Y1 
        1 0  # 3 times electric charge
        2 3  # number of spin states (2S+1)
        3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 0  # Particle/Antiparticle distinction (0=own anti)
""".format(mPhi = mPhi, mChi = mChi)


import sys
import os
import glob
import shutil
import tarfile

if __name__ == "__main__":
    args = sys.argv[1:]

    # 1 - mPhi, 2 - mChi
    assert len(args) == 2
    mphi = float(args[0])
    mchi = float(args[1])

    pathname = "DMsimp_s_spin1_{:d}_{:d}_801".format(int(mphi), int(mchi))
    if not os.path.exists(pathname):
        os.makedirs(pathname)

    # Create the parameter dependent files
    filename = os.path.join(pathname, "param_card.dat")
    with open(os.path.join(filename), 'w') as f:
        f.write(param_card(mphi, mchi))

    filename = os.path.join(pathname, "parameters.py")
    with open(os.path.join(filename), 'w') as f:
        f.write(parameters(mphi, mchi))

    filename = os.path.join(pathname, "restrict_test.dat")
    with open(os.path.join(filename), 'w') as f:
        f.write(parameters(mphi, mchi))

    # Copy other files
    sourcename = "DMsimp_s_spin1_template"
    for path in glob.glob(sourcename+"/*"):
        shutil.copyfile(path, pathname+"/"+os.path.basename(path))

    # Tarring time
    with tarfile.open(pathname+".tar.gz", mode='w:gz') as tar:
        tar.add(pathname)
