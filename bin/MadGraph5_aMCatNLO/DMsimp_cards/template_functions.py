def extramodels(mPhi, mChi):
    return """DMsimp_s_spin1.tar.gz"""

def customize_cards(mPhi, mChi):
    return """# change mass parameters
set param_card mass 55 {0:d}
set param_card mass 18 {1:d}
set param_card decay 55 auto""".format(int(mPhi), int(mChi))

def run_card(mPhi, mChi):
    return """#***********************************************************************
#                        MadGraph5_aMC@NLO                             *
#                                                                      *
#                      run_card.dat aMC@NLO                            *
#                                                                      *
#  This file is used to set the parameters of the run.                 *
#                                                                      *
#  Some notation/conventions:                                          *
#                                                                      *
#   Lines starting with a hash (#) are info or comments                *
#                                                                      *
#   mind the format:   value    = variable     ! comment               *
#                                                                      *
#   Some of the values of variables can be list. These can either be   *
#   comma or space separated.                                          *
#***********************************************************************
#
#*******************                                                 
# Running parameters
#*******************                                                 
#
#***********************************************************************
# Tag name for the run (one word)                                      *
#***********************************************************************
  tag_1	= run_tag ! name of the run 
#***********************************************************************
# Number of LHE events (and their normalization) and the required      *
# (relative) accuracy on the Xsec.                                     *
# These values are ignored for fixed order runs                        *
#***********************************************************************
  10	= nevents ! Number of unweighted events requested 
  -1.0	= req_acc ! Required accuracy (-1=auto determined from nevents)
  10	= nevt_job ! Max number of events per job in event generation. 
                 !  (-1= no split).
#***********************************************************************
# Normalize the weights of LHE events such that they sum or average to *
# the total cross section                                              *
#***********************************************************************
  average	= event_norm ! average or sum
#***********************************************************************
# Number of points per itegration channel (ignored for aMC@NLO runs)   *
#***********************************************************************
  0.01	= req_acc_fo ! Required accuracy (-1=ignored, and use the 
 	                   ! number of points and iter. below)
# These numbers are ignored except if req_acc_FO is equal to -1
  5000	= npoints_fo_grid ! number of points to setup grids
  4	= niters_fo_grid ! number of iter. to setup grids
  10000	= npoints_fo ! number of points to compute Xsec
  6	= niters_fo ! number of iter. to compute Xsec
#***********************************************************************
# Random number seed                                                   *
#***********************************************************************
  123456	= iseed ! rnd seed (0=assigned automatically=default))
#***********************************************************************
# Collider type and energy                                             *
#***********************************************************************
  1	= lpp1 ! beam 1 type (0 = no PDF)
  1	= lpp2 ! beam 2 type (0 = no PDF)
  6500.0	= ebeam1 ! beam 1 energy in GeV
  6500.0	= ebeam2 ! beam 2 energy in GeV
#***********************************************************************
# PDF choice: this automatically fixes also alpha_s(MZ) and its evol.  *
#***********************************************************************
  lhapdf	= pdlabel ! PDF set
  260000	= lhaid ! If pdlabel=lhapdf, this is the lhapdf number. Only 
              ! numbers for central PDF sets are allowed. Can be a list; 
              ! PDF sets beyond the first are included via reweighting.
#***********************************************************************
# Include the NLO Monte Carlo subtr. terms for the following parton    *
# shower (HERWIG6 | HERWIGPP | PYTHIA6Q | PYTHIA6PT | PYTHIA8)         *
# WARNING: PYTHIA6PT works only for processes without FSR!!!!          *
#***********************************************************************
  PYTHIA8	= parton_shower 
  1.0	= shower_scale_factor ! multiply default shower starting
                                  ! scale by this factor
#***********************************************************************
# Renormalization and factorization scales                             *
# (Default functional form for the non-fixed scales is the sum of      *
# the transverse masses divided by two of all final state particles    * 
# and partons. This can be changed in SubProcesses/set_scales.f or via *
# dynamical_scale_choice option)                                       *
#***********************************************************************
  False	= fixed_ren_scale ! if .true. use fixed ren scale
  False	= fixed_fac_scale ! if .true. use fixed fac scale
  91.118	= mur_ref_fixed ! fixed ren reference scale 
  91.118	= muf_ref_fixed ! fixed fact reference scale
  -1	= dynamical_scale_choice ! Choose one (or more) of the predefined
           ! dynamical choices. Can be a list; scale choices beyond the
           ! first are included via reweighting
  1.0	= mur_over_ref ! ratio of current muR over reference muR
  1.0	= muf_over_ref ! ratio of current muF over reference muF
#*********************************************************************** 
# Reweight variables for scale dependence and PDF uncertainty          *
#***********************************************************************
  1.0, 2.0, 0.5	= rw_rscale ! muR factors to be included by reweighting
  1.0, 2.0, 0.5	= rw_fscale ! muF factors to be included by reweighting
  True	= reweight_scale ! Reweight to get scale variation using the 
            ! rw_rscale and rw_fscale factors. Should be a list of 
            ! booleans of equal length to dynamical_scale_choice to
            ! specify for which choice to include scale dependence.
  False	= reweight_pdf ! Reweight to get PDF uncertainty. Should be a
            ! list booleans of equal length to lhaid to specify for
            !  which PDF set to include the uncertainties.
#***********************************************************************
# Store reweight information in the LHE file for off-line model-       *
# parameter reweighting at NLO+PS accuracy                             *
#***********************************************************************
  True	= store_rwgt_info ! Store info for reweighting in LHE file
#***********************************************************************
# ickkw parameter:                                                     *
#   0: No merging                                                      *
#   3: FxFx Merging - WARNING! Applies merging only at the hard-event  *
#      level. After showering an MLM-type merging should be applied as *
#      well. See http://amcatnlo.cern.ch/FxFx_merging.htm for details. *
#   4: UNLOPS merging (with pythia8 only). No interface from within    *
#      MG5_aMC available, but available in Pythia8.                    *
#  -1: NNLL+NLO jet-veto computation. See arxiv:1412.8408 [hep-ph].    *
#***********************************************************************
  3	= ickkw 
#***********************************************************************
#
#***********************************************************************
# BW cutoff (M+/-bwcutoff*Gamma). Determines which resonances are      *
# written in the LHE event file                                        *
#***********************************************************************
  15.0	= bwcutoff 
#***********************************************************************
# Cuts on the jets. Jet clustering is performed by FastJet.            *
#  - When matching to a parton shower, these generation cuts should be *
#    considerably softer than the analysis cuts.                       *
#  - More specific cuts can be specified in SubProcesses/cuts.f        *
#***********************************************************************
  1.0	= jetalgo ! FastJet jet algorithm (1=kT, 0=C/A, -1=anti-kT)
  0.7	= jetradius ! The radius parameter for the jet algorithm
  10.0	= ptj ! Min jet transverse momentum
  -1.0	= etaj ! Max jet abs(pseudo-rap) (a value .lt.0 means no cut)
#***********************************************************************
# Cuts on the charged leptons (e+, e-, mu+, mu-, tau+ and tau-)        *
# More specific cuts can be specified in SubProcesses/cuts.f           *
#***********************************************************************
  0.0	= ptl ! Min lepton transverse momentum
  -1.0	= etal ! Max lepton abs(pseudo-rap) (a value .lt.0 means no cut)
  0.0	= drll ! Min distance between opposite sign lepton pairs
  0.0	= drll_sf ! Min distance between opp. sign same-flavor lepton pairs
  0.0	= mll ! Min inv. mass of all opposite sign lepton pairs
  30.0	= mll_sf ! Min inv. mass of all opp. sign same-flavor lepton pairs
#***********************************************************************
# Photon-isolation cuts, according to hep-ph/9801442. When ptgmin=0,   *
# all the other parameters are ignored.                                *
# More specific cuts can be specified in SubProcesses/cuts.f           *
#***********************************************************************
  20.0	= ptgmin ! Min photon transverse momentum
  -1.0	= etagamma ! Max photon abs(pseudo-rap)
  0.4	= r0gamma ! Radius of isolation code
  1.0	= xn ! n parameter of eq.(3.4) in hep-ph/9801442
  1.0	= epsgamma ! epsilon_gamma parameter of eq.(3.4) in hep-ph/9801442
  True	= isoem ! isolate photons from EM energy (photons and leptons)
#***********************************************************************
# For aMCfast+APPLGRID use in PDF fitting (http://amcfast.hepforge.org)*
#***********************************************************************
  0	= iappl ! aMCfast switch (0=OFF, 1=prepare grids, 2=fill grids)
#***********************************************************************"""

def proc_card(mPhi, mChi):
    return """
#************************************************************
#*                     MadGraph5_aMC@NLO                    *
#*                                                          *
#*                *                       *                 *
#*                  *        * *        *                   *
#*                    * * * * 5 * * * *                     *
#*                  *        * *        *                   *
#*                *                       *                 *
#*                                                          *
#*                                                          *
#*         VERSION 2.4.3                 2016-08-01         *
#*                                                          *
#*    The MadGraph5_aMC@NLO Development Team - Find us at   *
#*    https://server06.fynu.ucl.ac.be/projects/madgraph     *
#*                                                          *
#************************************************************
#*                                                          *
#*               Command File for MadGraph5_aMC@NLO         *
#*                                                          *
#*     run as ./bin/mg5_aMC  filename                       *
#*                                                          *
#************************************************************
set group_subprocesses Auto
set ignore_six_quark_processes False
set loop_optimized_output True
set loop_color_flows False
set gauge unitary
set complex_mass_scheme False
set max_npoint_for_channel 0
import model DMsimp_s_spin1
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define p = 21 2 4 1 3 -2 -4 -1 -3 5 -5 # pass to 5 flavors
define j = p
generate       p p  > xd xd~ j   / h w+ w- z a [QCD]
add process    p p  > xd xd~ j j / h w+ w- z a [QCD]
output A12JToChiChi_MA-{0:d}_MChi-{1:d}
""".format(int(mPhi), int(mChi))

import sys
import os
import glob
import shutil

if __name__ == "__main__":
    args = sys.argv[1:]

    # 1 - mPhi, 2 - mChi
    assert len(args) == 2
    mphi = float(args[0])
    mchi = float(args[1])

    model = "A12JToChiChi_MA-{0:d}_MChi-{1:d}".format(int(mphi), int(mchi))
    if not os.path.exists(model):
        os.makedirs(model)

    # Create the parameter dependent files
    filename = os.path.join(model, "{}_customizecards.dat".format(model))
    with open(os.path.join(filename), 'w') as f:
        f.write(customize_cards(mphi, mchi))

    filename = os.path.join(model, "{}_extramodels.dat".format(model))
    with open(os.path.join(filename), 'w') as f:
        f.write(extramodels(mphi, mchi))

    filename = os.path.join(model, "{}_proc_card.dat".format(model))
    with open(os.path.join(filename), 'w') as f:
        f.write(proc_card(mphi, mchi))

    filename = os.path.join(model, "{}_run_card.dat".format(model))
    with open(os.path.join(filename), 'w') as f:
        f.write(run_card(mphi, mchi))