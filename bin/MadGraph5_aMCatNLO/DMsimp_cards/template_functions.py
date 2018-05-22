def extramodels(mPhi, mChi):
    return """DMsimp_s_spin1_{:d}_{:d}_801.tar.gz""".format(int(mPhi), int(mChi))

def customize_cards(mPhi, mChi):
    return """################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
#
# This File contains some configuration variable for MadGraph/MadEvent
#
# Line starting by #! are comment and should remain commented
# Line starting with # should be uncommented if you want to modify the default
#    value.
# Current value for all options can seen by typing "display options"
#    after either ./bin/mg5_aMC or ./bin/madevent 
#
# You can place this files in ~/.mg5/mg5_configuration.txt if you have more than
#    one version of MG5. 
#
################################################################################
#! Prefered Fortran Compiler
#! If None: try to find g77 or gfortran on the system
#!
# fortran_compiler = None
# f2py_compiler = None 


#! Prefered C++ Compiler
#! If None: try to find g++ or clang on the system
#!
# cpp_compiler = None

#! Prefered Text Editor
#!  Default: use the shell default Editor
#!           or try to find one available on the system
#!  Be careful: Only shell based editor are allowed
# text_editor = None

#! Prefered WebBrower
#! If None: try to find one available on the system
# web_browser = None

#! Prefered PS viewer
#!  If None: try to find one available on the system
# eps_viewer = None

#! Time allowed to answer question (if no answer takes default value)
#!  0: No time limit
# timeout = 60

#! Pythia8 path.
#!  Defines the path to the main pythia8 directory (i.e. that containing
#!  the pythia8 configure script) .
#!  If using a relative path, that starts from the mg5 directory
# pythia8_path =

#! Herwig++ paths
#!  specify here the paths also to HepMC ant ThePEG
#!  define the path to the herwig++, thepeg and hepmc directories.
#!  paths can be absolute or relative from mg5 directory
# hwpp_path = 
# thepeg_path = 
# hepmc_path = 

#! Control when MG5 checks if he is up-to-date.
#! Enter the number of day between two check (0 means never)
#! A question is always asked before any update
# auto_update = 7

################################################################################
#  INFO FOR MADEVENT / aMC@NLO 
################################################################################
# If this file is in a MADEVENT Template. 'main directory' is the directory
# containing the SubProcesses directory. Otherwise this is the MadGraph5_aMC@NLO main
# directory (containing the directories madgraph and Template)

#! Allow/Forbid the automatic opening of the web browser  (on the status page)
#!  when launching MadEvent [True/False]
# automatic_html_opening = True
#! allow notification of finished job in the notification center (Mac Only)
# notification_center = True


#! Default Running mode 
#!  0: single machine/ 1: cluster / 2: multicore
# run_mode = 2

#! Cluster Type [pbs|sge|condor|lsf|ge|slurm|htcaas|htcaas2] Use for cluster run only
#!  And cluster queue (or partition for slurm)
#!  And size of the cluster (some part of the code can adapt splitting accordingly)
# cluster_type = condor
# cluster_queue = madgraph
# cluster_size = 150 

#! Path to a node directory to avoid direct writing on the central disk
#!  Note that condor clusters avoid direct writing by default (therefore this
#!  options does not affect condor clusters)
# cluster_temp_path = None

#! path to a node directory where local file can be found (typically pdf)
#! to avoid to send them to the node (if cluster_temp_path is on True or condor)
# cluster_local_path =  None # example: /cvmfs/cp3.uclouvain.be/madgraph/

#! Cluster waiting time for status update 
#!  First number is when the number of waiting job is higher than the number 
#!  of running one (time in second). The second number is in the second case.
# cluster_status_update = 600 30

#! How to deal with failed submission (can occurs on cluster mode)
#!  0: crash, -1: print error, hangs the program up to manual instructions, N(>0) retry up to N times.
# cluster_nb_retry = 1

#! How much time to wait for the output file before resubmission/crash (filesystem can be very slow)
# cluster_retry_wait = 300

#! Nb_core to use (None = all) This is use only for multicore run
#!  This correspond also to the number core used for code compilation for cluster mode
# nb_core = None

#! Pythia-PGS Package
#!  relative path start from main directory
# pythia-pgs_path = ./pythia-pgs

#! Delphes Package
#!  relative path start from main directory
# delphes_path = ./Delphes

#! MadAnalysis Package [For Drawing output]
#!  relative path start from main directory
# madanalysis_path = ./MadAnalysis

#! ExRootAnalysis Package
#!  relative path start from main directory
# exrootanalysis_path = ./ExRootAnalysis

#! TOPDRAWER PATH
#!  Path to the directory containing td executables
#!  relative path start from main directory
# td_path = ./td

#! lhapdf-config
#!  If None: try to find one available on the system
lhapdf = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf6/6.1.5-cms/bin/lhapdf-config #  

#! fastjet-config
#!  If None: try to find one available on the system
# fastjet = fastjet-config

#! MCatNLO-utilities 
#!  relative path starting from main directory
# MCatNLO-utilities_path = ./MCatNLO-utilities

#! Set what OLP to use for the loop ME generation
# OLP = MadLoop

#! Set the PJFRy++ directory containing pjfry's library
#! if auto: try to find it automatically on the system (default)
#! if '' or None: disabling pjfry
#! if pjfry=/PATH/TO/pjfry/lib: use that specific installation path for PJFry++
# pjfry = auto

#! Set the Golem95 directory containing golem's library
#! It only supports version higher than 1.3.0
#! if auto: try to find it automatically on the system (default)
#! if '' or None: disabling Golem95
#! if golem=/PATH/TO/golem/lib: use that speficif installation path for Golem95
# golem = auto

#! Set the samurai directory containing samurai's library
#! It only supports version higher than 2.0.0
#! if auto: try to find it automatically on the system (default)
#! if '' or None: disabling samurai
#! if samurai=/PATH/TO/samurai/lib: use that specific installation path for samurai
# samurai = None

#! Set the samurai directory containing ninja's library
#! if auto: try to find it automatically on the system (default)
#! if '' or None: disabling ninja 
#! if ninja=/PATH/TO/ninja/lib: use that specific installation path for ninja 
# ninja = ./HEPTools/lib

#! Set how MadLoop dependencies (such as CutTools) should be handled
#!  > external : ML5 places a link to the MG5_aMC-wide libraries
#!  > internal : ML5 copies all dependencies in the output so that it is independent
#!  > environment_paths : ML5 searches for the dependencies in your environment path
# output_dependencies = external

#! SysCalc PATH
#! Path to the directory containing syscalc executables
#! relative path start from main directory
# syscalc_path = ./SysCalc

#! Absolute paths to config scripts in the bin directories for APPLgrid and aMCFast.
# applgrid = applgrid-config
# amcfast = amcfast-config



# MG5 MAIN DIRECTORY
mg5_path = /afs/cern.ch/work/p/pharris/public/bacon/gen/MadGraph5_aMCatNLO_Grid/Axial_MonoJR2_NLO_Mphi-MED_Mchi-XMASS_gSM-0p25_gDM-1p0_13TeV-madgraph_MG5_aMC_v2_4_3
mg5_path = ../mgbasedir
cluster_temp_path = None
lhapdf = /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/share/LHAPDF/../../bin/lhapdf-config
run_mode = 2
nb_core = 1
lhapdf = /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/share/LHAPDF/../../bin/lhapdf-config
run_mode = 2
nb_core = 1
lhapdf = /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/share/LHAPDF/../../bin/lhapdf-config
run_mode = 2
nb_core = 1"""

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
import model DMsimp_s_spin1_{0:d}_{1:d}_801
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
