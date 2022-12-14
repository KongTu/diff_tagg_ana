// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef DIFF_TAGG_ANA_H
#define DIFF_TAGG_ANA_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Reco.h>

#include <string>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <TFile.h>
#include <TNtuple.h>

#include <cassert>
#include <sstream>
#include <string>
#include <iostream>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

//PID includes
#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

using namespace std;

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;

class diff_tagg_ana : public SubsysReco
{
 public:

//  diff_tagg_ana(const std::string &name = "diff_tagg_ana");
  diff_tagg_ana(const std::string &name = "Diff_Tagg_ana", const std::string &fname = "MyNtuple.root");

  virtual ~diff_tagg_ana();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;
  
  /// Event tree
  void getEvent(PHCompositeNode *topNode);

  /// MC tree
  void getPHG4Truth(PHCompositeNode *topNode);
  
  /// REC tree
  void getTracks(PHCompositeNode *topNode);

  /// ZDC tree
  void getZDC(PHCompositeNode *topNode);

  /// RP tree
  void getRP(PHCompositeNode *topNode);

  /// OMD tree
  void getOMD(PHCompositeNode *topNode);

  /// B0 tree
  void getB0(PHCompositeNode *topNode);

// private:
 protected:

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TTree* m_eventtree;
  TTree *m_truthtree;
  TTree *m_tracktree;
  TNtuple *g4hitntuple;
  TNtuple *g4rphitntuple;
  TNtuple *g4omdhitntuple;
  TNtuple *g4b0hitntuple;

  SvtxEvalStack *m_svtxEvalStack = nullptr;

  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;

  int static_event_counter;

  // Beam parameter
  Float_t e_beam_energy;
  Float_t e_beam_pmag;
  Float_t ion_beam_energy;
  Float_t ion_beam_pmag;
  Float_t crossing_angle;

  Float_t mProt;
  Float_t mElec;

  TLorentzVector eBeam4Vect;
  TLorentzVector pBeam4Vect;
  TLorentzVector virtphoton4VectTruth;
  TLorentzVector e4VectTruth;

  //event level
  Float_t m_Q2_truth;
  Float_t m_x_truth;
  Float_t m_y_truth;
  Float_t m_e_eta_truth;
  Float_t m_e_phi_truth;
  Float_t m_e_pt_truth;

  //hepmc truth variables
  int m_nMCtracks;
  int m_mpi[200];
  int m_process_id[200];
  int m_partid1[200];
  int m_partid2[200];
  double m_x1[200];
  double m_x2[200];
  double m_truthenergy[200];
  double m_trutheta[200];
  double m_truthphi[200];
  double m_truthpx[200];
  double m_truthpy[200];
  double m_truthpz[200];
  double m_truthpt[200];
  double m_truthp[200];
  int m_numparticlesinevent[200];
  int m_truthpid[200];

  /// Track variables
  int m_nRECtracks;
  int m_nRECtracksMatch;
  int m_nRECtracksPID;
  double m_tr_px[200];
  double m_tr_py[200];
  double m_tr_pz[200];
  double m_tr_p[200];
  double m_tr_pt[200];
  double m_tr_phi[200];
  double m_tr_eta[200];
  int m_charge[200];
  double m_chisq[200];
  int m_ndf[200];
  double m_dca[200];
  double m_tr_x[200];
  double m_tr_y[200];
  double m_tr_z[200];
  float m_tr_pion_loglikelihood[200];
  float m_tr_kaon_loglikelihood[200];
  float m_tr_proton_loglikelihood[200];
  int m_truth_is_primary[200];
  double m_truthtrackpx[200];
  double m_truthtrackpy[200];
  double m_truthtrackpz[200];
  double m_truthtrackp[200];
  double m_truthtracke[200];
  double m_truthtrackpt[200];
  double m_truthtrackphi[200];
  double m_truthtracketa[200];
  int m_truthtrackpid[200];

  //FF detectors
  PHParameters Enclosure_params{"PHGEnclosure"};
  PHParameters ZDC_params{"PHG4RP"};
  PHParameters RP_1_params{"PHG4RP"};
  PHParameters RP2_params{"PHG4RP2"};
  PHParameters B0_params{"PHG4B0"};
  PHParameters BeamLineMagnet_params{"PHG4BeamLinMagnet"};

  PdbParameterMapContainer *encloseure_nodeparams; 
  PdbParameterMapContainer *zdc_nodeparams; 
  PdbParameterMapContainer *rp_nodeparams;
  PdbParameterMapContainer *rp2_nodeparams;
  PdbParameterMapContainer *b0_nodeparams;
  PdbParameterMapContainer *beamlinemagnet_nodeparams; 

};

#endif // DIFF_TAGG_ANA_H
