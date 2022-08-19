// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef DIFF_TAGG_ANA_H
#define DIFF_TAGG_ANA_H

#include <fun4all/SubsysReco.h>

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
#include <SvtxTrackMap.h>
#include <SvtxTrack.h>

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class SvtxTrackMap;
class SvtxTrackEval;
class SvtxEvalStack;


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
  
  void getHEPMCTruth(PHCompositeNode *topNode);
  
  void getTracks(PHCompositeNode *topNode);

// private:

 protected:

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TTree *m_hepmctree;
  TTree *m_tracktree;
  TNtuple *g4hitntuple;
  TNtuple *clusterntuple;

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

  Float_t Q2_truth;
  Float_t e_eta_truth;
  Float_t e_Phi_truth;
  Float_t e_E_truth;
  Float_t e_Polar_truth;

  Float_t mProt;
  Float_t mElec;

  TLorentzVector eBeam4Vect;
  TLorentzVector pBeam4Vect;
  TLorentzVector virtphoton4VectTruth;
  TLorentzVector e4VectTruth;

  //hepmc truth variables
  int m_mpi;
  int m_process_id;
  int m_partid1;
  int m_partid2;
  double m_x1;
  double m_x2;
  double m_truthenergy;
  double m_trutheta;
  double m_truthphi;
  double m_truthpx;
  double m_truthpy;
  double m_truthpz;
  double m_truthpt;
  double m_truthp;
  int m_numparticlesinevent;
  int m_truthpid;

  /// Track variables
  double m_tr_px;
  double m_tr_py;
  double m_tr_pz;
  double m_tr_p;
  double m_tr_pt;
  double m_tr_phi;
  double m_tr_eta;
  int m_charge;
  double m_chisq;
  int m_ndf;
  double m_dca;
  double m_tr_x;
  double m_tr_y;
  double m_tr_z;
  int m_truth_is_primary;
  double m_truthtrackpx;
  double m_truthtrackpy;
  double m_truthtrackpz;
  double m_truthtrackp;
  double m_truthtracke;
  double m_truthtrackpt;
  double m_truthtrackphi;
  double m_truthtracketa;
  int m_truthtrackpid;

};

#endif // DIFF_TAGG_ANA_H
