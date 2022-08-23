/*
**************************************************
Original Author: Bill Li (SBU)
Modification by: Kong Tu (BNL)
Purpose        : Look at diffractive physics
Date.          : Aug 18 2022

Descriptions:

* Run on DST trees;
* Save a few important trees for analysis.



**************************************************
*/
#include "diff_tagg_ana.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

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

#include <gsl/gsl_randist.h>

#include <gsl/gsl_rng.h>

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

#include <g4main/PHG4Reco.h>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

using namespace std;

diff_tagg_ana::diff_tagg_ana(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "Diff_Tagg_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
diff_tagg_ana::~diff_tagg_ana()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "diff_tagg_ana::~diff_tagg_ana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int diff_tagg_ana::Init(PHCompositeNode *topNode)
{

  static_event_counter = 0;
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  std::cout << "diff_tagg_ana::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  m_truthtree = new TTree("truthg4tree", "A tree with truth g4 particles");
  m_truthtree->Branch("m_nMCtracks", &m_nMCtracks, "m_nMCtracks/I");
  m_truthtree->Branch("m_truthenergy", m_truthenergy, "m_truthenergy[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthp", m_truthp, "m_truthp[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthpx", m_truthpx, "m_truthpx[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthpy", m_truthpy, "m_truthpy[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthpz", m_truthpz, "m_truthpz[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthpt", m_truthpt, "m_truthpt[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthphi", m_truthphi, "m_truthphi[m_nMCtracks]/D");
  m_truthtree->Branch("m_trutheta", m_trutheta, "m_trutheta[m_nMCtracks]/D");
  m_truthtree->Branch("m_truthpid", m_truthpid, "m_truthpid[m_nMCtracks]/I");

  m_tracktree = new TTree("tracktree", "A tree with svtx tracks");
  m_tracktree->Branch("m_nRECtracks", &m_nRECtracks, "m_nRECtracks/I");
  m_tracktree->Branch("m_nRECtracksMatch", &m_nRECtracksMatch, "m_nRECtracksMatch/I");
  m_tracktree->Branch("m_nRECtracksPID", &m_nRECtracksPID, "m_nRECtracksPID/I");
  m_tracktree->Branch("m_tr_px", m_tr_px, "m_tr_px[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_py", m_tr_py, "m_tr_py[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_pz", m_tr_pz, "m_tr_pz[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_p", m_tr_p, "m_tr_p[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_pt", m_tr_pt, "m_tr_pt[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_phi", m_tr_phi, "m_tr_phi[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_eta", m_tr_eta, "m_tr_eta[m_nRECtracks]/D");
  m_tracktree->Branch("m_charge", m_charge, "m_charge[m_nRECtracks]/I");
  m_tracktree->Branch("m_chisq", m_chisq, "m_chisq[m_nRECtracks]/D");
  m_tracktree->Branch("m_ndf", m_ndf, "m_ndf[m_nRECtracks]/I");
  m_tracktree->Branch("m_dca", m_dca, "m_dca[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_x", m_tr_x, "m_tr_x[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_y", m_tr_y, "m_tr_y[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_z", m_tr_z, "m_tr_z[m_nRECtracks]/D");
  m_tracktree->Branch("m_tr_pion_loglikelihood", m_tr_pion_loglikelihood, "m_tr_pion_loglikelihood[m_nRECtracks]/F");
  m_tracktree->Branch("m_tr_kaon_loglikelihood", m_tr_kaon_loglikelihood, "m_tr_kaon_loglikelihood[m_nRECtracks]/F");
  m_tracktree->Branch("m_tr_proton_loglikelihood", m_tr_proton_loglikelihood, "m_tr_proton_loglikelihood[m_nRECtracks]/F");
  m_tracktree->Branch("m_truth_is_primary", m_truth_is_primary, "m_truth_is_primary[m_nRECtracks]/I");
  m_tracktree->Branch("m_truthtrackpx", m_truthtrackpx, "m_truthtrackpx[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackpy", m_truthtrackpy, "m_truthtrackpy[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackpz", m_truthtrackpz, "m_truthtrackpz[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackp", m_truthtrackp, "m_truthtrackp[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtracke", m_truthtracke, "m_truthtracke[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackpt", m_truthtrackpt, "m_truthtrackpt[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackphi", m_truthtrackphi, "m_truthtrackphi[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtracketa", m_truthtracketa, "m_truthtracketa[m_nRECtracks]/D");
  m_tracktree->Branch("m_truthtrackpid", m_truthtrackpid, "m_truthtrackpid[m_nRECtracks]/I");

  ///**********************************/
  // Parameter definition

  mElec = 0.000510998950;
  mProt = 0.93827208816;

  // Define beam 4 vectors
  e_beam_energy = 18;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 275;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.025; 

  //Double_t Pi = TMath::ACos(-1);
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);

  /**********************************/

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::InitRun(PHCompositeNode *topNode)
{

  cout << " END initialization" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::process_event(PHCompositeNode *topNode)
{

  m_svtxEvalStack = new SvtxEvalStack(topNode);
  m_svtxEvalStack->set_verbosity(Verbosity());
  // Getting the MC information
  getPHG4Truth(topNode);
  // Getting the RECO information
  getTracks(topNode);

  event_itt++; 
 
  if(event_itt%100 == 0)
     std::cout << "Event Processing Counter: " << event_itt << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::ResetEvent(PHCompositeNode *topNode)
{
  //  std::cout << "diff_tagg_ana::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  //
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::EndRun(const int runnumber)
{
  std::cout << "diff_tagg_ana::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::End(PHCompositeNode *topNode)
{
  std::cout << "diff_tagg_ana::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  outfile->cd();
  m_truthtree->Write();
  m_tracktree->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::Reset(PHCompositeNode *topNode)
{
 std::cout << "diff_tagg_ana::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void diff_tagg_ana::Print(const std::string &what) const
{
  
  std::cout << "diff_tagg_ana::Print(const std::string &what) const Printing info for " << what << std::endl;
}

/**
 * This function collects the truth PHG4 stable particles that
 * are produced from PYTHIA, or some other event generator. These
 * are the stable particles, e.g. there are not any (for example)
 * partons here.
 */
void diff_tagg_ana::getPHG4Truth(PHCompositeNode *topNode)
{
  /// G4 truth particle node
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return;
  }

  /// Get the primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

  int itruth=0;
  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx[itruth] = truth->get_px();
    m_truthpy[itruth] = truth->get_py();
    m_truthpz[itruth] = truth->get_pz();
    m_truthp[itruth] = sqrt(m_truthpx[itruth] * m_truthpx[itruth]
                             + m_truthpy[itruth] * m_truthpy[itruth] 
                                  + m_truthpz[itruth] * m_truthpz[itruth]);
    m_truthenergy[itruth] = truth->get_e();

    m_truthpt[itruth] = sqrt(m_truthpx[itruth] * m_truthpx[itruth] + m_truthpy[itruth] * m_truthpy[itruth]);

    m_truthphi[itruth] = atan2(m_truthpy[itruth], m_truthpx[itruth]);

    m_trutheta[itruth] = atanh(m_truthpz[itruth] / m_truthenergy[itruth]);
    /// Check for nans
    if (m_trutheta[itruth] != m_trutheta[itruth])
      m_trutheta[itruth] = -99;
    m_truthpid[itruth] = truth->get_pid();

    itruth++;
  }
  m_nMCtracks=itruth;
  /// Fill the g4 truth tree
  m_truthtree->Fill();

}

/**
 * This method gets the tracks as reconstructed from the tracker. It also
 * compares the reconstructed track to its truth track.
 */
void diff_tagg_ana::getTracks(PHCompositeNode *topNode)
{
  //  Tracks node
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
  EICPIDParticleContainer *pidcontainer = findNode::getClass<EICPIDParticleContainer>(topNode, "EICPIDParticleMap");

  if (Verbosity() > 1 and pidcontainer == nullptr)
  {
    cout << "EICPIDParticleContainer named EICPIDParticleMap does not exist. Skip saving PID info" << endl;
  }

  if (!trackmap)
  {
    cout << PHWHERE
         << "TrackMap node is missing, can't collect tracks"
         << endl;
    return;
  }

  /// Get the range for primary tracks
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (Verbosity() > 1)
  {
    cout << "Get the tracks" << endl;
  }
  int ireco=0;
  int imatch=0;
  int ipid=0;
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    SvtxTrack *track = iter->second;

    /// Get the reconstructed track info
    m_tr_px[ireco] = track->get_px();
    m_tr_py[ireco] = track->get_py();
    m_tr_pz[ireco] = track->get_pz();
    m_tr_p[ireco] = sqrt(m_tr_px[ireco] * m_tr_px[ireco] 
                          + m_tr_py[ireco] * m_tr_py[ireco] 
                            + m_tr_pz[ireco] * m_tr_pz[ireco]);
    m_tr_pt[ireco] = sqrt(m_tr_px[ireco] * m_tr_px[ireco] 
                          + m_tr_py[ireco] * m_tr_py[ireco]);

    m_tr_phi[ireco] = track->get_phi();
    m_tr_eta[ireco] = track->get_eta();

    m_charge[ireco] = track->get_charge();
    m_chisq[ireco] = track->get_chisq();
    m_ndf[ireco] = track->get_ndf();
    m_dca[ireco] = track->get_dca();
    m_tr_x[ireco] = track->get_x();
    m_tr_y[ireco] = track->get_y();
    m_tr_z[ireco] = track->get_z();
    ireco++;

    /// Ensure that the reco track is a fast sim track
    SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(iter->second);
    if (!temp)
    {
      if (Verbosity() > 0)
        std::cout << "Skipping non fast track sim object..." << std::endl;
      continue;
    }

    /// Get truth track info that matches this reconstructed track
    PHG4Particle *truthtrack = truthinfo->GetParticle(temp->get_truth_track_id());
    if (truthtrack)
    {
      m_truth_is_primary[imatch] = truthinfo->is_primary(truthtrack);

      m_truthtrackpx[imatch] = truthtrack->get_px();
      m_truthtrackpy[imatch] = truthtrack->get_py();
      m_truthtrackpz[imatch] = truthtrack->get_pz();
      m_truthtrackp[imatch] = sqrt(m_truthtrackpx[imatch] * m_truthtrackpx[imatch] 
                                  + m_truthtrackpy[imatch] * m_truthtrackpy[imatch] 
                                  + m_truthtrackpz[imatch] * m_truthtrackpz[imatch]);

      m_truthtracke[imatch] = truthtrack->get_e();

      m_truthtrackpt[imatch] = sqrt(m_truthtrackpx[imatch] * m_truthtrackpx[imatch] 
                                    + m_truthtrackpy[imatch] * m_truthtrackpy[imatch]);
      m_truthtrackphi[imatch] = atan2(m_truthtrackpy[imatch], m_truthtrackpx[imatch]);
      m_truthtracketa[imatch] = atanh(m_truthtrackpz[imatch] / m_truthtrackp[imatch]);
      m_truthtrackpid[imatch] = truthtrack->get_pid();
      imatch++;
    }
    

    // match to PIDparticles
    if (pidcontainer)
    {
      // EICPIDParticle are index the same as the tracks
      const EICPIDParticle *pid_particle =
          pidcontainer->findEICPIDParticle(track->get_id());

      if (pid_particle)
      {
        // top level log likelihood sums.
        // More detailed per-detector information also available at  EICPIDParticle::get_LogLikelyhood(EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector)
        m_tr_pion_loglikelihood[ipid] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::PionCandiate);
        m_tr_kaon_loglikelihood[ipid] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::KaonCandiate);
        m_tr_proton_loglikelihood[ipid] = pid_particle->get_SumLogLikelyhood(EICPIDDefs::ProtonCandiate);
        ipid++;
      }
    }
   
  }
  m_nRECtracks=ireco;
  m_nRECtracksMatch=imatch;
  m_nRECtracksPID=ipid;

  m_tracktree->Fill();
}

