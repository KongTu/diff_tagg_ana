/*
**************************************************
Original Author: Bill Li (SBU)
Modification by: Kong Tu (BNL)
Purpose        : Look at diffractive physics
Date.          : Aug 18 2022

Descriptions:




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

#include <g4main/PHG4Reco.h>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <g4eval/SvtxEvalStack.h>

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

  m_hepmctree = new TTree("hepmctree", "A tree with hepmc truth particles");
  m_hepmctree->Branch("m_partid1", &m_partid1, "m_partid1/I");
  m_hepmctree->Branch("m_partid2", &m_partid2, "m_partid2/I");
  m_hepmctree->Branch("m_x1", &m_x1, "m_x1/D");
  m_hepmctree->Branch("m_x2", &m_x2, "m_x2/D");
  m_hepmctree->Branch("m_mpi", &m_mpi, "m_mpi/I");
  m_hepmctree->Branch("m_process_id", &m_process_id, "m_process_id/I");
  m_hepmctree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
  m_hepmctree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
  m_hepmctree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
  m_hepmctree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
  m_hepmctree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
  m_hepmctree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
  m_hepmctree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
  m_hepmctree->Branch("m_numparticlesinevent", &m_numparticlesinevent, "m_numparticlesinevent/I");
  m_hepmctree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

  m_partid1 = -99;
  m_partid2 = -99;
  m_x1 = -99;
  m_x2 = -99;
  m_mpi = -99;
  m_process_id = -99;
  m_truthenergy = -99;
  m_trutheta = -99;
  m_truthphi = -99;
  m_truthp = -99;
  m_truthpx = -99;
  m_truthpy = -99;
  m_truthpz = -99;
  m_truthpt = -99;
  m_numparticlesinevent = -99;
  m_truthpid = -99;

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

  // if( static_event_counter == 0) {
  
  // 	static_event_counter++;
  // }
  // cout << " END initialization" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::process_event(PHCompositeNode *topNode)
{

  std::cout << "diff_tagg_ana::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;

  SvtxEvalStack *_svtxEvalStack;

  _svtxEvalStack = new SvtxEvalStack(topNode);
  _svtxEvalStack->set_verbosity(Verbosity());

 /// Getting the Truth information
  getHEPMCTruth(topNode);

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
  m_hepmctree->Write();
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

void diff_tagg_ana::getHEPMCTruth(PHCompositeNode *topNode)
{
  /// Get the node from the node tree
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return;
  }

  /// Could have some print statements for debugging with verbosity
  if (Verbosity() > 1)
  {
    cout << "Getting HEPMC truth particles " << endl;
  }

  /// You can iterate over the number of events in a hepmc event
  /// for pile up events where you have multiple hard scatterings per bunch crossing
  for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
       eventIter != hepmceventmap->end();
       ++eventIter)
  {
    /// Get the event
    PHHepMCGenEvent *hepmcevent = eventIter->second;

    if (hepmcevent)
    {
      /// Get the event characteristics, inherited from HepMC classes
      HepMC::GenEvent *truthevent = hepmcevent->getEvent();
      if (!truthevent)
      {
        cout << PHWHERE
             << "no evt pointer under phhepmvgeneventmap found "
             << endl;
        return;
      }

      /// Get the parton info
      HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

      /// Get the parton info as determined from HEPMC
      m_partid1 = pdfinfo->id1();
      m_partid2 = pdfinfo->id2();
      m_x1 = pdfinfo->x1();
      m_x2 = pdfinfo->x2();

      /// Are there multiple partonic intercations in a p+p event
      m_mpi = truthevent->mpi();

      /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
      m_process_id = truthevent->signal_process_id();

      if (Verbosity() > 2)
      {
        cout << " Iterating over an event" << endl;
      }
      /// Loop over all the truth particles and get their information
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
           iter != truthevent->particles_end();
           ++iter)
      {
        /// Get each pythia particle characteristics
        m_truthenergy = (*iter)->momentum().e();
        m_truthpid = (*iter)->pdg_id();

        m_trutheta = (*iter)->momentum().pseudoRapidity();
        m_truthphi = (*iter)->momentum().phi();
        m_truthpx = (*iter)->momentum().px();
        m_truthpy = (*iter)->momentum().py();
        m_truthpz = (*iter)->momentum().pz();
        m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

        /// Fill the truth tree
        m_hepmctree->Fill();
        m_numparticlesinevent++;
      }
    }
  }
}


