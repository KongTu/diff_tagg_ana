//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in diff_tagg_ana.h.
//
// diff_tagg_ana(const std::string &name = "diff_tagg_ana")
// everything is keyed to diff_tagg_ana, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// diff_tagg_ana::~diff_tagg_ana()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int diff_tagg_ana::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int diff_tagg_ana::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int diff_tagg_ana::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int diff_tagg_ana::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int diff_tagg_ana::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int diff_tagg_ana::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int diff_tagg_ana::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void diff_tagg_ana::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

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

//____________________________________________________________________________..
//diff_tagg_ana::diff_tagg_ana(const std::string &name):
// SubsysReco(name)
//{
//  std::cout << "diff_tagg_ana::diff_tagg_ana(const std::string &name) Calling ctor" << std::endl;
//}


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

  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TNtuple("hitntup", "G4Hits", "x0:y0:z0:x1:y1:z1:edep");

  std::cout << "diff_tagg_ana::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  //**************
  // Truth

  gDirectory->mkdir("Truth");
  gDirectory->cd("Truth");

  m_truthtree = new TTree("truthg4tree", "A tree with truth g4 particles");
  m_truthtree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
  m_truthtree->Branch("m_truthp", &m_truthp, "m_truthp/D");
  m_truthtree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
  m_truthtree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
  m_truthtree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
  m_truthtree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
  m_truthtree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
  m_truthtree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
  m_truthtree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

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

  gDirectory->cd("/");

  m_mpi = -99;
  m_process_id = -99;
  m_partid1 = -99;
  m_partid2 = -99;
  m_x1 = -99;
  m_x2 = -99;
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
//  std::cout << "diff_tagg_ana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
//

  is_electron = "false";
  is_positron = "false";

  is_Jpsi = "false";

  if( static_event_counter == 0) {
  
  	encloseure_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_hFarFwdBeamLineEnclosure_0");

        encloseure_nodeparams->Print();

  	if (encloseure_nodeparams)
  	{
  	   Enclosure_params.FillFrom(encloseure_nodeparams, 0);
  	} else {
  	   cerr << "There is a issue finding the detector paramter node!" << endl;
  	}


  	beamlinemagnet_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_BEAMLINEMAGNET");

//	beamlinemagnet_nodeparams->print();

  	if (encloseure_nodeparams)
  	{
  	   BeamLineMagnet_params.FillFrom(beamlinemagnet_nodeparams, 0);
  	} else {
  	   cerr << "There is a issue finding the detector paramter node!" << endl;
  	}

//	exit(0);

    zdc_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_ZDCsurrogate");  	
    rp_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth");
    rp_nodeparams->print();
    rp2_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth2");  	
    b0_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_b0Truth_0");

  	static_event_counter++;

       /// Determining which IP design
	if (zdc_nodeparams) {
  	   if (rp2_nodeparams) {
		IP_design = "IP8";
	   } else {
		IP_design = "IP6";
           }
	} else {
	   IP_design = "UNKNOWN";
        }

  }

//  exit(0);

  cout << " END initialization" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int diff_tagg_ana::process_event(PHCompositeNode *topNode)
{

  is_electron = "false";
  is_positron = "false";
  is_Jpsi = "false";


  std::cout << "diff_tagg_ana::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;

  SvtxEvalStack *_svtxEvalStack;

  _svtxEvalStack = new SvtxEvalStack(topNode);
  _svtxEvalStack->set_verbosity(Verbosity());

 /// Getting the Truth information
  process_PHG4Truth(topNode);
  // getHEPMCTruth(topNode);

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
  // m_hepmctree->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

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

int diff_tagg_ana::process_PHG4Truth(PHCompositeNode* topNode) {

 PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  /// Get the primary particle range
  ///PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::Range range = truthinfo->GetParticleRange();

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);
    m_truthphi = atan(m_truthpy / m_truthpx);
    float m_trutheta = atanh(m_truthpz / m_truthenergy);
    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;
    float m_truthpid = truth->get_pid();

	m_truthpid = m_truthpid;

	if (m_truthpid == 11) {
	   is_electron = "true";
	}

	if (m_truthpid == -11) {
	   is_positron = "true";
	}

   // cout << "truth: " << m_truthpid << "  " << m_truthpx << "  " << m_truthpy 
   //      << "  " << m_truthpz << endl;

    /// Fill the g4 truth tree
   m_truthtree->Fill();
  }


  return Fun4AllReturnCodes::EVENT_OK;

}

// void diff_tagg_ana::getHEPMCTruth(PHCompositeNode *topNode)
// {
//   /// Get the node from the node tree
//   PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

//   /// If the node was not properly put on the tree, return
//   if (!hepmceventmap)
//   {
//     cout << PHWHERE
//          << "HEPMC event map node is missing, can't collected HEPMC truth particles"
//          << endl;
//     return;
//   }

//   /// Could have some print statements for debugging with verbosity
//   if (Verbosity() > 1)
//   {
//     cout << "Getting HEPMC truth particles " << endl;
//   }

//   /// You can iterate over the number of events in a hepmc event
//   /// for pile up events where you have multiple hard scatterings per bunch crossing
//   for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
//        eventIter != hepmceventmap->end();
//        ++eventIter)
//   {
//     /// Get the event
//     PHHepMCGenEvent *hepmcevent = eventIter->second;

//     if (hepmcevent)
//     {
//       /// Get the event characteristics, inherited from HepMC classes
//       HepMC::GenEvent *truthevent = hepmcevent->getEvent();
//       if (!truthevent)
//       {
//         cout << PHWHERE
//              << "no evt pointer under phhepmvgeneventmap found "
//              << endl;
//         return;
//       }

//       /// Get the parton info
//       HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

//       /// Get the parton info as determined from HEPMC
//       m_partid1 = pdfinfo->id1();
//       m_partid2 = pdfinfo->id2();
//       m_x1 = pdfinfo->x1();
//       m_x2 = pdfinfo->x2();

//       /// Are there multiple partonic intercations in a p+p event
//       m_mpi = truthevent->mpi();

//       /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
//       m_process_id = truthevent->signal_process_id();

//       if (Verbosity() > 2)
//       {
//         cout << " Iterating over an event" << endl;
//       }
//       /// Loop over all the truth particles and get their information
//       for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
//            iter != truthevent->particles_end();
//            ++iter)
//       {
//         /// Get each pythia particle characteristics
//         m_truthenergy = (*iter)->momentum().e();
//         m_truthpid = (*iter)->pdg_id();

//         m_trutheta = (*iter)->momentum().pseudoRapidity();
//         m_truthphi = (*iter)->momentum().phi();
//         m_truthpx = (*iter)->momentum().px();
//         m_truthpy = (*iter)->momentum().py();
//         m_truthpz = (*iter)->momentum().pz();
//         m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

//         /// Fill the truth tree
//         m_hepmctree->Fill();
//         m_numparticlesinevent++;
//       }
//     }
//   }
// }




