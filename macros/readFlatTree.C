#include "pleaseIncludeMe.h"
void readFlatTree(TString vm_name="jpsi"){

	TChain* eventtree = new TChain("eventtree");
	TChain* truthtree = new TChain("truthg4tree");
	TChain* tracktree = new TChain("tracktree");
	eventtree->Add("../flatTree/"+vm_name+".root");
	truthtree->Add("../flatTree/"+vm_name+".root");
	tracktree->Add("../flatTree/"+vm_name+".root");

	if(eventtree){
		eventtree->SetBranchAddress("m_Q2_truth",&m_Q2_truth);
		eventtree->SetBranchAddress("m_x_truth",&m_x_truth);
		eventtree->SetBranchAddress("m_y_truth",&m_y_truth);
		eventtree->SetBranchAddress("m_e_eta_truth",&m_e_eta_truth);
		eventtree->SetBranchAddress("m_e_phi_truth",&m_e_phi_truth);
		eventtree->SetBranchAddress("m_e_pt_truth",&m_e_pt_truth);
	}
	if(truthtree){
		truthtree->SetBranchAddress("m_nMCtracks",&m_nMCtracks);
		truthtree->SetBranchAddress("m_truthenergy",m_truthenergy);
		truthtree->SetBranchAddress("m_truthpt",m_truthpt);
		truthtree->SetBranchAddress("m_trutheta",m_trutheta);
		truthtree->SetBranchAddress("m_truthphi",m_truthphi);
		truthtree->SetBranchAddress("m_truthpid",m_truthpid);
	}
	if(tracktree){
		tracktree->SetBranchAddress("m_nRECtracks",&m_nRECtracks);
		tracktree->SetBranchAddress("m_nRECtracksMatch",&m_nRECtracksMatch);
		tracktree->SetBranchAddress("m_nRECtracksPID",&m_nRECtracksPID);
		tracktree->SetBranchAddress("m_tr_px",m_tr_px);
		tracktree->SetBranchAddress("m_tr_py",m_tr_py);
		tracktree->SetBranchAddress("m_tr_pz",m_tr_pz);
		tracktree->SetBranchAddress("m_tr_p",m_tr_p);
		tracktree->SetBranchAddress("m_tr_pt",m_tr_pt);
		tracktree->SetBranchAddress("m_tr_eta",m_tr_eta);
		tracktree->SetBranchAddress("m_tr_phi",m_tr_phi);
		tracktree->SetBranchAddress("m_charge",m_charge);
		tracktree->SetBranchAddress("m_truth_is_primary",m_truth_is_primary);
		tracktree->SetBranchAddress("m_tr_px",m_tr_px);
		tracktree->SetBranchAddress("m_truthtrackpt",m_truthtrackpt);
		tracktree->SetBranchAddress("m_truthtracketa",m_truthtracketa);
		tracktree->SetBranchAddress("m_truthtrackphi",m_truthtrackphi);
	}

	//Friend everyone.
	if(eventtree){
		eventtree->AddFriend(truthtree);
		eventtree->AddFriend(tracktree);
	}

	double vm_mass=0.;
	double dau_mass=0.;
	int dau_pid=0;
	double vm_mass_width=0.;
	if(vm_name=="jpsi"){
		vm_mass=3.09;
		dau_mass=MASS_ELECTRON;	
		dau_pid=11;
		vm_mass_width=0.007;
	}
	else if(vm_name=="phi"){
		vm_mass=1.02;
		dau_mass=MASS_KAON;
		dau_pid=321;
		vm_mass_width=0.02;
	}

	//histogram
	TH1D* h_t_MC=new TH1D("h_t_MC","-t (GeV^{2})",46,t_binning);
	TH1D* h_mass_MC=new TH1D("h_mass_MC","mass (GeV/c^{2})",100,0.4,4.0);
	TH1D* h_t_REC=new TH1D("h_t_REC","-t (GeV^{2})",46,t_binning);
	TH1D* h_mass_REC=new TH1D("h_mass_REC","mass (GeV/c^{2})",100,0.4,4.0);
	TH1D* h_pt_MC=new TH1D("h_pt_MC","p_{T} (GeV/c)",100,0.,4.0);
	TH1D* h_pt_match_RECO=new TH1D("h_pt_match_RECO","p_{T} (GeV/c)",100,0.,4.0);
	TH1D* h_pt_res[14];
	TH2D* h_pt_res_2D[14];
	for(int ieta=0;ieta<14;ieta++){
		h_pt_res[ieta]=new TH1D(Form("h_pt_res_%d",ieta),"p_{T} (GeV/c)",100,-0.3,0.3);
		h_pt_res_2D[ieta]=new TH2D(Form("h_pt_res_2D_%d",ieta),"p_{T} (GeV/c); res",100,0,10,100,-0.3,0.3);
	}
	//analysis
	int nEvents=eventtree->GetEntries();
	
	for(int i=0;i<nEvents;i++){
		eventtree->GetEntry(i);
		if(m_Q2_truth<1.||m_Q2_truth>10.) continue;
		//truth
		TLorentzVector scatElec,dau1,dau2,vmcand;
		vector< TLorentzVector> scatElecCand;
		//MC track loop
		for(int imc=0;imc<m_nMCtracks;imc++){
			TVector3 p1;
			p1.SetPtEtaPhi(m_truthpt[imc],m_trutheta[imc],m_truthphi[imc]);
			if( fabs(m_trutheta[imc])<0.5 ) h_pt_MC->Fill( m_truthpt[imc] );
			//find scat cand
			if(m_truthpid[imc]==+11){
				scatElec.SetPtEtaPhiM(m_truthpt[imc],m_trutheta[imc],m_truthphi[imc],MASS_ELECTRON);
				scatElecCand.push_back( scatElec );		
			}
			//find dau1
			if( m_truthpid[imc]==dau_pid){
				dau1.SetPtEtaPhiM(m_truthpt[imc],m_trutheta[imc],m_truthphi[imc],dau_mass);
			}
			//find dau2
			if(m_truthpid[imc]==-dau_pid){
				dau2.SetPtEtaPhiM(m_truthpt[imc],m_trutheta[imc],m_truthphi[imc],dau_mass);
			}
			for(int imatch=0;imatch<m_nRECtracksMatch;imatch++){
				TVector3 p2;
				p2.SetPtEtaPhi(m_tr_pt[imatch],m_tr_eta[imatch],m_tr_phi[imatch]);
				if( matchGen(p1,p2) 
					&& fabs(m_trutheta[imc])<0.5){
					h_pt_match_RECO->Fill( m_tr_pt[imatch] );
				}
			}//for efficiency.
		}

		//find scat
		scatElec=scatElecCand[0];
		double mineta=scatElec.Eta();
		double maxenergy=scatElec.E();
		for(unsigned jelec=1;jelec<scatElecCand.size();jelec++){
			if( scatElecCand[jelec].Eta()<mineta 
				&& scatElecCand[jelec].E()>maxenergy){
				scatElec=scatElecCand[jelec];
				mineta=scatElecCand[jelec].Eta();
				maxenergy=scatElecCand[jelec].E();
			}
		}

		//VM candidate.
		if(vm_name=="jpsi"){
			for(unsigned jelec=0;jelec<scatElecCand.size();jelec++){
				if(fabs(scatElecCand[jelec].Eta()-scatElec.Eta())<1e-5) continue;
				dau1 = scatElecCand[jelec];
			}
		}
		vmcand=dau1+dau2;
		if( fabs(vmcand.M()-vm_mass)<vm_mass_width 
				&& fabs(vmcand.Rapidity())<4. ) {//pick some real vmcand.
			h_t_MC->Fill( giveme_t_L(scatElec,vmcand) );
		}
		h_mass_MC->Fill( vmcand.M() );
		scatElecCand.clear();

		//REC track loop
		TLorentzVector scatElec_REC;
		vector<TLorentzVector> dau1_REC, dau2_REC;
		for(int ireco=0;ireco<m_nRECtracks;ireco++){
			//tracking resolution.
			for(int ieta=0;ieta<14;ieta++){
				if(m_truthtracketa[ireco]>eta_binning[ieta] 
					&& m_truthtracketa[ireco]<eta_binning[ieta+1])
				{
					double mc=m_truthtrackpt[ireco];
					double rec=m_tr_pt[ireco];
					h_pt_res[ieta]->Fill( (mc-rec)/mc );
					h_pt_res_2D[ieta]->Fill(mc, (mc-rec)/mc);
				}
			}

			//finding REC scat electron by God.
			TVector3 p;
			p.SetPtEtaPhi(m_truthtrackpt[ireco],m_truthtracketa[ireco],m_truthtrackphi[ireco]);
			if(fabs(scatElec.Vect().Pt()-p.Pt())<1e-4 
				&& fabs(scatElec.Vect().Eta()-p.Eta())<1e-4){//match GEN
				scatElec_REC.SetPtEtaPhiM(m_tr_pt[ireco],m_tr_eta[ireco],m_tr_phi[ireco], MASS_ELECTRON);
			}
			//charged particles		
			TLorentzVector d1tmp,d2tmp;
			if(m_charge[ireco]==+1){
				d1tmp.SetPtEtaPhiM(m_tr_pt[ireco],m_tr_eta[ireco],m_tr_phi[ireco], dau_mass);
				dau1_REC.push_back(d1tmp);
			}
			if(m_charge[ireco]==-1){
				d2tmp.SetPtEtaPhiM(m_tr_pt[ireco],m_tr_eta[ireco],m_tr_phi[ireco], dau_mass);
				dau2_REC.push_back(d2tmp);
			}
			
		}//end REC track loop

		//loop over all positive and negative cands.
		for(int id1=0;id1<dau1_REC.size();id1++){
			for(int id2=0;id2<dau2_REC.size();id2++){
				TLorentzVector vmcand_REC=dau1_REC[id1]+dau2_REC[id2];
				if( fabs(vmcand_REC.M()-vm_mass)<vm_mass_width 
						&& fabs(vmcand_REC.Rapidity())<4. ) {//pick some real vmcand_REC.
					h_t_REC->Fill( giveme_t_L(scatElec_REC,vmcand_REC) );
				}
				h_mass_REC->Fill( vmcand_REC.M() );
			}
		}
		dau1_REC.clear();
		dau2_REC.clear();

	}//event loop
	
	TFile* output = new TFile("output_"+vm_name+".root","RECREATE");
	//
	gDirectory->mkdir("Physics");
	gDirectory->cd("Physics");
	h_t_MC->Write();
	h_mass_MC->Write();
	h_t_REC->Write();
	h_mass_REC->Write();
	gDirectory->cd("/");
	gDirectory->mkdir("Tracking");
	gDirectory->cd("Tracking");
	h_pt_MC->Write();
	h_pt_match_RECO->Write();
	TH1D* h_Eff=(TH1D*)h_pt_match_RECO->Clone("h_Eff");
	h_Eff->Divide(h_pt_MC);
	h_Eff->Write();
	for(int ieta=0;ieta<14;ieta++){
		h_pt_res[ieta]->Write();
	}
	for(int ieta=0;ieta<14;ieta++){
		h_pt_res_2D[ieta]->Write();
	}
	output->Close();



}