#include "RiceStyle.h"
using namespace std;

#define PI 3.1415926
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

Float_t m_Q2_truth;
Float_t m_x_truth;
Float_t m_y_truth;
Float_t m_e_eta_truth;
Float_t m_e_phi_truth;
Float_t m_e_pt_truth;
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

// truth variables
int m_nMCtracks;
double m_truthenergy[200];
double m_trutheta[200];
double m_truthphi[200];
double m_truthpx[200];
double m_truthpy[200];
double m_truthpz[200];
double m_truthpt[200];
double m_truthp[200];
int m_truthpid[200];

double eta_binning[15]={-3.5,-3.0,-2.5,-2.,-1.5,-1,-0.5,0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5};
double t_binning[47]={0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.079,0.083,0.087,0.091,0.095,0.099,0.103,0.107,0.111,0.121,0.131,0.141,0.151,0.161,0.171,0.181};
double giveme_t_L(TLorentzVector eScat, TLorentzVector VM){

  TLorentzVector eIn(0,0,-18,18);
  TLorentzVector pIn(-2.749,0,109.996,110.034); //Crossing angle
  // TLorentzVector pIn(0,0,109.996,110.00);//no Crossing angle
  TLorentzVector eOut=eScat;
  TLorentzVector vmOut=VM;
  TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
  double method_L = -99.;
  TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
  double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
  double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
  double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
  TLorentzVector a_beam_scattered_corr; 
  a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
  method_L = (a_beam_scattered_corr-aInVec).Mag2();

  return -method_L;
}
double giveme_t_A(TLorentzVector eScat, TLorentzVector VM){

  double method_A = -99;
  TVector2 sum_pt(VM.Px()+eScat.Px(), VM.Py()+eScat.Py());
  method_A = sum_pt.Mod2();

  return method_A;
}
bool matchGen(TVector3 p1, TVector3 p2){

  if(p1.DeltaR(p2)<0.05 && fabs(p1.Pt()-p2.Pt())<0.1 ) return true;
  else return false;
}



