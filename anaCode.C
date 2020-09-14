
#include <TFile.h>
#include "ROOT/RVec.hxx"
#include <TLorentzVector.h>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>
using namespace ROOT;
using rvec_i = const RVec<int> &;
using rvec_f = const RVec<float> &;
using rvec_TLvector = const RVec<ROOT::Math::PtEtaPhiMVector> &;
using TLvector = const ROOT::Math::PtEtaPhiMVector &;

inline void vecSetPtEtaPhiM(ROOT::Math::PtEtaPhiMVector& vec, Float_t pt, Float_t eta, Float_t phi, Float_t m) {
    vec.SetPt(pt);
    vec.SetEta(eta);
    vec.SetPhi(phi);
    vec.SetM(m);
}
inline Float_t deltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
    Float_t deltaRSquare = (eta2-eta1)*(eta2-eta1)+(phi2-phi1)*(phi2-phi1);
    return TMath::Sqrt(deltaRSquare);
}
RVec<Int_t> doSelection(Int_t nThinJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, TLvector Hjet1, TLvector Hjet2) {
    const Float_t dRCut = 1.2;
    ROOT::Math::PtEtaPhiMVector v1, v2;
    RVec<Int_t> vbfj;
    for (Int_t i=0;i<nThinJet;i++) {
        if (Jet_pt[i]<50.) continue;
        if (Jet_eta[i]>5.) continue;
        if (deltaR(Hjet1.Eta(),Hjet1.Phi(),Jet_eta[i],Jet_phi[i])<dRCut) continue;
        if (deltaR(Hjet2.Eta(),Hjet2.Phi(),Jet_eta[i],Jet_phi[i])<dRCut) continue;
        vbfj.push_back(i);
    }
    for (UInt_t j=0;j<vbfj.size();j++) {
        Int_t jj = vbfj[j];
        for (UInt_t i=0;i<j;i++) {
            Int_t ii = vbfj[i];
            if (Jet_eta[ii]*Jet_eta[jj]>=0.) continue;
            if (abs(Jet_eta[ii]-Jet_eta[jj])<5.) continue;
            vecSetPtEtaPhiM(v1,Jet_pt[ii],Jet_eta[ii],Jet_phi[ii],Jet_mass[ii]);
            vecSetPtEtaPhiM(v2,Jet_pt[jj],Jet_eta[jj],Jet_phi[jj],Jet_mass[jj]);
            if ((v1+v2).M()<300.) continue;
            return RVec<Int_t> {ii,jj};
        }
    }
    return RVec<Int_t> {-1,-1};
}
void anaCode() {
    //ROOT::EnableImplicitMT(0);
    TFile *f = TFile::Open("VBF_HHTo4B_CV_1_5_C2V_1_C3_1_mc_10K_ok.root");
    TTree *t1 = (TTree *)f->Get("btagana/ttree");
    TTree *t2 = (TTree *)f->Get("btaganaFatJets/ttree");
    t2->AddFriend(t1,"t1");
    RDataFrame d0(*t2);

    auto d = d0;
    auto d1 = d.Filter("FatJetInfo.nJet>=2","nJet>=2");
    auto d2 = d1.Filter("abs(FatJetInfo.Jet_eta[0])<3. && abs(FatJetInfo.Jet_eta[1])<3.","eta cut");
    auto d3 = d2.Filter("abs(FatJetInfo.Jet_pt[0])>300. && abs(FatJetInfo.Jet_pt[1])>300.","pt cut")
        .Define("DeepAK8_jet0","FatJetInfo.Jet_DeepAK8_ZHbb[0]").Define("DeepAK8_jet1","FatJetInfo.Jet_DeepAK8_ZHbb[1]");
    auto d4 = d3.Filter("DeepAK8_jet0 <= 1. && DeepAK8_jet0 > 0.8 && DeepAK8_jet1 <= 1. && DeepAK8_jet1 > 0.8","DeepAK8 cut");
    auto d5 = d4.Filter("FatJetInfo.Jet_massSoftDrop[0]>90. && FatJetInfo.Jet_massSoftDrop[0]<140.","Fat Jet SDmass cut 1")
        .Filter("FatJetInfo.Jet_massSoftDrop[1]>90.&& FatJetInfo.Jet_massSoftDrop[1]<140.","Fat Jet SDmass cut 2");

    auto d6 = d5.Filter("FatJetInfo.Jet_tau2[0]/FatJetInfo.Jet_tau1[0]<0.6 && FatJetInfo.Jet_tau2[1]/FatJetInfo.Jet_tau1[1]<0.6","tau21 cut");
    auto b = d6.Define("lead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJetInfo.Jet_pt[0],FatJetInfo.Jet_eta[0],FatJetInfo.Jet_phi[0],FatJetInfo.Jet_mass[0])");
    auto bb = b.Define("sublead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJetInfo.Jet_pt[1],FatJetInfo.Jet_eta[1],FatJetInfo.Jet_phi[1],FatJetInfo.Jet_mass[1])");
    auto bbb = bb.Define("ind1"     ,"doSelection(t1.nJet,t1.Jet_pt,t1.Jet_eta,t1.Jet_phi,t1.Jet_mass,lead_fatjet,sublead_fatjet)[0]");
    auto bbbb = bbb.Define("ind2"   ,"doSelection(t1.nJet,t1.Jet_pt,t1.Jet_eta,t1.Jet_phi,t1.Jet_mass,lead_fatjet,sublead_fatjet)[1]");
    auto finalNode = bbbb.Filter("ind1>=0&&ind2>=0","vbf cut").Define("lead_fatjet_sdmass","FatJetInfo.Jet_massSoftDrop[0]").Define("sublead_fatjet_sdmass","FatJetInfo.Jet_massSoftDrop[1]");
    auto defineNode = finalNode.Define("vbf_jet0","ROOT::Math::PtEtaPhiMVector(t1.Jet_pt[ind1],t1.Jet_eta[ind1],t1.Jet_phi[ind1],t1.Jet_mass[ind1])")
                               .Define("vbf_jet1","ROOT::Math::PtEtaPhiMVector(t1.Jet_pt[ind2],t1.Jet_eta[ind2],t1.Jet_phi[ind2],t1.Jet_mass[ind2])")
                               .Define("HHinvMass","(lead_fatjet+sublead_fatjet).M()")
                               .Define("vbfjet_invmass","(vbf_jet0+vbf_jet1).M()");
    cout << "running" << endl;
    initializer_list< std::string > outputlist = {"lead_fatjet","sublead_fatjet","lead_fatjet_sdmass","sublead_fatjet_sdmass","DeepAK8_jet0","DeepAK8_jet1","nJet","Jet_pt","Jet_eta","Jet_phi","Jet_mass","ind1","ind2","vbf_jet0","vbf_jet1","HHinvMass","vbfjet_invmass"};
    //initializer_list< std::string > outputlist = {"lead_fatjet","sublead_fatjet","lead_fatjet_sdmass","sublead_fatjet_sdmass","nJet","Jet_pt","Jet_eta","Jet_phi","Jet_mass","ind1","ind2","vbf_jet0","vbf_jet1","HHinvMass","vbfjet_invmass"};
    defineNode.Snapshot("mytree","database.root",outputlist);
    defineNode.Report()->Print();

}
