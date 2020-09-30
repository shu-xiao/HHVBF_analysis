
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
        if (abs(Jet_eta[i])>5.) continue;
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
RVec<Int_t> vbfBasis(Int_t nThinJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, TLvector Hjet1, TLvector Hjet2) {
    RVec<Int_t> vbfj;
    const Float_t dRCut = 1.2;
    for (Int_t i=0;i<nThinJet;i++) {
        if (Jet_pt[i]<50.) continue;
        if (abs(Jet_eta[i])>5.) continue;
        vbfj.push_back(i);
    }
    if (vbfj.size()<2) return RVec<Int_t> {-1,-1};
    for (auto it = vbfj.begin(); it != vbfj.end(); ) {
        if (deltaR(Hjet1.Eta(),Hjet1.Phi(),Jet_eta[vbfj[*it]],Jet_phi[vbfj[*it]])>dRCut &&
            deltaR(Hjet2.Eta(),Hjet2.Phi(),Jet_eta[vbfj[*it]],Jet_phi[vbfj[*it]])>dRCut) ++it; 
        else it = vbfj.erase(it);
    }
    if (vbfj.size()<2) return RVec<Int_t> {-2,-2};
    else return vbfj;
}
RVec<Int_t> vbfCut(rvec_i vbfj, rvec_f Jet_eta) {
    if (vbfj[0]<0) return  RVec<Int_t> {-1}; 
    RVec<Int_t> vbfjCom;
    for (UInt_t i=0;i<vbfj.size();i++) {
        UInt_t ii = vbfj[i];
        for (UInt_t j=0;j<i;j++) {
            UInt_t jj = vbfj[j];
            if (Jet_eta[ii]*Jet_eta[jj]>=0.) continue;
            if (abs(Jet_eta[ii]-Jet_eta[jj])<5.) continue;
            vbfjCom.push_back(jj*100+ii);
        }
    }
    if (vbfjCom.size()>0) return vbfjCom;
    else return RVec<Int_t> {-1111};
}
Int_t vbfInv (rvec_i vbfjCom,rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, Float_t invMCut=300.) {
    if (vbfjCom.size()<1 || vbfjCom[0] < 0) return -999;
    ROOT::Math::PtEtaPhiMVector v1, v2;
    Int_t ii,jj;
    for (UInt_t i=0;i<vbfjCom.size();i++) {
        ii = vbfjCom[i]/100;
        jj = vbfjCom[i]%100;
        vecSetPtEtaPhiM(v1,Jet_pt[ii],Jet_eta[ii],Jet_phi[ii],Jet_mass[ii]);
        vecSetPtEtaPhiM(v2,Jet_pt[jj],Jet_eta[jj],Jet_phi[jj],Jet_mass[jj]);
        if ((v1+v2).M()<invMCut) continue;
        return vbfjCom[i];
    }
    return -1;

}
void nMinus1_anaCode(string inputFile="VBF_HHTo4B_CV_1_5_C2V_1_C3_1_mc_10K_ok.root",string outputFile="signalDataBase.root") {
    //ROOT::EnableImplicitMT(2);
    TFile *f = TFile::Open(inputFile.data());
    TTree *t1 = (TTree *)f->Get("btagana/ttree");
    TTree *t2 = (TTree *)f->Get("btaganaFatJets/ttree");
    t2->AddFriend(t1,"t1");
    RDataFrame d0(*t2);

    //auto d = d0.Range(500);
    auto d = d0;
    auto d1 = d.Filter("FatJetInfo.nJet>=2","nJet>=2");
    auto d2 = d1.Filter("abs(FatJetInfo.Jet_eta[0])<3. && abs(FatJetInfo.Jet_eta[1])<3.","eta cut");
    auto d3 = d2.Filter("abs(FatJetInfo.Jet_pt[0])>300. && abs(FatJetInfo.Jet_pt[1])>300.","pt cut")
        .Define("DeepAK8_jet0","FatJetInfo.Jet_DeepAK8_ZHbb[0]").Define("DeepAK8_jet1","FatJetInfo.Jet_DeepAK8_ZHbb[1]");
    auto d4 = d3.Define("tau21_0","FatJetInfo.Jet_tau2[0]/FatJetInfo.Jet_tau1[0]").Define("tau21_1","FatJetInfo.Jet_tau2[1]/FatJetInfo.Jet_tau1[1]").Filter("tau21_0 < 0.6 && tau21_1 < 0.6","tau21 cut");

    auto d5 = d4.Filter("FatJetInfo.Jet_massSoftDrop[0]>90. && FatJetInfo.Jet_massSoftDrop[0]<140.","Fat Jet SDmass cut 1")
        .Filter("FatJetInfo.Jet_massSoftDrop[1]>90.&& FatJetInfo.Jet_massSoftDrop[1]<140.","Fat Jet SDmass cut 2");
    auto d6 = d5.Filter("DeepAK8_jet0 <= 1. && DeepAK8_jet0 >= 0.0 && DeepAK8_jet1 <= 1. && DeepAK8_jet1 >= 0.0","DeepAK8 basis");

    auto d7 = d6.Define("lead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJetInfo.Jet_pt[0],FatJetInfo.Jet_eta[0],FatJetInfo.Jet_phi[0],FatJetInfo.Jet_mass[0])").Define("lead_fatjet_sdmass","FatJetInfo.Jet_massSoftDrop[0]");
    auto d8 = d7.Define("sublead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJetInfo.Jet_pt[1],FatJetInfo.Jet_eta[1],FatJetInfo.Jet_phi[1],FatJetInfo.Jet_mass[1])").Define("sublead_fatjet_sdmass","FatJetInfo.Jet_massSoftDrop[1]")
        .Define("HHinvMass","(lead_fatjet+sublead_fatjet).M()");

    auto b0 = d8.Filter("DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8","DeepAK8 cut");

    auto b1 = b0.Define("vbfjList" ,"vbfBasis(t1.nJet,t1.Jet_pt,t1.Jet_eta,t1.Jet_phi,t1.Jet_mass,lead_fatjet,sublead_fatjet)").Filter("vbfjList[0]>=0","AK4 jet basic cuts");
    auto b2 = b1.Define("vbfjKinList","vbfCut(vbfjList, t1.Jet_eta)").Filter("vbfjKinList[0]>0","AK4 eta cuts");
    auto b3 = b2.Define("vbfjCom","vbfInv(vbfjKinList,t1.Jet_pt,t1.Jet_eta,t1.Jet_phi,t1.Jet_mass,300)")
        .Filter("vbfjCom>=0","vbf jet invariance mass cut").Define("ind1","vbfjCom/100").Define("ind2","vbfjCom%100");
    /*
    auto finalNode = d8.Define("vbf_jet0","ROOT::Math::PtEtaPhiMVector(t1.Jet_pt[ind1],t1.Jet_eta[ind1],t1.Jet_phi[ind1],t1.Jet_mass[ind1])")
                               .Define("vbf_jet1","ROOT::Math::PtEtaPhiMVector(t1.Jet_pt[ind2],t1.Jet_eta[ind2],t1.Jet_phi[ind2],t1.Jet_mass[ind2])")
                               .Define("vbfjet_invmass","(vbf_jet0+vbf_jet1).M()");
    */
    auto finalNode = d8;
                               
    cout << "running" << endl;
    //initializer_list< std::string > outputlist = {"lead_fatjet","sublead_fatjet","lead_fatjet_sdmass","sublead_fatjet_sdmass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1","nJet","Jet_pt","Jet_eta","Jet_phi","Jet_mass","ind1","ind2","vbf_jet0","vbf_jet1","HHinvMass","vbfjet_invmass"};
    initializer_list< std::string > outputlist = {"lead_fatjet","sublead_fatjet","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1"};
    finalNode.Snapshot("mytree",outputFile.data(),outputlist);
    finalNode.Report()->Print();

}
