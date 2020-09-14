#include <TFile.h>
#include "ROOT/RVec.hxx"
#include <TLorentzVector.h>
#include "Math/Vector4D.h"

using namespace ROOT::VecOps;
using rvec_i = const RVec<int> &;
using rvec_f = const RVec<float> &;
using rvec_TLvector = const RVec<ROOT::Math::PtEtaPhiMVector> &;
using TLvector = const ROOT::Math::PtEtaPhiMVector &;
//using rvec_TLvector = const RVec<Math::PtEtaPhiMVector> &;
//using TLvector = const Math::PtEtaPhiMVector &;
using namespace ROOT::Math;
inline Float_t deltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
    Float_t deltaRSquare = (eta2-eta1)*(eta2-eta1)+(phi2-phi1)*(phi2-phi1);
    return TMath::Sqrt(deltaRSquare);
}
inline void vecSetPtEtaPhiM(ROOT::Math::PtEtaPhiMVector& vec, Float_t pt, Float_t eta, Float_t phi, Float_t m) {
    vec.SetPt(pt);
    vec.SetEta(eta);
    vec.SetPhi(phi);
    vec.SetM(m);
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

RVec<Int_t> doSelection(Int_t nThinJet, rvec_TLvector thinJetP4, TLvector Hjet1, TLvector Hjet2 ) {
    const Float_t dRCut = 1.2;
    TLorentzVector v1, v2;
    RVec<Int_t> vbfj;
    for (Int_t i=0;i<nThinJet;i++) {
        if (thinJetP4[i].Pt()<50.) continue;
        if (thinJetP4[i].Eta()>5.) continue;
        if (deltaR(Hjet1.Eta(),Hjet1.Phi(),thinJetP4[i].Eta(),thinJetP4[i].Phi())<dRCut) continue;
        if (deltaR(Hjet2.Eta(),Hjet2.Phi(),thinJetP4[i].Eta(),thinJetP4[i].Phi())<dRCut) continue;
        vbfj.push_back(i);
    }
    for (UInt_t j=0;j<vbfj.size();j++) {
        Int_t jj = vbfj[j];
        for (UInt_t i=0;i<j;i++) {
            Int_t ii = vbfj[i];
            if (thinJetP4[ii].Eta()*thinJetP4[jj].Eta()>=0.) continue;
            if (abs(thinJetP4[ii].Eta()-thinJetP4[jj].Eta())<5.) continue;
            if ((thinJetP4[ii]+thinJetP4[jj]).M()<300.) continue;
            return RVec<Int_t> {ii,jj};
        }
    }
    return RVec<Int_t> {-1,-1};
}

RVec<Int_t> doSelection(Int_t nThinJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, rvec_f FatJet_pt, rvec_f FatJet_eta, rvec_f FatJet_phi, rvec_f FatJet_mass ) {
    const Float_t dRCut = 1.2;
    ROOT::Math::PtEtaPhiMVector v1, v2;
    RVec<Int_t> vbfj;
    for (Int_t i=0;i<nThinJet;i++) {
        if (Jet_pt[i]<50.) continue;
        if (Jet_eta[i]>5.) continue;
        if (deltaR(FatJet_eta[0],FatJet_phi[0],Jet_eta[i],Jet_phi[i])<dRCut) continue;
        if (deltaR(FatJet_eta[1],FatJet_phi[1],Jet_eta[i],Jet_phi[i])<dRCut) continue;
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
