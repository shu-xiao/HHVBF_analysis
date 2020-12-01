
#include <TFile.h>
#include "ROOT/RVec.hxx"
#include <TLorentzVector.h>
#include "Math/Vector4D.h"
#include <ROOT/RDataFrame.hxx>

#define runShort 0
#define CONESIZE 0.4

using namespace ROOT;
using rvec_i = const RVec<int> &;
using rvec_f = const RVec<float> &;
using rvec_TLvector = const RVec<ROOT::Math::PtEtaPhiMVector> &;
using TLvector = const ROOT::Math::PtEtaPhiMVector &;
void setMax(TH1D &h1, TH1D &h2) {
    float a = h1.GetMaximum();
    float b = h2.GetMaximum();
    float max;
    if (a>b) max = a;
    else max = b;
    h1.SetMaximum(max*1.5);
    h2.SetMaximum(max*1.5);
}
inline void vecSetPtEtaPhiM(ROOT::Math::PtEtaPhiMVector& vec, Float_t pt, Float_t eta, Float_t phi, Float_t m) {
    vec.SetPt(pt);
    vec.SetEta(eta);
    vec.SetPhi(phi);
    vec.SetM(m);
}
Double_t invPt(Int_t nJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass) {
    ROOT::Math::PtEtaPhiMVector vec(0.,0.,0.,0.);
    //for (int i=0;i<nJet;i++) vec+=ROOT::Math::PtEtaPhiMVector(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
    for (int i=0;i<2;i++) vec+=ROOT::Math::PtEtaPhiMVector(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
    return vec.Pt();
}
RVec<Int_t> doSelection(Int_t nThinJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, TLvector Hjet1, TLvector Hjet2) ;
RVec<Int_t> vbfBasis(Int_t nThinJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, TLvector Hjet1, TLvector Hjet2) ;
RVec<Int_t> vbfCut(rvec_i vbfj, rvec_f Jet_eta) ;
Int_t vbfInv (rvec_i vbfjCom,rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, Float_t invMCut) ;

RVec<Int_t> matchInd(Int_t nGenJet, rvec_i PDGList) {
    RVec<Int_t> ind = {-1,-1};
    for (int i=0;i<nGenJet;i++) {
        if (PDGList[i]==25 && ind[0]<0) {ind[0] = i;continue;}
        else if (PDGList[i]==25 && ind[0]>=0) {ind[1] = i; break;}
    }
    return ind;
}
//RVec<Int_t> matchJet(Int_t ind1, Int_t ind2, rvec_f GenPruned_pT, rvec_f GenPruned_eta, rvec_f GenPruned_phi, rvec_f GenPruned_mass,Int_t nJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass, rvec_f Jet_sdMass) {
RVec<Int_t> matchJet(Int_t ind1, Int_t ind2, rvec_f GenPruned_pT, rvec_f GenPruned_eta, rvec_f GenPruned_phi, rvec_f GenPruned_mass,Int_t nJet, rvec_f Jet_pt, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_mass) {
    TLorentzVector vecJet, vecGen1, vecGen2;
    vecGen1.SetPtEtaPhiM(GenPruned_pT[ind1],GenPruned_eta[ind1],GenPruned_phi[ind1],GenPruned_mass[ind1]);
    vecGen2.SetPtEtaPhiM(GenPruned_pT[ind2],GenPruned_eta[ind2],GenPruned_phi[ind2],GenPruned_mass[ind2]);
    //ROOT::Math::PtEtaPhiMVector vecGen1(GenPruned_pT[ind1],GenPruned_eta[ind1],GenPruned_phi[ind1],GenPruned_mass[ind1]);
    //ROOT::Math::PtEtaPhiMVector vecGen2(GenPruned_pT[ind2],GenPruned_eta[ind2],GenPruned_phi[ind2],GenPruned_mass[ind2]);
    RVec<Int_t> ind = {-1,-1};
    for (int i=0;i<nJet;i++) {
        //ROOT::Math::PtEtaPhiMVector vecJet(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
        if (Jet_mass[i] <= 0) continue;
        //if (Jet_sdMass[i] <= 0) continue;
        vecJet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
        if (vecJet.DeltaR(vecGen1) < CONESIZE && ind[0]<0) ind[0] = i;
        if (vecJet.DeltaR(vecGen2) < CONESIZE && ind[1]<0) ind[1] = i;
        if (ind[0]>=0 && ind[1]>=0) {
            if (Jet_pt[ind[1]]>Jet_pt[ind[0]]) swap(ind[0],ind[1]);
            break;
        }
    }
    return ind;
}

void nMinus1_genMatch() {
    gStyle->SetOptStat(110011);

    TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
    //ROOT::EnableImplicitMT(2);  // error
    TH1::SetDefaultSumw2(1);
    string parList[] = {"CV_0_5_C2V_1_C3_1","CV_1_C2V_1_C3_0","CV_1_C2V_1_C3_2","CV_1_5_C2V_1_C3_1","CV_1_C2V_1_C3_1","CV_1_C2V_2_C3_1"};
    string puList[] = {"PU200","NoPU"};

    int spCase[] = {0,3,4};
    vector<vector<TFile*>> fList(2); 
    //TFile* fList[6][2];
    vector<TChain*> pu200List;
    vector<TH2D> pdgidList;
    vector<TH2D> AK8jetIDList;
    const int nHist = 4;
    vector<vector<TH1D>> histList(nHist);

    string pdfname = Form("matchInd_R%d.pdf",int(CONESIZE*10));
    c1->Print((pdfname+"[").data());

    c1->Divide(3,2);
    string xTitleList[nHist] = {"sdmass1","sdmass2","AK8mass1","AK8mass2"};
    //for (int i=11;i<12;i++) {
    for (int i=0;i<12;i++) {
        if (i!=0 && runShort) continue;
        // VBF_HHTo4B_CV_1_5_C2V_1_C3_1_nTuple_NoPU_1.root
        /*
        TFile *f = TFile::Open(Form("VBF_HHTo4B_%s_nTuple_PU200.root",parList[i].data()));
        TTree *t1 = (TTree *)f->Get("btagana/ttree");
        TTree *t2 = (TTree *)f->Get("btaganaFatJets/ttree");
        //t2->AddFriend(t1,"t1");
        */
        TChain ch1("btagana/ttree");
        TChain ch2("btaganaFatJets/ttree");
        ch1.SetAutoDelete(0);
        ch2.SetAutoDelete(0);
        if (i<6) {
            ch1.Add(Form("VBF_HHTo4B_%s_nTuple_PU200.root",parList[i].data()));
            ch2.Add(Form("VBF_HHTo4B_%s_nTuple_PU200.root",parList[i].data()));
        }
        else if (i==6||i==9||i==10) {
            ch1.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU_1.root",parList[i%6].data()));
            ch1.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU_2.root",parList[i%6].data()));
            ch2.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU_1.root",parList[i%6].data()));
            ch2.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU_2.root",parList[i%6].data()));
        }
        else {
            cout << Form("VBF_HHTo4B_%s_nTuple_NoPU.root",parList[i%6].data()) << endl;
            ch1.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU.root",parList[i%6].data()));
            ch2.Add(Form("VBF_HHTo4B_%s_nTuple_NoPU.root",parList[i%6].data()));
        }
        ch1.AddFriend(&ch2);

        //RDataFrame d0(*t2);
        RDataFrame d0(ch1);
        cout << "loading RDataframe " << i+1 << endl;

        auto d = d0;//.Range(100);
        auto d1 = d.Filter("FatJetInfo.nJet>=2","nJet>=2");
        auto d2 = d1.Define("Hind1","matchInd(nGenPruned,GenPruned_pdgID)[0]").Define("Hind2","matchInd(nGenPruned,GenPruned_pdgID)[1]").Filter("Hind2>=0 && Hind1>=0","Higgs Matching");
        //auto d3 = d2.Define("matchId1","matchJet(Hind1,Hind2,GenPruned_pT,GenPruned_eta,GenPruned_phi,GenPruned_mass,FatJetInfo.nJet,FatJetInfo.Jet_pt,FatJetInfo.Jet_eta,FatJetInfo.Jet_phi,FatJetInfo.Jet_mass,FatJetInfo.Jet_massSoftDrop)[0]")
        //            .Define("matchId2","matchJet(Hind1,Hind2,GenPruned_pT,GenPruned_eta,GenPruned_phi,GenPruned_mass,FatJetInfo.nJet,FatJetInfo.Jet_pt,FatJetInfo.Jet_eta,FatJetInfo.Jet_phi,FatJetInfo.Jet_mass,FatJetInfo.Jet_massSoftDrop)[1]")
        auto d3p = d2.Define("matchId1","matchJet(Hind1,Hind2,GenPruned_pT,GenPruned_eta,GenPruned_phi,GenPruned_mass,FatJetInfo.nJet,FatJetInfo.Jet_pt,FatJetInfo.Jet_eta,FatJetInfo.Jet_phi,FatJetInfo.Jet_mass)[0]")
                    .Define("matchId2","matchJet(Hind1,Hind2,GenPruned_pT,GenPruned_eta,GenPruned_phi,GenPruned_mass,FatJetInfo.nJet,FatJetInfo.Jet_pt,FatJetInfo.Jet_eta,FatJetInfo.Jet_phi,FatJetInfo.Jet_mass)[1]")
                    .Define("sdMass1","FatJetInfo.Jet_massSoftDrop[matchId1]")
                    .Define("sdMass2","FatJetInfo.Jet_massSoftDrop[matchId2]");
        auto d3 = d3p.Define("FatJetMass1","FatJetInfo.Jet_mass[matchId1]")
                    .Define("FatJetMass2","FatJetInfo.Jet_mass[matchId2]");
        d3.Report()->Print();
        auto h_pdgid = d3.Histo2D({"pdgId2D","pdgId2D",6,0,6,6,0,6},"Hind1","Hind2");
        pdgidList.push_back(*h_pdgid);        
        auto h_ind2D = d3.Histo2D({"matchInd2D","matchInd2D",5,0,5,5,0,5},"matchId1","matchId2");
        AK8jetIDList.push_back(*h_ind2D);
        //d3.Histo1D({"sdmass1",Form("SDmass1_%s",parList[i%6].data()),6,-1,5},"matchId1")->Draw("hist");
        //c1->Print("ggg.pdf");
        histList[0].push_back(*d3.Histo1D({"sdmass1",Form("SDmass1_%s",parList[i%6].data()),40,0,160},"sdMass1"));
        histList[1].push_back(*d3.Histo1D({"sdmass2",Form("SDmass2_%s",parList[i%6].data()),40,0,160},"sdMass2"));
        histList[2].push_back(*d3.Histo1D({"AK8mass1",Form("Jetmass1_%s",parList[i%6].data()),40,0,160},"FatJetMass1"));
        histList[3].push_back(*d3.Histo1D({"AK8mass2",Form("Jetmass2_%s",parList[i%6].data()),40,0,160},"FatJetMass2"));
        //c1->cd(i+1);
        //h_ind2D->Draw("colz");
    }
    //c1->Print(pdfname.data());
    //
    cout << "Drawing" << endl;
    
    for (int i=0;i<6;i++) {
        c1->cd(i+1);
        pdgidList[i].Draw("colz");
    }
    c1->Print(pdfname.data());

    for (int i=0;i<12;i++) {
        c1->cd(i%6+1);
        AK8jetIDList[i].SetTitle(Form("%s_%s;genmatch AK8 index1;genmatch AK8 index2",parList[i%6].data(),puList[i/6].data()));
        AK8jetIDList[i].Draw("colz text");
        if (i==5||i==11) c1->Print(pdfname.data());
    }
    //c1->Print(pdfname.data());

    gStyle->SetOptStat(0);
    TLegend *leg = 0;
    for (int n=0;n<nHist;n++) {

        for (int i=0;i<6;i++) {
            c1->cd(i+1);
            // bug
            //if (n<2) histList[n][i].SetTitle(Form("%s;%s",parList[i].data(),xTitleList[i].data()));
            histList[n][i+6].GetYaxis()->SetTitle("Normalize to 1");
            histList[n][i+6].SetLineColor(kBlack);
            histList[n][i].SetLineColor(kBlue);
            histList[n][i].Scale(1./histList[n][i].Integral());
            histList[n][i+6].Scale(1./histList[n][i+6].Integral());
            setMax(histList[n][i+6],histList[n][i]);
            if (!leg) {
                leg = new TLegend(0.7,0.7,0.9,0.9);
                leg->AddEntry(&histList[n][i+6],"NoPU");
                leg->AddEntry(&histList[n][i],"PU 200");
            }
            histList[n][i+6].DrawNormalized("hist");
            histList[n][i].DrawNormalized("histsame");
            leg->Draw();
            //gPad->BuildLegend(0.8,0.75,0.99,0.95,"","");
        }
        c1->Print(pdfname.data());
    }
    c1->Print((pdfname+"]").data());

}

