#include "setNCUStyle.C"
const double LPHASE2 = 3000*1000; // unit pb
inline double getWeight(double xs, double nEve) {
    return LPHASE2*xs/nEve;
}
void setMax(TH1D* h1, TH1D* h2, TH1D* h3, double scale=1.2) {
    double m1 = h1->GetMaximum();
    double m2 = h2->GetMaximum();
    double m3 = h3->GetMaximum();
    double hmax = max(max(m1,m2),m3)*scale;
    h1->SetMaximum(hmax);
    h2->SetMaximum(hmax);
    h3->SetMaximum(hmax);
}
void plotAll() {
    setNCUStyle(1);

    //gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();

    TFile *fSig[6];
    TTree *tSig[6];
    string parList[] = {"CV_0_5_C2V_1_C3_1","CV_1_C2V_1_C3_0","CV_1_C2V_1_C3_2","CV_1_5_C2V_1_C3_1","CV_1_C2V_1_C3_1","CV_1_C2V_2_C3_1"};
    for (int i=0;i<6;i++) {
        fSig[i] = TFile::Open(Form("mytree_VBF_HHTo4B_%s_nTuple_PU200.root",parList[i].data()));
        tSig[i] = (TTree *)fSig[i]->Get("mytree");
    }

    // Unit: pb
    const double xsTT = 984.50;
    const double xsSig = 1;
    //const double xsQCD = 0.000165445;
    const double xsQCD = 244.3;
    const int sigStr = 500;


    vector<string> varList = {"lead_fatjet.Pt()","sublead_fatjet.Pt()","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","DeepAK8_jet0","DeepAK8_jet1","vbfjet_invmass",\
    "Jet_pt[ind1]","Jet_pt[ind2]","Jet_eta[ind1]","Jet_eta[ind2]","abs(-Jet_eta[ind2]+Jet_eta[ind1])"};
    //vector<string> varList = {"lead_fatjet.Pt()","sublead_fatjet.Pt()","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1"};
    vector<string> titleList = {"LeadAK8_Pt","TrailAK8_Pt","Lead_AK8jet_sdmass","Trail_AK8jet_sdmass","HHinvMass","DeepAK8_leadAK8","DeepAK8_trailAK8","VBFJet_invmass","LeadVBFJet_Pt","TrailVBFJet_Pt","LeadVBFJet_eta","TrailVBFJet_eta","VBFJets_etadiff"};
    vector<vector<float>> histparams = {   
                                {25,300,800},
                                {25,300,800},
                                {35,80,150},
                                {35,80,150},
                                {30,600,2100},
                                //{60,0,0.6},
                                //{60,0,0.6},
                                {25,0.75,1},
                                {25,0.75,1},
                                {25,0,5000},
                                {50,0,500},
                                {50,0,500},
                                {25,-5,5},
                                {25,-5,5},
                                {40,5,9}};

    TCanvas *c1 = new TCanvas("c1","c1",1800,900);
    c1->Divide(3,2);
    string pdfName = "sigPlot.pdf";
    c1->Print((pdfName+"[").data());

    // TT >> QCD
   auto leg = new TLegend(0.6,0.7,0.9,0.9);

    for (int i=0;i<varList.size();i++) {

        for (int j=0;j<6;j++) {
            
            c1->cd(j+1);
            TH1D *hsig = new TH1D(Form("hsig_%d%d",i,j),parList[j].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            //TH1D *hsig = new TH1D(Form("hsig_%d%d",i,j),Form("%s_%s",titleList[i].data(),parList[j].data()),histparams[i][0],histparams[i][1],histparams[i][2]);
            
            //hsig->SetLineColor(kRed);

            //tSig[j]->Draw(Form("%s",varList[i].data()));
            gPad->SetLeftMargin(0.15);
            tSig[j]->Draw(Form("%s >> hsig_%d%d",varList[i].data(),i,j));
            hsig->GetYaxis()->SetTitle("Normalized to 1");
            hsig->GetXaxis()->SetTitle(titleList[i].data());
            hsig->DrawNormalized("hist");
            
            
            //hsig->Scale(wSig*sigStr);
            //hTT->Scale(wTT);
            //hQCD->Scale(wQCD);
            /*
            leg->Clear();
            leg->AddEntry(hQCD,"QCD","lf");
            leg->AddEntry(hTT,"TTbar","lf");
            leg->AddEntry(hsig,Form("signal X %d",sigStr),"lf");
            hQCD->Draw("hist");
            c1->Print(pdfName.data());
            hTT->Draw("hist");
            c1->Print(pdfName.data());
            */

            //hsig->Draw("hist");
            //leg->Draw();
            
        }
        c1->Print(pdfName.data());
    }
    c1->Print((pdfName+"]").data());
}
