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

    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();

    TFile *fSig = TFile::Open("myVBF_HHTo4B_CV_1_C2V_2_C3_1_Phase2tree.root");
    TFile *fTT = TFile::Open("myttbartree.root");
    TFile *fQCD = TFile::Open("myQCD_Pt_600oInf_mctree.root");

    TTree *tSig = (TTree *)fSig->Get("mytree");
    TTree *tTT = (TTree *)fTT->Get("mytree");
    TTree *tQCD = (TTree *)fQCD->Get("mytree");

    // Unit: pb
    const double xsTT = 984.50;
    const double xsSig = 1;
    //const double xsQCD = 0.000165445;
    const double xsQCD = 244.3;
    const int sigStr = 500;

    double wTT = getWeight(xsTT, tTT->GetEntries());
    double wSig = getWeight(xsSig, tSig->GetEntries());
    double wQCD = getWeight(xsQCD, tQCD->GetEntries());

    cout << "wTT:\t" << wTT << endl;
    cout << "wSig:\t" << wSig << endl;
    cout << "wQCD:\t" << wQCD << endl;

    vector<string> varList = {"lead_fatjet.Pt()","sublead_fatjet.Pt()","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1"};
    vector<vector<float>> histparams = {   
                                {50,300,800},
                                {50,300,800},
                                {35,80,150},
                                {35,80,150},
                                {30,600,2100},
                                {60,0,0.6},
                                {60,0,0.6},
                                {50,0,1},
                                {50,0,1} };

    TCanvas *c1 = new TCanvas("c1","c1",3);
    string pdfName = "sumPlot.pdf";
    c1->Print((pdfName+"[").data());

    // TT >> QCD
   auto leg = new TLegend(0.6,0.7,0.9,0.9);

    for (int i=0;i<varList.size();i++) {

        for (int j=1;j>=0;j--) {
            
            if (j && i>8) continue;

            TH1D *hsig = new TH1D(Form("hsig_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            TH1D *hTT  = new TH1D(Form("hTT_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            TH1D *hQCD = new TH1D(Form("hQCD_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            
            hTT->SetLineColor(55);
            hTT->SetFillColor(55);
            hQCD->SetLineColor(75);
            hQCD->SetFillColor(75);
            hsig->SetLineColor(kRed);

            tSig->Draw(Form("%s >> hsig_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tTT->Draw(Form("%s >> hTT_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tQCD->Draw(Form("%s >> hQCD_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            
            
            hsig->Scale(wSig*sigStr);
            hTT->Scale(wTT);
            hQCD->Scale(wQCD);
            
            leg->Clear();
            leg->AddEntry(hQCD,"QCD","lf");
            leg->AddEntry(hTT,"TTbar","lf");
            leg->AddEntry(hsig,Form("signal X %d",sigStr),"lf");
            /*
            hQCD->Draw("hist");
            c1->Print(pdfName.data());
            hTT->Draw("hist");
            c1->Print(pdfName.data());
            */

            hTT->Add(hQCD);
            if (!j) hTT->SetTitle((string(hTT->GetTitle())+"_DeepAK8_M").data());
            hTT->Draw("hist");
            hQCD->Draw("histsame");
            hsig->Draw("histsame");
            leg->Draw();
            c1->Print(pdfName.data());
            
            if (i>4 && j) {
                hTT->SetFillStyle(3001);
                hsig->SetFillStyle(3001);
                hQCD->SetFillStyle(3001);
                hTT->Add(hQCD,-1);
                hsig->Scale(1./hsig->Integral());
                hQCD->Scale(1./hQCD->Integral());
                hTT->Scale(1./hTT->Integral());
                setMax(hsig,hQCD,hTT);
                hTT->SetTitle((string(hTT->GetTitle())+"_NormToOne").data());
                hTT->DrawNormalized("hist");
                hQCD->DrawNormalized("histsame");
                hsig->DrawNormalized("histsame");
                leg->Draw();
                c1->Print(pdfName.data());
            }
        }
    }
    c1->Print((pdfName+"]").data());
}
