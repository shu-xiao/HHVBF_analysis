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

    TFile *fSig = TFile::Open("phase2An/myVBF_HHTo4B_CV_1_C2V_2_C3_1_Phase2tree.root");
    TFile *fTT = TFile::Open("phase2An/myttbartree.root");
    TFile *fQCD_600toInf = TFile::Open("mytree_QCD_Pt_600toInf_nTuple.root");
    TFile *fQCD_470to600 = TFile::Open("mytree_QCD_Pt_470to600_nTuple.root");
    TFile *fQCD_300to470 = TFile::Open("mytree_QCD_Pt_300to470_nTuple.root");

    TTree *tSig = (TTree *)fSig->Get("mytree");
    TTree *tTT = (TTree *)fTT->Get("mytree");
    TTree *tQCD_600toInf = (TTree *)fQCD_600toInf->Get("mytree");
    TTree *tQCD_470to600 = (TTree *)fQCD_470to600->Get("mytree");
    TTree *tQCD_300to470 = (TTree *)fQCD_300to470->Get("mytree");

    // Unit: pb
    const double xsTT = 984.50;
    const double xsSig = 1;
    //const double xsQCD = 0.000165445;
    const double xsQCD_600 = 244.3;
    const double xsQCD_470 = 648.2;
    const double xsQCD_300 = 7823;
    const int sigStr = 500;

    double wTT = getWeight(xsTT, tTT->GetEntries());
    double wSig = getWeight(xsSig, tSig->GetEntries());
    double wQCD_600 = getWeight(xsQCD_600, tQCD_600toInf->GetEntries());
    double wQCD_470 = getWeight(xsQCD_470, tQCD_470to600->GetEntries());
    double wQCD_300 = getWeight(xsQCD_300, tQCD_300to470->GetEntries());

    cout << "wTT:\t" << wTT << endl;
    cout << "wSig:\t" << wSig << endl;
    cout << "wQCD_600:\t" << wQCD_600 << endl;
    cout << "wQCD_470:\t" << wQCD_470 << endl;
    cout << "wQCD_300:\t" << wQCD_300 << endl;

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
            TH1D *hQCD600 = new TH1D(Form("hQCD600_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            TH1D *hQCD470 = new TH1D(Form("hQCD470_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            TH1D *hQCD300 = new TH1D(Form("hQCD300_%d%d",i,j),varList[i].data(),histparams[i][0],histparams[i][1],histparams[i][2]);
            
            hTT->SetLineColor(75);
            hTT->SetFillColor(75);
            hQCD300->SetLineColor(kBlue);
            hQCD300->SetFillColor(kBlue);
            hQCD470->SetLineColor(kBlue+2);
            hQCD470->SetFillColor(kBlue+2);
            hQCD600->SetLineColor(kBlue+4);
            hQCD600->SetFillColor(kBlue+4);
            hsig->SetLineColor(kRed);
            
            tSig->Draw(Form("%s >> hsig_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tTT->Draw(Form("%s >> hTT_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tQCD_600toInf->Draw(Form("%s >> hQCD600_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tQCD_470to600->Draw(Form("%s >> hQCD470_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));
            tQCD_300to470->Draw(Form("%s >> hQCD300_%d%d",varList[i].data(),i,j), Form("%d || (DeepAK8_jet0 > 0.8 && DeepAK8_jet1 > 0.8)",j));

            hsig->Scale(wSig*sigStr);
            hTT->Scale(wTT);
            hQCD600->Scale(wQCD_600);
            hQCD470->Scale(wQCD_470);
            hQCD300->Scale(wQCD_300);
            
            
            leg->Clear();
            leg->AddEntry(hQCD600,"QCD_Pt600toInf","lf");
            leg->AddEntry(hQCD470,"QCD_Pt470to600","lf");
            leg->AddEntry(hQCD300,"QCD_Pt300to470","lf");
            leg->AddEntry(hTT,"TTbar","lf");
            leg->AddEntry(hsig,Form("signal X %d",sigStr),"lf");
            
            /*
            hTT->Add(hQCD);
            hQCD->Draw("hist");
            c1->Print(pdfName.data());
            hTT->Draw("hist");
            c1->Print(pdfName.data());
            */
            //hTT->Add(hQCD);
            
            //THStack *h_bg = new THStack("h_bg","");
            THStack h_bg("h_bg",varList[i].data());
            h_bg.Add(hTT,"hist");
            h_bg.Add(hQCD300,"hist");
            h_bg.Add(hQCD470,"hist");
            h_bg.Add(hQCD600,"hist");

            //if (!j) hTT->SetTitle((string(hTT->GetTitle())+"_DeepAK8_M").data());
            h_bg.Draw();
            hsig->Draw("histsame");
            leg->Draw();
            c1->Print(pdfName.data());
           
            // draw normalized
            if (i>4 && j ) {

                auto h_sumQCD = (TH1F*)hQCD300->Clone(Form("h_sumQCD_%d%d",i,j));
                h_sumQCD->Add(hQCD470);
                h_sumQCD->Add(hQCD600);

                hTT->SetFillStyle(3001);
                hsig->SetFillStyle(3001);
                h_sumQCD->SetFillStyle(3001);
                
                leg->Clear();
                leg->AddEntry(h_sumQCD,"QCD","lf");
                leg->AddEntry(hTT,"TTbar","lf");
                leg->AddEntry(hsig,"signal","lf");
                h_sumQCD->SetTitle(("No Stack Normalized plot - " + (string)h_sumQCD->GetTitle()).data());
                h_sumQCD->DrawNormalized("hist");
                hTT->DrawNormalized("histsame");
                hsig->DrawNormalized("histsame");
                leg->Draw();
                c1->Print(pdfName.data());
                c1->SetLogy();
                c1->Update();
                c1->Print(pdfName.data());
                c1->SetLogy(0);
            }
        }
    }
    c1->Print((pdfName+"]").data());
}
