void effiSum_noTau21() {
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    //gStyle->SetLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetPaintTextFormat("2.0f");

    TCanvas *c1 = new TCanvas("c1","c1",1800,900);
    // data
    //           column meaning                       nTotal   nJet     eta       pt       SDMass  DeepAk8 AK4   AK4eta invMass
    vector <int> nEve_VBF_HHTo4B_CV_1_5_C2V_1_C3_1 = {500000  ,97910   ,97910    ,60592   ,14124  ,10536  ,3787  ,372  ,372};
    vector <int> nEve_VBF_HHTo4B_CV_0_5_C2V_1_C3_1 = {500000  ,178372  ,178372   ,117110  ,27343  ,20477  ,7318  ,765  ,765};
    vector <int> nEve_VBF_HHTo4B_CV_1_C2V_2_C3_1   = {499656  ,222699  ,222699   ,149427  ,34561  ,25973  ,9413  ,926  ,926};
    vector <int> nEve_VBF_HHTo4B_CV_1_C2V_1_C3_1   = {499652  ,13300   ,13300    ,4339    ,477    ,289    ,134    ,11   ,11};
    vector <int> nEve_VBF_HHTo4B_CV_1_C2V_1_C3_0   = {500000  ,10208   ,10208    ,2800    ,457    ,313    ,126    ,11   ,11};
    vector <int> nEve_VBF_HHTo4B_CV_1_C2V_1_C3_2   = {499308  ,26649   ,26649    ,9038    ,509    ,204    ,94     ,11   ,11};
    vector <int> nEve_QCD_600toInf                 = {50000   ,48487   ,48487    ,45994   ,8659   ,307    ,0};
    vector <int> nEve_TT                           = {300000  ,10270   ,2964     ,1778    ,239    ,22     ,4      ,0};
    
    TH1F * h1 = new TH1F("h1","",10,0,10);
    TH1F * h2 = new TH1F("h2","",10,0,10);
    THStack *hs_effi = new THStack("hffi","cut efficiency X 100");
    THStack *hssig_effi = new THStack("hffi_sig","Signal Efficiency of Each Selection");
    THStack *hs_cumeffi = new THStack("hffi","cumulative cut efficiency");
    THStack *hs_sigeffi = new THStack("hffi_sig","Signal Cumulative Efficiency of Each Selection");
    vector<string> namelist = {"QCD_600toInf","TT","VBF_HHTo4B_CV_1_5_C2V_1_C3_1","VBF_HHTo4B_CV_0_5_C2V_1_C3_1","VBF_HHTo4B_CV_1_C2V_2_C3_1",\
    "VBF_HHTo4B_CV_1_C2V_1_C3_1","VBF_HHTo4B_CV_1_C2V_1_C3_0","VBF_HHTo4B_CV_1_C2V_1_C3_2"};

    TH1F* heffi_QCD_600toInf = new TH1F("h_QCD_600toInf","",10,0,10);
    TH1F* heffi_TT = new TH1F("h_TT","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1 = new TH1F("Effi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1 = new TH1F("Effi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_1_C2V_2_C3_1 = new TH1F("Effi_VBF_HHTo4B_CV_1_C2V_2_C3_1","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_1_C2V_1_C3_1 = new TH1F("Effi_VBF_HHTo4B_CV_1_C2V_1_C3_1","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_1_C2V_1_C3_0 = new TH1F("Effi_VBF_HHTo4B_CV_1_C2V_1_C3_0","",10,0,10);
    TH1F* heffi_VBF_HHTo4B_CV_1_C2V_1_C3_2 = new TH1F("Effi_VBF_HHTo4B_CV_1_C2V_1_C3_2","",10,0,10);

    TH1F* hcumeffi_QCD_600toInf = new TH1F("hcum_QCD_600toInf","",10,0,10);
    TH1F* hcumeffi_TT = new TH1F("hcum_TT","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1 = new TH1F("CumEffi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1 = new TH1F("CumEffi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_1_C2V_2_C3_1 = new TH1F("CumEffi_VBF_HHTo4B_CV_1_C2V_2_C3_1","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_1 = new TH1F("CumEffi_VBF_HHTo4B_CV_1_C2V_1_C3_1","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_0 = new TH1F("CumEffi_VBF_HHTo4B_CV_1_C2V_1_C3_0","",10,0,10);
    TH1F* hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_2 = new TH1F("CumEffi_VBF_HHTo4B_CV_1_C2V_1_C3_2","",10,0,10);

    vector<vector<int>> dataList = {nEve_QCD_600toInf, nEve_TT,\
        nEve_VBF_HHTo4B_CV_1_C2V_2_C3_1,nEve_VBF_HHTo4B_CV_0_5_C2V_1_C3_1,\
        nEve_VBF_HHTo4B_CV_1_5_C2V_1_C3_1,nEve_VBF_HHTo4B_CV_1_C2V_1_C3_1,\
        nEve_VBF_HHTo4B_CV_1_C2V_1_C3_0,nEve_VBF_HHTo4B_CV_1_C2V_1_C3_2};

    vector<TH1F*> hlist = {heffi_QCD_600toInf, heffi_TT, \
        heffi_VBF_HHTo4B_CV_1_C2V_2_C3_1, heffi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1, \
        heffi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1,heffi_VBF_HHTo4B_CV_1_C2V_1_C3_1,\
        heffi_VBF_HHTo4B_CV_1_C2V_1_C3_0,heffi_VBF_HHTo4B_CV_1_C2V_1_C3_2};

    vector<TH1F*> hcumlist = {hcumeffi_QCD_600toInf, hcumeffi_TT, \
        hcumeffi_VBF_HHTo4B_CV_1_C2V_2_C3_1, hcumeffi_VBF_HHTo4B_CV_0_5_C2V_1_C3_1, \
        hcumeffi_VBF_HHTo4B_CV_1_5_C2V_1_C3_1,hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_1,\
        hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_0,hcumeffi_VBF_HHTo4B_CV_1_C2V_1_C3_2};
    
    vector<string> cutList = {"nTotal","nJet","eta","pt","SDMass","DeepAk8","AK4_matchH","AK4etaCorr","invMass"};
    for (int i=0;i<cutList.size();i++) {
        for (int j=0;j<hlist.size();j++) {
            if (i>0) hlist[j]->Fill(cutList[i].data(),0);
            hcumlist[j]->Fill(cutList[i].data(),0);
        }
    }
    for (int j=0;j<cutList.size();j++) h1->Fill(cutList[j].data(),0); 
    for (int j=0;j<cutList.size();j++) h2->Fill(cutList[j].data(),0); 
    
    c1->Print("hcumEffi_sig.pdf[");
    int colorList[] = {kBlack,kGray+1,kRed,kOrange,kPink-4,kViolet,kMagenta-9,kAzure};
    int markerList[] = {kFullCircle,kFullSquare,kFullTriangleUp,kFullStar,kFullDiamond,kFullCross};
    for (int i=2;i<hlist.size();i++) {
        h1->Reset();
        h2->Reset();
        // fill events
        for (int j=0;j<dataList[i].size()-1;j++) {
            for (int nn=0;nn<dataList[i][j];nn++) h1->Fill(j); 
            for (int nn=0;nn<dataList[i][j+1];nn++) h2->Fill(j); 
        }
        hlist[i]->Divide(h2,h1,1,1,"b");
        hlist[i]->Draw("hist text 0");
        c1->Print("hcumEffi_sig.pdf");


        for (int j=0;j<dataList[i].size();j++) {
            for (int nn=0;nn<dataList[i][j];nn++) hcumlist[i]->Fill(j); 
        }
        hcumlist[i]->Scale(1./dataList[i][0]);
        hcumlist[i]->SetLineColor(colorList[i]);
        hcumlist[i]->SetLineWidth(2);
        hs_cumeffi->Add(hcumlist[i],"");
        if (i>1) hs_sigeffi->Add(hcumlist[i],"");
    }
    TH1F* h_empty = (TH1F*)hlist[0]->Clone("");
    h_empty->Reset();
    auto mg = new TMultiGraph();
    mg->SetTitle("cumulative cut efficiency");
    vector<TGraph*> tgList;
    for (int i=0;i<hlist.size();i++) {
        hcumlist[i]->Draw("L");
        auto tgtemp = new TGraphErrors(hcumlist[i]);
        for (int j=0;j<tgtemp->GetN();j++) tgtemp->SetPointError(j,0,tgtemp->GetErrorY(j)); 
        tgList.push_back(tgtemp);
        //tgList.push_back(new TGraphErrors(hcumlist[i]));
        c1->Print("hcumEffi_sig.pdf");
    }
    for (int i=0;i<tgList.size();i++) mg->Add(tgList[i],"APL");
    for (int i=0;i<hlist.size();i++) {
        hlist[i]->SetFillColor(colorList[i]);
        hlist[i]->SetLineColor(colorList[i]);
        hlist[i]->SetLineWidth(0);
        hlist[i]->Scale(100);
        hs_effi->Add(hlist[i],"hist text 0");
        if (i>1) {
            hlist[i]->SetMarkerStyle(markerList[i-2]);
            hlist[i]->SetMarkerSize(1);
    }
        if(i>1) hssig_effi->Add(hlist[i],"hist text 0");
    }
    hs_effi->Draw("nostackb");
    gPad->BuildLegend(0.85,0.75,0.99,0.95,"","F");
    c1->Print("hcumEffi_sig.pdf");
    hssig_effi->SetMaximum(150);
    hssig_effi->Draw("nostackb");
    auto leg = gPad->BuildLegend(0.7,0.7,0.95,0.95,"","F");
    hssig_effi->Add(h_empty,"hist text 0");
    hssig_effi->Draw("nostackb");
    hssig_effi->GetYaxis()->SetTitle("Efficiency (%)");
    leg->Draw();
    c1->Print("hcumEffi_sig.pdf");


    //gStyle->SetHistLineWidth(2);
    hs_cumeffi->Draw("nostack");
    gPad->BuildLegend(0.7,0.70,0.99,0.95,"","");
    c1->SetLogy();
    c1->Print("hcumEffi_sig.pdf");
    hs_sigeffi->Draw("nostack");
    hs_sigeffi->GetYaxis()->SetTitle("Cumulative Efficiency");
    gPad->BuildLegend(0.67,0.65,0.95,0.9,"","");
    c1->SetLogy();
    c1->Print("hcumEffi_sig.pdf");
    mg->Draw("APL");
    gPad->BuildLegend(0.8,0.75,0.99,0.95,"","");
    c1->SetLogy();
    c1->Print("hcumEffi_sig.pdf");

    c1->Print("hcumEffi_sig.pdf]");
}
