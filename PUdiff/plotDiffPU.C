#include "setNCUStyle.C"
using namespace ROOT;
#define runShort 0
void setMax(TH1D* h1, TH1D* h2,double scale = 1.2) {
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    if (max1>max2) {
        h2->SetMaximum(max1*scale);
        h1->SetMaximum(max1*scale);
    }
    else {
        h1->SetMaximum(max2*scale);
        h2->SetMaximum(max2*scale);
    }
}
void plotDiffPU() {
    //setNCUStyle(1);
    TH1::SetDefaultSumw2(1);
    string parList[] = {"CV_0_5_C2V_1_C3_1","CV_1_C2V_1_C3_0","CV_1_C2V_1_C3_2","CV_1_5_C2V_1_C3_1","CV_1_C2V_1_C3_1","CV_1_C2V_2_C3_1"};
    string puList[] = {"PU200","NoPU"};
    //vector<std::string_view > varList = {"lead_fatjet.Pt()","sublead_fatjet.Pt()","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1"};
    //  AK4 information
    vector<std::string_view > varList = {"lead_fatjet_Pt","sublead_fatjet_Pt","lead_fatjet_sdmass","sublead_fatjet_sdmass","HHinvMass","tau21_0","tau21_1","DeepAK8_jet0","DeepAK8_jet1","vbfjet_invmass","vbfjet0_Pt","vbfjet1_Pt","vbfjet0_E","vbfjet1_E","vbfjet0_Eta","vbfjet1_Eta","vbfjet_deltaEta"};
    int spCase[] = {0,3,4};
    vector<vector<TFile*>> fList(2); 
    //TFile* fList[6][2];
    vector<TChain*> pu200List;
    vector<TChain*> nopuList;
    for (int i=0;i<6;i++) {
        pu200List.push_back(new TChain("mytree",parList[i].data())); 
        nopuList.push_back(new TChain("mytree",parList[i].data())); 
    }
    vector<vector<TH1D*>> histPU200List(6);
    vector<vector<TH1D*>> histNoPUList(6);
    vector<vector<TH2D*>> histPUtoSDMass(2);
    vector<vector<TH2D*>> histPVtoSDMass(2);
    for (int i=0;i<6;i++) {
        if (i!=0 && runShort) continue;
        pu200List[i]->Add(Form("mytree_VBF_HHTo4B_%s_nTuple_%s.root",parList[i].data(),puList[0].data()));
        if (i==0 || i==3 || i==4) {
            nopuList[i]->Add(Form("mytree_VBF_HHTo4B_%s_nTuple_%s_1.root",parList[i].data(),puList[1].data()));
            nopuList[i]->Add(Form("mytree_VBF_HHTo4B_%s_nTuple_%s_2.root",parList[i].data(),puList[1].data()));
        }
        else nopuList[i]->Add(Form("mytree_VBF_HHTo4B_%s_nTuple_%s.root",parList[i].data(),puList[1].data()));
        cout << "PU200 " << parList[i] << "\t" << pu200List[i]->GetEntries() << endl;
        cout << "NoPU  " << parList[i] << "\t" << nopuList[i]->GetEntries() << endl;
        
        //if (pu200List[i]->GetEntries() < 100) continue;
        
        RDataFrame d0(*pu200List[i]);
        RDataFrame s0(*nopuList[i]);
        auto d = d0.Define("lead_fatjet_Pt","lead_fatjet.Pt()").Define("sublead_fatjet_Pt","sublead_fatjet.Pt()")
            .Define("vbfjet0_Pt","vbf_jet0.Pt()").Define("vbfjet0_E","vbf_jet0.E()").Define("vbfjet0_Eta","vbf_jet0.Eta()")
            .Define("vbfjet1_Pt","vbf_jet1.Pt()").Define("vbfjet1_E","vbf_jet1.E()").Define("vbfjet1_Eta","vbf_jet1.Eta()")
            .Define("vbfjet_deltaEta","abs(vbf_jet0.Eta()-vbf_jet1.Eta())");

        auto s = s0.Define("lead_fatjet_Pt","lead_fatjet.Pt()").Define("sublead_fatjet_Pt","sublead_fatjet.Pt()")
            .Define("vbfjet0_Pt","vbf_jet0.Pt()").Define("vbfjet0_E","vbf_jet0.E()").Define("vbfjet0_Eta","vbf_jet0.Eta()")
            .Define("vbfjet1_Pt","vbf_jet1.Pt()").Define("vbfjet1_E","vbf_jet1.E()").Define("vbfjet1_Eta","vbf_jet1.Eta()")
            .Define("vbfjet_deltaEta","abs(vbf_jet0.Eta()-vbf_jet1.Eta())");
        for (int j=0;j<varList.size();j++) histPU200List[i].push_back((TH1D*)d.Histo1D(varList[j].data())->Clone(Form("PU200_%s_%s",parList[i].data(),varList[j].data())));
        for (int j=0;j<varList.size();j++) histNoPUList[i].push_back((TH1D*)s.Histo1D(varList[j].data())->Clone(Form("NoPU_%s_%s",parList[i].data(),varList[j].data())));
        string PUM2dName = Form("histPUtoSDMass_%s",parList[i].data());
        string PUM2dName2 = Form("histPUtoSDMass2_%s",parList[i].data());
        string PVM2dName = Form("histPVtoSDMass_%s",parList[i].data());
        string PVM2dName2 = Form("histPVtoSDMass2_%s",parList[i].data());
        auto temp2d_1 = d.Histo2D({PUM2dName.data(),PUM2dName.data(),10,150,250,10,90,140},"nPU","lead_fatjet_sdmass");
        auto temp2d_2 = d.Histo2D({PUM2dName2.data(),PUM2dName2.data(),10,150,250,10,90,140},"nPU","sublead_fatjet_sdmass");
        auto temp2d_3 = d.Histo2D({PVM2dName.data(),PVM2dName.data(),10,100,200,10,90,140},"nPV","lead_fatjet_sdmass");
        auto temp2d_4 = d.Histo2D({PVM2dName2.data(),PVM2dName2.data(),10,100,200,10,90,140},"nPV","sublead_fatjet_sdmass");
        //temp2d_1->Draw("colz");
        //temp2d_2->Draw("colz");
        histPUtoSDMass[0].push_back((TH2D*)temp2d_1->Clone(Form("PUtoFatJet1SDMass_%s",parList[i].data())));
        histPUtoSDMass[1].push_back((TH2D*)temp2d_2->Clone(Form("PUtoFatJet2SDMass_%s",parList[i].data())));
        histPVtoSDMass[0].push_back((TH2D*)temp2d_3->Clone(Form("PVtoFatJet1SDMass_%s",parList[i].data())));
        histPVtoSDMass[1].push_back((TH2D*)temp2d_4->Clone(Form("PVtoFatJet2SDMass_%s",parList[i].data())));
    }
    cout << "finish loading" << endl;

    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.03, "XYZ");
    
    TCanvas *c1 = new TCanvas("c1","",1800,1200);
    c1->Divide(3,2);
    string pdfName = "diffPU.pdf";
    c1->Print((pdfName+"[").data());
    TLegend* leg = 0;
    for (int j=0;j<varList.size();j++) {
        for (int i=0;i<6;i++) {
            if (i!=0 && runShort) continue;
            if (!leg) {
                leg = new TLegend(0.7,0.7,0.9,0.9);
                leg->AddEntry(histNoPUList[i][j],"No PU");
                leg->AddEntry(histPU200List[i][j],"PU200");
            }
            c1->cd(i+1);
            gPad->SetLeftMargin(0.15);
            histNoPUList[i][j]->SetLineColor(kBlack);
            histPU200List[i][j]->SetLineColor(kBlue);
            histNoPUList[i][j]->SetTitle(parList[i].data());
            histPU200List[i][j]->SetTitle(parList[i].data());
            histNoPUList[i][j]->GetYaxis()->SetTitle("A.U.");
            histNoPUList[i][j]->Scale(1./histNoPUList[i][j]->Integral());
            histPU200List[i][j]->Scale(1./histPU200List[i][j]->Integral());
            histNoPUList[i][j]->Rebin(2);
            histPU200List[i][j]->Rebin(2);
            setMax(histNoPUList[i][j],histPU200List[i][j]);
            histNoPUList[i][j]->Draw("hist");
            histPU200List[i][j]->Draw("histsame");
            //histNoPUList[i][j]->DrawNormalized("hist");
            //histPU200List[i][j]->DrawNormalized("histsame");
            //static auto leg = c1->BuildLegend();
            leg->Draw();
        }
        c1->Print(pdfName.data());
    }
    /*
    histPUtoSDMass[0]->Draw("colz");
    c1->Print(pdfName.data());
    histPVtoSDMass[0]->Draw("colz");
    c1->Print(pdfName.data());
    */
    cout << "2D plot" << endl;
    int pind[] = {0,3,5};
    for (int i=0;i<6;i++)  {c1->cd(i+1);histPUtoSDMass[0][i]->SetTitle(Form("%s;nPU;AK8jet1_SDMass",parList[i].data()));histPUtoSDMass[0][i]->Draw("colz"); }
    c1->Print(pdfName.data());
    
    for (int i=0;i<3;i++) {
        c1->cd(i+1);
        histPUtoSDMass[0][pind[i]]->ProfileX(((string)histPUtoSDMass[0][pind[i]]->GetName()+"_ProfileX").data())->Draw();
        c1->cd(i+4);
        histPUtoSDMass[0][pind[i]]->ProfileY(((string)histPUtoSDMass[0][pind[i]]->GetName()+"_ProfileY").data())->Draw();
    }
    c1->Print(pdfName.data());
    
    for (int i=0;i<6;i++)  {c1->cd(i+1);histPUtoSDMass[1][i]->SetTitle(Form("%s;nPU;AK8jet2_SDMass",parList[i].data()));histPUtoSDMass[1][i]->Draw("colz"); }
    c1->Print(pdfName.data());
    
    for (int i=0;i<3;i++) {
        c1->cd(i+1);
        histPUtoSDMass[1][pind[i]]->ProfileX(((string)histPUtoSDMass[1][pind[i]]->GetName()+"_ProfileX").data())->Draw();
        c1->cd(i+4);
        histPUtoSDMass[1][pind[i]]->ProfileY(((string)histPUtoSDMass[1][pind[i]]->GetName()+"_ProfileY").data())->Draw();
    }
    c1->Print(pdfName.data());
    
    for (int i=0;i<6;i++)  {c1->cd(i+1);histPVtoSDMass[0][i]->SetTitle(Form("%s;dof_PV;AK8jet1_SDMass",parList[i].data()));histPVtoSDMass[0][i]->Draw("colz"); }
    c1->Print(pdfName.data());
    
    for (int i=0;i<3;i++) {
        c1->cd(i+1);
        histPVtoSDMass[0][pind[i]]->ProfileX(((string)histPVtoSDMass[0][pind[i]]->GetName()+"_ProfileX").data())->Draw();
        c1->cd(i+4);
        histPVtoSDMass[0][pind[i]]->ProfileY(((string)histPVtoSDMass[0][pind[i]]->GetName()+"_ProfileY").data())->Draw();
    }
    c1->Print(pdfName.data());
    
    for (int i=0;i<6;i++)  {c1->cd(i+1);histPVtoSDMass[1][i]->SetTitle(Form("%s;dof_PV;AK8jet2_SDMass",parList[i].data()));histPVtoSDMass[1][i]->Draw("colz"); }
    c1->Print(pdfName.data());
    for (int i=0;i<3;i++) {
        c1->cd(i+1);
        histPVtoSDMass[1][pind[i]]->ProfileX(((string)histPVtoSDMass[1][pind[i]]->GetName()+"_ProfileX").data())->Draw();
        c1->cd(i+4);
        histPVtoSDMass[1][pind[i]]->ProfileY(((string)histPVtoSDMass[1][pind[i]]->GetName()+"_ProfileY").data())->Draw();
    }
    c1->Print(pdfName.data());
    cout << "profile" << endl;
        
    c1->Print((pdfName+"]").data());
}
