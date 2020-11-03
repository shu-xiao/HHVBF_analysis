#include "setNCUStyle.C"
using namespace ROOT;
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
    for (int i=0;i<6;i++) {
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
    c1->Print((pdfName+"]").data());
}
