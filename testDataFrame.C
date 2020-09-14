//#include <Math/GenVector/LorentzVector.h>
//#include <Math/GenVector/PtEtaPhiM4D.h>
//#include "Math/LorentzVector.h"
#include <TStopwatch.h>
#include "../myCut.cc"
#include <ROOT/RDataFrame.hxx>
#include "Math/Vector4D.h"

#define csvv2_M_102x  0.4184
using namespace ROOT;
bool doTest = 0;
void testDataFrame_base(RDataFrame& );
void drawhist(RDataFrame& );

void testDataFrame(int nFile=1) {
    TStopwatch timer;
    timer.Start();
    RDataFrame d(0);
    const int nLen = 1000;
    std::vector< std::string > list;
    if (doTest) {
         //d = RDataFrame("Events","root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/VBFHHTo4B_CV_1_C2V_1_C3_2_TuneCP5_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/D637D126-0B17-2446-A53A-3A9908A6FACF.root");
         d = RDataFrame("Events","../VBFHHto4B_nanov6.root ");
    }
    else {
         // dasgoclient --query="file dataset=/VBFHHTo4B_CV_1_C2V_2_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"
         ROOT::EnableImplicitMT();
         string prefix = "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/VBFHHTo4B_CV_1_C2V_2_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/";
         string bli[] = {"70000/FBADAE20-C7D3-F142-BEB4-EBFB7F61DBD9.root","260000/3735B467-315C-084F-9875-855850D606F9.root","100000/A1A504A6-45C9-9445-A08E-66292B57FDE8.root","100000/AE3ECB4E-189F-674A-9AF7-4FE23C7169E3.root","100000/7A8C9579-7129-324E-881C-3DA492BF333F.root","130000/88DBC34D-EB41-334E-A495-EE2C6E9E0F4F.root","110000/3E554CB0-0899-164B-86CF-8F2F4E330132.root","110000/9FCEC660-6953-894B-B7FC-25644AF2F0BB.root","130000/7075D526-BD7F-EF45-94AD-FACCBDFF5425.root","130000/8FC99EA2-5883-6240-8120-B164112F7447.root","70000/CEA5FEA8-43B3-004D-AC08-04AF053E11E6.root","130000/DEC02E22-8155-3F48-909A-5DEFE4B745A8.root"};
         // string prefix = "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv7/VBFHHTo4B_CV_1_C2V_1_C3_2_TuneCP5_PSWeights_13TeV-madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/";
         // string bli[13] = {"70000/D637D126-0B17-2446-A53A-3A9908A6FACF.root","70000/A10CD057-C41D-7641-897A-D31C7C3745ED.root","70000/840CEE7A-A350-A143-A2A9-7233124CF861.root","70000/427D1D0C-7FA3-4F4B-9068-E9D709223E5B.root","70000/2B95A36F-2E39-E74D-92C2-103D1D601C24.root","70000/1BBDA349-36CF-D24D-902A-807729453ED6.root","60000/DFB9FBCF-8858-F54D-9869-86BFB8F4BA63.root","60000/883A9438-CA06-674A-9B17-FBC71E4ACA29.root","60000/6DA6A77D-898A-4846-B6DA-D52FBD579673.root","270000/9310FBD7-EBC0-8047-BC36-CDF0D5BAC972.root","130000/726C033E-E402-EF46-B21D-64B549532EB5.root","120000/071CA18E-C0F3-EB45-8102-77AE6053C26A.root","100000/154922D2-ADD4-6946-BD32-5455BA08B872.root"};
         for (int i=0;i<nFile;i++) list.push_back(prefix+bli[i]);
         d = RDataFrame("Events",list);
         //return;
    }
    //drawhist(d);
    testDataFrame_base(d);
    timer.Stop();
     cout << "CpuTime\tRealTime" << endl;
     cout << timer.CpuTime()  << "\t" << timer.RealTime() << endl;
}


void testDataFrame_base(RDataFrame& d) {
     auto dd = d.Define("fatjet0_pt","FatJet_pt[0]").Define("fatjet1_pt","FatJet_pt[1]");
     auto g = dd.Filter("nFatJet>=2 && FatJet_pt[0] > 300 && FatJet_pt[1] > 300 && abs(FatJet_eta[0]) < 3.0 && abs(FatJet_eta[1]) < 3.0","FatJet pt&eta cut");
     auto gg = g.Filter("FatJet_tau2[0]/FatJet_tau1[0] < 0.6 && FatJet_tau2[1]/FatJet_tau1[1] < 0.6 && FatJet_msoftdrop[0] > 90. && FatJet_msoftdrop[0] < 140. && FatJet_msoftdrop[1] > 90. && FatJet_msoftdrop[1] < 140.","FatJet tau21 and sdmass cut");
     auto ggg = gg.Filter("nSubJet>=2").Define("nbsubjet1","(FatJet_subJetIdx1[0]>=0&&SubJet_btagCSVV2[FatJet_subJetIdx1[0]]>csvv2_M_102x)+(FatJet_subJetIdx2[0]>=0&&SubJet_btagCSVV2[FatJet_subJetIdx2[0]]>csvv2_M_102x)");
     auto gggg = ggg.Define("nbsubjet2","(FatJet_subJetIdx1[1]>=0&&SubJet_btagCSVV2[FatJet_subJetIdx1[1]]>csvv2_M_102x)+(FatJet_subJetIdx2[1]>=0&&SubJet_btagCSVV2[FatJet_subJetIdx2[1]]>csvv2_M_102x)");
     auto ggggg = gggg.Filter("nbsubjet1>=1&&nbsubjet2>=1","subjet cut");
     auto b = ggggg.Define("lead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJet_pt[0],FatJet_eta[0],FatJet_phi[0],FatJet_mass[0])");
     auto bb = b.Define("sublead_fatjet","ROOT::Math::PtEtaPhiMVector(FatJet_pt[1],FatJet_eta[1],FatJet_phi[1],FatJet_mass[1])");
     auto bbb = bb.Define("ind1","doSelection(nJet,Jet_pt,Jet_eta,Jet_phi,Jet_mass,lead_fatjet,sublead_fatjet)[0]");
     auto bbbb = bbb.Define("ind2","doSelection(nJet,Jet_pt,Jet_eta,Jet_phi,Jet_mass,lead_fatjet,sublead_fatjet)[1]");
     auto finalNode = bbbb.Filter("ind1>=0&&ind2>=0","vbf cut").Define("lead_fatjet_sdmass","FatJet_msoftdrop[0]").Define("sublead_fatjet_sdmass","FatJet_msoftdrop[1]");
     cout << "running" << endl;
     initializer_list< std::string > outputlist = {"lead_fatjet","sublead_fatjet","lead_fatjet_sdmass","sublead_fatjet_sdmass","nJet","Jet_pt","Jet_eta","Jet_phi","Jet_mass","ind1","ind2"};
     finalNode.Snapshot("mytree","database.root",outputlist);
     finalNode.Report()->Print();
     
}
void drawhist(RDataFrame& d) {
     auto dd = d.Filter("nFatJet>=2").Define("fatjet0_pt","FatJet_pt[0]").Define("fatjet1_pt","FatJet_pt[1]");
     auto h_leadjet_pt = dd.Histo1D("fatjet0_pt");
     auto h_nextleadjet_pt = dd.Histo1D("fatjet1_pt");
     h_leadjet_pt->DrawCopy("hist");
     h_nextleadjet_pt->DrawCopy("hist");
     TFile *hout = new TFile("histout.root","recreate");
     h_leadjet_pt->Write();
     h_nextleadjet_pt->Write();
     hout->Close();
}
