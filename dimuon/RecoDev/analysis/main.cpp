//
// Created by dinupa on 6/18/21.
//

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TEfficiency.h>
#include <iostream>

using namespace std;

// combine mutiple files
int com_file(){
    TChain* chain1 = new TChain("dim");
    TChain* chain2 = new TChain("trk");

    chain1->AddFile("data2.root");
    chain1->AddFile("data3.root");

    chain2->AddFile("data2.root");
    chain2->AddFile("data3.root");

    TFile* file = TFile::Open("data.root", "RECREATE");
    chain1->CloneTree(-1, "fast");
    chain2->CloneTree(-1, "fast");
    file -> Write();

    return 0;
}

// plots at vertex
int plots_vtx(){
    TFile* file = TFile::Open("data.root");
    auto tree = (TTree*)file->Get("dim");
    int n = tree->GetEntries();

    auto vtx = new TVector3(0.0, 0.0, 0.0);
    auto pmom = new TVector3(0.0, 0.0, 0.0);
    auto nmom = new TVector3(0.0, 0.0, 0.0);

    tree->SetBranchAddress("vtx", &vtx);
    tree->SetBranchAddress("pmom", &pmom);
    tree->SetBranchAddress("nmom", &nmom);

    auto hist1 = new TH1F("hist1", "x at vtx; x (cm); counts", 50, -1.0, +1.0);
    auto hist2 = new TH1F("hist2", "y at vtx; y (cm); counts", 50, -1.0, +1.0);
    auto hist3 = new TH1F("hist3", "z at vtx; z (cm); counts", 50, -305.0, -295.0);
    auto hist4 = new TH1F("hist4", "px at vtx; px (GeV/c); counts", 50, -5.0, 5.0);
    auto hist5 = new TH1F("hist5", "py at vtx; py (GeV/c); counts", 50, -5.0, 5.0);
    auto hist6 = new TH1F("hist6", "pz at vtx; pz (GeV/c); counts", 50, 20.0, 140.0);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        hist1->Fill(vtx->X());
        hist2->Fill(vtx->Y());
        hist3->Fill(vtx->Z());
        hist4->Fill(pmom->Px()+nmom->Px());
        hist5->Fill(pmom->Py()+nmom->Py());
        hist6->Fill(pmom->Pz()+nmom->Pz());
    }

    auto can1 = new TCanvas();can1->SetGrid();
    hist1->Draw();
    can1->SaveAs("pic1.png");

    auto can2 = new TCanvas();can2->SetGrid();
    hist2->Draw();
    can2->SaveAs("pic2.png");

    auto can3 = new TCanvas();can3->SetGrid();
    hist3->Draw();
    can3->SaveAs("pic3.png");

    auto can4 = new TCanvas();can4->SetGrid();
    hist4->Draw();
    can4->SaveAs("pic4.png");

    auto can5 = new TCanvas();can5->SetGrid();
    hist5->Draw();
    can5->SaveAs("pic5.png");

    auto can6 = new TCanvas();can6->SetGrid();
    hist6->Draw();
    can6->SaveAs("pic6.png");

    return 0;
}

// efficiency plots
int effi_plots(){
    TFile* file = TFile::Open("data.root");
    auto tree = (TTree*)file->Get("dim");
    int n = tree->GetEntries();

    double mass, rec_mass, mass_acc;

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("rec_mass", &rec_mass);
    tree->SetBranchAddress("mass_acc", &mass_acc);

    auto hist1 = new TH1F("hist1", "true mass; mass (GeV/c^{2}); counts", 20, 3.0, 8.0);
    hist1->SetLineColor(kRed);
    auto hist2 = new TH1F("hist2", "true mass (rec_mass > 0.0); mass (GeV/c^{2}); counts", 20, 3.0, 8.0);
    hist2->SetLineColor(kBlue);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        if(mass_acc > 0.0){
            hist1->Fill(mass);
            if(rec_mass > 0.0){
                hist2->Fill(mass);
            }
        }
    }

    auto effi = new TEfficiency(*hist2, *hist1);
    effi->SetTitle("mass vs. efficiency; mass (GeV/c^{2}); efficiency");
    effi->SetMarkerStyle(21);
    effi->SetMarkerColor(2);

    auto can1 = new TCanvas();can1->SetGrid();
    effi->Draw("APE1");
    can1->SaveAs("pic7.png");

    auto can2 = new TCanvas();can2->SetGrid();
    hist1->Draw();
    hist2->Draw("same");
    //gPad->BuildLegend(0.5,0.5,0.7,0.7,"");
    can2->SaveAs("pic8.png");

    return 0;
}

// mass distributions with different variables
int plot_mass(){
    TFile* file = TFile::Open("data.root");
    auto tree = (TTree*)file->Get("dim");
    int n = tree->GetEntries();

    double mass, mass_acc, rec_mass;
    auto pmom = new TVector3 (0.0, 0.0, 0.0);
    auto nmom = new TVector3 (0.0, 0.0, 0.0);

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("mass_acc", &mass_acc);
    tree->SetBranchAddress("rec_mass", &rec_mass);
    tree->SetBranchAddress("pmom", &pmom);
    tree->SetBranchAddress("nmom", &nmom);

    auto hist1 = new TH2F("hist1",
                          "acc. mass vs. momentum; acc. mass(GeV/c^{2}); momentum (GeV/c)",
                          20, 3.0, 9.0, 50, 20.0, 140.0);

    auto hist2 = new TH2F("hist2",
                          "acc. mass vs. rec. mass; acc. mass(GeV/c^{2}); rec. mass(GeV/c^{2})",
                          20, 1.0, 10.0, 20, 1.0, 10.0);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        if(mass_acc > 0.0){
            hist1->Fill(mass_acc, (*pmom+*nmom).Mag());
            if(rec_mass > 0.0){
                hist2->Fill(mass_acc, rec_mass);
            }
        }
    }

    auto can1 = new TCanvas();can1->SetGrid();
    hist1->Draw("COLZ TEXT");
    can1->SaveAs("pic9.png");

    auto can2 = new TCanvas();can2->SetGrid();
    hist2->Draw("COLZ TEXT");
    can2->SaveAs("pic10.png");

    return 0;
}

// properties in higher mass region
int high_mass(){
    TFile* file = TFile::Open("data.root");
    auto tree = (TTree*)file->Get("dim");
    int n = tree->GetEntries();

    double mass, mass_acc, rec_mass;
    auto pmom = new TVector3(0.0, 0.0, 0.0);
    auto nmom = new TVector3(0.0, 0.0, 0.0);

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("mass_acc", &mass_acc);
    tree->SetBranchAddress("rec_mass", &rec_mass);
    tree->SetBranchAddress("pmom", &pmom);
    tree->SetBranchAddress("nmom", &nmom);

    auto hist1 = new TH2F("hist1",
                          "acc. mass vs. momentum; acc. mass(GeV/c^{2}); mom (GeV/c)",
                          20, 7.0, 10.0, 50, 20.0, 140.0);

    auto hist2 = new TH2F("hist2",
                          "rec. mass vs. momentum; rec. mass(GeV/c^{2}); mom (GeV/c)",
                          20, 7.0, 10.0, 50, 20.0, 140.0);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        if(mass > 7.0){
            hist1->Fill(mass_acc, (*pmom+*nmom).Mag());
            if(rec_mass > 0.0){
                hist2->Fill(rec_mass, (*pmom+*nmom).Mag());
            }
        }
    }

    auto can1 = new TCanvas();can1->SetGrid();
    hist1->Draw("COLZ TEXT");
    can1->SaveAs("pic11.png");

    auto can2 = new TCanvas();can2->SetGrid();
    hist2->Draw("COLZ TEXT");
    can2->SaveAs("pic12.png");

    return 0;
}

// single muon recosntruction
int single_reco(){
    TFile* file = TFile::Open("data.root");
    auto tree = (TTree*)file->Get("trk");
    int n = tree->GetEntries();

    auto momvtx = new TVector3(0.0, 0.0, 0.0);
    auto rec_momvtx = new TVector3(0.0, 0.0, 0.0);

    tree->SetBranchAddress("momvtx", &momvtx);
    tree->SetBranchAddress("rec_momvtx", &rec_momvtx);

    auto hist1 = new TH1F("hist1",
                          "true momentum; p (GeV/c); counts",
                          50, 0.0, 100.0);

    auto hist2 = new TH1F("hist2",
                          "rec. momentum; p (GeV/c); counts",
                          50, 0.0, 100.0);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        hist1->Fill(momvtx->Mag());
        if(rec_momvtx->Px() > -999.){
            hist2->Fill(momvtx->Mag());
        }
    }

    auto effi = new TEfficiency(*hist2, *hist1);
    effi->SetTitle("momentum vs. efficiency; p (GeV/c); efficiency");
    effi->SetMarkerStyle(21);
    effi->SetMarkerColor(2);

    auto can1 = new TCanvas();can1->SetGrid();
    effi->Draw("APE1");
    can1->SaveAs("pic13.png");

    return 0;
}

int main(){

    // combine root files
    //com_file();

    //plots at vertex
    //plots_vtx();

    //plots of efficiency
    effi_plots();

    // mass distributions
    //plot_mass();

    // high mass ditributions
    //high_mass();

    // single track efficiency
    //single_reco();

    return 0;
}