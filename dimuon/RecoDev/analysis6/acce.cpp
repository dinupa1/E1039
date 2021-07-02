//
// Created by dinupa on 6/25/21.
//
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include <iostream>

using namespace std;

int acce(){
    TFile* file = TFile::Open("data4.root");
    auto tree = (TTree*)file->Get("dim"); // load dimuon tree
    int n = tree->GetEntries();

    double mass, mass_acc, rec_mass;

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("mass_acc", &mass_acc);
    tree->SetBranchAddress("rec_mass", &rec_mass);

    auto hist1 = new TH1F("hist1", "generated mass; mass(GeV/c^{2}); counts", 50, 0.0 , 11.0);
    auto hist2 = new TH1F("hist2", "accepted mass; mass(GeV/c^{2}); counts", 50, 0.0 , 11.0);

    for(int i = 0; i < n; i++){
        tree->GetEntry(i);
        hist1->Fill(mass);
        if(mass_acc > 0.0){
            hist2->Fill(mass);
        }
    }

    auto effi = new TEfficiency(*hist2, *hist1);
    effi->SetTitle("detector acceptance; mass(GeV/c^{2}); counts");
    effi->SetMarkerStyle(21);
    effi->SetMarkerColor(4);

    auto can1 = new TCanvas(); can1->SetGrid();
    hist1->Draw();
    can1->SaveAs("pic1.png");

    auto can2 = new TCanvas(); can2->SetGrid();
    hist2->Draw();
    can2->SaveAs("pic2.png");

    auto can3 = new TCanvas(); can3->SetGrid();
    effi->Draw("APE1");
    can3->SaveAs("pic3.png");

    return 0;
}
