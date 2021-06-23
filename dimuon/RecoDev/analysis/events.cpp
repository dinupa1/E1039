//
// Created by dinupa on 6/22/21.
//

#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <iostream>

using namespace std;

int events(){
    TFile* file = TFile::Open("data.root");
    auto tree1 = (TTree*)file->Get("trk"); // tree with track info
    auto tree2 = (TTree*)file->Get("dim"); // tree with dimuon info
    tree2->BuildIndex("eventID");
    tree1->AddFriend(tree2);
    int n = tree1->GetEntries();

    int eventID;
    double mass, rec_mass, mass_acc;
    auto pos1 = new TVector3(0.0, 0.0, 0.0);

    // branch info dimuon tree
    tree1->SetBranchAddress("eventID", &eventID);
    tree1->SetBranchAddress("mass_acc", &mass_acc);
    tree1->SetBranchAddress("mass", &mass);

    // branch info track tree
    tree1->SetBranchAddress("pos1", &pos1);

    auto hist1 = new TH1F("hist1",
                          "x postion station1",
                          50, -60.0, 60.0);

    for(int i = 0; i < n; i++){
        tree1->GetEntry(i);
        if(mass > 5.0 && mass_acc < 0.0){
            hist1->Fill(pos1->X());
        }
    }

    auto can1 = new TCanvas(); can1->SetGrid();
    hist1->Draw();
    can1->SaveAs("pic21.png");

    return 0;
}
