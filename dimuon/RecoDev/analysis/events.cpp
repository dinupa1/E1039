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
    int n = tree2->GetEntries();

    int eventID;
    double mass, rec_mass, mass_acc;

    tree2->SetBranchAddress("eventID", &eventID);
    tree2->SetBranchAddress("mass_acc", &mass_acc);
    tree2->SetBranchAddress("mass", &mass);

    for(int i = 0; i < n; i++){
        tree2->GetEntry(i);
        if(mass > 5.0 && mass_acc < 0.0){
            cout << "event ID : " << eventID << endl;
        }
    }

    return 0;
}
