#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>
#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>


#include "AnaModule.h"

AnaModule::AnaModule(const std::string& name): SubsysReco(name), legacyContainer(true)
{}

AnaModule::~AnaModule()
{
  delete pos1;
  delete pos2;
  delete pos3;
  delete posvtx;
  delete mom1;
  delete mom2;
  delete mom3;
  delete momvtx;
  delete rec_mom1;
  delete rec_momvtx;
  delete rec_posvtx;
  delete rec_momtgt;
  delete rec_postgt;

  delete pmom;
  delete nmom;
  delete rec_pmom;
  delete rec_nmom;
  delete rec_ppos;
  delete rec_npos;
  delete rec_vtx;
  delete vtx;
}

int AnaModule::Init(PHCompositeNode* topNode)
{
  pos1 = new TVector3();
  pos2 = new TVector3();
  pos3 = new TVector3();
  posvtx = new TVector3();
  mom1 = new TVector3();
  mom2 = new TVector3();
  mom3 = new TVector3();
  momvtx = new TVector3();
  rec_mom1 = new TVector3();
  rec_momvtx = new TVector3();
  rec_posvtx = new TVector3();
  rec_momtgt = new TVector3();
  rec_postgt = new TVector3();

  vtx  = new TVector3();
  pmom = new TVector3();
  nmom = new TVector3();
  rec_pmom = new TVector3();
  rec_nmom = new TVector3();
  rec_ppos = new TVector3();
  rec_npos = new TVector3();
  rec_vtx  = new TVector3();

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{
  int nTracks = trackVector->size();
  int nRecTracks = legacyContainer ? recEvent->getNTracks() : recTrackVector->size();
  for(int i = 0; i < nTracks; ++i)
  {
    SQTrack* track = trackVector->at(i);
    charge = track->get_charge();
    *pos1 = track->get_pos_st1();
    *mom1 = track->get_mom_st1().Vect();
    *pos3 = track->get_pos_st3();
    *mom3 = track->get_mom_st3().Vect();
    *posvtx = track->get_pos_vtx();
    *momvtx = track->get_mom_vtx().Vect();

    int recid = track->get_rec_track_id();
    if(recid >= 0 && recid < nRecTracks)
    {
      SRecTrack* recTrack = legacyContainer ? &(recEvent->getTrack(recid)) : dynamic_cast<SRecTrack*>(recTrackVector->at(recid));
      *rec_mom1 = recTrack->getMomentumVecSt1();
      *rec_momvtx = recTrack->getVertexMom();
      *rec_posvtx = recTrack->getVertexPos();
      *rec_momtgt = recTrack->getTargetMom();
      *rec_postgt = recTrack->getTargetPos();
    }
    else
    {
      rec_mom1->SetXYZ(-999., -999., -999.);
      rec_momvtx->SetXYZ(-999., -999., -999.);
      rec_posvtx->SetXYZ(-999., -999., -999.);
    }

    saveTree1->Fill();
  }

  int nDimuons = dimuonVector->size();
  int nRecDimuons = legacyContainer ? recEvent->getNDimuons() : (recDimuonVector ? recDimuonVector->size() : -1);
  for(int i = 0; i < nDimuons; ++i)
  {
    SQDimuon* dimuon = dimuonVector->at(i);
    mass = dimuon->get_mom().M();
    *vtx = dimuon->get_pos();
    *pmom = dimuon->get_mom_pos().Vect();
    *nmom = dimuon->get_mom_neg().Vect();
    if(AllChamberPlaneHits(dimuon->get_track_id_pos(), hitVector) && AllChamberPlaneHits(dimuon->get_track_id_neg(), hitVector)) mass_acc =  mass;
    else mass_acc = -999.;

    int recid = dimuon->get_rec_dimuon_id();
    if(recid >= 0 && recid < nRecDimuons)
    {
      SRecDimuon* recDimuon = legacyContainer ? &(recEvent->getDimuon(recid)) : dynamic_cast<SRecDimuon*>(recDimuonVector->at(recid));
      rec_mass = recDimuon->mass;
      *rec_pmom = recDimuon->p_pos.Vect();
      *rec_nmom = recDimuon->p_neg.Vect();
      *rec_ppos = recDimuon->vtx_pos;
      *rec_npos = recDimuon->vtx_neg;
      *rec_vtx  = recDimuon->vtx;
    }
    else
    {
      rec_mass = -100.;
    }

    saveTree2->Fill();
  }

  ++eventID;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveTree1->Write();
  saveTree2->Write();
  saveFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
  hitVector    = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  trackVector  = findNode::getClass<SQTrackVector>(topNode, "SQTruthTrackVector");
  dimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  if(!hitVector || !trackVector || !dimuonVector)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if(legacyContainer)
  {
    recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
    if(!recEvent)
    {
      recEvent = nullptr;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else
  {
    recTrackVector  = findNode::getClass<SQTrackVector>(topNode, "SQRecTrackVector");
    recDimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQRecDimuonVector");
    if(!recTrackVector)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");

  saveTree1 = new TTree("trk", "Track Tree Created by AnaModule");
  saveTree1->Branch("eventID", &eventID, "eventID/I");
  saveTree1->Branch("charge", &charge, "charge/I");
  saveTree1->Branch("pos1", &pos1, 256000, 99);
  saveTree1->Branch("pos2", &pos2, 256000, 99);
  saveTree1->Branch("pos3", &pos3, 256000, 99);
  saveTree1->Branch("posvtx", &posvtx, 256000, 99);
  saveTree1->Branch("mom1", &mom1, 256000, 99);
  saveTree1->Branch("mom2", &mom2, 256000, 99);
  saveTree1->Branch("mom3", &mom3, 256000, 99);
  saveTree1->Branch("momvtx", &momvtx, 256000, 99);
  saveTree1->Branch("rec_mom1", &rec_mom1, 256000, 99);
  saveTree1->Branch("rec_momvtx", &rec_momvtx, 256000, 99);
  saveTree1->Branch("rec_posvtx", &rec_posvtx, 256000, 99);
  saveTree1->Branch("rec_momtgt", &rec_momtgt, 256000, 99);
  saveTree1->Branch("rec_postgt", &rec_postgt, 256000, 99);

  saveTree2 = new TTree("dim", "Dimuon Tree Created by AnaModule");
  saveTree2->Branch("eventID", &eventID, "eventID/I");
  saveTree2->Branch("mass", &mass, "mass/D");
  saveTree2->Branch("mass_acc", &mass_acc, "mass_acc/D");
  saveTree2->Branch("rec_mass", &rec_mass, "rec_mass/D");
  saveTree2->Branch("vtx", &vtx, 256000, 99);
  saveTree2->Branch("pmom", &pmom, 256000, 99);
  saveTree2->Branch("nmom", &nmom, 256000, 99);
  saveTree2->Branch("rec_pmom", &rec_pmom, 256000, 99);
  saveTree2->Branch("rec_nmom", &rec_nmom, 256000, 99);
  saveTree2->Branch("rec_ppos", &rec_ppos, 256000, 99);
  saveTree2->Branch("rec_npos", &rec_npos, 256000, 99);
  saveTree2->Branch("rec_vtx", &rec_vtx, 256000, 99);
}

bool AnaModule::AllChamberPlaneHits(const int trk_id, SQHitVector *hitVector )
{
    bool D0Xhit, D0Xphit, D0Vhit, D0Vphit, D0Uhit, D0Uphit, D2Xhit, D2Xphit, D2Uhit, D2Uphit, D2Vhit, D2Vphit, D3pXhit, D3pXphit, D3pUhit, D3pUphit, D3pVhit, D3pVphit, D3mXhit, D3mXphit, D3mUhit, D3mUphit, D3mVhit, D3mVphit;

    D0Xhit = D0Xphit = D0Vhit =D0Vphit= D0Uhit= D0Uphit= D2Xhit= D2Xphit= D2Uhit= D2Uphit= D2Vhit= D2Vphit= D3pXhit= D3pXphit= D3pUhit= D3pUphit= D3pVhit= D3pVphit= D3mXhit= D3mXphit= D3mUhit= D3mUphit= D3mVhit= D3mVphit = false;


    for(int ihit=0; ihit<hitVector->size(); ++ihit) {
        SQHit *sqhit = hitVector->at(ihit);
        if(sqhit->get_track_id() == trk_id){
            if(sqhit->get_detector_id() == 1 ) D0Uhit=true;
            if(sqhit->get_detector_id() == 2 ) D0Uphit=true;
            if(sqhit->get_detector_id() == 3 ) D0Xphit=true;
            if(sqhit->get_detector_id() == 4 ) D0Xhit=true;
            if(sqhit->get_detector_id() == 5 ) D0Vhit=true;
            if(sqhit->get_detector_id() == 6 ) D0Vphit=true;

            if(sqhit->get_detector_id() == 13) D2Vhit=true;
            if(sqhit->get_detector_id() == 14) D2Vphit=true;
            if(sqhit->get_detector_id() == 15) D2Xphit=true;
            if(sqhit->get_detector_id() == 16) D2Xhit=true;
            if(sqhit->get_detector_id() == 17) D2Uhit=true;
            if(sqhit->get_detector_id() == 18) D2Uphit=true;

            if(sqhit->get_detector_id() == 19) D3pVphit=true;
            if(sqhit->get_detector_id() == 20) D3pVhit=true;
            if(sqhit->get_detector_id() == 21) D3pXphit=true;
            if(sqhit->get_detector_id() == 22) D3pXhit=true;
            if(sqhit->get_detector_id() == 23) D3pUphit=true;
            if(sqhit->get_detector_id() == 24) D3pUhit=true;

            if(sqhit->get_detector_id() == 25) D3mVphit=true;
            if(sqhit->get_detector_id() == 26) D3mVhit=true;
            if(sqhit->get_detector_id() == 27) D3mXphit=true;
            if(sqhit->get_detector_id() == 28) D3mXhit=true;
            if(sqhit->get_detector_id() == 29) D3mUphit=true;
            if(sqhit->get_detector_id() == 30) D3mUhit=true;
        }
    }

    st1hit = false;
    st2hit = false;
    st3phit = false;
    st3mhit = false;
    st3hit = false;

    st1hit = D0Uhit && D0Uphit && D0Xphit && D0Xhit && D0Vhit && D0Vphit;
    st2hit = D2Uhit && D2Uphit && D2Xphit && D2Xhit && D2Vhit && D2Vphit;
    st3phit = D3pUhit && D3pUphit && D3pXphit && D3pXhit && D3pVhit && D3pVphit;
    st3mhit = D3mUhit && D3mUphit && D3mXphit && D3mXhit && D3mVhit && D3mVphit;
    st3hit = st3phit || st3mhit;

    //if (st3phit && st3mhit) return false;
    if(st1hit && st2hit && st3hit) return true;
    else return false;
}
