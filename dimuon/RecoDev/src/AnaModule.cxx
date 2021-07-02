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
#include <interface_main/SQMCEvent_v1.h>


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
  delete acc_mom;
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

  delete mom_D0U;
  delete pos_D0U;
  delete mom_D0Up;
  delete pos_D0Up;
  delete mom_D0Xp;
  delete pos_D0Xp;
  delete mom_D0X;
  delete pos_D0X;
  delete mom_D0V;
  delete pos_D0V;
  delete mom_D0Vp;
  delete pos_D0Vp;

  delete mom_D2U;
  delete pos_D2U;
  delete mom_D2Up;
  delete pos_D2Up;
  delete mom_D2Xp;
  delete pos_D2Xp;
  delete mom_D2X;
  delete pos_D2X;
  delete mom_D2V;
  delete pos_D2V;
  delete mom_D2Vp;
  delete pos_D2Vp;

  delete mom_D3U;
  delete pos_D3U;
  delete mom_D3Up;
  delete pos_D3Up;
  delete mom_D3Xp;
  delete pos_D3Xp;
  delete mom_D3X;
  delete pos_D3X;
  delete mom_D3V;
  delete pos_D3V;
  delete mom_D3Vp;
  delete pos_D3Vp;

  delete pos_H1T;
  delete pos_H1B;
  delete pos_H1R;
  delete pos_H1L;

  delete pos_H2T;
  delete pos_H2B;
  delete pos_H2R;
  delete pos_H2L;

  delete pos_H3T;
  delete pos_H3B;
  delete pos_H4B;
  delete pos_H4T;

  delete pos_H4Y1L;
  delete pos_H4Y1R;
  delete pos_H4Y2L;
  delete pos_H4Y2R;

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
  acc_mom = new TVector3();
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

  mom_D0U = new TVector3();
  pos_D0U = new TVector3();
  mom_D0Up = new TVector3();
  pos_D0Up = new TVector3();
  mom_D0Xp = new TVector3();
  pos_D0Xp = new TVector3();
  mom_D0X = new TVector3();
  pos_D0X = new TVector3();
  mom_D0V = new TVector3();
  pos_D0V = new TVector3();
  mom_D0Vp = new TVector3();
  pos_D0Vp = new TVector3();

  mom_D2U = new TVector3();
  pos_D2U = new TVector3();
  mom_D2Up = new TVector3();
  pos_D2Up = new TVector3();
  mom_D2Xp = new TVector3();
  pos_D2Xp = new TVector3();
  mom_D2X = new TVector3();
  pos_D2X = new TVector3();
  mom_D2V = new TVector3();
  pos_D2V = new TVector3();
  mom_D2Vp = new TVector3();
  pos_D2Vp = new TVector3();

  mom_D3U = new TVector3();
  pos_D3U = new TVector3();
  mom_D3Up = new TVector3();
  pos_D3Up = new TVector3();
  mom_D3Xp = new TVector3();
  pos_D3Xp = new TVector3();
  mom_D3X = new TVector3();
  pos_D3X = new TVector3();
  mom_D3V = new TVector3();
  pos_D3V = new TVector3();
  mom_D3Vp = new TVector3();
  pos_D3Vp = new TVector3();

  pos_H1T = new TVector3();
  pos_H1B = new TVector3();
  pos_H1R = new TVector3();
  pos_H1L = new TVector3();

  pos_H2T = new TVector3();
  pos_H2B = new TVector3();
  pos_H2R = new TVector3();
  pos_H2L = new TVector3();

  pos_H3T = new TVector3();
  pos_H3B = new TVector3();
  pos_H4Y1R = new TVector3();
  pos_H4Y1L = new TVector3();
  pos_H4Y2L = new TVector3();
  pos_H4Y2R = new TVector3();
  pos_H4T = new TVector3();
  pos_H4B = new TVector3();

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
  int nHits = hitVector->size();
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

  for(int hit = 0; hit < nHits; hit++){
  SQHit *sqhit = hitVector->at(hit);
  
  // truth hit info in D0U
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 1){
     pos_D0U->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
     mom_D0U->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
   }
  
  // truth hit info in D0Up
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 2){
	pos_D0Up->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	mom_D0Up->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
	}

  // truth hit info in D0Xp
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 3){
	pos_D0Xp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	mom_D0Xp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
	}

  // truth hit info in D0X
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 4){
	pos_D0X->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	mom_D0X->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
	}
  
  // truth hit info in D0V
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 5){
	pos_D0V->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	mom_D0V->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
	}

  // truth hit info D0Vp
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 6){
	pos_D0Vp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	mom_D0Vp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
	}
 
  // truth hit info H1B
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 31){
	pos_H1B->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	}
  
  // truth hit info H1T
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 32){
	pos_H1B->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	}

  // truth hit info H1L
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 33){
	pos_H1L->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	}

  // truth hit info H1R
  if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 34){
	pos_H1R->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
	}

 // truth hit info in D2V
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 13){
    pos_D2V->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2V->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D2Vp
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 14){
    pos_D2Vp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2Vp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D2Xp
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 15){
    pos_D2Xp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2Xp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D2X
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 16){
    pos_D2X->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2X->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D2U
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 17){
    pos_D2U->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2U->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D2Up
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 18){
    pos_D2Up->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D2Up->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in H2L
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 35){
    pos_H2L->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H2R
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 36){
    pos_H2R->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H2B
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 37){
    pos_H2B->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H2T
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 38){
    pos_H2T->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in D3Vp
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 19 || sqhit->get_detector_id() == 25){
    pos_D3Vp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3Vp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D3V
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 20 || sqhit->get_detector_id() == 26){
    pos_D3V->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3V->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D3Xp
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 21 || sqhit->get_detector_id() == 27){
    pos_D3Xp->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3Xp->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D3X
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 22 || sqhit->get_detector_id() == 28){
    pos_D3X->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3X->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D3Up
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 23 || sqhit->get_detector_id() == 29){
    pos_D3Up->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3Up->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in D3U
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() ==  24 || sqhit->get_detector_id() == 30){
    pos_D3U->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    mom_D3U->SetXYZ(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
    }

 // truth hit info in H3B
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 39){
    pos_H3B->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H3T
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 40){
    pos_H3T->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4Y1L
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 41){
    pos_H4Y1L->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4Y1R
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 42){
    pos_H4Y1R->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4Y2R
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 43){
    pos_H4Y1R->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4Y2L
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 44){
    pos_H4Y2L->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4B
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 45){
    pos_H4B->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

 // truth hit info in H4T
 if(sqhit->get_track_id()==track->get_track_id() && sqhit->get_detector_id() == 46){
    pos_H4T->SetXYZ(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
    }

  }
    

if(AllChamberPlaneHits(track->get_track_id(), hitVector)) {*acc_mom = *momvtx;}
else {acc_mom->SetXYZ(-999., -999., -999.);}

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
/*
///abi add for legacy container
   if(legacyContainer && nRecTracks>0){
	for(int i = 0; i < nRecTracks; ++i)
	    {
		
	      SRecTrack* recTrack =  &recEvent->getTrack(i); 
	      *rec_mom1 = recTrack->getMomentumVecSt1();
	      *rec_momvtx = recTrack->getVertexMom();
	      *rec_posvtx = recTrack->getVertexPos();
	      *rec_momtgt = recTrack->getTargetMom();
	      *rec_postgt = recTrack->getTargetPos();
	    }
	}
    else
    {
      rec_mom1->SetXYZ(-999., -999., -999.);
      rec_momvtx->SetXYZ(-999., -999., -999.);
      rec_posvtx->SetXYZ(-999., -999., -999.);
    }

*/
    saveTree1->Fill();
  }

  int nDimuons = dimuonVector->size();
  int nRecDimuons = legacyContainer ? recEvent->getNDimuons() : (recDimuonVector ? recDimuonVector->size() : -1);

//std::cout << "recEvent->getNDimuons() : " << recEvent->getNDimuons() << std::endl;
//std::cout << "recDimuonVector->size() : " << recDimuonVector->size() << std::endl;

//std::cout << "nDimuons : " << nDimuons << std::endl;
//std::cout << "nRecDimuons : " << nRecDimuons << std::endl;

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
std::cout << "eventID : " << eventID << " nDimuons : " << nDimuons << " nRecDimuons : " << nRecDimuons << " recid : " << recid << std::endl;
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
      rec_mass = -999.;
    }
/*
///@temporary modification by Abi to work with lecayContainer
if(recid <0  && nDimuons==nRecDimuons)
    {
      SRecDimuon* recDimuon = &(recEvent->getDimuon(0));
      rec_mass = recDimuon->mass;
      *rec_pmom = recDimuon->p_pos.Vect();
      *rec_nmom = recDimuon->p_neg.Vect();
      *rec_ppos = recDimuon->vtx_pos;
      *rec_npos = recDimuon->vtx_neg;
      *rec_vtx  = recDimuon->vtx;
      // rec_opt_vz = recDimuon->vtx_opt_z;
      //rec_mass_pred = (recDimuon->p_neg_single+recDimuon->p_pos_single).M();
      //rec_mass_1fit = recDimuon->mass_filtered;
      //chisq_target=recDimuon->chisq_target;
      //chisq_dump = recDimuon->chisq_dump; 
      //chisq_upstream = recDimuon->chisq_upstream;
    }
    else
    {
      rec_mass = -999.;
    }
///@Abi
*/
    saveTree2->Fill();
  }

  ++eventID;

//std::cout << "eventID : " << eventID << std::endl;

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
  mcEvent = findNode::getClass<SQMCEvent>(topNode, "SQMCEvent");
  if(!hitVector || !trackVector || !dimuonVector ||!mcEvent)
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
  saveTree1->Branch("acc_mom", &acc_mom, 256000, 99);
  saveTree1->Branch("rec_mom1", &rec_mom1, 256000, 99);
  saveTree1->Branch("rec_momvtx", &rec_momvtx, 256000, 99);
  saveTree1->Branch("rec_posvtx", &rec_posvtx, 256000, 99);
  saveTree1->Branch("rec_momtgt", &rec_momtgt, 256000, 99);
  saveTree1->Branch("rec_postgt", &rec_postgt, 256000, 99);

  saveTree1->Branch("mom_D0U", &mom_D0U, 256000, 99);
  saveTree1->Branch("pos_D0U", &pos_D0U, 256000, 99);
  saveTree1->Branch("mom_D0Up", &mom_D0Up, 256000, 99);
  saveTree1->Branch("pos_D0Up", &pos_D0Up, 256000, 99);
  saveTree1->Branch("mom_D0X", &mom_D0X, 256000, 99);
  saveTree1->Branch("pos_D0X", &pos_D0X, 256000, 99);
  saveTree1->Branch("mom_D0Xp", &mom_D0Xp, 256000, 99);
  saveTree1->Branch("pos_D0Xp", &pos_D0Xp, 256000, 99);
  saveTree1->Branch("mom_D0V", &mom_D0V, 256000, 99);
  saveTree1->Branch("pos_D0V", &pos_D0V, 256000, 99);
  saveTree1->Branch("mom_D0Vp", &mom_D0Vp, 256000, 99);
  saveTree1->Branch("pos_D0Vp", &pos_D0Vp, 256000, 99);

  saveTree1->Branch("mom_D2U", &mom_D2U, 256000, 99);
  saveTree1->Branch("pos_D2U", &pos_D2U, 256000, 99);
  saveTree1->Branch("mom_D2Up", &mom_D2Up, 256000, 99);
  saveTree1->Branch("pos_D2Up", &pos_D2Up, 256000, 99);
  saveTree1->Branch("mom_D2X", &mom_D2X, 256000, 99);
  saveTree1->Branch("pos_D2X", &pos_D2X, 256000, 99);
  saveTree1->Branch("mom_D2Xp", &mom_D2Xp, 256000, 99);
  saveTree1->Branch("pos_D2Xp", &pos_D2Xp, 256000, 99);
  saveTree1->Branch("mom_D2V", &mom_D2V, 256000, 99);
  saveTree1->Branch("pos_D2V", &pos_D2V, 256000, 99);
  saveTree1->Branch("mom_D2Vp", &mom_D2Vp, 256000, 99);
  saveTree1->Branch("pos_D2Vp", &pos_D2Vp, 256000, 99);

  saveTree1->Branch("mom_D3U", &mom_D3U, 256000, 99);
  saveTree1->Branch("pos_D3U", &pos_D3U, 256000, 99);
  saveTree1->Branch("mom_D3Up", &mom_D3Up, 256000, 99);
  saveTree1->Branch("pos_D3Up", &pos_D3Up, 256000, 99);
  saveTree1->Branch("mom_D3X", &mom_D3X, 256000, 99);
  saveTree1->Branch("pos_D3X", &pos_D3X, 256000, 99);
  saveTree1->Branch("mom_D3Xp", &mom_D3Xp, 256000, 99);
  saveTree1->Branch("pos_D3Xp", &pos_D3Xp, 256000, 99);
  saveTree1->Branch("mom_D3V", &mom_D3V, 256000, 99);
  saveTree1->Branch("pos_D3V", &pos_D3V, 256000, 99);
  saveTree1->Branch("mom_D3Vp", &mom_D3Vp, 256000, 99);
  saveTree1->Branch("pos_D3Vp", &pos_D3Vp, 256000, 99);

  saveTree1->Branch("pos_H1T", &pos_H1T, 256000, 99);
  saveTree1->Branch("pos_H1B", &pos_H1B, 256000, 99);
  saveTree1->Branch("pos_H1L", &pos_H1L, 256000, 99);
  saveTree1->Branch("pos_H1R", &pos_H1R, 256000, 99);

  saveTree1->Branch("pos_H2T", &pos_H2T, 256000, 99);
  saveTree1->Branch("pos_H2B", &pos_H2B, 256000, 99);
  saveTree1->Branch("pos_H2L", &pos_H2L, 256000, 99);
  saveTree1->Branch("pos_H2R", &pos_H2R, 256000, 99);

  saveTree1->Branch("pos_H3T", &pos_H3T, 256000, 99);
  saveTree1->Branch("pos_H3B", &pos_H3B, 256000, 99);
  saveTree1->Branch("pos_H4T", &pos_H4T, 256000, 99);
  saveTree1->Branch("pos_H4B", &pos_H4B, 256000, 99);

  saveTree1->Branch("pos_H4Y1L", &pos_H4Y1L, 256000, 99);
  saveTree1->Branch("pos_H4Y1R", &pos_H4Y1R, 256000, 99);
  saveTree1->Branch("pos_H4Y2L", &pos_H4Y2L, 256000, 99);
  saveTree1->Branch("pos_H4Y2R", &pos_H4Y2R, 256000, 99);

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
