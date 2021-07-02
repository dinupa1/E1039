#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <ktracker/SRecEvent.h>

class TFile;
class TTree;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;
class SQMCEvent;

class AnaModule: public SubsysReco 
{
public:
  AnaModule(const std::string& name = "AnaModule");
  virtual ~AnaModule();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n) { saveName = n; }
  void set_legacy_rec_container(bool b) { legacyContainer = b; }

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();
  bool AllChamberPlaneHits(const int trk_id, SQHitVector *hitVector);
  TVector3 mom_detector(SQHitVector* hitVector, SQTrack* track, const int det_id);
  TVector3 pos_detector(SQHitVector* hitVector, SQTrack* track, const int det_id);

  bool legacyContainer;

  // Input
  SQHitVector* hitVector;
  SQTrackVector* trackVector;
  SQDimuonVector* dimuonVector;
  SQMCEvent* mcEvent;

  SRecEvent* recEvent;
  SQTrackVector*  recTrackVector;
  SQDimuonVector* recDimuonVector;

  // Output
  TString saveName;
  TFile* saveFile;
  int eventID;

  TTree* saveTree1;
  int charge;
  TVector3* pos1;
  TVector3* pos2;
  TVector3* pos3;
  TVector3* posvtx;
  TVector3* mom1;
  TVector3* mom2;
  TVector3* mom3;
  TVector3* momvtx;
  TVector3* acc_mom;
  TVector3* rec_mom1;
  TVector3* rec_momvtx;
  TVector3* rec_posvtx;
  TVector3* rec_momtgt;
  TVector3* rec_postgt;

  TVector3* mom_D0U;
  TVector3* pos_D0U;
  TVector3* mom_D0Up;
  TVector3* pos_D0Up;
  TVector3* mom_D0X;
  TVector3* pos_D0X;
  TVector3* mom_D0Xp;
  TVector3* pos_D0Xp;
  TVector3* mom_D0V;
  TVector3* pos_D0V;
  TVector3* mom_D0Vp;
  TVector3* pos_D0Vp;

  TVector3* mom_D2U;
  TVector3* pos_D2U;
  TVector3* mom_D2Up;
  TVector3* pos_D2Up;
  TVector3* mom_D2X;
  TVector3* pos_D2X;
  TVector3* mom_D2Xp;
  TVector3* pos_D2Xp;
  TVector3* mom_D2V;
  TVector3* pos_D2V;
  TVector3* mom_D2Vp;
  TVector3* pos_D2Vp;

  TVector3* mom_D3U;
  TVector3* pos_D3U;
  TVector3* mom_D3Up;
  TVector3* pos_D3Up;
  TVector3* mom_D3X;
  TVector3* pos_D3X;
  TVector3* mom_D3Xp;
  TVector3* pos_D3Xp;
  TVector3* mom_D3V;
  TVector3* pos_D3V;
  TVector3* mom_D3Vp;
  TVector3* pos_D3Vp;

  TVector3* pos_H1T;
  TVector3* pos_H1B;
  TVector3* pos_H1L;
  TVector3* pos_H1R;

  TVector3* pos_H2T;
  TVector3* pos_H2B;
  TVector3* pos_H2L;
  TVector3* pos_H2R;

  TVector3* pos_H3T;
  TVector3* pos_H3B;
  TVector3* pos_H4T;
  TVector3* pos_H4B;

  TVector3* pos_H4Y1L;
  TVector3* pos_H4Y1R;
  TVector3* pos_H4Y2L;
  TVector3* pos_H4Y2R;

  TTree* saveTree2;
  TVector3* pmom;
  TVector3* nmom;
  TVector3* rec_pmom;
  TVector3* rec_nmom;
  TVector3* rec_ppos;
  TVector3* rec_npos;
  TVector3* rec_vtx;
  TVector3* vtx;
  double mass;
  double rec_mass;
  double mass_acc;

  bool st1hit;
  bool st2hit;
  bool st3phit;
  bool st3mhit;
  bool st3hit;

};

#endif
