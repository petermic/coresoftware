#include "PHFullKFTrkProp.h"

#include "AssocInfoContainer.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <iostream>                            // for operator<<, basic_ostream

#define _DEBUG_

#ifdef _DEBUG_
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;

PHFullKFTrkProp::PHFullKFTrkProp(const std::string& name)
  : SubsysReco(name)
  , _cluster_map(nullptr)
  , _vertex_map(nullptr)
  , _track_map(nullptr)
  , _assoc_container(nullptr)
  , _track_map_name("SvtxTrackMap")
  , _aTrack()
  , _aProp()
  , _aField()
  , _max_sin_phi(0.999)
{
}

int PHFullKFTrkProp::InitRun(PHCompositeNode* topNode)
{
  return Setup(topNode);
}

int PHFullKFTrkProp::process_event(PHCompositeNode* topNode)
{
  return Process();
}

int PHFullKFTrkProp::Process()
{
  // get tracklets
  for(auto phtrk_iter = _track_map->begin();
    phtrk_iter != _track_map->end();
    ++phtrk_iter)
  {
    SvtxTrack* tracklet = phtrk_iter->second;
    // set initial alpha to be aligned with track head
    _aTrack.X() = sqrt(pow(tracklet->get_x(),2)+pow(tracklet->get_y(),2));
    _aTrack.Y() = 0.;
    _aTrack.Z() = tracklet->get_z();
    float alpha = atan(tracklet->get_y()/tracklet->get_x());
    _aTrack.SinPhi() = 0.;
    _aTrack.DzDs() = 0.;
    _aTrack.QPt() = tracklet->get_charge()/tracklet->get_pt();
    _aTrack.ResetCovariance();
    _aProp.SetMaterial(0.,1e100); // material corrections not needed, but these are fairly straightforward parameters for the vacuum
    _aProp.SetMaxSinPhi(_max_sin_phi);
    _aProp.SetToyMCEventsFlag(false);
    _aProp.SetSeedingErrors(true);
    _aProp.SetFitInProjections(true);
    _aField.SetFieldNominal(1.4);
    float BxCoeff[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float ByCoeff[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float BzCoeff[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    _aField.SetFieldTpc(BxCoeff,ByCoeff,BzCoeff);
    _aProp.SetPolynomialField(&_aField);
    _aProp.SetTrack(&_aTrack,alpha);
    
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHFullKFTrkProp::End(PHCompositeNode* topNode)
{
  End();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHFullKFTrkProp::Setup(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHFullKFTrkProp::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------


  //_cluster_map = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << _track_map_name << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
