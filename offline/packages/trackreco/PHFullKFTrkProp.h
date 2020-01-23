/*!
 *  \file		  PHFullKFTrkProp.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHFULLKFTRKPROP_H
#define TRACKRECO_PHFULLKFTRKPROP_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

#include "GPUTPCGMPropagator.h"
#include "GPUTPCGMTrackParam.h"
#include "GPUTPCGMPolynomialField.h"

// forward declarations
class PHCompositeNode;

class TrkrClusterContainer;
class SvtxVertexMap;
class SvtxTrackMap;
class AssocInfoContainer;

/// \class PHFullKFTrkProp
///
/// \brief Base class for track seeding
///
class PHFullKFTrkProp : public SubsysReco
{
 public:
  PHFullKFTrkProp(const std::string &name = "PHFullKFTrkProp");
  virtual ~PHFullKFTrkProp() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process() = 0;

  ///
  virtual int End() = 0;


  //SvtxClusterMap *_cluster_map;
  TrkrClusterContainer *_cluster_map;
  SvtxVertexMap *_vertex_map;
  SvtxTrackMap *_track_map;
  AssocInfoContainer *_assoc_container;
  std::string _track_map_name;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
  GPUTPCGMTrackParam _aTrack;
  GPUTPCGMPropagator _aProp;
  GPUTPCGMPolynomialField _aField;
  float _max_sin_phi;
};

#endif
