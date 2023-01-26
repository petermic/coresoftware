#ifndef ALICEKF_H
#define ALICEKF_H

#include "GPUTPCTrackParam.h"
#include <trackbase/ClusterErrorPara.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>
#include <phfield/PHFieldConfigv1.h>

#include <Acts/Definitions/Algebra.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <utility>

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
using TrackSeedAliceSeedMap = std::pair<std::vector<TrackSeed_v1>,std::vector<Eigen::Matrix<float,6,6>>>;

class ALICEKF
{
  public:
  ALICEKF(PHCompositeNode *topNode, 
          TrkrClusterContainer* cmap, 
	        PHField *B,
          float fieldDir,
          unsigned int min_clusters,
          float max_sin_phi,
          int verbosity)
  {
    if(!topNode) std::cout << "no topnode, too bad..." << std::endl;
    _B = B;
    _cluster_map = cmap;
    _fieldDir = fieldDir;
    _max_sin_phi = max_sin_phi;
    _v = verbosity;
    _min_clusters_per_track = min_clusters;
    _ClusErrPara = new ClusterErrorPara();
  }

  ~ALICEKF() {delete _ClusErrPara;}

  explicit ALICEKF(const ALICEKF&) = delete;
  ALICEKF& operator=(const ALICEKF&) = delete;

  // main KF function
  TrackSeedAliceSeedMap ALICEKalmanFilter(const std::vector<std::vector<TrkrDefs::cluskey>>& chains, bool use_nhits_limit, const PositionMap& globalPositions, std::vector<float>& trackChi2) const;

  // encapsulated methods used both here and in propagator
  bool InitializeSeed(const std::vector<TrkrDefs::cluskey> &clusters, GPUTPCTrackParam &trackSeed, const PositionMap &globalPositions, int nseeds) const;
  bool Transport(GPUTPCTrackParam &trackSeed, float nextAlice_x, float newPhi, float alpha, GPUTPCTrackLinearisation &trackLine, GPUTPCTrackParam::GPUTPCTrackFitParam &fp, int nseeds) const;
  bool ConvertToTrackSeedv1(GPUTPCTrackParam &trackSeed, const std::vector<TrkrDefs::cluskey> &trackKeyChain, TrackSeed_v1 &track, float phi, const PositionMap& globalPositions, int nseeds) const;

  // utilities mainly used internally, but might be publicly useful
  bool covIsPosDef(Eigen::Matrix<float,6,6>& cov) const;
  void repairCovariance(Eigen::Matrix<float,6,6>& cov) const;
  bool checknan(float val, const std::string &msg, int num) const;
  double get_Bz(double x, double y, double z) const;
  std::vector<double> GetCircleClusterResiduals(const std::vector<std::pair<double,double>>& pts, double R, double X0, double Y0) const;
  std::vector<double> GetLineClusterResiduals(const std::vector<std::pair<double,double>>& pts, double A, double B) const; 
  Eigen::Matrix<float,6,6> TransformCovarianceMatrix(GPUTPCTrackParam &trackSeed, float track_phi, int nseeds) const;
  float getClusterError(TrkrCluster* c, TrkrDefs::cluskey key, Acts::Vector3 global, int i, int j) const;

  // setters and getters
  void useConstBField(bool opt) {_use_const_field = opt;}
  void useFixedClusterError(bool opt) {_use_fixed_clus_error = opt;}
  void setFixedClusterError(int i,float val) {_fixed_clus_error.at(i)=val;}
  void set_cluster_version(int value) { m_cluster_version = value; }
  void set_n_transport_steps(int n) {_n_transport_steps = n;}
  float get_Bzconst() const { return _Bzconst; }
  int Verbosity() const { return _v; }

  private:

  PHField* _B = nullptr;
  size_t _min_clusters_per_track = 20;
  TrkrClusterContainer* _cluster_map = nullptr;
  ClusterErrorPara *_ClusErrPara;

  int _v = 0;
  float _Bzconst = 10*0.000299792458f;
  float _fieldDir = -1;
  float _max_sin_phi = 1.;
  bool _use_const_field = false;
  bool _use_fixed_clus_error = true;
  std::array<float,3> _fixed_clus_error = {.1,.1,.1};
  int _n_transport_steps = 2;

  int m_cluster_version = 4;
};

#endif
