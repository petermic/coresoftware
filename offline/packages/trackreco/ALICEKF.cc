#include "ALICEKF.h"

#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase/ClusterErrorPara.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <TMatrixFfwd.h>
#include <TMatrixT.h>   
#include <TMatrixTUtils.h>

using keylist = std::vector<TrkrDefs::cluskey>;

// anonymous namespace for local functions
namespace
{
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
}

bool ALICEKF::checknan(float val, const std::string &name, int num) const
{
  if(std::isnan(val))
  {
    if(Verbosity()>0) std::cout << "WARNING: " << name << " is NaN for seed " << num << ". Aborting this seed.\n";
  }
  return std::isnan(val);
}

float ALICEKF::get_Bz(float x, float y, float z) const
{
  if(_use_const_field) return 1.4;
  float p[4] = {x*cm,y*cm,z*cm,0.*cm};
  float bfield[3];
  _B->GetFieldValue(p,bfield);
  return bfield[2]/tesla;
}

float ALICEKF::getClusterError(TrkrCluster* c, TrkrDefs::cluskey key, Acts::Vector3 global, int i, int j) const
{
  if(_use_fixed_clus_error) 
  {
     if(i==j) return _fixed_clus_error.at(i)*_fixed_clus_error.at(i);
     else return 0.;
  }
  else 
    {
      TMatrixF localErr(3,3);
      if(m_cluster_version==3){
        localErr[0][0] = 0.;
        localErr[0][1] = 0.;
        localErr[0][2] = 0.;
        localErr[1][0] = 0.;
        localErr[1][1] = c->getActsLocalError(0,0);
        localErr[1][2] = c->getActsLocalError(0,1);
        localErr[2][0] = 0.;
        localErr[2][1] = c->getActsLocalError(1,0);
        localErr[2][2] = c->getActsLocalError(2,0);
      }else if(m_cluster_version==4){
        std::pair<float, float> para_errors = _ClusErrPara->get_fix_tpc_cluster_error(c,key);
        localErr[0][0] = 0.;
        localErr[0][1] = 0.;
        localErr[0][2] = 0.;
        localErr[1][0] = 0.;
        localErr[1][1] = para_errors.first;
        localErr[1][2] = 0.;
        localErr[2][0] = 0.;
        localErr[2][1] = 0.;
        localErr[2][2] = para_errors.second;
      }
      float clusphi = atan2(global(1), global(0));
      TMatrixF ROT(3,3);
      ROT[0][0] = cos(clusphi);
      ROT[0][1] = -sin(clusphi);
      ROT[0][2] = 0.0;
      ROT[1][0] = sin(clusphi);
      ROT[1][1] = cos(clusphi);
      ROT[1][2] = 0.0;
      ROT[2][0] = 0.0;
      ROT[2][1] = 0.0;
      ROT[2][2] = 1.0;
      TMatrixF ROT_T(3,3);
      ROT_T.Transpose(ROT);
  
      TMatrixF err(3,3);
      err = ROT * localErr * ROT_T;
      
      return err[i][j];
    }
}

bool ALICEKF::InitializeSeed(const keylist& clusters, GPUTPCTrackParam& trackSeed, const PositionMap& globalPositions, int nseeds)
{
  // Transform sPHENIX coordinates into ALICE-compatible coordinates
  const auto& globalpos = globalPositions.at(trackKeyChain.at(0));
  const float x0 = globalpos(0);
  const float y0 = globalpos(1);
  const float z0 = globalpos(2);
  if(Verbosity()>0) std::cout << "Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << std::endl;
  // ALICE x coordinate = radius
  const float alice_x0 = sqrt(x0*x0+y0*y0);
  // ALICE y coordinate = perpendicular to x and z
  // We intially align our coordinates with the first cluster, so this is 0 by definition
  const float alice_y0 = 0;
  const float alice_z0 = z0;

  // Initialize track and linearisation
  trackSeed.InitParam();
  trackSeed.SetX(alice_x0);
  trackSeed.SetY(alice_y0);
  trackSeed.SetZ(alice_z0);
  const float x = x0;
  const float y = y0;
  const float z = z0;
  const float first_phi = atan2(y0,x0);
  const float alice_x = sqrt(x0*x0+y0*y0);
  const float trackCartesian_x = 0.;
  const float trackCartesian_y = 0.;
  const float trackCartesian_z = 0.;
  
  // Pre-set momentum-based parameters to improve numerical stability
  const auto& secondpos = globalPositions.at(trackKeyChain.at(1));
  const float second_x = secondpos(0);
  const float second_y = secondpos(1);
  const float second_z = secondpos(2);
  const float second_phi = atan2(second_y,second_x);
  const float second_alice_x = second_x*cos(first_phi)+second_y*sin(first_phi);
  const float delta_alice_x = second_alice_x - alice_x0;
  const float second_alice_y = -second_x*sin(first_phi)+second_y*cos(first_phi);
  const float init_SinPhi = second_alice_y / sqrt(square(delta_alice_x) + square(second_alice_y));
  const float delta_z = second_z - z0;
  // first alice y is 0 by definition
  const float init_DzDs = -delta_z / sqrt(square(delta_alice_x) + square(second_alice_y));
  trackSeed.SetSinPhi(init_SinPhi);
  trackSeed.SetDzDs(init_DzDs);
  if(Verbosity()>0)
  {
    std::cout << "Set initial SinPhi to " << init_SinPhi << std::endl;
    std::cout << "Set initial DzDs to " << init_DzDs << std::endl;
  }
  
  // get initial pt estimate
  std::vector<std::pair<float,float>> pts;
  std::transform( trackKeyChain.begin(), trackKeyChain.end(), std::back_inserter( pts ), [&globalPositions]( const TrkrDefs::cluskey& key )
  {
    const auto& clpos = globalPositions.at(key);
    return std::make_pair(clpos(0),clpos(1));
  });
  
  const auto [R, x_center, y_center] = TrackFitUtils::circle_fit_by_taubin( pts );
  if(Verbosity()>1) std::cout << "circle fit parameters: R=" << R << ", X0=" << x_center << ", Y0=" << y_center << std::endl;
  
  // check circle fit success
  /* failed fit will result in infinite momentum for the track, which in turn will break the kalman filter */
  if( std::isnan(R) ) return false;
  
  float init_QPt = 1./(0.3*R/100.*get_Bz(x0,y0,z0));
  // determine charge
  if(Verbosity()>2) std::cout << "phi_first: " << first_phi << std::endl;
  if(Verbosity()>2) std::cout << "phi_second: " << second_phi << std::endl;
  float dphi = second_phi - first_phi;
  if(Verbosity()>2) std::cout << "dphi: " << dphi << std::endl;
  if(dphi>M_PI) dphi = 2*M_PI - dphi;
  if(dphi<-M_PI) dphi = 2*M_PI + dphi;
  if(Verbosity()>2) std::cout << "corrected dphi: " << dphi << std::endl;
  if(dphi<0) init_QPt = -1*init_QPt;
  if(Verbosity()>0) std::cout << "initial QPt: " << init_QPt << std::endl;
  trackSeed.SetQPt(init_QPt);
  return true;
}

bool Transport(GPUTPCTrackParam &trackSeed, float X, float alpha, GPUTPCTrackLinearisation &trackLine, GPUTPCTrackFitParam &fp, int nseeds)
{

  if(!trackSeed.Rotate(alpha/2.,trackLine,_max_sin_phi))
  {
    if(Verbosity()>0) std::cout << "WARNING: Rotate failed! Aborting for this seed..." << std::endl;
    return false;
  }

  if(Verbosity()>1) std::cout << "track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << std::endl;
  if(Verbosity()>1) std::cout << "Transporting from " << alice_x << " to " << nextAlice_x << "..." << std::endl;
  float track_x = trackSeed.GetX()*cos(newPhi)-trackSeed.GetY()*sin(newPhi);
  float track_y = trackSeed.GetX()*sin(newPhi)+trackSeed.GetY()*cos(newPhi);
  float track_z = trackSeed.GetZ();
  if(!trackSeed.TransportToXWithMaterial((nextAlice_x+trackSeed.GetX())/2.,trackLine,fp,_Bzconst*get_Bz(track_x,track_y,track_z),_max_sin_phi))
  {
    if(Verbosity()>0) std::cout << "WARNING: Transport failed! Aborting for this seed..." << std::endl;
    return false;
    }
  if(!trackSeed.Rotate(alpha/2.,trackLine,_max_sin_phi))
  {
    if(Verbosity()>0) std::cout << "WARNING: Rotate failed! Aborting for this seed..." << std::endl;
    return false;
  }
  if(!trackSeed.TransportToXWithMaterial(nextAlice_x,trackLine,fp,_Bzconst*get_Bz(track_x,track_y,track_z),_max_sin_phi)) 
  {
    if(Verbosity()>0) std::cout << "WARNING: Rotate failed! Aborting for this seed..." << std::endl;
    return false;
  }
}

bool ALICEKF::ConvertToTrackSeedv1(GPUTPCTrackParam &trackSeed, TrackSeed_v1 &track, float phi, int nseeds)
{
  if(checknan(track_pt,"pT",nseeds)) return false;
  if(checknan(track_pterr,"pT err",nseeds)) return false;
  //float track_x = trackSeed.GetX()*cos(track_phi)-trackSeed.GetY()*sin(track_phi);
  //float track_y = trackSeed.GetX()*sin(track_phi)+trackSeed.GetY()*cos(track_phi);
  float track_z = trackSeed.GetZ();
  if(checknan(track_z,"z",nseeds)) return false;
  float track_zerr = sqrt(trackSeed.GetErr2Z());
  if(checknan(track_zerr,"zerr",nseeds)) return false;

  // get last cluster phi error 
  auto lcluster = _cluster_map->findCluster(trackKeyChain.back());
  const auto& lclusterglob = globalPositions.at(trackKeyChain.back());
  const float lclusterrad = sqrt(lclusterglob(0)*lclusterglob(0) + lclusterglob(1)*lclusterglob(1));
  float last_cluster_phierr = 0;
  if(m_cluster_version==3)
  {
    last_cluster_phierr = lcluster->getRPhiError() / lclusterrad;
  }
  else if(m_cluster_version==4)
  {
    auto para_errors = _ClusErrPara->get_fix_tpc_cluster_error(lcluster,trackKeyChain.back());
    last_cluster_phierr  = sqrt(para_errors.first);
  }

  // phi error assuming error in track radial coordinate is zero
  float track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
    pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
  if(checknan(track_phierr,"phierr",nseeds)) return false;
  if(Verbosity()>0)
  {
    std::cout << "Track phi = " << atan2(track_py,track_px) << std::endl;
    std::cout << "Track phierr = " << track_phierr << std::endl;
  }
  float track_curvature = trackSeed.GetKappa(_Bzconst*get_Bz(track_x,track_y,track_z));
  if(checknan(track_curvature,"curvature",nseeds)) return false;
  float track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bzconst*get_Bz(track_x,track_y,track_z);
  if(checknan(track_curverr,"curvature error",nseeds)) return false;

  //track.set_vertex_id(_vertex_ids[best_vtx]);
  for (unsigned int j = 0; j < trackKeyChain.size(); ++j)
  {
    track.insert_cluster_key(trackKeyChain.at(j));
  }

  int track_charge = 0;
  if(trackSeed.GetQPt()<0) track_charge = -1 * _fieldDir;
  else track_charge = 1 * _fieldDir;
  
  float s = sin(track_phi);
  float c = cos(track_phi);
  float p = trackSeed.GetSinPhi();
  if(checknan(p,"ALICE sinPhi",nseeds)) return false;
  float d = trackSeed.GetDzDs();
  if(checknan(d,"ALICE dz/ds",nseeds)) return false;
  
  const float* cov = trackSeed.GetCov();
  bool cov_nan = false;
  for(int i=0;i<15;i++)
  {
    if(checknan(cov[i],"covariance element "+std::to_string(i),nseeds)) cov_nan = true;
  }
  if(cov_nan) return false;
  std::vector<std::pair<float,float>> pts;
  std::transform( trackKeyChain.begin(), trackKeyChain.end(), std::back_inserter( pts ), [&globalPositions]( const TrkrDefs::cluskey& key )
  {
    const auto& clpos = globalPositions.at(key);
    return std::make_pair(clpos(0),clpos(1));
  });
  std::vector<std::pair<float,float>> rz_pts;
  std::transform( trackKeyChain.begin(), trackKeyChain.end(), std::back_inserter( pts ), [&globalPositions]( const TrkrDefs::cluskey& key )
  {
    const auto& clpos = globalPositions.at(key);
    return std::make_pair(sqrt(clpos(0)*clpos(0)+clpos(1)*clpos(1)),clpos(2));
  });
  const auto [s, Z0] = TrackFitUtils::line_fit( rz_pts );
  
  /// We set the qoverR to get the good charge estimate from the KF
  /// which helps the Acts fit
  track.set_qOverR(trackSeed.GetQPt()*(0.3*1.4)/100.);
  track.set_X0(x_center);
  track.set_Y0(y_center);
  track.set_Z0(Z0);
  track.set_slope(d);
  return true;
}

std::tuple<TrackSeed_v1,Eigen::Matrix<float,6,6>,float,bool> FailedFit = std::make_tuple(TrackSeed_v1(),Eigen::Matrix<float,6,6>(),-1.,false);

TrackSeedAliceSeedMap ALICEKF::ALICEKalmanFilter(const std::vector<keylist>& trackSeedKeyLists,bool use_nhits_limit, const PositionMap& globalPositions, std::vector<float>& trackChi2) const
{
  std::unordered_set<std::tuple<TrackSeed_v1,Eigen::Matrix<float,6,6>,float>> seedmap;
  std::vector<std::pair<int,keylist>> id_keylistvec;

  for(int i=0;i<trackSeedKeyLists.size();i++)
  {
    id_keylist.push_back(std::make_pair(i,trackSeedKeyLists[i]));
  }

  if(Verbosity()>0) std::cout << "min clusters per track: " << _min_clusters_per_track << "\n";
  
  std::for_each( std::execution::par, id_keylistvec.begin(), id_keylistvec.end(), [&](std::pair<int,keylist> id_keylist)
  {
    int nseeds = id_keylist.first;
    keylist trackKeyChain = id_keylist.second;
    // skip if too small
    if(trackKeyChain.size()<2) return;
    if(use_nhits_limit && trackKeyChain.size() < _min_clusters_per_track) return;
    // ensure consistent cluster order
    if(TrkrDefs::getLayer(trackKeyChain.front())<TrkrDefs::getLayer(trackKeyChain.back())) std::reverse(trackKeyChain.begin(),trackKeyChain.end());

    // initialize seed and related components
    GPUTPCTrackParam trackSeed;
    if(!InitializeSeed(trackKeyChain,trackSeed,globalPositions,nseeds)) return;
    GPUTPCTrackLinearisation trackLine(trackSeed);
    GPUTPCTrackParam::GPUTPCTrackFitParam fp;
    trackSeed.CalculateFitParameters(fp);

    if(Verbosity()>0) std::cout << std::endl << std::endl << "------------------------" << std::endl << "seed size: " << trackKeyChain.size() << std::endl << std::endl << std::endl;
    int cluster_ctr = 1;
    // bool aborted = false;
    
    // starting at second cluster
    for(auto clusterkeyIter = std::next(trackKeyChain.begin()); clusterkeyIter != trackKeyChain.end(); ++clusterkeyIter)
    {
      if(std::isnan(trackSeed.GetX()) ||
         std::isnan(trackSeed.GetY()) ||
         std::isnan(trackSeed.GetZ())) return;

      TrkrDefs::cluskey clusterkey = *clusterkeyIter;

      if(Verbosity()>2) 
      {
        std::cout << "-------------------------------------------------------------" << std::endl;
        std::cout << "cluster " << cluster_ctr << " -> " << cluster_ctr + 1 << std::endl;
        std::cout << "this cluster (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
        std::cout << "layer " << (int)TrkrDefs::getLayer(clusterkey) << std::endl;
      }

      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->find(clusterkey);
      const auto& nextpos = globalPositions.at(clusterkey);

      // find ALICE x-coordinate
      float nextCluster_x = nextpos(0);
      float nextCluster_xerr = sqrt(getClusterError(nextCluster,clusterkey,nextpos,0,0));
      float nextCluster_y = nextpos(1);
      float nextCluster_yerr = sqrt(getClusterError(nextCluster,clusterkey,nextpos,1,1));
      float nextCluster_z = nextpos(2);
      float nextCluster_zerr = sqrt(getClusterError(nextCluster,clusterkey,nextpos,2,2));

      float newPhi = atan2(nextCluster_y,nextCluster_x);
      float oldPhi = atan2(y,x);
      float alpha = newPhi - oldPhi;
      float nextAlice_x = nextCluster_x*cos(newPhi)+nextCluster_y*sin(newPhi);

      if(Verbosity()>1)
      {
        LogDebug("new phi = " << newPhi << std::endl);
        LogDebug("old phi = " << oldPhi << std::endl);
        LogDebug("alpha = " << alpha << std::endl);
      }

      if(!Transport(trackSeed,nextAlice_x,alpha,trackLine,fp,nseeds)) return;

      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging

      float predicted_alice_x = trackSeed.GetX();
      float predicted_alice_y = trackSeed.GetY();
      float predicted_z = trackSeed.GetZ();
      float cos_phi = x/sqrt(x*x+y*y);
      float sin_phi = y/sqrt(x*x+y*y);
      
      if(Verbosity()>0)
      {
        std::cout << "new track ALICE x = " << trackSeed.GetX() << std::endl;
        std::cout << "new track ALICE y = " << trackSeed.GetY() << std::endl;
        std::cout << "new track z = " << trackSeed.GetZ() << std::endl;
        std::cout << "cos_phi = " << cos_phi << std::endl;
        std::cout << "sin phi = " << sin_phi << std::endl;
      }
      
      trackCartesian_x = predicted_alice_x*cos_phi-predicted_alice_y*sin_phi;
      trackCartesian_y = predicted_alice_x*sin_phi+predicted_alice_y*cos_phi;
      trackCartesian_z = predicted_z;

      if(Verbosity()>0)
      {
        std::cout << "Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << std::endl;
        std::cout << "Track position ALICE Y error: " << sqrt(trackSeed.GetCov(0)) << std::endl;
        std::cout << "Track position x error: " << sqrt(trackSeed.GetCov(0))*sin_phi << std::endl;
        std::cout << "Track position y error: " << sqrt(trackSeed.GetCov(0))*cos_phi << std::endl;
        std::cout << "Track position z error: " << sqrt(trackSeed.GetCov(5)) << std::endl;
        std::cout << "Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << std::endl;
        std::cout << "Cluster errors: (" << nextCluster_xerr << ", " << nextCluster_yerr << ", " << nextCluster_zerr << ")" << std::endl;
        std::cout << "track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << std::endl;
      }

      float nextCluster_alice_y = -nextCluster_x*sin(newPhi)+nextCluster_y*cos(newPhi);
      float y2_error = getClusterError(nextCluster,*clusterkey,nextpos,0,0)*sin(newPhi)*sin(newPhi)+2*getClusterError(nextCluster,*clusterkey,nextpos,0,1)*cos(newPhi)*sin(newPhi)+getClusterError(nextCluster,*clusterkey,nextpos,1,1)*cos(newPhi)*cos(newPhi);
      float z2_error = getClusterError(nextCluster,*clusterkey,nextpos,2,2);

      if(Verbosity()>0)
      {
        std::cout << "track ALICE SinPhi = " << trackSeed.GetSinPhi() << std::endl;
        std::cout << "track DzDs = " << trackSeed.GetDzDs() << std::endl;
        std::cout << "chi2 = " << trackSeed.GetChi2() << std::endl;
        std::cout << "NDF = " << trackSeed.GetNDF() << std::endl;
        std::cout << "chi2 / NDF = " << trackSeed.GetChi2()/trackSeed.GetNDF() << std::endl;
      }
  
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,_max_sin_phi))
      {
	      if(Verbosity()>0) std::cout << "Kalman filter failed for seed " << nseeds << "! Aborting for this seed..." << std::endl;
        break;
      }

      if(Verbosity()>1)
      {
        float track_pt = 1./trackSeed.GetQPt();
        float track_pY = track_pt*trackSeed.GetSinPhi();
        float track_pX = sqrt(track_pt*track_pt-track_pY*track_pY);
        float track_px = track_pX*cos(newPhi)-track_pY*sin(newPhi);
        float track_py = track_pX*sin(newPhi)+track_pY*cos(newPhi);
        float track_pz = -track_pt*trackSeed.GetDzDs();
        float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
        std::cout << "track pt = " << track_pt << " +- " << track_pterr << std::endl;
        std::cout << "track ALICE p = (" << track_pX << ", " << track_pY << ", " << track_pz << ")" << std::endl;
        std::cout << "track p = (" << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;
      }
      x = nextCluster_x;
      y = nextCluster_y;
      z = nextCluster_z;
      alice_x = nextAlice_x;
      ++cluster_ctr;
	    float nextclusrad = sqrt(nextCluster_x*nextCluster_x + nextCluster_y*nextCluster_y);
	    float nextclusphierr = 0;
	    if(m_cluster_version==3)
      {
	      nextclusphierr = nextCluster->getRPhiError() / nextclusrad;
	    }
      else if(m_cluster_version==4)
      {
	      auto para_errors = _ClusErrPara->get_fix_tpc_cluster_error(nextCluster,*clusterkey);
	      nextclusphierr = sqrt(para_errors.first);
	    }   
    }

    float track_phi = atan2(y,x);
    float track_pt = fabs(1./trackSeed.GetQPt());
    float track_pY = track_pt*trackSeed.GetSinPhi();
    float track_pX = sqrt(track_pt*track_pt-track_pY*track_pY);
    float track_px = track_pX*cos(track_phi)-track_pY*sin(track_phi);
    float track_py = track_pX*sin(track_phi)+track_pY*cos(track_phi);
    float track_pz = track_pt*trackSeed.GetDzDs();
    float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());

    // If Kalman filter doesn't do its job (happens often with short seeds), use the circle-fit estimate as the central value
    // if(trackKeyChain.size()<10) track_pt = fabs(1./init_QPt);

    if(Verbosity()>0)
    {
      std::cout << "track pt = " << track_pt << " +- " << track_pterr << std::endl;
      std::cout << "track ALICE p = (" << track_pX << ", " << track_pY << ", " << track_pz << ")" << std::endl;
      std::cout << "track p = (" << track_px << ", " << track_py << ", " << track_pz << ")" << std::endl;
      std::cout << "Track pterr = " << track_pterr << std::endl;
    }

/*    
    if(cluster_ctr!=1 && !trackSeed.CheckNumericalQuality())
    {
      std::cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
      aborted = true;
      continue;
    } 
*/    
    TrackSeed_v1 track;
    if(!ConvertToTrackSeedv1(trackSeed,track,track_phi,nseeds)) return;
    Eigen::Matrix<float,6,6> scov = TransformCovarianceMatrix(trackSeed,track_phi);

    seedmap.insert(std::make_tuple(track,scov,trackSeed.GetChi2()/trackSeed.GetNDF()));
  });

  if(Verbosity()>0) std::cout << "number of seeds: " << seedmap.size() << "\n";

  std::vector<TrackSeed_v1> seeds_vector;
  std::vector<Eigen::Matrix<float,6,6>> alice_seeds_vector;

  for(auto seedtriple : seedmap)
  {
    seeds.push_back(std::get<0>(seedtriple));
    alice_seeds_vector.push_back(std::get<1>(seedtriple));
    trackChi2.push_back(std::get<2>(seedtriple));
  }

  return std::make_pair(seeds_vector, alice_seeds_vector);

}

bool ALICEKF::covIsPosDef(Eigen::Matrix<float,6,6>& cov) const
{
  // attempt Cholesky decomposition
  Eigen::LLT<Eigen::Matrix<float,6,6>> chDec(cov);
  // if Cholesky decomposition does not exist, matrix is not positive definite
  return (chDec.info() != Eigen::NumericalIssue);
}

void ALICEKF::repairCovariance(Eigen::Matrix<float,6,6>& cov) const
{
  Eigen::Matrix<float,6,6> repaircov = cov;
  // find closest positive definite matrix
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float,6,6>> solver(repaircov);
  Eigen::Matrix<float,6,1> D = solver.eigenvalues();
  Eigen::Matrix<float,6,6> Q = solver.eigenvectors();
  Eigen::Matrix<float,6,1> Dp = D.cwiseMax(1e-15);
  Eigen::Matrix<float,6,6> Z = Q*Dp.asDiagonal()*Q.transpose();
  // updates covariance matrix
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      cov(i,j) = Z(i,j);
    }
  }
  
}

std::vector<float> ALICEKF::GetCircleClusterResiduals(const std::vector<std::pair<float,float>>& points, float R, float X0, float Y0) const
{
  std::vector<float> residues;
  std::transform( points.begin(), points.end(), std::back_inserter( residues ), [R,X0,Y0]( const std::pair<float,float>& point )
  {
    float x = point.first;
    float y = point.second;

    // The shortest distance of a point from a circle is along the radial line from the circle center to the point
    return sqrt( square(x-X0) + square(y-Y0) )  -  R;  
  } );
  return residues;  
}

std::vector<float> ALICEKF::GetLineClusterResiduals(const std::vector<std::pair<float,float>>& points, float A, float B) const
{
  std::vector<float> residues;
  // calculate cluster residuals from the fitted circle
  std::transform( points.begin(), points.end(), std::back_inserter( residues ), [A,B]( const std::pair<float,float>& point )
  {
    float r = point.first;
    float z = point.second;
    
    // The shortest distance of a point from a circle is along the radial line from the circle center to the point
    
    float a = -A;
    float b = 1.0;
    float c = -B;
    return std::abs(a*r+b*z+c)/sqrt(square(a)+square(b));
  });
  return residues;  
}

Eigen::Matrix<float,6,6> ALICEKF::TransformCovarianceMatrix(GPUTPCTrackParam &trackSeed, float track_phi)
{
  float track_pt = 1./trackSeed.GetQPt();
  int track_charge = 0;
  if(trackSeed.GetQPt()<0) track_charge = -1 * _fieldDir;
  else track_charge = 1 * _fieldDir;
  float s = sin(track_phi);
  float c = cos(track_phi);
  float p = trackSeed.GetSinPhi();
  float d = trackSeed.GetDzDs();
  // make this into an Eigen matrix
  Eigen::Matrix<float,5,5> ecov;
  ecov(0,0)=cov[0];
  ecov(0,1)=cov[1];
  ecov(0,2)=cov[2];
  ecov(0,3)=cov[3];
  ecov(0,4)=cov[4];
  ecov(1,1)=cov[5];
  ecov(1,2)=cov[6];
  ecov(1,3)=cov[7];
  ecov(1,4)=cov[8];
  ecov(2,2)=cov[9];
  ecov(2,3)=cov[10];
  ecov(2,4)=cov[11];
  ecov(3,3)=cov[12];
  ecov(3,4)=cov[13];
  ecov(4,4)=cov[14];
  // symmetrize
  ecov(1,0)=ecov(0,1);
  ecov(2,0)=ecov(0,2);
  ecov(3,0)=ecov(0,3);
  ecov(4,0)=ecov(0,4);
  ecov(2,1)=ecov(1,2);
  ecov(3,1)=ecov(1,3);
  ecov(4,1)=ecov(1,4);
  ecov(3,2)=ecov(2,3);
  ecov(4,2)=ecov(2,4);
  ecov(4,3)=ecov(3,4);
  // make rotation matrix based on the following:
  // x = X*cos(track_phi) - Y*sin(track_phi)
  // y = X*sin(track_phi) + Y*cos(track_phi)
  // z = Z
  // pY = pt*sinphi
  // pX = sqrt(pt**2 - pY**2)
  // px = pX*cos(track_phi) - pY*sin(track_phi)
  // py = pX*sin(track_phi) + pY*cos(track_phi)
  // pz = pt*(dz/ds)
  Eigen::Matrix<float,6,5> J;
  J(0,0) = -s; // dx/dY
  J(0,1) = 0.; // dx/dZ
  J(0,2) = 0.; // dx/d(sinphi)
  J(0,3) = 0.; // dx/d(dz/ds)
  J(0,4) = 0.; // dx/d(Q/pt)

  J(1,0) = c;  // dy/dY
  J(1,1) = 0.; // dy/dZ
  J(1,2) = 0.; // dy/d(sinphi)
  J(1,3) = 0.; // dy/d(dz/ds)
  J(1,4) = 0.; // dy/d(Q/pt)

  J(2,0) = 0.; // dz/dY
  J(2,1) = 1.; // dz/dZ
  J(2,2) = 0.; // dz/d(sinphi)
  J(2,3) = 0.; // dz/d(dz/ds)
  J(2,4) = 0.; // dz/d(Q/pt)

  J(3,0) = 0.; // dpx/dY
  J(3,1) = 0.; // dpx/dZ
  J(3,2) = -track_pt*(p*c/sqrt(1-p*p)+s); // dpx/d(sinphi)
  J(3,3) = 0.; // dpx/d(dz/ds)
  J(3,4) = track_pt*track_pt*track_charge*(p*s-c*sqrt(1-p*p)); // dpx/d(Q/pt)

  J(4,0) = 0.; // dpy/dY
  J(4,1) = 0.; // dpy/dZ
  J(4,2) = track_pt*(c-p*s/sqrt(1-p*p)); // dpy/d(sinphi)
  J(4,3) = 0.; // dpy/d(dz/ds)
  J(4,4) = -track_pt*track_pt*track_charge*(p*c+s*sqrt(1-p*p)); // dpy/d(Q/pt)

  J(5,0) = 0.; // dpz/dY
  J(5,1) = 0.; // dpz/dZ
  J(5,2) = 0.; // dpz/d(sinphi)
  J(5,3) = track_pt; // dpz/d(dz/ds)
  J(5,4) = -track_pt*track_pt*track_charge*d; // dpz/d(Q/pt)
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<5;j++)
    {
      checknan(J(i,j),"covariance rotator element ("+std::to_string(i)+","+std::to_string(j)+")",nseeds);
    }
  }

  // the heavy lifting happens here
  Eigen::Matrix<float,6,6> scov = J*ecov*J.transpose();
  if(!covIsPosDef(scov))
  {
    repairCovariance(scov);
  }
  return scov;
  /*
  // Derived from:
  // 1) Taking the Jacobian of the conversion from (Y,Z,SinPhi,DzDs,Q/Pt) to (x,y,z,px,py,pz)
  // 2) Computing (Jacobian)*(ALICE covariance matrix)*(transpose of Jacobian)
  track.set_error(0, 0, cov[0]*s*s);
  track.set_error(0, 1, -cov[0]*c*s);
  track.set_error(0, 2, -cov[1]*s);
  track.set_error(0, 3, cov[2]*s*s/q-cov[4]*s*(-c/(q*q)+p*s/(q*q)));
  track.set_error(0, 4, -cov[2]*c*s/q-cov[4]*s*(-c*p/(q*q)-s/(q*q)));
  track.set_error(0, 5, cov[4]*d*s/(q*q)-cov[3]*s/q);
  track.set_error(1, 1, cov[0]*c*c);
  track.set_error(1, 2, cov[1]*c);
  track.set_error(1, 3, -cov[2]*c*s/q+cov[4]*c*(-c/(q*q)+p*s/(q*q)));
  track.set_error(1, 4, cov[2]*c*c/q+cov[4]*c*(-c*p/(q*q)-s/(q*q)));
  track.set_error(1, 5, cov[4]*d*c/(q*q)+cov[3]*c/q);
  track.set_error(2, 2, cov[5]);
  track.set_error(2, 3, -cov[6]*s/q+cov[8]*(-c/(q*q)+p*s/(q*q)));
  track.set_error(2, 4, cov[6]*c/q+cov[8]*(-c*p/(q*q)-s/(q*q)));
  track.set_error(2, 5, -cov[8]*d/(q*q)+cov[7]/q);
  track.set_error(3, 3, cov[9]*s*s/(q*q)-cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c/(q*q)+p*s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
  track.set_error(3, 4, -cov[9]*c*s/(q*q)+cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c*p/(q*q)-s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
  track.set_error(3, 5, -cov[10]*s/(q*q)+cov[13]/q*(-c/(q*q)+p*s/(q*q))-d/(q*q)*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
  track.set_error(4, 4, c/q*(c/q*cov[9]+cov[11]*(-c*p/(q*q)-s/(q*q)))+(-c*p/(q*q)-s/(q*q))*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
  track.set_error(4, 5, cov[10]*c/(q*q)+cov[13]/q*(-c*p/(q*q)-s/(q*q))-d/(q*q)*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
  track.set_error(5, 5, -d/(q*q)*(-d*cov[14]/(q*q)+cov[13]/q)-d*cov[13]/(q*q*q)+cov[12]/(q*q));
  // symmetrize covariance
  track.set_error(1, 0, track.get_error(0, 1));
  track.set_error(2, 0, track.get_error(0, 2));
  track.set_error(3, 0, track.get_error(0, 3));
  track.set_error(4, 0, track.get_error(0, 4));
  track.set_error(5, 0, track.get_error(0, 5));
  track.set_error(2, 1, track.get_error(1, 2));
  track.set_error(3, 1, track.get_error(1, 3));
  track.set_error(4, 1, track.get_error(1, 4));
  track.set_error(5, 1, track.get_error(1, 5));
  track.set_error(3, 2, track.get_error(2, 3));
  track.set_error(4, 2, track.get_error(2, 4));
  track.set_error(5, 2, track.get_error(2, 5));
  track.set_error(4, 3, track.get_error(3, 4));
  track.set_error(5, 3, track.get_error(3, 5));
  track.set_error(5, 4, track.get_error(4, 5));
  */
}
