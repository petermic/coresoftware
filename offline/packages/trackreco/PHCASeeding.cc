
/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */

#include "PHCASeeding.h"
#include "ALICEKF.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

// trackbase_historic includes
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v1.h>

//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>

//BOOST for combi seeding
//#define BOOST_NO_CXX98_FUNCTION_BASE
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>  // for pair, make_pair
#include <vector>
#include <algorithm> // for find
#include <set>
#include <execution>
#include <memory>
#include <mutex>

//end

typedef bg::model::point<double, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
typedef std::pair<std::array<double,3>,TrkrDefs::cluskey> coordKey;
typedef std::array<coordKey,2> keylink;
typedef std::array<coordKey,3> triplet;
typedef std::vector<TrkrDefs::cluskey> keylist;

// apparently there is no builtin STL hash function for a std::array
// so to use std::unordered_set (essentially a hash table), we have to make our own hasher

namespace std
{
  template<typename T,size_t N>
  struct hash<std::array<T,N>>
  {
    typedef std::array<T,N> argument_type;
    typedef size_t result_type;

    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for(result_type i = 0; i < N; ++i)
      {
        h = h * 31 + hasher(a[i]);
      }
      return h;
    }
  };
  template<typename A,typename B>
  struct hash<pair<A,B>>
  {
    typedef pair<A,B> argument_type;
    typedef size_t result_type;
    
    result_type operator()(const argument_type& a) const
    {
      hash<A> hashA;
      hash<B> hashB;
      return (hashA(a.first)*31+hashB(a.second));
    }
  }; 
}


// anonymous namespace for local functions
namespace
{
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  /// phi angle of Acts::Vector3
  inline double get_phi( const Acts::Vector3& position )
  {
    double phi = std::atan2( position.y(), position.x() );
    //if( phi < 0 ) phi += 2.*M_PI; //THIS WAS A MAJOR BUG???
    return phi;
  } 
  
  /// pseudo rapidity of Acts::Vector3
  inline double get_eta( const Acts::Vector3& position )
  {
    const double norm = std::sqrt( square(position.x()) + square(position.y()) + square(position.z()) );
    return std::log((norm+position.z())/(norm-position.z()))/2;
  }
  
  ///@name utility convertion functions
  //@{

  coordKey fromPointKey(const pointKey& p)
  { return std::make_pair(std::array<double,3>({p.first.get<0>(),p.first.get<1>(),p.first.get<2>()}),p.second); }

  std::vector<coordKey> fromPointKey(const std::vector<pointKey>& p)
  {
    std::vector<coordKey> output;
    output.resize(p.size());
    std::transform( p.begin(), p.end(), std::back_inserter( output ), []( const pointKey& point )
      { return fromPointKey(point); } );
    return output;
  }
  //@}
  
  double breaking_angle(double x1, double y1, double x2, double y2)
  {
/*
    double l1 = sqrt(x1*x1+y1*y1+z1*z1);
    double l2 = sqrt(x2*x2+y2*y2+z2*z2);
    double sx = (x1/l1+x2/l2);
    double sy = (y1/l1+y2/l2);
    double sz = (z1/l1+z2/l2);
    double dx = (x1/l1-x2/l2);
    double dy = (y1/l1-y2/l2);
    double dz = (z1/l1-z2/l2);
    return 2*atan2(sqrt(dx*dx+dy*dy+dz*dz),sqrt(sx*sx+sy*sy+sz*sz));
*/
    return (x1*x2+y1*y2)/(sqrt(x1*x1+y1*y1)*sqrt(x2*x2+y2*y2));
  }

  std::mutex mtx;
  std::mutex mtxn;
}

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const std::string &name,
    unsigned int start_layer,
    unsigned int end_layer,
    unsigned int min_nhits_per_cluster,
    unsigned int min_clusters_per_track,
    double neighbor_phi_width,
    double neighbor_eta_width,
    double maxSinPhi,
    double cosTheta_limit)
  : PHTrackSeeding(name)
  , _start_layer(start_layer)
  , _end_layer(end_layer)
  , _min_nhits_per_cluster(min_nhits_per_cluster)
  , _min_clusters_per_track(min_clusters_per_track)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_eta_width(neighbor_eta_width)
  , _max_sin_phi(maxSinPhi)
  , _cosTheta_limit(cosTheta_limit)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode *topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!tGeometry)
    {
      std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHCASeeding::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  auto globalpos = tGeometry->getGlobalPosition(key, cluster);

  // check if TPC distortion correction are in place and apply
  if( m_dcc ) { globalpos = m_distortionCorrection.get_corrected_position( globalpos, m_dcc ); }

  return globalpos;
}

void PHCASeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, 
  double phimin, double zmin, double lmin, 
  double phimax, double zmax, double lmax, 
  std::vector<pointKey> &returned_values) const
{
//  double phimin_2pi = phimin;
//  double phimax_2pi = phimax;
//  if (phimin < 0) phimin_2pi = 2*M_PI+phimin;
//  if (phimax > 2*M_PI) phimax_2pi = phimax-2*M_PI;
  rtree.query(bgi::intersects(
      box(point(phimin, zmin, lmin), 
      point(phimax, zmax, lmax))), 
      std::back_inserter(returned_values));
}

PositionMap PHCASeeding::FillTree()
{ 
  t_fill->stop();
  int n_dupli = 0;
  int nlayer[60];

  PositionMap cachedPositions;
  t_fill->restart();
  for (int j = 0; j < 60; ++j) nlayer[j] = 0;
  const auto hitsetkeys = _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
  std::for_each( std::execution::par, hitsetkeys.begin(), hitsetkeys.end(), [&](auto hitsetkey)
  { 
    mtx.lock();
    auto range = _cluster_map->getClusters(hitsetkey);
    std::vector<std::pair<TrkrDefs::cluskey,TrkrCluster*>> clist;
    for(auto clusIter = range.first; clusIter != range.second; ++clusIter)
    {
      clist.push_back(std::make_pair(clusIter->first,clusIter->second));
    } 
    mtx.unlock();
    for(auto clusIter = clist.begin(); clusIter != clist.end(); ++clusIter)
    {
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster *cluster = clusIter->second;
      unsigned int layer = TrkrDefs::getLayer(ckey);
      if (layer < _start_layer || layer >= _end_layer)
      {
        if(Verbosity()>0) std::cout << "invalid layer: " << layer << " for cluster key " << ckey << std::endl;
        continue;
      }
      if(_iteration_map != NULL && _n_iteration >0)
      {
	if( _iteration_map->getIteration(ckey) > 0) continue; // skip hits used in a previous iteration
      }

      // get global position, convert to Acts::Vector3 and store in map
      const Acts::Vector3 globalpos_d = getGlobalPosition(ckey, cluster);

      if(Verbosity() > 3)
      {
        auto global_before = tGeometry->getGlobalPosition(ckey, cluster);
        std::cout << "CA Seeder: Cluster: " << ckey << std::endl;
        std::cout << " Global before: " << global_before[0] << "  " << global_before[1] << "  " << global_before[2] << std::endl;
        std::cout << " Global after   : " << globalpos_d[0] << "  " << globalpos_d[1] << "  " << globalpos_d[2] << std::endl;
      }
      const Acts::Vector3 globalpos = { globalpos_d.x(), globalpos_d.y(), globalpos_d.z() };
      mtxn.lock();
      cachedPositions.insert(std::make_pair(ckey, globalpos));
      mtxn.unlock();
      const double clus_phi = get_phi( globalpos );      
      //const double clus_eta = get_eta( globalpos );
      const double clus_l = layer;

      if(Verbosity() > 2) std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
      
      std::vector<pointKey> testduplicate;
      QueryTree(_rtree, clus_phi - 0.00001, globalpos(2) - 0.00001, layer - 0.5, clus_phi + 0.00001, globalpos(2) + 0.00001, layer + 0.5, testduplicate);
      mtx.lock();
      if (!testduplicate.empty())
      {
        ++n_dupli;
        continue;
      }
      // mirrors the edges of the phi range so queries there still see all available clusters
      if(fabs(M_PI-clus_phi)<_neighbor_phi_width)
      {
        ++nlayer[layer];
        _rtree.insert(std::make_pair(point(clus_phi-2*M_PI,globalpos(2),clus_l),ckey));
      }
      if(fabs(-M_PI-clus_phi)<_neighbor_phi_width)
      {
        ++nlayer[layer];
        _rtree.insert(std::make_pair(point(clus_phi+2*M_PI,globalpos(2),clus_l),ckey));
      }
      ++nlayer[layer];
      //if(clus_l != 22 && clus_l != 23 && clus_l != 39 && clus_l != 40) 
        _rtree.insert(std::make_pair(point(clus_phi, globalpos(2), clus_l), ckey));
      mtx.unlock();
    }
  });
  t_fill->stop();
  if(Verbosity()>1) for (int j = 0; j < 60; ++j) std::cout << "nhits in layer " << j << ":  " << nlayer[j] << std::endl;
  if(Verbosity()>0) std::cout << "fill time: " << t_fill->elapsed() / 1000. << " sec" << std::endl;
  if(Verbosity()>0) std::cout << "number of duplicates : " << n_dupli << std::endl;
  return cachedPositions;
}

int PHCASeeding::Process(PHCompositeNode */*topNode*/)
{
  if(_n_iteration>0){
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  t_seed->restart();

  _rtree.clear();
  PositionMap globalClusPositions = FillTree();
  t_seed->stop();
  if(Verbosity()>0) std::cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  int numberofseeds = 0;
  numberofseeds += FindSeeds(globalClusPositions);
  t_seed->stop();
  if(Verbosity()>0) std::cout << "number of seeds " << numberofseeds << std::endl;
  if(Verbosity()>0) std::cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::FindSeeds(const PositionMap& globalPositions)
{
  std::vector<pointKey> allClusters;
  QueryTree(_rtree,
            -2*M_PI, // phi
            -200, // eta
            _start_layer-0.5, // layer 
            2*M_PI, // phi
            200, // eta
            _end_layer+0.5, // layer
            allClusters);
  t_seed->stop();
  if(Verbosity()>0) std::cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  if(Verbosity()>1) std::cout << " number of clusters: " << allClusters.size() << std::endl;
  t_seed->restart();

  std::shared_ptr<TripletGraph> triplets = CreateTriplets(fromPointKey(allClusters), globalPositions);
  std::vector<keylist> trackSeedKeyLists = ConnectTriplets(triplets,globalPositions);
  std::vector<TrackSeed_v1> seeds = ConvertToSeeds(trackSeedKeyLists, globalPositions);
   
  publishSeeds(seeds);
  return seeds.size();
}

std::shared_ptr<TripletGraph> PHCASeeding::CreateTriplets(const std::vector<coordKey>& clusters, const PositionMap& globalPositions) const
{
  std::shared_ptr<TripletGraph> triplets = std::make_shared<TripletGraph>();
  for(int i = 0; i<clusters.size(); i++)
  {
    triplets->cluster_to_node.insert({clusters[i].second,i});
    triplets->clusters.push_back(clusters[i].second);
  }
  triplets->nodes.resize(clusters.size());

  double triplet_time = 0.;
  t_seed->stop();
  t_seed->restart();

  std::for_each (std::execution::par, triplets->cluster_to_node.begin(), triplets->cluster_to_node.end(), [&](std::pair<TrkrDefs::cluskey,size_t> entry)
  {
    size_t node_id = entry.second;
    coordKey StartCluster = clusters[node_id];
    // get clusters near this one in adjacent layers
    TrkrDefs::cluskey startKey = StartCluster.second;
    const double StartPhi = StartCluster.first[0];
    const double StartZ = StartCluster.first[1];
    unsigned int StartLayer = StartCluster.first[2];
    if(StartLayer < _start_layer || StartLayer > _end_layer) return;
    const auto& globalpos = globalPositions.at(startKey);
    const double StartX = globalpos(0);
    const double StartY = globalpos(1);

    if(Verbosity()>2)
    {
      std::cout << " starting cluster:" << std::endl;
      std::cout << " z: " << StartZ << std::endl;
      std::cout << " phi: " << StartPhi << std::endl;
      std::cout << " layer: " << StartLayer << std::endl;
    }

    std::vector<pointKey> ClustersAbove;
    std::vector<pointKey> ClustersBelow;
    //if(StartLayer==22 || StartLayer==23 || StartLayer==39 || StartLayer==40) return;
    int lower_offset = 1;
    int upper_offset = 1;
    if(StartLayer==21 || StartLayer==38) upper_offset = 2; 
    if(StartLayer==23 || StartLayer==40) lower_offset = 2;
    // gets clusters one layer below within a (phi,z) neighborhood
    QueryTree(_rtree,
              StartPhi - _neighbor_phi_width,
              StartZ - _neighbor_eta_width,
              (double) StartLayer - 0.5 - lower_offset,
              StartPhi + _neighbor_phi_width,
              StartZ + _neighbor_eta_width,
              (double) StartLayer - 0.5,
              ClustersBelow);
    // gets clusters one layer above within a (phi,z) neighborhood
    QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartZ-_neighbor_eta_width,
              (double) StartLayer + 0.5,
              StartPhi+_neighbor_phi_width,
              StartZ+_neighbor_eta_width,
              (double) StartLayer + 0.5 + upper_offset,
              ClustersAbove);

    if(lower_offset==1 && ClustersBelow.size()==0)
    {
          ClustersBelow.clear();
          QueryTree(_rtree,
              StartPhi - _neighbor_phi_width,
              StartZ - _neighbor_eta_width,
              (double) StartLayer - 0.5 - lower_offset - 1,
              StartPhi + _neighbor_phi_width,
              StartZ + _neighbor_eta_width,
              (double) StartLayer - 0.5,
              ClustersBelow);
    }
    else if(upper_offset==1 && ClustersAbove.size()==0)
    {
          ClustersAbove.clear();
          QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartZ-_neighbor_eta_width,
              (double) StartLayer + 0.5,
              StartPhi+_neighbor_phi_width,
              StartZ+_neighbor_eta_width,
              (double) StartLayer + 0.5 + upper_offset + 1,
              ClustersAbove);
    }

    if(Verbosity()>2) 
    {
      std::cout << " entries in below layer: " << ClustersBelow.size() << std::endl;
      std::cout << " entries in above layer: " << ClustersAbove.size() << std::endl;
    }

    // calculate (dx, dy, dz) vector for each neighboring cluster
/*
    std::vector<std::array<double,3>> delta_below;
    std::vector<std::array<double,3>> delta_above;
*/

    std::vector<std::array<double,2>> deltaUV_below;
    std::vector<std::array<double,2>> deltaSZ_below;
    std::vector<std::array<double,2>> deltaUV_above;
    std::vector<std::array<double,2>> deltaSZ_above;

/*
    delta_below.clear();
    delta_above.clear();
    delta_below.resize(ClustersBelow.size());
    delta_above.resize(ClustersAbove.size());
*/
    const double StartS2 = (StartX*StartX+StartY*StartY);
    const double StartU = (StartX)/StartS2;
    const double StartV = (StartY)/StartS2;
    const double StartS = sqrt(StartS2);
/*
    std::transform(ClustersBelow.begin(), ClustersBelow.end(), delta_below.begin(), [&](pointKey BelowCandidate)
    {
      const auto& belowpos = globalPositions.at(BelowCandidate.second);
      const double belowS2 = belowpos(0)*belowpos(0)+belowpos(1)*belowpos(1);
      const double belowU = belowpos(0)/belowS2;
      const double belowV = belowpos(1)/belowS2;
      return std::array<double,3> { belowU-StartU, belowV-StartV, belowpos(2)-StartZ };
    });

    std::transform(ClustersAbove.begin(), ClustersAbove.end(), delta_above.begin(), [&](pointKey AboveCandidate)
    {
      const auto& abovepos = globalPositions.at(AboveCandidate.second);
      const double aboveS2 = abovepos(0)*abovepos(0)+abovepos(1)*abovepos(1);
      const double aboveU = abovepos(0)/aboveS2;
      const double aboveV = abovepos(1)/aboveS2;
      return std::array<double,3> { aboveU-StartU, aboveV-StartV, abovepos(2)-StartZ };
    });
*/

    for(pointKey cluster : ClustersBelow)
    {
      const auto& globalpos = globalPositions.at(cluster.second);
      const double cx = globalpos(0);
      const double cy = globalpos(1);
      const double S2 = cx*cx+cy*cy;
      const double U = cx/S2;
      const double V = cy/S2;
      deltaUV_below.push_back(std::array<double,2>{U-StartU,V-StartV});
      deltaSZ_below.push_back(std::array<double,2>{globalpos(2)-StartZ,sqrt(S2)-StartS});
    }


    for(pointKey cluster : ClustersAbove)
    {
      const auto& globalpos = globalPositions.at(cluster.second);
      const double cx = globalpos(0);
      const double cy = globalpos(1);
      const double S2 = cx*cx+cy*cy;
      const double U = cx/S2;
      const double V = cy/S2;
      deltaUV_above.push_back(std::array<double,2>{U-StartU,V-StartV});
      deltaSZ_above.push_back(std::array<double,2>{globalpos(2)-StartZ,sqrt(S2)-StartS});
    }

    // find all triplets that are "straight enough"
    // (in other words, all triplets for which the angle between the (dx, dy, dz) vectors is close enough to 180 degrees)

    double maxCosUVPlaneAngle = -0.95;
    double maxCosSZPlaneAngle = -0.95;

    std::vector<coordKey> bestBelowClusters;
    std::vector<coordKey> bestAboveClusters;
    //bestBelowClusters.resize(1);
    //bestAboveClusters.resize(1);

    for(size_t iAbove = 0; iAbove<deltaUV_above.size(); ++iAbove)
    {
      for(size_t iBelow = 0; iBelow<deltaUV_below.size(); ++iBelow)
      {
/*
        double angle = breaking_angle(
          delta_below[iBelow][0],
          delta_below[iBelow][1],
          delta_below[iBelow][2],
          delta_above[iAbove][0],
          delta_above[iAbove][1],
          delta_above[iAbove][2]);
*/

        double UVangle = breaking_angle(
          deltaUV_below[iBelow][0],
          deltaUV_below[iBelow][1],
          deltaUV_above[iAbove][0],
          deltaUV_above[iAbove][1]);
        double SZangle = breaking_angle(
          deltaSZ_below[iBelow][0],
          deltaSZ_below[iBelow][1],
          deltaSZ_above[iAbove][0],
          deltaSZ_above[iAbove][1]);

        if(UVangle < maxCosUVPlaneAngle && SZangle < maxCosSZPlaneAngle)
        {
          //maxCosPlaneAngle = cos(angle);
          //bestBelowClusters[0] = fromPointKey(ClustersBelow[iBelow]);
          //bestAboveClusters[0] = fromPointKey(ClustersAbove[iAbove]);
          //int layer_index = StartLayer - (_nlayers_intt + _nlayers_maps);
          
          triplets->nodes[node_id].push_back({ 
            triplets->cluster_to_node[ClustersBelow[iBelow].second],
            node_id,
            triplets->cluster_to_node[ClustersAbove[iAbove].second],
            startKey,
            false
          });
          
          //triplets[layer_index].push_back({ fromPointKey(ClustersBelow[iBelow]), StartCluster, fromPointKey(ClustersAbove[iAbove]) });
        }
      }
    }
/*    
    if(ClustersBelow.size()>0 && ClustersAbove.size()>0) triplets->nodes[node_id].push_back({
      triplets->cluster_to_node[bestBelowClusters[0].second],
      node_id,
      triplets->cluster_to_node[bestAboveClusters[0].second],
      startKey,
      false
    });
*/    
//    if(Verbosity()>2) std::cout << " max cos(vector angle): " << maxCosPlaneAngle << std::endl;
    
  });
  t_seed->stop();
  triplet_time += t_seed->elapsed();
  if(Verbosity()>0)
  {
    std::cout << "triplet forming time: " << triplet_time / 1000 << " s" << std::endl;
  }

  return triplets;
}
/*
double PHCASeeding::getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PositionMap& globalPositions) const
{
  // Menger curvature = 1/R for circumcircle of triangle formed by most recent three clusters
  // We use here 1/R = 2*sin(breaking angle)/(hypotenuse of triangle)
  auto& a_pos = globalPositions.at(a);
  auto& b_pos = globalPositions.at(b);
  auto& c_pos = globalPositions.at(c);
  double hypot_length = sqrt(square<double>(c_pos.x()-a_pos.x())+square<double>(c_pos.y()-a_pos.y())+square<double>(c_pos.z()-a_pos.z()));
  double break_angle = breaking_angle(
    a_pos.x()-b_pos.x(),
    a_pos.y()-b_pos.y(),
    a_pos.z()-b_pos.z(),
    c_pos.x()-b_pos.x(),
    c_pos.y()-b_pos.y(),
    c_pos.z()-b_pos.z());
  return 2*sin(break_angle)/hypot_length;
}
*/
std::vector<keylist> PHCASeeding::ConnectTriplets(const std::shared_ptr<TripletGraph> triplets, const PositionMap& globalPositions) const
{
/*
  // follow bidirectional links to form lists of cluster keys
  // (to be fitted for track seed parameters)
  std::vector<keylist> trackSeedPairs;
  // get starting cluster keys, create a keylist for each
  // (only check last element of each pair because we start from the outer layers and go inward)
  for(unsigned int layer = 0; layer < _nlayers_tpc-1; ++layer)
  {
    for(auto startCand = bidirectionalLinks[layer].begin(); startCand != bidirectionalLinks[layer].end(); ++startCand)
    {
      bool has_above_link = false;
      unsigned int imax = 1;
      if(layer==_nlayers_tpc-2) imax = 1;
      for(unsigned int i=1;i<=imax;i++)
      {
        has_above_link = has_above_link || std::any_of(bidirectionalLinks[layer+i].begin(),bidirectionalLinks[layer+i].end(),[&](keylink k){return (*startCand)[0]==k[1];});
      }
//      for(std::vector<keylink>::iterator testlink = bidirectionalLinks[layer+1].begin(); !has_above_link && (testlink != bidirectionalLinks[layer+1].end()); ++testlink)
//      {
//        if((*startCand) == (*testlink)) continue;
//        if((*startCand)[0] == (*testlink)[1]) has_above_link = true;
//      } 
      if(!has_above_link)
      {
        trackSeedPairs.push_back({(*startCand)[0].second,(*startCand)[1].second});
      }
    }
  }

  // form all possible starting 3-cluster tracks (we need that to calculate curvature)
  std::vector<keylist> trackSeedKeyLists;
  for(auto& trackKeyChain : trackSeedPairs)
  {
    TrkrDefs::cluskey trackHead = trackKeyChain.back();
    unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt+_nlayers_maps);
    for(auto& testlink : bidirectionalLinks[trackHead_layer])
    {
      if(testlink[0].second==trackHead)
      {
        keylist trackSeedTriplet;
        trackSeedTriplet.push_back(trackKeyChain[0]);
        trackSeedTriplet.push_back(trackKeyChain[1]);
        trackSeedTriplet.push_back(testlink[1].second);
        trackSeedKeyLists.push_back(trackSeedTriplet);
      }
    }
  }

  t_seed->stop();
  if(Verbosity()>0) std::cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)

  // std::cout << "STARTING SEED ASSEMBLY" << std::endl;
*/
  t_seed->stop();
  t_seed->restart();
  std::vector<keylist> seedEnds;

  for(std::vector<Triplet> node : triplets->nodes)
  {
    for(Triplet t : node)
    {
      std::vector<Triplet> lower_linked_node = triplets->nodes[t.lower_node];
      if(lower_linked_node.size()==0)
      {
        TrkrDefs::cluskey lower_key = triplets->clusters[t.lower_node];
        TrkrDefs::cluskey upper_key = triplets->clusters[t.upper_node];
        seedEnds.push_back({
          lower_key,
          t.center_cluster,
          upper_key});
        //t.used = false;
      }
    }
  }

  std::vector<keylist> pendingSeeds;
  for(auto& first_triplet : seedEnds)
  {
    pendingSeeds.push_back(first_triplet);
  }
  
  std::vector<keylist> trackSeedKeyLists;
  while(pendingSeeds.size()>0)
  {
    //std::cout << "pending seed size: " << pendingSeeds.size() << std::endl;
    std::vector<keylist> newPendingSeeds;
    newPendingSeeds.clear();
    for(auto& seed : pendingSeeds)
    {
      bool no_next_link = true;
      TrkrDefs::cluskey lastCluster = seed.back();
      //std::cout << "last cluster layer: " << std::to_string(TrkrDefs::getLayer(lastCluster)) << std::endl;
      TrkrDefs::cluskey secondLastCluster = seed.rbegin()[1];
      //std::cout << "second last cluster layer: " << std::to_string(TrkrDefs::getLayer(secondLastCluster)) << std::endl;
      size_t last_nodeid = triplets->cluster_to_node[lastCluster];
      //std::cout << "last node id: " << last_nodeid << std::endl;
      size_t secondlast_nodeid = triplets->cluster_to_node[secondLastCluster];
      //std::cout << "second last node id: " << secondlast_nodeid << std::endl;
      for(auto& trp : triplets->nodes[last_nodeid])
      {
        if(trp.lower_node == secondlast_nodeid && trp.center_node == last_nodeid)
        {
/*
          auto last_pos = globalPositions.at(lastCluster);
          auto secondlast_pos = globalPositions.at(secondLastCluster);
          double x1 = last_pos.x();
          double y1 = last_pos.y();
          double z1 = last_pos.z();
          double x2 = secondlast_pos.x();
          double y2 = secondlast_pos.y();
          double z2 = secondlast_pos.z();
          //std::cout << "last cluster position: (" << x1 << ", " << y1 << ", " << z1 << ")" << std::endl;
          //std::cout << "second last cluster position: (" << x2 << ", " << y2 << ", " << z2 << ")" << std::endl;
          double dr_12 = sqrt(x1*x1+y1*y1)-sqrt(x2*x2+y2*y2);
          auto& test_pos = globalPositions.at(triplets->clusters[trp.upper_node]);
          //std::cout << "test node id: " << trp.lower_node << std::endl;
          //std::cout << "test layer: " << std::to_string(TrkrDefs::getLayer(triplets->clusters[trp.lower_node])) << std::endl;
          double xt = test_pos.x();
          double yt = test_pos.y();
          double zt = test_pos.z();
          //std::cout << "test position: (" << xt << ", " << yt << ", " << zt << ")" << std::endl;
          double new_dr = sqrt(xt*xt+yt*yt)-sqrt(x1*x1+y1*y1);
          //std::cout << "dz/dr current: " << (z1-z2)/dr_12 << std::endl;
          //std::cout << "dz/dr new: " << (zt-z1)/new_dr << std::endl;
          
          //if(fabs( (z1-z2)/dr_12 - (zt-z1)/new_dr )>0.5) continue;
          
          //auto& third_pos = globalPositions.at(seed.rbegin()[2]);
          //double x3 = third_pos.x();
          //double y3 = third_pos.y();
          //std::cout << "third last cluster layer: " << std::to_string(TrkrDefs::getLayer(seed.rbegin()[2])) << std::endl;
          //std::cout << "third last cluster position: (" << x3 << ", " << y3 << ", " << third_pos.z() << ")" << std::endl;
          //double dr_23 = sqrt(x2*x2+y2*y2)-sqrt(x3*x3+y3*y3);
          //double phi1 = atan2(y1,x1);
          //double phi2 = atan2(y2,x2);
          //double phi3 = atan2(y3,x3);
          //double dphi12 = phi1-phi2;
          //double dphi23 = phi2-phi3;
          //double d2phidr2 = dphi12/(dr_12*dr_12)-dphi23/(dr_23*dr_23);
          //std::cout << "d2phi/dr2 current: " << d2phidr2 << std::endl;
          //double new_dphi = atan2(yt,xt)-atan2(y1,x1);
          //double new_d2phidr2 = new_dphi/(new_dr*new_dr)-dphi12/(dr_12*dr_12);
          //std::cout << "d2phi/dr2 new: " << new_d2phidr2 << std::endl;
          const double s2_1 = x1*x1+y1*y1;
          const double s2_2 = x2*x2+y2*y2;
          const double s2_t = xt*xt+yt*yt;
          const double u1 = x1/s2_1;
          const double u2 = x2/s2_2;
          const double ut = xt/s2_t;
          const double v1 = y1/s2_1;
          const double v2 = y2/s2_2;
          const double vt = yt/s2_t;
          //std::cout << "dv/du " << (vt-v1)/(ut-u1) << " -> " << (v1-v2)/(u1-u2) << std::endl;
          //if(fabs((vt-v1)/(ut-u1)-(v1-v2)/(u1-u2))<0.5)
*/
          {
            //std::cout << "added" << std::endl;
            no_next_link = false;
            keylist newseed = seed;
            newseed.push_back(triplets->clusters[trp.upper_node]);
            newPendingSeeds.push_back(newseed);
            trp.used = true;
          }
        }
      }
      if(no_next_link && seed.size()>_min_clusters_per_track)
      {
        trackSeedKeyLists.push_back(seed);
      }
      //std::cout << "----------" << std::endl;
    }
    //std::cout << "new pending seed size: " << newPendingSeeds.size() << std::endl;
    //std::cout << "final seed size: " << trackSeedKeyLists.size() << std::endl;
    pendingSeeds = newPendingSeeds;
  }

/*  for(auto& layer : triplets)
  {
    for(auto& trp : layer)
    {
      tempSeedKeyLists.push_back({trp[0].second, trp[1].second, trp[2].second});
    }
  }
*/
//  std::vector<keylist> tempSeedKeyLists = trackSeedKeyLists;
//  trackSeedKeyLists.clear();

/*  for(int layer_idx = 0; layer_idx<triplets.size(); layer_idx++)
  {
    if(Verbosity()>0) std::cout << "temp size: " << tempSeedKeyLists.size() << std::endl;
    if(Verbosity()>0) std::cout << "final size: " << trackSeedKeyLists.size() << std::endl;
    std::vector<keylist> newTempSeedKeyLists;
    std::vector<bool> usedTriplet;
    std::vector<triplet> triplets_thislayer = triplets[layer_idx];

    usedTriplet.resize(triplets_thislayer.size());
    std::fill(usedTriplet.begin(),usedTriplet.end(),false);

    std::for_each( std::execution::par, tempSeedKeyLists.begin(), tempSeedKeyLists.end(), [&](keylist seed)
    {
      bool no_next_link = true;
      TrkrDefs::cluskey lastCluster = seed.back();
      TrkrDefs::cluskey secondLastCluster = seed.rbegin()[1];
      for(int trp_idx = 0; trp_idx<triplets_thislayer.size(); trp_idx++)
      {
        triplet trp = triplets_thislayer[trp_idx];
        if(trp[1].second != lastCluster || trp[0].second != secondLastCluster) continue;
        auto last_pos = globalPositions.at(lastCluster);
        auto secondlast_pos = globalPositions.at(secondLastCluster);
        double x1 = last_pos.x();
        double y1 = last_pos.y();
        double z1 = last_pos.z();
        double x2 = secondlast_pos.x();
        double y2 = secondlast_pos.y();
        double z2 = secondlast_pos.z();
        double dr_12 = sqrt(x1*x1+y1*y1)-sqrt(x2*x2+y2*y2);
        auto& test_pos = globalPositions.at(trp[2].second);
        double xt = test_pos.x();
        double yt = test_pos.y();
        double zt = test_pos.z();
        double new_dr = sqrt(xt*xt+yt*yt)-sqrt(x1*x1+y1*y1);
        if(fabs( (z1-z2)/dr_12 - (zt-z1)/new_dr )>1e9) continue;
        auto& third_pos = globalPositions.at(seed.rbegin()[2]);
        double x3 = third_pos.x();
        double y3 = third_pos.y();
        double dr_23 = sqrt(x2*x2+y2*y2)-sqrt(x3*x3+y3*y3);
        double phi1 = atan2(y1,x1);
        double phi2 = atan2(y2,x2);
        double phi3 = atan2(y3,x3);
        double dphi12 = std::fmod(phi1-phi2,M_PI);
        double dphi23 = std::fmod(phi2-phi3,M_PI);
        double d2phidr2 = dphi12/(dr_12*dr_12)-dphi23/(dr_23*dr_23);
        double new_dphi = std::fmod(atan2(yt,xt)-atan2(y1,x1),M_PI);
        double new_d2phidr2 = new_dphi/(new_dr*new_dr)-dphi12/(dr_12*dr_12);
        if(fabs(d2phidr2-new_d2phidr2)<1e9)
        {
          no_next_link = false;
          keylist newseed = seed;
          newseed.push_back(trp[2].second);
          mtxn.lock();
          newTempSeedKeyLists.push_back(newseed);
          mtxn.unlock();
          usedTriplet[trp_idx] = true;
        }
      }
      if(no_next_link && seed.size()>=_min_clusters_per_track)
      {
        mtx.lock();
        trackSeedKeyLists.push_back(seed);
        mtx.unlock();
      }
         
    });

    for(int trp_idx = 0; trp_idx<triplets_thislayer.size(); trp_idx++)
    {
      if(!usedTriplet[trp_idx])
      {
        triplet trp = triplets_thislayer[trp_idx];
        keylist tkeys;
        tkeys.push_back(trp[0].second);
        tkeys.push_back(trp[1].second);
        tkeys.push_back(trp[2].second);
        newTempSeedKeyLists.push_back(tkeys);
      }
    }
    if(Verbosity()>0) std::cout << "new temp size: " << newTempSeedKeyLists.size() << std::endl;
    tempSeedKeyLists = newTempSeedKeyLists;
  }
*/
/*
  while(tempSeedKeyLists.size()>0)
  {
    if(Verbosity()>0) std::cout << "temp size: " << tempSeedKeyLists.size() << std::endl;
    if(Verbosity()>0) std::cout << "final size: " << trackSeedKeyLists.size() << std::endl;
    std::vector<keylist> newtempSeedKeyLists;
    std::for_each( std::execution::par, tempSeedKeyLists.begin(), tempSeedKeyLists.end(), [&](keylist seed)
    {
      TrkrDefs::cluskey trackHead = seed.back();
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead)-(_nlayers_intt+_nlayers_maps);
      if(trackHead_layer==_nlayers_tpc-1) return;
      bool no_next_link = true;
      for(auto& trp : triplets[trackHead_layer])
      {
        if(trp[1].second != trackHead) continue;
        if(trp[0].second != seed.rbegin()[1]) continue;

        auto& head_pos = globalPositions.at(trackHead);
        auto& prev_pos = globalPositions.at(seed.rbegin()[1]);
        double x1 = head_pos.x();
        double y1 = head_pos.y();
        double z1 = head_pos.z();
        double x2 = prev_pos.x();
        double y2 = prev_pos.y();
        double z2 = prev_pos.z();
        double dr_12 = sqrt(x1*x1+y1*y1)-sqrt(x2*x2+y2*y2);
        TrkrDefs::cluskey testCluster = trp[2].second;
        auto& test_pos = globalPositions.at(testCluster);
        double xt = test_pos.x();
        double yt = test_pos.y();
        double zt = test_pos.z();
        double new_dr = sqrt(xt*xt+yt*yt)-sqrt(x1*x1+y1*y1);
        if(fabs( (z1-z2)/dr_12 - (zt-z1)/new_dr )>0.5) continue;
        auto& third_pos = globalPositions.at(seed.rbegin()[2]);
        double x3 = third_pos.x();
        double y3 = third_pos.y();
        double dr_23 = sqrt(x2*x2+y2*y2)-sqrt(x3*x3+y3*y3);
        double phi1 = atan2(y1,x1);
        double phi2 = atan2(y2,x2);
        double phi3 = atan2(y3,x3);
        double dphi12 = std::fmod(phi1-phi2,M_PI);
        double dphi23 = std::fmod(phi2-phi3,M_PI);
        double d2phidr2 = dphi12/(dr_12*dr_12)-dphi23/(dr_23*dr_23);
        double new_dphi = std::fmod(atan2(yt,xt)-atan2(y1,x1),M_PI);
        double new_d2phidr2 = new_dphi/(new_dr*new_dr)-dphi12/(dr_12*dr_12);

        if(fabs(d2phidr2-new_d2phidr2)<.005)
        {
          no_next_link = false;
          keylist newseed = seed;
          newseed.push_back(trp[2].second);
          mtxn.lock();
          newtempSeedKeyLists.push_back(newseed);
          mtxn.unlock();
        }
      }
      if(no_next_link && seed.size()>5)
      {
        mtx.lock();
        trackSeedKeyLists.push_back(seed);
        mtx.unlock();
      }
    });
    if(Verbosity()>0) std::cout << "new temp size: " << newtempSeedKeyLists.size() << std::endl;
    tempSeedKeyLists = newtempSeedKeyLists;
  }


//  trackSeedKeyLists = tempSeedKeyLists;

  for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    bool reached_end = false;
    while(!reached_end)
    {
      TrkrDefs::cluskey trackHead = trackKeyChain->back();
      TrkrDefs::cluskey secondToLast = trackKeyChain->rbegin()[1];
      TrkrDefs::cluskey thirdToLast = trackKeyChain->rbegin()[2];
      auto& head_pos = globalPositions.at(trackHead);
      auto& sec_pos = globalPositions.at(secondToLast);
      auto& third_pos = globalPositions.at(thirdToLast);
      double dz_avg = ((head_pos.z()-sec_pos.z())+(sec_pos.z()-third_pos.z()))/2.;
      double dx1 = head_pos.x()-sec_pos.x();
      double dy1 = head_pos.y()-sec_pos.y();
      double dx2 = sec_pos.x()-third_pos.x();
      double dy2 = sec_pos.y()-third_pos.y();
      double ddx = dx1-dx2;
      double ddy = dy1-dy2;
      double new_dx = dx1+ddx;
      double new_dy = dy1+ddy;
      double new_x = head_pos.x()+new_dx;
      double new_y = head_pos.y()+new_dy;
      double new_z = head_pos.z()+dz_avg;
      std::cout << "(x,y,z) = (" << head_pos.x() << ", " << head_pos.y() << ", " << head_pos.z() << ")" << std::endl;
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
      std::cout << "layer " << trackHead_layer << std::endl;
      std::cout << "projected: (" << new_x << ", " << new_y << ", " << new_z << ")" << std::endl;
      TrkrDefs::cluskey nextCluster;
      double bestDist = 1e9;
      bool no_next_link = true;
      for(auto testlink = bidirectionalLinks[trackHead_layer].begin(); testlink != bidirectionalLinks[trackHead_layer].end(); ++testlink)
      {
        if((*testlink)[0].second==trackHead)
        {
          TrkrDefs::cluskey testCluster = (*testlink)[1].second;
          auto& test_pos = globalPositions.at(testCluster);
          std::cout << "test cluster: (" << test_pos.x() << ", " << test_pos.y() << ", " << test_pos.z() << ")" << std::endl;
          double distToNew = sqrt(square<double>(test_pos.x()-new_x)+square<double>(test_pos.y()-new_y)+square<double>(test_pos.z()-new_z));
          if(distToNew<bestDist)
          {
            std::cout << "current best" << std::endl;
            nextCluster = testCluster;
            bestDist = distToNew;
          }
          no_next_link = false;
        }
      }
      if(!no_next_link) trackKeyChain->push_back(nextCluster);
      if(no_next_link) reached_end = true;
    }
  }
*/
  t_seed->stop();
  if(Verbosity()>0) std::cout << "keychain assembly time: " << t_seed->elapsed() / 1000 << " s" << std::endl;
  t_seed->restart();
  if(Verbosity()>0) std::cout << "track key chains assembled: " << trackSeedKeyLists.size() << std::endl;
  if(Verbosity()>2)
  {
    std::cout << "track key chain lengths: " << std::endl;
    for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
    {
      std::cout << " " << trackKeyChain->size() << std::endl;
    }
  }
/*
  int jumpcount = 0;
  if(Verbosity()>1)
  {
    std::cout << " track key associations:" << std::endl;
    for(size_t i=0;i<trackSeedKeyLists.size();++i)
    {
      std::cout << " seed " << i << ":" << std::endl;

      double lasteta = -100;
      double lastphi = -100;
      for(size_t j=0;j<trackSeedKeyLists[i].size();++j)
      {
        const auto& globalpos = globalPositions.at(trackSeedKeyLists[i][j]);           
        const double clus_phi = get_phi( globalpos );
        const double clus_eta = get_eta( globalpos );
        const double etajump = clus_eta-lasteta;
        const double phijump = clus_phi-lastphi;
        unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j].second);
        if((fabs(etajump)>0.1 && lasteta!=-100) || (fabs(phijump)>1 && lastphi!=-100))
	      {
           std::cout << " Eta or Phi jump too large! " << std::endl;
           ++jumpcount;
        }
        std::cout << " (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
          " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl;
      }
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
    std::cout << " Total large jumps: " << jumpcount << std::endl;
  }
  t_seed->stop();
  if(Verbosity()>0) std::cout << "eta-phi sanity check time: " << t_seed->elapsed() / 1000 << " s" << std::endl;
  t_seed->restart();
*/
  std::vector<keylist> final_keylists;
  if(Verbosity()>0) std::cout << "SEED LIST" << std::endl;
  for(auto& kl : trackSeedKeyLists)
  {
    if(Verbosity()>0)
    {
      for(auto& c : kl)
      {
        const auto& globalpos = globalPositions.at(c);
        const unsigned int layer = TrkrDefs::getLayer(c);
        std::cout << "layer " << layer << ", (x, y, z)=(" << globalpos(0) << ", " << globalpos(1) << ", " << globalpos(2) << ")" << std::endl;
      }
      std::cout << "----------------" << std::endl;
    }
    final_keylists.push_back(kl);
  }
  return final_keylists;
}

std::vector<TrackSeed_v1> PHCASeeding::ConvertToSeeds(const std::vector<keylist>& chains, const PositionMap& globalPositions) const
{
  if(Verbosity()>0) std::cout << "removing bad clusters" << std::endl;
  std::vector<TrackSeed_v1> clean_chains;

  for(const auto& chain : chains)
  {
    if(chain.size()<3) continue;
    if(Verbosity()>0) std::cout << "chain size: " << chain.size() << std::endl;

    TrackFitUtils::position_vector_t xy_pts;
    for( const auto& ckey:chain )
    {
      const auto &global = globalPositions.at(ckey);
      xy_pts.emplace_back( global.x(), global.y() );
    } 

    // fit a circle through x,y coordinates
    //const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin( xy_pts );

    // skip chain entirely if fit fails
    //if( std::isnan( R ) ) continue;

    // calculate residuals
    //const std::vector<double> xy_resid = fitter->GetCircleClusterResiduals(xy_pts,R,X0,Y0);
    
    // assign clusters to seed
    TrackSeed_v1 trackseed;
    for(size_t i=0;i<chain.size();i++)
    {
      //if(xy_resid[i]>_xy_outlier_threshold) continue;
      trackseed.insert_cluster_key(chain.at(i));
    }

    trackseed.circleFitByTaubin(_cluster_map,tGeometry,7,55);
    trackseed.lineFit(_cluster_map,tGeometry,7,55);

    clean_chains.push_back(trackseed);
    if(Verbosity()>0) std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
  }

//  std::vector<double> chi2;
//  clean_chains = fitter->ALICEKalmanFilter(chains,false,globalPositions,chi2);
  return clean_chains;
}


void PHCASeeding::publishSeeds(const std::vector<TrackSeed_v1>& seeds)
{
 
  for( const auto&  seed:seeds )
  {
    auto pseed = std::make_unique<TrackSeed_v1>(seed);
    if(Verbosity() > 4)
      { pseed->identify(); }
    _track_map->insert(pseed.get());
  }
}

int PHCASeeding::Setup(PHCompositeNode *topNode)
{
  if(Verbosity()>0) std::cout << "Called Setup" << std::endl;
  if(Verbosity()>0) std::cout << "topNode:" << topNode << std::endl;
  PHTrackSeeding::Setup(topNode);
  
  // geometry initialization
  int ret = InitializeGeometry(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    { return ret; }
    
  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( m_dcc )
  { std::cout << "PHCASeeding::Setup - found static TPC distortion correction container" << std::endl; }
  
  t_fill = std::make_unique<PHTimer>("t_fill");
  t_seed = std::make_unique<PHTimer>("t_seed");
  t_fill->stop();
  t_seed->stop();
  PHFieldConfigv1 fcfg;
  fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);
  char *calibrationsroot = getenv("CALIBRATIONROOT");
  assert(calibrationsroot);
  auto magField = std::string(calibrationsroot) +
    std::string("/Field/Map/sphenix3dtrackingmapxyz.root"); 
  fcfg.set_filename(magField);
  //  fcfg.set_rescale(1);
  std::unique_ptr<PHField> field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));

  fitter = std::make_unique<ALICEKF>(topNode,_cluster_map,field_map.get(),_fieldDir,_min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if(Verbosity()>0) std::cout << "Called End " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
