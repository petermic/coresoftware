#ifndef __TrkrClusterv1_H__
#define __TrkrClusterv1_H__

#include "TrkrCluster.h"
#include "TrkrDefUtil.h"

#include <limits.h>
#include <cmath>
#include <iostream>
#include <set>

class TrkrClusterv1 : public TrkrCluster
{
 public:
  TrkrClusterv1();
  virtual ~TrkrClusterv1() {}
  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const;
  void Reset() {}
  int isValid() const;
  TrkrCluster* Clone() const { return new TrkrClusterv1(*this); }


  TrkrDefs::cluskey get_id() const { return _id; }
  void set_id(TrkrDefs::cluskey id) { _id = id; }
  //
  // cluster position
  //
  float get_x() const { return _pos[0]; }
  void set_x(float x) { _pos[0] = x; }
  float get_y() const { return _pos[1]; }
  void set_y(float y) { _pos[1] = y; }
  float get_z() const { return _pos[2]; }
  void set_z(float z) { _pos[2] = z; }
  float get_position(int coor) const { return _pos[coor]; }
  void set_position(int coor, float xi) { _pos[coor] = xi; }
  void set_global() { _is_global = true; }
  void set_local() { _is_global = false; }
  bool is_global() { return _is_global; }
  // 
  // cluster info
  //
  unsigned int get_adc() const { return _adc; }
  void set_adc(unsigned int adc) { _adc = adc; }
  float get_size(unsigned int i, unsigned int j) const;        //< get cluster dimension covar
  void set_size(unsigned int i, unsigned int j, float value);  //< set cluster dimension covar

  float get_error(unsigned int i, unsigned int j) const;        //< get cluster error covar
  void set_error(unsigned int i, unsigned int j, float value);  //< set cluster error covar

  //
  // clustered hit ids methods
  //
  void clear_hits() { _hit_ids.clear(); }
  bool empty_hits() { return _hit_ids.empty(); }
  size_t size_hits() { return _hit_ids.size(); }
  void insert_hit(TrkrDefs::hitsetkey hit_id) { _hit_ids.insert(hit_id); }
  size_t erase_hit(TrkrDefs::hitsetkey hit_id) { return _hit_ids.erase(hit_id); }
  ConstHitIter begin_hits() const { return _hit_ids.begin(); }
  ConstHitIter find_hit(TrkrDefs::hitsetkey hitid) const { return _hit_ids.find(hitid); }
  ConstHitIter end_hits() const { return _hit_ids.end(); }
  HitIter begin_hits() { return _hit_ids.begin(); }
  HitIter find_hit(TrkrDefs::hitsetkey hitid) { return _hit_ids.find(hitid); }
  HitIter end_hits() { return _hit_ids.end(); }
  // convenience interface

  float get_phi_size() const;
  float get_z_size() const;

  float get_rphi_error() const;
  float get_phi_error() const;
  float get_z_error() const;

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  TrkrDefs::cluskey _id; //< unique identifier within container
  float _pos[3];         //< mean position x,y,z
  bool _is_global;       //< flag for coord sys (true = global)
  unsigned int _adc;     //< cluster sum adc (D. McGlinchey - Do we need this)
  float _size[6];        //< size covariance matrix (packed storage) (+/- cm^2)
  float _err[6];         //< covariance matrix: rad, arc and z
  HitSet _hit_ids;       //< list of cell hit ids

  ClassDef(TrkrClusterv1, 1);
};

#endif
