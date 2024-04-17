/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Modified by : Krishnakanth B, Theoretical Biophysics Laboratory, 
   Molecular Biophysics Unit, Indian Institute of Science, Bangalore - 560012  
   
   Last Modified :  29 March 2023
------------------------------------------------------------------------- */   

#ifdef PAIR_CLASS
// clang-format off
PairStyle(list/sop,PairList_SOP);
// clang-format on
#else

#ifndef LMP_PAIR_LIST_SOP_H
#define LMP_PAIR_LIST_SOP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairList_SOP : public Pair {
 public:
  PairList_SOP(class LAMMPS *);
  ~PairList_SOP() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

 protected:
  void allocate();

  // potential specific parameters
  struct harm_p {
    double k, r0;
  };
  struct morse_p {
    double d0, alpha, r0;
  };
  struct ljsop_p {
    double epsilon, sigma;
  };
  struct repsop_p {
    double epsilon, sigma;
  };  
 

  union param_u {
    harm_p harm;
    morse_p morse;
    ljsop_p ljsop;
    repsop_p repsop;
  };

  struct list_param {
    int style;              // potential style indicator
    tagint id1, id2;        // global atom ids
    double cutsq;           // cutoff**2 for this pair
    double offset;          // energy offset
    union param_u param;    // parameters for style
  };

 protected:
  double cut_global;     // global cutoff distance
  list_param *params;    // lisf of pair interaction parameters
  int npairs;            // # of atom pairs in global list
  int check_flag;        // 1 if checking for missing pairs
};

}    // namespace LAMMPS_NS

#endif
#endif
