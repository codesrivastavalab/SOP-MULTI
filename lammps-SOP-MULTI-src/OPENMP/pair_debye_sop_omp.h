/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)

   Modified by : Krishnakanth B, Theoretical Biophysics Laboratory, 
   Molecular Biophysics Unit, Indian Institute of Science, Bangalore - 560012
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(debye/sop/omp,PairDebye_SOP_OMP);
// clang-format on
#else

#ifndef LMP_PAIR_DEBYE_SOP_OMP_H
#define LMP_PAIR_DEBYE_SOP_OMP_H

#include "pair_debye_sop.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairDebye_SOP_OMP : public PairDebye_SOP, public ThrOMP {

 public:
  PairDebye_SOP_OMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
