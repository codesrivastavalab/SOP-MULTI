// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "pair_debye_sop_omp.h"
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
//#include <iostream>


#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDebye_SOP_OMP::PairDebye_SOP_OMP(LAMMPS *lmp) :
  PairDebye_SOP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
//  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairDebye_SOP_OMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairDebye_SOP_OMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,rinv,forcecoul,factor_coul, screening;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int i_dom_id, j_dom_id, i_mol_id, j_mol_id;  

  ecoul = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  const int * _noalias const domain_id = atom->domain_id; // domain tag of the atom
  const int * _noalias const molecule_id = atom->molecule; // molecule tag of the atom  

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    i_dom_id  = domain_id[i];
    i_mol_id = molecule_id[i];     
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      j_dom_id = domain_id[j];
      j_mol_id = molecule_id[j];

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq[itype][jtype]) {
      
        if((i_dom_id == 0) || (j_dom_id == 0) || (i_dom_id != j_dom_id) || (i_mol_id != j_mol_id)){
          r2inv = 1.0/rsq;
          r = sqrt(rsq);
          rinv = 1.0/r;
        screening = exp(-kappa[itype][jtype]*r);
        forcecoul = qqrd2e * qtmp*q[j] * epsilon_coul[itype][jtype] *screening * (kappa[itype][jtype] + rinv);
        forcecoul *= factor_coul;
        } else {forcecoul = 0.0;}
        fpair = forcecoul * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG){
        if((i_dom_id == 0) || (j_dom_id == 0) || (i_dom_id != j_dom_id) || (i_mol_id != j_mol_id)){
          ecoul = factor_coul * qqrd2e * epsilon_coul[itype][jtype] * qtmp*q[j]*rinv*screening;
//           std::cout << screening << " "<< epsilon_coul[itype][jtype]  << " "<< factor_coul <<std::endl;
          }else{ecoul = 0.0 ;}
        }
        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 0.0,ecoul,fpair,delx,dely,delz,thr);
      }
      }

   f[i].x += fxtmp;
   f[i].y += fytmp;
   f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDebye_SOP_OMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDebye_SOP::memory_usage();

  return bytes;
}
