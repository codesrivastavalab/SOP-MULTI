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

#include "pair_lj_sop_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"
#include <iostream>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJ_SOP_OMP::PairLJ_SOP_OMP(LAMMPS *lmp) :
  PairLJ_SOP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
//  respa_enable = 0;
//  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJ_SOP_OMP::compute(int eflag, int vflag)
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

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJ_SOP_OMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const double * _noalias const special_lj = force->special_lj;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int * _noalias const domain_id = atom->domain_id; // domain tag of the atom
  const int * _noalias const molecule_id = atom->molecule; // molecule tag of the atom
  

  double xtmp,ytmp,ztmp,delx,dely,delz,fxtmp,fytmp,fztmp;
  double rsq,r2inv,r6inv,forcelj,factor_lj,evdwl,fpair;
  int connect_degree,i_dom_id, j_dom_id, i_mol_id, j_mol_id;
  
  const int nlocal = atom->nlocal;
  int j,jj,jnum,jtype;

  evdwl = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int    * _noalias const jlist = firstneigh[i];
    const double * _noalias const cutsqi = cutsq[itype];
    const double * _noalias const offseti = offset[itype];
    const double * _noalias const lj1i = lj1[itype];
    const double * _noalias const lj2i = lj2[itype];
    const double * _noalias const lj3i = lj3[itype];
    const double * _noalias const lj4i = lj4[itype];
    const double * _noalias const rep1i = rep1[itype];    
    const double * _noalias const rep2i = rep2[itype];    
    
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    jnum = numneigh[i];
    i_dom_id  = domain_id[i];
    i_mol_id = molecule_id[i];    
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      connect_degree = sbmask(j); // identify connectivity between beads i and j
      
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      j_dom_id = domain_id[j];
      j_mol_id = molecule_id[j];      

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
//          std::cout << rep1i[jtype]<< " " << ;        
        if((connect_degree == 1) || (connect_degree == 2)){
          forcelj =  r6inv*rep1i[jtype];
          }
          else{
            if ((i_dom_id == 0) || (j_dom_id == 0) || (i_mol_id != j_mol_id)||(i_dom_id != j_dom_id)){
	      forcelj = factor_lj*r6inv * (lj1i[jtype] * r6inv - lj2i[jtype]);
            }else{forcelj =0.0;} // end of lj non-bonded calculation          
          }
          fpair = factor_lj*forcelj*r2inv;
        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }
          
        if (EFLAG) {
        if((connect_degree == 1) || (connect_degree == 2))
          {
            evdwl = r6inv*rep2i[jtype];
          }else{
          if ((i_dom_id == 0) || (j_dom_id == 0) || (i_mol_id != j_mol_id)||(i_dom_id != j_dom_id)){
          evdwl =r6inv * (lj3i[jtype]*r6inv - lj4i[jtype]) - offseti[jtype];  // SOP lj potential 
          evdwl *= factor_lj;
          }else{ evdwl = 0.0;}
          }
        
        } //end of e flag 
        
        
         if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
        
     
        } // end of lj non-bonded calculation
      } // end of j loop
     f[i].x += fxtmp;
     f[i].y += fytmp;
     f[i].z += fztmp;     
    } // end of i loop

} // end of eval

/* ---------------------------------------------------------------------- */

double PairLJ_SOP_OMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJ_SOP::memory_usage();

  return bytes;
}
