/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_sop.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJ_SOP::PairLJ_SOP(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJ_SOP::~PairLJ_SOP()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon_lj);
    memory->destroy(epsilon_rep);    
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(rep1);
    memory->destroy(rep2);        
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJ_SOP::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj;
  double r, rshift, rshiftsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *domain_id = atom->domain_id;       // domain tag of the atom
  int *molecule_id = atom->molecule;     // molecule tag of the atom  
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  int connect_degree,i_dom_id, j_dom_id, i_mol_id, j_mol_id;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    i_dom_id  = domain_id[i];
    i_mol_id = molecule_id[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      connect_degree = sbmask(j); // identify connectivity between beads i and j
      //std::cout << domain_id[j]<<"\t" ;
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
    //  std::cout << domain_id[i]<<"\n" ;
      j_dom_id = domain_id[j];
      j_mol_id = molecule_id[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
       

      
      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;
        if((connect_degree == 1) || (connect_degree == 2)){
	  forcelj =  r6inv*rep1[itype][jtype];
	        
        } // end of excluded force calculation
         else{        
        if ((i_dom_id == 0) || (j_dom_id == 0) || (i_mol_id != j_mol_id)||(i_dom_id != j_dom_id)){
	   forcelj = factor_lj*r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
      }else{forcelj =0.0;} // end of lj non-bonded calculation
    }
     fpair = factor_lj * forcelj * r2inv;
     f[i][0] += delx * fpair;
     f[i][1] += dely * fpair;
     f[i][2] += delz * fpair;
     if (newton_pair || j < nlocal) {
       f[j][0] -= delx * fpair;
       f[j][1] -= dely * fpair;
       f[j][2] -= delz * fpair;
      }
    
    if (eflag)
    {

      if((connect_degree == 1) || (connect_degree == 2))
      {
         evdwl = r6inv*rep2[itype][jtype];
      }else
        {
          if ((i_dom_id == 0) || (j_dom_id == 0) || (i_mol_id != j_mol_id)||(i_dom_id != j_dom_id))
            {
             evdwl =r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]) - offset[itype][jtype];  // SOP lj potential 
             evdwl *= factor_lj;
            }else
              {
               evdwl = 0.0;
              }
        }

     } // end of eflag
    if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);

    } // within cutoff
    } // end of loop over j bead
  } // end of loop over i bead

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJ_SOP::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");

  memory->create(cut, np1, np1, "pair:cut");
  memory->create(epsilon_rep, np1, np1, "pair:epsilon_rep");
  memory->create(epsilon_lj, np1, np1, "pair:epsilon_lj");  
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(lj1, np1, np1, "pair:lj1");
  memory->create(lj2, np1, np1, "pair:lj2");
  memory->create(lj3, np1, np1, "pair:lj3");
  memory->create(lj4, np1, np1, "pair:lj4");
  memory->create(rep1, np1, np1, "pair:rep1");  
  memory->create(rep2, np1, np1, "pair:rep2");    
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJ_SOP::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJ_SOP::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_rep_one = utils::numeric(FLERR, arg[2], false, lmp);  
  double epsilon_lj_one = utils::numeric(FLERR, arg[3], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[4], false, lmp);

  double cut_one = cut_global;
  if (narg == 6) cut_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon_rep[i][j] = epsilon_rep_one;
      epsilon_lj[i][j] = epsilon_lj_one;      
      sigma[i][j] = sigma_one;

      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJ_SOP::init_one(int i, int j)
{


  if (setflag[i][j] == 0)  error->all(FLERR, "All pair coeffs are not set");

  lj1[i][j] = 12.0 * epsilon_lj[i][j] * pow(sigma[i][j],12.0); // LJ SOP force related
  lj2[i][j] = 12.0 * epsilon_lj[i][j] * pow(sigma[i][j],6.0); // LJ SOP force related
  lj3[i][j] =  epsilon_lj[i][j] * pow(sigma[i][j],12.0); // LJ SOP energy related
  lj4[i][j] = 2.0 * epsilon_lj[i][j] * pow(sigma[i][j],6.0); // LJ SOP energy related
 // std::cout << "lj_rep: " << lj2[i][j] <<"\n";  
  rep1[i][j] = 6.0*epsilon_rep[i][j] * pow(sigma[i][j],6.0); // excluded volume force related
  rep2[i][j] = epsilon_rep[i][j] * pow(sigma[i][j],6.0); // excluded volume energy related  


  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / (cut[i][j] );
    offset[i][j] = epsilon_lj[i][j] * (pow(ratio, 12.0) - 2* pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;


  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  rep1[j][i] = rep1[i][j];
  rep2[j][i] = rep2[i][j];  
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce
//std::cout <<"tail_flag : "<<tail_flag<<"\n";
  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;
//    double prefactor = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 / (9.0 * rc9);
//    etail_ij = prefactor * (sig6 - 3.0 * rc6);

    double prefactor = 2.0 * MY_PI * all[0] * all[1] * epsilon_lj[i][j] * sig6 / (9.0 * rc9);
    etail_ij = prefactor * (sig6 - 6.0 * rc6);    
//	std::cout <<"etail : "<<etail_ij << "\n";    
    ptail_ij = 2.0 * prefactor * (2.0 * sig6 - 3.0 * rc6);
  }

//  return cut[i][j] + shift[i][j];
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_rep[i][j], sizeof(double), 1, fp);
        fwrite(&epsilon_lj[i][j], sizeof(double), 1, fp);        
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ_SOP::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon_rep[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &epsilon_lj[i][j], sizeof(double), 1, fp, nullptr, error);          
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon_rep[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&epsilon_lj[i][j], 1, MPI_DOUBLE, 0, world);        
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ_SOP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, epsilon_rep[i][i], epsilon_lj[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, epsilon_rep[i][j], epsilon_lj[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */



void *PairLJ_SOP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon_rep") == 0) return (void *) epsilon_rep;
  if (strcmp(str, "epsilon_lj") == 0) return (void *) epsilon_lj;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;

  return nullptr;
}
