/* ----------------------------------------------------------------------
 *   This is the MLIP-LAMMPS interface
 *   MLIP is a software for Machine Learning Interatomic Potentials
 *   distributed by A. Shapeev, Skoltech (Moscow)
 *   Contributors: Evgeny Podryabinkin

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   LAMMPS is distributed under a GNU General Public License
   and is not a part of MLIP.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Evgeny Podryabinkin (Skoltech)
------------------------------------------------------------------------- */

/* This is a modified version of the interface wrapper of MLIP - MTP based on the 
   correction about the many body potential heat flux definition stated in Zheyong Fan et al.
   (2015), " Force and heat currentt formulas for many-body potentials in molecular dynamics 
   simulation with applications to thermal conductivity calculations", Phys. Rev. B 92, 094301.
   The effect of this version to the calculation result is discussed in the article "Revisit 
   Many-body Interaction Heat Current and Thermal Conductivity Calculation in Moment Tensor 
   Potential/LAMMPS Interface",2024,arXiv:2411.01255
   Modified by Siu Ting Tai, HKU MECH
*/

#ifdef PAIR_CLASS

PairStyle(mlip,PairMLIP)

#else

#ifndef LMP_PAIR_MLIP_H
#define LMP_PAIR_MLIP_H

#include <stdio.h>
#include "pair.h"

extern void MLIP_init(const char*, const char*, int, double&, int&);
extern void MLIP_calc_cfg(int, double*, double**, int*, int*, double&, double**, double*);
//add extra pointer for nbh evaluation for total virial
extern void MLIP_calc_nbh(int, int*, int*, int**, int, int, double**, int*, double**, double&, double*, double**, double*);
extern void MLIP_finalize();

namespace LAMMPS_NS {

class PairMLIP : public Pair {
 public:
  double cutoff;

  PairMLIP(class LAMMPS *);
  virtual ~PairMLIP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  int mode; // 0 - nbh mode (can't learn on the fly), 1 - cfg mode (typically for non-parallel lammps)
  bool inited;
  char MLIPsettings_filename[1000];
  char MLIPlog_filename[1000];
  double cutoffsq;
  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
