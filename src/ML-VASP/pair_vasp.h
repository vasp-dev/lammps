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
   Contributing author: VASP Software GmbH (www.vasp.at)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(vasp,PairVASP);
// clang-format on
#else

#ifndef LMP_PAIR_VASP_H
#define LMP_PAIR_VASP_H

#include "pair.h"

namespace vaspml {
  class InterfaceLAMMPS;
}

namespace LAMMPS_NS {

class PairVASP : public Pair {

 public:
  PairVASP(class LAMMPS *);
  ~PairVASP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  virtual void allocate();

  double cutoff;
  vaspml::InterfaceLAMMPS *vasp;

};

} // namespace LAMMPS_NS

#endif
#endif
