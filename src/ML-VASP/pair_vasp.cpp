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

/* ----------------------------------------------------------------------
   Contributing author: VASP Software GmbH (www.vasp.at)
------------------------------------------------------------------------- */

#include "pair_vasp.h"
#include "InterfaceLAMMPS.hpp"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace std;

using VaspError = vaspml::InterfaceLAMMPS::Exception;

/* ----------------------------------------------------------------------
   report an error of VASPml through LAMMPS, all_ranks: every rank of
   world is known to raise the same error (else only this rank does)
------------------------------------------------------------------------- */

[[noreturn]] static void report_vasp_error(Error *error, const VaspError &e, bool all_ranks)
{
  if (e.isInternal())
    error->one(FLERR, "Internal error in VASPml, please report it: " + std::string(e.what()));
  if (all_ranks) error->all(FLERR, std::string(e.what()));
  error->one(FLERR, std::string(e.what()));
}

static const char cite_user_vaspml_package[] =
    "VASP package: 10.1103/PhysRevB.100.014105\n\n"
    "@Article{JinnouchiKarsaiEtAl2019PRB,\n"
    "author = {Jinnouchi, Ryosuke and Karsai, Ferenc and Kresse, Georg},\n"
    "title = {On-the-Fly Machine Learning Force Field Generation: {{Application}} to Melting "
    "Points},\n"
    "year = {2019},\n"
    "month = jul,\n"
    "journal = {Physical Review B},\n"
    "volume = {100},\n"
    "number = {1},\n"
    "pages = {014105},\n"
    "publisher = {{American Physical Society}},\n"
    "doi = {10.1103/PhysRevB.100.014105}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairVASP::PairVASP(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_vaspml_package);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  //unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  vasp = new vaspml::InterfaceLAMMPS();

  // route alerts, warnings and advice of VASPml to LAMMPS warnings and its
  // regular output to screen and log file, both on rank 0 only

  vasp->setMessageHandlers(
      [this](const std::string &text) {
        if (comm->me == 0) error->warning(FLERR, text);
      },
      [this](const std::string &text) {
        if (comm->me == 0) utils::logmesg(this->lmp, text + "\n");
      });
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairVASP::~PairVASP()
{
  delete vasp;
  memory->destroy(setflag);
  memory->destroy(cutsq);
}

/* ---------------------------------------------------------------------- */

void PairVASP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  double localEnergy = 0.0;
  try {
    vasp->setTags(atom->tag);
    vasp->computeNeighborLists(atom->x,
                               atom->type,
                               list->inum,
                               list->ilist,
                               list->numneigh,
                               list->firstneigh,
                               NEIGHMASK);
    vasp->predict();
    localEnergy = vasp->getLocalEnergy();
    vasp->getForces(atom->f);
  } catch (VaspError &e) {
    report_vasp_error(error, e, false);
  }
  if (eflag_global) ev_tally(0,0,atom->nlocal,1,localEnergy,0.0,0.0,0.0,0.0,0.0);
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairVASP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairVASP::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR, "Pair style vasp takes no arguments, found {}", arg[0]);
  if (strcmp("metal", update->unit_style) != 0)
    error->all(FLERR, "VASP potentials require 'metal' units");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairVASP::coeff(int narg, char **arg)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "pair_coeff", error);
  if (!allocated) allocate();

  if (comm->me == 0) utils::logmesg(lmp, "pair_coeff vasp:\n");
  auto mlff_file_name = utils::get_potential_file_path(arg[2]);
  if (mlff_file_name.empty()) error->all(FLERR, "Cannot open VASP force field file {}", arg[2]);
  if (comm->me == 0)
  {
    utils::logmesg(lmp, "  VASP ML force field file   = {}\n", mlff_file_name);
    utils::logmesg(lmp, "  Setting up force field ...\n");
  }
  try {
    cutoff = vasp->setupForceField(mlff_file_name, world);
  } catch (VaspError &e) {
    report_vasp_error(error, e, true);
  }
  if (comm->me == 0)
    utils::logmesg(lmp, "  VASP ML force field cutoff = {:.2f} Angstrom\n", cutoff);
  map_element2type(narg - 3, arg + 3);
  std::vector<std::string> typesToNames;
  const int n = atom->ntypes;
  typesToNames.push_back(""); // Type 0 is not used in LAMMPS.
  for (int i = 1; i <= n; ++i) {
    if (map[i] > -1) typesToNames.push_back(elements[map[i]]);
    else typesToNames.push_back("");
  }
  std::string typeReport;
  try {
    typeReport = vasp->setupTypes(typesToNames);
  } catch (VaspError &e) {
    report_vasp_error(error, e, true);
  }
  if (comm->me == 0)
  {
    utils::logmesg(lmp, "  Please review the following mapping between LAMMPS and VASP types:\n\n");
    utils::logmesg(lmp, typeReport);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairVASP::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR,"Pair style vasp requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR,"Pair style vasp requires newton pair on");

  // need a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairVASP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  return cutoff;
}
