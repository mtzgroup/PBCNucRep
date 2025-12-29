#ifndef PERIODIC_LATTICE_EXPORT_H_
#define PERIODIC_LATTICE_EXPORT_H_

#include "periodic_types.h"

#define DSOGLOBAL __attribute__((visibility("default")))

#ifdef __cplusplus
extern "C" {
#endif
  namespace PeriodicBox
  {
    namespace LatticeInfoMethods
    {
      extern DSOGLOBAL void calculate_lattice_vector(LatticeInfo* lattice);
      extern DSOGLOBAL void set_lattice_with_vector(LatticeInfo* lattice);
      extern DSOGLOBAL void set_trash_default(LatticeInfo* lattice);
      extern DSOGLOBAL void print(const LatticeInfo& lattice);
    }

    extern DSOGLOBAL void grid_density_to_n_grid_each_dimension(const double grid_density, const LatticeInfo& lattice, int n_grid_each_dimension[3]);

    extern DSOGLOBAL void get_absolute_coord_from_fractional_coord(double absolute[3], const double fractional[3], const double unit_cell[9]);
  }
#ifdef __cplusplus
}  // End extern C
#endif

#undef DSOGLOBAL

#endif
