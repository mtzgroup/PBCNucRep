#include "periodic_lattice_export.h"
#include "periodic_lattice_internal.h"

#include "helper.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define FABS_MAX(X, Y) ((fabs(X) > fabs(Y)) ? fabs(X) : fabs(Y))
#define NORM2(X) (((X)[0]) * ((X)[0]) + ((X)[1]) * ((X)[1]) + ((X)[2]) * ((X)[2]))

namespace PeriodicBox
{
  static double get_min_nonzero_lattice_vector_norm(const double unit_cell[9])
  {
    double P[9];
    memcpy(P, unit_cell, 9 * sizeof(double));
    const double P2_on_P1 = (P[0] * P[3] + P[1] * P[4] + P[2] * P[5]) / NORM2(P + 0);
    const int P2_on_P1_n_copy = (int) round(P2_on_P1);
    for (int i_xyz = 0; i_xyz < 3; i_xyz++)
      P[3 + i_xyz] -= P2_on_P1_n_copy * P[0 + i_xyz];

    const double P3_on_P1 = (P[0] * P[6] + P[1] * P[7] + P[2] * P[8]) / NORM2(P + 0);
    const int P3_on_P1_n_copy = (int) round(P3_on_P1);
    const double P3_on_P2 = (P[3] * P[6] + P[4] * P[7] + P[5] * P[8]) / NORM2(P + 3);
    const int P3_on_P2_n_copy = (int) round(P3_on_P2);
    for (int i_xyz = 0; i_xyz < 3; i_xyz++)
    {
      P[6 + i_xyz] -= P3_on_P1_n_copy * P[0 + i_xyz];
      P[6 + i_xyz] -= P3_on_P2_n_copy * P[3 + i_xyz];
    }

    const double P1_norm2 = NORM2(P + 0), P2_norm2 = NORM2(P + 3), P3_norm2 = NORM2(P + 6);
    const double min_norm = sqrt(MIN(MIN(P1_norm2, P2_norm2), P3_norm2));

    if (min_norm < 1e-6)
      DIE("Lattice vectors are integer linear combination of each other!\n");

    return min_norm;
  }

  // Some explanation:
  // unit_cell.real_space contains R matrix, each row of R matrix is a lattice vector
  // unit_cell.reciprocal_space contains G matrix, G * R^T = 2pi I, (R^T)^-1 = 1/2pi G
  // The lattice index can be computed as n = (R * R^T)^-1 (R * v) = 1/2pi G * v
  static double get_min_norm_lattice_image_real(const LatticeVector unit_cell, const double v[3], double v_min[3])
  {
    const double* G = unit_cell.reciprocal_space;
    const double lattice_index_with_fraction[3] {
      0.5 / PI * (G[0] * v[0] + G[1] * v[1] + G[2] * v[2]),
      0.5 / PI * (G[3] * v[0] + G[4] * v[1] + G[5] * v[2]),
      0.5 / PI * (G[6] * v[0] + G[7] * v[1] + G[8] * v[2]),
    };
    const double lattice_index[3] {
      round(lattice_index_with_fraction[0]),
      round(lattice_index_with_fraction[1]),
      round(lattice_index_with_fraction[2]),
    };

    const double* R = unit_cell.real_space;
    v_min[0] = v[0] - (R[0] * lattice_index[0] + R[3] * lattice_index[1] + R[6] * lattice_index[2]);
    v_min[1] = v[1] - (R[1] * lattice_index[0] + R[4] * lattice_index[1] + R[7] * lattice_index[2]);
    v_min[2] = v[2] - (R[2] * lattice_index[0] + R[5] * lattice_index[1] + R[8] * lattice_index[2]);

    return sqrt(NORM2(v_min));
  }

  void LatticeInfoMethods::calculate_lattice_vector(LatticeInfo* lattice)
  {
    if (lattice->dimension <= 0 || lattice->dimension > 3 ||
        lattice->a <= 0.0 || lattice->b <= 0.0 || lattice->c <= 0.0 ||
        lattice->alpha <= 0.0 || lattice->beta <= 0.0 || lattice->gamma <= 0.0)
      DIE("PeriodicBox::LatticeInfoMethods::CalculateLatticeVector() Lattice parameters not properly set!");

    double* R = lattice->unit_cell.real_space;
    double* K = lattice->unit_cell.reciprocal_space;

    const double a = lattice->a, b = lattice->b, c = lattice->c;
    const double alpha = lattice->alpha, beta = lattice->beta, gamma = lattice->gamma;
    // https://en.wikipedia.org/wiki/Fractional_coordinates#In_crystallography
    R[0] = a;
    R[1] = 0.0;
    R[2] = 0.0;
    R[3] = b * cos(gamma);
    R[4] = b * sin(gamma);
    R[5] = 0.0;
    R[6] = c * cos(beta);
    R[7] = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma);
    lattice->V_real = a * b * c
                    * sqrt(1 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma));
    R[8] = lattice->V_real / (a * b * sin(gamma));

    K[0] = 2 * PI / lattice->V_real * (R[4] * R[8] - R[5] * R[7]);
    K[1] = 2 * PI / lattice->V_real * (R[5] * R[6] - R[3] * R[8]);
    K[2] = 2 * PI / lattice->V_real * (R[3] * R[7] - R[4] * R[6]);
    K[3] = 2 * PI / lattice->V_real * (R[7] * R[2] - R[8] * R[1]);
    K[4] = 2 * PI / lattice->V_real * (R[8] * R[0] - R[6] * R[2]);
    K[5] = 2 * PI / lattice->V_real * (R[6] * R[1] - R[7] * R[0]);
    K[6] = 2 * PI / lattice->V_real * (R[1] * R[5] - R[2] * R[4]);
    K[7] = 2 * PI / lattice->V_real * (R[2] * R[3] - R[0] * R[5]);
    K[8] = 2 * PI / lattice->V_real * (R[0] * R[4] - R[1] * R[3]);

    lattice->K_min_norm = get_min_nonzero_lattice_vector_norm(K);
  }

  void LatticeInfoMethods::set_lattice_with_vector(LatticeInfo* lattice)
  {
    if (lattice->dimension <= 0 || lattice->dimension > 3)
      DIE("PeriodicBox::LatticeInfoMethods::SetLatticeWithVector() Lattice parameters not properly set!");

    const double* R = lattice->unit_cell.real_space;
    double* K = lattice->unit_cell.reciprocal_space;

    lattice->V_real = R[0] * R[4] * R[8] + R[1] * R[5] * R[6] + R[2] * R[3] * R[7]
                    - R[2] * R[4] * R[6] - R[0] * R[5] * R[7] - R[1] * R[3] * R[8];

    K[0] = 2 * PI / lattice->V_real * (R[4] * R[8] - R[5] * R[7]);
    K[1] = 2 * PI / lattice->V_real * (R[5] * R[6] - R[3] * R[8]);
    K[2] = 2 * PI / lattice->V_real * (R[3] * R[7] - R[4] * R[6]);
    K[3] = 2 * PI / lattice->V_real * (R[7] * R[2] - R[8] * R[1]);
    K[4] = 2 * PI / lattice->V_real * (R[8] * R[0] - R[6] * R[2]);
    K[5] = 2 * PI / lattice->V_real * (R[6] * R[1] - R[7] * R[0]);
    K[6] = 2 * PI / lattice->V_real * (R[1] * R[5] - R[2] * R[4]);
    K[7] = 2 * PI / lattice->V_real * (R[2] * R[3] - R[0] * R[5]);
    K[8] = 2 * PI / lattice->V_real * (R[0] * R[4] - R[1] * R[3]);

    lattice->K_min_norm = get_min_nonzero_lattice_vector_norm(K);

    lattice->a = sqrt(NORM2(R + 0));
    lattice->b = sqrt(NORM2(R + 3));
    lattice->c = sqrt(NORM2(R + 6));
    lattice->alpha = acos((R[3] * R[6] + R[4] * R[7] + R[5] * R[8]) / (lattice->b * lattice->c));
    lattice->beta = acos((R[0] * R[6] + R[1] * R[7] + R[2] * R[8]) / (lattice->a * lattice->c));
    lattice->gamma = acos((R[0] * R[3] + R[1] * R[4] + R[2] * R[5]) / (lattice->a * lattice->b));
  }

  void LatticeInfoMethods::set_trash_default(LatticeInfo* lattice)
  {
    lattice->dimension = -1;
    lattice->V_real = NAN;
    lattice->K_min_norm = NAN;
    lattice->a = 0.0;
    lattice->b = 0.0;
    lattice->c = 0.0;
    lattice->alpha = 0.0;
    lattice->beta  = 0.0;
    lattice->gamma = 0.0;
    for (int i = 0; i < 9; i++)
    {
      lattice->unit_cell.real_space[i] = NAN;
      lattice->unit_cell.reciprocal_space[i] = NAN;
    }
  }

  void LatticeInfoMethods::print(const LatticeInfo& lattice)
  {
    printf("%dD periodic system\n", lattice.dimension);
    printf("a = %.10f A, b = %.10f A, c = %.10f A\n",
      lattice.a * BohrToAng, lattice.b * BohrToAng, lattice.c * BohrToAng);
    printf("alpha = %.10f deg, beta = %.10f deg, gamma = %.10f deg\n",
      lattice.alpha / PI * 180.0, lattice.beta / PI * 180.0, lattice.gamma / PI * 180.0);

    printf("Real space lattice vector in a.u.:\n");
    for (size_t i_dimension = 0; i_dimension < 3; i_dimension++)
      printf("(%.10f, %.10f, %.10f), ",
        lattice.unit_cell.real_space[i_dimension * 3 + 0],
        lattice.unit_cell.real_space[i_dimension * 3 + 1],
        lattice.unit_cell.real_space[i_dimension * 3 + 2]);
    printf("\n"); fflush(stdout);

    printf("Reciprocal space lattice vector in a.u.:\n");
    for (size_t i_dimension = 0; i_dimension < 3; i_dimension++)
      printf("(%.10f, %.10f, %.10f), ",
        lattice.unit_cell.reciprocal_space[i_dimension * 3 + 0],
        lattice.unit_cell.reciprocal_space[i_dimension * 3 + 1],
        lattice.unit_cell.reciprocal_space[i_dimension * 3 + 2]);

    printf("\n\n"); fflush(stdout);
  }

  void grid_density_to_n_grid_each_dimension(const double grid_density, const LatticeInfo& lattice, int n_grid_each_dimension[3])
  {
    // The following logic assumes the lattice vector matrix is triangular
    const double lenght_xyz[3] { lattice.unit_cell.real_space[0], lattice.unit_cell.real_space[4], lattice.unit_cell.real_space[8], };
    const double n_grid_total = grid_density * lattice.V_real;
    n_grid_each_dimension[0] = (int)ceil( pow(n_grid_total * lenght_xyz[0] * lenght_xyz[0] / lenght_xyz[1] / lenght_xyz[2], 1.0 / 3.0) );
    n_grid_each_dimension[1] = (int)ceil( pow(n_grid_total * lenght_xyz[1] * lenght_xyz[1] / lenght_xyz[0] / lenght_xyz[2], 1.0 / 3.0) );
    n_grid_each_dimension[2] = (int)ceil( pow(n_grid_total * lenght_xyz[2] * lenght_xyz[2] / lenght_xyz[0] / lenght_xyz[1], 1.0 / 3.0) );
  }

  void get_absolute_coord_from_fractional_coord(double absolute[3], const double fractional[3], const double unit_cell[9])
  {
    dgemv('N', 3, 3, 1.0, unit_cell, 3, fractional, 1, 0.0, absolute, 1);
  }

  static void get_cube_bound(int positive_bound[3], int negative_bound[3],
                            const LatticeVector unit_cell, const double origin_offset_absolute[3], const double radius,
                            const bool if_real)
  {
    const double positive_direction_standard_basis[9] { radius - origin_offset_absolute[0], 0.0, 0.0,
                                                        0.0, radius - origin_offset_absolute[1], 0.0,
                                                        0.0, 0.0, radius - origin_offset_absolute[2] };
    const double negative_direction_standard_basis[9] { -radius - origin_offset_absolute[0], 0.0, 0.0,
                                                        0.0, -radius - origin_offset_absolute[1], 0.0,
                                                        0.0, 0.0, -radius - origin_offset_absolute[2] };
    double P_inverse_T[9];
    memcpy(P_inverse_T, if_real ? unit_cell.reciprocal_space : unit_cell.real_space, 9 * sizeof(double));
    for (int i = 0; i < 9; i++)
      P_inverse_T[i] /= 2.0 * PI;

    double positive_direction_lattice_basis[9];
    double negative_direction_lattice_basis[9];
    dgemm('T', 'N', 3, 3, 3, 1.0, P_inverse_T, 3, positive_direction_standard_basis, 3, 0.0, positive_direction_lattice_basis, 3);
    dgemm('T', 'N', 3, 3, 3, 1.0, P_inverse_T, 3, negative_direction_standard_basis, 3, 0.0, negative_direction_lattice_basis, 3);

    for (int i = 0; i < 3; i++)
    {
      positive_bound[i] = (int) ceil(FABS_MAX(FABS_MAX(positive_direction_lattice_basis[i + 0],
                                                       positive_direction_lattice_basis[i + 3]),
                                                       positive_direction_lattice_basis[i + 6]));
      negative_bound[i] = (int) ceil(FABS_MAX(FABS_MAX(negative_direction_lattice_basis[i + 0],
                                                       negative_direction_lattice_basis[i + 3]),
                                                       negative_direction_lattice_basis[i + 6]));
      // Notice the algorithm above will include unnecessary ceils
      // if origin_offset_absolute is outside of the sphere centered at zero with radius radius.
    }
  }

  std::vector<double3> compute_lattice_vectors_within_sphere_real(const LatticeVector unit_cell,
                                                                  const double origin_offset_absolute[3],
                                                                  const double radius,
                                                                  const bool remove_origin)
  {
    int positive_direction_lattice_bound[3] { 0, 0, 0 };
    int negative_direction_lattice_bound[3] { 0, 0, 0 };

    double min_norm_origin[3];
    get_min_norm_lattice_image_real(unit_cell, origin_offset_absolute, min_norm_origin);

    get_cube_bound(positive_direction_lattice_bound, negative_direction_lattice_bound,
                   unit_cell, min_norm_origin, radius, true);

    const double r2 = radius * radius;

    std::vector<double3> lattice_vectors;
    for (int i_x = -negative_direction_lattice_bound[0]; i_x <= positive_direction_lattice_bound[0]; i_x++)
      for (int i_y = -negative_direction_lattice_bound[1]; i_y <= positive_direction_lattice_bound[1]; i_y++)
        for (int i_z = -negative_direction_lattice_bound[2]; i_z <= positive_direction_lattice_bound[2]; i_z++)
        {
          if (remove_origin && i_x == 0 && i_y == 0 && i_z == 0)
            continue;

          const double fractional_v[3] { (double)i_x, (double)i_y, (double)i_z };
          double absolute_v[3];
          get_absolute_coord_from_fractional_coord(absolute_v, fractional_v, unit_cell.real_space);
          double absolute_v_plus_C[3];
          for (int i = 0; i < 3; i++)
            absolute_v_plus_C[i] = absolute_v[i] + min_norm_origin[i];
          for (int i = 0; i < 3; i++)
            absolute_v[i] -= origin_offset_absolute[i] - min_norm_origin[i];
          if (NORM2(absolute_v_plus_C) < r2)
          {
            const double3 lattice_vector { absolute_v[0], absolute_v[1], absolute_v[2] };
            lattice_vectors.push_back(lattice_vector);
          }
        }

    return lattice_vectors;
  }

  std::vector<double3> compute_lattice_vectors_within_sphere_reciprocal(const LatticeVector unit_cell,
                                                                        const double radius)
  {
    int positive_direction_lattice_bound[3] { 0, 0, 0 };
    int negative_direction_lattice_bound[3] { 0, 0, 0 };

    const double origin[3] { 0.0, 0.0, 0.0 };
    get_cube_bound(positive_direction_lattice_bound, negative_direction_lattice_bound,
                   unit_cell, origin, radius, false);

    const double r2 = radius * radius;

    std::vector<double3> lattice_vectors;
    for (int i_x = -negative_direction_lattice_bound[0]; i_x <= positive_direction_lattice_bound[0]; i_x++)
      for (int i_y = -negative_direction_lattice_bound[1]; i_y <= positive_direction_lattice_bound[1]; i_y++)
        for (int i_z = -negative_direction_lattice_bound[2]; i_z <= positive_direction_lattice_bound[2]; i_z++)
        {
          if (i_x == 0 && i_y == 0 && i_z == 0)
            continue;

          const double fractional_v[3] { (double)i_x, (double)i_y, (double)i_z };
          double absolute_v[3];
          get_absolute_coord_from_fractional_coord(absolute_v, fractional_v, unit_cell.reciprocal_space);
          if (NORM2(absolute_v) < r2)
          {
            const double3 lattice_vector { absolute_v[0], absolute_v[1], absolute_v[2] };
            lattice_vectors.push_back(lattice_vector);
          }
        }

    return lattice_vectors;
  }

  double3 get_image_in_origin_cell(const LatticeVector unit_cell, const double3 point)
  {
    const double* G = unit_cell.reciprocal_space;
    const double lattice_index_with_fraction[3] {
      0.5 / PI * (G[0] * point.x + G[1] * point.y + G[2] * point.z),
      0.5 / PI * (G[3] * point.x + G[4] * point.y + G[5] * point.z),
      0.5 / PI * (G[6] * point.x + G[7] * point.y + G[8] * point.z),
    };
    const double lattice_index[3] {
      floor(lattice_index_with_fraction[0]),
      floor(lattice_index_with_fraction[1]),
      floor(lattice_index_with_fraction[2]),
    };

    const double* R = unit_cell.real_space;
    const double3 shifted_point {
      point.x - (R[0] * lattice_index[0] + R[3] * lattice_index[1] + R[6] * lattice_index[2]),
      point.y - (R[1] * lattice_index[0] + R[4] * lattice_index[1] + R[7] * lattice_index[2]),
      point.z - (R[2] * lattice_index[0] + R[5] * lattice_index[1] + R[8] * lattice_index[2]),
    };
    return shifted_point;
  }
}
