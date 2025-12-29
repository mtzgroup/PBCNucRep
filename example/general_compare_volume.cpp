#include "../src/periodic_lattice_export.h"
#include "../src/helper.h"
#include "../src/periodic_grid.h"
#include "../src/periodic_becke/periodic_becke.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <deque>
#include <vector>

int main()
{
    const char* filename = "compare_volume.txt";

    std::ifstream infile(filename);
    if (!infile) {
        printf("Cannot open file %s\n", filename);
        return 1;
    }

    std::string line;
    std::deque<std::vector<double>> data;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);

        if (line.empty())
            continue;

        std::vector<double> data_one_line;
        double value;

        while (iss >> value) {
            data_one_line.push_back(value);
        }

        data.push_back(data_one_line);
    }

    const std::vector<double> line1 = data.front();
    data.pop_front();
    if (line1.size() != 9) {
        printf("The first line should contain 9 of the lattice vector values");
        return 1;
    }

    PeriodicBox::LatticeInfo lattice;
    lattice.dimension = 3;
    for (int i = 0; i < 9; i++)
        lattice.unit_cell.real_space[i] = line1[i];
    PeriodicBox::LatticeInfoMethods::set_lattice_with_vector(&lattice);
    PeriodicBox::LatticeInfoMethods::print(lattice);

    PeriodicBox::PeriodicParameter periodic_parameter;
    periodic_parameter.min_primitive_exponent = NAN;
    periodic_parameter.max_primitive_exponent = NAN;
    periodic_parameter.lattice = lattice;
    periodic_parameter.thresholds.periodic_orbtial_cutoff = NAN;
    periodic_parameter.thresholds.periodic_charge_cutoff_real = 1e-14;
    periodic_parameter.thresholds.periodic_charge_cutoff_reciprocal = 1e-14;
    periodic_parameter.omega = NAN;

    std::vector<double> atom_xyz;
    std::vector<double> atom_radius;

    while (!data.empty()) {
        const std::vector<double> line_xyz = data.front();
        if (line_xyz.size() != 4)
            break;

        data.pop_front();
        atom_radius.push_back(line_xyz[0]);
        atom_xyz.push_back(line_xyz[1]);
        atom_xyz.push_back(line_xyz[2]);
        atom_xyz.push_back(line_xyz[3]);
    }
    const int n_atom = atom_radius.size();

    std::vector<PeriodicBox::GridPoint> grid_points;

    while (!data.empty()) {
        const std::vector<double> line_grid = data.front();
        if (line_grid.size() != 5)
            break;

        data.pop_front();
        PeriodicBox::GridPoint grid_point;
        grid_point.x = line_grid[0];
        grid_point.y = line_grid[1];
        grid_point.z = line_grid[2];
        grid_point.w_fixed = line_grid[3];
        grid_point.w_total = NAN;
        grid_point.i_atom = (int) round(line_grid[4]);

        grid_points.push_back(grid_point);
    }
    const int n_grid_point = grid_points.size();

    printf("Finish loading, n_atom = %d, n_grid_point = %d\n", n_atom, n_grid_point);
    fflush(stdout);

    BeckeWeights(n_grid_point, grid_points.data(), n_atom, atom_xyz.data(), atom_radius.data(), periodic_parameter.lattice.unit_cell);

    double weight_sum = 0.0;
    for (int i_grid = 0; i_grid < n_grid_point; i_grid++)
        weight_sum += grid_points[i_grid].w_total;

    printf("Weight sum = %.5f, cell volume = %.5f, difference = %.5f\n", weight_sum, lattice.V_real, weight_sum - lattice.V_real);

    return 0;
}
