#ifndef FINAL_PROJ_INPUT_H
#define FINAL_PROJ_INPUT_H

#include <string>
#include <vector>

class Input
{
public:
    Input();

    ~Input();

    void read_input(const std::string &fn);

    std::vector<std::vector<double>> get_points();

    std::vector<double> get_distributions();

    void read_v(double **p);

    int isHexahedral = 0, lx = 1000, ly = 1000, lz = 1000, p_count = 0;
    int thetaxy = 0, thetayz = 0, thetaxz = 0;
    int support_SH = 0, support_Periodic_Boundary = 0, multi_parallel_strategies = 0;
    std::string diago_lib{"lapack"};
    std::string points_path, v_path, distribution_path;
    int nx, ny, nz;
    double cutoff, dr;
    int mesh;
};

#endif // FINAL_PROJ_INPUT_H
