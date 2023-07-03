#include <iostream>
#include "Input.h"
#include <vector>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "omp.h"
#include <fstream>

#ifdef __APPLE__

#include <vecLib/clapack.h>

#else
#include "lapacke.h"
#endif

int main(int argc, char **argv)
{

    Input input;
    input.read_input("INPUT.txt");
    std::vector<std::vector<double>> pt = input.get_points();
    std::vector<double> dis = input.get_distributions();
    boost::math::cubic_b_spline<double> spline(dis.begin(), dis.end(), 0, input.dr);

    double *v;
    input.read_v(&v);

    double x = input.lx / (double)input.nx;
    double y = input.ly / (double)input.ny;
    double z = input.lz / (double)input.nz;
    double h[pt.size() * pt.size()];

    for (int m = 0; m < pt.size(); ++m)
    {
        for (int n = 0; n < pt.size(); ++n)
        {
            double sum = 0;
            std::vector<double> pt1 = pt[m], pt2 = pt[n];
#pragma omp parallel for reduction(+ : sum)
            for (int i = 0; i < input.nx; ++i)
            {
                for (int j = 0; j < input.ny; ++j)
                {
                    for (int k = 0; k < input.nz; ++k)
                    {
                        double dis1 = sqrt((i * x - pt1[0]) * (i * x - pt1[0]) + (j * y - pt1[1]) * (j * y - pt1[1]) +
                                           (k * z - pt1[2]) * (k * z - pt1[2]));
                        double dis2 = sqrt((i * x - pt2[0]) * (i * x - pt2[0]) + (j * y - pt2[1]) * (j * y - pt2[1]) +
                                           (k * z - pt2[2]) * (k * z - pt2[2]));
                        if (dis1 > input.cutoff or dis2 > input.cutoff)
                        {
                            continue;
                        }
                        else
                        {
                            sum += v[i * input.ny * input.nz + j * input.nz + k] * spline(dis1) * spline(dis2) * x * y * z;
                        }
                    }
                }
            }
            h[m * pt.size() + n] = sum;
        }
    }
    free(v);

    int info, lwork;
    double wkopt;
    double w[pt.size()];
    int r = pt.size();
#ifdef __APPLE__
    lwork = -1;
    dsyev_("V", "U", &r, h, &r, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    double work[(int)wkopt];
    /* Solve eigenproblem */
    dsyev_("V", "U", &r, h, &r, w, work, &lwork, &info);
    /* Check for convergence */
    if (info > 0)
    {
        std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
        exit(-1);
    }

#else
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', r, h, r, w);
    if (info > 0)
    {
        std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
        exit(-1);
    }
#endif
    std::ofstream ofs1{"eigenvalues.log"}, ofs2{"eigenvectors.log"};

    for (auto i : w)
    {
        ofs2 << std::setprecision(6) << std::fixed << i << std::endl;
    }
    for (int i = 1; i <= pt.size() * pt.size(); ++i)
    {
        ofs1 << h[i] << ' ';
        if (i % pt.size() == 0)
        {
            ofs1 << std::endl;
        }
    }

    return 0;
}
