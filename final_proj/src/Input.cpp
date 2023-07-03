#include <fstream>
#include <iostream>
#include <sstream>
#include "Input.h"

Input::Input()
{
}

Input::~Input()
{
}

void Input::read_input(const std::string &fn)
{

    std::ifstream ifs;
    if (fn.empty())
        ifs.open("INPUT.txt");
    else
        ifs.open(fn);
    std::cout << "read from : " + fn << std::endl;
    if (!ifs.good())
    {
        std::cout << "Cannot Find Input File!" << std::endl;
        exit(0);
    }
    std::string line;

    while (getline(ifs, line))
    {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        std::string arg1;

        iss >> arg1;
        if (arg1 == std::string("isHexahedral"))
        {
            iss >> isHexahedral;
        }
        else if (arg1 == std::string("lx"))
        {
            iss >> lx;
        }
        else if (arg1 == std::string("ly"))
        {
            iss >> ly;
        }
        else if (arg1 == std::string("lz"))
        {
            iss >> lz;
        }
        else if (arg1 == std::string("thetaxy"))
        {
            iss >> thetaxy;
        }
        else if (arg1 == std::string("thetayz"))
        {
            iss >> thetayz;
        }
        else if (arg1 == std::string("thetaxz"))
        {
            iss >> thetaxz;
        }
        else if (arg1 == std::string("support_SH"))
        {
            iss >> support_SH;
        }
        else if (arg1 == std::string("diago_lib"))
        {
            iss >> diago_lib;
        }
        else if (arg1 == std::string("support_Periodic_Boundary"))
        {
            iss >> support_Periodic_Boundary;
        }
        else if (arg1 == std::string("multi_parallel_strategies"))
        {
            iss >> multi_parallel_strategies;
        }
        else if (arg1 == std::string("points_path"))
        {
            iss >> points_path;
        }
        else if (arg1 == std::string("v_path"))
        {
            iss >> v_path;
        }
        else if (arg1 == std::string("distribution_path"))
        {
            iss >> distribution_path;
        }
    }
    ifs.close();
}

std::vector<std::vector<double>> Input::get_points()
{
    std::ifstream ifs;
    std::vector<std::vector<double>> vec{};
    if (points_path.empty())
    {
        std::cout << "Cannot Find points!" << std::endl;
        exit(-1);
    }
    ifs.open(points_path);

    std::string line;

    while (getline(ifs, line))
    {
        if (line.empty())
            continue;
        double a, b, c;
        char d;
        std::istringstream iss(line);
        std::vector<double> vec1{};
        iss >> d >> a >> d >> b >> d >> c >> d;
        vec1.push_back(a);
        vec1.push_back(b);
        vec1.push_back(c);
        vec.push_back(vec1);
    }
    p_count = (int)vec.size(); // 点的数量不太多，可以类型转换不溢出
    ifs.close();
    return vec;
}

std::vector<double> Input::get_distributions()
{
    std::ifstream ifs;
    std::vector<double> vec{};

    if (points_path.empty())
    {
        std::cout << "Cannot Find distributions!";
        exit(-1);
    }
    ifs.open(distribution_path);

    std::string line;
    while (getline(ifs, line))
    {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        std::string arg1;
        iss >> arg1;
        if (arg1 == std::string("cutoff"))
        {
            iss >> cutoff;
        }
        else if (arg1 == std::string("dr"))
        {
            iss >> dr;
        }
        else if (arg1 == std::string("mesh"))
        {
            iss >> mesh;
        }
        else if (arg1 == std::string("l"))
        {
            continue;
        }
        if (arg1 == std::string("f:"))
        {
            break;
        }
    }
    getline(ifs, line);
    std::istringstream iss(line);
    char c;
    double a;
    for (int i = 0; i < mesh - 1; ++i)
    {
        iss >> a;
        iss >> c;
        vec.push_back(a);
    }
    iss >> a;
    vec.push_back(a);
    ifs.close();
    return vec;
}

void Input::read_v(double **p)
{
    std::ifstream ifs;
    ifs.open(v_path);
    std::string line;

    while (getline(ifs, line))
    {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        std::string arg1;

        iss >> arg1;
        if (arg1 == std::string("nx"))
        {
            iss >> nx;
        }
        else if (arg1 == std::string("ny"))
        {
            iss >> ny;
        }
        else if (arg1 == std::string("nz"))
        {
            iss >> nz;
        }
        else if (arg1 == std::string("V:"))
        {
            break;
        }
    }
    *p = (double *)malloc(sizeof(double) * (nx * ny * nz));
    double *pp = *p;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int l = 0; l < nz / 6; ++l)
            {
                getline(ifs, line);
                std::istringstream iss(line);
                double a;
                for (int k = 0; k < 6; ++k)
                {
                    iss >> a;
                    pp[i * ny * nz + j * nz + l * 6 + k] = a;
                }
            }
            if (nz - nz / 6 * 6 != 0)
            {
                getline(ifs, line);
                std::istringstream iss(line);
                double a;
                for (int k = 0; k < nz - nz / 6 * 6; ++k)
                {
                    iss >> a;
                    pp[i * ny * nz + j * nz + nz / 6 * 6 + k] = a;
                }
            }
        }
    }
}
