/*

This file contains all the functions used to calculated the time optimal 
gradient waveforms.

minTimeGradientRIV: Computes the rotationally invariant solution
    RungeKutte_riv: Used to solve the ODE using RK4
    beta: calculates sqrt (gamma^2 * smax^2 - k^2 * st^4) in the ODE. 
          Used in RungeKutte_riv

minTimeGradientRV: Computes the rotationally variant solution
    RungeKutte_rv: Used to solve the ODE using RK4
     sdotdot: calculates the maximum possible value for d^2s/dt^t, 
              used in RugeKutte_rv 
              
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mtg.h"

namespace py = pybind11;


/* 
function: riv
=============

Finds the time optimal gradient waveforms for the rotationally invariant 
constraints case.

Input
-----
g0: Initial gradient amplitude.
gfin: Gradient value at the end of the trajectory. 
      If given value is not possible, the result would be the largest possible 
      amplitude. Enter -1 for default.
gmax: Maximum gradient [G/cm] (4 default)
smax: Maximum slew [G/cm/ms] (15 default)
T: Sampling time intervale [ms] (4e-3 default)
ds: Step size for ODE integration. Enter -1 to use default.

Output
------
[gx gy gz]: gradient waveforms [G/cm]

Others
------
[Cx Cy Cz]: reparameterized curve, sampled at T[ms]
[sx sy sz]: slew rate [G/cm/ms]
[kx ky kz]: exact k-space corresponding to gradient g
sta: Solution for the forward ODE
stb: Solution for the backward ODE
phi: Geometry constrains on amplitude vs. arclength
 */
std::tuple< py::array_t<double>, py::array_t<double>, py::array_t<double> > 
riv(py::array_t<double> py_x, py::array_t<double> py_y, 
    py::array_t<double> py_z, int Lp, double g0, double gfin, double gmax, 
    double smax, double T, double ds) 
{
    /* from python to c */

    py::buffer_info buf_x = py_x.request();
    py::buffer_info buf_y = py_y.request();
    py::buffer_info buf_z = py_z.request();

    double *x = static_cast<double *>(buf_x.ptr);
    double *y = static_cast<double *>(buf_y.ptr);
    double *z = static_cast<double *>(buf_z.ptr);

    double *Cx, *Cy, *Cz;
    int l_t = riv_c(&Cx, &Cy, &Cz, x, y, z, Lp, g0, gfin, gmax, smax, T, ds);
    
    /* Final gradient waveforms to be returned */

    auto py_gx = py::array_t<double>(l_t);
    auto py_gy = py::array_t<double>(l_t);
    auto py_gz = py::array_t<double>(l_t);

    py::buffer_info buf_gx = py_gx.request();
    py::buffer_info buf_gy = py_gy.request();
    py::buffer_info buf_gz = py_gz.request();

    double *gx = static_cast<double *>(buf_gx.ptr);
    double *gy = static_cast<double *>(buf_gy.ptr);
    double *gz = static_cast<double *>(buf_gz.ptr);
    
    double dt = T;
    double gamma = 4.257;

    for (int i = 0 ; i < l_t - 1 ; i++) {
        gx[i] = (Cx[i+1] - Cx[i]) / (gamma * dt);
        gy[i] = (Cy[i+1] - Cy[i]) / (gamma * dt);
        gz[i] = (Cz[i+1] - Cz[i]) / (gamma * dt);
    }
    
    gx[l_t-1] = gx[l_t-2] + gx[l_t-2] - gx[l_t-3];
    gy[l_t-1] = gy[l_t-2] + gy[l_t-2] - gy[l_t-3];
    gz[l_t-1] = gz[l_t-2] + gz[l_t-2] - gz[l_t-3];
    
    free(Cx); free(Cy); free(Cz);
    
    std::tuple< py::array_t<double>, 
                py::array_t<double>, 
                py::array_t<double> > result(py_gx, py_gy, py_gz);

    return result;
}


/* 
function: rv
=============

Finds the time optimal gradient waveforms for the rotationally variant 
constraints case.

Input
-----
g0: Initial gradient amplitude.
gfin: Gradient value at the end of the trajectory. 
      If given value is not possible, the result would be the largest possible 
      amplitude. Enter -1 for default.
gmax: Maximum gradient [G/cm] (4 default)
smax: Maximum slew [G/cm/ms] (15 default)
T: Sampling time intervale [ms] (4e-3 default)
ds: Step size for ODE integration. Enter -1 to use default.

Output
------
[gx gy gz]: gradient waveforms [G/cm]

Others
------
[Cx Cy Cz]: reparameterized curve, sampled at T[ms]
[sx sy sz]: slew rate [G/cm/ms]
[kx ky kz]: exact k-space corresponding to gradient g
sta: Solution for the forward ODE
stb: Solution for the backward ODE
phi: Geometry constrains on amplitude vs. arclength
 */
std::tuple< py::array_t<double>, py::array_t<double>, py::array_t<double> > 
rv(py::array_t<double> py_x, py::array_t<double> py_y, 
    py::array_t<double> py_z, int Lp, double g0, double gfin, double gmax, 
    double smax, double T, double ds) 

{
    /* from python to c */

    py::buffer_info buf_x = py_x.request();
    py::buffer_info buf_y = py_y.request();
    py::buffer_info buf_z = py_z.request();

    double *x = static_cast<double *>(buf_x.ptr);
    double *y = static_cast<double *>(buf_y.ptr);
    double *z = static_cast<double *>(buf_z.ptr);

    double *Cx, *Cy, *Cz;
    int l_t = rv_c(&Cx, &Cy, &Cz, x, y, z, Lp, g0, gfin, gmax, smax, T, ds);
    
    /* Final gradient waveforms to be returned */

    auto py_gx = py::array_t<double>(l_t);
    auto py_gy = py::array_t<double>(l_t);
    auto py_gz = py::array_t<double>(l_t);

    py::buffer_info buf_gx = py_gx.request();
    py::buffer_info buf_gy = py_gy.request();
    py::buffer_info buf_gz = py_gz.request();

    double *gx = static_cast<double *>(buf_gx.ptr);
    double *gy = static_cast<double *>(buf_gy.ptr);
    double *gz = static_cast<double *>(buf_gz.ptr);
    
    double dt = T;
    double gamma = 4.257;

    for (int i = 0 ; i < l_t - 1 ; i++) {
        gx[i] = (Cx[i+1] - Cx[i]) / (gamma * dt);
        gy[i] = (Cy[i+1] - Cy[i]) / (gamma * dt);
        gz[i] = (Cz[i+1] - Cz[i]) / (gamma * dt);
    }
    
    gx[l_t-1] = gx[l_t-2] + gx[l_t-2] - gx[l_t-3];
    gy[l_t-1] = gy[l_t-2] + gy[l_t-2] - gy[l_t-3];
    gz[l_t-1] = gz[l_t-2] + gz[l_t-2] - gz[l_t-3];
        
    free(Cx); free(Cy); free(Cz);
    
    std::tuple< py::array_t<double>, 
                py::array_t<double>, 
                py::array_t<double> > result(py_gx, py_gy, py_gz);

    return result;
}

PYBIND11_MODULE(mintimegrad, m) {
    m.def("riv", &riv, 
          "Return time optimal gradient waveforms under rotationally invariant constraints.", 
          py::return_value_policy::move);

    m.def("rv", &rv, 
          "Return time optimal gradient waveforms under rotationally variant constraints.", 
          py::return_value_policy::move);
}