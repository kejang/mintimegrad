#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

#include "spline.h"
#include "mtg_c.h"

double *double_malloc_check(int L)
{
    /* Double array malloc */
    double *rtn = (double *)malloc(L * sizeof(double));
    if (!rtn)
    {
        printf("malloc failed");
        exit(1);
    }
    return rtn;
}

double sdotdot(double xpp, double ypp, double zpp, double xp, double yp,
               double zp, double st, double smax)
{
    /* Computes maximum possible value for sdotdot, the second time derivative
     * of s, the arc-length.
     *
     * xp, yp, zp are the first derivatives of the curve derivatives of the
     * curve in the arc-length(s) parameterization.
     *
     * xpp, ypp, and zpp are the second derivatives.
     */
    double gamma = 4.257;
    double sx, sy, sz; /* constraints for sdotdot in x, y, and z directions */

    /* maximum sdotdot will be limited by the smallest of these */
    sx = (-xpp * st * st + gamma * smax) / xp;
    sy = (-ypp * st * st + gamma * smax) / yp;
    sz = (-zpp * st * st + gamma * smax) / zp;
    double rtn = sx;
    if (sy < rtn)
    {
        rtn = sy;
    }
    if (sz < rtn)
    {
        rtn = sz;
    }
    return rtn;
}

double beta(double k, double st, double smax)
{
    /* calculates sqrt (gamma^2 * smax^2 - k^2 * st^4) used in RK4 method
     * for rotationally invariant ODE solver.
     */
    double gamma = 4.257;
    return sqrt(sqrt((gamma * gamma * smax * smax - k * k * st * st * st * st) * (gamma * gamma * smax * smax - k * k * st * st * st * st)));
}

double RungeKutte_riv(double ds, double st, double k[], double smax)
{
    /* Solves ODE for rotationally invariant solution using Runge-Kutte */
    double k1 = ds * (1 / st) * beta(k[0], st, smax);
    double k2 = ds * 1 / (st + k1 / 2) * beta(k[1], st + k1 / 2, smax);
    double k3 = ds * 1 / (st + k2 / 2) * beta(k[1], st + k2 / 2, smax);
    double k4 = ds * 1 / (st + k3 / 2) * beta(k[2], st + k3 / 2, smax);
    double rtn = k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
    return rtn;
}

double RungeKutte_rv(double ds, double st, double xpp[], double ypp[],
                     double zpp[], double xp[], double yp[], double zp[],
                     double smax)
{
    /*  Solves ODE for rotationally variant solution using Runge-Kutte */
    double k1, k2, k3, k4;
    k1 = ds * (1 / st) * sdotdot(xpp[0], ypp[0], zpp[0], xp[0], yp[0], zp[0], st, smax);

    k2 = ds * (1 / (st + k1 / 2)) * sdotdot(xpp[1], ypp[1], zpp[1], xp[1], yp[1], zp[1], st + k1 / 2, smax);

    k3 = ds * (1 / (st + k2 / 2)) * sdotdot(xpp[1], ypp[1], zpp[1], xp[1], yp[1], zp[1], st + k2 / 2, smax);

    k4 = ds * (1 / (st + k3 / 2)) * sdotdot(xpp[2], ypp[2], zpp[2], xp[2], yp[2], zp[2], st + k3 / 2, smax);

    double rtn = k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
    return rtn;
}

int riv_c(double **Cx, double **Cy, double **Cz, double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax,
          double smax, double T, double ds)
{
    /* iflag used in spline method to signal error */

    int *iflag;
    int iflagp;
    iflag = &iflagp;

    /* Representing the curve with parameter p */

    double p[Lp];

    for (int i = 0; i < Lp; i++)
    {
        p[i] = i;
    }

    double dt = T;
    double gamma = 4.257;

    /* Interpolation of curve for gradient accuracy, using cubic spline
       interpolation */

    double *c1x, *c2x, *c3x,
        *c1y, *c2y, *c3y,
        *c1z, *c2z, *c3z; /* storage for spline coefficients */

    c1z = double_malloc_check(Lp * sizeof(double));
    c2z = double_malloc_check(Lp * sizeof(double));
    c3z = double_malloc_check(Lp * sizeof(double));

    c1x = double_malloc_check(Lp * sizeof(double));
    c2x = double_malloc_check(Lp * sizeof(double));
    c3x = double_malloc_check(Lp * sizeof(double));

    c1y = double_malloc_check(Lp * sizeof(double));
    c2y = double_malloc_check(Lp * sizeof(double));
    c3y = double_malloc_check(Lp * sizeof(double));

    spline(Lp, 0, 0, 1, 1, p, x, c1x, c2x, c3x, iflag);
    spline(Lp, 0, 0, 1, 1, p, y, c1y, c2y, c3y, iflag);
    spline(Lp, 0, 0, 1, 1, p, z, c1z, c2z, c3z, iflag);

    double dp = 0.1;
    int num_evals = (int)floor((Lp - 1) / dp) + 1;

    double *CCx, *CCy, *CCz;
    CCx = double_malloc_check(num_evals * sizeof(double));
    CCy = double_malloc_check(num_evals * sizeof(double));
    CCz = double_malloc_check(num_evals * sizeof(double));

    double toeval = 0;

    int *last;
    int holder = 0;
    last = &holder;

    /* interpolated curve in p-parameterization */

    double *Cpx, *Cpy, *Cp_abs, *Cpz;
    Cpx = double_malloc_check(num_evals * sizeof(double));
    Cpy = double_malloc_check(num_evals * sizeof(double));
    Cpz = double_malloc_check(num_evals * sizeof(double));
    Cp_abs = double_malloc_check(num_evals * sizeof(double));

    for (int i = 0; i < num_evals; i++)
    {
        toeval = (double)i * dp;
        CCx[i] = seval(Lp, toeval, p, x, c1x, c2x, c3x, last);
        CCy[i] = seval(Lp, toeval, p, y, c1y, c2y, c3y, last);
        CCz[i] = seval(Lp, toeval, p, z, c1z, c2z, c3z, last);
        Cpx[i] = deriv(Lp, toeval, p, c1x, c2x, c3x, last);
        Cpy[i] = deriv(Lp, toeval, p, c1y, c2y, c3y, last);
        Cpz[i] = deriv(Lp, toeval, p, c1z, c2z, c3z, last);
        Cp_abs[i] = sqrt(Cpx[i] * Cpx[i] + Cpy[i] * Cpy[i] + Cpz[i] * Cpz[i]);
    }
    free(Cpx);
    free(Cpy);
    free(Cpz);

    /* converting to arc-length parameterization from p, using trapezoidal
       integration */

    double *s_of_p;
    s_of_p = double_malloc_check(num_evals * sizeof(double));
    s_of_p[0] = 0;

    double sofar = 0;

    for (int i = 0; i < num_evals; i++)
    {
        sofar += (Cp_abs[i] + Cp_abs[i - 1]) / 2;
        s_of_p[i] = dp * sofar;
    }

    free(Cp_abs);

    /* length of the curve */
    double L = s_of_p[num_evals - 1];

    /* decide ds and compute st for the first point */
    double stt0 = gamma * smax;   /* always assumes first point is max slew */
    double st0 = (stt0 * dt) / 2; /* start at half the gradient for accuracy
                                     close to g=0 */
    double s0 = st0 * dt;
    if (ds < 0)
    {
        ds = s0 / 1.5; /* smaller step size for numerical accuracy */
    }

    int length_of_s = (int)floor(L / ds);
    int half_ls = (int)floor(L / (ds / 2));

    double *s;
    s = double_malloc_check(length_of_s * sizeof(double));

    double *sta;
    sta = double_malloc_check(length_of_s * sizeof(double));

    double *stb;
    stb = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        s[i] = i * ds;
        sta[i] = 0;
        stb[i] = 0;
    }

    double *s_half;
    s_half = double_malloc_check(half_ls * sizeof(double));

    for (int i = 0; i < half_ls; i++)
    {
        s_half[i] = (double)i * (ds / 2);
    }

    double *p_of_s_half;
    p_of_s_half = double_malloc_check(half_ls * sizeof(double));

    /* Convert from s(p) to p(s) and interpolate for accuracy */

    double *a1x, *a2x, *a3x;
    a1x = double_malloc_check(num_evals * sizeof(double));
    a2x = double_malloc_check(num_evals * sizeof(double));
    a3x = double_malloc_check(num_evals * sizeof(double));

    double sop_num[num_evals];
    for (int i = 0; i < num_evals; i++)
    {
        sop_num[i] = i * dp;
    }

    spline(num_evals, 0, 0, 1, 1, s_of_p, sop_num, a1x, a2x, a3x, iflag);

    for (int i = 0; i < half_ls; i++)
    {
        p_of_s_half[i] = seval(num_evals, s_half[i], s_of_p, sop_num, a1x, a2x,
                               a3x, last);
    }

    free(a1x);
    free(a2x);
    free(a3x);
    free(s_of_p);

    int size_p_of_s = half_ls / 2;

    double *p_of_s;
    p_of_s = double_malloc_check(size_p_of_s * sizeof(double));

    for (int i = 0; i < size_p_of_s; i++)
    {
        p_of_s[i] = p_of_s_half[2 * i];
    }

    free(p_of_s_half);

    double *k; /* curvature along the curve */
    k = double_malloc_check(half_ls * sizeof(double));

    double *Cspx, *Cspy, *Cspz;

    /* Csp is C(s(p)) = [Cx(p(s)) Cy(p(s)) Cz(p(s))] */
    Cspx = double_malloc_check(length_of_s * sizeof(double));
    Cspy = double_malloc_check(length_of_s * sizeof(double));
    Cspz = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        Cspx[i] = seval(Lp, p_of_s[i], p, x, c1x, c2x, c3x, last);
        Cspy[i] = seval(Lp, p_of_s[i], p, y, c1y, c2y, c3y, last);
        Cspz[i] = seval(Lp, p_of_s[i], p, z, c1z, c2z, c3z, last);
    }

    double *Csp1x, *Csp2x, *Csp3x,
        *Csp1y, *Csp2y, *Csp3y,
        *Csp1z, *Csp2z, *Csp3z; /* storage for spline coefficients. */

    Csp1x = double_malloc_check(length_of_s * sizeof(double));
    Csp2x = double_malloc_check(length_of_s * sizeof(double));
    Csp3x = double_malloc_check(length_of_s * sizeof(double));
    Csp1y = double_malloc_check(length_of_s * sizeof(double));
    Csp2y = double_malloc_check(length_of_s * sizeof(double));
    Csp3y = double_malloc_check(length_of_s * sizeof(double));
    Csp1z = double_malloc_check(length_of_s * sizeof(double));
    Csp2z = double_malloc_check(length_of_s * sizeof(double));
    Csp3z = double_malloc_check(length_of_s * sizeof(double));

    spline(length_of_s, 0, 0, 1, 1, s, Cspx, Csp1x, Csp2x, Csp3x, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspy, Csp1y, Csp2y, Csp3y, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspz, Csp1z, Csp2z, Csp3z, iflag);

    free(Cspx);
    free(Cspy);
    free(Cspz);

    for (int i = 0; i < half_ls; i++)
    {
        double kx = deriv2(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x,
                           last);

        double ky = deriv2(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y,
                           last);

        double kz = deriv2(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z,
                           last);

        /* the curvature, magnitude of the second derivative of the curve
           in arc-length parameterization */
        k[i] = sqrt(kx * kx + ky * ky + kz * kz);
    }

    free(s_half);
    free(CCx);
    free(CCy);
    free(CCz);

    free(Csp1x);
    free(Csp2x);
    free(Csp3x);
    free(Csp1y);
    free(Csp2y);
    free(Csp3y);
    free(Csp1z);
    free(Csp2z);
    free(Csp3z);

    /* computing geomtry dependent constraints (forbidden line curve) */

    double *sdot1, *sdot2, *sdot;

    sdot1 = double_malloc_check(half_ls * sizeof(double));
    sdot2 = double_malloc_check(half_ls * sizeof(double));
    sdot = double_malloc_check(half_ls * sizeof(double));

    /* Calculating the upper bound for the time parametrization */
    /* sdot (which is a non scaled max gradient constaint) as a function of s. */
    /* sdot is the minimum of gamma*gmax and sqrt(gamma*gmax / k) */

    for (int i = 0; i < half_ls; i++)
    {
        sdot1[i] = gamma * gmax;
        sdot2[i] = sqrt((gamma * smax) / (fabs(k[i] + (DBL_EPSILON))));
        if (sdot1[i] < sdot2[i])
        {
            sdot[i] = sdot1[i];
        }
        else
        {
            sdot[i] = sdot2[i];
        }
    }

    free(sdot1);
    free(sdot2);

    int size_k2 = half_ls + 2; /* extend of k for RK4 */

    double *k2;
    k2 = double_malloc_check(size_k2 * sizeof(double));
    for (int i = 0; i < half_ls; i++)
    {
        k2[i] = k[i];
    }

    k2[size_k2 - 2] = k2[size_k2 - 3];
    k2[size_k2 - 1] = k2[size_k2 - 3];

    double g0gamma = g0 * gamma + st0;
    double gammagmax = gamma * gmax;

    if (g0gamma < gammagmax)
    {
        sta[0] = g0gamma;
    }
    else
    {
        sta[0] = gammagmax;
    }

    /* Solving ODE Forward */

    for (int i = 1; i < length_of_s; i++)
    {
        double k_rk[3];
        k_rk[0] = k2[2 * i - 2];
        k_rk[1] = k2[2 * i - 1];
        k_rk[2] = k2[2 * i];

        double dstds = RungeKutte_riv(ds, sta[i - 1], k_rk, smax);
        double tmpst = sta[i - 1] + dstds;
        if (sdot[2 * i + 1] < tmpst)
        {
            sta[i] = sdot[2 * i + 1];
        }
        else
        {
            sta[i] = tmpst;
        }
    }

    free(k2);

    /* Solving ODE Backwards */

    double max;
    if (gfin < 0)
    { /*if gfin is not provided */
        stb[length_of_s - 1] = sta[length_of_s - 1];
    }
    else
    {
        if (gfin * gamma > st0)
        {
            max = gfin * gamma;
        }
        else
        {
            max = st0;
        }

        if (gamma * gmax < max)
        {
            stb[length_of_s - 1] = gamma * gmax;
        }
        else
        {
            stb[length_of_s - 1] = max;
        }
    }

    for (int i = length_of_s - 2; i > -1; i--)
    {
        double k_rk[3];
        k_rk[0] = k[2 * i + 2];
        k_rk[1] = k[2 * i + 1];
        k_rk[2] = k[2 * i];

        double dstds = RungeKutte_riv(ds, stb[i + 1], k_rk, smax);
        double tmpst = stb[i + 1] + dstds;

        if (sdot[2 * i] < tmpst)
        {
            stb[i] = sdot[2 * i];
        }
        else
        {
            stb[i] = tmpst;
        }
    }

    free(k);
    free(sdot);

    /* take st(s) to be the minimum of the curves sta and stb */
    double *st_of_s, *st_ds_i;
    st_of_s = double_malloc_check(length_of_s * sizeof(double));
    st_ds_i = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        if (sta[i] < stb[i])
        {
            st_of_s[i] = sta[i];
        }
        else
        {
            st_of_s[i] = stb[i];
        }

        /* ds * 1/st(s) used in below calculation of t(s) */
        st_ds_i[i] = ds * (1 / st_of_s[i]);
    }

    free(st_of_s);
    free(sta);
    free(stb);

    /* Final interpolation */

    /* Converting to the time parameterization, t(s) using trapezoidal
       integration. t(s) = integral (1/st) ds */
    double *t_of_s;
    t_of_s = double_malloc_check(length_of_s * sizeof(double));

    t_of_s[0] = 0;
    for (int i = 1; i < length_of_s; i++)
    {
        t_of_s[i] = t_of_s[i - 1] + (st_ds_i[i] + st_ds_i[i - 1]) / 2;
    }

    free(st_ds_i);

    int l_t = (int)floor(t_of_s[length_of_s - 1] / dt);

    double t[l_t];
    for (int i = 0; i < l_t; i++)
    {
        t[i] = i * dt; /* time array */
    }

    double *t1x, *t2x, *t3x; /* coefficient arrays for spline interpolation
                                of t(s) to get s(t) */

    t1x = double_malloc_check(length_of_s * sizeof(double));
    t2x = double_malloc_check(length_of_s * sizeof(double));
    t3x = double_malloc_check(length_of_s * sizeof(double));

    double *s_of_t;
    s_of_t = double_malloc_check(l_t * sizeof(double));

    spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);

    for (int i = 0; i < l_t; i++)
    {
        s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
    }

    free(t1x);
    free(t2x);
    free(t3x);
    free(t_of_s);

    double *p1x, *p2x, *p3x; /* coefficient arrays for spline interpolation
                                of p(s) with s(t) to get p(s(t)) = p(t) */
    p1x = double_malloc_check(length_of_s * sizeof(double));
    p2x = double_malloc_check(length_of_s * sizeof(double));
    p3x = double_malloc_check(length_of_s * sizeof(double));

    spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);

    double *p_of_t;
    p_of_t = double_malloc_check(l_t * sizeof(double));

    for (int i = 0; i < l_t; i++)
    {
        p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x,
                          last);
    }

    free(s);
    free(p_of_s);
    free(p1x);
    free(p2x);
    free(p3x);
    free(s_of_t);

    /* interpolated k-space trajectory */

    *Cx = double_malloc_check(l_t * sizeof(double));
    *Cy = double_malloc_check(l_t * sizeof(double));
    *Cz = double_malloc_check(l_t * sizeof(double));

    for (int i = 0; i < l_t; i++)
    {
        (*Cx)[i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
        (*Cy)[i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
        (*Cz)[i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
    }

    free(p_of_t);
    free(c1x);
    free(c2x);
    free(c3x);
    free(c1y);
    free(c2y);
    free(c3y);
    free(c1z);
    free(c2z);
    free(c3z);

    return l_t;
}

int rv_c(double **Cx, double **Cy, double **Cz, double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax,
         double smax, double T, double ds)
{
    /* iflag used in spline method to signal error */

    int *iflag;
    int iflagp;
    iflag = &iflagp;

    /* Representing the curve with parameter p */

    double p[Lp];
    for (int i = 0; i < Lp; i++)
    {
        p[i] = i;
    }

    double dt = T;
    double gamma = 4.257;

    /* Interpolation of curve in p-parameterization for gradient accuracy,
       using cubic spline interpolation */
    double *c1x, *c2x, *c3x,
        *c1y, *c2y, *c3y,
        *c1z, *c2z, *c3z; /* storage for spline coefficients */

    c1z = double_malloc_check(Lp * sizeof(double));
    c2z = double_malloc_check(Lp * sizeof(double));
    c3z = double_malloc_check(Lp * sizeof(double));
    c1x = double_malloc_check(Lp * sizeof(double));
    c2x = double_malloc_check(Lp * sizeof(double));
    c3x = double_malloc_check(Lp * sizeof(double));
    c1y = double_malloc_check(Lp * sizeof(double));
    c2y = double_malloc_check(Lp * sizeof(double));
    c3y = double_malloc_check(Lp * sizeof(double));

    spline(Lp, 0, 0, 1, 1, p, x, c1x, c2x, c3x, iflag);
    spline(Lp, 0, 0, 1, 1, p, y, c1y, c2y, c3y, iflag);
    spline(Lp, 0, 0, 1, 1, p, z, c1z, c2z, c3z, iflag);

    double dp = 0.1;
    int num_evals = (int)floor((Lp - 1) / dp) + 1;

    double *CCx, *CCy, *CCz;
    CCx = double_malloc_check(num_evals * sizeof(double));
    CCy = double_malloc_check(num_evals * sizeof(double));
    CCz = double_malloc_check(num_evals * sizeof(double));

    double toeval = 0; /* used by spline eval function, seval */

    int *last;
    int holder = 0;
    last = &holder;

    /* interpolated curve in p-parameterization */

    double *Cpx, *Cpy, *Cp_abs, *Cpz;
    Cpx = double_malloc_check(num_evals * sizeof(double));
    Cpy = double_malloc_check(num_evals * sizeof(double));
    Cpz = double_malloc_check(num_evals * sizeof(double));
    Cp_abs = double_malloc_check(num_evals * sizeof(double));

    for (int i = 0; i < num_evals; i++)
    {
        toeval = (double)i * dp;
        CCx[i] = seval(Lp, toeval, p, x, c1x, c2x, c3x, last);
        CCy[i] = seval(Lp, toeval, p, y, c1y, c2y, c3y, last);
        CCz[i] = seval(Lp, toeval, p, z, c1z, c2z, c3z, last);
        Cpx[i] = deriv(Lp, toeval, p, c1x, c2x, c3x, last);
        Cpy[i] = deriv(Lp, toeval, p, c1y, c2y, c3y, last);
        Cpz[i] = deriv(Lp, toeval, p, c1z, c2z, c3z, last);
        Cp_abs[i] = sqrt(Cpx[i] * Cpx[i] + Cpy[i] * Cpy[i] + Cpz[i] * Cpz[i]);
    }
    free(Cpx);
    free(Cpy);
    free(Cpz);
    free(CCx);
    free(CCy);
    free(CCz);

    /* converting to arc-length parameterization from p, using trapezoidal
       integration */

    double *s_of_p;
    s_of_p = double_malloc_check(num_evals * sizeof(double));
    s_of_p[0] = 0;

    double sofar = 0;

    for (int i = 0; i < num_evals; i++)
    {
        sofar += (Cp_abs[i] + Cp_abs[i - 1]) / 2;
        s_of_p[i] = dp * sofar;
    }

    free(Cp_abs);

    /* length of the curve */
    double L = s_of_p[num_evals - 1];

    /* decide ds and compute st for the first point */
    double stt0 = gamma * smax;   /* always assumes first point is max slew */
    double st0 = (stt0 * dt) / 2; /* start at half the gradient for accuracy
                                     close to g=0 */
    double s0 = st0 * dt;
    if (ds < 0)
    {
        ds = s0 / 1.5; /* smaller step size for numerical accuracy */
    }

    int length_of_s = (int)floor(L / ds);
    int half_ls = (int)floor(L / (ds / 2));

    double *s;
    s = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        s[i] = i * ds;
    }

    double *s_half;
    s_half = double_malloc_check(half_ls * sizeof(double));

    for (int i = 0; i < half_ls; i++)
    {
        s_half[i] = (double)i * (ds / 2);
    }

    double *p_of_s_half;
    p_of_s_half = double_malloc_check(half_ls * sizeof(double));

    /* Convert from s(p) to p(s) and interpolate for accuracy */

    double *a1x, *a2x, *a3x;
    a1x = double_malloc_check(num_evals * sizeof(double));
    a2x = double_malloc_check(num_evals * sizeof(double));
    a3x = double_malloc_check(num_evals * sizeof(double));

    double sop_num[num_evals];
    for (int i = 0; i < num_evals; i++)
    {
        sop_num[i] = i * dp;
    }
    spline(num_evals, 0, 0, 1, 1, s_of_p, sop_num, a1x, a2x, a3x, iflag);

    for (int i = 0; i < half_ls; i++)
    {
        p_of_s_half[i] = seval(num_evals, s_half[i], s_of_p, sop_num, a1x, a2x,
                               a3x, last);
    }

    free(a1x);
    free(a2x);
    free(a3x);
    free(s_of_p);

    int size_p_of_s = half_ls / 2;

    double *p_of_s;
    p_of_s = double_malloc_check(size_p_of_s * sizeof(double));

    for (int i = 0; i < size_p_of_s; i++)
    {
        p_of_s[i] = p_of_s_half[2 * i];
    }
    free(p_of_s_half);

    double *k; /* k is the curvature along the curve */
    k = double_malloc_check(half_ls * sizeof(double));

    double *Cspx, *Cspy, *Cspz;

    /* Csp is C(s(p)) = [Cx(p(s)) Cy(p(s)) Cz(p(s))]  */
    Cspx = double_malloc_check(length_of_s * sizeof(double));
    Cspy = double_malloc_check(length_of_s * sizeof(double));
    Cspz = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        Cspx[i] = seval(Lp, p_of_s[i], p, x, c1x, c2x, c3x, last);
        Cspy[i] = seval(Lp, p_of_s[i], p, y, c1y, c2y, c3y, last);
        Cspz[i] = seval(Lp, p_of_s[i], p, z, c1z, c2z, c3z, last);
    }

    double *Csp1x, *Csp2x, *Csp3x,
        *Csp1y, *Csp2y, *Csp3y,
        *Csp1z, *Csp2z, *Csp3z; /* spline coefficients */

    Csp1x = double_malloc_check(length_of_s * sizeof(double));
    Csp2x = double_malloc_check(length_of_s * sizeof(double));
    Csp3x = double_malloc_check(length_of_s * sizeof(double));
    Csp1y = double_malloc_check(length_of_s * sizeof(double));
    Csp2y = double_malloc_check(length_of_s * sizeof(double));
    Csp3y = double_malloc_check(length_of_s * sizeof(double));
    Csp1z = double_malloc_check(length_of_s * sizeof(double));
    Csp2z = double_malloc_check(length_of_s * sizeof(double));
    Csp3z = double_malloc_check(length_of_s * sizeof(double));

    /* interpolation for C(s(p)) */
    spline(length_of_s, 0, 0, 1, 1, s, Cspx, Csp1x, Csp2x, Csp3x, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspy, Csp1y, Csp2y, Csp3y, iflag);
    spline(length_of_s, 0, 0, 1, 1, s, Cspz, Csp1z, Csp2z, Csp3z, iflag);

    free(Cspx);
    free(Cspy);
    free(Cspz);

    double *xpp, *ypp, *zpp; /* xpp, ypp, zpp are d^2/ds^2
                                (second derivative in s parameterization) */

    xpp = double_malloc_check(half_ls * sizeof(double));
    ypp = double_malloc_check(half_ls * sizeof(double));
    zpp = double_malloc_check(half_ls * sizeof(double));

    double *xp, *yp, *zp; /* d/ds - first derivative in s parameterization */
    xp = double_malloc_check(half_ls * sizeof(double));
    yp = double_malloc_check(half_ls * sizeof(double));
    zp = double_malloc_check(half_ls * sizeof(double));

    /* Computing the curvature along the curve */
    for (int i = 0; i < half_ls; i++)
    {
        double kx = deriv2(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x,
                           last);

        double ky = deriv2(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y,
                           last);

        double kz = deriv2(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z,
                           last);

        zpp[i] = sqrt(kz * kz);
        xpp[i] = sqrt(kx * kx);
        ypp[i] = sqrt(ky * ky);

        xp[i] = sqrt(
            deriv(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last) * deriv(length_of_s, s_half[i], s, Csp1x, Csp2x, Csp3x, last));

        yp[i] = sqrt(
            deriv(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last) * deriv(length_of_s, s_half[i], s, Csp1y, Csp2y, Csp3y, last));

        zp[i] = sqrt(
            deriv(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last) * deriv(length_of_s, s_half[i], s, Csp1z, Csp2z, Csp3z, last));

        /* the curvature, magnitude of the second derivative of the curve in
           arc-length parameterization */
        k[i] = sqrt(kx * kx + ky * ky + kz * kz);
    }

    free(s_half);
    free(Csp1x);
    free(Csp2x);
    free(Csp3x);
    free(Csp1y);
    free(Csp2y);
    free(Csp3y);
    free(Csp1z);
    free(Csp2z);
    free(Csp3z);

    /* Compute the geometry dependant constraints */
    double *sdot1, *sdot2, *sdot3;
    sdot1 = double_malloc_check(half_ls * sizeof(double));
    sdot2 = double_malloc_check(half_ls * sizeof(double));
    sdot3 = double_malloc_check(half_ls * sizeof(double));

    double *sdot;
    sdot = double_malloc_check(half_ls * sizeof(double));

    for (int i = 0; i < half_ls; i++)
    {
        sdot1[i] = gamma * gmax / xp[i];
        sdot2[i] = gamma * gmax / yp[i];
        sdot3[i] = gamma * gmax / zp[i];
        sdot[i] = sdot1[i];

        if (sdot2[i] < sdot[i])
        {
            sdot[i] = sdot2[i];
        }
        if (sdot3[i] < sdot[i])
        {
            sdot[i] = sdot3[i];
        }
    }

    free(sdot1);
    free(sdot2);
    free(sdot3);

    /* Extend k, xp, yp, zp, xpp, ypp, zpp for RK4 end points */
    int size_k2 = half_ls + 2;
    double *k2;
    k2 = double_malloc_check(size_k2 * sizeof(double));

    for (int i = 0; i < half_ls; i++)
    {
        k2[i] = k[i];
    }
    double *xpp2, *ypp2, *zpp2, *xp2, *yp2, *zp2;
    xpp2 = double_malloc_check(size_k2 * sizeof(double));
    ypp2 = double_malloc_check(size_k2 * sizeof(double));
    zpp2 = double_malloc_check(size_k2 * sizeof(double));
    xp2 = double_malloc_check(size_k2 * sizeof(double));
    yp2 = double_malloc_check(size_k2 * sizeof(double));
    zp2 = double_malloc_check(size_k2 * sizeof(double));

    for (int i = 0; i < half_ls; i++)
    {
        xpp2[i] = xpp[i];
        ypp2[i] = ypp[i];
        zpp2[i] = zpp[i];
        xp2[i] = xp[i];
        yp2[i] = yp[i];
        zp2[i] = zp[i];
    }

    free(xpp);
    free(ypp);
    free(zpp);
    free(xp);
    free(yp);
    free(zp);

    xpp2[size_k2 - 2] = xpp2[size_k2 - 3];
    xpp2[size_k2 - 1] = xpp2[size_k2 - 3];

    ypp2[size_k2 - 2] = ypp2[size_k2 - 3];
    ypp2[size_k2 - 1] = ypp2[size_k2 - 3];

    zpp2[size_k2 - 2] = zpp2[size_k2 - 3];
    zpp2[size_k2 - 1] = zpp2[size_k2 - 3];

    xp2[size_k2 - 2] = xp2[size_k2 - 3];
    xp2[size_k2 - 1] = xp2[size_k2 - 3];

    yp2[size_k2 - 2] = yp2[size_k2 - 3];
    yp2[size_k2 - 1] = yp2[size_k2 - 3];

    zp2[size_k2 - 2] = zp2[size_k2 - 3];
    zp2[size_k2 - 1] = zp2[size_k2 - 3];

    k2[size_k2 - 2] = k2[size_k2 - 3];
    k2[size_k2 - 1] = k2[size_k2 - 3];

    /* Solving the ODE */

    double *sta; /* forward solution */
    sta = double_malloc_check(length_of_s * sizeof(double));

    double *stb; /* backward solution */
    stb = double_malloc_check(length_of_s * sizeof(double));

    /* Forward solution  */
    /* Initial cond.  */
    double g0gamma = g0 * gamma + st0;
    double gammagmax = gamma * gmax;
    if (g0gamma < gammagmax)
    {
        sta[0] = g0gamma;
    }
    else
    {
        sta[0] = gammagmax;
    }

    for (int i = 1; i < length_of_s; i++)
    {
        double xpp_rk[3], ypp_rk[3], zpp_rk[3], xp_rk[3], yp_rk[3], zp_rk[3];

        xpp_rk[0] = xpp2[2 * i - 2];
        xpp_rk[1] = xpp2[2 * i - 1];
        xpp_rk[2] = xpp2[2 * i];

        ypp_rk[0] = ypp2[2 * i - 2];
        ypp_rk[1] = ypp2[2 * i - 1];
        ypp_rk[2] = ypp2[2 * i];

        zpp_rk[0] = zpp2[2 * i - 2];
        zpp_rk[1] = zpp2[2 * i - 1];
        zpp_rk[2] = zpp2[2 * i];

        xp_rk[0] = xp2[2 * i - 2];
        xp_rk[1] = xp2[2 * i - 1];
        xp_rk[2] = xp2[2 * i];

        yp_rk[0] = yp2[2 * i - 2];
        yp_rk[1] = yp2[2 * i - 1];
        yp_rk[2] = yp2[2 * i];

        zp_rk[0] = zp2[2 * i - 2];
        zp_rk[1] = zp2[2 * i - 1];
        zp_rk[2] = zp2[2 * i];

        double dstds = RungeKutte_rv(ds, sta[i - 1], xpp_rk, ypp_rk, zpp_rk,
                                     xp_rk, yp_rk, zp_rk, smax);

        double tmpst = sta[i - 1] + dstds;
        if (sdot[2 * i] < tmpst)
        {
            sta[i] = sdot[2 * i];
        }
        else
        {
            sta[i] = tmpst;
        }
    }

    free(k2);

    if (gfin < 0)
    { /*if gfin is not provided */
        stb[length_of_s - 1] = sta[length_of_s - 1];
    }
    else
    {
        double max;
        if (gfin * gamma > st0)
        {
            max = gfin * gamma;
        }
        else
        {
            max = st0;
        }

        if (gamma * gmax < max)
        {
            stb[length_of_s - 1] = gamma * gmax;
        }
        else
        {
            stb[length_of_s - 1] = max;
        }
    }
    for (int i = length_of_s - 2; i > -1; i--)
    {
        double xpp_rk[3], ypp_rk[3], zpp_rk[3], xp_rk[3], yp_rk[3], zp_rk[3];

        xpp_rk[0] = xpp2[2 * i + 2];
        xpp_rk[1] = xpp2[2 * i + 1];
        xpp_rk[2] = xpp2[2 * i];

        ypp_rk[0] = ypp2[2 * i + 2];
        ypp_rk[1] = ypp2[2 * i + 1];
        ypp_rk[2] = ypp2[2 * i];

        zpp_rk[0] = zpp2[2 * i + 2];
        zpp_rk[1] = zpp2[2 * i + 1];
        zpp_rk[2] = zpp2[2 * i];

        xp_rk[0] = xp2[2 * i + 2];
        xp_rk[1] = xp2[2 * i + 1];
        xp_rk[2] = xp2[2 * i];

        yp_rk[0] = yp2[2 * i + 2];
        yp_rk[1] = yp2[2 * i + 1];
        yp_rk[2] = yp2[2 * i];

        zp_rk[0] = zp2[2 * i + 2];
        zp_rk[1] = zp2[2 * i + 1];
        zp_rk[2] = zp2[2 * i];

        double dstds = RungeKutte_rv(ds, stb[i + 1], xpp_rk, ypp_rk, zpp_rk,
                                     xp_rk, yp_rk, zp_rk, smax);

        double tmpst = stb[i + 1] + dstds;

        if (sdot[2 * i] < tmpst)
        {
            stb[i] = sdot[2 * i];
        }
        else
        {
            stb[i] = tmpst;
        }
    }

    free(k);
    free(sdot);
    free(xp2);
    free(yp2);
    free(zp2);
    free(xpp2);
    free(ypp2);
    free(zpp2);

    double *st_of_s, *st_ds_i;
    st_of_s = double_malloc_check(length_of_s * sizeof(double));
    st_ds_i = double_malloc_check(length_of_s * sizeof(double));

    for (int i = 0; i < length_of_s; i++)
    {
        if (sta[i] < stb[i])
        {
            /* use the minimum of the two solutions, sta and stb */
            st_of_s[i] = sta[i];
        }
        else
        {
            st_of_s[i] = stb[i];
        }
        /* ds * 1/st(s) used in below calculation of t(s) */
        st_ds_i[i] = ds * (1 / st_of_s[i]);
    }

    free(sta);
    free(stb);
    free(st_of_s);

    /* st is v(s) */

    /* Converting to the time parameterization, t(s) using trapezoidal
       integration. t(s) = integral (1/st) ds  */
    double *t_of_s = double_malloc_check(length_of_s * sizeof(double));
    t_of_s[0] = 0;
    for (int i = 1; i < length_of_s; i++)
    {
        t_of_s[i] = t_of_s[i - 1] + (st_ds_i[i] + st_ds_i[i - 1]) / 2;
    }

    free(st_ds_i);

    /* size of the interpolated trajectory  */
    int l_t = (int)floor(t_of_s[length_of_s - 1] / dt);

    double t[l_t];
    for (int i = 0; i < l_t; i++)
    {
        t[i] = i * dt; /* time array  */
    }

    double *t1x, *t2x, *t3x; /* coefficient arrays for spline interpolation
                                of t(s) to get s(t)  */
    t1x = double_malloc_check(length_of_s * sizeof(double));
    t2x = double_malloc_check(length_of_s * sizeof(double));
    t3x = double_malloc_check(length_of_s * sizeof(double));

    double *s_of_t;
    s_of_t = double_malloc_check(l_t * sizeof(double));
    spline(length_of_s, 0, 0, 1, 1, t_of_s, s, t1x, t2x, t3x, iflag);

    for (int i = 0; i < l_t; i++)
    {
        s_of_t[i] = seval(length_of_s, t[i], t_of_s, s, t1x, t2x, t3x, last);
    }
    free(t1x);
    free(t2x);
    free(t3x);
    free(t_of_s);

    double *p1x, *p2x, *p3x; /* coefficient arrays for spline interpolation
                                of p(s) with s(t) to get p(s(t)) = p(t) */

    p1x = double_malloc_check(length_of_s * sizeof(double));
    p2x = double_malloc_check(length_of_s * sizeof(double));
    p3x = double_malloc_check(length_of_s * sizeof(double));

    spline(length_of_s, 0, 0, 1, 1, s, p_of_s, p1x, p2x, p3x, iflag);

    double *p_of_t;
    p_of_t = double_malloc_check(l_t * sizeof(double));

    for (int i = 0; i < l_t; i++)
    {
        p_of_t[i] = seval(length_of_s, s_of_t[i], s, p_of_s, p1x, p2x, p3x,
                          last);
    }

    free(s);
    free(p_of_s);
    free(p1x);
    free(p2x);
    free(p3x);
    free(s_of_t);

    /*  Reparameterized curve C, sampled at T */

    *Cx = double_malloc_check(l_t * sizeof(double));
    *Cy = double_malloc_check(l_t * sizeof(double));
    *Cz = double_malloc_check(l_t * sizeof(double));

    for (int i = 0; i < l_t; i++)
    {
        (*Cx)[i] = seval(Lp, p_of_t[i], p, x, c1x, c2x, c3x, last);
        (*Cy)[i] = seval(Lp, p_of_t[i], p, y, c1y, c2y, c3y, last);
        (*Cz)[i] = seval(Lp, p_of_t[i], p, z, c1z, c2z, c3z, last);
    }
    free(p_of_t);

    free(c1x);
    free(c2x);
    free(c3x);
    free(c1y);
    free(c2y);
    free(c3y);
    free(c1z);
    free(c2z);
    free(c3z);

    return l_t;
}