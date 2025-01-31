#ifndef MTG_H
#define MTG_H

/**
 * @brief Generates a rotationally invariant k-space trajectory.
 * 
 * @param Cx Pointer to store the dynamically allocated x-coordinates of the trajectory.
 * @param Cy Pointer to store the dynamically allocated y-coordinates of the trajectory.
 * @param Cz Pointer to store the dynamically allocated z-coordinates of the trajectory.
 * @param x Array containing the x-coordinates of the input curve.
 * @param y Array containing the y-coordinates of the input curve.
 * @param z Array containing the z-coordinates of the input curve.
 * @param Lp Length of the input arrays x, y, and z.
 * @param g0 Initial gradient magnitude.
 * @param gfin Final gradient magnitude.
 * @param gmax Maximum allowable gradient magnitude.
 * @param smax Maximum allowable slew rate.
 * @param T Time interval between samples.
 * @param ds Arc-length step size for ODE solving.
 * 
 * @return The number of points in the generated trajectory.
 */
int riv_c(double **Cx, double **Cy, double **Cz, double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax,
          double smax, double T, double ds);

/**
 * @brief Generates a rotationally variant k-space trajectory.
 * 
 * @param Cx Pointer to store the dynamically allocated x-coordinates of the trajectory.
 * @param Cy Pointer to store the dynamically allocated y-coordinates of the trajectory.
 * @param Cz Pointer to store the dynamically allocated z-coordinates of the trajectory.
 * @param x Array containing the x-coordinates of the input curve.
 * @param y Array containing the y-coordinates of the input curve.
 * @param z Array containing the z-coordinates of the input curve.
 * @param Lp Length of the input arrays x, y, and z.
 * @param g0 Initial gradient magnitude.
 * @param gfin Final gradient magnitude.
 * @param gmax Maximum allowable gradient magnitude.
 * @param smax Maximum allowable slew rate.
 * @param T Time interval between samples.
 * @param ds Arc-length step size for ODE solving.
 * 
 * @return The number of points in the generated trajectory.
 */
int rv_c(double **Cx, double **Cy, double **Cz, double *x, double *y, double *z, int Lp, double g0, double gfin, double gmax,
         double smax, double T, double ds);

#endif /* MTG_H */