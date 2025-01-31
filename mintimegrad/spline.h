#ifndef SPLINE_H
#define SPLINE_H

/* 
 * spline.h - Header for cubic spline interpolation functions.
 * 
 * Functions:
 *  - spline(): Compute the cubic spline coefficients.
 *  - seval(): Evaluate the cubic spline at a given point.
 *  - deriv(): Compute the first derivative of the spline at a point.
 *  - deriv2(): Compute the second derivative of the spline at a point.
 *  - sinteg(): Compute the integral of the spline up to a given point.
 */

/* Function prototypes */

/**
 * Compute the cubic spline coefficients.
 * 
 * @param n      Number of data points.
 * @param end1   Boundary condition for the first endpoint.
 * @param end2   Boundary condition for the last endpoint.
 * @param slope1 Slope at the first endpoint.
 * @param slope2 Slope at the last endpoint.
 * @param x      Array of x-coordinates of data points.
 * @param y      Array of y-coordinates of data points.
 * @param b      Output array for spline coefficients.
 * @param c      Output array for spline coefficients.
 * @param d      Output array for spline coefficients.
 * @param iflag  Status flag (0 = success, nonzero = error).
 * 
 * @return Status code (0 on success, nonzero on failure).
 */
int spline(int n, int end1, int end2, double slope1, double slope2,
           double x[], double y[], double b[], double c[], double d[], int *iflag);

/**
 * Evaluate the cubic spline at a given point.
 * 
 * @param n     Number of data points.
 * @param u     Point at which to evaluate the spline.
 * @param x     Array of x-coordinates of data points.
 * @param y     Array of y-coordinates of data points.
 * @param b     Array of spline coefficients.
 * @param c     Array of spline coefficients.
 * @param d     Array of spline coefficients.
 * @param last  Segment in which the last evaluation occurred.
 * 
 * @return Value of the spline at u.
 */
double seval(int n, double u, double x[], double y[], double b[], double c[], double d[], int *last);

/**
 * Evaluate the first derivative of the cubic spline at a given point.
 * 
 * @param n     Number of data points.
 * @param u     Point at which to evaluate the derivative.
 * @param x     Array of x-coordinates of data points.
 * @param b     Array of spline coefficients.
 * @param c     Array of spline coefficients.
 * @param d     Array of spline coefficients.
 * @param last  Segment in which the last evaluation occurred.
 * 
 * @return Value of the first derivative at u.
 */
double deriv(int n, double u, double x[], double b[], double c[], double d[], int *last);

/**
 * Evaluate the second derivative of the cubic spline at a given point.
 * 
 * @param n     Number of data points.
 * @param u     Point at which to evaluate the second derivative.
 * @param x     Array of x-coordinates of data points.
 * @param b     Array of spline coefficients.
 * @param c     Array of spline coefficients.
 * @param d     Array of spline coefficients.
 * @param last  Segment in which the last evaluation occurred.
 * 
 * @return Value of the second derivative at u.
 */
double deriv2(int n, double u, double x[], double b[], double c[], double d[], int *last);

/**
 * Compute the integral of the spline up to a given point.
 * 
 * @param n     Number of data points.
 * @param u     Point up to which to compute the integral.
 * @param x     Array of x-coordinates of data points.
 * @param y     Array of y-coordinates of data points.
 * @param b     Array of spline coefficients.
 * @param c     Array of spline coefficients.
 * @param d     Array of spline coefficients.
 * @param last  Segment in which the last evaluation occurred.
 * 
 * @return Integral of the spline up to u.
 */
double sinteg(int n, double u, double x[], double y[], double b[], double c[], double d[], int *last);

#endif /* SPLINE_H */