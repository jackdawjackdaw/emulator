//! A C implementation of the BFGS algorithm as sketched out in mathematica

/** 
 * This should provide an additional maximisation method for libEmu
 */

/**
 * function pointers
 * float fn (float a, float b)
 * 
 * void foo( float (*ptr)(float, float)){
 * result = ptr(a, b);
 * }
 *
 * foo(&fn);
 */
