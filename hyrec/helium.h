/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         helium.h: all functions related to helium recombination                               */
/*                                                                                               */
/*         Version: January 2011  (first released November 2010)                                 */
/*************************************************************************************************/ 

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#define M_PI_2 1.57079632679489661923132169164
#endif

double rec_sahaHeII(double nH0, double Tr0, double fHe, double z, double *xHeIII);
double rec_sahaHeI(double nH0, double Tr0, double fHe, double z);
double rec_helium_dxedt(double xe, double nH0, double Tr0, double fHe, double H, double z);
double rec_saha_xe_H(double nH0, double T0, double z);
double xe_PostSahaHe(double nH0, double Tr0, double fHe, double H, double z, double *Delta_xe);
