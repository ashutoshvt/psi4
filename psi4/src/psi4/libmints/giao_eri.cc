/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/physconst.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <stdexcept>
#include <string>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
;
using namespace psi;

namespace {

/**
     * @brief Fills the primitive data structure used by libint/libderiv with information from the ShellPairs
     * @param PrimQuartet The structure to hold the data.
     * @param fjt Object used to compute the fundamental integrals.
     * @param p12 ShellPair data structure for the left
     * @param p34 ShellPair data structure for the right
     * @param am Total angular momentum of this quartet
     * @param nprim1 Number of primitives on center 1
     * @param nprim2 Number of primitives on center 2
     * @param nprim3 Number of primitives on center 3
     * @param nprim4 Number of primitives on center 4
     * @param sh1eqsh2 Is the shell on center 1 identical to that on center 2?
     * @param sh3eqsh4 Is the shell on center 3 identical to that on center 4?
     * @param deriv_lvl Derivitive level of the integral
     * @return The total number of primitive combinations found. This is passed to libint/libderiv.
     */
static size_t fill_primitive_data(prim_data *PrimQuartet, Fjt *fjt,
                                  const ShellPair *p12, const ShellPair *p34,
                                  int am,
                                  int nprim1, int nprim2, int nprim3, int nprim4,
                                  bool sh1eqsh2, bool sh3eqsh4, int deriv_lvl)
{
    double zeta, eta, ooze, rho, poz, coef1, PQx, PQy, PQz, PQ2, Wx, Wy, Wz, o12, o34, T, *F;
    double a1, a2, a3, a4;
    int p1, p2, p3, p4, i;
    size_t nprim = 0L;
    double *pai = p12->ai;
    double *pgamma12 = p12->gamma[0];
    double *poverlap12 = p12->overlap[0];
    for (p1 = 0; p1 < nprim1; ++p1) {
        a1 = *pai;
        ++pai;
        double *paj = p12->aj;
        for (p2 = 0; p2 < nprim2; ++p2) {
            a2 = *paj;
            zeta = *pgamma12;
            o12 = *poverlap12;
            ++paj;
            ++pgamma12;
            ++poverlap12;
            double PAx = p12->PA[p1][p2][0];
            double PAy = p12->PA[p1][p2][1];
            double PAz = p12->PA[p1][p2][2];
            double PBx = p12->PB[p1][p2][0];
            double PBy = p12->PB[p1][p2][1];
            double PBz = p12->PB[p1][p2][2];
            double PABx = p12->P[p1][p2][0];
            double PABy = p12->P[p1][p2][1];
            double PABz = p12->P[p1][p2][2];

            double *pak = p34->ai;
            double *pgamma34 = p34->gamma[0];
            double *poverlap34 = p34->overlap[0];
            for (p3 = 0; p3 < nprim3; ++p3) {
                a3 = *pak;
                ++pak;
                double *pal = p34->aj;
                for (p4 = 0; p4 < nprim4; ++p4) {
                    a4 = *pal;
                    eta = *pgamma34;
                    o34 = *poverlap34;
                    ++pal;
                    ++pgamma34;
                    ++poverlap34;

                    double PCx = p34->PA[p3][p4][0];
                    double PCy = p34->PA[p3][p4][1];
                    double PCz = p34->PA[p3][p4][2];
                    double PDx = p34->PB[p3][p4][0];
                    double PDy = p34->PB[p3][p4][1];
                    double PDz = p34->PB[p3][p4][2];
                    double PCDx = p34->P[p3][p4][0];
                    double PCDy = p34->P[p3][p4][1];
                    double PCDz = p34->P[p3][p4][2];

                    ooze = 1.0 / (zeta + eta);
                    poz = eta * ooze;
                    rho = zeta * poz;
                    coef1 = 2.0 * sqrt(rho * M_1_PI) * o12 * o34;

                    PrimQuartet[nprim].poz = poz;
                    PrimQuartet[nprim].oo2zn = 0.5 * ooze;
                    PrimQuartet[nprim].pon = zeta * ooze;
                    PrimQuartet[nprim].oo2z = 0.5 / zeta;
                    PrimQuartet[nprim].oo2n = 0.5 / eta;
                    PrimQuartet[nprim].twozeta_a = 2.0 * a1;
                    PrimQuartet[nprim].twozeta_b = 2.0 * a2;
                    PrimQuartet[nprim].twozeta_c = 2.0 * a3;
                    PrimQuartet[nprim].twozeta_d = 2.0 * a4;

                    PQx = PABx - PCDx;
                    PQy = PABy - PCDy;
                    PQz = PABz - PCDz;
                    PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

                    Wx = (PABx * zeta + PCDx * eta) * ooze;
                    Wy = (PABy * zeta + PCDy * eta) * ooze;
                    Wz = (PABz * zeta + PCDz * eta) * ooze;

                    // PA
                    PrimQuartet[nprim].U[0][0] = PAx;
                    PrimQuartet[nprim].U[0][1] = PAy;
                    PrimQuartet[nprim].U[0][2] = PAz;
                    // PB
                    PrimQuartet[nprim].U[1][0] = PBx;
                    PrimQuartet[nprim].U[1][1] = PBy;
                    PrimQuartet[nprim].U[1][2] = PBz;
                    // QC
                    PrimQuartet[nprim].U[2][0] = PCx;
                    PrimQuartet[nprim].U[2][1] = PCy;
                    PrimQuartet[nprim].U[2][2] = PCz;
                    // QD
                    PrimQuartet[nprim].U[3][0] = PDx;
                    PrimQuartet[nprim].U[3][1] = PDy;
                    PrimQuartet[nprim].U[3][2] = PDz;
                    // WP
                    PrimQuartet[nprim].U[4][0] = Wx - PABx;
                    PrimQuartet[nprim].U[4][1] = Wy - PABy;
                    PrimQuartet[nprim].U[4][2] = Wz - PABz;
                    // WQ
                    PrimQuartet[nprim].U[5][0] = Wx - PCDx;
                    PrimQuartet[nprim].U[5][1] = Wy - PCDy;
                    PrimQuartet[nprim].U[5][2] = Wz - PCDz;

                    T = rho * PQ2;
                    fjt->set_rho(rho);
                    F = fjt->values(am + deriv_lvl, T);

                    for (i = 0; i <= am + deriv_lvl; ++i)
                        PrimQuartet[nprim].F[i] = F[i] * coef1;

                    nprim++;
                }
            }
        }
    }
    return nprim;
}

} // end namespace

size_t GiaoERI::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
#ifdef MINTS_TIMER
    timer_on("setup");
#endif

    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();
    int am = am1 + am2 + am3 + am4; // total am
    int nprim1;
    int nprim2;
    int nprim3;
    int nprim4;
    double A[3], B[3], C[3], D[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];
    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];
    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libint_.AB[0] = A[0] - B[0];
    libint_.AB[1] = A[1] - B[1];
    libint_.AB[2] = A[2] - B[2];
    libint_.CD[0] = C[0] - D[0];
    libint_.CD[1] = C[1] - D[1];
    libint_.CD[2] = C[2] - D[2];

#ifdef MINTS_TIMER
    timer_off("setup");
#endif

#ifdef MINTS_TIMER
    timer_on("Primitive setup");
#endif

    // Prepare all the data needed by libint
    size_t nprim = 0;
    nprim1 = s1.nprimitive();
    nprim2 = s2.nprimitive();
    nprim3 = s3.nprimitive();
    nprim4 = s4.nprimitive();

    // If we can, use the precomputed values found in ShellPair.
    if (use_shell_pairs_) {
        ShellPair *p12, *p34;
        // 1234 -> 1234 no change
        p12 = &(pairs12_[sh1][sh2]);
        p34 = &(pairs34_[sh3][sh4]);

        nprim = fill_primitive_data(libint_.PrimQuartet, fjt_, p12, p34, am, nprim1, nprim2, nprim3, nprim4, sh1 == sh2, sh3 == sh4, 0);
    } else {
        const double *a1s = s1.exps();
        const double *a2s = s2.exps();
        const double *a3s = s3.exps();
        const double *a4s = s4.exps();
        const double *c1s = s1.coefs();
        const double *c2s = s2.coefs();
        const double *c3s = s3.coefs();
        const double *c4s = s4.coefs();

        // Old version - without ShellPair - STILL USED BY RI CODES
        for (int p1 = 0; p1 < nprim1; ++p1) {
            double a1 = a1s[p1];
            double c1 = c1s[p1];
            for (int p2 = 0; p2 < nprim2; ++p2) {
                double a2 = a2s[p2];
                double c2 = c2s[p2];
                double zeta = a1 + a2;
                double ooz = 1.0 / zeta;
                double oo2z = 1.0 / (2.0 * zeta);

                double PA[3], PB[3];
                double P[3];

                P[0] = (a1 * A[0] + a2 * B[0]) * ooz;
                P[1] = (a1 * A[1] + a2 * B[1]) * ooz;
                P[2] = (a1 * A[2] + a2 * B[2]) * ooz;
                PA[0] = P[0] - A[0];
                PA[1] = P[1] - A[1];
                PA[2] = P[2] - A[2];
                PB[0] = P[0] - B[0];
                PB[1] = P[1] - B[1];
                PB[2] = P[2] - B[2];

                double Sab = pow(M_PI * ooz, 3.0 / 2.0) * exp(-a1 * a2 * ooz * AB2) * c1 * c2;

                for (int p3 = 0; p3 < nprim3; ++p3) {
                    double a3 = a3s[p3];
                    double c3 = c3s[p3];
                    for (int p4 = 0; p4 < nprim4; ++p4) {
                        double a4 = a4s[p4];
                        double c4 = c4s[p4];
                        double nu = a3 + a4;
                        double oon = 1.0 / nu;
                        double oo2n = 1.0 / (2.0 * nu);
                        double oo2zn = 1.0 / (2.0 * (zeta + nu));
                        double rho = (zeta * nu) / (zeta + nu);
                        double oo2rho = 1.0 / (2.0 * rho);

                        double QC[3], QD[3], WP[3], WQ[3], PQ[3];
                        double Q[3], W[3], a3C[3], a4D[3];

                        a3C[0] = a3 * C[0];
                        a3C[1] = a3 * C[1];
                        a3C[2] = a3 * C[2];

                        a4D[0] = a4 * D[0];
                        a4D[1] = a4 * D[1];
                        a4D[2] = a4 * D[2];

                        Q[0] = (a3C[0] + a4D[0]) * oon;
                        Q[1] = (a3C[1] + a4D[1]) * oon;
                        Q[2] = (a3C[2] + a4D[2]) * oon;

                        QC[0] = Q[0] - C[0];
                        QC[1] = Q[1] - C[1];
                        QC[2] = Q[2] - C[2];
                        QD[0] = Q[0] - D[0];
                        QD[1] = Q[1] - D[1];
                        QD[2] = Q[2] - D[2];
                        PQ[0] = P[0] - Q[0];
                        PQ[1] = P[1] - Q[1];
                        PQ[2] = P[2] - Q[2];

                        double PQ2 = 0.0;
                        PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                        PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                        PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                        W[0] = (zeta * P[0] + nu * Q[0]) / (zeta + nu);
                        W[1] = (zeta * P[1] + nu * Q[1]) / (zeta + nu);
                        W[2] = (zeta * P[2] + nu * Q[2]) / (zeta + nu);
                        WP[0] = W[0] - P[0];
                        WP[1] = W[1] - P[1];
                        WP[2] = W[2] - P[2];
                        WQ[0] = W[0] - Q[0];
                        WQ[1] = W[1] - Q[1];
                        WQ[2] = W[2] - Q[2];

                        for (int i = 0; i < 3; ++i) {
                            libint_.PrimQuartet[nprim].U[0][i] = PA[i];
                            libint_.PrimQuartet[nprim].U[2][i] = QC[i];
                            libint_.PrimQuartet[nprim].U[4][i] = WP[i];
                            libint_.PrimQuartet[nprim].U[5][i] = WQ[i];
                        }
                        libint_.PrimQuartet[nprim].oo2z = oo2z;
                        libint_.PrimQuartet[nprim].oo2n = oo2n;
                        libint_.PrimQuartet[nprim].oo2zn = oo2zn;
                        libint_.PrimQuartet[nprim].poz = rho * ooz;
                        libint_.PrimQuartet[nprim].pon = rho * oon;
                        libint_.PrimQuartet[nprim].oo2p = oo2rho;

                        double T = rho * PQ2;
                        fjt_->set_rho(rho);
                        double *F = fjt_->values(am, T);

                        // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                        double Scd = pow(M_PI * oon, 3.0 / 2.0) * exp(-a3 * a4 * oon * CD2) * c3 * c4;
                        double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;
                        for (int i = 0; i <= am; ++i) {
                            libint_.PrimQuartet[nprim].F[i] = F[i] * val;
                        }
                        nprim++;
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Primitive setup");
#endif

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

#ifdef MINTS_TIMER
    timer_on("libint overhead");
#endif

    // Compute the integral
    if (am) {
        double *target_ints;

        target_ints = build_eri[am1][am2][am3][am4](&libint_, nprim);

        memcpy(source_, target_ints, sizeof(double) * size);
    } else {
        // Handle (ss|ss)
        double temp = 0.0;
        for (size_t i = 0; i < nprim; ++i)
            temp += (double) libint_.PrimQuartet[i].F[0];
        source_[0] = temp;
//        outfile->Printf( "s-functions = %8.5f\n", temp);
    }

#ifdef MINTS_TIMER
    timer_off("libint overhead");
#endif

    // The following two functions time themselves.

    // Normalize the integrals for angular momentum
    //normalize_am(s1, s2, s3, s4);

    // Transform the integrals into pure angular momentum
    if (!force_cartesian_)
        pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
    return size;
}
