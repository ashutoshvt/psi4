/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <psi4/libciomr/libciomr.h>
#include <psi4/libdpd/dpd.h>
#include <psi4/libqt/qt.h>
#include "Params.h"
#define EXTERN
//#include "globals.h"
#include "MOInfo.h"
#include "Local.h"

#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/vector.h>
#include <psi4/libpsio/psio.h>
#include "New_Local.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/// Constructor
New_Local::New_Local(MOInfo & moinfo_, Local & local_) {
    /* Set up basic info */
    nocc_    = moinfo_.clsdpi[0];
    nvir_    = moinfo_.nmo - nocc_ - moinfo_.frdocc[0] - moinfo_.fruocc[0];
    npairs_  = nocc_ * nocc_;
    local_type_ = "PNO"; 
    pno_cut_    = local_.pno_cutoff;
    pair_cut_   = local_.pair_energy_cutoff;
    singles_cut_ = local_.singles_cutoff;
    //pno_cut_ = options.get_double("PNO_CUTOFF");

    /* Grab Local Occ. Orb. Energies */
    get_occ_epsilons();

    /* Initialize PNO */
    init_pno();
}

New_Local::New_Local(MOInfo & moinfo_, Local & local_, bool init) {
    /* Set up basic info */
    //  Getting nocc_ and vir_ from local struct doesn't work, b/c
    //  b/c they get set in local_init(), which doesn't get called for
    //  the PNO approach.  However, local.method and local.pno_cutoff
    //  both get set in get_params().
    //nocc_    = local.nocc;
    //nvir_    = local.nvir;
    nocc_    = moinfo_.clsdpi[0];
    nvir_    = moinfo_.nmo - nocc_ - moinfo_.frdocc[0] - moinfo_.fruocc[0];
    npairs_  = nocc_ * nocc_;
    pno_cut_ = local_.pno_cutoff;
    pair_cut_ = local_.pair_energy_cutoff;
    singles_cut_ = local_.singles_cutoff;
    local_type_ = "PNO"; 
    //pno_cut_ = options.get_double("PNO_CUTOFF");

    /* Grab Local Occ. Orb. Energies */
    get_occ_epsilons();

    /* Initialize PNO */
    if(init) init_pno();
}

New_Local::New_Local(MOInfo & moinfo_, Local & local_, std::string local_type, bool init) {
    /* Set up basic info */
    nocc_    = moinfo_.clsdpi[0];
    nvir_    = moinfo_.nmo - nocc_  - moinfo_.frdocc[0] - moinfo_.fruocc[0];
    npairs_  = nocc_ * nocc_;
    pair_cut_ = local_.pair_energy_cutoff;
    local_type_  = local_type;
    singles_cut_ = local_.singles_cutoff;
    /*
    fprintf(outfile, "\tDimension Check:\n");
    fprintf(outfile, "\tNum. Active Occ. = %d\n", nocc_);
    fprintf(outfile, "\tNum. Frozen Occ. = %d\n", moinfo_.frdocc[0]);
    fprintf(outfile, "\tNum. Active Vir. = %d\n", nvir_);
    fprintf(outfile, "\tNum. Frozen Vir. = %d\n", moinfo_.fruocc[0]);
    fprintf(outfile, "\tNum. NMO    = %d\n", moinfo_.nmo);
    fprintf(outfile, "\tNum. OCCPI  = %d\n", moinfo_.occpi[0]);
    fprintf(outfile, "\tNum. CLSDPI = %d\n", moinfo_.clsdpi[0]);
    fflush(outfile);
    */

    /* Grab Local Occ. Orb. Energies */
    get_occ_epsilons();

    /* Initialize PNO */
    if(local_type_ == "PNO") {
        //pno_cut_ = options.get_double("PNO_CUTOFF");
        pno_cut_ = local_.pno_cutoff;
        if(init) init_pno();
    }
    if(local_type_ == "OSV") {
        //osv_cut_ = options.get_double("OSV_CUTOFF");
        osv_cut_ = local_.osv_cutoff;
        osv_type_ = local_.osv_type;
        if(init) init_osv();
    }
}

void New_Local::init_pno(void) {
    //fprintf(outfile, "\tINIT PNO\n");
    //fflush(outfile);
    dpdbuf4 T2; // T
    dpdbuf4 N;  // T~
    dpdbuf4 D;  // Energy denominator
  
    // DENSITIES //
    //    1) T   = tIjAb  = <ij|ab>/(fii+fjj-faa-fbb)
    //    2) T~  = tIjAb~ = (2<ij|ab> - <ij|ba>)/(fii+fjj-faa-fbb)
    //    3) Nrm = 1/(1 + i==j)
    //    4) Dij = 2.0 * [T(T~)* + T*(T~)] * Nrm
    
    /* Put T2, T2~ into vectors of IJ-matrices */
    std::vector<SharedMatrix> Tij; 
    std::vector<SharedMatrix> Ttij;
    std::vector<SharedMatrix> dints;
    SharedMatrix temp(new Matrix(nvir_, nvir_));
  
    /* Initialize T Amps */
    // This happens in mp2_energy().

    /* Initialize T~ */
    global_dpd_->buf4_init(&N, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    get_matvec(&N, &dints); // Grab dints here for MP2 pair energy computation.
    global_dpd_->buf4_copy(&N, PSIF_CC_TMP0, "tIjAb ~");
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_init(&N, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb ~");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&D, &N);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&N);
  
    /* "Vectorize" T2, T2~ */
    global_dpd_->buf4_init(&N, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb ~");
    get_matvec(&N, &Ttij);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
/*
// New start
    global_dpd_->buf4_init(&T2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "tIjAb");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
// New end
*/
    get_matvec(&T2, &Tij);
    global_dpd_->buf4_close(&T2);
    // Print checking
    /*
    fprintf(outfile, "\t***** T Matrices *****\n");
    for(int ij=0; ij<Tij.size(); ++ij)
        Tij[ij]->print();
    fflush(outfile);
    fprintf(outfile, "\t***** T~ Matrices *****\n");
    for(int ij=0; ij<Ttij.size(); ++ij)
        Ttij[ij]->print();
    fflush(outfile);
    */

    /* Compute Pair Energies, Check for Weak Pairs */
    //if(weak_pairs_) {
        lmp2_pair_energy(Tij,dints);
    //}
  
    /* Form Pair Densities */
    std::vector<SharedMatrix> Dij;
    temp->zero();
    for(int ij=0; ij < npairs_; ++ij) {
        temp->zero();
        int i = ij/nocc_;
        int j = ij%nocc_;
        temp->gemm(0, 1, 1, Tij[ij], Ttij[ij], 0);   
        temp->gemm(1, 0, 1, Tij[ij], Ttij[ij], 1);   
        temp->scale(2.0 * 1/(1+(i==j)));
        Dij.push_back(temp->clone());
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Dij Matrics *****\n");
    for(int ij=0; ij < npairs_; ++ij) {
        int i = ij/nocc_;
        int j = ij%nocc_;
        fprintf(outfile, "ij = %2d, i = %2d, j = %2d", ij, i, j);
        Dij[ij]->print();
    }
    fflush(outfile);
    */
  
    /* Get PNO Transforms, Occupation numbers (aka diagonalize each pair density) */
    std::vector<SharedMatrix> Q_full(npairs_);
    occ_num_.resize(npairs_);
    for(int ij=0; ij < npairs_; ++ij) {
        occ_num_[ij] = SharedVector(new Vector("PNO Occupation #'s", nvir_));
        Q_full[ij]   = SharedMatrix(new Matrix("Full PNO Transform", nvir_, nvir_));
    }
  
    for(int ij=0; ij < npairs_; ++ij) {
        Dij[ij]->diagonalize(Q_full[ij], occ_num_[ij], descending);
    }
    /*
    for(int ij=0; ij < npairs_; ++ij) {
        fprintf(outfile, "Pair: %d\n", ij);
        for(int a=0; a < nvir_; ++a) {
            fprintf(outfile, "%20.12lf\n", occ_num_[ij]->get(a));
        }
    }
    */
  
    /* Truncate the Transforms */
    if(singles_cut_) {
        init_sep_singles(Q_full);
    }

    // List of PNO dimensions.
    double n_pno;
    int survivor; 
    dim_ = SharedVector(new Vector(npairs_));
    for(int ij=0; ij < npairs_; ++ij) {
        survivor = 0;
        for(int a=0; a < nvir_; ++a) {
            n_pno = fabs(occ_num_[ij]->get(a));
            if(n_pno >= pno_cut_) {
                ++survivor;
            }
        }
        dim_->set(ij,survivor);
    }

    /* Print some stats */
    //fprintf(outfile, "\nNumber of Pairs: %d\n", npairs_);
    double avg_pno = 0.0;
    for(int ij=0; ij < npairs_; ++ij)
        if(!weak_pairs_[ij]) avg_pno += dim_->get(ij);
    avg_pno /= npairs_;
    //fprintf(outfile, "\tAverage # PNO's: %10.3lf\n", avg_pno);
    //fprintf(outfile, "\nPNO Dimensions:\n");
    //fflush(outfile);
    dim_->print("outfile");
    //fflush(outfile);

    /* T1,T2 Length for direct comparison to PAO approach */
    int t1_length  = 0;
    int t2_length  = 0;
    for(int ij=0; ij < npairs_; ++ij) {
            int i = ij/(nocc_);
            int j = ij%(nocc_);
            if(i==j) t1_length += (int) dim_->get(ij);
            if(!weak_pairs_[ij])
                t2_length += (int) (dim_->get(ij) * dim_->get(ij));
    }
    if(singles_cut_) {
        t1_length = 0;
        for(int i=0; i < nocc_; ++i)
            t1_length += (int) t1dim_->get(i);
    }
    //fprintf(outfile, "\n\tT1 Length = %d (local), %d (canonical)\n",
    //        t1_length, nocc_*nvir_);
    //fprintf(outfile, "\tT2 Length = %d (local), %d (canonical)\n\n",
    //        t2_length, nocc_*nocc_*nvir_*nvir_);

    // Make list of surviving pairs.
    get_pair_list(dim_);
    int spairs = pair_list_.size();
    /*
    fprintf(outfile, "\tSurviving Pairs:\n");
    //fflush(outfile);
    for(int x=0; x < spairs; ++x) {
        fprintf(outfile, "\tPair: %d\n", pair_list_[x] + 1);
        fflush(outfile);
    }
    */

    // Form truncated transforms.
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        int pno = dim_->get(pair_idx);
        SharedMatrix qtemp(new Matrix(nvir_, pno));
        for(int a=0; a < nvir_; ++a)
            for(int aij=0; aij < pno; ++aij)
                qtemp->set(a,aij, Q_full[pair_idx]->get(a,aij));
        Q_.push_back(qtemp->clone());
        qtemp->zero();
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Truncated Q Matrices ****\n");
    for(int ij=0; ij < spairs; ++ij)
        Q_[ij]->print();
    fflush(outfile);
    */
  
    /* Get semicanonical transforms from virtual fock block */
     get_semicanonical_transform();
  
    /* Write the PNO dimensions and transforms to CC_INFO */
    psio_address next;
    psio_write_entry(PSIF_CC_INFO, "Local Dimensions", (char *) dim_->pointer(), npairs_ * sizeof(double));
    next = PSIO_ZERO;
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        psio_write(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) Q_[ij]->pointer()[0],
            nvir_ * dim_->get(pair_idx) * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        psio_write(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) L_[ij]->pointer()[0],
            dim_->get(pair_idx) * dim_->get(pair_idx) * sizeof(double), next, &next);
    }

    /* Now, can I read them back in?? */
    //  The test below worked!!  Setting "next = PSIO_ZERO" started the psio_read over, so I 
    //  could have done it like seven times, and would still get L[0] rather than L[6].
    //  Consecutive reads like below, though, grab the appropriate sets of data.  Yay!
    /*
    SharedMatrix Ltest(new Matrix("Ltest of PSIO",nvir_,nvir_));
    next = PSIO_ZERO;
    psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) Ltest->pointer()[0],
        nvir_ * nvir_ * sizeof(double), next, &next);
    psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) Ltest->pointer()[0],
        nvir_ * nvir_ * sizeof(double), next, &next);
    psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) Ltest->pointer()[0],
        nvir_ * nvir_ * sizeof(double), next, &next);
    L_[2]->print();
    Ltest->print();
    */
    
    free(weak_pairs_);
}


/*
 *  OSV Initializer - Same algorithm as for PNO's, but in this case
 *  we get the transforms from the diagonal (ii, jj) pairs only.
 *  The ij-pair transforms are the ii & jj transforms collated together.
 */
void New_Local::init_osv(void) {
    //fprintf(outfile, "\tINIT OSV\n");
    //fflush(outfile);
    dpdbuf4 T2; // T
    dpdbuf4 N;  // T~
    dpdbuf4 D;  // Energy denominator
  
    // DENSITIES //
    //    1) T   = tIjAb  = <ij|ab>/(fii+fjj-faa-fbb)
    //    2) T~  = tIjAb~ = (2<ij|ab> - <ij|ba>)/(fii+fjj-faa-fbb)
    //    3) Nrm = 1/(1 + i==j)
    //    4) Dij = 2.0 * [T(T~)* + T*(T~)] * Nrm
    
    /* Initialize T Amps */
    // This happens in mp2_energy().

    /* Matrix Vectors */
    std::vector<SharedMatrix> Tij;
    std::vector<SharedMatrix> Ttij;
    std::vector<SharedMatrix> dints;
    SharedMatrix temp(new Matrix(nvir_, nvir_));
  
    /* Initialize T~ */
    global_dpd_->buf4_init(&N, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    get_matvec(&N, &dints); // Grab dints here for MP2 pair energy computation.
    global_dpd_->buf4_copy(&N, PSIF_CC_TMP0, "tIjAb ~");
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_init(&N, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb ~");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&D, &N);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&N);
  
    /* Put T2, T2~ into vectors of IJ-matrices */
    global_dpd_->buf4_init(&N, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb ~");
    get_matvec(&N, &Ttij);
    global_dpd_->buf4_close(&N);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
/*
// New start
    global_dpd_->buf4_init(&T2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "tIjAb");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
// New end
*/
    get_matvec(&T2, &Tij);
    global_dpd_->buf4_close(&T2);
    // Print checking
    /*
    fprintf(outfile, "\t***** T Matrices *****\n");
    for(int ij=0; ij<Tij.size(); ++ij)
        Tij[ij]->print();
    fflush(outfile);
    fprintf(outfile, "\t***** T~ Matrices *****\n");
    for(int ij=0; ij<Ttij.size(); ++ij)
        Ttij[ij]->print();
    fflush(outfile);
    */

    /* Compute MP2 pair energies, Check for Weak Pairs */
    //if(weak_pairs_) {
        lmp2_pair_energy(Tij,dints);
    //}
  
    /* Form Diagonal Pair Densities */
    std::vector<SharedMatrix> Dii;
    temp->zero();
    if(osv_type_ == "DENSITY") {
        for(int i=0; i < nocc_; ++i) {
            temp->zero();
            int ii = i*nocc_ + i;
            temp->gemm(0, 1, 1, Tij[ii], Ttij[ii], 0);   
            temp->gemm(1, 0, 1, Tij[ii], Ttij[ii], 1);   
            //temp->scale(2.0 * 1/(1+(i==j)));
            Dii.push_back(temp->clone());
        }
    }
    else if(osv_type_ == "AMPLITUDE") {
        for(int i=0; i < nocc_; ++i) {
            int ii = i*nocc_ + i;
            Dii.push_back(Tij[ii]->clone());
        }
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Dij Matrics *****\n");
    for(int i=0; i < nocc_; ++i) {
        int ii = i*nocc_ + i;
        fprintf(outfile, "ii = %2d, i = %2d\n", ij, i, j);
        Dii[i]->print();
    }
    fflush(outfile);
    */
  
    /* Get OSV Transforms, Occupation numbers (aka diagonalize each pair density) */
    std::vector<SharedMatrix> Q_full(nocc_);
    occ_num_.resize(nocc_);
    for(int ii=0; ii < nocc_; ++ii) {
        occ_num_[ii] = SharedVector(new Vector("OSV Occupation #'s", nvir_));
        Q_full[ii]   = SharedMatrix(new Matrix("Full OSV Transform", nvir_, nvir_));
    }
  
    for(int ii=0; ii < nocc_; ++ii) {
        Dii[ii]->diagonalize(Q_full[ii], occ_num_[ii], descending);
    }
    /*
    for(int ii=0; ii < nocc_; ++ii) {
        fprintf(outfile, "Pair: %d (%d orig.)", ii, ii*nocc_ + ii);
        for(int a=0; a < nvir_; ++a) {
            fprintf(outfile, "%20.12lf\n", occ_num_[ii]->get(a));
        }
    }
    */
  
    /* Truncate the Transforms */
    if(singles_cut_) {
        init_sep_singles(Q_full);
    }
    // List of OSV dimensions.
    double n_osv;
    int survivor; 
    //dim_ = SharedVector(new Vector(nocc_));
    SharedVector dim(new Vector(nocc_));
    for(int ii=0; ii < nocc_; ++ii) {
        survivor = 0;
        for(int a=0; a < nvir_; ++a) {
            n_osv = fabs(occ_num_[ii]->get(a));
            if(n_osv >= osv_cut_) {
                ++survivor;
            }
        }
        dim->set(ii,survivor);
    }

    /* Print some stats */
    //fprintf(outfile, "\nNumber of Diagonal Pairs: %d\n", nocc_);
    double avg_osv = 0.0;
    for(int ii=0; ii < nocc_; ++ii)
        avg_osv += dim->get(ii);
    avg_osv /= nocc_;
    //fprintf(outfile, "\tDiagonal Pair Average # OSV's: %10.3lf\n", avg_osv);
    /*
    fprintf(outfile, "\nOSV Dimensions:\n");
    fflush(outfile);
    dim->print(outfile);
    fflush(outfile);
    */

    // Make list of surviving pairs.
    get_pair_list(dim);
    int spairs = pair_list_.size();
    /*
    fprintf(outfile, "\tSurviving Pairs:\n");
    fflush(outfile);
    for(int x=0; x < spairs; ++x) {
        fprintf(outfile, "\tPair: %d\n", pair_list_[x] + 1);
        fflush(outfile);
    }
    */

    // Form truncated transforms.
    std::vector<SharedMatrix> Qii;
    for(int ii=0; ii < spairs; ++ii) {
        int pair_idx = pair_list_[ii];
        int osv = dim->get(pair_idx);
        SharedMatrix qtemp(new Matrix(nvir_, osv));
        for(int a=0; a < nvir_; ++a)
            for(int aii=0; aii < osv; ++aii)
                qtemp->set(a,aii, Q_full[pair_idx]->get(a,aii));
        Qii.push_back(qtemp->clone());
        qtemp->zero();
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Truncated Q Matrices ****\n");
    for(int ii=0; ii < spairs; ++ii)
        Qii[ii]->print();
    fflush(outfile);
    */

    /* Built IJ-pair Densities (horizonally concatenate Qii and Qjj) */
    std::vector<SharedMatrix> ijvec;
    for(int ii=0; ii < spairs; ++ii) {    
        for(int jj=0; jj < spairs; ++jj) {    
            SharedMatrix ijmat;
            ijvec.push_back(Qii[ii]);
            ijvec.push_back(Qii[jj]);
            ijmat = Matrix::horzcat(ijvec);
            Q_.push_back(ijmat->clone());
            ijvec.clear();
        }
    }

    // Keep Track of Original Indices - NOT REALLY NECESSARY
    /*
    SharedMatrix idx_lookup(new Matrix(nocc_, nocc_));
    for(int ii=0; ii < spairs; ++ii) {    
        int i = pair_list_[ii];
        for(int jj=0; jj < spairs; ++jj) {    
            int j = pair_list_[jj];
            idx_lookup->set(i,j, i*nocc_ + j);
        }
    }
    */

    /* Reset Q, Dimensions, Significant Pair List */
    // Print checking
    /*
    fprintf(outfile, "\tMade new Q vector...\n");
    fflush(outfile);
    for(int ij=0; ij < npairs_; ++ij) {
        Q_[ij]->print();
    }
    */
    dim_ = SharedVector(new Vector("OSV Pair Dimensions", npairs_));
    dim_->zero();
    int ijx;     // Index for Q matrices
    int ij_orig; // Index for original I,J,IJ
    for(int ii=0, ijx = 0; ii < spairs; ++ii) {
        int i = pair_list_[ii];
        for(int jj=0; jj < spairs; ++jj, ++ijx) {    
            int j = pair_list_[jj];
            //ij_orig = idx_lookup(i,j);
            ij_orig = i*nocc_ + j;
            dim_->set(ij_orig, Q_[ijx]->ncol());
        }
    }
    //dim_->print();

    get_pair_list(dim_);
    spairs = pair_list_.size();
    /// At this point, everything should be essentially the same as for PNO.
     // -> Q's are set up for all surviviing IJ-pairs.
     // -> The dim_ vector contains dimensions for all IJ-pairs.
  
    /* Get semicanonical transforms from virtual fock block */
    get_semicanonical_transform();

    /* Redo the Stats since removed linear dependencies */
    //fprintf(outfile, "\nNumber of Pairs: %d\n", npairs_);
    avg_osv = 0.0;
    for(int ij=0; ij < npairs_; ++ij)
        avg_osv += dim_->get(ij);
    avg_osv /= npairs_;
    //fprintf(outfile, "\tActual Average # OSV's: %10.3lf\n", avg_osv);
    //fflush(outfile);
    dim_->print("outfile");
    //fflush(outfile);

    /* T1,T2 Length for direct comparison to PAO approach */
    int t1_length  = 0;
    int t2_length  = 0;
    for(int ij=0; ij < npairs_; ++ij) {
        int i = ij/(nocc_);
        int j = ij%(nocc_);
        if(i==j) t1_length += (int) dim_->get(ij);
        if(!weak_pairs_[ij]) t2_length += (int) (dim_->get(ij) * dim_->get(ij));
    }
    if(singles_cut_) {
        t1_length = 0;
        for(int i=0; i < nocc_; ++i)
            t1_length += (int) t1dim_->get(i);
    }
    //fprintf(outfile, "\n\tT1 Length = %d (local), %d (canonical)\n",
    //        t1_length, nocc_*nvir_);
    //fprintf(outfile, "\tT2 Length = %d (local), %d (canonical)\n\n",
    //        t2_length, nocc_*nocc_*nvir_*nvir_);
  
    /* Write the PNO dimensions and transforms to CC_INFO */
    psio_address next;
    psio_write_entry(PSIF_CC_INFO, "Local Dimensions", (char *) dim_->pointer(), npairs_ * sizeof(double));
    next = PSIO_ZERO;
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        psio_write(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) Q_[ij]->pointer()[0],
            nvir_ * dim_->get(pair_idx) * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        psio_write(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) L_[ij]->pointer()[0],
            dim_->get(pair_idx) * dim_->get(pair_idx) * sizeof(double), next, &next);
    }
    /* Now, can I read them back in?? YES!! */
    /*
    SharedMatrix Qtest(new Matrix("Qtest of PSIO",nvir_,nvir_*2));
    fprintf(outfile, "Qtest:\n");
    fflush(outfile);
    next = PSIO_ZERO;
    psio_read(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) Qtest->pointer()[0],
        2 * nvir_ * nvir_ * sizeof(double), next, &next);
    psio_read(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) Qtest->pointer()[0],
        2 * nvir_ * nvir_ * sizeof(double), next, &next);
    Q_[1]->print();
    Qtest->print();
    fflush(outfile);
    */
}


/* Get set of PNO (or OSV) -> Semicanoncial Transforms */
   // 1) Transform Virtual Fock block to each PNO basis
   // 2) Diagonalize each PNO virtual Fock
   // BTW, assumes RHF for now.
void New_Local::get_semicanonical_transform(void) {
    int ij, a, b;
    int pno;
    dpdfile2 fab;
    SharedMatrix fvir(new Matrix("Virtual Fock Block",nvir_,nvir_));
    std::vector<SharedMatrix> fpno;  // Set of ij PNO-vir Fock blocks

    /* Grab virtual Fock block */
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);
    for(a=0; a < nvir_; ++a)
        for(b=0; b < nvir_; ++b)
            fvir->set(a,b, fab.matrix[0][a][b]);
    // Print checking
    /*
    fvir->print();
    global_dpd_->file2_mat_print(&fab, outfile);
    */
    global_dpd_->file2_close(&fab);

    /* Get Orthogonal Transform, Remove Linear Dependencies (OSV only) */
    if(local_type_ == "OSV") {
        int spairs = pair_list_.size();
        int pair_idx;
        int osv;

        /* Build Pair OSV Overlap */
        for(ij=0; ij < spairs; ++ij) {
            pair_idx = pair_list_[ij];
            osv = dim_->get(pair_idx); 
            SharedMatrix xtemp(new Matrix("OSV Overlap",osv,osv));
            SharedMatrix evec(new Matrix(osv,osv));
            SharedVector eval(new Vector(osv));
            xtemp->gemm(1, 0, 1, Q_[ij], Q_[ij], 0);
            X_.push_back(xtemp->clone());
        }
        /*
        fprintf(outfile, "\t FIRST TWO OSV PAIR OVERLAPS\n");
        X_[0]->print();
        X_[1]->print();
        */

        /* Diagonalize and Adjust Eigenvectors */
        std::vector<SharedVector> eval_list;
        for(ij=0; ij < spairs; ++ij) {
            pair_idx = pair_list_[ij];
            osv = dim_->get(pair_idx); 
            SharedMatrix evecs(new Matrix("X Matrix",osv,osv));
            SharedVector evals(new Vector("X Eigenvalues",osv));
            X_[ij]->diagonalize(evecs,evals,descending);

            eval_list.push_back(evals);
            X_[ij] = evecs->clone();
            /*
            evecs->print();
            evals->print();
            */
        }
        /*
        fprintf(outfile, "\t FIRST TWO NONADJUSTED X Matrices, eigenvectors\n");
        X_[0]->print();
        eval_list[0]->print();
        X_[1]->print();
        eval_list[1]->print();
        */

        /* Adjust Eigenvectors */
        int n_osv;
        for(ij=0; ij < spairs; ++ij) {
            pair_idx = pair_list_[ij];
            osv = dim_->get(pair_idx); 
            n_osv = 0;

            // Count Dimension!
            for(a=0; a < osv; ++a)
                if(eval_list[ij]->get(a) >= 1e-6) n_osv++;
            nr_dim_.push_back(n_osv);

            // Set up new Transform
            SharedMatrix xtemp(new Matrix("Reduced X Matrix", osv, n_osv));
            SharedVector xcol(new Vector(osv));
            for(a=0; a < n_osv; ++a) {
                xcol = X_[ij]->get_column(0, a);
                xcol->scale(1.0/sqrt(eval_list[ij]->get(a)));
                xtemp->set_column(0, a, xcol); 
            }
            X_[ij] = xtemp->clone();
        }
        /*
        fprintf(outfile, "\t FIRST TWO ADJUSTED X Matrices\n");
        X_[0]->print();
        X_[1]->print();
        */
        
        /* Combine Q and X, and track new dimensions */
        for(ij=0; ij < spairs; ++ij) {
            pair_idx = pair_list_[ij];
            osv = nr_dim_[ij];
            dim_->set(pair_idx, osv);
            SharedMatrix qprime(new Matrix(nvir_, osv));
            qprime->gemm(0, 0, 1, Q_[ij], X_[ij], 0);
            Q_[ij] = qprime->clone();
        }
        
    }

    /* Transform Fvir to PNO basis */
    int spairs = pair_list_.size();
    for(ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        pno = dim_->get(pair_idx);
        SharedMatrix atemp(new Matrix(nvir_,pno));
        SharedMatrix btemp(new Matrix(pno,pno));
        atemp->gemm(0, 0, 1, fvir, Q_[ij], 0);
        btemp->gemm(1, 0, 1, Q_[ij], atemp, 1);
        fpno.push_back(btemp->clone());
    }
    
    /* Diagonalize Fpno */
    // For non-truncated case, Evals should be MO-basis fock values (they are),
    // but we'll grab 'em here to be safe b/c with truncation (when Fpno is smaller
    // than Fvir) the energies may differ slightly.
    eps_vir_.resize(spairs);
    for(ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        pno = dim_->get(pair_idx);
        SharedMatrix evecs(new Matrix(pno,pno));
        eps_vir_[ij] = SharedVector(new Vector("Semicanonical Virtual Orbital Energies",pno));
        fpno[ij]->diagonalize(evecs, eps_vir_[ij], ascending);
        L_.push_back(evecs->clone());
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Ortho. Transforms *****\n");
    for(ij=0; ij<spairs; ++ij) 
        L_[ij]->print();
    */

    /* Write Virtual Orbital Energies */
    psio_address next;
    next = PSIO_ZERO;
    for(int ij=0; ij < spairs; ++ij) {
        int pair_idx = pair_list_[ij];
        pno = dim_->get(pair_idx);
        psio_write(PSIF_CC_INFO, "IJ - Virtual Orbital Energies", (char *) eps_vir_[ij]->pointer(),
            pno * sizeof(double), next, &next);
    }
}


/* Separate Singles Cutoff Setup */
void New_Local::init_sep_singles(std::vector<SharedMatrix> Q_full) {
    // List of PNO dimensions.
    double n_pno;
    int survivor; 
    t1dim_ = SharedVector(new Vector(nocc_));
    for(int i=0; i < nocc_; ++i) {
        int ii;
        if(occ_num_.size() == nocc_) ii = i;
        if(occ_num_.size() == npairs_) ii = i*nocc_ + i;

        survivor = 0;
        for(int a=0; a < nvir_; ++a) {
            n_pno = fabs(occ_num_[ii]->get(a));
            if(n_pno >= singles_cut_) {
                ++survivor;
            }
        }
        t1dim_->set(i,survivor);
    }

    /* Print some stats */
    double avg_pno = 0.0;
    for(int i=0; i < nocc_; ++i)
        avg_pno += t1dim_->get(i);
    avg_pno /= nocc_;
    //fprintf(outfile, "\tAverage # Virtuals (singles): %10.3lf\n", avg_pno);
    //fprintf(outfile, "\nLocal Singles Dimensions:\n");
    t1dim_->print("outfile");
    //fflush(outfile);

    // Form truncated transforms.
    for(int i=0; i < nocc_; ++i) {
        int ii;
        if(occ_num_.size() == nocc_) ii = i;
        if(occ_num_.size() == npairs_) ii = i*nocc_ + i;
        int pno = t1dim_->get(i);
        SharedMatrix qtemp(new Matrix(nvir_, pno));
        for(int a=0; a < nvir_; ++a)
            for(int aii=0; aii < pno; ++aii)
                qtemp->set(a,aii, Q_full[ii]->get(a,aii));
        Qs_.push_back(qtemp->clone());
        qtemp->zero();
    }
    // Print checking
    /*
    fprintf(outfile, "\t***** Truncated Q Matrices (Singles) ****\n");
    for(int i=0; i < nocc_; ++i)
        Qs_[i]->print();
    fflush(outfile);
    */
  
    /* Get semicanonical transforms from virtual fock block */
    sep_singles_semicanonical_transform();
    //fprintf(outfile, "\tFinished sep semi trans.\n");
    //fflush(outfile);
  
    /* Write the PNO dimensions and transforms to CC_INFO */
    psio_address next;
    psio_write_entry(PSIF_CC_INFO, "Local Singles Dimensions", (char *) t1dim_->pointer(), nocc_ * sizeof(double));
    next = PSIO_ZERO;
    for(int i=0; i < nocc_; ++i) {
        int pno = t1dim_->get(i);
        psio_write(PSIF_CC_INFO, "Local Singles Transformation Matrix (Q)", (char *) Qs_[i]->pointer()[0],
            nvir_ * pno * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int i=0; i < nocc_; ++i) {
        int pno = t1dim_->get(i);
        psio_write(PSIF_CC_INFO, "Semicanonical Singles Transformation Matrix (L)", (char *) Ls_[i]->pointer()[0],
            pno * pno * sizeof(double), next, &next);
    }
}


/* Sep. Singles Semicanonicalization */
void New_Local::sep_singles_semicanonical_transform(void) {
    int ii, i, a, b;
    int pno;
    dpdfile2 fab;
    SharedMatrix fvir(new Matrix("Virtual Fock Block",nvir_,nvir_));
    std::vector<SharedMatrix> fpno;  // Set of ij PNO-vir Fock blocks

    /* Grab virtual Fock block */
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);
    for(a=0; a < nvir_; ++a)
        for(b=0; b < nvir_; ++b)
            fvir->set(a,b, fab.matrix[0][a][b]);
    global_dpd_->file2_close(&fab);

    /* Get Orthogonal Transform, Remove Linear Dependencies (OSV only) */
//    if(local_type_ == "OSV") {
//        int osv;
//
//        /* Build Pair OSV Overlap */
//        for(i=0; i < nocc_; ++i) {
//            osv = t1dim_->get(i); 
//            SharedMatrix xtemp(new Matrix("OSV Overlap",osv,osv));
//            SharedMatrix evec(new Matrix(osv,osv));
//            SharedVector eval(new Vector(osv));
//            xtemp->gemm(1, 0, 1, Qs_[i], Qs_[i], 0);
//            X_.push_back(xtemp->clone());
//        }
//        /*
//        fprintf(outfile, "\t FIRST TWO OSV PAIR OVERLAPS\n");
//        X_[0]->print();
//        X_[1]->print();
//        */
//
//        /* Diagonalize and Adjust Eigenvectors */
//        std::vector<SharedVector> eval_list;
//        for(i=0; i < nocc_; ++i) {
//            osv = t1dim_->get(i); 
//            SharedMatrix evecs(new Matrix("X Matrix",osv,osv));
//            SharedVector evals(new Vector("X Eigenvalues",osv));
//            X_[i]->diagonalize(evecs,evals,descending);
//
//            eval_list.push_back(evals);
//            X_[i] = evecs->clone();
//            /*
//            evecs->print();
//            evals->print();
//            */
//        }
//
//        /* Adjust Eigenvectors */
//        int n_osv;
//        for(i=0; i < nocc_; ++i) {
//            osv = t1dim_->get(i); 
//            n_osv = 0;
//
//            // Count Dimension!
//            for(a=0; a < osv; ++a)
//                if(eval_list[i]->get(a) >= 1e-6) n_osv++;
//            nr_dim_.push_back(n_osv);
//
//            // Set up new Transform
//            SharedMatrix xtemp(new Matrix("Reduced X Matrix", osv, n_osv));
//            SharedVector xcol(new Vector(osv));
//            for(a=0; a < n_osv; ++a) {
//                xcol = X_[i]->get_column(0, a);
//                xcol->scale(1.0/sqrt(eval_list[i]->get(a)));
//                xtemp->set_column(0, a, xcol); 
//            }
//            X_[i] = xtemp->clone();
//        }
//        
//        /* Combine Q and X, and track new dimensions */
//        for(i=0; i < nocc_; ++i) {
//            osv = nr_dim_[i];
//            t1dim_->set(i, osv); 
//            SharedMatrix qprime(new Matrix(nvir_, osv));
//            qprime->gemm(0, 0, 1, Qs_[i], X_[i], 0);
//            Qs_[i] = qprime->clone();
//        }
//        X_.clear();
//    }   // end if OSV loop

    /* Transform Fvir to PNO basis */
    for(i=0; i < nocc_; ++i) {
        pno = t1dim_->get(i);
        SharedMatrix atemp(new Matrix(nvir_,pno));
        SharedMatrix btemp(new Matrix(pno,pno));
        atemp->gemm(0, 0, 1, fvir, Qs_[i], 0);
        btemp->gemm(1, 0, 1, Qs_[i], atemp, 0);
        fpno.push_back(btemp->clone());
    }
    
    /* Diagonalize Fpno */
    // For non-truncated case, Evals should be MO-basis fock values (they are),
    // but we'll grab 'em here to be safe b/c with truncation (when Fpno is smaller
    // than Fvir) the energies may differ slightly.
    sep_eps_vir_.resize(nocc_);
    for(i=0; i < nocc_; ++i) {
        pno = t1dim_->get(i);
        SharedMatrix evecs(new Matrix(pno,pno));
        sep_eps_vir_[i] = SharedVector(new Vector("Semicanonical Singles Virtual Orbital Energies",pno));
        fpno[i]->diagonalize(evecs, sep_eps_vir_[i], ascending);
        Ls_.push_back(evecs->clone());
    }

    /* Write Virtual Orbital Energies */
    psio_address next;
    next = PSIO_ZERO;
    for(int i=0; i < nocc_; ++i) {
        pno = t1dim_->get(i);
        psio_write(PSIF_CC_INFO, "II - Virtual Orbital Energies", (char *) sep_eps_vir_[i]->pointer(),
            pno * sizeof(double), next, &next);
    }
}

/* T2 Update (stepping through PNO (and Ortho.) basis) */
void New_Local::update_pno(dpdbuf4 * T2) {
    int ij, i, j, ab, a, b;
    int pno, pair_idx, spairs;
    
    // This function will use transforms stored on disk
    // instead of the member variables, Q_, L_.
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_vir;

    /* Grab PNO Dimensions */
    SharedVector dim(new Vector("Local Pair-Virtual Dimensions",npairs_));
    psio_read_entry(PSIF_CC_INFO, "Local Dimensions", (char *) dim->pointer(), npairs_ * sizeof(double));
    /*
    if(local_type_ == "PNO") fprintf(outfile, "\tGrabbed PNO dimensions...\n");
    if(local_type_ == "OSV") fprintf(outfile, "\tGrabbed OSV dimensions...\n");
    fflush(outfile);
    */

    /* Get Siginficant-Pair List */
    get_pair_list(dim);
    spairs = pair_list_.size();
    /*
    fprintf(outfile, "\tSurviving Pairs:\n");
    fflush(outfile);
    for(int x=0; x < spairs; ++x) {
        fprintf(outfile, "\tPair: %d\n", pair_list_[x] + 1);
        fflush(outfile);
    }
    fprintf(outfile, "\tNon-Surviving Pairs:\n");
    fflush(outfile);
    for(int x=0; x < zero_list_.size(); ++x) {
        fprintf(outfile, "\tPair: %d\n", zero_list_[x] + 1);
        fflush(outfile);
    }
    */

    /* Grab MO-Basis T2's */
    global_dpd_->buf4_mat_irrep_init(T2, 0);
    global_dpd_->buf4_mat_irrep_rd(T2, 0);
    /*
    fprintf(outfile, "\tGrabbed T2's...\n");
    fflush(outfile);
    */

    /* Grab PNO, Semicanonical Transforms */
    psio_address next;
    next = PSIO_ZERO;
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedMatrix temp(new Matrix(nvir_, pno));
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) temp->pointer()[0],
            nvir_ * pno * sizeof(double), next, &next);
        Q.push_back(temp->clone());
    }
    /*
    fprintf(outfile, "\tGrabbed Q's...\n");
    fflush(outfile);
    */
    next = PSIO_ZERO;
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedMatrix temp(new Matrix(pno, pno));
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) temp->pointer()[0],
            pno * pno * sizeof(double), next, &next);
        L.push_back(temp->clone());
    }
    /*
    fprintf(outfile, "\tGrabbed L's...\n");
    fflush(outfile);
    */

    /* Grab Virtual Orbital Energies */
    next = PSIO_ZERO;
    eps_vir.resize(spairs);
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedVector temp(new Vector(pno));
        psio_read(PSIF_CC_INFO, "IJ - Virtual Orbital Energies", (char *) temp->pointer(),
            pno * sizeof(double), next, &next);
        eps_vir[ij] = temp;
    }
    /*
    fprintf(outfile, "\tGrabbed virtual epsilons...\n");
    fflush(outfile);
    */
    // Print checking
    /*
    FILE * myout;
    myout = fopen("epsilons.dat", "w");
    fprintf(myout, "\tVirtual Orbital Energies:\n");
    for(ij=0; ij < spairs; ++ij) 
        eps_vir[ij]->print(myout);
    fclose(myout);
    */

    /* Grab Weak Pair List */
    int *weak_pairs = init_int_array(npairs_);
    psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *) weak_pairs, npairs_ * sizeof(int));

    /* Transform from MO to Local to Semicanonical Basis */
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];

        if(!weak_pairs[pair_idx]) {
            i = pair_idx/nocc_;
            j = pair_idx%nocc_;
            pno = dim->get(pair_idx);
            SharedMatrix atemp(new Matrix(nvir_,pno));
            SharedMatrix btemp(new Matrix(pno,pno));
            SharedMatrix t2bar(new Matrix(pno,pno));
            SharedMatrix t2temp(new Matrix(nvir_,nvir_));

            for(ab=0; ab<nvir_*nvir_; ++ab) {
                a = ab/nvir_;
                b = ab%nvir_;
                t2temp->set(a,b, T2->matrix[0][pair_idx][ab]);
            }
            /*
            fprintf(outfile, "\tT2 Matrix\n");
            t2temp->print();
            fflush(outfile);
            */

            // To PNO basis
            /*
            C_DGEMM('n', 'n', nvir_, pno, nvir_, 1.0, &(T2->matrix[0][pair_idx][0]), nvir_, Q[ij]->pointer()[0],
                nvir_, 0.0, atemp->pointer()[0], nvir_); // result is vir by pno
            */
            atemp->gemm(0, 0, 1, t2temp, Q[ij], 0);
            t2bar->gemm(1, 0, 1, Q[ij], atemp, 0);       // result is pno by pno
            /*
            fprintf(outfile, "\tResidual in PNO basis...\n");
            fflush(outfile);
            fprintf(outfile, "\tResidual in PNO Basis:\n");
            t2bar->print();
            */

            // To semi-can. basis
            btemp->gemm(0, 0, 1, t2bar, L[ij], 0);       // result is pno by pno
            t2bar->gemm(1, 0, 1, L[ij], btemp, 0);       // result is pno by pno
            /*
            fprintf(outfile, "\tResidual in S.C. basis...\n");
            fflush(outfile);
            fprintf(outfile, "\tResidual in Semicanonical Basis:\n");
            t2bar->print();
            */

            // Apply energy denominators 
            btemp->copy(t2bar);
            for(a=0; a < pno; ++a)
                for(b=0; b < pno; ++b) {
                    double val = btemp->get(a,b);
                    t2bar->set(a,b, val * 1.0/(eps_occ_->get(i) + eps_occ_->get(j)
                                             - eps_vir[ij]->get(a) - eps_vir[ij]->get(b)));
                }
            /*
            fprintf(outfile, "\tApplied denominator...\n");
            fprintf(outfile, "\tUpdate in Semicanonical Basis:\n");
            t2bar->print();
            fflush(outfile);
            */

            // Back to PNO basis
            btemp->gemm(0, 0, 1, L[ij], t2bar, 0);         // result is pno by pno
            t2bar->gemm(0, 1, 1, btemp, L[ij], 0);         // result is pno by pno
            /*
            fprintf(outfile, "\tUpdate in PNO basis...\n");
            fflush(outfile);
            */
            // Back to MO basis
            t2temp->zero();
            atemp->gemm(0, 0, 1, Q[ij], t2bar, 0);         // result is vir by pno
            t2temp->gemm(0, 1, 1, atemp, Q[ij], 0);
            /*
            C_DGEMM('n', 't', nvir_, nvir_, pno, 1.0, atemp->pointer()[0], nvir_, Q[ij]->pointer()[0],
                pno, 0.0, &(T2->matrix[0][pair_idx][0]), nvir_); // result is vir by vir
            */
            /*
            fprintf(outfile, "\tUpdate in canonical basis...\n");
            fflush(outfile);
            */
            for(ab=0; ab<nvir_*nvir_; ++ab) {
                a = ab/nvir_;
                b = ab%nvir_;
                T2->matrix[0][pair_idx][ab] = t2temp->get(a,b);
            }
        }
        else {
            for(ab=0; ab<nvir_*nvir_; ++ab) {
                a = ab/nvir_;
                b = ab%nvir_;
                T2->matrix[0][pair_idx][ab] = 0.0;
            }
        }
    }

    /* Make sure T2 data for non-significant pairs are zero */
    for(int ij=0; ij < zero_list_.size(); ++ij) {
        pair_idx = zero_list_[ij];
        for(int ab=0; ab < nvir_*nvir_; ++ab)
            T2->matrix[0][pair_idx][ab] = 0.0;
        //memset((void *) T2->matrix[0][pair_idx], 0, nvir_*nvir_*sizeof(double));
    }

    /* Write New T2's to disk */
    global_dpd_->buf4_mat_irrep_wrt(T2, 0);
    global_dpd_->buf4_mat_irrep_close(T2, 0);
    
}


/* T1 Update (stepping through PNO (and Ortho.) basis) */
void New_Local::update_singles(dpdfile2 * T1) {
    //fprintf(outfile, "\n\t*** Local Singles Update ***\n");
    //fflush(outfile);
    int ij, ii, i, j, ab, a, b;
    int pno, pair_idx, spairs;
    
    // This function will use transforms stored on disk
    // instead of the member variables, Q_, L_.
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_vir;

    /* Grab PNO Dimensions */
    SharedVector dim(new Vector("Local Pair-Virtual Dimensions",npairs_));
    psio_read_entry(PSIF_CC_INFO, "Local Dimensions", (char *) dim->pointer(), npairs_ * sizeof(double));

    /* Get Siginficant-Pair List */
    get_pair_list(dim);
    spairs = pair_list_.size();
    /*
    fprintf(outfile, "\tSurviving Pairs:\n");
    fflush(outfile);
    for(int x=0; x < spairs; ++x) {
        fprintf(outfile, "\tPair: %d\n", pair_list_[x] + 1);
        fflush(outfile);
    }
    fprintf(outfile, "\tNon-Surviving Pairs:\n");
    fflush(outfile);
    for(int x=0; x < zero_list_.size(); ++x) {
        fprintf(outfile, "\tPair: %d\n", zero_list_[x] + 1);
        fflush(outfile);
    }
    */

    /* Grab MO-Basis T1's */
    global_dpd_->file2_mat_init(T1);
    global_dpd_->file2_mat_rd(T1);

    /* Grab PNO (or OSV), Semicanonical Transforms */
    psio_address next;
    next = PSIO_ZERO;
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedMatrix temp(new Matrix(nvir_, pno));
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) temp->pointer()[0],
            nvir_ * pno * sizeof(double), next, &next);
        Q.push_back(temp->clone());
    }
    next = PSIO_ZERO;
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedMatrix temp(new Matrix(pno, pno));
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) temp->pointer()[0],
            pno * pno * sizeof(double), next, &next);
        L.push_back(temp->clone());
    }

    /* Grab Virtual Orbital Energies */
    next = PSIO_ZERO;
    eps_vir.resize(spairs);
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        pno = dim->get(pair_idx);
        SharedVector temp(new Vector(pno));
        psio_read(PSIF_CC_INFO, "IJ - Virtual Orbital Energies", (char *) temp->pointer(),
            pno * sizeof(double), next, &next);
        eps_vir[ij] = temp;
    }

    /* Grab Weak Pair List */
    int *weak_pairs = init_int_array(npairs_);
    psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *) weak_pairs, npairs_ * sizeof(int));

    /* Transform from MO to Local to Semicanonical Basis */
    int ldim; // Local dimension, generally (either PNO or OSV)
    double *T1tilde, *T1bar;
    double **Qtemp, **Ltemp;

    /* Make list of significant diagonal pairs */
    std::vector<int> diag_pairs;
    for(ij=0; ij < spairs; ++ij) {
        pair_idx = pair_list_[ij];
        i = pair_idx/nocc_;
        j = pair_idx%nocc_;
        if(i==j) {
            diag_pairs.push_back(ij);
        }
    }
    /*
    for(ii=0; ii<diag_pairs.size(); ++ii)
        fprintf(outfile, "Relative pair (diagonal) index: %d\n", diag_pairs[ii]);
    */

    // diag_idx -> relative index for accessing correct Q, L, eps_vir data.
    // pair_idx -> absolute, original pair index.
    // i        -> absolute occupied index.
    for(ii=0; ii < diag_pairs.size(); ++ii) {
        int diag_idx = diag_pairs[ii];
        pair_idx = pair_list_[diag_idx];
        i = pair_idx/nocc_;
        /*
        fprintf(outfile, "Absolute Diagonal Pair Index: %d\n", pair_idx);
        fprintf(outfile, "Relative Diagonal Pair Index: %d\n", diag_idx);
        fprintf(outfile, "Occ. Index: %d\n", i);
        fflush(outfile);
        */

        // This was replaced by tracking "diag_pairs" and keeping up with
        // relative and absolute indices in the code above.
        /*
        // Only do this if the diagonal pair dimension is not zero.
        if(pair_list_[ii]) 
            pair_idx = pair_list_[ii];
            ldim = dim->get(pair_idx);
            T1tilde = init_array(ldim);
            T1bar   = init_array(ldim);
            Qtemp = Q[pair_idx]->to_block_matrix();
            Ltemp = L[pair_idx]->to_block_matrix();
        */

        // Update T1 for non-weak diagonal pairs
        /* Hey, Weak Pairs applies to DOUBLES ONLY! */
        ldim = dim->get(pair_idx);
        T1tilde = init_array(ldim);
        T1bar   = init_array(ldim);
        Qtemp = Q[diag_idx]->to_block_matrix();
        Ltemp = L[diag_idx]->to_block_matrix();
    
        // To PNO/OSV basis
        C_DGEMV('t', nvir_, ldim, 1.0, Qtemp[0], ldim, 
        &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

        // To semi-can. basis
        C_DGEMV('t', ldim, ldim, 1.0, Ltemp[0], ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

        // Apply energy denominators 
        for(a=0; a < ldim; ++a)
            T1bar[a] /= eps_occ_->get(i) - eps_vir[diag_idx]->get(a);
            //T1bar[a] /= eps_occ_->get(i) - eps_vir[ii]->get(a);

        // Back to PNO basis
        C_DGEMV('n', ldim, ldim, 1.0, Ltemp[0], ldim,
        &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

        // Back to MO basis
        C_DGEMV('n', nvir_, ldim, 1.0, Qtemp[0], ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1->matrix[0][i][0]), 1);

/*
        // To PNO Basis
        C_DGEMV('t', nvir_, ldim, 1.0, &(Q[ii]->pointer()[0][0]), ldim, 
        &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

        // To semi-can. basis
        C_DGEMV('t', ldim, ldim, 1.0, &(L[ii]->pointer()[0][0]), ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

        // Apply energy denominators 
        for(a=0; a < ldim; ++a)
            T1bar[a] /= eps_occ_->get(i) - eps_vir[ii]->get(a);

        // Back to PNO basis
        C_DGEMV('n', ldim, ldim, 1.0, &(L[ii]->pointer()[0][0]), ldim,
        &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

        // Back to MO basis
        C_DGEMV('n', nvir_, ldim, 1.0, &(Q[ii]->pointer()[0][0]), ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1->matrix[0][i][0]), 1);
*/

        free(T1bar);
        free(T1tilde);
    }      // End ii-loop.
    free_block(Qtemp);
    free_block(Ltemp);
    
    /* Write out the updated T1's */
    global_dpd_->file2_mat_wrt(T1);
    global_dpd_->file2_mat_close(T1);
    
}

/* Singles Update for Separate Singles Cutoff */
void New_Local::sep_update_singles(dpdfile2 * T1) {
    /* I'm going to assume that all diagonal pairs are significant/kept. */
    int ij, ii, i, j, ab, a, b;
    int pno, pair_idx, spairs;
    
    // This function will use transforms stored on disk
    // instead of the member variables, Q_, L_.
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_vir;

    /* Grab PNO Dimensions */
    SharedVector dim(new Vector("Local Singles-Virtual Dimensions",nocc_));
    psio_read_entry(PSIF_CC_INFO, "Local Singles Dimensions", (char *) dim->pointer(), nocc_ * sizeof(double));

    /* Grab MO-Basis T1's */
    global_dpd_->file2_mat_init(T1);
    global_dpd_->file2_mat_rd(T1);

    /* Grab PNO (or OSV), Semicanonical Transforms */
    psio_address next;
    next = PSIO_ZERO;
    for(i=0; i < nocc_; ++i) {
        pno = dim->get(i);
        SharedMatrix temp(new Matrix(nvir_, pno));
        psio_read(PSIF_CC_INFO, "Local Singles Transformation Matrix (Q)", (char *) temp->pointer()[0],
            nvir_ * pno * sizeof(double), next, &next);
        Q.push_back(temp->clone());
    }
    next = PSIO_ZERO;
    for(i=0; i < nocc_; ++i) {
        pno = dim->get(i);
        SharedMatrix temp(new Matrix(pno, pno));
        psio_read(PSIF_CC_INFO, "Semicanonical Singles Transformation Matrix (L)", (char *) temp->pointer()[0],
            pno * pno * sizeof(double), next, &next);
        L.push_back(temp->clone());
    }

    /* Grab Virtual Orbital Energies */
    next = PSIO_ZERO;
    eps_vir.resize(nocc_);
    for(i=0; i < nocc_; ++i) {
        pno = dim->get(i);
        SharedVector temp(new Vector(pno));
        psio_read(PSIF_CC_INFO, "II - Virtual Orbital Energies", (char *) temp->pointer(),
            pno * sizeof(double), next, &next);
        eps_vir[i] = temp;
    }

    /* Grab Weak Pair List */
    // Diagonal pairs aren't/shouldn't be weak!
    /*
    int *weak_pairs = init_int_array(npairs_);
    psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *) weak_pairs, npairs_ * sizeof(int));
    */

    /* Transform from MO to Local to Semicanonical Basis */
    int ldim; // Local dimension, generally (either PNO or OSV)
    double *T1tilde, *T1bar;
    double **Qtemp, **Ltemp;

    for(i=0; i < nocc_; ++i) {

        /* Hey, Weak Pairs applies to DOUBLES ONLY! */
        ldim = dim->get(i);
        T1tilde = init_array(ldim);
        T1bar   = init_array(ldim);
        Qtemp = Q[i]->to_block_matrix();
        Ltemp = L[i]->to_block_matrix();
    
        // To PNO/OSV basis
        C_DGEMV('t', nvir_, ldim, 1.0, Qtemp[0], ldim, 
        &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

        // To semi-can. basis
        C_DGEMV('t', ldim, ldim, 1.0, Ltemp[0], ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

        // Apply energy denominators 
        for(a=0; a < ldim; ++a)
            T1bar[a] /= eps_occ_->get(i) - eps_vir[i]->get(a);

        // Back to PNO basis
        C_DGEMV('n', ldim, ldim, 1.0, Ltemp[0], ldim,
        &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

        // Back to MO basis
        C_DGEMV('n', nvir_, ldim, 1.0, Qtemp[0], ldim, 
        &(T1tilde[0]), 1, 0.0, &(T1->matrix[0][i][0]), 1);

        free(T1bar);
        free(T1tilde);
    }      // End i-loop.
    free_block(Qtemp);
    free_block(Ltemp);
    
    /* Write out the updated T1's */
    global_dpd_->file2_mat_wrt(T1);
    global_dpd_->file2_mat_close(T1);
}

/// IGNORING THIS VERSION FOR NOW.
 // -> If you want to update this one to be like the above function, probably best to just
 //    recopy the function above and change it to be in keeping with the process of this one.
void New_Local::alt_update_pno(dpdbuf4 * T2) {
    int ij, i, j, ab, a, b;
    int pno;

    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_vir(npairs_);
    std::vector<SharedMatrix> t2can;

    /* Grab MO-Basis T2's */
    get_matvec(T2, &t2can);

    /* Grab PNO Dimensions */
    SharedVector dim(new Vector("PNO Dimensions",npairs_));
    psio_read_entry(PSIF_CC_INFO, "Local Dimensions", (char *) dim->pointer(), npairs_ * sizeof(double));

    /* Grab PNO, Semicanonical Transforms and Virtual Orbital Energies */
    psio_address next;
    next = PSIO_ZERO;
    for(ij=0; ij < npairs_; ++ij) {
        pno = dim->get(ij);
        SharedMatrix temp(new Matrix(pno, pno));
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (Q)", (char *) temp->pointer()[0],
            pno * pno * sizeof(double), next, &next);
        Q.push_back(temp->clone());
    }
    next = PSIO_ZERO;
    for(ij=0; ij < npairs_; ++ij) {
        pno = dim->get(ij);
        SharedMatrix temp(new Matrix(pno, pno));
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix (L)", (char *) temp->pointer()[0],
            pno * pno * sizeof(double), next, &next);
        L.push_back(temp->clone());
    }
    next = PSIO_ZERO;
    for(ij=0; ij < npairs_; ++ij) {
        pno = dim->get(ij);
        SharedVector temp(new Vector(pno));
        psio_read(PSIF_CC_INFO, "IJ - Virtual Orbital Energies", (char *) temp->pointer(),
            pno * sizeof(double), next, &next);
        eps_vir[ij] = temp;
    }

    global_dpd_->buf4_mat_irrep_init(T2, 0);
    global_dpd_->buf4_mat_irrep_rd(T2, 0);

    /* Transform from MO to PNO, Semicanonical Basis */
    for(ij=0; ij < npairs_; ++ij) {
        pno = dim->get(ij);
        SharedMatrix atemp(new Matrix(nvir_,pno));
        SharedMatrix btemp(new Matrix(pno,pno));
        SharedMatrix t2bar(new Matrix(pno,pno));

        // To PNO basis
        atemp->gemm(0, 0, 1, t2can[ij], Q[ij], 0);  // result is vir by pno
        t2bar->gemm(1, 0, 1, Q[ij], atemp, 0);      // result is pno by pno
        // To semi-can. basis
        btemp->gemm(0, 0, 1, t2bar, L[ij], 0);      // result is pno by pno
        t2bar->gemm(1, 0, 1, L[ij], btemp, 0);      // result is pno by pno
    
        // Apply energy denominators */
        for(a=0; a < pno; ++a)
            for(b=0; b < pno; ++b) {
                double val = t2bar->get(a,b);
                t2bar->set(a,b, val * 1.0/(eps_occ_->get(i) + eps_occ_->get(j)
                                         - eps_vir[ij]->get(a) - eps_vir[ij]->get(b)));
            }

        // Back to PNO basis
        btemp->gemm(0, 0, 1, L[ij], t2bar, 0);      // result is pno by pno
        t2bar->gemm(0, 1, 1, btemp, L[ij], 0);      // result is pno by pno
        // Back to MO basis
        atemp->gemm(0, 0, 1, Q[ij], t2bar, 0);      // result is vir by pno
        t2can[ij]->gemm(0, 1, 1, atemp, Q[ij], 0);  // result is vir by vir

        for(a=0,ab=0; a<nvir_; ++a)
            for(b=0; b<nvir_; ++b, ++ab)
                T2->matrix[0][ij][ab] = t2can[ij]->get(a,b);
    }

    global_dpd_->buf4_mat_irrep_wrt(T2, 0);
    global_dpd_->buf4_mat_irrep_close(T2, 0);
    
}


/* Copy from dpdbuf4 to vector of matrices */
   // Assumes virtual by virtual matrix for each occupied pair,
   // BUT needs to be generalized.
void New_Local::get_matvec(dpdbuf4 *myfile, std::vector<SharedMatrix> *matvec) {
    int i,j,ij;
    int a,b,ab;
 
    global_dpd_->buf4_mat_irrep_init(myfile, 0);
    global_dpd_->buf4_mat_irrep_rd(myfile, 0);
    SharedMatrix mat(new Matrix(nvir_,nvir_));
 
    for(ij=0; ij < npairs_; ++ij) {
        for(ab=0; ab < nvir_*nvir_; ++ab) {
            a = ab/(nvir_);
            b = ab%(nvir_);
            mat->set(a,b, myfile->matrix[0][ij][ab]);
        }
        matvec->push_back(mat->clone());
    }

    global_dpd_->buf4_mat_irrep_close(myfile, 0);

    // Print checking
    /*
    fprintf(outfile, "\t***** IJ Matrices *****\n");
    for(ij=0; ij<matvec->size(); ++ij)
        matvec->at(ij)->print();
    */
}

/* Only Keep up with significant pairs */
   // Currently only determined by whether or not PNO/OSV dimension is zero.
void New_Local::get_pair_list(SharedVector dim) {
    int ij;
    if(pair_list_.size()) pair_list_.clear();
    if(zero_list_.size()) zero_list_.clear();

    //for(ij=0; ij < npairs_; ++ij) {
    for(ij=0; ij < dim->dim(); ++ij) {
        if(dim->get(ij)) {
            pair_list_.push_back(ij);
        }
        else zero_list_.push_back(ij);
    }

}


/* Takes list of MP2 pair-energy-screened pairs and removes pairs with zero PNO/OSV */
void New_Local::update_pair_list(SharedVector dim) {
    int ij, pair_idx;
    std::vector<int> new_list;

    for(ij=0; ij < pair_list_.size(); ++ij) {
        pair_idx = pair_list_[ij];
        if(dim->get(pair_idx)) {
            new_list.push_back(pair_idx);
        }
    }
    pair_list_.clear();
    pair_list_ = new_list;
    
}

void New_Local::get_occ_epsilons(void) {
    // Unnecessary, b/c PSIF_CC_OEI only has active orbital energies.
    /*
    int nfzc = moinfo_.frdocc[0];
    fprintf(outfile, "\tGetting occupied orbital energies...\n");
    fprintf(outfile, "\tNocc = %d, Nfzc = %d\n", nocc_, nfzc);
    fflush(outfile);
    */
    eps_occ_ = SharedVector(new Vector("Local Occ. Fock Diagonals",nocc_));

    dpdfile2 fij;
    SharedMatrix focc(new Matrix(nocc_,nocc_));
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);
    for(int i=0; i < nocc_; ++i)
        eps_occ_->set(i, fij.matrix[0][i][i]);
        //eps_occ_->set(i, fij.matrix[0][i+nfzc][i+nfzc]);
    global_dpd_->file2_close(&fij);
}

void New_Local::lmp2_pair_energy(std::vector<SharedMatrix> T2mats, std::vector<SharedMatrix> Dints) {
    int ij;
    double eij = 0.0;
    //pair_cut_ = 0.0001;
    SharedVector pair_energies(new Vector("LMP2 Pair Energies", npairs_));
    //int *weak_pairs = init_int_array(npairs_); /* 0 for strong, 1 for weak */
    weak_pairs_ = init_int_array(npairs_); /* 0 for strong, 1 for weak */

    int nwp = 0;
    for(ij=0; ij < npairs_; ++ij) {
        eij = T2mats[ij]->vector_dot(Dints[ij]);        
        pair_energies->set(ij, eij);
        if(fabs(eij) < pair_cut_) {
            weak_pairs_[ij] = 1;
            nwp += 1;
        }
    }
    /*
    T2mats[1]->print();
    Dints[1]->print();
    */
    pair_energies->print();
    //fprintf(outfile, "\t|Pair Energy| cutoff: %3.1e\n", pair_cut_);
    //fprintf(outfile, "\tNumber of Weak Pairs: %d\n", nwp);
    if(nwp) {
        //fprintf(outfile, "\t*** Weak Pairs ***\n");
        for(ij=0; ij < npairs_; ++ij) {
            if(weak_pairs_[ij]) 
                {//fprintf(outfile, "\t\t %2d\n", ij+1);
                }
        }
    }
    //fprintf(outfile, "\n");

    psio_address next;
    psio_write_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *) weak_pairs_, npairs_ * sizeof(int));
    //free(weak_pairs_);  // Wait to free until after computing T2 Length, # Average PNO's, etc.
}

New_Local::~New_Local() {
}

}} // namespace psi::ccenergy


