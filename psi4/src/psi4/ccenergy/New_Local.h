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

#ifndef _psi_src_bin_ccenergy_new_local_h
#define _psi_src_bin_ccenergy_new_local_h

namespace psi { namespace ccenergy {

class New_Local {
    /// Basic Variables
    int nocc_;   // Active double occupied.
    int nvir_;   // Active virtuals.
    int npairs_; 
    double pno_cut_;
    double osv_cut_;
    double pair_cut_;
    std::string local_type_;
    std::string osv_type_;
    /// Different threshold for singles
    bool separate_singles_;
    double singles_cut_;

    /// Transforms
    std::vector<SharedMatrix> Q_;  // PNO or OSV Transforms
    std::vector<SharedMatrix> L_;  // Semicanonical Transforms
    std::vector<SharedMatrix> X_;  // Orthogonal Transforms
    std::vector<SharedMatrix> Qs_;  // PNO or OSV Transforms (singles specific)
    std::vector<SharedMatrix> Ls_;  // Semicanonical Transforms (singles specific)

    /// Miscellaneous Essentials
    std::vector<SharedVector> occ_num_; // Array of PNO,OSV occupation numbers.
    SharedVector dim_;                  // Array of truncated virtual dimensions per each ij-pair.
    std::vector<int>  nr_dim_;          // Array of trunc., lin. independ. virt. dim. per each ij-pair.
    SharedVector eps_occ_;              // Virtual orbital energies.
    std::vector<SharedVector> eps_vir_; // Occupied orbital energies.
    std::vector<int> pair_list_;        // List of significant/surviving pairs.
    std::vector<int> zero_list_;        // List of eliminated pairs.
    int *weak_pairs_;                   // List of weak pairs to be neglected.
    /// Sep. Singles Misc.
    SharedVector t1dim_;
    std::vector<SharedVector> sep_eps_vir_;

    /// Utility functions
    void get_occ_epsilons();
    void get_matvec(dpdbuf4 *, std::vector<SharedMatrix> *);
    void lmp2_pair_energy(std::vector<SharedMatrix>, std::vector<SharedMatrix>);
    void get_semicanonical_transform();
    void get_pair_list(SharedVector);
    void update_pair_list(SharedVector);
    /// Sep. Singles
    void init_sep_singles(std::vector<SharedMatrix>);
    void sep_singles_semicanonical_transform();

public:
    /// Constructors
     // -> Sets basic variables.
     // -> Grabs eps_occ_. 
     // -> If bool is true, calls init_pno().
    New_Local();
    New_Local(bool);
    New_Local(std::string);
    New_Local(std::string, bool);

    /// Initial PNO
     // -> Computes and writes Q_, L_, dim_, and eps_vir_ to CC_INFO.
    void init_pno();
    /// PNO update procedure
     // -> Simulates local calculation (a la local_filter_T2()).
     // -> Reads Q, L, dim, and eps_vir from CC_INFO.
    void update_pno(dpdbuf4 *);       // Handles dpdbuf4 T2 directly
    void alt_update_pno(dpdbuf4 *);   // Everything in terms of SharedMatrix
    void update_singles(dpdfile2 *);  // "Filter" the singles, handles T1 directly..

    /// Initial OSV
     // -> Computes and writes Q_, L_, dim_, and eps_vir_ to CC_INFO.
    void init_osv();
    /// PNO update procedure
     // -> Simulates local calculation (a la local_filter_T2()).
     // -> Reads Q, L, dim, and eps_vir from CC_INFO.
    void update_osv(dpdbuf4 *);     // Handles dpdbuf4 T2 directly

    /// Singles Alternative Version
    void sep_update_singles(dpdfile2 *);

    /// Destructor
     // Shouldn't need to do anything...
    ~New_Local();
    
};

}} // namespace psi::ccenergy

#endif // _psi_src_bin_ccenergy_new_local_h
