list(APPEND sources
  atom_data.cc
  bend.cc
  cart.cc
  combo_coordinates.cc
  fb_frag.cc
  frag.cc
  frag_H_guess.cc
  frag_apply_frozen_constraints.cc
  frag_disp.cc
  frag_natural.cc
  frag_print.cc
  geom_gradients_io.cc
  getIntcoFileName.cc
  interfrag.cc
  interfrag_orient.cc
  lindh_guess.cc
  linear_algebra.cc
  mem.cc
  molecule.cc
  molecule_backstep.cc
  molecule_fragments.cc
  molecule_irc_step.cc
  molecule_linesearch_step.cc
  molecule_nr_step.cc
  molecule_prfo_step.cc
  molecule_print.cc
  molecule_read_coords.cc
  molecule_rfo_step.cc
  molecule_sd_step.cc
  molecule_tests.cc
  oofp.cc
  opt_data.cc
  opt_data_io.cc
  optking.cc
  print.cc
  print.cc
  set_params.cc
  stre.cc
  tors.cc
  v3d.cc
  )
psi4_add_module(bin optking sources)
