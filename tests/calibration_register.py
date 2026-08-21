"""The register of checking instruments and the case that calibrates each.

Split from ``test_calibration.py`` at the 500-line limit and along the seam that
file already had: the tests there are the mechanism — every public callable in a
checking module must appear here, and every named test must exist — and this is
the register itself, which is the part a reader checks by eye.

An instrument that has never been run on a case whose answer is known
independently is worse than none: it manufactures findings. Every entry maps
``module.callable`` to the test that drives it against such a case.
"""

from __future__ import annotations

__all__ = ["CALIBRATED"]

#: instrument -> the test that drives it against a known answer.
CALIBRATED = {
    "homology.reviewed_family":
        "test_homology.py::test_nine_are_reviewed_and_the_tenth_is_a_real_paralogue",
    "homology.unreviewed_family":
        "test_homology.py::test_nine_are_reviewed_and_the_tenth_is_a_real_paralogue",

    # --- family_constraint: replication on boundaries the census did not choose
    "family_constraint.domain_constraint":
        "test_family_constraint.py::test_the_pore_machinery_is_the_most_constrained_on_our_boundaries",
    "family_constraint.paralogue_asymmetry":
        "test_family_constraint.py::test_our_identities_reproduce_the_census_where_it_matters",
    "family_constraint.blade_gradient":
        "test_family_constraint.py::test_the_distal_blade_finding_is_a_property_of_the_bands_not_the_blades",
    "family_constraint.constraint_on_structure":
        "test_family_constraint.py::test_the_track_is_refused_on_an_entry_it_cannot_be_read_on",
    "family_constraint.compare_with_own_conservation":
        "test_family_constraint.py::test_our_own_conservation_and_the_census_agree",
    "family_constraint.census_domain_constraint": "test_family_findings.py",
    "family_constraint.selection_track": "test_family_findings.py",

    # --- constraint_mechanics: a correlation, with the null that makes it one
    "constraint_mechanics.circular_shift_null":
        "test_constraint_mechanics.py::test_the_shift_null_is_centred_on_zero_for_an_unrelated_feature",
    "constraint_mechanics.partial_spearman":
        "test_constraint_mechanics.py::test_the_partial_correlation_removes_a_correlation_that_is_only_burial",
    "constraint_mechanics.couple":
        "test_constraint_mechanics.py::test_the_verdict_says_no_when_nothing_survives",
    "constraint_mechanics.align_track_to_features":
        "test_constraint_mechanics.py::test_the_track_lines_up_with_the_feature_table_through_the_alignment",

    # --- core_periphery: a ratio, with the case where its denominator fails
    "core_periphery.compare":
        "test_family_structure.py::test_a_pair_whose_pore_modules_do_not_superpose_is_refused_a_splay_ratio",
    "core_periphery.core_residues":
        "test_family_structure.py::test_the_core_and_periphery_sets_do_not_overlap",
    "core_periphery.periphery_residues":
        "test_family_structure.py::test_the_core_and_periphery_sets_do_not_overlap",
    "core_periphery.correspondence":
        "test_family_structure.py::test_an_entry_against_itself_gives_a_zero_core_and_no_splay",
    # The fit with its transform kept. Calibrated on a known answer the
    # instrument cannot fake: the same pair must fit the same either way
    # round, because a rigid fit on one residue set is symmetric — and it did
    # not, until the core stopped being selected in whichever numbering the
    # mobile happened to be in.
    "core_periphery.core_fit":
        "test_family_structure.py::test_the_same_pair_fits_the_same_either_way_round",

    # --- equivalent_positions: the register, not the distance
    "equivalent_positions.locate":
        "test_family_structure.py::test_the_control_shows_the_whole_core_superposes_not_just_these_two",
    "equivalent_positions.alignment_agrees":
        "test_family_structure.py::test_this_projects_own_alignment_pairs_the_same_residues_as_the_census",
    "equivalent_positions.map_position":
        "test_family_structure.py::test_the_claimed_pairs_land_within_one_residue_after_a_core_fit",

    # --- disease_geography: a Fisher test, against a hand-computable table
    "disease_geography.fisher_exact_greater":
        "test_family_constraint.py::test_fisher_exact_matches_a_table_computed_by_hand",
    "disease_geography.pore_module_enrichment":
        "test_family_constraint.py::test_the_two_partitions_disagree_and_the_disputed_band_holds_disease",
    "disease_geography.both_partitions":
        "test_family_constraint.py::test_the_two_partitions_disagree_and_the_disputed_band_holds_disease",
    "disease_geography.boundary_disagreement":
        "test_family_constraint.py::test_the_boundary_disagreement_is_confined_to_the_anchor_end",
    "disease_geography.constraint_classifier":
        "test_family_constraint.py::test_the_constraint_score_still_classifies_against_a_different_negative_set",
    "disease_geography.pathogenic_positions":
        "test_family_constraint.py::test_pathogenic_positions_exclude_truncating_variants",
    "disease_geography.pore_module_residues":
        "test_family_constraint.py::test_the_two_partitions_disagree_and_the_disputed_band_holds_disease",
    "disease_geography.population_positions": "test_family_constraint.py",
    "disease_geography.census_comparison": "test_family_findings.py",

    # --- family_motifs: an absence, with a positive control beside it
    "family_motifs.motif_scan":
        "test_family_findings.py::test_the_absent_motif_search_has_a_control_that_finds_something",
    "family_motifs.control_motif":
        "test_family_findings.py::test_the_absent_motif_search_has_a_control_that_finds_something",
    "family_motifs.deep_windows":
        "test_family_findings.py::test_what_is_conserved_to_family_depth_is_the_pore_machinery",

    # --- piezo3: a paralogue with no experiment behind it
    "piezo3.census_to_model":
        "test_family_structure.py::test_the_two_piezo3_records_are_mapped_by_alignment_not_by_an_offset",
    "piezo3.model_to_census":
        "test_family_structure.py::test_the_two_piezo3_records_are_mapped_by_alignment_not_by_an_offset",
    "piezo3.kept_positions":
        "test_family_structure.py::test_piezo3_keeps_the_human_residue_at_all_fourteen_pore_positions",
    "piezo3.fold_comparison":
        "test_family_structure.py::test_the_piezo3_fold_agrees_with_piezo1_at_the_core_and_not_the_blades",
    "piezo3.template_survey":
        "test_family_structure.py::test_the_worm_template_is_measurably_worse_than_a_paralogue_one",
    "piezo3.best_paralogue_template":
        "test_family_structure.py::test_the_worm_template_is_measurably_worse_than_a_paralogue_one",
    "piezo3.load_model": "test_family_structure.py",
    "piezo3_channel.build_channel":
        "test_family_structure.py::test_the_assembled_channel_says_how_much_of_itself_it_borrowed",
    "piezo3_channel.measure_channel":
        "test_family_structure.py::test_the_comparison_entry_reproduces_this_projects_own_dome_number",

    # crosscheck: re-derivations of physics results
    "crosscheck.dome_curvature_by_cap_geometry":
        "test_crosscheck.py::test_the_exact_cap_route_is_also_exact_at_every_slope",
    "crosscheck.dome_curvature_by_parabola":
        "test_crosscheck.py::test_the_parabola_route_degrades_with_slope",
    "crosscheck.gating_overlap_by_distances":
        "test_crosscheck.py::test_distance_route_recovers_a_planted_mode",
    "crosscheck.t50_by_ode_integration":
        "test_crosscheck.py::test_ode_integration_reproduces_the_matrix_exponential",
    "crosscheck.t50_by_steady_state":
        "test_crosscheck.py::test_the_steady_state_is_a_different_quantity",
    "crosscheck.compare": "test_crosscheck.py",

    # crosscheck_methods: re-derivations of algorithms
    "crosscheck_methods.pore_radius_by_random_search":
        "test_crosscheck_methods.py (slab with a known gap)",
    "crosscheck_methods.sasa_by_monte_carlo":
        "test_crosscheck_methods.py (isolated sphere, 4 pi r^2)",
    "crosscheck_methods.conservation_by_kmer_anchoring":
        "test_crosscheck_methods.py (random sequences isolate the bias)",
    "crosscheck_methods.pc1_by_power_iteration":
        "test_crosscheck_methods.py (planted principal component)",
    "crosscheck_methods.check_pore_radius": "test_crosscheck_methods.py",
    "crosscheck_methods.check_sasa": "test_crosscheck_methods.py",
    "crosscheck_methods.check_conservation": "test_crosscheck_methods.py",
    "crosscheck_methods.check_pca": "test_crosscheck_methods.py",

    # model_error
    "model_error.fit_spheroid":
        "test_model_error.py (known semi-axes, incl. the original failure)",
    "model_error.spring_model_error":
        "test_calibration.py::test_spring_model_error_recovers_a_planted_mode",
    "model_error.pore_convention_error": "test_model_error.py (the null's mechanism)",
    "model_error.dome_model_error": "test_model_error.py",
    "model_error.compare_with_sampling": "test_model_error.py",

    # uncertainty
    "uncertainty.bootstrap":
        "test_uncertainty.py::test_bootstrap_recovers_a_known_interval",
    "uncertainty.sensitivity":
        "test_calibration.py::test_sensitivity_range_is_exactly_the_known_extremes",
    "uncertainty.parameter_range":
        "test_uncertainty.py::test_parameter_range_restores_even_when_the_statistic_raises",
    "uncertainty.format_with_interval": "test_uncertainty.py",

    # validation
    "validation.permutation_test":
        "test_calibration.py::test_permutation_p_matches_exhaustive_enumeration",
    "validation.cliffs_delta": "test_validation.py::test_cliffs_delta_extremes",
    "validation.bootstrap_cliffs_delta": "test_validation.py",
    "validation.auroc": "test_validation.py::test_auroc_known_cases",
    "validation.interpret_delta": "test_validation.py::test_delta_interpretation_thresholds",

    # homology: the family, and whether a percentage may be believed
    "homology.align_pair":
        "test_homology.py::test_the_null_is_calibrated_on_two_known_answers",
    "homology.shuffled_null":
        "test_homology.py::test_the_null_can_say_no_when_there_is_no_homology",
    "homology.relationship":
        "test_homology.py::test_at_least_one_pair_has_an_identity_indistinguishable_from_chance",
    "homology.family_matrix": "test_homology.py (36 pairs, each against its null)",
    "homology.family": "test_homology.py::test_the_family_is_ten_and_every_member_is_committed",
    "homology.member": "test_homology.py (family lookup)",
    "homology.group_of": "test_homology.py::test_no_two_members_share_a_length_or_a_helix_count_architecture",

    # assembly: is a trimer built from one protomer a model of anything?
    "assembly.assemble_trimer":
        "test_assembly.py::test_a_protomer_rebuilt_on_its_own_trimer_reproduces_it_exactly",
    "assembly.borrowed_fraction":
        "test_assembly.py::test_most_of_an_assembly_s_shape_is_the_template_s_and_it_says_so",
    "assembly.best_template":
        "test_assembly.py::test_the_template_is_chosen_by_protein_then_by_coverage",
    "clashes.count_clashes":
        "test_assembly.py::test_the_clash_counter_is_near_zero_on_real_trimers",
    "assembly.is_monomer":
        "test_assembly.py::test_a_real_trimer_is_refused_rather_than_replaced_by_a_model_of_itself",

    # homology_structure: is one entry pair a measurement of two proteins?
    "homology_structure.mode_overlap_spread":
        "test_homology_structure.py::test_the_spread_reports_instability_and_piezo2_is_the_positive_control",
    "homology_structure.compare_structures":
        "test_homology_structure.py::test_a_single_entry_pair_is_not_a_measurement_of_the_two_proteins",
    "homology_structure.index_pairing_valid":
        "test_homology_structure.py::test_index_pairing_is_refused_where_the_helix_counts_differ",
    "homology_structure.comparable_entries":
        "test_homology_structure.py::test_comparable_entries_are_best_resolved_first",
    "homology_structure.helix_counts":
        "test_homology.py::test_no_two_members_share_a_length_or_a_helix_count_architecture",

    # alignment_windows: is the alignment in register here?
    "alignment_windows.window_score":
        "test_homology.py::test_the_window_finds_a_planted_block_that_whole_sequence_identity_misses",
    "alignment_windows.window_identity":
        "test_homology.py::test_the_window_is_calibrated_on_a_sequence_against_itself",
    "alignment_windows.alignment_windows":
        "test_homology.py::test_the_window_refuses_a_composition_matched_shuffle",
    "alignment_windows.column_scores":
        "test_homology.py (self-alignment scores above its own null)",
    "alignment_windows.window_null_distribution":
        "test_homology.py::test_window_width_was_chosen_by_a_power_scan_not_by_taste",

    # design
    "design.power_curve":
        "test_design.py::test_null_effect_gives_the_nominal_false_positive_rate",
    "design.minimum_detectable_effect":
        "test_calibration.py::test_minimum_detectable_effect_delivers_its_target_power",
    "design.sample_size_for": "test_design.py",
    "design.shift_for_delta": "test_design.py::test_simulated_effect_matches_the_requested_effect",
    "design.delta_for_shift": "test_design.py",
    "design.benjamini_hochberg": "test_design.py (known p-value sets)",
    "design.leave_one_out": "test_design.py::test_leave_one_out_recovers_a_real_signal",

    # provenance_chain
    "provenance_chain.number_in_document":
        "test_provenance_chain.py::test_it_finds_a_number_written_in_any_reasonable_form",
    "provenance_chain.resolved_keys":
        "test_provenance_chain.py::test_a_registered_parameter_that_nothing_reads_is_detected",
    "provenance_chain.unwired_parameters": "test_provenance_chain.py",
    "provenance_chain.record_sources":
        "test_provenance_chain.py::test_it_records_the_parameters_a_computation_actually_reads",
    "provenance_chain.trace": "test_provenance_chain.py",
    "provenance_chain.walk": "test_provenance_chain.py",
    "provenance_chain.git_state": "test_provenance_chain.py",

    # dead_code
    "dead_code.audit":
        "test_dead_code.py::test_the_detector_is_calibrated",
    "dead_code.calibration":
        "test_dead_code.py::test_the_detector_is_calibrated",
    "dead_code.reference_counts":
        "test_dead_code.py::test_same_file_references_count",

    # fluctuations: the elastic network against the deposited B-factors.
    # A checking instrument in the strict sense — it is the standard test of
    # whether the network describes the molecule — so it is registered here
    # rather than treated as an analysis that happens to produce a number.
    "fluctuations.pearson":
        "test_fluctuations.py::test_the_correlations_return_their_analytic_limits",
    "fluctuations.spearman":
        "test_fluctuations.py::test_ranks_average_ties_rather_than_inventing_an_order",
    "fluctuations.predicted_msf":
        "test_fluctuations.py::"
        "test_the_comparison_recovers_a_planted_fluctuation_and_a_bad_network_does_not",
    "fluctuations.contact_number":
        "test_fluctuations.py::"
        "test_contact_number_is_a_control_and_not_a_copy_of_the_prediction",
    "fluctuations.assess_b_factors":
        "test_fluctuations.py::"
        "test_a_predicted_model_is_refused_and_the_gate_points_the_right_way",
    "fluctuations.observed_b_factors":
        "test_fluctuations.py::"
        "test_observed_b_factors_are_read_per_residue_not_per_atom",
    "fluctuations.compare_fluctuations":
        "test_fluctuations.py::test_a_shuffled_observation_correlates_with_nothing",
    "fluctuations.survey_fluctuations":
        "test_fluctuations.py::test_the_survey_of_every_downloaded_entry",

    # hydropathy: the 4-TM repeat test is a checking instrument — it is the
    # evidence for the nine-unit architecture the whole project is built on,
    # including the distal blade the full-length model grafts.
    "hydropathy.repeat_periodicity":
        "test_hydropathy.py::test_the_repeat_test_finds_a_planted_period "
        "and ::test_the_repeat_test_says_no_to_a_true_null",
    "hydropathy.hydropathy_profile":
        "test_hydropathy.py::test_a_uniform_sequence_gives_its_own_value_everywhere",
    "hydropathy.predict_segments":
        "test_hydropathy.py::"
        "test_a_planted_hydrophobic_block_is_found_where_it_was_planted",
    "hydropathy.compare_with_reference":
        "test_hydropathy.py::"
        "test_the_default_threshold_recovers_few_helices_and_that_is_reported",
    "hydropathy.annotated_hydropathy":
        "test_hydropathy.py::test_piezo1_helices_sit_below_the_conventional_membrane_cut",
    "hydropathy.threshold_scan":
        "test_hydropathy.py::test_the_threshold_scan_trades_recall_for_specificity_monotonically",
    "hydropathy.load_reference":
        "test_hydropathy.py::test_a_missing_reference_names_what_it_looked_for",

    # planarity: Figure 4a's claim as residuals. A checking instrument because
    # it decides whether a published statement about shape holds, and because
    # its first version reported a tautology (points replicated three times)
    # as a control.
    "planarity.fit_plane":
        "test_planarity.py::test_a_known_out_of_plane_displacement_is_recovered",
    "planarity.planarity":
        "test_planarity.py::"
        "test_the_flattened_structure_is_the_control_that_makes_this_mean_something",
    "planarity.blade_dependence":
        "test_planarity.py::test_coverage_decides_the_answer_and_the_module_says_so",
    "planarity.beam_angle":
        "test_planarity.py::"
        "test_the_beam_angle_opens_towards_90_when_the_channel_flattens",

    # architecture: the helix detector, calibrated on analytic helices of
    # known rise, radius and turn before it is used on coordinates.
    "architecture.helical_windows":
        "test_cross_helices.py::test_the_estimator_is_unbiased_on_an_ideal_helix",
    "architecture.helical_segments":
        "test_cross_helices.py::test_the_beam_is_found_as_one_helix_not_two",
    "architecture.cross_helices":
        "test_cross_helices.py::test_the_beam_is_not_reported_as_a_cross_helix",
    "architecture.cross_helix_scan":
        "test_cross_helices.py::"
        "test_the_threshold_is_reported_against_a_scan_not_asserted",
    "architecture.ideal_helix":
        "test_cross_helices.py::test_the_ideal_helix_generator_is_what_it_claims",

    # micelle: the modelled Figure 4b envelope. A measuring instrument rather
    # than a checking one, but it produces a curvature that gets compared with
    # a published number, so it is registered on the same terms as the rest.
    "micelle.distance_field":
        "test_ui_membrane_views.py::"
        "test_the_envelope_of_one_point_is_a_sphere_of_the_offset_radius",
    "micelle.build_micelle":
        "test_ui_membrane_views.py::"
        "test_the_micelle_encloses_the_belt_and_says_it_is_a_model",
    "micelle.belt_atoms":
        "test_ui_membrane_views.py::"
        "test_the_belt_is_apolar_transmembrane_side_chains",

    # parameter_effect
    "parameter_effect.measure_effect":
        "test_parameter_effect.py::test_a_parameter_with_no_effect_is_reported_as_such",
    "parameter_effect.probe_effects": "test_parameter_effect.py",
}
