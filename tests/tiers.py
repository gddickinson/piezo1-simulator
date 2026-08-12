"""Which situation each test file belongs to, so short runs are selections
with nothing falling between them.

The full suite has grown to 129 files and tens of minutes, which is the right
cost for an occasional check and the wrong one for "did this edit break
anything obvious". But a shorter suite is a **selection**, and a selection can
silently drop things — a test no situation runs would decay without anyone
choosing that. So the tiers are a *partition*: every test file belongs to
exactly one, ``test_tiers.py`` fails on any file left out or listed twice,
and ``make test`` runs everything with no filter at all, exactly as before.

Selection is by file, through ``--suite`` (see ``conftest.py``), because a
file is the unit this suite is organised in: each one covers a module and
carries its own skips for missing data.

The tiers are situations, not speeds — though ``quick`` is also the fast one:

- ``quick``    — sanity on every edit. No Qt, no GL, no heavy computation on
                 real coordinates; under a minute together.
- ``science``  — the physics, structure and analysis suites that compute on
                 downloaded coordinates. Run when the science changes.
- ``ui``       — the Qt suites that run on the offscreen platform. Run when
                 panels, controllers or menus change.
- ``render``   — the suites that need a real OpenGL context and judge the
                 screen in pixels. Run when the renderer or a drawing path
                 changes; they exist because everything else once passed
                 while nothing was drawn at all.
- ``records``  — the documentation, record and reproducibility guards: the
                 claims registry, the READMEs, the notebooks, the frozen
                 validation rounds. Run before anything is committed that
                 touches prose or a recorded number.
"""

from __future__ import annotations

__all__ = ["TIERS", "all_files", "files_for"]

TIERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "quick": (
        "sanity on every edit: imports, registries, resources, pure logic",
        (
            "test_architecture.py",
            "test_calibration.py",
            "test_cif_reader.py",
            "test_data_routes.py",
            "test_dead_code.py",
            "test_dome_idealised.py",
            "test_engineered.py",
            "test_export.py",
            "test_external.py",
            "test_external_md.py",
            "test_feasibility.py",
            "test_fetch_content.py",
            "test_harvest.py",
            "test_harvest_curation.py",
            "test_imports.py",
            "test_kinetics.py",
            "test_measurement_set.py",
            "test_membrane.py",
            "test_paired_feasibility.py",
            "test_parameters.py",
            "test_sequence_and_resources.py",
            "test_sequences.py",
            "test_superpose.py",
            "test_tiers.py",
            "test_tour.py",
            "test_validation.py",
            "test_variant_sets.py",
            "test_workflow.py",
        ),
    ),
    "science": (
        "physics, structure and analysis on real coordinates",
        (
            "test_allostery.py",
            "test_anm.py",
            "test_annotation_coverage.py",
            "test_assembly.py",
            "test_conduction_path.py",
            "test_conduction_verdict.py",
            "test_conservation.py",
            "test_cross_helices.py",
            "test_crosscheck.py",
            "test_crosscheck_methods.py",
            "test_design.py",
            "test_elastica.py",
            "test_electrostatics.py",
            "test_ensemble.py",
            "test_entities.py",
            "test_features.py",
            "test_fluctuations.py",
            "test_frame.py",
            "test_full_length.py",
            "test_fusion.py",
            "test_fusion_pose.py",
            "test_gating_area.py",
            "test_geometry.py",
            "test_guo2017.py",
            "test_homology.py",
            "test_homology_structure.py",
            "test_hybrid.py",
            "test_hydration.py",
            "test_hydropathy.py",
            "test_interactions.py",
            "test_labelling.py",
            "test_ligands.py",
            "test_liu2025.py",
            "test_martini.py",
            "test_measure.py",
            "test_model_error.py",
            "test_morph.py",
            "test_nanodomain.py",
            "test_paired_variant.py",
            "test_paralogue.py",
            "test_parameter_effect.py",
            "test_performance.py",
            "test_permeation.py",
            "test_planarity.py",
            "test_pockets.py",
            "test_pore.py",
            "test_pore_charge.py",
            "test_pore_regions.py",
            "test_prediction_confidence.py",
            "test_provenance_chain.py",
            "test_renumbering.py",
            "test_reproduce_young2023.py",
            "test_selectivity.py",
            "test_substitution.py",
            "test_uncertainty.py",
            "test_variant_impact.py",
            "test_variant_structures.py",
        ),
    ),
    "ui": (
        "the Qt suites, on the offscreen platform",
        (
            "test_hazards.py",
            "test_ui_analysis.py",
            "test_ui_companions.py",
            "test_ui_components.py",
            "test_ui_context_menu.py",
            "test_ui_display.py",
            "test_ui_dome_surface.py",
            "test_ui_feature_styles.py",
            "test_ui_fluctuation_color.py",
            "test_ui_fusion.py",
            "test_ui_hybrid.py",
            "test_ui_interactions.py",
            "test_ui_measure.py",
            "test_ui_membrane_views.py",
            "test_ui_morph.py",
            "test_ui_nanodomain.py",
            "test_ui_path.py",
            "test_ui_pockets.py",
            "test_ui_pore_surface.py",
            "test_ui_shell.py",
            "test_ui_status_bar.py",
            "test_ui_structure_panel.py",
            "test_ui_topology.py",
        ),
    ),
    "render": (
        "needs a real OpenGL context; judges the screen in pixels",
        (
            "test_ion_flux.py",
            "test_render_impostors.py",
            "test_ui_controls.py",
            "test_ui_morph_fusion.py",
            "test_ui_overlays_render.py",
        ),
    ),
    "records": (
        "documentation, frozen records and reproducibility guards",
        (
            "test_cold_clone_check.py",
            "test_conclusion.py",
            "test_entry_points.py",
            "test_methods_note.py",
            "test_not_preregistered_round64.py",
            "test_notebooks.py",
            "test_prediction_record.py",
            "test_predictor_retired.py",
            "test_published_interval.py",
            "test_readme.py",
            "test_reproduce_from_clone.py",
            "test_reproducibility.py",
            "test_roadmap.py",
            "test_validation_round22.py",
            "test_validation_round36.py",
            "test_validation_round41.py",
            "test_validation_round48.py",
        ),
    ),
}


def all_files() -> set[str]:
    """Every file the tier map accounts for."""
    return {name for _desc, files in TIERS.values() for name in files}


def files_for(suites) -> set[str]:
    """The files the named tiers select. Unknown names raise, loudly:
    a typo that silently selected nothing would report a green run."""
    unknown = [s for s in suites if s not in TIERS]
    if unknown:
        raise KeyError(
            f"unknown suite(s) {unknown}; choose from {sorted(TIERS)}")
    return {name for s in suites for name in TIERS[s][1]}
