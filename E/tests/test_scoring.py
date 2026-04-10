from scoring import (
    build_pair_profile,
    build_pair_source_profiles,
    build_rankings,
)


def test_build_pair_profile_basic():
    d_result = {
        "pair_index": 1,
        "forward": {"sequence": "ATGGAAGGATCATCAGGT"},
        "reverse": {"sequence": "CTAGCCTTGTAGGTCAAG"},
        "tm_diff": 0.8,
        "heterodimer_dG": -5.0,
        "product_size": 200,
        "target_amplicons": [{"subject_id": "XM_001"}],
        "off_targets": {
            "near": [],
            "far": [{"subject_id": "XM_far1"}],
        },
    }

    profile = build_pair_profile(d_result)

    assert profile["pair_index"] == 1
    assert profile["forward"] == "ATGGAAGGATCATCAGGT"
    assert profile["reverse"] == "CTAGCCTTGTAGGTCAAG"

    assert "thermo" in profile
    assert "specificity" in profile

    assert profile["thermo"]["label"] in {"excellent", "good", "usable", "warning"}
    assert profile["specificity"]["label"] in {
        "clean",
        "far_off_target_warning",
        "near_off_target_warning",
        "failed_target",
    }


def test_build_pair_source_profiles_basic():
    d_result = {
        "pair_index": 2,
    }

    e_eval = {
        "pair_index": 2,
        "target_expr_by_tissue": {
            "root": 20.0,
            "leaf": 0.5,
        },
        "near_expr_by_tissue": {
            "root": 1.0,
            "leaf": 0.2,
        },
        "far_expr_by_tissue": {
            "root": 0.5,
            "leaf": 0.1,
        },
        "target_vs_near_ratio": {
            "root": 20.0,
            "leaf": 2.5,
        },
        "target_vs_far_ratio": {
            "root": 40.0,
            "leaf": 5.0,
        },
        "dominant_target_tissues": ["root"],
    }

    profiles = build_pair_source_profiles(d_result, e_eval)

    assert len(profiles) == 2

    root_profile = [p for p in profiles if p["source"] == "root"][0]
    leaf_profile = [p for p in profiles if p["source"] == "leaf"][0]

    assert root_profile["pair_index"] == 2
    assert root_profile["expression_score"] >= leaf_profile["expression_score"]
    assert "source_tags" in root_profile
    assert "source_warnings" in root_profile


def test_low_expression_tissue_is_kept_with_warning():
    d_result = {
        "pair_index": 3,
    }

    e_eval = {
        "pair_index": 3,
        "target_expr_by_tissue": {
            "root": 12.0,
            "seed": 0.2,
        },
        "near_expr_by_tissue": {
            "root": 0.3,
            "seed": 0.0,
        },
        "far_expr_by_tissue": {
            "root": 0.1,
            "seed": 0.0,
        },
        "target_vs_near_ratio": {
            "root": 40.0,
            "seed": 10.0,
        },
        "target_vs_far_ratio": {
            "root": 120.0,
            "seed": 10.0,
        },
        "dominant_target_tissues": [],
    }

    profiles = build_pair_source_profiles(d_result, e_eval)
    seed_profile = [p for p in profiles if p["source"] == "seed"][0]

    assert seed_profile["source"] == "seed"
    assert "very_low_expression" in seed_profile["source_tags"]
    assert any("very low" in w for w in seed_profile["source_warnings"])


def test_build_rankings_basic():
    d_results = [
        {
            "pair_index": 1,
            "forward": {"sequence": "AAA"},
            "reverse": {"sequence": "TTT"},
            "tm_diff": 0.5,
            "heterodimer_dG": -5.0,
            "product_size": 180,
            "target_amplicons": [{"subject_id": "XM_target"}],
            "off_targets": {"near": [], "far": []},
        },
        {
            "pair_index": 2,
            "forward": {"sequence": "CCC"},
            "reverse": {"sequence": "GGG"},
            "tm_diff": 2.5,
            "heterodimer_dG": -8.0,
            "product_size": 250,
            "target_amplicons": [{"subject_id": "XM_target2"}],
            "off_targets": {"near": [{"subject_id": "XM_near"}], "far": []},
        },
    ]

    e_evals = [
        {
            "pair_index": 1,
            "target_expr_by_tissue": {"root": 15.0, "leaf": 3.0},
            "near_expr_by_tissue": {"root": 0.2, "leaf": 0.5},
            "far_expr_by_tissue": {"root": 0.1, "leaf": 0.2},
            "target_vs_near_ratio": {"root": 75.0, "leaf": 6.0},
            "target_vs_far_ratio": {"root": 150.0, "leaf": 15.0},
            "dominant_target_tissues": ["root"],
        },
        {
            "pair_index": 2,
            "target_expr_by_tissue": {"root": 10.0, "leaf": 4.0},
            "near_expr_by_tissue": {"root": 5.0, "leaf": 1.0},
            "far_expr_by_tissue": {"root": 0.5, "leaf": 0.3},
            "target_vs_near_ratio": {"root": 2.0, "leaf": 4.0},
            "target_vs_far_ratio": {"root": 20.0, "leaf": 13.3},
            "dominant_target_tissues": ["root"],
        },
    ]

    rankings = build_rankings(d_results, e_evals)

    assert "thermo_ranking" in rankings
    assert "expression_ranking" in rankings
    assert "combined_ranking" in rankings

    assert len(rankings["thermo_ranking"]) == 2
    assert len(rankings["expression_ranking"]) == 4
    assert len(rankings["combined_ranking"]) == 4

    assert rankings["expression_ranking"][0]["expression_score"] >= rankings["expression_ranking"][-1]["expression_score"]
    assert rankings["combined_ranking"][0]["combined_score"] >= rankings["combined_ranking"][-1]["combined_score"]