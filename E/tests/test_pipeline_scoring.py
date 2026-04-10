from scoring import build_rankings


def test_pipeline_style_scoring_integration():
    """
    模拟 pipeline 最后拿到 d_results 和 e_evals 后，
    能正常生成 rankings。
    """
    d_results = [
        {
            "pair_index": 0,
            "forward": {"sequence": "ATGGAAGGATCATCAGGT"},
            "reverse": {"sequence": "CTAGCCTTGTAGGTCAAG"},
            "tm_diff": 0.7,
            "heterodimer_dG": -5.2,
            "product_size": 220,
            "target_amplicons": [{"subject_id": "XM_target"}],
            "off_targets": {
                "near": [],
                "far": [{"subject_id": "XM_far1"}],
            },
        }
    ]

    e_evals = [
        {
            "pair_index": 0,
            "target_expr_by_tissue": {
                "root": 18.0,
                "leaf": 2.0,
                "seed": 0.4,
            },
            "near_expr_by_tissue": {
                "root": 0.5,
                "leaf": 0.7,
                "seed": 0.0,
            },
            "far_expr_by_tissue": {
                "root": 0.3,
                "leaf": 0.5,
                "seed": 0.1,
            },
            "target_vs_near_ratio": {
                "root": 36.0,
                "leaf": 2.86,
                "seed": 4.0,
            },
            "target_vs_far_ratio": {
                "root": 60.0,
                "leaf": 4.0,
                "seed": 4.0,
            },
            "dominant_target_tissues": ["root"],
        }
    ]

    rankings = build_rankings(d_results, e_evals)

    assert isinstance(rankings, dict)
    assert len(rankings["thermo_ranking"]) == 1
    assert len(rankings["expression_ranking"]) == 3
    assert len(rankings["combined_ranking"]) == 3

    best_expr = rankings["expression_ranking"][0]
    assert "pair_index" in best_expr
    assert "source" in best_expr
    assert "expression_score" in best_expr

    best_combined = rankings["combined_ranking"][0]
    assert "combined_score" in best_combined
    assert "final_recommendation" in best_combined
    assert "warnings" in best_combined