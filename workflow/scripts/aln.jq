(
    [
        "db", "id",
        "primer1.aln", "primer1.sub", "primer1.ind",
        "probe.aln", "probe.sub", "probe.ind",
        "primer2.aln", "primer2.sub", "primer2.ind",
        "TP", "TN", "FP", "FN"
    ] | @tsv
),
.db as $db |
.assay.id as $id |
.hits |
    map(
        [
            $db,
            $id,
            "\"qry:\(.evals.primer1.qaln)|sbj:\(.evals.primer1.saln)\"",
            .evals.primer1.atypes.TRV + .evals.primer1.atypes.TRS,
            .evals.primer1.atypes.INS + .evals.primer1.atypes.DEL,
            "\"qry:\(.evals.probe.qaln)|sbj:\(.evals.probe.saln)\"",
            .evals.probe.atypes.TRV + .evals.probe.atypes.TRS,
            .evals.probe.atypes.INS + .evals.probe.atypes.DEL,
            "\"qry:\(.evals.primer2.qaln)|sbj:\(.evals.primer2.saln)\"",
            .evals.primer2.atypes.TRV + .evals.primer2.atypes.TRS,
            .evals.primer2.atypes.INS + .evals.primer2.atypes.DEL,
            (
                [
                    .calls.TP + .calls.TPN,
                    .calls.TN + .calls.TNN,
                    .calls.FP + .calls.FPN,
                    .calls.FN + .calls.FNN
                ] | map(length)[]
            )
        ] | map(. // 0)
    )[] | @tsv
