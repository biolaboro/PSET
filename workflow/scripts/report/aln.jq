.assay.id as $id |
(if .hits[0] | length > 0 then .hits[0].evals | keys_unsorted else empty end) as $com |
(["id", ($com | map("\(.).aln", "\(.).sub", "\(.).ind")[]), "TP", "TN", "FP", "FN"] | @tsv),
(
    .hits |
    map(
        .evals as $evals |
        .calls as $calls |
        [
            $id,
            (
                $com | map(
                    "\"qry:\($evals[.].qaln)|sbj:\($evals[.].saln)\"",
                    $evals[.].atypes.TRV + $evals[.].atypes.TRS,
                    $evals[.].atypes.INS + $evals[.].atypes.DEL
                )[]
            ),
            (
                [
                    $calls.TP + $calls.TPN,
                    $calls.TN + $calls.TNN,
                    $calls.FP + $calls.FPN,
                    $calls.FN + $calls.FNN
                ] | map(length)[]
            )
        ] | map(.//0) | @tsv
    )[]
)
