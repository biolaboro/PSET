(
    .assay.id as $id |
    .db as $db |
    .hits |
    map(.calls | to_entries) |
    reduce .[][] as $ele ({}; .[$ele.key[:2]]+=($ele.value|length)) |
    [$db, $id, .TP, .FN, .FP, .TN] |
    map(. // 0) |
    @tsv
)