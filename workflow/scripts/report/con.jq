(
    .assay.id as $id |
    .hits |
    map(.calls | to_entries) |
    reduce .[][] as $ele ({}; .[$ele.key[:2]]+=($ele.value|length)) |
    [$id, .TP, .FN, .FP, .TN] |
    map(. // 0) |
    @tsv
)