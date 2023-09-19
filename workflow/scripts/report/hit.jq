(
    ["acc", "call", "heat"] | @tsv
),
(
    .hits |
    map(
        [
            (.evals | map(.qsim) | length as $n | if $n > 0 then add / $n else 0 end),
            (.calls | to_entries[] | .key as $call | .value | map([., $call])[])
        ] | .[0] as $heat | .[1:] | map(. + [$heat])[]
    )[] | @tsv
)
