(select(.results != []) | .results[0].definition | [match("\\["; "g")] | length) as $nbrackets |
(
    [
        "id", "targets", "definition", "penalty", "key", "pos1", "pos2", "len", "Tm", "dG5", "dG3", "gc", "seq",
        if $nbrackets == 6 then "fip", "bip" else empty end
    ] | @tsv
),
(
    select(.results != []) |
    .results |
    map(
        .id as $id |
        .definition as $df |
        .fip as $fip |
        .bip as $bip |
        .penalty as $pn |
        (
            .sequences |
            to_entries |
            map(
                [
                    $id, $targets, $df, $pn, .key, .value.pos[], [.value[]][1:][],
                    if $nbrackets == 6 then $fip, $bip else empty end
                ] | @tsv
            )[]
        )
    )[]
)
