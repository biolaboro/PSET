{
  "results": map(
        (.args.conf | match("(\\w+)\\.json$").captures[0].string | gsub("_"; "-")) as $conf |
        select(.results != []) |
        .results |
        group_by(.id | match("_(\\d+-\\d+)_").captures[0].string) |
        map(.[:$n_per_region][] | .id += "_\($conf)")[]
    ) |
    sort_by(.penalty)[:$n_total]
}