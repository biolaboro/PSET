(["id", "definition", "targets"] | @tsv), (select(.results != []) | .results | map([.id, .definition, $targets] | @tsv)[])
