(["id", "definition", "targets"] | @tsv),
map(
    (.SEQUENCE_IDENTIFIER | sub("[^a-zA-Z0-9_. ]+"; "_")) as $id |
    to_entries |
    map(select(.key | contains("PSET_DEFINITION_")))[:($max | tonumber)] |
    map([$id + "_" + (.key | split("_") | last), (.value | ascii_upcase), $targets] | @tsv)[]
)[]
