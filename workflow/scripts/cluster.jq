(input_filename | split("/")) as $name |
.ref as $ref |
(.accs | length) as $nacc |
.coor as $coor |
.hits as $hits |
.seqs | to_entries |
map([
    $name[-3],
    ($name[-1] | split(".")[0]),
    $ref, $coor[.key][],
    $nacc,
    ($hits[.key] | add / length),
    (.value | length),
    if .value | length <= 5000 then .value else "..." end
] | @tsv)[]
