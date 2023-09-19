(["id", "molecule", "organism", "seq", "taxid"] | @tsv),
(.Data |
    map(
        [
            (.["Epitope ID - IEDB IRI"] | split("/")[-1]),
            (
                [
                    .["Epitope - Molecule Parent"],
                    .["Related Object - Molecule Parent"]
                ] | map(select(. != "")) | first // ""
            ),
            (
                [
                    .["Epitope - Source Organism"],
                    .["Related Object - Source Organism"]
                ] | map(select(. != "")) | first // ""
            ),
            (.["Epitope - Name"] | split(" ")[0]),
            (
                [
                    .["Epitope - Source Organism IRI"],
                    .["Related Object - Source Organism IRI"]
                ] | map(select(contains("NCBITaxon_"))) | first // "" | split("_")[-1]
            )
        ]
    ) |
    map(
        select((.[-2] | match("^[ACDEFGHIKLMNPQRSTVWY]+$")) and (.[-1] != null)) | @tsv
    )[]
)