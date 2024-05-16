import json
import os
import sqlite3
from pathlib import Path
from subprocess import PIPE, Popen, check_call, check_output

import pandas as pd
from Bio import SeqIO
from shiny import App
from shiny import experimental as ex
from shiny import reactive, render, ui

from pset.assay import parse_assays

BLAST_DIR = Path("resources") / "blast"
PATH_TAXA = "resources/taxa/taxa.db"
DBURL = f"sqlite:///{PATH_TAXA}"
RESULTS = Path("results")
CPU_COUNT = os.cpu_count()


def blastdbcmd_info():
    fields = "\t".join(("%f", "%p", "%t", "%d", "%l", "%n", "%U", "%v"))
    cmd = ("blastdbcmd", "-recursive", "-remove_redundant_dbs", "-list", BLAST_DIR, "-list_outfmt", fields)
    with Popen(cmd, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdout as file:
            data = pd.read_table(file, names=("path", "type", "title", "date", "bases", "sequences", "bytes", "version"))
    return data


def nucl_db_v5_choices():
    data = blastdbcmd_info()
    data = data[(data["type"] == "Nucleotide") & (data["version"] == 5)]
    return dict(zip(data["path"], data["title"]))


def dict_factory(cursor, row):
    fields = [column[0] for column in cursor.description]
    return {key: value for key, value in zip(fields, row)}


def ancestors(curs, tax_id):
    ROOT = 1
    curs.execute("SELECT new_tax_id FROM tax_merged WHERE old_tax_id = ?;", (tax_id,))
    row = curs.fetchone()
    tax_id = row["new_tax_id"] if row else tax_id
    query = """
        SELECT
            tax_node.tax_id, tax_node.parent_tax_id, tax_node.rank,
            tax_name.name_txt, tax_name.unique_name, tax_name.name_class
        FROM tax_node
        LEFT JOIN tax_name ON tax_node.tax_id == tax_name.tax_id
        WHERE
            tax_node.tax_id == ? AND
            tax_name.name_class == 'scientific name'
        ;
    """
    while tax_id != ROOT:
        curs.execute(query, (tax_id,))
        row = curs.fetchone()
        yield row
        tax_id = row["parent_tax_id"]


def monitor_snakemake(session, cmd):
    print(cmd)
    with ui.Progress(min=0, max=100, session=session) as prog:
        prog.set(message="running...", detail=" ".join(cmd))
        with Popen(cmd, universal_newlines=True, bufsize=1, stderr=PIPE) as proc:
            with proc.stderr as file:
                for line in file:
                    line = line.strip()
                    print(line)
                    if line.endswith("%) done"):
                        prog.set(value=float(line.split(" ")[-2][1:-2]), message=line)


app_ui = ui.page_fluid(
    ui.navset_pill_card(
        ui.nav(
            "PSET",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_action_button("run_pset", "run"),
                    ui.hr(),
                    ex.ui.accordion(
                        ex.ui.accordion_panel(
                            "input",
                            ui.input_file("info", "assay"),
                            ui.input_select("database", "database", choices=[""], selected=""),
                        ),
                        ex.ui.accordion_panel(
                            "parameters",
                            ui.input_numeric("context", "context", value=6, min=0),
                            ui.input_switch("flank", "use flank mode", value=False),
                            ui.input_slider("simlcl", "local threshold", min=0, max=1, value=0.85),
                            ui.input_slider("simglc", "glocal threshold", min=0, max=1, value=0.90),
                            ui.input_text_area("confb", "BLAST+ config", value="-task=blastn -num_alignments=10000 -max_hsps=1 -subject_besthit"),
                            ui.input_text_area("confg", "glsearch36 config", value="-E 10000"),
                            ui.input_slider("dFR", "min/max forward/reverse primer distance", min=1, max=10000, value=(1, 1000)),
                            ui.input_slider("dF3F2", "min/max F3/F2 primer distance", min=1, max=1000, value=(20, 80)),
                            ui.input_slider("dF2F1c", "min/max F2/F1c primer distance", min=1, max=1000, value=(20, 80)),
                            ui.input_slider("dF1cB1c", "min/max F1c/B1c primer distance", min=1, max=1000, value=(1, 100)),
                            ui.input_text_area("xtaxa", "exclude taxa", value="81077"),
                        ),
                        ex.ui.accordion_panel(
                            "threading",
                            ui.input_numeric("max_threads", label="max threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT),
                            ui.input_numeric("lcl_threads", label="local threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT),
                            ui.input_numeric("glc_threads", label="glocal threads", value=min(2, CPU_COUNT), min=0, max=CPU_COUNT),
                        ),
                    ),
                ),
                ui.panel_main(
                    ui.navset_pill_card(
                        ui.nav(
                            "input",
                            ui.navset_tab(ui.nav("assay", ui.output_table("assay_table")), ui.nav("databases", ui.output_table("database_table"))),
                        ),
                        ui.nav(
                            "output",
                            ui.navset_tab(
                                ui.nav("confusion", ui.output_table("result_table")),
                                ui.nav("report", ui.output_ui("result_report")),
                            ),
                        ),
                    )
                ),
            ),
        ),
        ui.nav(
            "CUSTOM DB",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    "Input FASTA file and accession-taxon mapping. "
                    "The mapping file must not contain column names. "
                    "The first column of the mapping file is the accession and the second is the taxon identifier. ",
                    "If the mapping file is an Excel file, then the first sheet is assumed.",
                    ui.hr(),
                    ui.input_file("info_fasta", "DNA sequences (FASTA)"),
                    ui.input_file("info_taxon", "accession-taxon mapping"),
                    ui.input_text("db_title", "title"),
                    ui.input_action_button("run_build", "build"),
                ),
                ui.panel_main(
                    ui.output_table("result_mapping"),
                ),
            ),
        ),
        ui.nav(
            "DOWNLOAD",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_action_button("run_database_listing", "⓵ download listing"),
                    ui.hr(),
                    ui.input_select("select_db", "Select BLAST+ Database", choices=[], selected=""),
                    ui.input_action_button("run_database_download", "⓶ download database", enabled=False),
                ),
                ui.panel_main(ui.output_table("remote_databases")),
            ),
        ),
        ui.nav(
            "GENERATE",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_action_button("agen_run", "run"),
                    ui.hr(),
                    ex.ui.accordion(
                        ex.ui.accordion_panel("input", ui.input_file("agen_fasta", "FASTA")),
                        ex.ui.accordion_panel(
                            "parameters",
                            ui.input_text(
                                "agen_cstr", "config override", placeholder="format='section:key=val,...', example: 'GLOBAL:PRIMER_NUM_RETURN=500'"
                            ),
                            ui.input_radio_buttons("agen_mode", "mode", choices=("PCR", "LAMP"), inline=True),
                            ui.input_switch("agen_loop", "optional loop"),
                            ui.input_numeric("agen_limit", "limit", value=10, min=1, max=1000),
                        ),
                        ex.ui.accordion_panel(
                            "threading", ui.input_numeric("agen_threads", label="max threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT)
                        ),
                    ),
                ),
                ui.panel_main(ui.output_table("agen_output")),
            ),
        ),
        ui.nav(
            "TAXA",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_action_button("run_taxa", "run"),
                    ui.hr(),
                    ui.input_numeric("tax_id", "NCBI Taxonomy", value=666),
                ),
                ui.panel_main(ui.output_table("taxa_table")),
            ),
        ),
        id="tabs",
    )
)


def server(input, output, session):
    root = reactive.Value()
    agen_expected_output = reactive.Value()
    df_local_db = reactive.Value()
    df_remote_db = reactive.Value()
    df_taxa = reactive.Value()
    df_sequences = reactive.Value()
    df_mapping = reactive.Value()

    @reactive.Effect
    def _():
        # populate local databases
        ui.update_select("database", choices=nucl_db_v5_choices(), session=session)
        data = blastdbcmd_info()
        df_local_db.set(data[(data["type"] == "Nucleotide") & (data["version"] == 5)].iloc[:, 2:-1])

    @output(suspend_when_hidden=False)
    @render.table
    def assay_table():
        info = input.info()
        if not info:
            return
        with open(info[0]["datapath"]) as file:
            data = []
            for record in parse_assays(file):
                data.append(dict(id=record.id, target=record.targets, definition=record.definition))
            return pd.DataFrame(data)

    @output(suspend_when_hidden=False)
    @render.table
    def database_table():
        return df_local_db()

    @reactive.Effect
    @reactive.event(input.run_pset)
    def run_pset_workflow():
        info = input.info()
        if not info:
            return
        label = Path(info[0]["name"]).stem
        path_out = RESULTS / "pset" / label
        os.makedirs(path_out, exist_ok=True)
        print(path_out)
        path_config = path_out.joinpath("config.json")
        with path_config.open("w") as file:
            json.dump(
                dict(
                    file=info[0]["datapath"],
                    db=input.database(),
                    out=str(path_out),
                    flank=input.flank(),
                    context=f"{input.context()},{input.context()}",
                    dburl=DBURL,
                    confb=input.confb(),
                    confg=input.confg(),
                    simlcl=input.simlcl(),
                    simglc=input.simglc(),
                    dFR=",".join(map(str, input.dFR())),
                    dF3F2=",".join(map(str, input.dF3F2())),
                    dF2F1c=",".join(map(str, input.dF2F1c())),
                    dF1cB1c=",".join(map(str, input.dF1cB1c())),
                    xtaxa=input.xtaxa(),
                ),
                fp=file,
                indent=True,
            )
        for ele in ("tsv", "report"):
            cmd = (
                "snakemake",
                "--cores",
                str(input.max_threads()),
                "--set-threads",
                f"local={input.lcl_threads()}",
                f"glocal={input.glc_threads()}",
                "--configfile",
                str(path_config),
                "--",
                f"target_{ele}",
            )
            cmd = list(filter(len, cmd))
            monitor_snakemake(session, cmd)
        path_root = path_out.joinpath(Path(input.database()).stem)
        if path_root.joinpath("con.tsv").exists():
            root.set(path_root)

    @output(suspend_when_hidden=False)
    @render.table
    def result_table():
        return pd.read_table(root().joinpath("con.tsv"), delimiter="\t")

    @output(suspend_when_hidden=False)
    @render.ui
    def result_report():
        with ui.Progress() as prog:
            with root().joinpath("report.html").open() as file:
                return ui.HTML(file.read())

    @reactive.Effect
    @reactive.event(input.info_fasta)
    def load_fasta():
        info = input.info_fasta()
        if not info:
            return
        df_sequences.set(pd.DataFrame(dict(accession=[ele.id for ele in SeqIO.parse(info[0]["datapath"], "fasta")])))

    @reactive.Effect
    @reactive.event(input.info_taxon)
    def load_mapping():
        info = input.info_taxon()
        if not info:
            return
        path = Path(info[0]["datapath"])
        df = (
            pd.read_excel(path, names=["accession", "taxon"], dtype=str)
            if path.suffix == ".xlsx" else
            pd.read_table(path, names=["accession", "taxon"], dtype=str, sep="\s+")
        )
        pd.set_option("display.max_rows", None)
        conn = sqlite3.connect(PATH_TAXA)
        curs = conn.cursor()
        curs.row_factory = dict_factory
        df_right = pd.DataFrame([next(ancestors(curs, ele)) for ele in df["taxon"].unique()])
        df_right.rename(columns=dict(tax_id="taxon", parent_tax_id="parent_taxon"), inplace=True)
        df_right["parent_taxon"] = df_right["parent_taxon"].astype(str)
        df_right["taxon"] = df_right["taxon"].astype(str)
        df = df.merge(df_right, on="taxon", how="left")
        curs.close()
        conn.close()
        df_taxa.set(df)

    @reactive.Effect
    @output
    @render.table
    def result_mapping():
        df = df_sequences().merge(df_taxa(), on="accession", how="left")
        df_mapping.set(df)
        return df

    @reactive.Effect
    @reactive.event(input.run_build)
    def run_build():
        title = input.db_title()
        path_dir = BLAST_DIR.joinpath(title)
        out = str(path_dir.joinpath(title))
        df = df_mapping()
        if df["taxon"].isnull().any():
            modal = ui.modal(
                "Accession(s) missing taxon identifiers!",
                title="Error",
                easy_close=True,
                footer=None,
            )
        else:
            os.makedirs(path_dir, exist_ok=True)
            path_map = path_dir.joinpath(title).with_suffix(".ssv")
            with path_map.open("w") as file:
                df.to_csv(file, sep=" ", columns=["accession", "taxon"], header=False, index=False)
            cmd = (
                "makeblastdb",
                "-in", input.info_fasta()[0]["datapath"],
                "-input_type", "fasta",
                "-dbtype", "nucl",
                "-title", title,
                "-parse_seqids",
                "-hash_index",
                "-out", str(path_dir.joinpath(title)),
                "-blastdb_version", "5",
                "-logfile", str(path_dir.joinpath(title).with_suffix(".log")),
                "-taxid_map", str(path_map)
            )
            with ui.Progress(session=session) as prog:
                prog.set(message=f"building {title} BLAST+ database...")
                check_call(cmd)
                print(" ".join(cmd))
                prog.set(message="getting info...")

            modal = ui.modal(
                f"Database '{title}' built!",
                title="Complete",
                easy_close=True,
                footer=None,
            )
            ui.update_select("database", choices=nucl_db_v5_choices(), session=session)
            data = blastdbcmd_info()
            df_local_db.set(data[(data["type"] == "Nucleotide") & (data["version"] == 5)].iloc[:, 2:-1])
        ui.modal_show(modal)

    @output(suspend_when_hidden=False)
    @render.table
    @reactive.event(input.run_taxa)
    def taxa_table():
        tax_id = input.tax_id()
        conn = sqlite3.connect(PATH_TAXA)
        curs = conn.cursor()
        curs.row_factory = dict_factory
        result = pd.DataFrame(ancestors(curs, tax_id))
        curs.close()
        conn.close()
        return result

    @reactive.Effect
    @reactive.event(input.run_database_listing)
    def run_download_listing():
        with ui.Progress(session=session) as prog:
            prog.set(message="downloading NCBI BLAST+ database listing...")
            cmd = ("update_blastdb.pl", "--showall", "tsv")
            data = check_output(cmd, universal_newlines=True).split("\n")[1:]  # remove "Connected to NCBI"
            data = sorted(line.strip().split("\t") for line in data)
            data = (line for line in data if line)
            # name, description, size in gigabytes, date of last update (YYYY-MM-DD format)
            data = pd.DataFrame(data, columns=("name", "description", "size_gb", "date"))
            df_remote_db.set(data)

    @output(suspend_when_hidden=False)
    @render.table
    def remote_databases():
        data = df_remote_db()
        choices = {ele: ele for ele in data["name"]}
        selected = "nt" if "nt" in data["name"].values else data["name"][0]
        ui.update_select("select_db", choices=choices, selected=selected, session=session)
        return data

    @reactive.Effect
    @reactive.event(input.run_database_download)
    def download():
        if not input.select_db():
            return
        cmd = (
            "snakemake",
            "--printshellcmds",
            "--cores",
            str(input.max_threads()),
            "--set-threads",
            f"download={str(input.max_threads())}",
            "-s",
            "workflow/rules/setup.smk",
            "--config",
            f"db={input.select_db()}",
            "--",
            "download",
        )
        monitor_snakemake(session, cmd)
        ui.update_select("database", choices=nucl_db_v5_choices(), session=session)

    @reactive.Effect
    @reactive.event(input.agen_run)
    def run_agen_workflow():
        info = input.agen_fasta()
        if not info:
            return
        label = Path(info[0]["name"]).stem
        path_out = RESULTS / "agen" / label
        agen_expected_output.set(
            [(path_out / ele.id /input.agen_mode()).with_suffix(".tsv")  for ele in SeqIO.parse(info[0]["datapath"], "fasta")]
        )
        os.makedirs(path_out, exist_ok=True)
        print(path_out)
        path_config = path_out.joinpath("config.json")
        with path_config.open("w") as file:
            obj = dict(
                file=info[0]["datapath"],
                out=str(path_out),
                mode=input.agen_mode(),
                limit=input.agen_limit(),
            )
            if input.agen_cstr():
                obj["cstr"] = input.agen_cstr()
            if input.agen_loop():
                obj["optional_loop"] = input.agen_loop()
            json.dump(obj, fp=file, indent=True)
        cmd = (
            "snakemake",
            "-s",
            "./workflow/rules/agen.smk",
            "--cores",
            str(input.agen_threads()),
            "--set-threads",
            f"agen={input.agen_threads()}",
            "--configfile",
            str(path_config),
            "--",
            "target",
        )
        print(cmd)
        monitor_snakemake(session, cmd)

    @output(suspend_when_hidden=False)
    @render.table
    def agen_output():
        dfs = []
        for ele in agen_expected_output():
            temp = pd.read_csv(ele, sep="\t")
            if not temp.empty:
                temp["subject"] = ele.parent.name
                dfs.append(temp)
        return pd.concat(dfs)


app = App(app_ui, server)
