import json
from pathlib import Path
from subprocess import check_call, check_output
from itertools import chain
from datetime import datetime

import pandas as pd
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from shiny import App, reactive, render, ui
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from taxa.taxa import ancestors, descendants, session_scope
from pset.assay import AlignType, parse_assays
from shared import *
from ui import *


def server(input, output, session):
    agen_expected_output = reactive.Value()
    max_comp_len = reactive.Value()
    df_conf = reactive.Value()
    df_heat_1 = reactive.Value()
    df_heat_2 = reactive.Value()
    df_muts = reactive.Value()
    df_local_db = reactive.Value()
    df_remote_db = reactive.Value()
    df_taxa = reactive.Value()
    df_sequences = reactive.Value()
    df_mapping = reactive.Value()
    plot_params = reactive.Value()

    def plot_heat(width, height, exts=("png", )):
        df = df_heat_2()
        with NamedTemporaryFile() as t1, NamedTemporaryFile(delete=False) as t2:
            df.to_csv(t1.name, sep="\t", index=False)
            dims = (str(ele * DPI_PLOT / DPI_APP) for ele in (width, height))
            cmd = ("Rscript", Path(__file__).parent.resolve() / "heat.R", t1.name, t2.name, *dims, "px", str(DPI_PLOT), *exts)
            return check_output(cmd, universal_newlines=True).split("\n")

    def plot_muts(width, height, exts=("png", )):
        df = df_muts()
        with NamedTemporaryFile() as t1, NamedTemporaryFile(delete=False) as t2:
            df.to_csv(t1.name, sep="\t", index=False)
            df.to_csv("~/Desktop/df2.tsv", sep="\t", index=False)
            dims = (str(ele * DPI_PLOT / DPI_APP) for ele in (width, height))
            cmd = ("Rscript", Path(__file__).parent.resolve() / "muts.R", t1.name, t2.name, *dims, "px", str(DPI_PLOT), *exts)
            return check_output(cmd, universal_newlines=True).split("\n")

    @reactive.Effect
    def _():
        # populate local databases
        choices = nucl_db_v5_choices()
        ui.update_select("database", choices=choices, session=session)
        ui.update_select("taxa_blastdb", choices=choices, session=session)
        data = blastdbcmd_info()
        df_local_db.set(data[(data["type"] == "Nucleotide") & (data["version"] == 5)].iloc[:, 2:-1])
        # plot parameters
        plot_params.set(
            dict(
                heat=dict(width=800, height=800, keyify=True, call=CALLS[:-1], comp=COMPONENTS, taxa=dict()),
                muts=dict(width=800, height=800, keyify=True, call=(CALLS[0], CALLS[3]), comp=COMPONENTS, taxa=dict()),
            )
        )

    @reactive.Effect
    def _():
        obj = plot_params()[input.report_navset_plots()]
        for ele in ("keyify", ):
            ui.update_switch(f"report_{ele}", value=obj[ele])
        for ele in ("width", "height"):
            ui.update_numeric(f"report_{ele}", value=obj[ele])
        for ele in ("call", "comp"):
            ui.update_checkbox_group(f"report_{ele}", selected=obj[ele])
        for ele in ("taxa", ):
            ui.update_select(f"report_{ele}", selected=obj[ele])

    @reactive.Effect
    def _():
        key = input.report_navset_plots()
        obj = plot_params()
        obj[key] = {ele: getattr(input, f"report_{ele}")() for ele in ("width", "height", "keyify", "call", "comp", "taxa")}
        plot_params.set(obj)

    @render.data_frame
    def report_runs():
        data = []
        for ele in sorted(PATH_RESULTS.rglob("pset/*/*/*/assay.json")):
            if ele.with_name("call.json").exists():
                with ele.open() as file:
                    obj = json.load(file)
                    data.append(dict(zip(
                        ("batch", "db", "id", "type", "targets"),
                        (ele.parts[-4], ele.parts[-3], ele.parts[-2], obj["type"], obj["targets"])
                    )))
        return render.DataTable(pd.DataFrame(data), selection_mode="rows", filters=True)

    @reactive.Effect
    def _():
        ui.update_action_button("report_load", disabled=not report_runs.cell_selection()["rows"])
        ui.update_action_button("report_save", disabled=not report_runs.cell_selection()["rows"])

    @reactive.Effect
    @reactive.event(input.report_save, ignore_init=True)
    def _():
        with ui.Progress(session=session) as prog:
            df = df_conf()
            name = datetime.now().isoformat().replace(":", "_")
            root = PATH_RESULTS / "app"
            os.makedirs(root, exist_ok=True)
            prog.set(message="output confusion matrix...")
            df.to_csv(root.joinpath(f"{name}_conf.tsv"), sep="\t")
            obj = plot_params()
            for key, plot_fn in zip(("heat", "muts"), (plot_heat, plot_muts)):
                width, height = obj[key]["width"], obj[key]["height"]
                prog.set(message=f"output {key} plot...")
                for path in map(Path, plot_fn(width, height, exts=("png", "pdf"))):
                    print(path, path.exists())
                    path.rename(root / path.with_stem(f"{name}_{key}").name)

    @reactive.Effect
    @reactive.event(input.report_load, ignore_init=True)
    def _():
        paths = [
            PATH_RESULTS / "pset" / Path(ele.batch, ele.db, ele.id, "call.json")
            for ele in report_runs.data().loc[list(report_runs.cell_selection()["rows"])].itertuples()
        ]
        df_hits = pd.DataFrame(chain.from_iterable(map(read_hits, paths)))
        if len(df_hits):
            with ui.Progress(session=session) as prog:
                prog.set(message="load assay hits...")
                df_hits = add_key(
                    df_hits.merge(
                        pd.concat(map(lambda x: read_func(x.with_name("lib.tax.tsv"), pd.read_table, header=None), paths)).set_axis(["acc", "tax"], axis=1).drop_duplicates(),
                        on="acc", how="left"
                    )
                )
                df_heat_1.set(df_hits)

                prog.set(message="record max component length...")
                max_comp_len.set(max(map(len, df_hits.astr)))

                prog.set(message="calculate subject taxonomy counts...")
                df_conf_temp = df_hits.drop(columns=["com", "psim", "astr"]).drop_duplicates()
                df_taxa = (
                    pd.DataFrame(
                        (
                            (db, *ele.split("\t")) for db in df_conf_temp.db.unique()
                            for ele in blastdb_taxidlist(BLAST_DIR / db / db, df_conf_temp.tax)
                        ),
                        columns=(columns := ["db", "tax"])
                    ).
                    astype(dict(tax="Int64")).
                    groupby(columns, observed=True).
                    size().
                    reset_index(name="count")
                )
                df_taxa = df_taxa if len(df_taxa) else df_conf_temp[columns].drop_duplicates()

                prog.set(message="query scientific names...")
                with session_scope(sessionmaker(bind=create_engine(DBURL))) as curs:
                    rows = []
                    for ele in df_taxa[columns].to_dict(orient="records"):
                        lineage = ancestors(curs, ele["tax"])
                        genus = next((ele for ele in lineage if ele["rank"] == "genus"), {})
                        species = next((ele for ele in lineage if ele["rank"] == "species"), {})
                        rows.append(dict(zip(
                            ("db", "tax", "sci", "rank", "genus", "genus_tax", "species"),
                            (ele["db"], ele["tax"], lineage[-1]["name_txt"], lineage[-1]["rank"], genus.get("name_txt", ""), genus.get("tax_id", ""), species.get("name_txt", "")
                             ))))
                    df_sci = pd.DataFrame(rows if rows else dict(db=pd.Series(dtype=str), tax=pd.Series(dtype=int), sci=pd.Series(dtype=str)))

                prog.set(message="set confusion matrix...")
                df_conf.set(df_conf_temp.merge(df_taxa, on=columns, how="left").merge(df_sci, on=columns))

                prog.set(message="set taxa dropdowns...")
                choices = defaultdict(dict)
                for ele in rows:
                    parens = f" ({val})" if (val := ele["genus_tax"]) else ""
                    key = (val + parens) if (val := ele["genus"]) else "?"
                    val = (val + " (" + str(ele["tax"]) + ")") if (val := ele["sci"]) else ele["tax"]
                    choices[key][ele["tax"]] = val
                selected = list(chain.from_iterable(choices.values()))
                ui.update_select(f"report_taxa", choices=choices, selected=selected, session=session)
                obj = plot_params()
                obj["heat"]["taxa"] = selected
                obj["muts"]["taxa"] = selected
                plot_params.set(obj)
        else:
            ui.notification_show("no results to plot!")

    @reactive.Effect
    @reactive.event(input.report_plot_heat)
    def _():
        if len(df := df_heat_1().copy()):
            df["assay"] = df.key if input.report_keyify() else df.id
            df = df[df.call.isin(calls := input.report_call()) & df.com.isin(coms := input.report_comp()) & df.tax.isin(set(map(int, input.report_taxa())))]
            df.com = pd.Categorical(df.com, categories=(ele for ele in coms if ele in df.com.unique()))
            df.call = pd.Categorical(df.call, categories=(ele for ele in calls if ele in df.call.unique()))
            df_heat_2.set(df)

    @reactive.Effect
    @reactive.event(input.report_plot_muts)
    def _():
        if len(df := df_heat_1().copy()):
            df = df[df.call.isin(calls := input.report_call()) & df.com.isin(coms := input.report_comp()) & df.tax.isin(set(map(int, input.report_taxa())))]
            if len(df):
                columns = ["id", "key", "com", "call", "astr"]
                records = (
                    df[~df.astr.str.contains(f"^[{AlignType.IDN.value}{AlignType.SIM.value}]+$", regex=True)].
                    groupby(columns, observed=True).
                    size().
                    reset_index(name="count").
                    to_dict(orient="records")
                )
                columns = ["id", "key", "com", "call", "pos", "mut", "count"]
                df = pd.DataFrame(
                    dict(zip(columns, (obj["id"], obj["key"], obj["com"], obj["call"], idx, ele, obj["count"])))
                    for obj in records
                    for idx, ele in enumerate(map(int, obj["astr"]), start=1)
                    if ele not in (AlignType.IDN.value, AlignType.SIM.value)
                )
                if len(df):
                    df = df.groupby(columns[:-1], observed=True).sum().reset_index()
                    muts = [None, *(ele.name for ele in AlignType)]
                    df.mut = pd.Categorical([muts[ele] for ele in df.mut], categories=(ele.name for ele in AlignType if ele.name not in ("IDN", "SIM")))
                    df.com = pd.Categorical(df.com, categories=(ele for ele in coms if ele in df.com.unique()))
                    df.call = pd.Categorical(df.call, categories=(ele for ele in calls if ele in df.call.unique()))
                    df = df.groupby(df.columns.tolist()[:-1], observed=True).sum().reset_index()
                df["assay"] = df.key if input.report_keyify() else df.id
                df_muts.set(df)

    @render.ui
    @reactive.event(input.report_plot_heat)
    def report_heat_ui():
        obj = plot_params()[key := input.report_navset_plots()]
        return ui.TagList(ui.output_image(f"report_{key}", width=obj["width"], height=obj["height"]))

    @render.ui
    @reactive.event(input.report_plot_muts)
    def report_muts_ui():
        obj = plot_params()[key := input.report_navset_plots()]
        return ui.TagList(ui.output_image(f"report_{key}", width=obj["width"], height=obj["height"]))

    @render.data_frame
    def report_confusion():
        columns = ["id", "key", *input.report_aggregate(), "call"]
        df = df_conf()
        return (
            render.DataTable(
                df.
                groupby(columns, group_keys=False, as_index=False, dropna=False, observed=True).
                size().
                pivot(index=columns[:-1], columns=columns[-1], values="size").reset_index().rename_axis(None, axis=1).
                merge(df[["id", "key"]].drop_duplicates(), on=["id", "key"], how="left"),
                width="100%"
            )
        )

    @render.image(delete_file=True)
    @reactive.event(input.report_plot_heat)
    def report_heat():
        obj = plot_params()[input.report_navset_plots()]
        width, height = obj["width"], obj["height"]
        return dict(src=plot_heat(width, height)[0], width=width, height=height)

    @render.image(delete_file=True)
    @reactive.event(input.report_plot_muts)
    def report_muts():
        obj = plot_params()[input.report_navset_plots()]
        width, height = obj["width"], obj["height"]
        return dict(src=plot_muts(width, height)[0], width=width, height=height)

    @output(suspend_when_hidden=False)
    @render.data_frame
    def assay_table():
        info = input.info()
        if not info:
            return
        with open(info[0]["datapath"]) as file:
            data = []
            nprob = 0
            for rec, row in parse_assays(file):
                nprob += (prob := rec is None)
                data.append(
                    dict(ok=False, id=row.get("id", "?"), target=row.get("targets", "?"), definition=row.get("definition", "?"))
                    if prob else
                    dict(ok=True, id=rec.id, target=";".join(map(str, rec.targets)), definition=rec.definition)
                )
            ui.update_action_button("run_pset", label="run", disabled=bool(nprob))
            ui.notification_show(
                f"Found {nprob} assay problem{'s' * (nprob > 1)}, check the 'ok' column..."
                if nprob else
                f"Loaded {len(data)} assay{'s' * (len(data) > 1)}!"
            )
            return render.DataTable(pd.DataFrame(data), width="100%")

    @output(suspend_when_hidden=False)
    @render.data_frame
    def database_table():
        return render.DataTable(df_local_db(), width="100%")

    @reactive.Effect
    @reactive.event(input.run_pset)
    def run_pset_workflow():
        info = input.info()
        if not info:
            return
        label = Path(info[0]["name"]).stem
        path_out = PATH_RESULTS / "pset" / label
        os.makedirs(path_out, exist_ok=True)
        for _, db in enumerate(input.database(), start=1):
            path_config = path_out.joinpath("config.json")
            with path_config.open("w") as file:
                json.dump(
                    dict(
                        file=info[0]["datapath"],
                        db=db,
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
            cmd = (
                "snakemake",
                "--forceall" * (input.forceall()),
                "--rerun-incomplete",
                "--cores",
                str(input.max_threads()),
                "--set-threads",
                f"local={input.lcl_threads()}",
                f"glocal={input.glc_threads()}",
                "--configfile",
                str(path_config),
                "--",
                "target_tsv",
            )
            cmd = list(filter(len, cmd))
            monitor_snakemake(session, cmd, db)

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
        with session_scope(sessionmaker(bind=create_engine(DBURL))) as curs:
            df_right = pd.DataFrame([next(ancestors(curs, ele)) for ele in df["taxon"].unique()])
            df_right.rename(columns=dict(tax_id="taxon", parent_tax_id="parent_taxon"), inplace=True)
            df_right["parent_taxon"] = df_right["parent_taxon"].astype(str)
            df_right["taxon"] = df_right["taxon"].astype(str)
            df = df.merge(df_right, on="taxon", how="left")
            df_taxa.set(df)

    # @reactive.Effect
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
    @render.data_frame
    @reactive.event(input.run_taxa)
    def taxa_table():
        with ui.Progress(session=session) as prog:
            prog.set(message=f"calculating {(mode := input.lineage_mode())}...")
            with session_scope(sessionmaker(bind=create_engine(DBURL))) as curs:
                if mode == "ancestors":
                    return pd.DataFrame(ancestors(curs, input.tax_id()))
                elif mode == "descendants":
                    return pd.DataFrame(descendants(curs, input.tax_id()))
                elif mode == "count":
                    qry = "count >= 0" if input.missing_taxa() else "count > 0"
                    if len(df := pd.DataFrame(count_nntaxa_in_blastdb(curs, input.taxa_blastdb(), input.tax_id(), input.near_neighbors()).values())):
                        return render.DataTable(df.drop(columns=["id"]).query(qry).sort_values(["count"], ascending=False), width="100%")
                else:
                    raise ValueError("incorrect mode for TAXA page...")

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
        path_out = PATH_RESULTS / "agen" / label
        agen_expected_output.set(
            [(path_out / ele.id / input.agen_mode()).with_suffix(".tsv") for ele in SeqIO.parse(info[0]["datapath"], "fasta")]
        )
        os.makedirs(path_out, exist_ok=True)
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
