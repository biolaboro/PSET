from shiny import ui

from shared import *

app_ui = ui.page_fluid(
    ui.navset_card_pill(
        ui.nav_panel(
            "PSET",
            ui.navset_card_pill(
                ui.nav_panel(
                    "input",
                    ui.layout_sidebar(
                        ui.sidebar(
                            ui.input_action_button("run_pset", "run"),
                            ui.accordion(
                                ui.accordion_panel(
                                    "input",
                                    ui.input_file("info", "assay"),
                                    ui.input_select("database", "database", choices=[""], multiple=True, selectize=True),
                                ),
                                ui.accordion_panel(
                                    "parameters",
                                    ui.input_switch("forceall", "force rerun of analysis", value=False),
                                    ui.input_switch("flank", "use flank mode", value=False),
                                    ui.input_numeric("context", "context", value=6, min=0),
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
                                ui.accordion_panel(
                                    "threading",
                                    ui.input_numeric("max_threads", label="max threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT),
                                    ui.input_numeric("lcl_threads", label="local threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT),
                                    ui.input_numeric("glc_threads", label="glocal threads", value=min(2, CPU_COUNT), min=0, max=CPU_COUNT),
                                ),
                            ),
                            width="400px"
                        ),
                        ui.navset_tab(
                            ui.nav_panel("assay", ui.card(ui.output_data_frame("assay_table"))),
                            ui.nav_panel("databases", ui.card(ui.output_data_frame("database_table")))
                        ),
                    ),
                ),
                ui.nav_panel(
                    "output",
                    ui.layout_sidebar(
                        ui.sidebar(
                            ui.input_action_button("report_save", "save"),
                            ui.accordion(
                                ui.accordion_panel(
                                    "data",
                                    ui.card(ui.input_action_button("report_load", "load")),
                                    ui.card(ui.output_data_frame("report_runs")),
                                ),
                                ui.accordion_panel(
                                    "dimensions",
                                    ui.card(
                                        ui.layout_columns(
                                            ui.input_numeric(f"report_width", "plot width (px)", value=800, min=0, step=100),
                                            ui.input_numeric(f"report_height", "plot height (px)", value=800, min=0, step=100),
                                        ),
                                    ),
                                ),
                                ui.accordion_panel(
                                    "components",
                                    ui.card(ui.input_switch(f"report_keyify", "use assay id key", value=False)),
                                    ui.card(
                                        ui.layout_columns(
                                            ui.input_checkbox_group(f"report_call", "calls", choices=CALLS, selected=CALLS[:-1], inline=False),
                                            ui.input_checkbox_group(f"report_comp", "components", choices=COMPONENTS, selected=COMPONENTS, inline=False)
                                        )
                                    ),
                                ),
                                ui.accordion_panel(
                                    "taxa",
                                    ui.card(ui.input_select(f"report_taxa", "", choices=[], multiple=True, size=20, width="100%"))
                                ),
                            ),
                            id="report_sidebar",
                            width="800px",
                        ),
                        ui.navset_tab(
                            ui.nav_panel(
                                "tables",
                                ui.card(
                                    ui.input_checkbox_group(
                                        "report_aggregate",
                                        "aggregate by assay id, including",
                                        choices=CONFUSION_AGGREGATE_KEYS,
                                        selected=CONFUSION_AGGREGATE_KEYS[1],
                                        inline=True
                                    )
                                ),
                                ui.card(ui.output_data_frame("report_confusion")),
                            ),
                            ui.nav_panel(
                                "plots",
                                ui.navset_tab(
                                    *(
                                        ui.nav_panel(
                                            ele,
                                            ui.card(
                                                ui.card(ui.input_action_button(f"report_plot_{ele}", "plot")),
                                                ui.card(ui.output_ui(f"report_{ele}_ui"), open=True)
                                            )
                                        ) for ele in ("heat", "muts")
                                    ),
                                    id="report_navset_plots"
                                )
                            )
                        ),
                    ),
                ),
            ),
        ),
        ui.nav_panel(
            "CUSTOM DB",
            ui.layout_sidebar(
                ui.sidebar(
                    "Input FASTA file and accession-taxon mapping. "
                    "The mapping file must not contain column names. "
                    "The first column of the mapping file is the accession and the second is the taxon identifier. ",
                    "If the mapping file is an Excel file, then the first sheet is assumed.",
                    ui.input_file("info_fasta", "DNA sequences (FASTA)"),
                    ui.input_file("info_taxon", "accession-taxon mapping"),
                    ui.input_text("db_title", "title"),
                    ui.input_action_button("run_build", "build"),
                    width="400px",
                ),
                ui.output_table("result_mapping"),
            ),
        ),
        ui.nav_panel(
            "DOWNLOAD",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_action_button("run_database_listing", "⓵ download listing"),
                    ui.input_select("select_db", "Select BLAST+ Database", choices=[], selected=""),
                    ui.input_action_button("run_database_download", "⓶ download database", enabled=False),
                    width="400px",
                ),
                ui.output_data_frame("remote_databases"),
            ),
        ),
        ui.nav_panel(
            "GENERATE",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_action_button("agen_run", "run"),
                    ui.accordion(
                        ui.accordion_panel("input", ui.input_file("agen_fasta", "FASTA")),
                        ui.accordion_panel(
                            "parameters",
                            ui.input_text(
                                "agen_cstr", "config override", placeholder="format='section:key=val,...', example: 'GLOBAL:PRIMER_NUM_RETURN=500'"
                            ),
                            ui.input_radio_buttons("agen_mode", "mode", choices=("PCR", "LAMP"), inline=True),
                            ui.input_switch("agen_loop", "optional loop"),
                            ui.input_numeric("agen_limit", "limit", value=10, min=1, max=1000),
                        ),
                        ui.accordion_panel(
                            "threading", ui.input_numeric("agen_threads", label="max threads", value=min(8, CPU_COUNT), min=0, max=CPU_COUNT)
                        ),
                    ),
                    width="400px",
                ),
                ui.output_table("agen_output"),
            ),
        ),
        ui.nav_panel(
            "TAXA",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.input_action_button("run_taxa", "run"),
                    ui.input_numeric("tax_id", "NCBI Taxonomy", value=666),
                    ui.input_radio_buttons("lineage_mode", "mode", choices=("ancestors", "descendants", "count"), inline=True),
                    ui.HTML('<hr style="color: purple;">'),
                    ui.HTML("<b>parameters for count mode</b>"),
                    ui.input_select("taxa_blastdb", "database", choices=[""], multiple=False, selectize=True),
                    ui.input_switch("near_neighbors", "include near neighbors", value=True),
                    ui.input_switch("missing_taxa", "include missing taxa", value=True),
                    width="400px",
                ),
                ui.output_data_frame("taxa_table"),
            ),
        ),
        id="tabs",
    )
)
