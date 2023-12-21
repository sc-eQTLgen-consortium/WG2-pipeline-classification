#!/usr/bin/env python


rule split:
    input:
        file = config["inputs"]["query_rds"] if config["inputs"]["query_rds"] is not None else []
    output:
        files = temp(expand(config["outputs"]["output_dir"] + "split/split_object_{pool}.RDS", pool=POOLS))
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["preprocess"]["split_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["preprocess"]["split_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["preprocess"]["split_time"]]
    threads: config["preprocess"]["split_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/split.R",
        batch = "Pool",
        out = "split_object",
        path = config["outputs"]["output_dir"] + "split/"
    log: config["outputs"]["output_dir"] + "log/split.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {input.file} \
            --batch {params.batch} \
            --out {params.out} \
            --path {params.path}
        """


rule make_seurat:
    input:
        poolsheet = config["inputs"]["poolsheet_path"]
    output:
        files = temp(config["outputs"]["output_dir"] + "make_seurat/query_{pool}.RDS")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["preprocess"]["make_seurat_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["preprocess"]["make_seurat_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["preprocess"]["make_seurat_time"]]
    threads: config["preprocess"]["make_seurat_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/make_seurat.R",
        batch = lambda wildcards: "--batch " + wildcards.pool if config["settings"]["split"] else "",
        out = "query",
        path = config["outputs"]["output_dir"] + "make_seurat/"
    log: config["outputs"]["output_dir"] + "log/make_seurat.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            {params.batch} \
            --out {params.out} \
            --path {params.path}
        """


def get_map_input(wildcards):
    if config["inputs"]["query_rds"] is not None and config["settings"]["split"]:
            return config["outputs"]["output_dir"] + "split/split_object_{pool}.RDS"
    elif config["inputs"]["query_rds"] is not None:
            return config["inputs"]["query_rds"]
    else:
            return config["outputs"]["output_dir"] + "make_seurat/query_{pool}.RDS"


rule map_azimuth:
    input:
        file = get_map_input,
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["azimuth_reference"] if config["settings"]["run_azimuth"] else []
    output:
        umap = report(config["outputs"]["output_dir"] + "map/azimuth_{pool}_ref_umap.png", category="Azimuth", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        pca = report(config["outputs"]["output_dir"] + ("map/azimuth_{pool}_ref_spca.png" if config["settings"]["cite_seq"] else "map/azimuth_{pool}_ref_pca.png"), category="Azimuth", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        metadata = config["outputs"]["output_dir"] + "map/azimuth_{pool}.metadata.tsv.gz"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["map"]["azimuth_memory"] * config["map"]["azimuth_threads"] - config["settings_extra"]["memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["azimuth_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["azimuth_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["map"]["azimuth_time"]]
    threads: config["map"]["azimuth_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/" + ("map_azimuth_cite_seq.R" if config["settings"]["cite_seq"] else "map_azimuth.R"),
        batch = "" if config["settings"]["split"] else "--batch Pool",
        refdata = config["settings"]["refdata"],
        plan = config["map_extra"]["azimuth_plan"],
        palette = "--palette '{}'".format(config["settings"]["palette"]) if config["settings"]["palette"] is not None else "",
        out = "azimuth_{pool}",
        path = config["outputs"]["output_dir"] + "map/"
    log: config["outputs"]["output_dir"] + "log/map_azimuth.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {input.file} \
            {params.batch} \
            --reference {input.reference} \
            --refdata '{params.refdata}' \
            --plan {params.plan} \
            --workers {threads} \
            --mem {resources.mem} \
            {params.palette} \
            --out {params.out} \
            --path {params.path}
        """


rule map_hierscpred:
    input:
        file = get_map_input,
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["hierscpred_reference"] if config["settings"]["run_hierscpred"] else []
    output:
        metadata = config["outputs"]["output_dir"] + "map/hier_scpred_{pool}.metadata.tsv.gz"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["map"]["hierscpred_memory"] * config["map"]["hierscpred_threads"] - config["settings_extra"]["memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["hierscpred_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["hierscpred_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["map"]["hierscpred_time"]]
    threads: config["map"]["hierscpred_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/map_hierscpred.R",
        batch = "" if config["settings"]["split"] else "--batch Pool",
        thr = config["map_extra"]["hierscpred_thr"],
        iter = config["map_extra"]["hierscpred_iter"],
        plan = config["map_extra"]["hierscpred_plan"],
        out = "hier_scpred_{pool}",
        path = config["outputs"]["output_dir"] + "map/"
    log: config["outputs"]["output_dir"] + "log/map_hierscpred.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {input.file} \
            {params.batch} \
            --reference {input.reference} \
            --thr {params.thr} \
            --iter {params.iter} \
            --plan {params.plan} \
            --workers {threads} \
            --mem {resources.mem} \
            --out {params.out} \
            --path {params.path}
        """


rule reduce:
    input:
        files = lambda wildcards: expand(config["outputs"]["output_dir"] + "map/{method}_{pool}.metadata.tsv.gz", method=wildcards.method, pool=POOLS),
        poolsheet = config["inputs"]["poolsheet_path"]
    output:
        file = config["outputs"]["output_dir"] + "map/{method}.metadata.tsv.gz"
    wildcard_constraints:
        method="[a-zA-Z]+"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["reduce"]["reduce_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["reduce"]["reduce_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["reduce"]["reduce_time"]]
    threads: config["reduce"]["reduce_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/reduce.R",
        indir = config["outputs"]["output_dir"] + "map/",
        path = config["outputs"]["output_dir"] + "map/"
    log: config["outputs"]["output_dir"] + "log/reduce.{method}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --indir {params.indir} \
            --method {wildcards.method} \
            --path {params.path}
        """


rule visualise_azimuth:
    input:
        metadata = config["outputs"]["output_dir"] + "map/azimuth.metadata.tsv.gz" if len(POOLS) > 1 else config["outputs"]["output_dir"] + "map/azimuth_{pool}.metadata.tsv.gz".format(pool=POOLS[0])
    output:
        proportions = report(expand(config["outputs"]["output_dir"] + "visualise/azimuth_predicted.{ct_col}_proportion.png", ct_col=[refdata.split("=")[0] for refdata in config["settings"]["refdata"].split(";")]), category="Azimuth combined", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        scores = report(expand(config["outputs"]["output_dir"] + "visualise/azimuth_predicted.{ct_col}_score.png", ct_col=[refdata.split("=")[0] for refdata in config["settings"]["refdata"].split(";")]), category="Azimuth combined", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        done = config["outputs"]["output_dir"] + "visualise/visualise.azimuth.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["visualise"]["visualise_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["visualise"]["visualise_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["visualise"]["visualise_time"]]
    threads: config["visualise"]["visualise_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/visualise.R",
        columns = config["settings"]["refdata"],
        palette = "--palette '{}'".format(config["settings"]["palette"]) if config["settings"]["palette"] is not None else "",
        out = "azimuth",
        path = config["outputs"]["output_dir"] + "visualise/"
    log: config["outputs"]["output_dir"] + "log/visualise.azimuth.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --metadata {input.metadata} \
            --columns '{params.columns}' \
            {params.palette} \
            --out {params.out} \
            --path {params.path}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """


rule visualise_hierscpred:
    input:
        metadata = config["outputs"]["output_dir"] + "map/hier_scpred.metadata.tsv.gz" if len(POOLS) > 1 else config["outputs"]["output_dir"] + "map/hier_scpred_{pool}.metadata.tsv.gz".format(pool=POOLS[0])
    output:
        proportions = report(config["outputs"]["output_dir"] + "visualise/hier_scpred_scpred_prediction_proportion.png", category="Hierarchical scPred combined", caption=config["inputs"]["repo_dir"] + "report_captions/hier_scpred.rst"),
        done = config["outputs"]["output_dir"] + "visualise/visualise.hier_scpred.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["visualise"]["visualise_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["visualise"]["visualise_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["visualise"]["visualise_time"]]
    threads: config["visualise"]["visualise_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/visualise.R",
        columns = "scpred_prediction",
        palette = "--palette '{}'".format(config["settings"]["palette"]) if config["settings"]["palette"] is not None else "",
        out = "hier_scpred",
        path = config["outputs"]["output_dir"] + "visualise/"
    log: config["outputs"]["output_dir"] + "log/visualise.hier_scpred.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --metadata {input.metadata} \
            --columns '{params.columns}' \
            {params.palette} \
            --out {params.out} \
            --path {params.path}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """


rule compare:
    input:
        metadata1 = config["outputs"]["output_dir"] + "map/azimuth.metadata.tsv.gz" if len(POOLS) > 1 else config["outputs"]["output_dir"] + "map/azimuth_{pool}.metadata.tsv.gz".format(pool=POOLS[0]),
        metadata2 = config["outputs"]["output_dir"] + "map/hier_scpred.metadata.tsv.gz" if len(POOLS) > 1 else config["outputs"]["output_dir"] + "map/hier_scpred_{pool}.metadata.tsv.gz".format(pool=POOLS[0])
    output:
        counts_heatmap = report(config["outputs"]["output_dir"] + "compare/comp_heatmap_counts.pdf", category="Compare", caption=config["inputs"]["repo_dir"] + "report_captions/compare.rst"),
        prop_heatmap = report(config["outputs"]["output_dir"] + "compare/comp_heatmap_prop.pdf", category="Compare", caption=config["inputs"]["repo_dir"] + "report_captions/compare.rst"),
        contingency_table = config["outputs"]["output_dir"] + "compare/comp_contingency_table.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["compare"]["compare_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["compare"]["compare_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["compare"]["compare_time"]]
    threads: config["compare"]["compare_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/compare.R",
        xaxis = config["compare_extra"]["xaxis"],
        yaxis = config["compare_extra"]["yaxis"],
        sort = config["compare_extra"]["sort"],
        out = "comp",
        path = config["outputs"]["output_dir"] + "compare/"
    log: config["outputs"]["output_dir"] + "log/compare.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --metadata1 {input.metadata1} \
            --metadata2 {input.metadata2} \
            --xaxis {params.xaxis} \
            --yaxis {params.yaxis} \
            --sort {params.sort} \
            --out {params.out} \
            --path {params.path}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """