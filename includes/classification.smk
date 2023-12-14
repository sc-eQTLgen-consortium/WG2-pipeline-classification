#!/usr/bin/env python


rule split:
    input:
        file = config["inputs"]["query_rds"]
    output:
        files = temp(expand(config["outputs"]["output_dir"] + "split/split_object_{batch}.RDS", batch=BATCHES))
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["split"]["split_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["split"]["split_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["split"]["split_time"]]
    threads: config["split"]["split_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/split.R",
        batch = BATCH,
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


rule map_azimuth:
    input:
        file = config["outputs"]["output_dir"] + "split/split_object_{batch}.RDS" if len(BATCHES) > 1 else config["inputs"]["query_rds"],
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["azimuth"]
    output:
        umap = report(config["outputs"]["output_dir"] + "map/azimuth_{batch}_ref_umap.png", category="Azimuth", subcategory="{batch}", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        pca = report(config["outputs"]["output_dir"] + ("map/azimuth_{batch}_ref_spca.png" if config["settings"]["citeseq"] else "map/azimuth_{batch}_ref_pca.png"), category="Azimuth", subcategory="{batch}", caption=config["inputs"]["repo_dir"] + "report_captions/azimuth.rst"),
        file = temp(config["outputs"]["output_dir"] + "map/azimuth_{batch}.RDS") if len(BATCHES) > 1 else config["outputs"]["output_dir"] + "map/azimuth_{batch}.RDS"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["map"]["azimuth_memory"] * config["map"]["azimuth_threads"] - config["settings_extra"]["memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["azimuth_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["map"]["azimuth_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["map"]["azimuth_time"]]
    threads: config["map"]["azimuth_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/" + ("map_azimuth_citeseq.R" if config["settings"]["citeseq"] else "map_azimuth.R"),
        refdata = config["settings"]["refdata"],
        plan = config["map_extra"]["azimuth_plan"],
        out = "azimuth_{batch}",
        path = config["outputs"]["output_dir"] + "map/"
    log: config["outputs"]["output_dir"] + "log/map_azimuth.{batch}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {input.file} \
            --reference {input.reference} \
            --refdata {params.refdata} \
            --plan {params.plan} \
            --workers {threads} \
            --mem {resources.mem} \
            --out {params.out} \
            --path {params.path}
        """


rule map_hierscpred:
    input:
        file = config["outputs"]["output_dir"] + "split/split_object_{batch}.RDS" if not config["settings"]["azimuth"] else config["outputs"]["output_dir"] + "map/azimuth_{batch}.RDS",
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["hierscpred"]
    output:
        file = temp(config["outputs"]["output_dir"] + "map/hier_scpred_{batch}.RDS") if len(BATCHES) > 1 else config["outputs"]["output_dir"] + "map/hier_scpred_{batch}.RDS"
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
        thr = config["map_extra"]["hierscpred_thr"],
        iter = config["map_extra"]["hierscpred_iter"],
        plan = config["map_extra"]["hierscpred_plan"],
        out = "hier_scpred_{batch}",
        path = config["outputs"]["output_dir"] + "map/"
    log: config["outputs"]["output_dir"] + "log/map_hierscpred.{batch}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {input.file} \
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
        files = expand(config["outputs"]["output_dir"] + "map/hier_scpred_{batch}.RDS", batch=BATCHES) if config["settings"]["hierscpred"] else expand(config["outputs"]["output_dir"] + "map/azimuth_{batch}.RDS", batch=BATCHES)
    output:
        file = config["outputs"]["output_dir"] + "reduce/reduced_data.RDS"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["reduce"]["reduce_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["reduce"]["reduce_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["reduce"]["reduce_time"]]
    threads: config["reduce"]["reduce_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/reduce.R",
        file = config["outputs"]["output_dir"] + "map/",
        out = "reduced_data",
        path = config["outputs"]["output_dir"] + "reduce/"
    log: config["outputs"]["output_dir"] + "log/reduce.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --file {params.file} \
            --out {params.out} \
            --path {params.path}
        """


rule compare:
    input:
        file = config["outputs"]["output_dir"] + "reduce/reduced_data.RDS" if len(BATCHES) > 1 else config["outputs"]["output_dir"] + "map/hier_scpred_{batch}.RDS".format(batch=BATCHES[0])
    output:
        counts_heatmap = report(config["outputs"]["output_dir"] + "compare/comp_heatmap_counts.pdf", category="Compare", caption=config["inputs"]["repo_dir"] + "report_captions/compare.rst"),
        prop_heatmap = report(config["outputs"]["output_dir"] + "compare/comp_heatmap_prop.pdf", category="Compare", caption=config["inputs"]["repo_dir"] + "report_captions/compare.rst"),
        contingency_table = config["outputs"]["output_dir"] + "compare/comp_contingency_table.tsv"
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
            --file {input.file} \
            --xaxis {params.xaxis} \
            --yaxis {params.yaxis} \
            --sort {params.sort} \
            --out {params.out} \
            --path {params.path}
        """