"""
This rule runs the notebook
"""
rule notebook_graphs:
    input:
        count = "aggregated_salmon_counts/NumReads.sf.tsv"
    output:
        "logs/notebook/processed_notebook.ipynb"
    message:
        "Testing notebooks"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(
                attempt * len(sample_id_list) * 250, 10240
            )
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 60)
        )
    conda:
        "../envs/py37.yaml"
    log:
        notebook = "logs/notebook/processed_notebook.ipynb"
    notebook:
        "../notebook/Explore_results.ipynb"


"""
This rule converts notebook into a stad-alone code-free html file.
"""
rule notebook_html:
    input:
        "logs/notebook/processed_notebook.ipynb"
    output:
        report(
            "notebook/notebook.html",
            caption="../report/notebook.rst",
            category="Notebook"
        )
    message:
        "Converting notebook as HTML"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(
                attempt * len(sample_id_list) * 250, 10240
            )
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 60)
        )
    conda:
        "../envs/py37.yaml"
    params:
    log:
        "logs/notebook/html_conversion.log"
    shell:
        "jupyter nbconvert "
        "--to html "
        "--TemplateExporter.exclude_input=True "
        "--output notebook.html "
        "--output-dir notebook/ "
        "{input} "
        "> {log} 2>&1"
