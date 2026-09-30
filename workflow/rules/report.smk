rule report:
    input: rules.coloc.output
    output:
        f"{RUN_DIR}/report/report.html",
        f"{RUN_DIR}/report/report.md",
        f"{RUN_DIR}/report/summary.tsv",
    shell: "python -m gwas2mechanism.workflow_stage report --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
