rule discover:
    output: f"{RUN_DIR}/.stages/discover.done"
    shell: "python -m gwas2mechanism.workflow_stage discover --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
