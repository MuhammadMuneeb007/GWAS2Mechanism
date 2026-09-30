rule harmonize:
    input: rules.download.output
    output: f"{RUN_DIR}/.stages/harmonize.done"
    shell: "python -m gwas2mechanism.workflow_stage harmonize --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
