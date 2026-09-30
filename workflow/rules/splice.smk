rule splice:
    input: rules.annotate.output
    output: f"{RUN_DIR}/.stages/splice.done"
    shell: "python -m gwas2mechanism.workflow_stage splice --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
