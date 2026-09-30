rule qtl:
    input: rules.annotate.output
    output: f"{RUN_DIR}/.stages/qtl.done"
    shell: "python -m gwas2mechanism.workflow_stage qtl --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
