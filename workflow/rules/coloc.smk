rule coloc:
    input: rules.qtl.output, rules.splice.output
    output: f"{RUN_DIR}/.stages/coloc.done"
    shell: "python -m gwas2mechanism.workflow_stage coloc --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
