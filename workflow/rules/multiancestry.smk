rule multiancestry:
    input: rules.finemap.output
    output: f"{RUN_DIR}/.stages/multiancestry.done"
    shell: "python -m gwas2mechanism.workflow_stage multiancestry --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
