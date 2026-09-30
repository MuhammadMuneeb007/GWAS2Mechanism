rule clump:
    input: rules.references.output
    output: f"{RUN_DIR}/.stages/clump.done"
    shell: "python -m gwas2mechanism.workflow_stage clump --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
