rule references:
    input: rules.harmonize.output
    output: f"{RUN_DIR}/.stages/references.done"
    shell: "python -m gwas2mechanism.workflow_stage references --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
