rule annotate:
    input: rules.multiancestry.output
    output: f"{RUN_DIR}/.stages/annotate.done"
    shell: "python -m gwas2mechanism.workflow_stage annotate --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
