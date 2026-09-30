rule finemap:
    input: rules.clump.output
    output: f"{RUN_DIR}/.stages/finemap.done"
    shell: "python -m gwas2mechanism.workflow_stage finemap --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
