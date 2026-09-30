rule download:
    input: rules.discover.output
    output: f"{RUN_DIR}/.stages/download.done"
    shell: "python -m gwas2mechanism.workflow_stage download --run-dir {RUN_DIR:q} --phenotype {PHENOTYPE:q}"
