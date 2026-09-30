# Legacy migration map

The files beside `GWAS2Mechanism/` are read-only reference implementations. They were not copied or modified. This map records every legacy script and every declared function/class, the useful scientific behavior retained, and the phenotype-specific behavior removed.

## Step01–Step02: discovery and summary statistics

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `Step01_Find_CAD_Trait.py` (top-level Catalog query) | `discovery/phenotype.py`, `discovery/gwas_catalog.py`, `discovery/service.py` | Replaced free-text-only searching with ontology/synonym-aware resolution, match reasons, canonical provenance, and backend fallback. |
| `Step02_EUR_CAD_GWAS.py::normalise_text`, `is_exact_cad` | `discovery/phenotype.py` | Generalized to normalized exact mapped-trait, reported-trait, synonym, and controlled fuzzy matching. |
| `estimate_sample_size` | `discovery/ancestry.py::parse_sample_size`, `apply_sample_size` | Retains parsed values and sources; never invents counts. |
| `build_gwas_ftp_directory`, `has_harmonised_sumstats` | `discovery/base.py`, `discovery/gwas_catalog.py` | Generalized GCST directory enumeration, bounded header probing, and raw/harmonised classification. |
| `download_sumstats` | `sumstats/download.py` | Resumable verified download with GWASLab preference and robust HTTP fallback. |
| remaining Step02 top-level ranking/selection | `discovery/ranking.py`, `discovery/service.py` | Replaced first-hit/EUR-only selection with an exported, component-wise ranking and strongest usable study per actual population. |

## Step03.0: 1000 Genomes download

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `safe_print`, `banner`, `human_size` | Rich CLI and standard logging | Presentation-only helpers replaced. |
| `chromosome_filename`, `build_download_list` | `reference/manager.py::vcf_name`, `download` | Release/template driven; all autosomes configurable. |
| `create_session`, `get_remote_size`, `download_file` | `utils/download.py::Downloader` | Retains retries, progress/resume intent, adds `.part` atomic rename, checksum/size stamps, aria2c acceleration, and cache reuse. |
| `main` | `setup_resources.py`, `gwas2m setup --reference` | Downloads once and records a resource manifest. |

## Step03.1: population reference preparation

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `run_command`, `check_program`, `human_size` | `utils/proc.py`, `utils/tools.py`, CLI logging | Centralized safe subprocess execution and tool resolution. |
| `raw_vcf_path`, `check_raw_reference`, `find_panel` | `reference/manager.py` | Versioned cache paths and completeness markers. |
| `create_european_sample_list` | `ReferencePanel.write_sample_list` | Generalized to EUR/AFR/EAS/SAS/AMR using 1000 Genomes `super_pop`; no pooled panel. |
| `subset_european_chromosome`, `convert_to_plink` | `ReferencePanel.prepare` | Generalized population-by-population PLINK2 preparation with normalized IDs. |
| `get_psam_sample_count`, `get_pvar_variant_count` | `ReferencePanel.sample_size`, cached PVAR Parquet | Validation and counts retained without repeated text scans. |
| `process_chromosome`, `main` | `ReferencePanel.download/prepare`, Snakemake `references.smk` | Restartable per population/chromosome. |

## Step03.2: clumping and loci

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `check_program`, `run_command`, `reference_prefix` | shared tool/process utilities and `ReferencePanel.prefix` | Centralized and population-aware. |
| `check_reference` | `ReferencePanel.is_ready` | Verifies the requested population rather than assuming EUR. |
| `clean_rsid`, `load_gwas` | `sumstats/harmonize.py`, `utils/variants.py` | Canonical schema, lazy Parquet scans, normalized IDs. |
| `match_chromosome` | `finemap/region.py::match_to_reference` | Retains unordered allele-pair matching and effect orientation; vectorized in Polars. |
| `create_clump_input`, `run_clumping`, `read_clumps` | `clumping/loci.py::run_plink_clump` and workflow parser | PLINK2 is invoked separately against each population's panel. |
| `clean_sp2_id`, `extract_leads`, `build_clump_members` | `clumping/loci.py::greedy_ld_clump` (fixtures) and PLINK output normalization | Member/lead provenance retained. |
| `construct_loci`, `annotate_member_table` | `clumping/loci.py::construct_loci`, `union_loci` | Configurable windows and union across real population results; never pooled-LD clumping. |
| `main` | `workflow/rules/clump.smk` | Parallel/restartable by population/chromosome. |

## Step04: fine-mapping

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `run_command`, `check_program`, `reference_prefix`, `check_inputs` | shared utilities, reference manager, workflow validation | Generalized and centrally validated. |
| `clean_rsid`, `allele_key` | `utils/variants.py` | One normalized allele-key implementation. |
| `check_susie` | `utils/tools.py`, `gwas2m doctor` | Reports installed R/susieR availability. |
| `load_loci`, `regional_gwas_path`, `extract_all_regional_gwas` | `finemap/region.py::extract_region` | Lazy Parquet region extraction; no fixed locus count. |
| `create_regional_reference_pvar`, `read_pvar` | `ReferencePanel.region_variants/cache_pvar` | Cached per-population Parquet variant indices. |
| `match_gwas_to_reference` | `finemap/region.py::match_to_reference` | Preserves allele/frequency/effect orientation logic. |
| `calculate_ld` | `finemap/ld.py::LDProvider.compute`, `stabilize` | Signed per-population LD, size verification, mmap binary reads, diagnostics and ridge stabilization. |
| `write_susie_r_script`, `run_susie` | `finemap/susie.py::run_susier` | Production `susieR::susie_rss`; explicit sample size and diagnostics. |
| `process_locus`, `combine_results` | `finemap/susie.py::FineMapResult`, workflow `finemap.smk` | Standard Parquet outputs for variants, credible sets, summaries, diagnostics, and failures. |
| `main` | `workflow/rules/finemap.smk` | Population/locus jobs; no accession, sample count, phenotype, or locus-number constants. |

Cross-ancestry capability did not exist in Step04. `finemap/multisusie.py` adds separate `z`, `n`, and LD inputs per population; `finemap/susiex.py` is an optional fallback. `meta/fixed_effect.py` adds allele-aligned fixed-effect meta-analysis with Q, Q p-value, I², direction, and contributing populations.

## Step05: VEP annotation

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `run_command`, `check_program`, `check_step04`, `check_vep_cache` | shared utilities, `doctor`, workflow prerequisites, `setup_resources.py` | Centralized checks and automated cache preparation. |
| `load_finemapped_variants`, `select_variants` | `annotation/vep.py::select_prioritised_variants` | Retains credible sets, PIP threshold, and per-locus top variant; unions/deduplicates across populations and cross-ancestry results. |
| `create_vep_vcf`, `run_vep`, `read_vep_output` | `write_vep_vcf`, `run_vep`, `parse_vep_tsv` | GRCh38 local-cache VEP; one run for the union. |
| `functional_category` | `annotation/vep.py::functional_category` | Retains SPLICING/CODING/UTR/REGULATORY/NONCODING_RNA/INTRONIC/UPSTREAM/DOWNSTREAM/INTERGENIC/OTHER. |
| `merge_annotations`, `select_best_consequence`, `collapse_values`, `make_variant_summary`, `create_summary_tables` | VEP parser, evidence integration, report | Consequence rows and one best consequence retained in separate Parquet outputs. |
| `main` | `workflow/rules/annotate.smk` | Restartable prioritized-only annotation. |

## Step06: GTEx integration

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `chunks`, `create_session`, `api_get`, `api_get_paginated` | `qtl/gtex.py`, `utils/download.py` | Compact fine-mapping archives and columnar scans replace high-volume row/API loops in fast mode. |
| `check_input`, `load_variants` | prioritized variant union from annotation/workflow | One deduplicated cross-population input. |
| `get_tissue_metadata`, `select_tissues` | `qtl/gtex.py::requested_tissues`, `discover_release` | Default is every tissue found in the stable release; no phenotype/heart prefilter. |
| `map_variants_to_gtex` | `parse_gtex_variant`, `_normalise_qtl`, `integrate_local_archive` | Exact normalized `CHR:POS:REF:ALT` matching. |
| `query_qtl`, `add_variant_context`, `collapse`, `build_variant_gene_evidence`, `top_association`, `build_variant_summary`, `build_locus_summary` | QTL integration outputs, evidence layer, report | Gene/tissue/QTL type/PIP/context retained in long form and summaries. |
| `main` | `workflow/rules/qtl.smk` | All tissues, eQTL+sQTL; optional apaQTL. |

## Step07: SpliceAI

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `check_program`, `run_command`, `check_software`, `check_input` | shared tools/process utilities and workflow prerequisites | Isolated pinned environment support. |
| `load_variants` | `annotation.vep.select_prioritised_variants` | Union/deduplication across populations before prediction. |
| `find_vep_fasta`, `prepare_reference_fasta` | `setup_resources.py` | Shared cached GRCh38 FASTA plus index. |
| `create_raw_vcf`, `normalize_vcf` | `annotation.vep.write_vep_vcf`; bcftools workflow | Canonical variant VCF generation; production workflow normalizes before scoring. |
| `run_spliceai`, `parse_info_field`, `safe_float`, `safe_int`, `parse_spliceai` | `splicing/spliceai.py::run_spliceai`, `parse_spliceai_info`, `parse_spliceai_vcf` | Independent execution and gene-level/all-prediction parsing. |
| `spliceai_class`, `boolean_series` | `splicing/spliceai.py::classify_score`, native Polars booleans | Threshold meanings retained. |
| `create_variant_summary`, `create_prioritized_table`, `create_locus_summary`, `create_category_summary` | splice summary, evidence and report layers | Scientific context is merged after prediction. |
| `main` | `workflow/rules/splice.smk` | One execution per unique variant union. |

## Step07B: Pangolin

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `require`, `open_text`, `count_vcf`, `get_info`, `to_float`, `parse_pos_score`, `as_bool`, `run_and_tee`, `check_inputs` | shared tools/process helpers and `splicing/pangolin.py` parser | Common infrastructure replaces script-local utilities. |
| `prepare_pangolin_csv`, `run_pangolin`, `read_input_variants`, `parse_pangolin` | `splicing/pangolin.py::run_pangolin`, `parse_pangolin_table` | Uses the same deduplicated union but remains independent of SpliceAI. |
| `make_variant_summary`, `integrate`, `write_outputs`, `final_summary` | `integrate_splicing`, evidence and report layers | Predictor-specific values remain separate columns. |
| `main` | `workflow/rules/splice.smk` | Pinned isolated environment. |

## Step08: fine-mapping-based colocalisation screen

| Legacy symbol | New location | Migration decision |
|---|---|---|
| `banner`, `ensure_dirs`, `require`, `check_pyarrow`, `download_with_progress`, `extract_tar` | shared I/O/download/setup infrastructure | Atomic/cached resource handling. Archive extraction is controlled by setup. |
| `first_existing`, `norm_chr`, `norm_allele`, `key_from_parts`, `parse_gtex_variant` | canonical schemas, `utils/variants.py`, `qtl/gtex.py` | Exact normalized variant identity. |
| `load_cad`, `load_spliceai` | phenotype-independent evidence inputs | No fixed phenotype paths. |
| `find_parquets`, `tissue_from_filename`, `should_process_tissue`, `normalize_qtl_parquet` | `qtl/gtex.py::integrate_local_archive`, `_normalise_qtl` | Every discovered tissue is processed by default. |
| `multi_causal_lclpp` | `coloc/screen.py::locus_clpp` | Vectorized and numerically stable `1 - product(1-vCLPP)`. |
| `summarize_pairs` | `coloc/screen.py::screen` | Produces variant and gene/tissue summaries per GWAS population and cross-ancestry result. |
| `add_spliceai`, `gene_summary` | evidence integration and report | Prediction context added after the coloc screen. |
| `main` | `workflow/rules/coloc.smk` | Explicitly labeled a fine-mapping-based screen, not formal coloc. |

## Optional Matt/lncRNA experiments and reporting

These are not mandatory core stages. The future optional `modules/lncrna/` route is disabled by default and must preserve the difference between “outside protein-coding CDS” and “inside an lncRNA gene body.”

| Legacy symbols | New destination/status |
|---|---|
| `Matt_Project_Map_CAD_lncRNA.py`: `banner`, `first_existing`, `norm_chr`, `chrom_sort_key`, `parse_attributes`, `broad_gene_class`, `safe_bool_series`, `variant_key`, `parse_variant_key`, `download_gtf_if_needed`, `load_cad_variants`, `merge_intervals`, `parse_gtf`, `ActiveIntervals`, `names_string`, `ids_string`, `map_variants`, `summarize_primary_classes`, `subset_masks`, `summarize_evidence_subsets`, `rank_lncrnas`, `summarize_loci`, `fisher_enrichment`, `yesno`, `build_conclusion`, `main` | Optional `modules/lncrna`; generic variant/GTF helpers reuse core normalization and cached GENCODE. Statistical summaries are deferred rather than silently included in core. |
| `Matt_Project_GENCODE_Validation.py`: `banner`, `first_existing`, `norm_chr`, `chrom_sort_key`, `parse_attributes`, `variant_key`, `parse_variant_key`, `safe_bool_series`, `download`, `merge_intervals`, `load_cad`, `prepare_gencode`, `parse_full_gencode`, `parse_lncrna_gencode`, `ActiveIntervals`, `join_names`, `join_ids`, `map_variants`, `subset_masks`, `subset_summary`, `class_summary`, `rank_lncrnas`, `locus_summary`, `enrichment`, `compare_ensembl`, `yn`, `make_conclusion`, `main` | Optional validation module/documentation; GENCODE download is centralized in setup. Ensembl-vs-GENCODE validation remains an opt-in experiment. |
| `SpliceHeartAI_Matt_Project_Findings.py`: `banner`, `first_existing`, `read_tsv`, `fmt_int`, `fmt_float`, `pct`, `add_result`, `summarize_step01`, `summarize_step02`, `summarize_step03`, `summarize_step04`, `summarize_step05`, `summarize_step06`, `summarize_step07`, `summarize_step08`, `summarize_matt_ensembl`, `summarize_matt_gencode`, `build_markdown`, `main` | Replaced by `reporting/report.py`, which reads run manifests and population-agnostic Parquet outputs. Optional lncRNA sections appear only when enabled. |

## Shell/job wrappers

| Legacy script | New location |
|---|---|
| `bioprs.sh` | `scripts/install.sh`, `scripts/bootstrap.sh`, `scripts/download_resources.sh`, `environment.yml` |
| `GPUJob.sh` | Snakemake execution plus `profiles/slurm/config.yaml`; no partition/account is hard-coded. |

## Deliberately removed assumptions

The new core contains no fixed phenotype, accession, ancestry, heart-tissue requirement, sample counts, or expected locus count. Published multi-ancestry results remain `COMBINED_PUBLISHED`; only genuinely stratified inputs can feed population-specific analyses or `COMBINED_META`. All expensive variant-level tools run once on a cross-population union.
