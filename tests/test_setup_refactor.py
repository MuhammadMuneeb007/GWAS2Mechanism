from __future__ import annotations

import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import ClassVar

import pytest

from gwas2mechanism.config import load_config
from gwas2mechanism.reference import manager as reference_module
from gwas2mechanism.reference.manager import ReferencePanel
from gwas2mechanism.resource_status import inspect_resources
from gwas2mechanism.setup_resources import setup_gtex
from gwas2mechanism.slurm_setup import generate_slurm_jobs
from gwas2mechanism.utils import download as download_module
from gwas2mechanism.utils.download import Downloader, DownloadOptions, DownloadResult
from gwas2mechanism.utils.proc import ToolError
from gwas2mechanism.utils.progress import is_batch, progress
from gwas2mechanism.utils.tools import ToolResolver


def _config(tmp_path: Path):
    return load_config(overrides={"resources_dir": str(tmp_path / "resources"), "genome": {"chromosomes": ["1", "2"]}})


def _dataset(panel: ReferencePanel, population: str | None, chrom: str) -> None:
    directory = panel.all_dir if population is None else panel.pop_dir(population)
    prefix = panel.all_prefix(chrom) if population is None else panel.prefix(population, chrom)
    parquet = panel.all_pvar_parquet(chrom) if population is None else panel.pvar_parquet(population, chrom)
    directory.mkdir(parents=True, exist_ok=True)
    outputs = [Path(f"{prefix}.{suffix}") for suffix in ("pgen", "pvar", "psam")] + [parquet]
    for path in outputs:
        path.write_bytes(b"x")
    marker = panel.chromosome_marker(directory, chrom)
    marker.write_text(json.dumps({"status": "COMPLETE", "n_samples": 1, "n_variants": 1, "output_sizes": {path.name: 1 for path in outputs}}), encoding="utf-8")


def test_http_downloader_retries() -> None:
    options = DownloadOptions(retry_attempts=3, retry_initial_seconds=0, retry_max_seconds=0)
    downloader = Downloader(options=options)
    calls = 0

    def flaky(_client, _url, part, _expected):
        nonlocal calls
        calls += 1
        if calls < 3:
            raise OSError("temporary")
        part.write_bytes(b"ok")

    downloader._httpx = flaky  # type: ignore[method-assign]
    downloader._httpx_with_retry(object(), "https://example.invalid/a", Path("unused.part"), 2)  # type: ignore[arg-type]
    Path("unused.part").unlink()
    assert calls == 3


def test_aria2_failure_falls_back_to_http(tmp_path: Path, monkeypatch) -> None:
    downloader = Downloader(use_aria2=["aria2c"], options=DownloadOptions(retry_attempts=1))
    monkeypatch.setattr(downloader, "_aria2_with_retry", lambda *_: (_ for _ in ()).throw(ToolError("failed")))
    monkeypatch.setattr(downloader, "_httpx_with_retry", lambda _client, _url, part, _size: part.write_bytes(b"payload"))
    result = downloader.fetch("https://example.invalid/file", tmp_path / "file", expected_size=7)
    assert result.status == "DOWNLOADED"
    assert result.path.read_bytes() == b"payload"


def test_http_resume_initializes_from_part(tmp_path: Path, monkeypatch) -> None:
    part = tmp_path / "file.part"
    part.write_bytes(b"abc")
    seen: dict[str, object] = {}

    class Response:
        status_code = 206
        headers: ClassVar[dict[str, str]] = {"content-length": "3"}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def raise_for_status(self):
            return None

        def iter_bytes(self, _chunk):
            yield b"def"

    class Client:
        def stream(self, _method, _url, headers):
            seen["range"] = headers["Range"]
            return Response()

    class Bar:
        n = 3

        def update(self, value):
            seen["updated"] = value

        def close(self):
            pass

    monkeypatch.setattr(download_module, "progress", lambda **kwargs: seen.update(kwargs) or Bar())
    Downloader()._httpx(Client(), "https://example.invalid/file", part, 6)  # type: ignore[arg-type]
    assert part.read_bytes() == b"abcdef"
    assert seen["initial"] == 3
    assert seen["range"] == "bytes=3-"


def test_verification_stamp_reuse(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.write_bytes(b"same")
    destination = tmp_path / "destination"
    first = Downloader().fetch(str(source), destination)
    first_mtime = destination.stat().st_mtime_ns
    second = Downloader().fetch(str(source), destination)
    assert first.status == "COPIED"
    assert second.status == "EXISTS"
    assert destination.stat().st_mtime_ns == first_mtime


def test_gtex_downloads_are_parallel(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    monkeypatch.setattr("gwas2mechanism.setup_resources.discover_release", lambda *_: "v11")
    objects = [{"name": f"v11/tissue{i}_eqtl_susie.parquet", "size": 1, "md5Hash": None} for i in range(4)]
    monkeypatch.setattr("gwas2mechanism.setup_resources._gtex_objects", lambda *_: objects)
    active = 0
    maximum = 0
    lock = threading.Lock()

    def fake_fetch(self, url, path, **_kwargs):
        nonlocal active, maximum
        with lock:
            active += 1
            maximum = max(maximum, active)
        time.sleep(0.03)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"x")
        with lock:
            active -= 1
        return DownloadResult(url, path, 1, None, "DOWNLOADED")

    monkeypatch.setattr(Downloader, "fetch", fake_fetch)
    result = setup_gtex(cfg)
    assert result["gtex_objects"] == 4
    assert maximum > 1


def test_resource_status_is_read_only(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    monkeypatch.setattr("httpx.get", lambda *_a, **_k: pytest.fail("status must not use network"))
    report = inspect_resources(cfg)
    assert report["missing"] > 0
    assert report["resource_cache"] == str(cfg.resources_path.resolve())


def test_chromosome_completion_requires_every_output(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    panel = ReferencePanel(cfg.reference, cfg.resources_path, ToolResolver(cfg.tools))
    _dataset(panel, None, "1")
    assert panel.all_chromosome_complete("1")
    panel.all_pvar_parquet("1").unlink()
    assert not panel.all_chromosome_complete("1")


def test_population_sample_lists_generated_once(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    panel = ReferencePanel(cfg.reference, cfg.resources_path, ToolResolver(cfg.tools))
    panel.raw.mkdir(parents=True)
    (panel.raw / reference_module.PANEL_FILE).write_text("sample\tpop\tsuper_pop\tgender\nA\tCEU\tEUR\t1\nB\tYRI\tAFR\t2\n", encoding="utf-8")
    paths = panel.write_sample_lists(["EUR", "AFR"])
    before = paths["EUR"].stat().st_mtime_ns
    assert paths["EUR"].read_text(encoding="utf-8") == "#IID\nA\n"
    panel.write_sample_lists(["EUR"])
    assert paths["EUR"].stat().st_mtime_ns == before


def test_slurm_scripts_and_population_scope(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("GWAS2M_PROJECT_ROOT", str(tmp_path))
    cfg = _config(tmp_path)
    plan = generate_slurm_jobs(cfg, reference=True, gtex=False, vep=False, splice=False, populations=["EUR"], partition="ascher")
    scripts = plan["scripts"]
    assert Path(scripts["reference_download"]).exists()
    population = Path(scripts["reference_population"]).read_text(encoding="utf-8")
    assert "POPS=(EUR)" in population
    assert "AFR" not in population
    assert "#SBATCH --partition=ascher" in population


def test_slurm_submit_uses_parsable_dependencies(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("GWAS2M_PROJECT_ROOT", str(tmp_path))
    plan = generate_slurm_jobs(_config(tmp_path), reference=True, gtex=True, vep=True, splice=True, populations=["EUR"])
    submit = Path(plan["submit_script"]).read_text(encoding="utf-8")
    assert "sbatch --parsable" in submit
    assert "--dependency=afterok:$RAW_JOB" in submit
    assert "GTEX_JOB=$(submit --dependency=afterok:$ENV_JOB" in submit
    assert "$ENV_JOB:$GENOME_JOB" in submit


def test_completed_conversion_is_skipped(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    panel = ReferencePanel(cfg.reference, cfg.resources_path, ToolResolver(cfg.tools))
    _dataset(panel, None, "1")
    monkeypatch.setattr(reference_module, "run", lambda *_a, **_k: pytest.fail("completed chromosome reran"))
    panel.convert_all_chromosome("1", 1)


def test_partial_conversion_reruns_and_completes(tmp_path: Path, monkeypatch) -> None:
    cfg = _config(tmp_path)
    cfg.tools["plink2"] = ["plink2"]
    panel = ReferencePanel(cfg.reference, cfg.resources_path, ToolResolver(cfg.tools))
    panel.raw.mkdir(parents=True)
    for name in (panel.vcf_name("1"), f"{panel.vcf_name('1')}.tbi"):
        path = panel.raw / name
        path.write_bytes(b"raw")
        path.with_name(path.name + ".verified.json").write_text("{}", encoding="utf-8")
    calls = 0

    def fake_run(command, **_kwargs):
        nonlocal calls
        calls += 1
        prefix = Path(command[command.index("--out") + 1])
        Path(f"{prefix}.pgen").write_bytes(b"pgen")
        Path(f"{prefix}.pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n1\t1\t1:1:A:G\tA\tG\n", encoding="utf-8")
        Path(f"{prefix}.psam").write_text("#IID\nA\n", encoding="utf-8")

    monkeypatch.setattr(reference_module, "run", fake_run)
    panel.convert_all_chromosome("1", 1)
    assert calls == 1
    assert panel.all_chromosome_complete("1")


def test_two_workers_do_not_corrupt_same_download(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.write_bytes(b"safe" * 1000)
    destination = tmp_path / "destination"
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: Downloader().fetch(str(source), destination), range(2)))
    assert destination.read_bytes() == source.read_bytes()
    assert sorted(result.status for result in results) == ["COPIED", "EXISTS"]


def test_no_pooled_cross_ancestry_panel(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    panel = ReferencePanel(cfg.reference, cfg.resources_path, ToolResolver(cfg.tools))
    with pytest.raises(ValueError, match="pooled"):
        panel.pop_dir("ALL")
    with pytest.raises(ValueError, match="pooled"):
        panel.subset_population_chromosome("MULTI", "1", 1)


def test_tqdm_batch_and_redirected_modes(monkeypatch) -> None:
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    assert is_batch()
    bar = progress(total=1, desc="test", disable=True)
    bar.update(1)
    bar.close()
