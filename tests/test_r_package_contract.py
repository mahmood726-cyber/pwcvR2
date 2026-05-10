import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_package_metadata_and_exports_are_present():
    description = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")
    namespace = (ROOT / "NAMESPACE").read_text(encoding="utf-8")
    source = (ROOT / "R" / "pwcv_r2.R").read_text(encoding="utf-8")

    assert "Package: pwcvR2" in description
    assert "Precision-Weighted Cross-Validation" in description
    assert "export(" in namespace
    assert "pwcv_r2 <- function" in source


def test_testthat_suite_is_wired():
    runner = (ROOT / "tests" / "testthat.R").read_text(encoding="utf-8")
    assert 'test_check("pwcvR2")' in runner
    assert (ROOT / "tests" / "testthat" / "test-pwcv_r2.R").exists()


def test_r_testthat_suite_passes_when_rscript_is_available():
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed in this environment")
    subprocess.run([rscript, "tests/testthat.R"], cwd=ROOT, check=True)
