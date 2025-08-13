import subprocess, sys, tempfile, csv, pathlib

def test_by_ihw_example_runs_and_writes_csv():
    tmp = tempfile.TemporaryDirectory()
    tmpdir = pathlib.Path(tmp.name)
    csv_path = tmpdir / "data.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pval", "cov"])
        for p, c in [(0.001, 0.9), (0.02, 0.8), (0.2, 0.5), (0.6, 0.2), (0.9, 0.1)]:
            w.writerow([p, c])
    out_csv = tmpdir / "ihw.csv"
    cmd = [
        sys.executable,
        "-m",
        "mc_cc0.cli",
        "ihw-example",
        "--csv",
        str(csv_path),
        "--p",
        "pval",
        "--covariate",
        "cov",
        "--alpha",
        "0.05",
        "--bins",
        "3",
        "--out",
        str(out_csv),
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert out_csv.exists(), "IHW results CSV not created"
    assert out_csv.stat().st_size > 0, "IHW results CSV empty"
