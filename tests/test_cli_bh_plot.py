import subprocess, sys, tempfile, csv, pathlib

def test_bh_plot_cli_runs_and_writes_png():
    tmp = tempfile.TemporaryDirectory()
    tmpdir = pathlib.Path(tmp.name)
    csv_path = tmpdir / "p.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pval"])
        w.writerows([[0.005], [0.012], [0.03], [0.045], [0.055], [0.2], [0.4], [0.6], [0.8], [0.9]])
    out_csv = tmpdir / "out.csv"
    out_png = tmpdir / "bh.png"
    cmd = [
        sys.executable,
        "-m",
        "mc_cc0.cli",
        "correct",
        "--csv",
        str(csv_path),
        "--col",
        "pval",
        "--alpha",
        "0.05",
        "--method",
        "fdr_bh",
        "--bh-plot",
        "--plot-out",
        str(out_png),
        "--out",
        str(out_csv),
    ]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert out_csv.exists(), "results CSV not created"
    assert out_png.exists(), "BH plot PNG not created"
