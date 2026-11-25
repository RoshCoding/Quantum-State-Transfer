import numpy as np

from qst.cli import main

def test_cli_perfect_peak(capsys, tmp_path):
    out = tmp_path / "p.txt"
    main(["--chain", "perfect", "--dt", "0.01", "--steps", "400", "-o", str(out)])
    assert "peak" in capsys.readouterr().out
    P = np.loadtxt(out)
    assert P.shape == (400, 64)
    assert P[:, -1].max() > 0.9999
