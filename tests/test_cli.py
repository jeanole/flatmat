import json
import os

from parts import Part, make_solid, write_step

from flatmat.cli import main


def test_cli_unfold(tmp_path, capsys):
    step = tmp_path / "u.step"
    write_step(make_solid(Part("u", 1.5, 80, [30, 60, 30], [(90, 1.5), (90, 1.5)])), str(step))
    rc = main(["unfold", str(step), "-k", "0.4", "-o", str(tmp_path / "u.dxf"), "--table", str(tmp_path / "u.json")])
    assert rc == 0
    assert "2 bends" in capsys.readouterr().out
    assert os.path.exists(tmp_path / "u.dxf")
    assert json.loads((tmp_path / "u.json").read_text())["k_model"] == "K=0.4"


def test_cli_refusal(tmp_path, capsys):
    p = tmp_path / "bad.step"
    p.write_text("nope")
    assert main(["unfold", str(p)]) == 2
    assert json.loads(capsys.readouterr().err.strip().splitlines()[-1])["code"] == "READ_FAILED"
