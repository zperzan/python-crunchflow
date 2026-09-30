import pytest

from crunchflow.cli.main import main


def test_clear_output_command(tmp_path):
    (tmp_path / "conc1.tec").write_text("")

    assert main(["clear-output", "--folder", str(tmp_path)]) == 0
    assert not (tmp_path / "conc1.tec").exists()


def test_clear_output_command_dryrun(tmp_path):
    (tmp_path / "conc1.tec").write_text("")

    assert main(["clear-output", "--folder", str(tmp_path), "--dry-run"]) == 0
    assert (tmp_path / "conc1.tec").exists()


def test_clear_output_command_suffixes(tmp_path):
    (tmp_path / "conc1.tec").write_text("")
    (tmp_path / "conc1.out").write_text("")

    main(["clear-output", "--folder", str(tmp_path), "--suffixes", ".out"])

    assert sorted(p.name for p in tmp_path.iterdir()) == ["conc1.tec"]


def test_clear_output_command_quiet(tmp_path, capsys):
    (tmp_path / "conc1.tec").write_text("")

    main(["clear-output", "--folder", str(tmp_path)])

    assert capsys.readouterr().out == ""


def test_no_command_prints_help(capsys):
    assert main([]) == 1
    assert "clear-output" in capsys.readouterr().out


def test_version():
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])

    assert excinfo.value.code == 0


if __name__ == "__main__":
    pytest.main()
