import stat
import sys

import pytest

from phu.core import exec as phu_exec


def _fake_tool(directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    tool = directory / name
    tool.write_text("#!/bin/sh\n")
    tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    return tool


def test_executable_found_beside_python_when_not_on_path(tmp_path, monkeypatch):
    env_bin = tmp_path / "tool-env" / "bin"
    tool = _fake_tool(env_bin, "vclust")
    monkeypatch.setattr(sys, "executable", str(env_bin / "python"))
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))

    assert phu_exec._executable(["vclust", "vclust.py"]) == str(tool)


def test_executable_falls_back_to_path(tmp_path, monkeypatch):
    path_dir = tmp_path / "on-path"
    tool = _fake_tool(path_dir, "vclust.py")
    monkeypatch.setattr(sys, "executable", str(tmp_path / "other" / "python"))
    monkeypatch.setenv("PATH", str(path_dir))

    assert phu_exec._executable(["vclust", "vclust.py"]) == str(tool)


def test_executable_missing_everywhere(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python"))
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(phu_exec.CmdNotFound, match="beside Python or in PATH"):
        phu_exec._executable(["vclust"])
