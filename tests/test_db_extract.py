import hashlib
import io
import tarfile
from dataclasses import replace
from pathlib import Path

import pytest

from phu.db import bundle as db_bundle
from phu.db.bundle import BundleMismatch, extract_bundle
from phu.db.fetch import DatabaseError
from phu.db.manifest import DATABASE_NAMES, MANIFEST_MEMBER
from phu.db.archive import ArchiveError, UnsafeArchive, check_tar_members, extract_member

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "phu_db"
PFAM_HMM = "phu-db/pfam/profiles.hmm"
AVG_MEMBERS = tuple(
    f"phu-db/avg/{name}"
    for name in (
        "filters_amg.tsv",
        "filters_apg.tsv",
        "filters_areg.tsv",
        "positive_amg.tsv",
        "positive_apg.tsv",
        "positive_areg.tsv",
        "v_scores.tsv",
    )
)


@pytest.fixture
def writes(monkeypatch):
    """Members extract_bundle wrote, recorded through its one write path."""
    written = []
    real = db_bundle.extract_member

    def recording(tar, member, dest, **kwargs):
        written.append(member.name)
        real(tar, member, dest, **kwargs)

    monkeypatch.setattr(db_bundle, "extract_member", recording)
    return written


@pytest.fixture
def extract(tmp_path, fixture_pins, fixture_manifest, fixture_manifest_bytes):
    db_folder = tmp_path / "db"

    def run(tarball, *, pins=fixture_pins):
        return extract_bundle(
            tarball, db_folder, pins, fixture_manifest, fixture_manifest_bytes
        )

    run.db_folder = db_folder
    return run


def _temp_dirs(db_folder: Path) -> list[Path]:
    return sorted(db_folder.glob(".phu-db-tmp-*")) if db_folder.exists() else []


def _assert_matches_fixture(root: Path, members) -> None:
    files = {path.relative_to(root.parent).as_posix() for path in root.rglob("*") if path.is_file()}
    assert files == set(members)
    for name in members:
        if name == MANIFEST_MEMBER:
            expected = (FIXTURE_DIR / "phu-db.json").read_bytes()
        else:
            expected = (FIXTURE_DIR / "prepared" / name.removeprefix("phu-db/")).read_bytes()
        assert (root.parent / name).read_bytes() == expected, name


def test_fixture_pins_cover_the_full_bundle(fixture_pins, fixture_manifest, fixture_bundle):
    assert len(fixture_pins.expected_members) == 16
    assert fixture_pins.expected_members == fixture_manifest.expected_members()
    assert fixture_pins.bundle_version == "phu-db-0.2.0-fixture"
    data = fixture_bundle.read_bytes()
    assert (fixture_pins.size, fixture_pins.sha256) == (len(data), hashlib.sha256(data).hexdigest())


def test_fixture_bundle_is_deterministic(fixture_bundle, make_bundle):
    assert make_bundle().read_bytes() == fixture_bundle.read_bytes()
    with tarfile.open(fixture_bundle) as tar:
        members = tar.getmembers()
    assert [member.name for member in members] == sorted(member.name for member in members)
    assert all(member.isreg() and member.mtime == 0 for member in members)
    assert {(member.uid, member.gid) for member in members} == {(0, 0)}


def test_valid_bundle_extracts(extract, fixture_bundle, fixture_pins, fixture_manifest_bytes):
    result = extract(fixture_bundle)

    assert result.databases == DATABASE_NAMES
    assert result.tempdir.parent == extract.db_folder
    assert result.tempdir.name.startswith(".phu-db-tmp-")
    assert result.root == result.tempdir / "phu-db"
    assert (result.root / "phu-db.json").read_bytes() == fixture_manifest_bytes
    _assert_matches_fixture(result.root, fixture_pins.expected_members)


def test_extracted_files_get_current_times_not_the_archive_epoch(extract, fixture_bundle):
    result = extract(fixture_bundle)
    assert (result.root / "pfam" / "profiles.hmm").stat().st_mtime > 0


def test_bundle_missing_a_whole_database_is_accepted(extract, make_bundle, fixture_pins):
    result = extract(make_bundle(drop=AVG_MEMBERS))

    assert result.databases == DATABASE_NAMES - {"avg"}
    assert not (result.root / "avg").exists()
    _assert_matches_fixture(result.root, fixture_pins.expected_members - set(AVG_MEMBERS))


def test_bundle_with_only_the_manifest_is_accepted(extract, make_bundle):
    everything = set(make_bundle.files()) - {MANIFEST_MEMBER}
    assert extract(make_bundle(drop=everything)).databases == frozenset()


def _symlink(make_bundle):
    return make_bundle(
        drop=[PFAM_HMM],
        extra=[(make_bundle.member(PFAM_HMM, tarfile.SYMTYPE, linkname="/etc/passwd"), None)],
    )


def _hardlink(make_bundle):
    link = make_bundle.member(PFAM_HMM, tarfile.LNKTYPE, linkname="phu-db/kofam/profiles.hmm")
    return make_bundle(drop=[PFAM_HMM], extra=[(link, None)])


def _device(make_bundle):
    return make_bundle(
        drop=[PFAM_HMM], extra=[(make_bundle.member(PFAM_HMM, tarfile.CHRTYPE), None)]
    )


def _fifo(make_bundle):
    return make_bundle(
        drop=[PFAM_HMM], extra=[(make_bundle.member(PFAM_HMM, tarfile.FIFOTYPE), None)]
    )


def _directory(make_bundle):
    return make_bundle(extra=[(make_bundle.member("phu-db/pfam", tarfile.DIRTYPE), None)])


def _member_with(name):
    return lambda make_bundle: make_bundle(extra=[(make_bundle.member(name), b"x\n")])


def _manifest_byte_flipped(make_bundle):
    data = bytearray(make_bundle.files()[MANIFEST_MEMBER])
    data[-2] ^= 0x01
    return make_bundle(replace={MANIFEST_MEMBER: bytes(data)})


@pytest.mark.parametrize(
    "build, error, message",
    [
        pytest.param(
            lambda b: b(drop=["phu-db/kofam/ko_list.tsv"]), UnsafeArchive, "part of the kofam",
            id="partial-database",
        ),
        pytest.param(
            lambda b: b(drop=[MANIFEST_MEMBER]), UnsafeArchive, "no phu-db/phu-db.json",
            id="missing-manifest",
        ),
        pytest.param(
            _member_with("phu-db/pfam/extra.hmm"), UnsafeArchive, "unexpected archive members",
            id="extra-member",
        ),
        pytest.param(
            lambda b: b(extra=[(b.member(PFAM_HMM), b"HMMER3/f\n")]), UnsafeArchive,
            "more than once", id="duplicate-member",
        ),
        pytest.param(_symlink, UnsafeArchive, "symbolic link", id="symlink"),
        pytest.param(_hardlink, UnsafeArchive, "hard link", id="hardlink"),
        pytest.param(_device, UnsafeArchive, "device", id="device"),
        pytest.param(_fifo, UnsafeArchive, "fifo", id="fifo"),
        pytest.param(_directory, UnsafeArchive, "directory", id="directory-member"),
        pytest.param(_member_with("phu-db/../evil.tsv"), UnsafeArchive, "unsafe", id="dotdot"),
        pytest.param(_member_with("../evil.tsv"), UnsafeArchive, "unsafe", id="leading-dotdot"),
        pytest.param(_member_with("/tmp/evil.tsv"), UnsafeArchive, "unsafe", id="absolute"),
        pytest.param(_member_with("phu-db\\evil.tsv"), UnsafeArchive, "unsafe", id="backslash"),
        pytest.param(_member_with("phu-db//evil.tsv"), UnsafeArchive, "unsafe", id="empty-segment"),
        pytest.param(_manifest_byte_flipped, BundleMismatch, "differs", id="manifest-byte"),
        pytest.param(
            lambda b: b(replace={MANIFEST_MEMBER: b.files()[MANIFEST_MEMBER] + b"\n"}),
            BundleMismatch, "differs", id="manifest-longer",
        ),
        pytest.param(
            lambda b: b(replace={PFAM_HMM: b"HMMER3/f\n"}), BundleMismatch, "bytes",
            id="prepared-size",
        ),
    ],
)
def test_bad_bundle_is_rejected_before_writing(extract, make_bundle, writes, build, error, message):
    tarball = build(make_bundle)

    with pytest.raises(error, match=message):
        extract(tarball)

    assert writes == []
    assert _temp_dirs(extract.db_folder) == []


def test_errors_are_database_errors():
    assert issubclass(UnsafeArchive, DatabaseError)
    assert issubclass(ArchiveError, DatabaseError)
    assert issubclass(BundleMismatch, DatabaseError)


@pytest.mark.parametrize(
    "content",
    [b"", b"not a tarball", b"\x1f\x8b\x08\x00garbage"],
    ids=["empty", "text", "bad-gzip"],
)
def test_unreadable_bundle_raises_archive_error(extract, tmp_path, writes, content):
    tarball = tmp_path / "broken.tar.gz"
    tarball.write_bytes(content)

    with pytest.raises(ArchiveError):
        extract(tarball)
    assert writes == []
    assert _temp_dirs(extract.db_folder) == []


def test_truncated_bundle_raises_archive_error(extract, fixture_bundle, tmp_path):
    tarball = tmp_path / "truncated.tar.gz"
    tarball.write_bytes(fixture_bundle.read_bytes()[:-200])

    with pytest.raises(ArchiveError):
        extract(tarball)
    assert _temp_dirs(extract.db_folder) == []


@pytest.mark.parametrize(
    "change",
    [
        {"manifest_sha256": "0" * 64},
        {"bundle_version": "phu-db-0.0.0"},
        {"expected_members": frozenset({MANIFEST_MEMBER})},
    ],
    ids=["manifest-sha256", "bundle-version", "members"],
)
def test_pins_that_disagree_with_the_vendored_manifest_are_refused(
    extract, fixture_bundle, fixture_pins, writes, change
):
    with pytest.raises(BundleMismatch):
        extract(fixture_bundle, pins=replace(fixture_pins, **change))
    assert writes == []
    assert _temp_dirs(extract.db_folder) == []


def test_failure_during_extraction_removes_the_temp_dir(extract, fixture_bundle, monkeypatch):
    calls = []
    real = db_bundle.extract_member

    def failing(tar, member, dest, **kwargs):
        calls.append(member.name)
        if len(calls) == 3:
            raise OSError("disk full")
        real(tar, member, dest, **kwargs)

    monkeypatch.setattr(db_bundle, "extract_member", failing)

    with pytest.raises(ArchiveError, match="disk full"):
        extract(fixture_bundle)
    assert len(calls) == 3
    assert _temp_dirs(extract.db_folder) == []


def test_bundle_swapped_between_scan_and_extraction_is_refused(
    extract, fixture_bundle, make_bundle, monkeypatch
):
    swapped = make_bundle(drop=AVG_MEMBERS)
    real_scan = db_bundle._scan

    def scan_then_swap(tarball, *args):
        result = real_scan(tarball, *args)
        tarball.write_bytes(swapped.read_bytes())
        return result

    tarball = make_bundle()
    monkeypatch.setattr(db_bundle, "_scan", scan_then_swap)

    with pytest.raises(UnsafeArchive, match="changed while"):
        extract(tarball)
    assert _temp_dirs(extract.db_folder) == []


@pytest.mark.skipif(not hasattr(tarfile, "data_filter"), reason="tarfile has no data filter")
def test_extraction_uses_the_data_filter(extract, fixture_bundle, monkeypatch):
    filters = []
    real = tarfile.TarFile.extract

    def recording(self, member, *args, **kwargs):
        filters.append(kwargs.get("filter"))
        return real(self, member, *args, **kwargs)

    monkeypatch.setattr(tarfile.TarFile, "extract", recording)

    extract(fixture_bundle)
    assert filters == ["data"] * 16


def test_extraction_without_the_data_filter(extract, fixture_bundle, fixture_pins, monkeypatch):
    monkeypatch.delattr(tarfile, "data_filter")

    def no_tarfile_extract(*args, **kwargs):
        raise AssertionError("the fallback must not call TarFile.extract")

    monkeypatch.setattr(tarfile.TarFile, "extract", no_tarfile_extract)

    result = extract(fixture_bundle)
    assert result.databases == DATABASE_NAMES
    _assert_matches_fixture(result.root, fixture_pins.expected_members)


def test_unsafe_bundle_is_still_rejected_without_the_data_filter(
    extract, make_bundle, writes, monkeypatch
):
    monkeypatch.delattr(tarfile, "data_filter")
    with pytest.raises(UnsafeArchive, match="symbolic link"):
        extract(_symlink(make_bundle))
    assert writes == []
    assert _temp_dirs(extract.db_folder) == []


# check_tar_members without an expected set: the contract's extract verb.


def _tar(entries) -> tarfile.TarFile:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for info, data in entries:
            if data is not None:
                info.size = len(data)
            archive.addfile(info, io.BytesIO(data) if data is not None else None)
    buffer.seek(0)
    return tarfile.open(fileobj=buffer, mode="r")


def _info(name, type=tarfile.REGTYPE, linkname=""):
    info = tarfile.TarInfo(name)
    info.type = type
    info.linkname = linkname
    return info


def test_upstream_archive_may_hold_directories():
    tar = _tar(
        [(_info("profiles", tarfile.DIRTYPE), None), (_info("profiles/K00001.hmm"), b"HMMER3/f\n")]
    )
    assert [member.name for member in check_tar_members(tar, None)] == [
        "profiles",
        "profiles/K00001.hmm",
    ]


@pytest.mark.parametrize(
    "entry",
    [
        (_info("profiles/K1.hmm", tarfile.SYMTYPE, "/etc/passwd"), None),
        (_info("profiles/K1.hmm", tarfile.LNKTYPE, "profiles/K0.hmm"), None),
        (_info("profiles/dev", tarfile.BLKTYPE), None),
        (_info("profiles/pipe", tarfile.FIFOTYPE), None),
        (_info("../profiles/K1.hmm"), b"x"),
        (_info("/profiles/K1.hmm"), b"x"),
    ],
    ids=["symlink", "hardlink", "block-device", "fifo", "dotdot", "absolute"],
)
def test_upstream_archive_rejects_unsafe_members(entry):
    tar = _tar([(_info("profiles/K0.hmm"), b"HMMER3/f\n"), entry])
    with pytest.raises(UnsafeArchive):
        check_tar_members(tar, None)


# extract_member refuses unsafe members itself, on both write paths.


@pytest.fixture(params=["data-filter", "fallback"])
def write_path(request, monkeypatch):
    if request.param == "data-filter":
        if not hasattr(tarfile, "data_filter"):
            pytest.skip("tarfile has no data filter")
    else:
        monkeypatch.delattr(tarfile, "data_filter", raising=False)
    return request.param


def _files_below(path: Path) -> set[str]:
    return {item.relative_to(path).as_posix() for item in path.rglob("*")}


@pytest.mark.parametrize(
    "entry, message",
    [
        ((_info("a/link", tarfile.SYMTYPE, "../b"), None), "symbolic link"),
        ((_info("a/link", tarfile.LNKTYPE, "a/x"), None), "hard link"),
        ((_info("a/pipe", tarfile.FIFOTYPE), None), "fifo"),
        ((_info("../x"), b"x"), "unsafe"),
        ((_info("a/x/../../y"), b"x"), "unsafe"),
        ((_info("/x"), b"x"), "unsafe"),
    ],
    ids=["symlink", "hardlink", "fifo", "leading-dotdot", "inner-dotdot", "absolute"],
)
def test_extract_member_refuses_unchecked_unsafe_members(tmp_path, write_path, entry, message):
    dest = tmp_path / "dest"
    dest.mkdir()
    tar = _tar([entry])

    with pytest.raises(UnsafeArchive, match=message):
        extract_member(tar, tar.getmembers()[0], dest)
    assert _files_below(tmp_path) == {"dest"}


def test_extract_member_regular_only_refuses_directories(tmp_path, write_path):
    tar = _tar([(_info("profiles", tarfile.DIRTYPE), None)])
    with pytest.raises(UnsafeArchive, match="directory"):
        extract_member(tar, tar.getmembers()[0], tmp_path, regular_only=True)
    assert _files_below(tmp_path) == set()


def test_extract_member_writes_directories_and_files(tmp_path, write_path):
    tar = _tar(
        [(_info("profiles", tarfile.DIRTYPE), None), (_info("profiles/K1.hmm"), b"HMMER3/f\n")]
    )
    for member in tar.getmembers():
        extract_member(tar, member, tmp_path)
    assert (tmp_path / "profiles" / "K1.hmm").read_bytes() == b"HMMER3/f\n"


def test_extract_member_does_not_write_through_an_existing_link(tmp_path, write_path):
    outside = tmp_path / "outside.txt"
    outside.write_bytes(b"keep\n")
    dest = tmp_path / "dest"
    dest.mkdir()
    (dest / "x").symlink_to(outside)
    tar = _tar([(_info("x"), b"evil\n")])

    # The data filter refuses the resolved path; the fallback's O_NOFOLLOW fails with ELOOP.
    with pytest.raises((tarfile.TarError, OSError)):
        extract_member(tar, tar.getmembers()[0], dest)
    assert outside.read_bytes() == b"keep\n"
