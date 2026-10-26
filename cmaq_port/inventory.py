"""Read-only inventory of a local CMAQ Git checkout."""

import hashlib
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "--no-pager", "-C", str(root), *arguments],
        check=True, capture_output=True, text=True, encoding="utf-8", timeout=60,
    ).stdout


def git_blob(root: Path, blob: str) -> bytes:
    return subprocess.run(
        ["git", "--no-pager", "-C", str(root), "cat-file", "blob", blob],
        capture_output=True, check=True, timeout=60,
    ).stdout


def inventory(root: Path) -> dict[str, object]:
    root = root.resolve(strict=True)
    if not (root / "CCTM" / "src" / "driver" / "cmaq_main.F").is_file():
        raise ValueError(f"Not a recognized CMAQ source checkout: {root}")
    commit = git(root, "rev-parse", "HEAD").strip()
    status = git(root, "status", "--porcelain=v1", "--untracked-files=all")
    index = git(root, "ls-files", "--stage", "-t", "-z")
    entries = []
    for entry in index.split("\0"):
        if not entry:
            continue
        metadata, relative = entry.split("\t", 1)
        flag, mode, blob, stage = metadata.split()
        if stage != "0":
            raise ValueError(f"Unmerged Git index entry: {relative}")
        entries.append((relative, flag, mode, blob))
    files: list[dict[str, object]] = []
    groups: Counter[str] = Counter()
    excluded: list[str] = []
    links: list[str] = []
    for relative, flag, mode, blob in sorted(entries):
        path = root / relative
        groups[relative.split("/")[0]] += 1
        if mode == "120000":
            stored = git_blob(root, blob)
            files.append({
                "path": relative, "bytes": len(stored),
                "sha256": hashlib.sha256(stored).hexdigest(), "git_blob": blob,
                "content_origin": "git_index_link_target", "link_target": stored.decode("utf-8"),
                "available_in_worktree": path.exists() or path.is_symlink(),
                "materialized_as_symlink": path.is_symlink(),
            })
            links.append(relative)
            continue
        if mode not in ("100644", "100755") or path.is_symlink():
            raise ValueError(f"Inventory does not follow tracked symlinks: {path}")
        if not path.is_file():
            if flag != "S":
                raise ValueError(f"Tracked file missing or not a regular file: {path}")
            stored = git_blob(root, blob)
            files.append({
                "path": relative, "bytes": len(stored),
                "sha256": hashlib.sha256(stored).hexdigest(), "git_blob": blob,
                "content_origin": "git_index_blob", "available_in_worktree": False,
            })
            excluded.append(relative)
            continue
        resolved = path.resolve(strict=True)
        if not resolved.is_relative_to(root):
            raise ValueError(f"Tracked file resolves outside the source checkout: {path}")
        files.append({
            "path": relative, "bytes": path.stat().st_size, "sha256": sha256(path),
            "git_blob": blob, "content_origin": "working_tree", "available_in_worktree": True,
        })
    if commit != git(root, "rev-parse", "HEAD").strip() or status != git(
        root, "status", "--porcelain=v1", "--untracked-files=all",
    ):
        raise ValueError("Checkout changed during inventory; retry when edits have stopped.")
    if index != git(root, "ls-files", "--stage", "-t", "-z"):
        raise ValueError("Git index changed during inventory; retry when edits have stopped.")
    return {
        "schema_version": 1, "kind": "source_inventory", "source_root": str(root),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "commit": commit, "branch": git(root, "branch", "--show-current").strip(),
        "git_status": status, "clean": status == "",
        "readme_identity": (root / "README.md").read_text(encoding="utf-8").splitlines()[0],
        "tracked_file_count": len(files), "top_level_counts": dict(sorted(groups.items())),
        "files": files,
        "excluded_worktree_paths": excluded, "worktree_complete": not excluded,
        "tracked_link_paths": links,
        "warnings": (
            [f"{len(excluded)} skip-worktree files unavailable on disk; hashed Git index blobs instead."]
            if excluded else []
        ) + ([f"{len(links)} tracked symbolic links inventoried as Git link targets, not followed."] if links else []),
        "scope": "Tracked file hashes with explicit content origin; not a model build or numerical validation.",
    }


def save_report(report: dict[str, object], kind: str, *, project: Path = PROJECT_ROOT) -> Path:
    folder = project / "runs" / (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + f"-{kind}-" + uuid4().hex[:8]
    )
    folder.mkdir(parents=True)
    destination = folder / f"{kind}.json"
    destination.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return destination
