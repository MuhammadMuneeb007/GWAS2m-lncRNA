"""Resumable HTTP(S) downloads with retries. No wget/curl/aria2 required."""

from __future__ import annotations

import os
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

USER_AGENT = "gwas2m-lncrna/0.1 (+https://github.com/MuhammadMuneeb007/GWAS2m-lncRNA)"
CHUNK = 4 * 1024 * 1024


def http_session(retries: int = 6) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=retries, connect=retries, read=retries, status=retries,
        backoff_factor=1.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "HEAD"],
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.mount("http://", HTTPAdapter(max_retries=retry))
    session.headers.update({"User-Agent": USER_AGENT})
    return session


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def download(url: str, target: Path, session: requests.Session | None = None,
             attempts: int = 10, label: str | None = None) -> str:
    """Download url to target via target.part, resuming partial files.

    Returns "SKIPPED_ALREADY_COMPLETE" if target already exists, else "COMPLETE".
    The final file only appears after the transfer finished, so its existence
    means it is complete.
    """
    target = Path(target)
    if target.is_file() and target.stat().st_size > 0:
        return "SKIPPED_ALREADY_COMPLETE"
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    session = session or http_session()
    name = label or target.name

    for attempt in range(1, attempts + 1):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            with session.get(url, headers=headers, stream=True, timeout=(30, 120)) as response:
                if response.status_code == 416:  # already fully downloaded
                    break
                response.raise_for_status()
                if offset and response.status_code != 206:
                    offset = 0  # server ignored Range: start again
                total = response.headers.get("Content-Length")
                total = int(total) + offset if total and total.isdigit() else None
                mode = "ab" if offset else "wb"
                done = offset
                last_report = time.monotonic()
                print(f"[download] {name}: {'resuming at ' + _human(offset) if offset else 'starting'}"
                      + (f" of {_human(total)}" if total else ""), flush=True)
                with open(partial, mode) as handle:
                    for chunk in response.iter_content(chunk_size=CHUNK):
                        if not chunk:
                            continue
                        handle.write(chunk)
                        done += len(chunk)
                        if time.monotonic() - last_report > 30:
                            share = f" ({100 * done / total:.0f}%)" if total else ""
                            print(f"[download] {name}: {_human(done)}{share}", flush=True)
                            last_report = time.monotonic()
                if total is not None and done < total:
                    raise IOError(f"incomplete transfer ({done} of {total} bytes)")
            break
        except (requests.RequestException, IOError) as exc:
            if attempt == attempts:
                raise RuntimeError(f"Download failed after {attempts} attempts: {url}\n{exc}") from exc
            wait = min(60, 5 * attempt)
            print(f"[download] {name}: {exc}; retrying in {wait}s ({attempt}/{attempts})", flush=True)
            time.sleep(wait)

    if not partial.exists() or partial.stat().st_size == 0:
        raise RuntimeError(f"Download produced no data: {url}")
    os.replace(partial, target)
    print(f"[download] {name}: done ({_human(target.stat().st_size)})", flush=True)
    return "COMPLETE"
