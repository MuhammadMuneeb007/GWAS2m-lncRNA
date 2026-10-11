"""Run one step for many studies: in this process, never through a scheduler.

A study that fails is recorded as FAILED with its reason and the run moves on
to the next study; downstream steps simply skip it.
"""

from __future__ import annotations

import contextlib
import sys
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable

from gwas2m_lncrna.common import DONE, banner
from gwas2m_lncrna.project import Project

# A study function returns (status, reason, extra_dict).
StudyFunction = Callable[[Project, dict, bool], tuple]


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def _call(func: StudyFunction, project: Project, study: dict, force: bool) -> tuple:
    try:
        return func(project, study, force)
    except Exception as exc:  # one study must never stop the others
        traceback.print_exc()
        return "FAILED", f"{type(exc).__name__}: {exc}"[:500], {}


def _execute(project: Project, step: str, func: StudyFunction, study: dict, force: bool,
             echo: bool, capture: bool = True) -> tuple[str, str, str, dict]:
    accession = study["STUDY_ACCESSION"]
    log_path = project.log_file(step, accession)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if not capture:  # thread pools share sys.stdout, so output goes to the console only
        status, reason, extra = _call(func, project, study, force)
    else:
        with open(log_path, "a", encoding="utf-8") as log:
            out = _Tee(sys.__stdout__, log) if echo else log
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
                status, reason, extra = _call(func, project, study, force)
    project.set_status(accession, step, status, reason, log=str(log_path), **(extra or {}))
    return accession, status, reason, extra or {}


def run_for_studies(project: Project, step: str, title: str, func: StudyFunction,
                    studies: list[dict], force: bool = False, jobs: int = 1,
                    threads_not_processes: bool = False) -> Counter:
    banner(f"{title} | {project.phenotype} | {project.ancestry_code} | {len(studies)} studies")
    tally = Counter()
    if not studies:
        print("No studies are ready for this step.")
        return tally

    def report(result):
        accession, status, reason, _ = result
        tally[status] += 1
        mark = "OK  " if status in DONE else status[:4]
        print(f"[{mark}] {accession:<15} {status}" + (f" - {reason}" if reason and status not in DONE else ""),
              flush=True)

    if jobs <= 1:
        for study in studies:
            report(_execute(project, step, func, study, force, echo=True))
    else:
        pool_type = ThreadPoolExecutor if threads_not_processes else ProcessPoolExecutor
        with pool_type(max_workers=jobs) as pool:
            futures = [pool.submit(_execute, project, step, func, study, force, False,
                                   not threads_not_processes) for study in studies]
            for future in as_completed(futures):
                report(future.result())
        print(f"Per-study logs: {project.root / 'logs' / step}")
    print("Summary: " + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    return tally


def ready_for(project: Project, studies: list[dict], upstream: str,
              accept: set[str] | None = None) -> list[dict]:
    """Studies whose upstream step finished (COMPLETE / SKIPPED_ALREADY_COMPLETE)."""
    accept = accept or DONE
    return [s for s in studies if project.step_status(s["STUDY_ACCESSION"], upstream) in accept]
