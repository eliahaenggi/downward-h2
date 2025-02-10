#! /usr/bin/env python

import os

import custom_parser
import project

from downward import suites
from downward.cached_revision import CachedFastDownwardRevision
from downward.experiment import FastDownwardAlgorithm, FastDownwardRun
from lab.experiment import Experiment

REPO = project.get_repo_base()
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
REVISION_CACHE = (
    os.environ.get("DOWNWARD_REVISION_CACHE") or project.DIR / "data" / "revision-cache"
)
SUITE = ["depot:p02.pddl", "driverlog:p02.pddl", "gripper:prob01.pddl", "gripper:prob02.pddl", "satellite:p01-pfile1.pddl", "trucks:p01.pddl", "logistics00:probLOGISTICS-4-0.pddl",  "logistics00:probLOGISTICS-5-0.pddl", "movie:prob01.pddl", "satellite:p01-pfile1.pddl"]

ENV = project.LocalEnvironment(processes=2)

CONFIGS = [
    ("astar-h2", ["--search", "astar(h2())"]),
    ("astar-hmax-pi-m", ["--search", "astar(hmax(pi_m_compilation = True))"]),
]

BUILD_OPTIONS = []
DRIVER_OPTIONS = [
    "--validate",
    "--overall-time-limit",
    "5m",
    "--overall-memory-limit",
    "3584M",
]
# Pairs of revision identifier and optional revision nick.
REV_NICKS = [
    ("0075e6bf47f610a35e9d6303aed87e69bf4df4e7", "pi_m_compilation")
]

ATTRIBUTES = [
    "error",
    "run_dir",
    "search_start_time",
    "search_start_memory",
    "total_time",
    "h_values",
    "coverage",
    "expansions",
    "memory",
    project.EVALUATIONS_PER_TIME,
]

exp = Experiment(environment=ENV)
for rev, rev_nick in REV_NICKS:
    cached_rev = CachedFastDownwardRevision(REVISION_CACHE, REPO, rev, BUILD_OPTIONS)
    cached_rev.cache()
    exp.add_resource("", cached_rev.path, cached_rev.get_relative_exp_path())
    for config_nick, config in CONFIGS:
        algo_name = f"{rev_nick}-{config_nick}" if rev_nick else config_nick

        for task in suites.build_suite(BENCHMARKS_DIR, SUITE):
            algo = FastDownwardAlgorithm(
                algo_name,
                cached_rev,
                DRIVER_OPTIONS,
                config,
            )
            run = FastDownwardRun(exp, algo, task)
            exp.add_run(run)

exp.add_parser(project.FastDownwardExperiment.EXITCODE_PARSER)
exp.add_parser(project.FastDownwardExperiment.TRANSLATOR_PARSER)
exp.add_parser(project.FastDownwardExperiment.SINGLE_SEARCH_PARSER)
exp.add_parser(custom_parser.get_parser())
exp.add_parser(project.FastDownwardExperiment.PLANNER_PARSER)

exp.add_step("build", exp.build)
exp.add_step("start", exp.start_runs)
exp.add_step("parse", exp.parse)
exp.add_fetcher(name="fetch")

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES,
    filter=[project.add_evaluations_per_time, project.group_domains],
)

exp.run_steps()
