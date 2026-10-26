# Track AgentCore sessions

See [SKILL.md](../SKILL.md) for concepts and the shared skill-version check — this
covers what's specific to LaminHub's AgentCore **Python** chat mode.

The runtime supplies the instance,
authenticated user or bot, chat ID, and an in-memory snapshot of the saved chat
plus the current turn's public messages. Do not look for Claude/Codex log files,
run `lamin connect` or `login`, install packages, or discover a dev-dir.

## Start

An explicit request to track is consent. Otherwise ask once: **Track this session
in LaminDB?** Offer **Track** and **Do not track** without a default. Wait for the
user's choice; plain chat text is appropriate when no interactive question tool
exists. If they decline, do not track.

After they agree, run the shared version check through the existing in-process
CLI helper in `execute_python`:

```python
lamin("--skill-version")
```

Read `SKILL.md` with `read_lamindb_guide(section="skill")` and compare the result
with its `metadata.version`, as the shared check describes. If they differ,
stop tracking and report that the runtime's bundled
skill and package need to be updated together; do not install desktop skills or
upgrade packages inside the worker. Then start tracking:

```python
lamin("track", "agentcore", "--name", "A short description of the task")
```

The fixed Transform is `__agentcore__`. One Run represents this chat in this
instance. On a follow-up, call track again before the task to resume the same Run.
Its identity survives a Python restart or runtime replacement; no local marker
files are required. Commands use the execution user's normal write permissions.
If tracking fails, report it and continue the user's task without claiming it was
tracked.

## Inputs, outputs and lineage

Scripts/notebooks remain Transforms and should self-track with `ln.track()` and
`ln.finish()`. Track sets `LAMIN_INITIATED_BY_RUN_UID` in the Python worker so
synchronous child processes inherit it. In the hosted worker, run scripts with
`runpy.run_path(path, run_name="__main__")` to reuse its authenticated interpreter;
a new shell process does not inherit the in-memory connection. Do not start
background processes. Finish links these child Transforms to the
agent Run using the same helper as the desktop harnesses.

For direct Python work, link inputs and outputs to the agent Run explicitly:

```python
from lamin_cli.agents.agentcore import get_agentcore_run
run = get_agentcore_run()
source = ln.Artifact.get("<input-uid>", is_run_input=run)
df = source.load(is_run_input=False)  # already linked to the agent Run above
# Compute result from df, then save it as this Run's output.
ln.Artifact.from_dataframe(result, key="results/summary.csv", run=run).save()
```

`display_table`, `display_figure` and `attach_file` publish chat results; they do
not automatically register LaminDB Artifacts. Only register outputs requested by
the user. Tracking does not create a branch or open a change request; use the
existing selected instance/branch, and make branch changes only as requested.

## Finish

**Write the user-facing summary as public assistant text before the finish tool
executes**, alongside the structured `execute_python` call. Do not put the summary
inside `print()`: Python output is not part of the saved public conversation.
Then execute:

```python
lamin("finish")
```

Finish immediately saves `index.html` and `trace.jsonl` using the shared report
renderer and writer. It includes saved user/assistant messages and the current
turn's public commentary, including the summary. It excludes this finish call,
later replies, hidden reasoning, and tool traces that Hub does not save. No wait
for desktop transcript files is needed. A later finish replaces the same Run's
report with the new snapshot. A repeated finish is safe.

If a turn is cancelled or crashes before finish, it has no new completed report.
On a later turn, resume tracking and finish explicitly. A report save failure
must be reported; do not claim success or close the Run until saving succeeds.
