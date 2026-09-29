---
name: lamindb
description: "MUST invoke before responding to ANY message — including greetings, small talk, trivial math, anything that looks unrelated. Not a judgment call: never skip it, never decide a message is too trivial. Tracks this session in LaminDB as a Transform + Run. If you are about to respond without invoking this first, stop — that is already a mistake."
metadata:
  version: "1.1"
---

# LaminDB

Official LaminDB skill to write code with best practices, keeping up to date with new versions and features.

Run `lamin --skill-version` only after the user agrees to track, from the same environment that provides `lamin`, and compare the printed value with this file's `metadata.version`. If they differ, stop and tell the user this skill is stale: remove `.agents/skills/lamindb` (and `.claude/skills/lamindb` if present), then run `uvx library-skills` (add `--claude` for Claude Code). Do not continue tracking on the old skill.

> Agent tracking requires lamindb >= 2.11.

## Concepts

- **Transform**: code, not data. Each coding agent/harness has its own fixed Transform for the whole project representing that harness's sessions (the "agent run") — internally keyed `__claudecode__` for Claude Code, `__copilot__` for Copilot, and `__cursor__` for Cursor. **These are internal database identifiers only, never CLI arguments or command names** — the actual commands are `lamin track claude`, `lamin track copilot`, and `lamin track cursor`; see your harness's reference file for the exact syntax. **Any script you write to accomplish the user's task (`.py`/`.ipynb`/`.R`/`.Rmd`/`.qmd`) is its own separate Transform, tracked automatically the moment it runs** — never save a script as a plain Artifact. Getting this backwards destroys the lineage from script to the data it produced, which is the entire point of LaminDB.
- **Run**: an execution. Each harness conversation gets one Run of your harness's fixed Transform (the **agent run**) per LaminDB instance. If the user sends a follow-up after you already ran `lamin finish`, run `lamin track <agent>` again: it resumes that same agent Run using the harness session ID (Cursor uses its generated session marker to resolve the conversation in its local chat database) instead of creating a duplicate, and the next `lamin finish` replaces its report with the complete conversation. Every script you write self-tracks its _own_ Run the instant it executes, linked back to the agent run via `initiated_by_run` — see "Self-tracking scripts" below. You never construct the script's Transform/Run by hand from outside.
- **Two distinct link fields — do not conflate them**: `Run.initiated_by_run` (on the _Run_ model) says "this execution was triggered by that other run" — it only exists once a script actually executes, and renders in its own "This run initiated" panel in the UI, not as an output. `Transform.run` (on the _Transform_ model, separate field) says "this piece of code was authored/produced during that run" — it's what makes a script show up in the agent run's **Output** column (alongside artifacts), the way a plain output file does. `ln.track()` never sets `Transform.run` on its own — `lamin finish` stamps it explicitly at session close, so a script counts as a session output even if it's the _only_ thing produced.
- **Never save a script as a plain Artifact.** Scripts (`.py`/`.ipynb`/`.R`/`.Rmd`/`.qmd`) must use `ln.track()` inside them. If you call `ln.Artifact("script.py").save()` you destroy the lineage between the code and the data it produced — that is the entire point of LaminDB and must never happen.
- **run.report**: rendered HTML of the transcript, saved as an Artifact and linked to the agent run.
- **Artifact**: data only — output files (csv, txt, images, fasta, etc.). A script's own `ln.Artifact(path).save()` calls (no `run=` needed) auto-attach to that script's own run. Only files you create directly, with no script involved, get attached to the agent run manually.
- **Always pass a meaningful, identity-specific `key`** when saving an Artifact (e.g. `key="datasets/ataqseq_counts.csv"`). Only reuse a key when saving a new version of the same dataset.
- **Create artifacts from in-memory objects when possible**: Prefer `ln.Artifact.from_*()` over writing objects to disk and then calling `ln.Artifact(path).save()`. For example, use `ln.Artifact.from_anndata()` for `AnnData` and `ln.Artifact.from_dataframe()` for pandas `DataFrame`. **Do not use `df.to_csv(...)` (or similar) as an intermediate step if a matching `from_*` constructor exists**. Use path-based `ln.Artifact(path).save()` only when the output is genuinely file-native (e.g. image, FASTA, binary export, or a format without a `from_*` helper). If you're unsure whether a `from_*` helper exists for an object type, run `help(ln.Artifact)` to inspect supported constructors.
- **When a script needs data that an earlier script in this workflow already produced, retrieve it from LaminDB — never read the local file path directly.** Use `ln.Artifact.get(key="...")` (the same key it was saved under) followed by `.load()`; this is what registers that artifact as this run's input and forms the lineage edge between the two scripts. Reading the file straight off disk produces the same result but leaves LaminDB with no record that the two scripts are connected, silently breaking the workflow's lineage graph.

## Self-tracking scripts and notebooks

Every script or notebook you write to do the user's actual task must instrument itself — this is what gives each output file a real lineage back to the exact code that produced it:

```python
import lamindb as ln
ln.track()
# if this step consumes an earlier step's output, retrieve it — never read the local file path directly:
input_artifact = ln.Artifact.get(key="<key used when it was saved>")
df = input_artifact.load()  # registers it as this run's input, forming the lineage edge
# ... the actual task ... e.g. process_data(df)
# for in-memory objects, prefer from_* constructors (no to_csv/to_parquet/any other intermediate write)
artifact = ln.Artifact.from_dataframe(df, key="<meaningful/folder/path>/output.csv", description="...")
artifact.save()
# use path-based save only for genuinely file-native outputs:
ln.Artifact("figure.png", key="<meaningful/folder/path>/figure.png", description="...").save()
ln.finish()
```

For Jupyter notebooks, `ln.track()` cannot auto-detect the notebook path when run via `nbconvert`, and without `new_run=True` it loads the latest run instead of creating a new one. Fix both by passing `path` and `new_run=True` explicitly in the first cell:

```python
import lamindb as ln
ln.track(path="notebook.ipynb", new_run=True)
```

When you actually **run** such a script or notebook, always set `LAMIN_INITIATED_BY_RUN_UID` and use the same project environment and LaminDB installation that created the agent Run and selected its branch. If execution fails because a task dependency such as `pandas` is missing, add only that dependency to the same environment and retry with the same interpreter. **Never switch to a temporary `uv run --with lamindb` environment**: another LaminDB installation may resolve different branch state and save the child Run and artifacts to the wrong branch. If the matching interpreter cannot import LaminDB, stop and ask the user to repair that environment. This links the script's/notebook's self-created Run back to the agent run automatically. On its own this only sets `Run.initiated_by_run`; `lamin finish` separately stamps `Transform.run` so the file also shows up as a session **output**, not just an "initiated" run. How you find your own `LAMIN_INITIATED_BY_RUN_UID` value differs by coding agent — see your harness's reference file (linked below) for the exact command, and run it exactly as shown. Writing your own simplified tracking (e.g. calling `ln.track()` without `LAMIN_INITIATED_BY_RUN_UID`, or skipping this entirely) breaks the lineage back to the agent run and defeats the point of tracking at all.

## Step 1 — Start of session (before the user's actual task)

The first action of a session is one question, before any `lamin` command, shell search, or version check.

Ask one blocking interactive question using this exact sentence: **"Track this session in LaminDB?"** Use exactly these two labels, in this order, without descriptions or recommendation text. Do not mark either answer as recommended or default, and never add "(Recommended)" to an option label or description.

1. **Track**
2. **Do not track**

**If your harness has a dedicated clarifying-question or ask-user tool, you must use it.** Only ask directly in response text if no such tool exists. Do not show a command or additional explanation. Stop and wait for the user's actual selection; do not assume one or continue in the same turn. Show this dialogue once at the start and never repeat it on a normal follow-up.

If the user selects **Do not track**, do the task with no LaminDB commands and no branch switch for the rest of the conversation.

### After Track

Run the version check described above. If the skill is stale, stop there.

Then resolve the session working directory. A development directory is one working tree for one instance, branch, and space. It is not the coding agent's own sandbox folder. Run these commands from the workspace directory.

1. `lamin settings dev-dir get`. A path means that directory already contains the session. Use it.
2. If it prints `None`, run `lamin settings dev-dir find` with the workspace as `PATH`. Do not pass `$HOME` and do not scan the machine.
3. One result: that directory is the session working directory. Run later LaminDB commands there.
4. Several: use the one that contains the files being edited. If none contains them, show the list and ask. Do not guess, and do not pick `$HOME`.
5. None: ask the user to run `lamin connect <account/name> --here` in the project directory. `lamin init` also creates a dev-dir in the working directory. Do not run either command unless they ask.
6. From that directory, choose a concise branch name in the form `<meaningful-task-slug>-<session-id-suffix>`. The slug must describe the user's actual task; never use a generic or timestamp-only name. Derive the suffix as specified in your harness reference (Cursor uses a unique agent-chosen suffix because it does not expose its session ID to shell commands); do not print it separately. Use only letters, digits, hyphens, or underscores, and never `/`. Then run:

```bash
lamin switch -c <branch-name>
```

`lamin switch -c` switches the branch of that directory. It does not create a child directory. Then `lamin track <agent>` as in the harness reference.

Escalate to the fallback below only if `lamin settings dev-dir get` errors (non-zero exit status). A command error here usually means `lamin` is only installed in a project-local virtualenv rather than on `PATH`:

```bash
LAMIN_BIN=$(find . -maxdepth 6 -type f -name lamin 2>/dev/null | head -1)
[ -z "$LAMIN_BIN" ] && LAMIN_BIN=$(command -v lamin 2>/dev/null)
if [ -z "$LAMIN_BIN" ]; then
  echo "NOT_FOUND: lamin"
else
  "$LAMIN_BIN" settings dev-dir get
fi
```

Determine which coding agent you are running as and follow the matching file under Quick reference below. Those references assume the session directory is already chosen; they do not run `dev-dir get` or `find` themselves.

The session working directory is immutable after it is resolved. `lamin switch` cannot change the parent agent process's working directory. The harness session folder may stay outside the dev-dir; that is fine. Therefore, **run every later LaminDB command and every task command from the session working directory**, using the execution tool's working-directory option when available or an explicit `cd "<session-working-dir>" &&` prefix otherwise. This includes `lamin track`, lineage verification, scripts and notebooks, direct-output attachment, tests, and `lamin finish`. Never set dev-dir to the harness session folder.

### Command hygiene

Keep every prescribed command free of diagnostic shell noise. Required working-directory setup (`cd` or the execution tool's working-directory option) and required environment configuration such as `LAMIN_SETTINGS_DIR` are allowed. Do not add any of the following:

- status headings or separators such as `echo "--- dev-dir ---"`;
- manual exit-code output such as `echo "exit: $?"` or `echo "switch exit: $?"` — rely on the execution tool's reported exit status;
- commands that inspect or print harness session IDs, including `echo`, `printenv`, `env`, or a Python command — except for the single Cursor marker command prescribed in its reference file, expand the environment variable only inside the branch name or required state-file path;
- convenience aliases such as `LAMIN=...` or `PYBIN=...` — invoke the required `lamin` or matching Python executable directly.

Run `lamin track` and `lamin finish` as standalone substantive commands: do not place another diagnostic or task command before or after either one in the same tool call.

For a follow-up prompt in the same harness conversation after Step 3 completed, do not ask again or create/switch another branch. Return to the same session working directory and run the same `lamin track <agent>` command before doing the follow-up work. The CLI uses the harness session ID to resume the existing agent Run and restores the active UID file needed for child-run lineage.

Each tracked mode starts with `lamin track <agent>`, which creates (or reuses) that harness's fixed Transform and opens a Run — see your reference file for the exact command and what it writes. **Run the exact command shown in your reference file from the session working directory as its own tool call — do not write your own tracking logic, add another command alongside it, or skip straight to the user's task.** If tracking isn't available (`lamin` not found, or the command errors — e.g. no lamindb instance connected), tell the user and proceed with their actual task untracked. Do not attempt Step 2/3 for the rest of the conversation because there is no Run to attach anything to. Do not undo the branch choice merely because tracking failed.

## Step 2 — During the session

Every script you write to do the task — the first one and every later one, on any message in this session — gets the `ln.track()`/`ln.finish()` instrumentation from "Self-tracking scripts" above, and gets run with `LAMIN_INITIATED_BY_RUN_UID` set per your harness's reference file — this already handles lineage for anything a script produces. Having done Step 1 once already does not exempt a later script from this wrapper. **This applies even to a script that just generates, fetches, or prepares input data before the "real" pipeline runs — there is no "just setup," "just a quick test," or "just dummy data" exemption.** If a script exists and gets executed, it gets the wrapper — full stop, regardless of how small or preliminary it feels. A local file that never went through `ln.track()`/`ln.Artifact(...).save()` cannot be turned into a real input later; the lineage gap it leaves is permanent.

For any script that's written, its inputs should not be local; they must be lamindb artifacts. If an input is available locally, it must be uploaded to lamindb as an artifact before running the script and adding the run to lamin. If any script was run or code was generated within the session to create that input, this script or code should also be added to lamin.

For EVERY python script that you plan to run, whether it was created/modified by you in the session or already existed, whether its run directly or indirectly, we want to verify that it insures proper lineage in lamindb before running it — if it does not, we will need to modify it to do so before running it. User confirmation is not required before doing this. This needs to be run without fail. If you can't run this, quit the session. **Hard gate: never run a script unless the most recent `verify_lineage(script_path)` result for that exact file prints `Lineage is fully tracked`. If it prints any `Missing lineage: ...`, do not run the script; modify it and repeat `verify_lineage` until it passes based on the 'Lineage remediation retry policy' section below. Never stop after the first failed fix attempt unless blocked by tool/runtime failure; you must complete up to 2 remediation iterations first.**

You can verify lineage by using the code below. The `script_path` variable should be set to the path of the script you are verifying. Run it with the same interpreter and LaminDB installation that will execute the script. A missing unrelated task dependency is not a reason to switch LaminDB environments.

```python
from lamindb.core import verify_lineage
result = verify_lineage(script_path)

if result.is_fully_tracked:
  print(f"Lineage is fully tracked")
else:
  print(f"Missing lineage: {result.missing_lineage}")

```

Mentally note only the files you create or modify **directly, with no script involved** — those need manual attachment before finishing (see Step 3, and your harness's reference file for how to resolve your run).

After a self-tracking script or notebook exits successfully, accept that result. Do not query or load its Transform, Run, or artifacts merely to reconfirm that they were saved, and do not search the filesystem for tracking state. Perform additional inspection only when the execution itself reports an error or when verification is part of the user's actual task rather than a tracking check.

Make sure you always do Step 3 at the end of the session, even if the user doesn't ask.

### Lineage remediation retry policy (mandatory)

When `verify_lineage(script_path)` fails for a script that is intended to run:

1. Attempt to fix the script and rerun `verify_lineage(script_path)`.
2. If it still fails, attempt one more fix and rerun `verify_lineage(script_path)` again.
3. Maximum remediation attempts: **2**.

After 2 failed remediation attempts:

- **Do not run the script.**
- Ask the user for guidance or a manual fix using the interactive ask-user tool (when available).
- Report both failed verify outputs and the exact remaining `missing_lineage` items.

Hard gate remains: only run when the **most recent** verify result for that exact file is `Lineage is fully tracked`.

### Lineage remediation guardrails

When fixing a script after `verify_lineage(script_path)` reports missing lineage, preserve script behavior and only add lineage tracking.

Non-negotiable rule:

- **Do not delete, comment out, or bypass file/folder path usage just to make verification pass. Only do it if overall script behavior can be preserved.**

Allowed direction:

- Add or adjust lineage instrumentation (`ln.track`, `ln.finish`, `ln.Artifact.get(...).load()`, `ln.Artifact(...).save()` or `ln.Artifact.from_*().save()`). Run `help(ln.Artifact)` for help with finding other methods for tracking artifacts in lamindb.

If lineage cannot be fixed without changing what the script does, stop and ask the user for guidance.

## Step 3 — End of session

User confirmation is not required. Always do Step 3. **Run the commands below exactly as shown — do not skip this step, and do not consider the task done until `lamin finish` has actually been run.**

If you created output files directly (no script involved), attach them first — see your harness's reference file for the exact command to resolve your run and attach files to it.

Then close the session from the immutable session working directory as its own tool call. Use the required environment prefix from the harness reference when one is prescribed; otherwise run:

```bash
lamin finish
```

Escalate to the fallback below only if this command errors (non-zero exit status) — under no other circumstance (not a lamindb warning, not wanting to double-check, not comparing against a local virtualenv version) should you run any additional command before or instead of accepting this result. A command error here usually means `lamin` is likely only installed in a project-local virtualenv rather than on `PATH`:

```bash
LAMIN_BIN=$(find . -maxdepth 6 -type f -name lamin 2>/dev/null | head -1)
[ -z "$LAMIN_BIN" ] && LAMIN_BIN=$(command -v lamin 2>/dev/null)
if [ -z "$LAMIN_BIN" ]; then
  echo "NOT_FOUND: lamin"
else
  "$LAMIN_BIN" finish
fi
```

This is the same command regardless of harness — it resolves whichever session is currently active on its own, renders the transcript as HTML, saves it as a report artifact, stamps all child scripts as session outputs (`Transform.run`), and closes the active tracking cycle. If the same harness conversation receives a later follow-up, `lamin track <agent>` resumes this Run and the next finish replaces its report and cumulative metrics with the complete conversation rather than creating a duplicate Run.

If Step 1 printed `NOT_FOUND`, there is no run to close — skip Step 3 entirely. If this command prints `NOT_FOUND`, or the binary itself errors (e.g. no lamindb instance connected): tell the user, skip the rest of tracking, and proceed with their actual task anyway — tracking infrastructure should never block the user's real request.

## Quick reference

- [Track Claude Code sessions](references/track_claude.md).
- [Track Copilot sessions](references/track_copilot.md). If this Copilot chat spawned a child, read that file first and skip dev-dir / branch / track.
- [Track Cursor IDE sessions](references/track_cursor.md).
- [Curate datasets](references/curate_datasets.md).
