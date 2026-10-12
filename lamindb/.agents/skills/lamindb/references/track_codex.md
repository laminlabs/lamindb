# Track Codex session in LaminDB

See [SKILL.md](../SKILL.md) for concepts and the shared steps — this covers only what's specific to Codex.

When shared Step 1 chooses a branch, use only the task slug from SKILL.md, for example `lamin switch -c curate-cell-types`. Other agent sessions may write to that same branch. If it already exists and its status is `merged`, ask as SKILL.md describes before reusing it. Never run `echo`, `printenv`, `env`, Python, or any other command to inspect or print the session ID, and never use a generic or timestamp-only branch name.

## Step 1 — Start of session

Run this now from the session working directory resolved in [SKILL.md](../SKILL.md)'s Step 1. **`--name` is mandatory — never omit it, never run this command without it:**

```bash
lamin track codex --name "<one sentence describing this session's task>"
```

Codex keeps `.codex` read-only while `sandbox_mode` is `workspace-write`. Read this turn's permissions instructions. When they allow `sandbox_permissions=require_escalated`, run this command, the fallback below, and the later `lamin finish` with `sandbox_permissions=require_escalated` on the first call. Do not run any of them in the default sandbox first. When `sandbox_mode` is full access, or those instructions say escalated permissions must not be requested, run the same commands with no `sandbox_permissions` argument. If a default-sandbox run still prints `Operation not permitted` for `.codex` and escalation is allowed, rerun that same command once with `sandbox_permissions=require_escalated`.

Escalate to the fallback below only if this command errors (non-zero exit status). The `.codex` permission warning above is the one warning that must be rerun, and only when escalation is allowed. Under no other circumstance (not another lamindb warning, not wanting to double-check, not comparing against a local virtualenv version) should you run any additional command before or instead of accepting this result. A command error here usually means `lamin` is likely only installed in a project-local virtualenv rather than on `PATH`:

```bash
LAMIN_BIN=$(find . -maxdepth 6 -type f -name lamin 2>/dev/null | head -1)
[ -z "$LAMIN_BIN" ] && LAMIN_BIN=$(command -v lamin 2>/dev/null)
if [ -z "$LAMIN_BIN" ]; then
  echo "NOT_FOUND: lamin"
else
  "$LAMIN_BIN" track codex --name "<one sentence describing this session's task>"
fi
```

This writes `.codex/.lamindb_run_uid_codex_${CODEX_THREAD_ID}` under the session working directory, keyed by Codex's own thread id. Codex sets `$CODEX_THREAD_ID` in every tool call. This keeps tracking state inside that directory. Repeating this command from the same session working directory for a follow-up in the same Codex conversation resumes its existing Run; it does not create another Run.

## Running self-tracking scripts and notebooks

`$CODEX_THREAD_ID` is already set in every subprocess Codex spawns, so finding your own run is a plain `cat`. Run the script or notebook from the session working directory with `LAMIN_INITIATED_BY_RUN_UID` set first. **Use the Python interpreter executable from the exact environment that provided the `lamin` executable used to start tracking: `/path/to/env/bin/lamin` requires `/path/to/env/bin/python`. Invoke that Python executable directly.** Do not add flags, error-suppression (`2>/dev/null`, `|| true`), or any other modification to the `cat` command — run it exactly as shown. If the file doesn't exist, let `cat` fail visibly rather than silently substituting an empty value:

```bash
printf 'y\n' | LAMIN_INITIATED_BY_RUN_UID=$(cat .codex/.lamindb_run_uid_codex_${CODEX_THREAD_ID}) <however you'd normally run this file>
```

The leading `printf 'y\n' |` auto-answers the "overwrite existing source code?" prompt `ln.track()` shows when a previously-tracked script's content has changed — normal when iterating — otherwise it hangs/crashes waiting for input that will never come.

If this command fails because a non-LaminDB task dependency is missing, add only that dependency to the same active project environment and retry with the same interpreter. If LaminDB itself cannot be imported by that interpreter, stop and ask the user to repair the environment. Do not create or switch to a temporary environment: it may resolve different LaminDB settings and save lineage to another branch.

## Step 3 — Attaching direct output files

If you created output files directly (no script involved), run this from the session working directory. Build the state-file path from that same directory. **Invoke the Python interpreter executable directly from the exact environment that provided the `lamin` executable used to start tracking. If tracking used `/path/to/env/bin/lamin`, this command must use `/path/to/env/bin/python`; do not substitute another Python executable or place `lamin` before the Python arguments.**

```bash
<matching-python-executable> -c "
import lamindb as ln
from pathlib import Path
run = ln.Run.get(uid=Path('.codex/.lamindb_run_uid_codex_${CODEX_THREAD_ID}').read_text().strip())
ln.Artifact('output.csv', key='<meaningful/folder/path>/output.csv', description='<what it is>', run=run).save()
# repeat for each direct file
"
```

Replace `<matching-python-executable>` with the exact interpreter described above. If it cannot import LaminDB, stop and ask the user to repair that environment; do not attach through another environment.

A zero exit status means the direct files were attached successfully, even if the command prints no explicit artifact confirmation. Do not search the filesystem for tracking state or query LaminDB merely to reconfirm success. If the command fails, handle its reported error directly.

Then close the session per [SKILL.md](../SKILL.md) Step 3 (`lamin finish`), with `sandbox_permissions=require_escalated` as stated in Step 1. A later finish in the same Codex conversation updates this Run's report and cumulative metrics.
