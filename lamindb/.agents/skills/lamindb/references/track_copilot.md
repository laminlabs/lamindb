# Track Copilot session in LaminDB

See [SKILL.md](../SKILL.md) for concepts and the shared steps — this covers only what's specific to Copilot.

## Parent vs child chats

If you created another chat with `create_session`, you are the **parent**:

- Do not run `lamin settings dev-dir get`, `lamin settings dev-dir find`, `lamin switch`, or `lamin track copilot`.
- Do not do the user's work in this chat.
- You may run `lamin finish` at the end. That updates the same report the child already wrote (parent chat, then child chat).

If this chat was created as a child, follow [SKILL.md](../SKILL.md) and the rest of this file as written: dev-dir, branch, track, artifacts, finish.

When shared Step 1 needs a branch, use only the task slug from SKILL.md, for example `lamin switch -c favorite-protein-fasta`. Other agent sessions may write to that same branch. If it already exists and its status is `merged`, ask as SKILL.md describes before reusing it. Never run `echo`, `printenv`, `env`, Python, or any other command to inspect or print the session ID, and never use a generic or timestamp-only branch name.

A Copilot session folder under `~/.copilot/repos/copilot-worktrees/` is not a dev-dir and is not the session working directory. Do not set dev-dir to that folder. The session working directory is the one resolved in [SKILL.md](../SKILL.md).

**Do not write your own tracking logic.** Run every command below exactly as shown, as its own tool call, in order — [SKILL.md](../SKILL.md)'s Step 1 first, including session-working-directory resolution, then Step 1 here, then the script/notebook command each time you run one, and Step 3 at the end. Don't skip a step because the task seems simple. If the user selected **Do not track**, stop here — there's nothing further to run, including Step 3. Otherwise, don't consider tracking finished until Step 3's `lamin finish` has actually run.

## Step 1 — Start of session

Run this now from the session working directory resolved in [SKILL.md](../SKILL.md)'s Step 1. **`--name` is mandatory — never omit it, never run this command without it:**

```bash
lamin track copilot --name "<one sentence describing this session's task>"
```

Escalate to the fallback below only if this command errors (non-zero exit status) — under no other circumstance (not a lamindb warning, not wanting to double-check, not comparing against a local virtualenv version) should you run any additional command before or instead of accepting this result. A command error here usually means `lamin` is likely only installed in a project-local virtualenv rather than on `PATH`:

```bash
LAMIN_BIN=$(find . -maxdepth 6 -type f -name lamin 2>/dev/null | head -1)
[ -z "$LAMIN_BIN" ] && LAMIN_BIN=$(command -v lamin 2>/dev/null)
if [ -z "$LAMIN_BIN" ]; then
  echo "NOT_FOUND: lamin"
else
  "$LAMIN_BIN" track copilot --name "<one sentence describing this session's task>"
fi
```

This reads Copilot's own `$COPILOT_AGENT_SESSION_ID` and writes `.copilot/.lamindb_run_uid_copilot_${COPILOT_AGENT_SESSION_ID}` under the session working directory. This keeps tracking state inside that directory. Repeating this command from the same session working directory for a follow-up in the same Copilot conversation resumes its existing Run; it does not create another Run.

You don't need to remember anything from this command's output — every later command below reads `$COPILOT_AGENT_SESSION_ID` from its own environment directly, the same way this one did.

## Running self-tracking scripts and notebooks

Run this exact pattern from the session working directory every time you execute a script or notebook — never without this wrapper, and never a hand-rolled `ln.track()` call without it either. **Do not add flags, error-suppression (`2>/dev/null`, `|| true`), or any other modification to the `cat` command — run it exactly as shown.** If the file doesn't exist, let `cat` fail visibly rather than silently substituting an empty value.

```bash
printf 'y\n' | LAMIN_INITIATED_BY_RUN_UID=$(cat ".copilot/.lamindb_run_uid_copilot_${COPILOT_AGENT_SESSION_ID}") <however you'd normally run this file>
```

The leading `printf 'y\n' |` auto-answers the "overwrite existing source code?" prompt `ln.track()` shows when a previously-tracked script's content has changed — normal when iterating — otherwise it hangs/crashes waiting for input that will never come.

Run the file with `LAMIN_INITIATED_BY_RUN_UID` set first. **Use the Python interpreter executable from the exact environment that provided the `lamin` executable used to start tracking: `/path/to/env/bin/lamin` requires `/path/to/env/bin/python`. Invoke that Python executable directly.**

If this command fails because a non-LaminDB task dependency is missing, add only that dependency to the same active project environment and retry with the same interpreter. If LaminDB itself cannot be imported by that interpreter, stop and ask the user to repair the environment. Do not create or switch to a temporary environment: it may resolve different LaminDB settings and save lineage to another branch.

## Step 3 — Attaching direct output files

If you created output files directly (no script involved), attach them from the session working directory. Build the state-file path from that same directory; do not use the base dev-dir when working in a branch-dir. **Invoke the Python interpreter executable directly from the exact environment that provided the `lamin` executable used to start tracking. If tracking used `/path/to/env/bin/lamin`, this command must use `/path/to/env/bin/python`; do not substitute another Python executable or place `lamin` before the Python arguments.**

```bash
<matching-python-executable> -c "
import lamindb as ln
from pathlib import Path
run = ln.Run.get(uid=Path('.copilot/.lamindb_run_uid_copilot_${COPILOT_AGENT_SESSION_ID}').read_text().strip())
ln.Artifact('output.csv', key='<meaningful/folder/path>/output.csv', description='<what it is>', run=run).save()
# repeat for each direct file
"
```

Replace `<matching-python-executable>` with the exact interpreter described above. If it cannot import LaminDB, stop and ask the user to repair that environment; do not attach through another environment.

A zero exit status means the direct files were attached successfully, even if the command prints no explicit artifact confirmation. Do not search the filesystem for tracking state or query LaminDB merely to reconfirm success. If the command fails, handle its reported error directly.

Then run [SKILL.md](../SKILL.md)'s Step 3 closing command (`lamin finish`) as its own tool call — it reads `$COPILOT_AGENT_SESSION_ID` from its own environment the same way Step 1 did. A later finish in the same Copilot conversation updates this Run's report and cumulative metrics. Don't stop after just writing/running the user's script.
