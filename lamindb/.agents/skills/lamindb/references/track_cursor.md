# Track Cursor IDE session in LaminDB

See [SKILL.md](../SKILL.md) for the shared steps. This reference covers only Cursor IDE Agent chats; it does not apply to Cursor's standalone CLI or cloud agents.

Cursor does not expose its conversation ID to commands run by the IDE Agent. Generate one random marker per Cursor conversation and pass it to the CLI, which uses the marker to identify the conversation in Cursor's local chat database. Parallel conversations must use different markers.

When shared Step 1 chooses a new branch, use a task-specific name with a unique agent-chosen suffix of at least eight hexadecimal characters, for example `favorite-protein-fasta-a1b2c3d4`. Keep this branch dedicated to the current Cursor conversation.

If the user chose **Do not track**, stop here. Otherwise complete [SKILL.md](../SKILL.md)'s Step 1, including session-working-directory resolution, before running the commands below. Stay in that directory. Do not write your own tracking logic.

## Step 1 — Start of session

Invent a 32-character lowercase hexadecimal marker yourself. Do not run Python, openssl, or any other program to create it, and do not import a library. Generate it once the session working directory is resolved. Print it once:

```bash
echo LAMIN_CURSOR_SESSION_ID=<32 lowercase hex characters>
```

Do not print or invent another marker on follow-ups; reuse the same value. Remember the complete value after `=` for this conversation. Do not ask the user to copy it. Then run the tracking command as its own tool call. `--name` is mandatory:

```bash
LAMIN_CURSOR_SESSION_ID=<generated value> lamin track cursor --name "<one sentence describing this session's task>"
```

Substitute the exact generated value without angle brackets. Only if this command errors, use the same `LAMIN_BIN` fallback described in [SKILL.md](../SKILL.md), preserving the `LAMIN_CURSOR_SESSION_ID` prefix and substituting `track cursor --name "<one sentence describing this session's task>"` for the command arguments.

Remember the LaminDB Run UID printed by this command for the execution wrappers below. Do not ask the user to copy it or run another command to retrieve it. A follow-up in the same Cursor conversation should repeat the prefixed tracking command with the same marker from the same session working directory, which resumes the existing Run.

## Running self-tracking scripts and notebooks

Use the matching Python executable from the environment that provided `lamin`, as required by [SKILL.md](../SKILL.md). Run each script or notebook with the agent Run UID linked:

```bash
printf 'y\n' | LAMIN_INITIATED_BY_RUN_UID=<Run UID printed by lamin track cursor> <however you'd normally run this file>
```

Substitute the exact Run UID printed earlier; do not include angle brackets. Do not run the file without the wrapper. Follow the shared lineage-verification requirement before execution.

## Step 3 — Attaching direct output files

For files created directly, with no script involved, attach them from the session working directory using the matching Python executable:

```bash
<matching-python-executable> -c "
import lamindb as ln
run = ln.Run.get(uid='<Run UID printed by lamin track cursor>')
ln.Artifact('output.csv', key='<meaningful/folder/path>/output.csv', description='<what it is>', run=run).save()
"
```

Substitute the exact Run UID printed earlier; do not include angle brackets. Then close the session as its own tool call using the same marker:

```bash
LAMIN_CURSOR_SESSION_ID=<generated value> lamin finish
```

Substitute the exact generated value without angle brackets. The CLI finds the corresponding conversation in Cursor's local chat database and renders its report. A later finish after a follow-up updates the same Run's report.
