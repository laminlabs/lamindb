# Lamin Skills

These skills teach AI coding agents how to correctly use LaminDB, curate datasets, and interact with the Lamin ecosystem. They follow the open [Agent Skills standard](https://agentskills.io).

## Installation

**Note:** If you have installed `lamindb` (e.g., via `pip install lamindb`), you **already have these skills**! They are bundled directly with the package. After installing, run `uvx library-skills --all` so your agent can read them (add `--claude` for Claude Code).

Install them manually only if you want the skills from `main` that have not shipped in a `lamindb` release yet, or if you are working without the `lamindb` package installed.

### 1. Via `npx skills`

```bash
npx skills add laminlabs/lamindb
```

### 2. Via GitHub CLI

```bash
gh skill install laminlabs/lamindb
```

## Compatibility with `library-skills`

Lamin Skills are bundled directly with the `lamindb` Python package in the `.agents/skills/` directory. This makes them fully consistent with [`library-skills`](https://github.com/tiangolo/library-skills), allowing you to use the official Lamin skills seamlessly alongside other agent skills you might have in your projects.
