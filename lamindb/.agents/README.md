# LaminDB skill

The contained skill teaches AI coding agents how to correctly use LaminDB and track their agent sessions. They follow the open [Agent Skills standard](https://agentskills.io).

## Installation

**Note:** If you have installed `lamindb` (e.g., via `pip install lamindb`), you **already have the skill**! It is bundled with the package. After installing, run `uvx library-skills --skill lamindb` so your agent can read it.

Install them manually only if you want the skills from `main` that have not shipped in a `lamindb` release yet, or if you are working without the `lamindb` package installed.

### 1. Via `npx skills`

```bash
npx skills add laminlabs/lamindb
```

### 2. Via GitHub CLI

```bash
gh skill install laminlabs/lamindb
```
