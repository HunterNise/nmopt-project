# Project conventions

These are action-specific instructions for working in this repository. Read
this index first, then follow the mandatory routing rules below.

| Category | What it governs | Convention |
| --- | --- | --- |
| Workflow and execution | Work units, review boundaries, turns, cost, context, and handoffs | [Workflow](workflow.md) |
| Git and worktrees | Status, branches, diffs, commits, remotes, and repository safety | [Git](git.md) |
| Code and interfaces | C++, headers, tests, layer boundaries, and implementation design | [Code](code.md) |
| Build and test | CMake, configuration, compilation, test execution, and generated artifacts | [Build and test](build.md) |
| Documentation | Markdown, LaTeX, rendered documentation, and mathematical notation | [Documentation](documentation.md) |
| Explanation and learning | Technical explanations, code walkthroughs, and project-guided learning | [Explanation](explanation.md) and [commit walkthrough](commit-walkthrough.md) |

If an action crosses categories, read every applicable convention file.

## Mandatory routing

Follow these rules for repository actions:

1. Before any non-read-only work, inspect `git status`, preserve existing user
   work, and read both [Workflow](workflow.md) and [Git](git.md). These two
   files are mandatory even when the change will not create a commit.
2. Read every additional convention file applicable to the action.
    - Read [Git](git.md) for read-only inspection of repository status,
      history, branches, diffs, commits, remotes, or worktrees.
    - Read [Code](code.md) before inspecting, designing, reviewing, or changing
      C++, headers, tests, or code-layer interfaces.
    - Read [Build and test](build.md) before configuring, compiling, running
      tests, or changing CMake, build, or test configuration.
    - Read [Documentation](documentation.md) before reading, writing,
      rendering, or reviewing project documentation, Markdown, or LaTeX
      content.
    - Read [Explanation](explanation.md) before explaining or teaching code,
      technical concepts, tooling, or project-specific engineering decisions.
    - Read [Commit walkthrough](commit-walkthrough.md) in addition to
      [Explanation](explanation.md) when explaining Git history or code changes
      commit by commit.
3. Use [docs/README.md](../docs/README.md) to select the relevant design and
   implementation documents.

Use these files for working practices and repository boundaries. Treat the
documents under [`docs/`](../docs/README.md) as the authority for the
mathematical and semantic model.
