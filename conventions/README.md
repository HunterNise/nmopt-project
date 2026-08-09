# Project conventions

These are action-specific instructions for working in this repository. Read
this index first, then follow the mandatory routing rules below.

| Category | What it governs | Convention |
| --- | --- | --- |
| Git and worktrees | Status, branches, diffs, commits, remotes, and repository safety | [Git](git.md) |
| Code and interfaces | C++, headers, tests, layer boundaries, and implementation design | [Code](code.md) |
| Build and test | CMake, configuration, compilation, test execution, and generated artifacts | [Build and test](build.md) |
| Documentation | Markdown, LaTeX, rendered documentation, and mathematical notation | [Documentation](documentation.md) |
| Explanation and learning | Technical explanations, code walkthroughs, and project-guided learning | [Explanation](explanation.md) and [commit walkthrough](commit-walkthrough.md) |

If an action crosses categories, read every applicable convention file.

## Mandatory routing

Follow these rules for repository actions:

1. Inspect `git status` before changing anything and preserve existing user work.
2. Read the appropriate convention file for the action you need to take.
    - Read [Git](git.md) before changing anything. Follow its branch, safety, diff-inspection, and commit-message rules, even when you are not creating a commit.
    - Read [Code](code.md) before inspecting, designing, reviewing, or changing C++, headers, tests, or code-layer interfaces.
    - Read [Build and test](build.md) before configuring, compiling, running tests, or changing CMake, build, or test configuration.
    - Read [Documentation](documentation.md) before reading, writing, rendering, or reviewing project documentation, Markdown, or LaTeX content.
    - Read [Explanation](explanation.md) before explaining or teaching code, technical concepts, tooling, or project-specific engineering decisions.
    - Read [Commit walkthrough](commit-walkthrough.md) in addition to [Explanation](explanation.md) when explaining Git history or code changes commit by commit.
3. Use [docs/README.md](../docs/README.md) to select the relevant design and implementation documents.


Use these files for working practices and repository boundaries. Treat the
documents under [`docs/`](../docs/README.md) as the authority for the
mathematical and semantic model.
