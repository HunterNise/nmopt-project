# Git conventions

## Safety and branch policy

- Work on the currently selected branch or worktree; do not switch branches
  for convenience.
- `main` is reserved for the user’s stable version. On `main`, or on any
  branch that does not match `codex/*`, do not edit, create, delete, or rename
  repository files unless the user explicitly requests that change.
- Codex may make in-scope changes on `codex/*` branches, but must not switch
  branches on its own.
- Do not merge, push, force-push, reset shared history, or rewrite existing
  authorship unless explicitly asked.
- Do not run `git add` or `git commit` unless the user explicitly requests it.
  When commits are not requested, leave each coherent work unit pending for
  review, keep multiple units distinguishable in the diff, inspect the result,
  and suggest the prospective commit messages and split in the final reply.

## Reviewable commits

- Shape every work unit as a prospective commit boundary, whether or not the
  user has authorized actual commits. This planning requirement does not grant
  permission to stage or commit.
- Map each actual commit to one review-sized work unit from the
  [workflow conventions](workflow.md).
  Keep directly required tests and contract updates with the behavior they
  verify unless they form an independently meaningful documentation or test
  unit.
- Prefer separate commits for independently meaningful prerequisites,
  mechanical migrations, interface changes, integrations, and documentation.
  Do not split merely to reduce line count when doing so creates an unbuildable
  or misleading intermediate state.
- Do not create a separate commit for a tiny follow-up with no independent
  review or revert value. Include it with the unit it supports; when it is
  independently meaningful, it may be a small separate commit.
- Treat roughly 400 changed non-generated lines or changes spanning multiple
  architectural layers as a prompt to revisit the decomposition before
  staging. This is a review heuristic, not a quota; exceeding it requires an
  explicit explanation of the coupling or cohesive generated/mechanical work.
  A genuinely coupled change may remain one detailed commit.
- Before every requested commit, inspect the staged diff and its statistics,
  run focused verification, and confirm that unrelated work is unstaged.

## Commit messages

Use this Conventional Commits format:

```text
<type>(<scope>): <imperative summary>
```

- Use a lowercase type such as `feat`, `fix`, `refactor`, `docs`, `test`,
  `build`, `ci`, `chore`, `perf`, or `revert`.
- Use a short, lowercase scope when it clarifies the affected component;
  omit it when it adds no information.
- Keep the summary imperative, specific, and preferably within 72 characters.
  Do not end it with a period.
- Add a body after a blank line when the reason, tradeoff, or compatibility
  impact is not clear from the summary.
- Mark breaking changes with `!` after the type or scope and explain them in
  the body or a `BREAKING CHANGE:` footer.
- Keep each commit focused on one coherent purpose.

Examples:

```text
docs(agents): define Codex edit and commit policy
docs(markdown): fix GitHub math rendering
refactor(spec): separate semantic contracts from backends
```
