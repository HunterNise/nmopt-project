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
- Do not run `git commit`; leave changes pending for user review. After each
  coherent change, inspect the diff and suggest a commit message in the final
  reply.

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
