# Git Commit Walkthrough

Use this workflow when the user wants to learn or understand code changes commit by commit.

Apply the general Explanation and Learning Conventions in addition to this workflow.

## 1. Establish the comparison

Determine:

* current branch;
* upstream or pushed baseline;
* relevant commit range;
* chronological commit order.

If the user specifies a baseline, use it.

If the user says they are familiar with previously pushed commits, treat the pushed state as known and focus on commits after that state.

Do not begin by reading every commit diff.

---

## 2. Reconnaissance

First gather lightweight metadata for the candidate commits.

Prefer information such as:

```text
hash
subject
commit type
changed files
diffstat
```

Use this information to determine which commits are relevant.

If the user asks to include only particular categories, such as `feat`, filter accordingly.

If the user asks to ignore documentation or chores, exclude them unless they contain code changes that materially affect understanding and therefore require explicit mention.

Do not deeply inspect commit contents during this stage unless metadata is insufficient to classify a commit.

---

## 3. Initial roadmap

Present a compact roadmap of the relevant commits.

For example:

```text
There are four relevant code commits:

1. abc123 — introduces mesh loading
2. def456 — adds boundary classification
3. 789abc — moves assembly into its own component
4. fed321 — connects the new component to the solver
```

This roadmap is provisional.

Do not pretend the commit subject alone establishes all implementation details.

If deeper inspection later changes the interpretation, update the roadmap.

---

## 4. Study one commit at a time

By default, deeply inspect only the first relevant unexplained commit.

Do not inspect later full diffs simply because they are available.

The default turn boundary is:

> one commit per turn.

A commit MAY require multiple turns if it contains several substantial conceptual units or the user asks detailed questions.

Do not proceed to the next commit automatically after a substantial explanation unless the user explicitly asks to continue or the next commit is required to complete the current explanation.

---

## 5. Context for the current commit

For the current commit, inspect:

* the full diff;
* enough surrounding code to understand changed sections;
* relevant declarations or definitions;
* relevant documentation if it records the reason for the change.

Retrieve additional files only when necessary to explain the current commit.

Do not recursively inspect unrelated implementation details.

Use the rule:

> Understand dependencies deeply enough to explain their role in this commit, but do not turn every dependency into a separate investigation.

---

## 6. Preserve the historical state

Explain the commit relative to the state immediately before it.

Think in terms of:

```text
previous known state
→ this commit
→ resulting state
```

Do not explain the current implementation primarily using concepts introduced in later commits.

If later work changes the design, leave that for the later commit unless a short forward reference prevents confusion.

Example:

> At this point the loader still owns the object directly. A later commit may change that ownership, but for this commit we should understand the direct-ownership design first.

---

## 7. Suggested explanation structure

Adapt the structure to the commit. Do not force empty sections.

A useful default is:

### Purpose

What capability or structural change does this commit introduce?

Distinguish documented intent from inferred intent.

### Files changed

Give a short map of which files matter and what role each plays.

Do not yet explain every line.

### Walkthrough

Inspect the important diff sections in a useful order.

This MAY follow:

* file order;
* execution order;
* dependency order;
* conceptual order.

Choose whichever makes the change easiest to understand.

### Syntax and APIs

Explain new or significant:

* C/C++ constructs;
* build-system syntax;
* deal.II APIs;
* standard-library APIs;
* language features;
* framework conventions.

Do not re-teach basic constructs already known.

### Behavior

Explain:

* control flow;
* data flow;
* state changes;
* lifetime and ownership;
* object relationships;
* invariants;
* effects on existing behavior.

### Design reasoning

Explain why the new structure makes sense when evidence supports it.

Discuss alternatives when doing so clarifies the chosen design.

Do not invent undocumented intent.

### Resulting mental model

End with a compact description of what is now true after this commit.

Do not automatically begin the next commit.

---

## 8. Previously learned concepts

Assume concepts from earlier explained commits remain known.

For example, if an earlier commit established:

* what `std::unique_ptr` means;
* how `Triangulation` is used;
* the purpose of a particular project class;

then later commits should refer to those concepts without repeating the full explanation.

Revisit an established concept only when:

* its role changes;
* a subtle new aspect appears;
* the user asks;
* comparison with a new construct is important.

Example:

> This is still the ownership pattern introduced in the previous commit, but the important difference here is that the second component now retains access beyond construction.

---

## 9. Generated and moved code

Treat generated, moved, renamed, and mechanically reformatted code carefully.

Do not spend substantial explanation time on a section merely because the diff is large.

Distinguish:

* semantic changes;
* structural moves;
* renames;
* formatting;
* generated boilerplate.

When code is moved unchanged, explain the architectural significance of the move rather than re-explaining its implementation unless necessary.

---

## 10. Commit filters

Honor user-provided filters strictly.

Examples:

```text
only feat commits
ignore docs
ignore chore
only commits affecting code
ignore formatting-only changes
```

If a supposedly excluded commit is required to understand the code state, mention this explicitly rather than silently treating it as a main learning unit.

For example:

> Commit X is labeled `chore`, so I am not treating it as a lesson unit. It does, however, rename the target used by the following feature commit, so I will account for that rename as context.

---

## 11. Updating the roadmap

The initial roadmap is based on lightweight inspection.

After studying a commit, new structure may become visible.

Update the remaining roadmap when useful.

Do not preserve an inaccurate initial classification merely for consistency.

Example:

```text
Initial view:
2. add mesh configuration

After inspecting commit 1:
2. separate mesh ownership from simulation configuration
```

Explain significant roadmap corrections briefly.

---

## 12. Completion state

After each commit, conceptually retain:

```text
commit understood
new behavior
new architecture
new terminology
new syntax/API concepts
changed assumptions
remaining open questions
```

Use this compact state as the baseline for the next commit.

Do not repeatedly reconstruct all earlier commits from scratch.

---

## 13. Default execution pattern

```text
identify baseline
→ list candidate commits
→ filter using metadata
→ show concise roadmap
→ inspect commit 1 deeply
→ explain commit 1
→ discuss questions
→ inspect commit 2 only when moving to it
→ ...
```

The key constraint is:

> Discover broadly, study narrowly.
