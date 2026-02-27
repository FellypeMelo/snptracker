# Project Workflow (XP + AI-Augmented Development Edition)

---

# 1. Core Philosophy

This project follows **Extreme Programming (XP)** enhanced with **AI-assisted development**, under strict engineering discipline.

AI is treated as:

> A fast, junior pair-programmer that must be supervised at all times.

The human developer acts as:

> Senior Engineer, Architect, Navigator, and Quality Gatekeeper.

---

## 1.1 Non-Negotiable Principles

### 1. The Plan is the Source of Truth

* All tasks live in `plan.md`
* No undocumented work
* No “quick hacks”
* No hidden scope expansion

---

### 2. Architecture First, Code Second

* The AI must never invent architecture.
* All architectural decisions must be:

  * Documented in `tech-stack.md`
  * Reviewed before implementation
  * Versioned with justification

---

### 3. Test-Driven Development (MANDATORY)

TDD is not optional.

Every feature must follow:

```
Red → Green → Refactor
```

* Tests define behavior
* AI is not allowed to write implementation before failing tests exist
* Every bug fix requires a regression test

No test = No feature.

---

### 4. AI Must Plan Before Coding

Before any implementation:

1. AI must propose:

   * File changes
   * Function signatures
   * Data structures
   * Edge cases
   * Risks
2. Human reviews and approves
3. Only then code is written

AI is never allowed to "just start coding."

---

### 5. Baby Steps Only

Never request:

> “Build full authentication system”

Instead:

* Create password validation function
* Test password validation
* Integrate hashing
* Test hashing
* Add login endpoint
* Test endpoint

Small, isolated, reversible changes.

---

### 6. Continuous Integration Discipline

All merges must pass:

* Unit tests
* Coverage threshold (>80%)
* Lint
* Static analysis
* Type checks
* Security checks

If CI fails → rollback.

---

### 7. High Code Coverage (>80%)

Coverage is:

* Enforced
* Verified before task completion
* Mandatory for new modules

Coverage gaps must be explained.

---

### 8. User Experience First

All decisions prioritize:

* Clarity
* Simplicity
* Responsiveness
* Accessibility
* Mobile usability (if applicable)

---

### 9. Non-Interactive & CI-Aware Execution

Always prefer:

```
CI=true npm test
pytest -q --disable-warnings
```

No watch modes.
No hanging processes.

---

# 2. AI Usage Protocol (Critical Section)

This section governs how AI is used.

---

## 2.1 AI Pair Programming Model

Human = Navigator
AI = Pilot

Human:

* Defines scope
* Defines constraints
* Reviews architecture
* Approves steps
* Validates output

AI:

* Generates drafts
* Implements small scoped changes
* Writes tests (when instructed)
* Refactors under supervision

AI never:

* Makes architectural decisions alone
* Refactors large areas unsupervised
* Modifies unrelated files
* Introduces new dependencies silently

---

## 2.2 AI Prompt Structure Standard

All implementation prompts must follow:

1. Context
2. Constraints
3. Files allowed to modify
4. Required tests
5. Acceptance criteria
6. Explicit instruction:

   > “Do not modify unrelated code.”

---

## 2.3 Anti-Spaghetti Safeguards

Before accepting AI output, verify:

* No duplicated logic
* No large monolithic functions
* No unexplained abstractions
* No new dependencies without approval
* No hidden global state
* No broken separation of concerns

If any are found → reject and refine prompt.

---

# 3. Standard Task Lifecycle

---

## Step 1 — Select Task

Choose next `[ ]` in `plan.md`.

Sequential execution required unless justified.

---

## Step 2 — Mark In Progress

Change:

```
[ ] → [~]
```

Before coding begins.

---

## Step 3 — Define Scope Explicitly

Before writing tests:

* What files will change?
* What behavior is expected?
* What are edge cases?
* What should NOT change?

Document this in task notes.

---

## Step 4 — Write Failing Tests (Red Phase)

* Create test file
* Cover:

  * Success case
  * Failure case
  * Edge case
* Run tests
* Confirm failure

No failure = invalid TDD cycle.

---

## Step 5 — Approve Minimal Implementation Plan

Before writing code:

AI must propose:

* Function signatures
* Flow description
* Affected files

Human approves.

---

## Step 6 — Implement Minimal Code (Green Phase)

* Write smallest code possible
* No premature optimization
* No speculative abstractions
* Run tests
* Confirm passing

---

## Step 7 — Refactor

Allowed only when:

* All tests pass
* Coverage remains ≥ 80%
* Behavior unchanged

Refactor goals:

* Remove duplication
* Improve naming
* Improve cohesion
* Improve clarity

---

## Step 8 — Coverage Verification

Example:

```bash
pytest --cov=app --cov-report=term-missing
```

If coverage < 80%:

* Add tests
* Do not reduce scope

---

## Step 9 — Tech Stack Deviation Protocol

If implementation requires:

* New dependency
* New framework
* Architectural change

STOP.

1. Update `tech-stack.md`
2. Add dated justification
3. Commit
4. Resume

---

## Step 10 — Commit Code

Format:

```
feat(auth): Add password validation with TDD
```

Commit must include only related changes.

---

## Step 11 — Attach Git Note

Must include:

* Task name
* Summary
* Files modified
* Architectural reasoning
* Testing summary

Attach via:

```bash
git notes add -m "<note>" <commit_hash>
```

---

## Step 12 — Update plan.md

```
[~] → [x] abc1234
```

Append first 7 characters of commit hash.

---

## Step 13 — Commit plan.md

```
conductor(plan): Mark task 'X' as complete
```

---

# 4. Phase Completion & Checkpoint Protocol

(Enhanced with stricter AI validation)

---

## 4.1 Phase Scope Detection

Identify changes since last checkpoint:

```bash
git diff --name-only <previous_checkpoint_sha> HEAD
```

---

## 4.2 Mandatory Test Validation

For every code file changed:

* Ensure test file exists
* Follow project naming convention
* Match testing style

If missing → create tests.

---

## 4.3 Automated Test Execution

Announce command first:

Example:

> Running full test suite.
> Command: `CI=true npm test`

If failing:

* Maximum 2 debugging attempts
* If still failing → stop and escalate

---

## 4.4 Manual Verification Plan

Must analyze:

* `product.md`
* `product-guidelines.md`
* `plan.md`

Provide step-by-step:

* Commands
* URLs
* Expected outputs

Then ask:

> Does this meet your expectations?

Wait for confirmation.

---

## 4.5 Checkpoint Commit

```
conductor(checkpoint): End of Phase X
```

Attach full verification report as Git note.

---

## 4.6 Record Checkpoint SHA

Append to phase header:

```
## Phase 2 – Authentication [checkpoint: abc1234]
```

Commit plan update.

---

# 5. Quality Gates (Hard Blockers)

Before marking task complete:

* [ ] All tests passing
* [ ] Coverage ≥ 80%
* [ ] Lint clean
* [ ] Static analysis clean
* [ ] Type checks pass
* [ ] Public APIs documented
* [ ] No new security risks
* [ ] No architectural violations
* [ ] Mobile validated (if applicable)
* [ ] No hidden AI-generated complexity

Failure in any = task incomplete.

---

# 6. Definition of Done (Strict XP Edition)

A task is complete only when:

1. Tests written first
2. Tests pass
3. Coverage ≥ 80%
4. Code refactored
5. No duplication introduced
6. CI passes
7. plan.md updated
8. Git note attached
9. Architectural integrity preserved

---

# 7. Emergency Procedures (AI-Aware)

---

## Critical Production Bug

1. Create hotfix branch
2. Write failing regression test
3. Implement minimal fix
4. Run full CI
5. Deploy
6. Postmortem in plan.md

AI must not refactor during hotfix.

---

## Security Breach

1. Rotate secrets
2. Audit commits
3. Patch vulnerability
4. Add security test
5. Document incident

---

# 8. Deployment Workflow

---

## Pre-Deployment Checklist

* [ ] All CI green
* [ ] Coverage ≥ 80%
* [ ] No TODOs
* [ ] No console logs
* [ ] Migrations reviewed
* [ ] Backups verified

---

## Deployment Steps

1. Merge to main
2. Tag release
3. Deploy
4. Run migrations
5. Smoke test
6. Monitor logs

---

# 9. Continuous Improvement

Weekly review:

* Did AI introduce unnecessary complexity?
* Were tasks too large?
* Were tests written first?
* Was architecture respected?

If not → refine workflow.

---

# Final Philosophy

AI can create a 10x speed increase.

But only if:

* You enforce XP rigor
* You write tests first
* You review architecture
* You move in baby steps
* You treat AI as junior

Without discipline, AI becomes a technical debt factory.

With XP discipline, AI becomes a force multiplier.
