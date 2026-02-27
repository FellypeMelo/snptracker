# Foundational Mandates (XP + AI-Augmented)

These instructions take absolute precedence over all other defaults. Adhere to them strictly.

## 1. Core Philosophy
- **AI as Pilot, Human as Navigator:** AI implements small, scoped changes under continuous human supervision. AI never makes architectural decisions alone.
- **Extreme Programming (XP) Foundation:** Work in baby steps, maintain high discipline, and prioritize simplicity.

## 2. Implementation Protocol
- **Test-Driven Development (MANDATORY):** Every feature and bug fix MUST follow the `Red -> Green -> Refactor` cycle. No implementation code is written without a preceding failing test.
- **Mandatory Planning:** Before ANY implementation, the AI must propose:
  - File changes
  - Function signatures
  - Data structures
  - Edge cases and risks
- **Baby Steps:** Break complex tasks into small, isolated, and reversible changes. Do not attempt monolithic implementations.

## 3. Architecture & Tech Stack
- **Architecture First:** Document all architectural decisions in `tech-stack.md` before implementation.
- **Plan as Source of Truth:** All work must be tracked in `plan.md`. No undocumented work or "quick hacks."

## 4. Quality & Verification
- **High Coverage (>80%):** Enforce and verify coverage before task completion.
- **Non-Interactive Execution:** Use non-interactive, CI-aware commands (e.g., `pytest -q`).
- **Standard Task Lifecycle:** Follow the 13-step lifecycle defined in `conductor/workflow.md` precisely.

## 5. Security & Safety
- **Credential Protection:** Never log, print, or commit secrets.
- **Fail Securely:** Ensure error handling does not expose sensitive information.

---
*Refer to `conductor/workflow.md` for the full operational protocol.*
