# Contributing to Khaos -- Semiconductor Analysis & Simulation Platform

Thank you for taking the time to look at Khaos. This project is an
education-and-research-grade real-time TCAD sandbox: every contribution
that improves accuracy, performance, robustness or pedagogical clarity
is welcome.

The goal of this document is to lower the friction of "I want to fix /
add something" so you can land your patch quickly without re-discovering
unwritten rules.

---

## Table of contents

1. [Code of conduct](#code-of-conduct)
2. [Where to start](#where-to-start)
3. [Build & run](#build--run)
4. [Reporting issues](#reporting-issues)
5. [Pull requests](#pull-requests)
6. [Coding style](#coding-style)
7. [Testing requirements](#testing-requirements)
8. [Physics / numerical changes](#physics--numerical-changes)
9. [Commit messages](#commit-messages)
10. [Licensing](#licensing)

---

## Code of conduct

Be kind, be precise. Disagreement is welcome; ad-hominem and personal
attacks are not. Maintainers reserve the right to lock or close threads
that drift away from the technical question.

If you are reporting a sensitive issue (security, data integrity), email
the maintainer directly rather than opening a public issue.

---

## Where to start

* Read the [README](README.md) and skim
  [`src/PhysicsEngine.hpp`](src/PhysicsEngine.hpp) and
  [`src/DriftDiffusion.hpp`](src/DriftDiffusion.hpp).
  Both headers carry the formulas and references that the code
  implements -- they are the source of truth, not the cpp files.
* Look at the open issues, especially anything labelled
  `good first issue` or `help wanted`.
* The [Roadmap](README.md#roadmap) lists ideas with agreed-upon designs
  but no implementation yet -- a good place to find a self-contained
  starter project.

---

## Build & run

The full instructions live in the [README](README.md#build). The short
version:

```bash
# Linux
sudo apt install build-essential cmake ninja-build libsfml-dev
cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

# Windows / MSYS2 UCRT64
pacman -S --needed mingw-w64-ucrt-x86_64-{gcc,cmake,ninja,sfml}
cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
cmake --build build -j
```

If your patch touches anything subtle (memory ownership, raw pointers,
multi-thread state), please also run the sanitizer build at least
once before pushing:

```bash
cmake -S . -B build-asan -G Ninja \
      -DCMAKE_BUILD_TYPE=Debug \
      -DBUILD_TESTS=ON \
      -DUSE_SANITIZER=ON
cmake --build build-asan -j
ASAN_OPTIONS=detect_leaks=1:abort_on_error=1 \
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 \
ctest --test-dir build-asan --output-on-failure
```

---

## Reporting issues

Open an issue **before** starting non-trivial work; it is much cheaper
to align on the design first than to throw away a 600-line PR.

### Bug reports

A high-quality bug report has:

1. **Title** -- one line, what is broken.
2. **Environment** -- OS + compiler + SFML / CMake versions, branch and
   commit hash (output of `git rev-parse --short HEAD`).
3. **Reproduction** -- the smallest sequence of UI clicks (or, even
   better, a CSV / preset) that shows the bug.
4. **Expected vs actual** -- what the physics / UI should do versus what
   it does.
5. **Logs / screenshots** -- console output, ImGui screenshots, or short
   GIFs are gold.

The repo ships an [issue template](.github/ISSUE_TEMPLATE) where
applicable; please fill in every section.

### Feature requests

Tell us:

1. The *physics or UX problem* the feature solves.
2. A reference (Sze / Pierret / Selberherr / IEEE paper) if the feature
   is a textbook model.
3. A sketch of the UI surface and the engine API change.

---

## Pull requests

We follow a fast-forward-friendly workflow:

1. **Fork** the repo and create a branch off `main`:
   `feat/heterojunction-poisson`, `fix/sg-bernoulli-overflow`,
   `docs/roadmap-update`, ...

   Branch name prefixes:

   | Prefix      | Use for                                               |
   | ----------- | ----------------------------------------------------- |
   | `feat/`     | new physics / UI features                             |
   | `fix/`      | bug fixes (must reference an issue if one exists)     |
   | `perf/`     | performance work without behavioural change          |
   | `refactor/` | internal cleanup with no observable change           |
   | `docs/`     | README / CONTRIBUTING / inline-comment changes        |
   | `ci/`       | GitHub Actions, CMake, sanitizers, build tooling      |
   | `test/`     | new or rewritten GoogleTest cases                     |

2. Keep the PR **focused**. Do not mix a refactor with a bugfix.

3. **Rebase, do not merge.** If `main` has moved, rebase your branch:
   `git fetch origin && git rebase origin/main`.

4. Open the PR, fill in the [pull request template](.github/PULL_REQUEST_TEMPLATE.md):

   * **What** -- one paragraph.
   * **Why** -- the user / physics motivation.
   * **How** -- the design choice (link to a paper if relevant).
   * **Tested how** -- which CTest targets, manual UI tests, sanitizer
     run.
   * **Risks** -- known limitations, follow-ups, behavioural changes.

5. CI must be green before review. The GitHub Actions matrix builds on
   Linux GCC, Linux GCC with ASan + UBSan, and Windows MSYS2 UCRT64 --
   if any one of those is red the PR will not be merged.

6. Address review comments by **adding new commits** (do not force-push
   during review except to fix CI). Once the PR is approved, the
   maintainer may squash-merge.

---

## Coding style

* **Formatter** -- `.clang-format` is the source of truth. Run
  `clang-format -i src/*.cpp src/*.hpp tests/*.cpp` before pushing.
  CI runs `clang-format --dry-run` and posts diff hints, but does not
  block merge -- maintaining the lint is the contributor's
  responsibility.

* **Standard** -- C++20. Prefer `std::span`, `std::optional`,
  `[[nodiscard]]`, `constexpr`, `concepts` and `if constexpr` where
  they replace older idioms cleanly.

* **Naming**

  | Kind                  | Convention            | Example                  |
  | --------------------- | --------------------- | ------------------------ |
  | Types / classes       | `UpperCamelCase`      | `DriftDiffusion`         |
  | Functions / methods   | `lowerCamelCase`      | `solveGummel`            |
  | Member variables      | `m_lowerCamelCase`    | `m_psi`, `m_T_ambient`   |
  | File-static / locals  | `lowerCamelCase`      | `cellPx`                 |
  | Enums                 | `UpperCamelCase` + `UpperCamelCase` values | `enum class CellRegion { Bulk, Emitter, Base }` |
  | `constexpr` constants | `kUpperCamelCase`     | `kWindowWidth`           |
  | Macros                | avoid; if needed `KHAOS_UPPER_SNAKE`                |

* **No exceptions in the hot path.** All per-frame solvers are marked
  `noexcept` and operate on pre-allocated buffers. New code in those
  paths must respect that contract.

* **Zero allocation** is a core invariant of every solver entry point.
  Use member buffers preallocated in the constructor; do **not** add
  `std::vector::resize` or `new` inside `solve*`, `step*` or `render*`
  functions on the hot path.

* **Header hygiene**

  * Public headers (`*.hpp`) keep `#include` minimal -- forward-declare
    where possible.
  * `<bits/stdc++.h>` is forbidden.
  * Each header includes its own dependencies (no transitive reliance
    on a `.cpp` file's earlier includes).

* **Comments** carry references to the formula and the source. Look at
  any block in `PhysicsEngine.hpp` or `DriftDiffusion.cpp` for the
  expected density of citations (Sze section / Pierret equation /
  Selberherr chapter / IEEE paper).

* **Do not commit** generated files: `build/`, `.vscode/settings.json`,
  IDE caches. The `.gitignore` already covers the common ones.

---

## Testing requirements

* Add or update **GoogleTest** cases in `tests/test_physics.cpp` for any
  new physics, mathematical helper, or numerical fix.

* Tests must:
  * Run in well under one second each.
  * Compare against an analytical limit, a tabulated value (Sze /
    Pierret), or a previous version's output -- not against a
    floating-point regurgitation of the implementation.
  * Be deterministic. No `std::rand`, no time-of-day dependence.

* The full suite must pass on every supported toolchain. CI enforces
  this; see [`.github/workflows/build.yml`](.github/workflows/build.yml).

* **Do not delete or weaken existing tests** to make a PR pass.
  If a test is genuinely wrong, document why in the PR body and have a
  maintainer review the change explicitly.

---

## Physics / numerical changes

Patches that touch the engine deserve extra care:

1. **Cite the source.** Every formula gets a reference: textbook with
   section / equation, or an IEEE / journal paper with year. The
   existing code is the template -- match its citation density.

2. **Document the limit.** State explicitly when the model breaks down
   (high field, low T, degenerate regime, ...). The simulator is used
   for teaching; misleading the user about the model's domain of
   validity is a worse bug than getting a number 5% wrong.

3. **Verify in equilibrium.** Most rates and currents must vanish under
   detailed balance (`np = n_i^2`). Add a test that pins this.

4. **Stability before speed.** A solver that diverges at 5 V is a
   non-starter, even if it is 2x faster. Show the convergence /
   stability margin in the PR body.

5. **Numerical safety.** Branch-clamp `exp` arguments at +/-40, gate
   `log` and `sqrt` on positive inputs, prefer `std::expm1` /
   `std::log1p` near zero. Look at `PhysicsEngine::bernoulli` for the
   in-house template.

---

## Commit messages

We use [Conventional Commits](https://www.conventionalcommits.org/) with
the scopes that match our directory layout.

```
<type>(<scope>): <imperative summary, <=72 chars>>

<optional body explaining *why*; wrap at 80 cols>

<optional footer: refs / breaking changes>
Refs: #123
```

Common types: `feat`, `fix`, `perf`, `refactor`, `docs`, `test`,
`build`, `ci`, `chore`.

Common scopes: `physics`, `dd` (drift-diffusion), `bandview`, `crystal`,
`ui`, `cmake`, `ci`.

Examples:

```
feat(dd): add Anderson chi-shift to Scharfetter-Gummel argument

The electron drift potential eta_n = psi + (chi - chi_ref) / q now
captures the conduction-band step at heterojunctions, so a Si/Ge
brush boundary blocks back-injected holes correctly (textbook
HBT result).

Refs: #42
```

```
fix(physics): clamp Bernoulli denominator to avoid IEEE inf at x ~ 700
```

```
docs(roadmap): mark Phase 5 (Wachutka + heterojunctions) as shipped
```

---

## Licensing

Khaos is MIT-licensed. By submitting a pull request you agree that your
contribution is licensed under the same MIT terms. If your employer
imposes additional licensing constraints, mention this in the PR body
and we will work with you to find a compatible path before merging.

We do not currently require a CLA / DCO sign-off, but adding
`Signed-off-by: Your Name <email>` to your commit messages
(`git commit -s`) is appreciated and forward-compatible should we
introduce one later.

---

## Thank you

Maintaining a real-time semiconductor sandbox is a niche pleasure.
Every issue, every benchmark, every typo fix makes the platform a
sharper teaching tool. Welcome aboard.
