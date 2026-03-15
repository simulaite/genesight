# Contributing to GeneSight

Thank you for your interest in contributing to GeneSight. This project aims to provide
privacy-first, local-only analysis of personal DNA raw data against public genome
databases. Every contribution helps make genetic self-knowledge more accessible and
trustworthy.

## Reporting Bugs

Please open a [GitHub Issue](https://github.com/simulaite/genesight/issues) and include:

- A clear, descriptive title
- Steps to reproduce the problem
- Expected behavior vs. actual behavior
- Your OS, Rust toolchain version (`rustc --version`), and GeneSight version
- Any relevant log output or error messages

**Important:** Never include real DNA data in bug reports. Use synthetic data or
describe the issue in general terms.

## Suggesting Features

Open a GitHub Issue with the `enhancement` label. Describe the use case, why it matters,
and how you envision it working. Discussion before implementation saves everyone time.

## Development Setup

```bash
git clone https://github.com/simulaite/genesight.git
cd genesight
cargo build
cargo test
cargo clippy
cargo fmt -- --check
```

All four commands must pass before submitting a pull request.

## Code Conventions

- **Error handling:** Use `thiserror` in library code (`genesight-core`), `anyhow` in binaries (`genesight-cli`, `genesight-server`).
- **No `unwrap()` in library code.** Use proper error propagation with `?` or explicit matching.
- **Doc-comments:** All public functions and types must have English doc-comments.
- **Confidence tiers:** Every analysis result must carry a `ConfidenceTier` (`Tier1Reliable`, `Tier2Probable`, or `Tier3Speculative`).
- **No real DNA data:** The repository must only contain synthetic test fixtures. Never commit real genetic data.
- **Data source attribution:** Reports must cite the source databases (ClinVar, GWAS Catalog, etc.).
- **Medical disclaimer:** All user-facing reports must include a medical disclaimer. GeneSight does not provide medical advice.

## Architecture Constraints

GeneSight is a Rust workspace with four crates:

| Crate | Role | Constraints |
|---|---|---|
| `genesight-core` | Library | **No filesystem IO, no network calls.** Accepts `&str`, `&[u8]`, and `rusqlite::Connection` as inputs. |
| `genesight-cli` | CLI binary | Handles file IO, opens databases, calls into core. |
| `genesight-gui` | Desktop GUI | Native egui application, calls into core directly. |
| `genesight-server` | Web API binary | Axum-based HTTP layer (planned). |

If your change touches `genesight-core`, ensure it does not introduce any filesystem
access or network calls. All IO belongs in the binary crates.

## Pull Request Process

1. Fork the repository and create a feature branch from `main`.
2. Make your changes in focused, reviewable commits.
3. Ensure all checks pass: `cargo test`, `cargo clippy`, `cargo fmt -- --check`.
4. Write or update tests as appropriate.
5. Submit a pull request with a clear description of what changed and why.

A maintainer will review your PR. Small, focused PRs are reviewed faster than large ones.

## Testing

- **Unit tests** live alongside the code they test (in `#[cfg(test)]` modules).
- **Integration tests** use synthetic fixtures in `tests/fixtures/`.
- For pharmacogenomics and annotation testing, see `docs/PGP_TEST_DATA.md` for guidance on using public PGP (Personal Genome Project) data.

Run the full test suite before submitting:

```bash
cargo test
```

## Commit Conventions

This project follows [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` -- a new feature
- `fix:` -- a bug fix
- `docs:` -- documentation changes
- `refactor:` -- code restructuring without behavior change
- `test:` -- adding or updating tests
- `chore:` -- build, CI, or tooling changes

Example: `feat(pgx): add CYP2D6 metabolizer status lookup`

## Privacy and Safety

GeneSight processes sensitive genetic data. Keep these principles in mind:

- **Privacy first.** All processing happens locally. No data leaves the user's machine.
- **No real DNA data in the repository.** Use only synthetic or publicly consented (PGP) data for testing.
- **No medical claims.** GeneSight is an informational tool, not a diagnostic device. Never frame outputs as medical advice.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/).
We are committed to providing a welcoming and inclusive experience for everyone.

## License

By contributing to GeneSight, you agree that your contributions will be licensed under
the [GNU General Public License v3.0 or later](LICENSE) (GPL-3.0-or-later).
