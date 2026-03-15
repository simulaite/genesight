# Security Policy

GeneSight processes sensitive personal genetic data. All analysis runs locally on the
user's machine, and no data is transmitted externally. We treat any violation of this
principle as a security issue.

## Supported Versions

| Version | Supported |
|---|---|
| 0.1.x | Yes |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly. **Do not open a
public GitHub issue for security vulnerabilities.**

Please contact us directly:

- **Email:** [info@simulaite.ai](mailto:info@simulaite.ai)
- **GitHub Security Advisory:** Use the [private security advisory](https://github.com/simulaite/genesight/security/advisories/new) feature on the repository.

Include as much of the following as you can:

- A description of the vulnerability and its potential impact
- Steps to reproduce the issue
- Affected versions
- Any suggested fix, if you have one

## What Qualifies as a Security Issue

Because GeneSight handles genetic data, the bar for what constitutes a security issue is
broader than a typical CLI tool. The following are all treated as security vulnerabilities:

- **Data leakage:** Any code path that transmits DNA data, analysis results, or user identifiers off the local machine.
- **Unintended network calls:** Any network activity from `genesight-core` or unexpected outbound connections from the CLI or GUI.
- **DNA data exposure:** Writing genetic data to unintended locations (temp files, logs, crash reports) where it could be accessed by other software.
- **Database integrity:** Corruption or unauthorized modification of the local annotation databases.
- **Privacy violations:** Any behavior that undermines the guarantee that all processing stays local.

General software security issues (e.g., buffer overflows, dependency vulnerabilities,
path traversal) are also in scope.

## Response Timeline

- **Acknowledgment:** Within 48 hours of your report.
- **Initial assessment:** Within 5 business days.
- **Fix or mitigation:** Depending on severity, but we aim to resolve critical issues as quickly as possible with a patch release.

We will coordinate with you on disclosure timing. We ask that you give us reasonable time
to address the issue before any public disclosure.

## Scope Exclusions

The following are generally not treated as security issues:

- Bugs in synthetic test fixtures
- Issues requiring physical access to the user's machine (beyond normal OS-level access)
- Vulnerabilities in upstream dependencies with no demonstrated impact on GeneSight

If you are unsure whether something qualifies, please report it anyway. We would rather
review a non-issue than miss a real vulnerability.

## Contact

**STONKS GmbH**
Buber-Neumann-Weg 68
60439 Frankfurt am Main

Managing Director: Farschad Hoshiar

Tel.: +49 69 4080 9836
E-Mail: [info@simulaite.ai](mailto:info@simulaite.ai)
