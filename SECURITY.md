# Security Policy

## Reporting a security issue

Please do not disclose a potential security vulnerability in a public issue when the report could expose credentials, private data, or a reproducible exploit.

Contact a project maintainer privately through the contact mechanism associated with the GitHub repository.

## Scope

Potential issues include unsafe handling of local files, command execution, temporary-file handling, credential exposure, or behavior that could compromise a user's system or data.

baGWAS executes external command-line programs (PLINK and GEMMA). Install these programs from trusted sources and verify executable paths before running analyses.
