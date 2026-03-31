# Security Policy

## Supported Versions

We currently provide security updates for the following versions of LISFLOOD-FP:

| Version | Supported          |
| ------- | ------------------ |
| 8.2.x   | :white_check_mark: |
| 8.1.x   | :white_check_mark: |
| 8.0.x   | :x:                |
| < 8.0   | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability in LISFLOOD-FP, please report it by opening a new issue with the label "security".

For severe vulnerabilities that should not be disclosed publicly, please contact the maintainers directly.

We will acknowledge receipt of your vulnerability report as soon as possible and will send you regular updates about our progress. If you have a solution, we would welcome a pull request.

## Security Best Practices

When using LISFLOOD-FP:

1. Always verify input files before processing, especially if they come from untrusted sources
2. Be cautious when running simulations with extremely large datasets, as they may cause resource exhaustion
3. Keep your installation updated with the latest security patches
4. When using LISFLOOD-FP in production environments, consider running it with appropriate resource limitations
5. Review any custom scripts or integrations for potential security issues