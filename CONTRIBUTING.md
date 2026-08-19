# Contributing to CompGen-lnc

Thank you for your interest in contributing to CompGen-lnc! You can follow the guidelines below when contributing to this project.

## Summary

- [Contributing to CompGen-lnc](#contributing-to-compgen-lnc)
  - [Summary](#summary)
  - [Requirements](#requirements)
  - [Preparing the Environment](#preparing-the-environment)
  - [Code Convention](#code-convention)
  - [Commit Messages](#commit-messages)
  - [Reporting Issues](#reporting-issues)
  - [Suggestions](#suggestions)

## Requirements

- Python 3.12
- Git
- WSL (if Windows)

## Preparing the Environment

1. Make a fork of the repository (if not administrator) and clone:

    ```bash
    git clone https://github.com/CompGen-lnc/CompGenlnc.git
    cd CompGenlnc
    ```

2. Create and activate and virtual environment:

    ```bash
    python3.12 -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```
    
3. Install the dependencies:

    ```bash
    # Manually by the time being
    pip install <package>
    ```

## Code Convention

- Follow the PEP8 guide for names and format

## Commit Messages

It is recommended to follow the [Conventional Commits](https://www.conventionalcommits.org/) format:

- `feat:` feature addition, modification or removal
- `fix:` bug correction
- `docs:` documentation changes
- `test:` tests addition or corrections
- `refactor:` features changes without altering API
- `chore:` maintenance tasks

## Reporting Issues

Open an issue that includes:

- Description
- Steps to reproduce
- Expected behavior
- Python and OS version
- Relevant traceback or logs

## Suggestions

Issues also can be opened for making suggestions for improving or changing
