# Contributing to viva-equilibrator

## Development setup

uv is required (`brew install uv` or `pip install uv`).

    uv venv .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest

The real-engine tests skip automatically when the `equilibrator_api` offline
cache is unavailable, so the suite stays green on a machine without it.

## Releasing to PyPI

Cut a GitHub Release tagged `v<VERSION>`; `.github/workflows/release.yml`
publishes to PyPI via OIDC trusted publishing (no tokens). The one-time PyPI
pending-publisher setup is described in that workflow's header.
