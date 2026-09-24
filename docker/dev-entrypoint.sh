#!/bin/bash
#
# Entrypoint for the SCM development container
#

# The source tree is bind-mounted from the host, so its owner may not match
# the container user; let git operate on it regardless.
git config --global --add safe.directory '*' 2>/dev/null || true

exec "$@"
