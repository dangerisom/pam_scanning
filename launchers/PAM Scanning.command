#!/bin/bash
# PAM Scanning — double-click launcher (macOS).
#
# First run: creates the 'pam_scanning' conda environment (Python + the app +
# NCBI BLAST+) from environment.yml, which can take a few minutes.
# Every run after that: just opens the app.
#
# Requirement: install Miniforge once (a normal clickable installer):
#   https://conda-forge.org/download/

ENV_NAME="pam_scanning"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"   # this script lives in <repo>/launchers/

pause_and_exit() {
    echo
    read -r -p "Press Return to close this window…" _
    exit "${1:-1}"
}

find_conda() {
    local c
    for c in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" \
             "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
             "/opt/homebrew/Caskroom/miniforge/base/bin/conda" \
             "$(command -v conda 2>/dev/null)"; do
        if [ -n "$c" ] && [ -x "$c" ]; then
            echo "$c"
            return 0
        fi
    done
    return 1
}

CONDA="$(find_conda)" || {
    echo "Could not find conda on this Mac."
    echo "Please install Miniforge first (clickable installer):"
    echo "    https://conda-forge.org/download/"
    pause_and_exit 1
}

# Create the environment the first time only.
if ! "$CONDA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "First-time setup: creating the '$ENV_NAME' environment."
    echo "This downloads NCBI BLAST+ and the app, and may take a few minutes…"
    echo
    ( cd "$REPO" && "$CONDA" env create -f environment.yml ) || {
        echo
        echo "Setup did not finish. Please see the messages above."
        pause_and_exit 1
    }
    echo
    echo "Setup complete."
fi

echo "Starting PAM Scanning…"
exec "$CONDA" run --no-capture-output -n "$ENV_NAME" pam-scan-gui
