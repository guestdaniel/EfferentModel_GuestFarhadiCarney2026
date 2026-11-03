#!/usr/bin/env bash
# Build the model's C shared library (deps/model/libgfc2026.{dylib,so}).
#
#   deps/build_clib.sh            Rebuild only if a source file is newer than the library
#   deps/build_clib.sh --force    Rebuild unconditionally
#   deps/build_clib.sh --all      Also build the 2014 reference library used by the tests
#                                 (external/julia/libzbc2014debug)
#
# Compiler output is shown only on failure; on success, one summary line is printed.
# Exits nonzero if a build fails.
#
# A running Julia session picks up a rebuilt model library on its next model call: the
# package keeps the library open but reopens it when the file's modification time changes.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "$(uname -s)" in Darwin*) EXT=dylib ;; *) EXT=so ;; esac

force=0; all=0
for arg in "$@"; do
    case "$arg" in
        --force) force=1 ;;
        --all) all=1 ;;
        -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "build_clib.sh: unknown argument: $arg" >&2; exit 2 ;;
    esac
done

# build NAME DIR LIB SCRIPT: run SCRIPT in DIR if LIB is missing or older than any source
build() {
    local name="$1" dir="$2" lib="$3" script="$4"
    if [[ $force -eq 0 && -f "$dir/$lib" ]] &&
       [[ -z "$(find "$dir" -maxdepth 1 \( -name '*.c' -o -name '*.h' -o -name '*.hpp' -o -name "$script" \) -newer "$dir/$lib")" ]]; then
        echo "$name: up to date ($dir/$lib)"
        return 0
    fi
    local log t0=$SECONDS
    if log="$(cd "$dir" && bash "$script" 2>&1)"; then
        rm -f "$dir"/*.o
        echo "$name: built $lib in $((SECONDS - t0)) s ($(grep -c 'warning:' <<<"$log" || true) compiler warnings)"
    else
        rm -f "$dir"/*.o
        echo "$log" >&2
        echo "$name: BUILD FAILED" >&2
        return 1
    fi
}

build "model" "$ROOT/deps/model" "libgfc2026.$EXT" compile.sh
if [[ $all -eq 1 ]]; then
    build "2014 reference" "$ROOT/external/julia" "libzbc2014debug.$EXT" shared_library.sh
    # shared_library.sh always writes a .so as well; on macOS it is an unused duplicate
    [[ "$EXT" == dylib ]] && rm -f "$ROOT/external/julia/libzbc2014debug.so"
fi
