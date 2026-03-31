# Plan: Integrate LISFLOOD-FP Source Directly with Auto-Compile

## Context
The autocalib framework previously required a git submodule for LISFLOOD-FP 8.2 source code. To simplify user experience, we've integrated the source directly into the repository so users can just clone and run self-calibration without submodule initialization.

## Approach

### 1. Integrated LISFLOOD-FP 8.2 source directly
```
src/lisflood-fp/ → Complete LISFLOOD-FP 8.2 source tree
```

### 2. Updated compilation helper (`scripts/build_lisflood.py`)
- Detect platform (Linux/macOS/Windows)
- Run cmake + make in `src/lisflood-fp/build/`
- Output binary at `src/lisflood-fp/build/lisflood`
- Support `--cuda` flag for GPU build
- Skip if binary already exists (unless `--force` or `--clean`)

### 3. Updated adapter auto-compile functionality
- Modified `models/lisflood82/adapter.py` `_resolve_executable()` method
- Now resolves `executable: "auto"` to `src/lisflood-fp/build/lisflood`
- Automatically compiles from integrated source if binary missing
- Removed submodule-specific paths and error messages

### 4. Updated main.py configuration loading
- No changes needed - adapter handles executable resolution
- Maintains backward compatibility with explicit executable paths

### 5. Updated config examples
- Default `executable` remains `"auto"` (now points to integrated source)
- Documentation updated to reflect direct source integration

### 6. Updated `.gitignore`
- Changed from ignoring `vendor/lisflood-fp/build/` to `src/lisflood-fp/build/`

### 7. Updated documentation
- CLAUDE.md: Updated to reflect direct source integration approach
- Other documentation files maintained as-is (general build process unchanged)

## Files Modified

| File | Action |
|------|--------|
| `src/lisflood-fp/` | New — integrated LISFLOOD-FP 8.2 source |
| `scripts/build_lisflood.py` | Modified — updated paths from vendor/ to src/ |
| `models/lisflood82/adapter.py` | Modified — updated `_resolve_executable()` for direct source |
| `.gitignore` | Modified — changed build artifact path |
| `plan.md` | Updated — reflects completed integration approach |

## Verification
1. `git clone` repository (no submodule initialization needed)
2. `python scripts/build_lisflood.py` compiles successfully from integrated source
3. `python main.py --config config/project_lisflood82.yaml` with `executable: "auto"` auto-compiles and runs
4. `python main.py --config config/project_lisflood82.yaml` with explicit path still works (backward compatible)
5. No `.gitmodules` file or submodule references remain
