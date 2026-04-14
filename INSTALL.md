# Installation

## Requirements

- GRASS GIS ≥ 8.0 (tested with 8.6)
- Python ≥ 3.8 (bundled with GRASS on most platforms)
- No additional Python packages are required — only the standard GRASS
  Python library (`grass.script`) is used.

---

## Method 1 — `g.extension` (recommended)

Once the addon is published to the GRASS Addons repository, install it with
a single command from inside any GRASS session:

```sh
g.extension extension=i.ocean
```

To uninstall:

```sh
g.extension -d extension=i.ocean
```

---

## Method 2 — Manual install from source

### 1. Clone or download the source

```sh
git clone https://github.com/<your-org>/i.ocean.git
cd i.ocean
```

Or copy the `i.ocean/` directory to any convenient location.

### 2. Install via `g.extension` pointing at the local source

Start a GRASS session, then run:

```sh
g.extension extension=i.ocean url=/path/to/i.ocean
```

GRASS copies the script to its addon directory and makes it available as
any other module.

### 3. Alternative: copy manually

```sh
# Locate the user addon script directory
ADDON_DIR=$(grass --config path)/scripts   # system-wide, needs write access
# or the user-level directory:
ADDON_DIR=$HOME/.grass8/addons/scripts

cp i.ocean.py "$ADDON_DIR/i.ocean"
chmod +x "$ADDON_DIR/i.ocean"
```

Make sure `$ADDON_DIR` is on the `PATH` recognised by GRASS
(set `GRASS_ADDON_PATH` if needed).

---

## Method 3 — Build with CMake (grass-addons tree)

Place `i.ocean/` inside a local clone of the
[grass-addons](https://github.com/OSGeo/grass-addons) repository under
`src/imagery/` (or any `src/` subdirectory), then build with CMake:

```sh
cd grass-addons
cmake -B build -S . -DGRASS_INSTALLATION=$(grass --config path)
cmake --build build --target i.ocean
cmake --install build --component i.ocean
```

---

## Verify the installation

Inside a GRASS session:

```sh
i.ocean --help
```

You should see the full parameter list.  Run the self-tests with:

```sh
cd /path/to/i.ocean
python testsuite/test_i_ocean.py
```
