# Dependencies and environment

This document records the software needed to configure, build, test, run, and
post-process `nmopt-project`. It is a human installation reference; it is not
an additional build manifest or a lockfile.

The version output below is a snapshot of the reference environment collected
on 2026-08-19. It documents a verified environment, not exact version pins.
The repository currently enforces only CMake 3.20 or newer and the C++17
language level. The checked-in CMake presets additionally use Ninja.

> [!WARNING]
> The project is developed and tested in a Linux environment. Windows users
> should install WSL2 with an Ubuntu distribution and run the project inside
> the Linux environment, preferably from the WSL filesystem rather than a
> mounted Windows directory.

## Dependency summary

### Linux environment and shell

| Dependency | Repository requirement or verified version | Required for |
| --- | --- | --- |
| Linux environment | Required platform | All supported workflows |
| Linux distribution | Verified Ubuntu `24.04.1 LTS` | Reference environment |
| Bash and standard Unix tools | Standard Linux userland; verified Bash `5.2.21` | `tools/run_chapter6.sh` |

### Git and development tools

| Dependency | Repository requirement or verified version | Required for |
| --- | --- | --- |
| Git | Verified `2.43.0` | Cloning the repository and recording run provenance |

### Build system

| Dependency | Repository requirement or verified version | Required for |
| --- | --- | --- |
| CMake | Minimum `3.20`; verified `3.28.3` | Configuring and building |
| Ninja | Required by the checked-in presets; verified `1.11.1` | Configuring and building with the presets |

### C++ compiler and deal.II

| Dependency | Repository requirement or verified version | Required for |
| --- | --- | --- |
| C++ compiler | C++17; verified GCC/G++ `13.3.0` | Compiling the library, runner, and tests |
| deal.II | Optional for neutral builds; verified `9.5.1` | `debug-dealii`, `nmopt_runner`, and Chapter 6 execution |

### Python and post-processing libraries

| Dependency | Repository requirement or verified version | Required for |
| --- | --- | --- |
| Python 3 | Verified `3.12.3` | Post-processing and reporting tools |
| NumPy | Verified `1.26.4` | Python field and geometry processing |
| Matplotlib | Verified `3.6.3` | Rendering PNG and SVG field plots |
| meshio | Verified `5.3.5` | Reading and writing native mesh/field files |

The backend-neutral profile does not discover deal.II and does not require the
Python packages. The deal.II profile is needed for the runner and Chapter 6
application tests. Python is needed when generating post-processing plots or
reports; the C++ build itself does not import Python.

The project does not directly request MPI, BLAS, or other numerical libraries
from CMake. Any transitive libraries required by deal.II belong to the
selected deal.II installation method.

## Official sources

- [Microsoft WSL installation](https://learn.microsoft.com/en-us/windows/wsl/install)
- [Git downloads](https://git-scm.com/downloads)
- [CMake downloads](https://cmake.org/download/)
- [Ninja](https://ninja-build.org/)
- [GCC](https://gcc.gnu.org/)
- [deal.II](https://www.dealii.org/), [deal.II CMake documentation](https://www.dealii.org/developer/users/cmake_dealii.html)
- [Python downloads](https://www.python.org/downloads/)
- [NumPy installation](https://numpy.org/install/)
- [Matplotlib installation](https://matplotlib.org/stable/install/index.html)
- [meshio on PyPI](https://pypi.org/project/meshio/)
- [uv installation](https://docs.astral.sh/uv/getting-started/installation/), [uv Python environments](https://docs.astral.sh/uv/pip/environments/)

## Installation on Ubuntu or WSL

The following installs the base repository toolchain on Ubuntu, including the
optional deal.II backend:

```bash
sudo apt update
sudo apt install \
  git \
  build-essential \
  cmake \
  ninja-build \
  libdeal.ii-dev \
  python3
```

For backend-neutral development only, `libdeal.ii-dev` can be omitted. The
Ubuntu package supplies the deal.II CMake package expected by the
`debug-dealii` preset. If deal.II is installed elsewhere, pass its CMake
package directory when configuring:

```bash
cmake --preset debug-dealii \
  -Ddeal.II_DIR=<path-to-deal.II-cmake-package>
```

The post-processing packages can either be installed into the system Python
environment with Ubuntu packages or kept separate in a user-owned `uv`
virtual environment.

### Option A: Ubuntu Python packages

This is the simplest option when the project is being used only on Ubuntu or
WSL:

```bash
sudo apt install python3-numpy python3-matplotlib python3-meshio
```

These packages are managed by Ubuntu and are available to the system
interpreter.

### Option B: uv-managed virtual environment

Install [`uv`](https://docs.astral.sh/uv/getting-started/installation/) using
its official instructions, then create and activate a user-owned environment
outside the repository:

```bash
export NMOPT_VENV="${NMOPT_VENV:-$HOME/.venvs/nmopt}"
mkdir -p "$(dirname "$NMOPT_VENV")"
uv venv "$NMOPT_VENV"
source "$NMOPT_VENV/bin/activate"
uv pip install numpy matplotlib meshio
```

The `source` command enters the environment in the current shell. While it is
active, `python` and `python3` refer to the environment's interpreter:

```bash
tools/run_chapter6.sh \
  --benchmark b1 \
  --refinement 1
```

Use `deactivate` to leave the environment. `uv run` is useful for one-off
commands, but it is not required here: `tools/run_chapter6.sh` uses the
activated `python3` by default. If you prefer not to activate the environment,
pass its interpreter explicitly with `--python "$NMOPT_VENV/bin/python"`.

## Checking the environment

Run these commands from any directory. The paths intentionally are not part
of the recorded output so that the checks remain useful on different Linux
installations. The `$` prefix denotes a command entered in a POSIX shell.

### Linux environment and shell

```text
$ uname -s
Linux

$ . /etc/os-release && printf '%s\n' "$PRETTY_NAME"
Ubuntu 24.04.1 LTS

$ bash --version | head -1
GNU bash, version 5.2.21(1)-release (x86_64-pc-linux-gnu)
```

### Git and development tools

```text
$ git --version
git version 2.43.0
```

### Build system

```text
$ cmake --version
cmake version 3.28.3

$ ninja --version
1.11.1
```

### C++ compiler and deal.II

```text
$ c++ --version
c++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0

$ cmake --find-package -DNAME=deal.II -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST
deal.II found.

$ dpkg-query -W -f='${Version}\n' libdeal.ii-dev
9.5.1-2build3

$ deal_ii_config="$(find /usr/local /usr -type f -name deal.IIConfig.cmake -print -quit)"
$ sed -n 's/^set(DEAL_II_VERSION "\(.*\)")/\1/p' "$deal_ii_config"
9.5.1
```

The `dpkg-query` command reports the Ubuntu package version, including its
distribution revision. The `find`/`sed` command reads the upstream deal.II
version from the CMake package and searches the common `/usr/local` and `/usr`
prefixes; adjust those search roots for a custom installation. The repository
does not currently encode an exact deal.II version constraint; compatibility
is established by successfully configuring and building the deal.II preset.

### Python and post-processing libraries

```text
$ python3 --version
Python 3.12.3

$ python3 -c '
from importlib.metadata import version

for package in ("numpy", "matplotlib", "meshio"):
    print(package, version(package))
'
numpy 1.26.4
matplotlib 3.6.3
meshio 5.3.5
```

## Build and test

Run the commands from the repository root. The backend-neutral profile is the
fastest way to check the repository without deal.II:

```bash
cmake --preset debug-neutral
cmake --build --preset debug-neutral
ctest --preset debug-neutral
```

The full deal.II profile builds the runner and deal.II-backed tests:

```bash
cmake --preset debug-dealii
cmake --build --preset debug-dealii
ctest --preset debug-dealii
```

On memory-constrained machines, limit the deal.II build to one job:

```bash
cmake --build --preset debug-dealii --parallel 1
```

The sanitizer profile is backend-neutral:

```bash
cmake --preset sanitize-neutral
cmake --build --preset sanitize-neutral
ctest --preset sanitize-neutral
```

The `release-dealii` profile is reserved for optimized verification and
source-scale Chapter 6 timing or reproduction runs. It is not required for a
normal installation check.

## Run Chapter 6 applications

After building `debug-dealii`, list the registered application entries:

```bash
build/debug-dealii/bin/nmopt_runner --list
```

The repository helper runs a development refinement, post-processes the native
fields, and writes a report under the ignored `runs/` directory:

```bash
tools/run_chapter6.sh \
  --benchmark b1 \
  --refinement 1
```

The script requires the built `nmopt_runner`, Bash, Git, Python, NumPy,
Matplotlib, and meshio. Use `--python <path>` when the post-processing
packages are installed in a separate virtual environment.

For run schemas, artifact paths, report contents, and direct post-processing
commands, see the [application execution reference](docs/reference/application-execution.md).
