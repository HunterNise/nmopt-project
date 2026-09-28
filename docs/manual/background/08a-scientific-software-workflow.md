# 08a · Scientific software workflow: shell, compilation, CMake, and tests

**Background navigation:** [Index](README.md) \
Main chapter: [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) \
Previous: [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md) \
Next: [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md)

## Purpose

Knowing the C++ language is not enough to work comfortably with a scientific C++ project. Source files must be found, compiled, linked with libraries, configured for the local machine, tested, and run from the correct working directory. A large fraction of early frustration comes from not knowing **which layer failed**.

This companion develops the practical software-engineering background needed to build and inspect `nmopt` and deal.II programs. It is intentionally narrower than a general Linux, shell, CMake, or software-engineering course. The aim is to establish a reliable mental model for commands you will repeatedly encounter:

```text
shell
  ↓
source tree
  ↓
CMake configuration
  ↓
generated build system
  ↓
compiler + linker
  ↓
executable / library
  ↓
CTest or application run
```

Once these layers are separated, an error such as “command not found,” “CMake could not find deal.II,” “undefined reference,” or “test failed” stops being one undifferentiated build problem.

Wolfgang Bangerth's MATH 676 video series contains short interactive treatments of the command line, compiling programs, CMake setup, debug versus optimized builds, and scientific-software workflows. Those recordings are useful companions to this document. Some demonstrations are old enough that current deal.II APIs or setup details may differ; the current deal.II tutorial and the current project's build documentation remain the authority for exact commands.

## Before you start

You should have basic programming familiarity and enough C++ to understand source files at the level of [08 · Modern C++ for scientific code](08-modern-cpp-for-scientific-code.md).

No shell expertise is assumed. The examples use a Unix-like command line and Bash conventions because that is the environment in which deal.II and this project are commonly built. The conceptual distinction between source files, object files, libraries, executables, configuration, and tests applies more broadly.

## What you will be able to do

After this companion, you should be able to:

- navigate a source tree from a shell;
- distinguish the current working directory from a path written in a command;
- use basic file-inspection commands without treating the shell as a black box;
- understand quoting, globbing, pipes, redirection, environment variables, and exit status at a practical level;
- explain the compile and link stages of a C++ build;
- distinguish a compiler error from a linker error and both from a run-time failure;
- explain what Make and Ninja do;
- explain what CMake does and why it is not itself a C++ compiler;
- distinguish a source directory from an out-of-source build directory;
- read basic target-oriented CMake code such as `add_executable`, `add_library`, `target_link_libraries`, and `target_compile_features`;
- understand `find_package` and the role of package-discovery hints such as `deal.II_DIR`;
- distinguish Debug and optimized builds;
- run tests through CTest and interpret a failing test as a different layer from a failing build;
- understand CMake presets at a practical level;
- follow the current `nmopt` build helper and know when generic CMake commands are only a fallback;
- locate the layer responsible for common scientific-build failures.

## Roadmap

We first treat the shell as a process launcher with a current directory rather than as a collection of magic commands. We then follow one C++ source file through compilation and linking. That gives enough context to understand why low-level build tools such as Make and Ninja exist.

CMake comes next. We will separate its configure/generate phase from the actual build phase, then introduce target-based dependency descriptions, external-package discovery, build types, and CTest. Finally, we map these general ideas onto the current `nmopt` build workflow and collect common failure patterns.

## 1. The shell launches programs in a working directory

A shell reads a command line, interprets shell syntax, and launches programs. For example,

```bash
ls -l docs
```

contains:

- the program name `ls`;
- the option `-l`;
- the argument `docs`.

The shell locates `ls`, starts it, and passes the arguments to it.

This simple model is more useful than memorizing commands individually.

### 1.1 The current working directory is part of every relative path

The command

```bash
pwd
```

prints the current working directory.

If it prints

```text
/home/user/nmopt-project
```

then the relative path

```text
docs/manual
```

means

```text
/home/user/nmopt-project/docs/manual
```

A path beginning with `/` is absolute:

```text
/tmp/example
```

A path without an initial `/` is interpreted relative to the current working directory unless a tool defines some other path convention explicitly.

### 1.2 `.` and `..` are path components

In Unix-like paths:

```text
.
```

means the current directory, while

```text
..
```

means the parent directory.

Therefore

```bash
cd ..
```

moves to the parent directory, and

```bash
./build.sh
```

means “run the file named `build.sh` from the current directory.”

The explicit `./` matters because the shell normally searches commands through `PATH` rather than automatically treating every executable file in the current directory as a command.

## 2. A small command-line vocabulary goes a long way

The following commands cover much of the read-only navigation needed for a scientific source tree:

```bash
pwd                 # print working directory
ls                  # list directory entries
ls -la               # include hidden entries and details
cd path              # change directory
cat file             # print a small text file
less file            # page through a larger text file
grep pattern file    # find matching lines
find . -name '*.cc'  # recursively find paths by name pattern
```

For directory/file creation and manipulation:

```bash
mkdir directory
cp source destination
mv source destination
rm file
```

`rm` removes files rather than moving them to a desktop trash folder in the usual shell environment. Do not use recursive removal commands casually, especially around source and build trees.

### 2.1 Commands have their own option syntax

The shell does not define what `-l`, `--parallel`, or `--preset` means. The invoked program does.

For example:

```bash
cmake --build build
```

passes `--build` and `build` to CMake, while

```bash
ctest --output-on-failure
```

passes a different option to CTest.

When in doubt, use the program's help:

```bash
cmake --help
ctest --help
```

or its authoritative manual.

## 3. Quoting controls how the shell turns text into arguments

Whitespace normally separates arguments:

```bash
mkdir my results
```

asks `mkdir` to create two directories, `my` and `results`.

To pass a space inside one argument:

```bash
mkdir "my results"
```

### 3.1 Single and double quotes are not identical

Single quotes suppress most shell interpretation:

```bash
echo '$HOME'
```

prints the literal text `$HOME`.

Double quotes preserve spaces while still allowing expansions such as variables:

```bash
echo "$HOME"
```

prints the value of the environment variable `HOME`.

For ordinary project paths without spaces or shell-special characters, no quoting is needed. The important point is that quoting belongs to the shell layer, before the invoked program sees its arguments.

### 3.2 Globs are expanded by the shell

A command such as

```bash
ls *.md
```

uses the glob `*.md`. The shell expands it to matching path names before launching `ls`.

This differs from a regular expression used by tools such as `grep`. A shell glob and a regex are two different pattern languages.

## 4. Redirection and pipes connect process input and output

Most command-line programs have standard input, standard output, and standard error streams.

The shell can redirect standard output:

```bash
program > output.txt
```

This overwrites `output.txt`. To append instead:

```bash
program >> output.txt
```

Standard error can be redirected separately:

```bash
program 2> errors.txt
```

### 4.1 A pipe feeds one program's output into another's input

```bash
cmake --build build 2>&1 | less
```

is a more advanced example: it combines error and ordinary output and feeds the result to `less`.

A simpler and more common pattern is

```bash
find include -name '*.hpp' | less
```

The pipe `|` connects the left process's standard output to the right process's standard input.

The useful mental model is a data flow between processes, not a special feature of the individual programs.

## 5. Exit status lets commands report success or failure

A Unix-like process returns an integer exit status. Conventionally:

```text
0       success
nonzero failure or another non-success condition
```

The shell makes the last status available as `$?`:

```bash
cmake --build build
echo $?
```

The operator `&&` runs the second command only if the first succeeds:

```bash
cmake --build build && ctest --test-dir build
```

This is useful because testing an executable that failed to build would usually be pointless.

In project automation, exit status is what allows scripts and continuous-integration systems to stop when a step fails.

## 6. Environment variables are inherited process configuration

An environment variable is a named string made available to a process and normally inherited by child processes.

Inspect one with

```bash
echo "$PATH"
```

and set/export one with

```bash
export EXAMPLE=value
```

### 6.1 `PATH` controls command lookup

When you type

```bash
cmake
```

without a slash, the shell searches directories listed in `PATH`.

You can ask what executable would be used:

```bash
command -v cmake
command -v c++
```

This is often more useful than assuming which compiler or CMake installation is active.

### 6.2 Project and package variables have narrower meanings

Variables such as

```text
DEAL_II_DIR
CMAKE_PREFIX_PATH
```

may be used by CMake/package-discovery logic, depending on the package and configuration.

They are not built-in shell knowledge about finite elements. They are simply values made available to programs such as CMake, which then interpret them according to their own package-discovery rules.

## 7. A C++ source file becomes machine code in stages

Consider a small program:

```cpp
// main.cc
#include <iostream>

int main()
{
  std::cout << "hello\n";
}
```

A simplified build pipeline is

```text
source + included declarations
          ↓
       compiler
          ↓
      object file
          ↓
        linker
          ↓
      executable
```

Real compiler drivers may internally perform preprocessing and other stages, but this model is enough to diagnose most project-level issues.

### 7.1 Compilation translates one translation unit

A command such as

```bash
c++ -std=c++17 -c main.cc -o main.o
```

compiles `main.cc` into an object file without producing the final executable.

The compiler checks syntax, names, types, template instantiations, and other language rules for that translation unit.

Typical compile-time failures include:

```text
header not found
unknown type
no matching function
invalid conversion
template instantiation error
syntax error
```

### 7.2 Linking combines object files and libraries

A second command could be

```bash
c++ main.o -o example
```

The linker resolves references among compiled object files and libraries and produces the executable.

A classic linker error is

```text
undefined reference to ...
```

This often means that a declaration was visible during compilation, but the corresponding compiled definition was not supplied to the linker.

This distinction is crucial:

```text
compiler error
    problem understanding one translation unit

linker error
    compiled pieces exist, but required definitions cannot be connected
```

### 7.3 Running is a third layer

If linking succeeds, the program may still fail when executed:

```bash
./example
```

Possible run-time failures include invalid input, failed assertions, exceptions, numerical nonconvergence, segmentation faults, or missing run-time data files.

Do not treat these as compile errors merely because they occur during a “build-and-run” workflow.

## 8. Libraries package compiled code for reuse

Instead of linking every implementation directly into every executable, compiled code can be packaged as libraries.

At a high level:

```text
application object files
       +
library code
       ↓
     linker
       ↓
 executable
```

Libraries may be static or shared/dynamic. The low-level differences matter in some deployment situations but are not necessary for the current background path.

The important conceptual distinction is between:

- **headers**, which make declarations and template definitions visible during compilation;
- **libraries**, which provide compiled definitions to the linker;
- **executables**, which are runnable program products.

A large library such as deal.II supplies all three kinds of material through its installed development package: headers, compiled libraries, and CMake metadata describing how a client target should be configured and linked.

## 9. Why build tools exist

Manually compiling a multi-file project quickly becomes tedious. Suppose an executable depends on ten source files and several libraries. If one source file changes, we would like to rebuild only what is affected rather than compile everything from scratch.

A build tool represents dependencies such as

```text
main.cc ─────► main.o ──┐
solver.cc ───► solver.o ├─► executable
model.cc ─────► model.o ─┘
```

It determines which outputs are out of date and runs the necessary commands.

This dependency-graph role is the core idea behind Make and Ninja.

## 10. Make describes targets, prerequisites, and recipes

A minimal Makefile could contain

```make
example: main.o model.o
	c++ main.o model.o -o example

main.o: main.cc
	c++ -std=c++17 -c main.cc -o main.o

model.o: model.cc model.hpp
	c++ -std=c++17 -c model.cc -o model.o
```

Each rule has:

```text
target: prerequisites
    recipe
```

If a prerequisite is newer than its target, or the target does not exist, Make can run the recipe.

### 10.1 The concept matters more here than handwritten Make syntax

Real C++ projects need platform-specific compiler flags, include directories, library dependencies, generated files, tests, and external packages. Maintaining all of this by hand in Makefiles is possible but quickly becomes a project of its own.

For `nmopt` and deal.II applications, CMake is the higher-level build-system description. CMake may generate a Make-based or Ninja-based build system underneath.

Therefore you should understand what Make *does* without assuming that the project expects you to edit Makefiles directly.

## 11. Ninja is another low-level build executor

Ninja also consumes a build dependency graph and executes commands needed to bring targets up to date. It is deliberately a relatively low-level tool and is commonly driven by a generator such as CMake rather than by manually maintained build files.

The conceptual relation is:

```text
CMake project description
       ↓ generate
Makefiles or build.ninja
       ↓ execute
compiler / linker commands
```

You normally edit `CMakeLists.txt`, not generated `build.ninja` files or generated Makefiles.

The choice between Make and Ninja does not change the C++ mathematics or finite-element method. It changes how the generated build graph is executed.

## 12. CMake configures and generates a build system

A small CMake project might contain

```cmake
cmake_minimum_required(VERSION 3.20)
project(example LANGUAGES CXX)

add_executable(example main.cc model.cc)
target_compile_features(example PRIVATE cxx_std_17)
```

CMake reads this project description and generates build files for a selected generator.

A modern out-of-source configuration looks like

```bash
cmake -S . -B build
```

where

- `-S .` selects the source tree;
- `-B build` selects a separate build tree.

The actual build is then

```bash
cmake --build build
```

These are different phases.

### 12.1 Configure/generate versus build

During configuration/generation, CMake may:

- detect the compiler;
- inspect platform properties;
- find external packages;
- evaluate project options;
- generate build rules.

During the build, the generated build tool runs compiler/linker commands.

Therefore an error such as

```text
Could not find a package configuration file provided by "deal.II"
```

is a **CMake configuration/discovery failure**, not a C++ compiler error.

## 13. Keep source and build trees separate

An out-of-source build keeps generated material away from tracked source files:

```text
project/
├── CMakeLists.txt
├── include/
├── src/
├── tests/
└── build/
    ├── CMakeCache.txt
    ├── build.ninja or Makefiles
    ├── object files
    └── executables
```

This separation has practical advantages:

- generated files are easier to remove or replace;
- several build configurations can coexist;
- source-control status stays cleaner;
- a broken cache can be discarded without touching source files.

A build directory is not merely a convenient output folder. It contains CMake's configured state for a particular build tree.

### 13.1 Do not casually move or reuse a build cache across incompatible configurations

CMake records absolute paths, compiler choices, generator choices, package locations, and options. Reusing a build directory created for a different source location or generator can lead to confusing failures.

When a build tree is disposable, recreating that generated directory is often clearer than trying to manually edit its internal cache files.

Project-specific rules may impose stricter guidance; `nmopt` does, as discussed below.

## 14. CMake is target-oriented

Modern CMake encourages thinking in terms of targets and relationships among targets rather than manually concatenating compiler flags.

For example:

```cmake
add_library(model model.cc)
target_include_directories(model PUBLIC include)
target_compile_features(model PUBLIC cxx_std_17)

add_executable(app main.cc)
target_link_libraries(app PRIVATE model)
```

The executable `app` depends on `model`.

CMake propagates usage requirements according to keywords such as `PRIVATE`, `PUBLIC`, and `INTERFACE`.

A useful first interpretation is:

```text
PRIVATE
    needed to build/use this target internally
    not part of what dependents inherit

PUBLIC
    needed by this target and by its consumers

INTERFACE
    usage requirement for consumers, not for compiling this target itself
```

The exact propagation rules can be learned when a concrete project needs them. The important shift is from “global compiler flags” to “requirements attached to targets.”

## 15. Interface libraries can carry requirements without compiled object code

CMake supports

```cmake
add_library(my_interface INTERFACE)
```

An interface target does not compile ordinary source files into a library artifact. It can still carry include directories, compile features, definitions, and transitive link requirements.

This is useful for header-only APIs or project-wide build policy.

The `nmopt` root build uses interface targets for precisely these kinds of responsibilities: one public contract target exposes headers and the C++17 requirement, while another carries project compiler-warning and sanitizer policy.

## 16. `find_package` connects a project to installed libraries

A client project should not normally guess all include directories and linker flags for a complex dependency by hand. CMake packages can provide configuration metadata and imported targets.

A project may contain

```cmake
find_package(deal.II REQUIRED)
```

or a library-specific variation. CMake searches for package metadata and, if successful, makes package-defined variables or targets available.

### 16.1 Discovery hints are not the dependency itself

A setting such as

```text
deal.II_DIR=/path/to/dealii/lib/cmake/deal.II
```

is a hint telling CMake where to find package configuration. It is not a substitute for an installed or built deal.II library.

Similarly, `CMAKE_PREFIX_PATH` can tell CMake about prefixes under which packages are installed.

If package discovery fails, diagnose the configuration layer:

1. Is the dependency installed?
2. Is the expected package configuration present?
3. Is CMake searching the right prefix/configuration directory?
4. Is the requested dependency enabled in this build profile?

Do not start changing C++ source code before answering those questions.

## 17. Debug and optimized builds serve different purposes

A Debug-oriented build typically favors diagnostics and debuggability. An optimized/Release-oriented build favors run-time performance.

The exact compiler flags are project-specific, but common consequences are:

```text
Debug
    easier source-level debugging
    assertions/checks may be enabled
    less aggressive optimization
    slower numerical execution

optimized build
    compiler optimizations enabled
    generated code may be much faster
    source-level debugging can be less direct
```

For large numerical simulations, performance differences can be dramatic. Timing a Debug build is therefore not reliable benchmark evidence.

At the same time, developing only in optimized mode may hide useful diagnostics or make failures harder to localize.

Bangerth's Lecture 18 is a useful interactive discussion of this distinction in a deal.II context.

## 18. Compiler warnings and sanitizers are additional diagnostic layers

A successful compilation only means the compiler accepted the program according to the selected language and diagnostic settings. Warnings can identify suspicious constructs that are legal C++ but often accidental.

Projects may intentionally treat warnings as errors in development builds to keep new warnings from accumulating.

Sanitizers instrument a program to detect certain classes of run-time errors, such as invalid memory accesses or undefined behavior. They are not proofs of correctness: they can only diagnose problems exercised by the run.

The important workflow idea is that build profiles may exist for different verification purposes rather than only “fast” versus “slow.”

## 19. CMake presets name reproducible configurations

A project can provide `CMakePresets.json` to give names to configure/build/test setups.

Instead of remembering a long list of cache variables, a user may run a command such as

```bash
cmake --preset debug
```

if the project defines that preset.

Presets can specify generator, binary directory, cache variables, environment settings, and related workflow information. Build and test presets can similarly name the matching later stages.

A preset is project configuration expressed declaratively. It is different from a shell alias that merely abbreviates a command string.

## 20. CTest runs tests registered by CMake

A CMake project can register tests. CTest is the companion tool used to execute and report them.

A generic workflow is

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

The distinction is:

```text
build succeeds
    executable/library targets were produced

CTest succeeds
    registered test processes returned success under the test configuration
```

A test failure does not imply that compilation failed. It means a built program or script did not satisfy the registered test.

### 20.1 Focused test selection is part of an efficient development loop

CTest can filter by names or labels. Large projects often use this to run a small relevant subset while developing and a broader suite before integration.

A short feedback loop is not the same thing as skipping the full verification gate. It means using the cheapest relevant check first, then expanding verification when the change warrants it.

## 21. The current `nmopt` build workflow wraps the lower-level tools

The concepts above explain what happens underneath `nmopt`, but the project deliberately provides a preferred entry point. For ordinary development and verification, use the root-level helper rather than reconstructing the raw CMake commands by hand.

The fast backend-neutral pipeline is

```bash
./build.sh pipeline debug-neutral
```

Changes that touch deal.II code, compiler/lowering behavior, or numerical deal.II functionality also use

```bash
./build.sh pipeline debug-dealii
```

The helper delegates to checked-in CMake presets and enforces project-local build policy.

### 21.1 Local machine configuration is intentionally separate

The helper expects an existing

```text
build.local.conf
```

for build/pipeline operations. That file carries machine-local choices such as allowed parallelism. If it is missing, initialize it explicitly with

```bash
./build.sh init-config
```

before using the build or pipeline actions. Keeping this step explicit prevents machine-local settings from being invented silently.

This illustrates a common scientific-software pattern:

```text
repository-owned configuration
    reproducible project choices

machine-local configuration
    installation paths / resources / local limits
```

Keeping them separate prevents one developer's machine assumptions from becoming repository-wide defaults.

### 21.2 Named build profiles live under `build/`

The current persistent profiles include

```text
build/debug-neutral/
build/debug-dealii/
build/sanitize-neutral/
build/release-dealii/
```

The root `build/` directory is only a container for these named configurations. The project explicitly avoids configuring CMake directly into `build/` itself.

This is a concrete instance of the out-of-source-build principle developed earlier, extended to multiple simultaneous configurations.

### 21.3 Manual CMake commands are a fallback, not the first project interface

When lower-level control is genuinely needed, a deal.II configuration can be driven manually through the checked-in presets:

```bash
cmake --preset debug-dealii
cmake --build --preset debug-dealii --parallel 1
ctest --preset debug-dealii
```

The `--parallel 1` restriction is a project-specific resource policy for manual deal.II builds. deal.II-heavy translation units can consume substantial memory, so the project does not permit uncontrolled parallel compilation in that fallback path.

This is a good example of why generic CMake knowledge and project instructions must both be respected. The generic tool tells you what is possible; the project tells you which usage is supported and safe here.

## 22. Reading the `nmopt` root CMake file with the target model

Several lines of the current root configuration become straightforward once the target mental model is established.

The project begins conceptually with

```cmake
project(nmopt LANGUAGES CXX)
```

and defines an interface target for the public contract:

```cmake
add_library(nmopt_contract INTERFACE)
target_include_directories(nmopt_contract INTERFACE include)
target_compile_features(nmopt_contract INTERFACE cxx_std_17)
```

This says, in effect:

```text
nmopt_contract
    compiled object code?       no, interface target
    public header location?     include/
    language requirement?       C++17
```

Executable helper functions then create a target and link it to the contract and build-policy targets. In the deal.II-enabled path, targets also receive deal.II's required setup.

You do not need to memorize the CMake language to understand the architecture. Read each command as adding a node or an edge to the build graph.

## 23. deal.II tutorials use the same layers, even when their commands are simpler

A tutorial program may show a compact sequence resembling

```bash
cmake .
make
make run
```

or an out-of-source variant. The teaching purpose is to expose the minimum path from an example's CMake project to an executable.

For current work, prefer the build procedure documented by the project you are actually using. In particular, do not copy an old tutorial's in-source build convention into `nmopt` simply because the mathematical example is similar.

Bangerth's Lecture 2.91 is useful precisely because it presents compilation and linking interactively before hiding them behind CMake. Lecture 8.01 then demonstrates direct CMake project setup. The current deal.II tutorial should be used to resolve any API or build-command drift in the older videos.

## 24. A diagnostic decision tree is more useful than random command changes

When something fails, first classify the layer.

### 24.1 “command not found”

Example:

```text
cmake: command not found
```

This is a shell/tool-installation or `PATH` problem. CMake has not even begun to configure the project.

Check:

```bash
command -v cmake
cmake --version
```

### 24.2 CMake configure failure

Example:

```text
Could not find deal.IIConfig.cmake
```

The compiler may be perfectly healthy. CMake cannot discover the external package.

Inspect the dependency installation and package-discovery configuration.

### 24.3 Compiler failure

Example:

```text
fatal error: some_header.hpp: No such file or directory
```

or a type/template diagnostic.

Now CMake successfully generated a compiler invocation, but that translation unit could not be compiled.

Read the first relevant project-owned error, the include paths or definitions involved, and the concrete template arguments if applicable.

### 24.4 Linker failure

Example:

```text
undefined reference to `Model::solve(...)'
```

Compilation succeeded far enough to produce object files. The linker cannot find a required compiled definition or compatible library symbol.

Inspect target linkage and whether the implementation source/library is part of the build.

### 24.5 Test failure

Example:

```text
90% tests passed, 1 tests failed
```

The build succeeded. A test process failed semantically or at run time.

Use CTest's failure output or run the focused test according to project instructions.

### 24.6 Application run failure

A built and tested executable may still reject an application configuration, fail to converge for particular numerical data, or be launched from a directory where relative input paths do not exist.

At this point the issue belongs to application/run semantics, not automatically to CMake.

## 25. Reproducibility requires recording configuration, not just source code

Two users can check out the same commit yet build materially different programs if they use different:

- compilers;
- library versions;
- CMake options;
- optimization modes;
- external package installations;
- environment variables.

A serious scientific workflow therefore records enough configuration to reconstruct what was built and run.

CMake presets, project build profiles, dependency metadata, parameter files, and generated run manifests all contribute to this goal at different layers.

This is why the build system should be treated as part of the scientific software rather than as an incidental command typed once before doing the “real” mathematics.

## 26. Scope frontier

This companion gives only the workflow needed to become productive around `nmopt` and deal.II. We stop before:

- shell scripting as a programming language;
- `awk`, `sed`, and advanced text-processing pipelines;
- package managers and Linux distribution administration;
- Git and version-control workflows;
- containers and environment managers;
- continuous-integration service configuration;
- CMake generator expressions and advanced install/export packaging;
- cross-compilation and toolchain files;
- ABI compatibility and shared-library deployment;
- debugger command languages;
- profiling tools;
- HPC batch schedulers;
- MPI launchers and parallel-runtime configuration.

Several of these become important in larger computational projects, but none is required to understand the current serial finite-element and optimization path.

## Where to go next

- **Default continuation:** [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) makes the source/build concepts concrete in a deal.II finite-element application.
- **When working directly in `nmopt`:** the repository's `.agents/build.md`, checked-in CMake presets, root `build.sh`, and project README remain authoritative for exact project commands. This companion explains the machinery underneath those instructions rather than replacing them.
- **Into the `nmopt` manual:** [manual 12 · Authoring and using compiled problems](../concepts/12-authoring-and-using-compiled-problems.md) and [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md) are the most direct continuations once you can configure, build, test, and run the project comfortably.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

- **Wolfgang Bangerth, MATH 676 video lectures.**
  - Lecture 2.9 for a brief command-line introduction.
  - Lecture 2.91 for compiling programs.
  - Lecture 8.01 for direct CMake project setup.
  - Lecture 18 for Debug versus optimized mode.
  - Lectures 42–43 for scientific-computing workflows and large-software development.
  - Some recordings predate current deal.II releases; use the current tutorial for present-day APIs and setup.
- **CMake official documentation.**
  - Getting Started and the CMake Tutorial for the introductory route.
  - `cmake-buildsystem(7)`, `cmake-presets(7)`, and the CTest manuals when targets, presets, and registered tests become concrete needs.
- **GNU Bash Reference Manual.**
  - Authoritative reference for shell parsing, quoting, expansion, redirection, variables, pipelines, and command execution.
- **GNU Make manual.**
  - Targets, prerequisites, recipes, dependency updating, and variables when Make is the generated build executor.
- **Ninja manual.**
  - Ninja's role as a low-level build executor, normally driven here by CMake.
- **deal.II tutorial programs.**
  - Current authority for building and running deal.II tutorial examples; [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) identifies the most relevant steps.
