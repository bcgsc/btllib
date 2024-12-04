[![GitHub release (latest by date)](https://img.shields.io/github/v/release/bcgsc/btllib?label=Release)](https://github.com/bcgsc/btllib/releases/latest)
[![Conda](https://img.shields.io/conda/dn/bioconda/btllib?label=Conda%20downloads)](https://anaconda.org/bioconda/btllib)
[![Build Status](https://dev.azure.com/bcgsc/btl_public/_apis/build/status/bcgsc.btllib)](https://dev.azure.com/bcgsc/btl_public/_build/latest?definitionId=1)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04720/status.svg)](https://doi.org/10.21105/joss.04720)

[Bioinformatics Technology Lab](http://www.birollab.ca/) common code library in C++ with Python wrappers.

Platforms
---
- Linux
- MacOS

Installation for users
---
The recommended way is to download using [Conda](https://docs.conda.io/en/latest/) package manager:  
`conda install -c bioconda -c conda-forge btllib`

Alternatively, you can compile the code from source. Download `btllib-$VERSION.tar.gz` from the GitHub [latest release](https://github.com/bcgsc/btllib/releases/latest) where `$VERSION` is the latest btllib version and do the following:
- `tar xzf btllib-$VERSION.tar.gz` to extract the source code.
- Have the dependencies ready:
  * GCC 6+ or Clang 5+ (with OpenMP and C++17 support)
  * Python 3.9+
  * Meson and Ninja Python3 packages, CMake (If not available, these will be automatically installed to a temporary directory.)
- Run `btllib/compile`
  * This will install btllib in the `btllib/install` directory. You can provide the `--prefix` parameter to change this.
  * The C++ compiler must be the same as the one used for compiling Python. E.g. if you installed Python using a package manager, you should use the C++ compiler from the same package manager. You can change the compiler by exporting the `CXX` environment variable to point to the compiler before running `btllib/compile`.
  * You can optionally run `python3 -m pip install $PREFIX/lib/btllib/python` afterwards to install the Python package. The Python wrappers are usable even without this step. `$PREFIX` is the path where btllib is installed.

Using the library
---
- Run time dependencies:
  * SAMtools for reading SAM, BAM, and CRAM files.
  * gzip, tar, pigz, bzip2, xz, lrzip, zip, and/or 7zip for compressing/decompressing files. Not all of these are necessary, only the ones whose compressions you'll be using.
    * Note that lrzip is not available on the btllib conda osx-arm64 build
  * wget for downloading sequences from a URL.
- Building C++ code (`$PREFIX` is the path where btllib is installed):
  * Link your code with `$PREFIX/lib/libbtllib.a` (pass `-L $PREFIX/lib -l btllib` flags to the compiler).
      * You can do so by typing the following in your console:
          * `export CPPFLAGS="-isystem /path/to/btllib/install/include $CPPFLAGS"`
          * `export LDFLAGS="-L/path/to/btllib/install//lib -lbtllib $LDFLAGS"`
  * `#include` any header from the `$PREFIX/include` directory (pass `-I $PREFIX/include` flag to the compiler).
  * `btllib` uses `C++11` features, so that standard should be enabled at a minimum.
- Running Python code:
  * The Python used to import btllib _must_ be the same as the one used to compile the library. Specifically, btllib uses `python3-config` to determine the flags used for compilation. Running `python3-config --exec-prefix` will give the path to the Python installation that needs to be used. The `python3` executable can be found at `$(python3-config --exec-prefix)/bin/python3`.
  * The wrappers correspond one-to-one with C++ code so any functions and classes can be used under the same name. The only exceptions are nested classes which are prefixed with outer class name (e.g. `btllib::SeqReader::Flag` in C++ versus `btllib.SeqReaderFlag` in Python), and (Kmer)CountingBloomFilter which provides `CountingBloomFilter8`, `CountingBloomFilter16`, `CountingBloomFilter32`, `KmerCountingBloomFilter8`, `KmerCountingBloomFilter16`, `CountingBloomFilter32` with counters 8, 16, and 32 bits wide.
  * If you compiled btllib from source code and didn't install the Python wrappers, you can use `PYTHONPATH` environment variable or `sys.path.append()` in your Python code to include `$PREFIX/lib/btllib/python/btllib` directory to make btllib available to the interpreter.
  * Include the library with `import btllib`
- Executables
  * btllib generated executables can be found in `$PREFIX/bin` directory. Append that path to the `PATH` environment variable to make it available to your shell.

Documentation
---
[Docs page](https://bcgsc.github.io/btllib/)

For btllib developers
---
- Initial setup:
  * `git clone --recurse-submodules https://github.com/bcgsc/btllib` in order to obtain all the code.
  * In `btllib` dir, run `meson build` to create a build directory.
- Every time you want to run tests, in the `build` dir:
  * `ninja wrap` to regenerate wrappers.
  * `ninja test` to build wrappers and tests, and run tests.
- Before making a pull request, in the `build` dir:
  * `ninja quality-assurance` to make sure all CI tests pass.
  * Make a commit after the above step, in case it has made any changes to wrappers or formatting. Don't commit the changes made to the `sdsl-lite` subproject. Meson config file adjusts the `sdsl-lite` config in order for it to work for `btllib`, but this is done ad hoc and is not necessary to be committed. By doing it ad hoc we keep a list of differences compared to the upstream repository.
- Before making a release, in the `build` dir:
  * Do the same as for a pull request and
  * `ninja docs` to regenerate docs to reflect the release and then commit the changes.
  * `meson dist --allow-dirty` to generate a self-contained package based on the last commit. `--allow-dirty` permits making a distributable with uncommited changes. This is necessary as `sdsl-lite` dependency has ad hoc changes made during the build process. The resulting distributable will be compressed with xz. For easier use, decompress it and then compress with gzip. Attach the resulting file to the release.

The following are all the available `ninja` commands which can be run within `build` directory:
- `ninja clang-format` formats the whitespace in code (requires clang-format 8+).
- `ninja wrap` wraps C++ code for Python (requires SWIG ≥4.0 and <4.3).
- `ninja clang-tidy` runs clang-tidy on C++ code and makes sure it passes (requires clang-tidy 8+).
- `ninja` builds the tests and wrapper libraries / makes sure they compile.
- `ninja test` runs the tests.
- `ninja code-coverage` assures code coverage threshold is satisfied. (requires gcovr 3.3+)
- `ninja sanitize-undefined` runs undefined sanitization.
- `ninja test-wrappers` tests whether wrappers work.
- `ninja docs` generates code documentation from comments (requires Doxygen).
- `ninja quality-assurance` runs `clang-format`, `wrap`, `clang-tidy`, `test`, `code-coverage`, `sanitize-undefined`, and `test-wrappers`. These are all checked at the CI test.

Credits
---
- Author: [Vladimir Nikolic](https://github.com/vlad0x00)
- Components:
  - [Hamid Mohamadi](https://github.com/mohamadi) and [Parham Kazemi](https://github.com/parham-k) for [ntHash](https://github.com/bcgsc/ntHash)
  - [Justin Chu](https://github.com/JustinChu) for [MIBloomFilter](https://github.com/bcgsc/btl_bloomfilter)
  - [Johnathan Wong](https://github.com/jwcodee) for [aaHash](https://github.com/bcgsc/btllib)
- Included dependencies:
  - [Chase Geigle](https://github.com/skystrife) for [cpptoml](https://github.com/skystrife/cpptoml)
  - Simon Gog, Timo Beller, Alistair Moffat, and Matthias Petri for [sdsl-lite](https://github.com/simongog/sdsl-lite)

Citing
---
If you use btllib in your research, please cite:

Nikolić et al., (2022). btllib: A C++ library with Python interface for efficient genomic sequence processing. Journal of Open Source Software, 7(79), 4720, https://doi.org/10.21105/joss.04720

If you use aaHash in your research, please cite:

Wong et al., (2023). aaHash: recursive amino acid sequence hashing. Bioinformatics Advances, vbad162, https://doi.org/10.1093/bioadv/vbad162.

