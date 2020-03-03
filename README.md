BTL common code library in C++ with Python and Java wrappers.

Platforms
---
- Linux
- MacOS

Documentation
---
Open `docs/index.html` with a browser and look up any classes/functions you need.

C++
---
- Dependencies
  * GCC 4.8.1+ or Clang 3.3.0+
- Copy the btllib directory in your project
- Use any header from the `btllib/include` directory

Python and Java
---
- Dependencies
  * GCC 4.8.1+ or Clang 3.3.0+
  * Python 3.5+
  * Meson and Ninja Python3 packages (optional - if they are missing, they will be automatically downloaded to a temporary directory)
- Copy the btllib directory in your project
- Run `btllib/compile`
- The wrappers correspond one-to-one with C++ code so any functions and classes can be used under the same name.
- Python
  * Use Python's `sys.path.append()` to include `btllib/python` directory
  * Include the library with `import btllib`
- Java
  * Add `btllib/java` to classpath
  * Include classes with `import btllib.*`

Contributing
---
If you want to contribute code to this repo, before making a pull request, make sure to:
- Create a build directory with `meson build` in `btllib`
- Run `ninja format` to format the whitespace in code
- Run `ninja wrap` to wrap any new C++ code for Python and Java (requires SWIG 4.0+)
- Run `ninja` to build the wrapper libraries / make sure they compile
- Run `ninja test` to run the tests
- Run `ninja tidycheck` to run clang-tidy on C++ code and make sure it passes
- Run `ninja cppcheck` to run cppcheck on C++ code and make sure it passes
- Run `ninja docs` to generate code documentation

Or simply run `ninja complete` to do all of the steps above after `meson build`