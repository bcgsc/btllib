project('btllib', 'cpp',
        version : '0.2.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++17' ])
add_global_arguments([ '-pthread', '-fopenmp' ], language : 'cpp')
add_global_link_arguments([ '-ldl', '-lpthread', '-fopenmp' ], language : 'cpp')

compiler = meson.get_compiler('cpp')
compiler_id = compiler.get_id()
compiler_version = compiler.version()
if compiler_id == 'gcc' and compiler_version.version_compare('<7.0.0')
  error('GCC ' + compiler_version + ' doesn\'t fully support C++17')
endif
if compiler_id == 'clang' and compiler_version.version_compare('<5.0.0')
  error('Clang ' + compiler_version + ' doesn\'t fully support C++17')
endif

btllib_include = include_directories('include')

clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target('tidycheck', command: 'scripts/tidycheck')
endif

clang_format = find_program('clang-format', required: false)
if clang_format.found()
    run_target('format', command: 'scripts/format')
endif

diff = find_program('diff', required : false)
if diff.found() and clang_format.found()
    run_target('formatcheck', command: 'scripts/formatcheck')
endif

cppcheck = find_program('cppcheck', required: false)
if cppcheck.found()
    run_target('cppcheck', command: 'scripts/cppcheck')
endif

run_target('wrap', command: 'scripts/wrap')

subdir('include')
subdir('tests')
subdir('extras')
