project('libbtl', 'cpp',
        version : '0.0.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++17' ])

thread_dep = dependency('threads')

btl_include = include_directories('include')

clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target('tidy', command: 'scripts/tidy.sh')
endif

diff = find_program('diff', required : false)
clang_format = find_program('clang-format', required: false)
if diff.found() and clang_format.found()
    run_target('format', command: 'scripts/format.sh')
endif

cppcheck = find_program('cppcheck', required: false)
if cppcheck.found()
    run_target('cppcheck', command: 'scripts/cppcheck.sh')
endif

run_target('wrap', command: 'scripts/wrap.sh')

subdir('include')
subdir('tests')
subdir('extras')
