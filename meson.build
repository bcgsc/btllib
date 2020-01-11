project('btllib', 'cpp',
        version : '0.2.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++17' ])
add_global_arguments([ '-pthread', '-fopenmp' ], language : 'cpp')
add_global_link_arguments([ '-ldl', '-lpthread', '-fopenmp', '-lboost_iostreams' ], language : 'cpp')

thread_dep = dependency('threads')

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
