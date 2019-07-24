project('libbtl', 'cpp',
        version : '0.0.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++17' ])

thread_dep = dependency('threads')

btl_include = include_directories('include')

tidy_files = files(['include/sequence.h'])
clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target('tidy', command: ['scripts/tidy.sh'] + tidy_files)
endif

format_files = files([])
cat = find_program('cat', required : false)
diff = find_program('diff', required : false)
clang_format = find_program('clang-format', required: false)
if cat.found() and diff.found() and clang_format.found()
    run_target('format', command: ['scripts/format.sh'] + format_files)
endif

subdir('include')
subdir('tests')
subdir('extras')