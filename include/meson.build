header_files = ['sequence.h']

install_headers(header_files, subdir : 'btl')

tidy_files = files(header_files)
clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target('tidy', command: ['scripts/tidy.sh'] + tidy_files)
endif

format_files = files(header_files)
cat = find_program('cat', required : false)
diff = find_program('diff', required : false)
clang_format = find_program('clang-format', required: false)
if cat.found() and diff.found() and clang_format.found()
    run_target('format', command: ['scripts/format.sh'] + format_files)
endif