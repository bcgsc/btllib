header_files = ['sequence.h']

install_headers(header_files, subdir : 'btl')

tidy_checks=',-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-hicpp-no-array-decay,' + \
    '-fuchsia-default-arguments,-fuchsia-overloaded-operator,-cppcoreguidelines-pro-bounds-pointer-arithmetic'
tidy_files = files(header_files)
clang_tidy = find_program('clang-tidy', required: false)
if clang_tidy.found()
    run_target(
        'tidy',
        command: [
            clang_tidy,
            '-checks=*' + tidy_checks,
            '-p', meson.build_root()
        ] + tidy_files)
endif