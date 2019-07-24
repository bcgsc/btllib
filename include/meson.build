header_files = ['sequence.h']

install_headers(header_files, subdir : 'btl')

tidy_checks=',-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-hicpp-no-array-decay,' + \
    '-fuchsia-default-arguments,-fuchsia-overloaded-operator,-cppcoreguidelines-pro-bounds-pointer-arithmetic,' + \
    '-hicpp-avoid-c-arrays,-google-readability-casting,-cppcoreguidelines-pro-bounds-constant-array-index,' + \
    '-modernize-avoid-c-arrays,-cppcoreguidelines-avoid-c-arrays,-modernize-use-equals-default,-hicpp-use-equals-default,' + \
    '-google-runtime-references,-google-explicit-constructor,-hicpp-explicit-conversions,' + \
    '-misc-non-private-member-variables-in-classes,-cppcoreguidelines-special-member-functions,-modernize-use-nodiscard,' + \
    '-hicpp-special-member-functions'
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