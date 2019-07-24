#!/bin/bash

cd "${MESON_BUILD_ROOT}"

tidy_checks="*,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-hicpp-no-array-decay, \
    -fuchsia-default-arguments,-fuchsia-overloaded-operator,-cppcoreguidelines-pro-bounds-pointer-arithmetic, \
    -hicpp-avoid-c-arrays,-google-readability-casting,-cppcoreguidelines-pro-bounds-constant-array-index, \
    -modernize-avoid-c-arrays,-cppcoreguidelines-avoid-c-arrays,-modernize-use-equals-default,-hicpp-use-equals-default, \
    -google-runtime-references,-google-explicit-constructor,-hicpp-explicit-conversions, \
    -misc-non-private-member-variables-in-classes,-cppcoreguidelines-special-member-functions,-modernize-use-nodiscard, \
    -hicpp-special-member-functions"

clang-tidy -checks="$tidy_checks" $@