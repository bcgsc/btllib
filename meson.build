project('btllib', 'cpp',
        version : '1.4.0',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++11', 'warning_level=3', 'werror=true' ])

# Configuration
# ===========================================================

conf = configuration_data()
conf.set('PROJECT_VERSION', meson.project_version(), description : 'btllib version. Useful to append to the help prompt of recipes.')
configure_file(input : 'include/btllib/config.hpp.in', output : 'config.hpp', configuration : conf)

# Modules
# ===========================================================

fs = import('fs')
cmake = import('cmake')

# Compiler and global flags
# ===========================================================

compiler = meson.get_compiler('cpp')
compiler_id = compiler.get_id()
compiler_version = compiler.version()

if compiler_id == 'gcc' and compiler_version.version_compare('<6.0.0')
  error('GCC ' + compiler_version + ' is unsupported.')
endif
if compiler_id == 'clang' and compiler_version.version_compare('<5.0.0')
  error('Clang ' + compiler_version + ' is unsupported.')
endif

global_args = []
global_link_args = [ '-ldl' ]

if compiler_id == 'clang'
  global_link_args += [ '-lstdc++', '-lm' ]
endif

add_global_arguments(global_args, language : 'cpp')
add_global_link_arguments(global_link_args, language : 'cpp')

# Dependencies and subprojects
# ===========================================================

threads_dep = dependency('threads')
openmp_dep = dependency('openmp', required : false)

cpptoml_subproject = cmake.subproject('cpptoml')
cpptomp_dep = cpptoml_subproject.dependency('cpptoml')

run_command(meson.source_root() + '/scripts/prepare-sdsl')

cmake_options = cmake.subproject_options()
cmake_options.set_override_option('werror', 'false')

sdsl_subproject = cmake.subproject('sdsl-lite', options : cmake_options)
sdsl_dep = sdsl_subproject.dependency('sdsl')

deps = [ threads_dep, openmp_dep, cpptomp_dep, sdsl_dep ]

# Custom run commands
# ===========================================================

rootpath = meson.source_root()

swig = find_program('swig', required : false)
if swig.found()
  run_target('wrap', command: join_paths(rootpath, 'scripts/wrap'))
endif

cppcheck = find_program('cppcheck', required : false)
if cppcheck.found()
  run_target('cppcheck', command: join_paths(rootpath, 'scripts/cppcheck'))
endif

run_target('sanitize-address', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'address'
])
run_target('sanitize-memory', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'memory'
])
run_target('sanitize-thread', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'thread'
])
run_target('sanitize-undefined', command: [
  join_paths(rootpath, 'scripts/sanitize'),
  'undefined'
])

doxygen = find_program('doxygen', required : false)
if doxygen.found()
  run_target('docs', command: join_paths(rootpath, 'scripts/docs'))
endif

run_target('quality-assurance', command: join_paths(rootpath, 'scripts/quality-assurance'))

run_target('test-wrappers', command: join_paths(rootpath, 'scripts/test-wrappers'))

# The library
# ===========================================================

btllib_include = [ include_directories('include'), include_directories('.') ]
btllib_sources = run_command('scripts/get-files', 'src').stdout().strip().split()

btllib_lib = static_library('btllib',
                           btllib_sources,
                           include_directories : btllib_include,
                           dependencies : deps,
                           install : true,
                           install_dir : 'lib')

btllib_dep = declare_dependency(
  link_with : btllib_lib,
  include_directories : btllib_include
)

subdir('wrappers')
subdir('recipes')
if get_option('buildtype') != 'release'
  subdir('tests')
endif
