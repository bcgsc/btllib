project('libbtl', 'cpp',
        version : '0.0.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++14' ])

thread_dep = dependency('threads')

inc = include_directories('include')

subdir('include')
subdir('tests')
subdir('extras')

java_subproject = subproject('java')