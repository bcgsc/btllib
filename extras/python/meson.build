boost_dep = dependency('boost')
cpp_extra_args = ['-I/usr/include/python3.7m']
link_extra_args = ['-lboost_python37', '-L/usr/lib', '-lpython3.7m', '-lcrypt',
                   '-lpthread', '-ldl', '-lutil', '-lm', '-Xlinker', '-export-dynamic']

py_sources = ['sequence.cpp']
libbtl_py = shared_library('btl_py',
                           py_sources,
                           include_directories : inc,
                           dependencies : boost_dep,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true)