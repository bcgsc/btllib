cpp_extra_args = ['-I/usr/include/python3.7m']
link_extra_args = ['-L/usr/lib', '-lpython3.7m', '-lcrypt',
                   '-lpthread', '-ldl', '-lutil', '-lm', '-Xlinker', '-export-dynamic']

py_sources = ['btl_wrap.cxx']
libbtl_py = shared_library('_btl',
                           py_sources,
                           include_directories : inc,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true,
                           install_dir : 'py',
                           name_prefix : '')

python_files = ['btl.py']
install_data(python_files, install_dir : 'py')