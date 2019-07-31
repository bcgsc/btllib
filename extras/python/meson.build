cpp_extra_args = run_command('python-config', '--includes').stdout().strip().split()
link_extra_args = run_command('python-config', '--libs').stdout().strip().split()

cpp_extra_args += ['-w']

py_sources = ['btl_wrap.cxx']
libbtl_py = shared_library('_btl',
                           py_sources,
                           include_directories : btl_include,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true,
                           install_dir : 'py',
                           name_prefix : '')

python_files = ['btl.py']
install_data(python_files, install_dir : 'py')