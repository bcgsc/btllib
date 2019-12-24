cpp_extra_args = run_command('python3-config', '--includes').stdout().strip().split()
link_extra_args = run_command('python3-config', '--libs').stdout().strip().split()

libpython = run_command('../../scripts/get_libpython').stdout().strip()
link_extra_args += [ libpython ]

cpp_extra_args += ['-w']

py_sources = ['btllib_wrap.cxx']
btllib_py = shared_library('_btllib',
                           py_sources,
                           include_directories : btllib_include,
                           cpp_args : cpp_extra_args,
                           link_args : link_extra_args,
                           install : true,
                           install_dir : 'python',
                           name_prefix : '')

python_files = ['btllib.py']
install_data(python_files, install_dir : 'python')