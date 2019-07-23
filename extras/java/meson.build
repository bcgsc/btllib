cpp_extra_args = ['-I/usr/lib/jvm/default/include', '-I/usr/lib/jvm/default/include/linux']

java_sources = ['btl_wrap.cxx']
libbtl_java = shared_library('btl',
                             java_sources,
                             include_directories : inc,
                             cpp_args : cpp_extra_args,
                             install : true,
                             install_dir : 'java')

#install_data('btl.py', install_dir : 'py')