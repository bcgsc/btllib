cpp_extra_args = ['-I/usr/lib/jvm/java-10-openjdk/include', '-I/usr/lib/jvm/java-10-openjdk/include/linux']

java_sources = ['hello.cpp']
libbtl_java = shared_library('btl_java',
                             java_sources,
                             include_directories : inc,
                             cpp_args : cpp_extra_args,
                             install : true)