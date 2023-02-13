recipe_files = run_command('../scripts/get-files', 'recipes').stdout().strip().split()

foreach recipe_file : recipe_files
  executable(fs.stem(recipe_file),
             recipe_file,
             include_directories : btllib_include,
             dependencies : deps + [ btllib_dep, argparse_dep ],
             install : true,
             install_dir : 'bin',
             override_options : ['cpp_std=c++17'])
endforeach