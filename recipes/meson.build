recipe_files = run_command('../scripts/get-files', 'recipes').stdout().strip().split()

executable('indexlr',
            meson.project_source_root() + '/recipes/indexlr.cpp',
            include_directories : btllib_include,
            dependencies : deps + [ btllib_dep ],
            install : true,
            install_dir : 'bin')

executable('seqgen',
            meson.project_source_root() + '/recipes/seqgen.cpp',
            include_directories : btllib_include,
            dependencies : deps + [ btllib_dep, argparse_dep ],
            install : true,
            install_dir : 'bin',
            override_options : ['cpp_std=c++17'])

executable('mi_bf_generate',
            meson.project_source_root() + '/recipes/mi_bloom_filter.cpp',
            include_directories : btllib_include,
            dependencies : deps + [ btllib_dep ],
            install : true,
            install_dir : 'bin')
