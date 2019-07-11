nt_sequence1 = executable('nt_sequence1', 'nt_sequence1.cpp',
                        include_directories : inc)
nt_sequence2 = executable('nt_sequence2', 'nt_sequence2.cpp',
                        include_directories : inc)
nt_sequence3 = executable('nt_sequence3', 'nt_sequence3.cpp',
                        include_directories : inc)

test('nt_sequence1', nt_sequence1)
test('nt_sequence2', nt_sequence2, should_fail : true)
test('nt_sequence3', nt_sequence3, should_fail : true)