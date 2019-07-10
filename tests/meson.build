sequence1 = executable('sequence1', 'sequence1.cpp',
                        include_directories : inc)
sequence2 = executable('sequence2', 'sequence2.cpp',
                        include_directories : inc)
sequence3 = executable('sequence3', 'sequence3.cpp',
                        include_directories : inc)

test('sequence1', sequence1)
test('sequence2', sequence2, should_fail : true)
test('sequence3', sequence3, should_fail : true)