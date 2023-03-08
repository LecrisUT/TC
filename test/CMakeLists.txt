# TODO: Add proper unit-tests
set(old_tests "")
# Note: Similar to source, all tests should be hard-coded in order for cmake to re-configure when a test is added
list(APPEND old_tests
		example1
		example2
		example3
		example4
		example5
		example6
		example7
		example8
		example9
		example10
		example11
		example12
		example13
		example14
		example15
		example16
		example17
		example18
		example19
		example20
		example21
		example22
		example23
		example24
		example25
		)

foreach (test IN LISTS old_tests)
	file(COPY ${test} DESTINATION ${test})
	# TODO: Correct executable format
	add_test(NAME ${test}
			COMMAND $<TARGET_FILE:TCpp_tc++>
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach ()
