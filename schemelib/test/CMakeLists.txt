
# shitty "data" files to compile faster
list( APPEND LIB_CC ../scheme/rosetta/score/EtableParams_init.lib.cc ../scheme/rosetta/atype_map.lib.cc ../scheme/chemical/ligand_factory.lib.cc )
SET_SOURCE_FILES_PROPERTIES( ${LIB_CC} PROPERTIES COMPILE_FLAGS -O0 )
add_library( scheme ${LIB_CC} )

FILE(GLOB L1 "../scheme/[0-9a-zA-Z_]*.cc")
FILE(GLOB L2 "../scheme/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L3 "../scheme/*/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L4 "../scheme/*/*/*/[0-9a-zA-Z_]*.cc")
FILE(GLOB L5 "../scheme/*/*/*/*/[0-9a-zA-Z_]*.cc")

list(APPEND EXTRA_LIBS boost_system  )

include_directories(".")

add_executable( test_libscheme main_test.cc ${L1} ${L2} ${L3} ${L4} ${L5} )
target_link_libraries( test_libscheme scheme ${EXTRA_LIBS} )
# install ( TARGETS test_libscheme RUNTIME DESTINATION bin )

add_executable(quick_test_libscheme quick_test.cc  )
target_link_libraries(quick_test_libscheme scheme ${EXTRA_LIBS})
