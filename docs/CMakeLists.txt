# TODO: Include sphinx build as well
find_package(Doxygen REQUIRED)
set(DOXYGEN_GENERATE_HTML NO)
set(DOXYGEN_GENERATE_XML YES)
doxygen_add_docs(TCpp_docs ${PROJECT_SOURCE_DIR}/src)
