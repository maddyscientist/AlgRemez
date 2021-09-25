add_executable(test main.C)
target_link_libraries(test algremez ${GMP_LIB} ${MPFR_LIBRARIES})
target_include_directories(test PUBLIC ${PROJECT_SOURCE_DIR}/include)
