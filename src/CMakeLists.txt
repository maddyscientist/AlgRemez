add_library(algremez STATIC alg_remez.C)
target_include_directories(algremez PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

if(ENABLE_OPENMP)
  target_link_libraries(algremez PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS algremez)
