set(EXAMPLES_PDES
    dirichlet
    acoustic_first_order_pde
    acoustic_1D_pde
)

#### BOOST LINKING ####
find_package(Boost 1.65.1 REQUIRED COMPONENTS
             date_time filesystem iostreams)

foreach(example ${EXAMPLES_PDES})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
    target_link_libraries(${example} Boost::date_time Boost::filesystem
                          Boost::iostreams)
endforeach()

add_custom_target(examples-pdes)
add_dependencies(examples-pdes ${EXAMPLES_PDES})
add_dependencies(examples ${EXAMPLES_PDES})
