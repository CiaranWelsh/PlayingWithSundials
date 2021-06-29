

set(
        SUNDIALS_LIBRARIES
        sundials_cvodes
)

set(
        SUNDIALS_INCLUDE_DIRS
        "/Users/ciaranwelsh/Documents/sundials/include"

)

set(tgt CWSimpleFluxFSAWithUserSuppliedFunctions)
add_executable(${tgt} SimpleFluxFSAWithUserSuppliedFunctions.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWSimpleFluxFSAWithFiniteDifferenceApprox)
add_executable(${tgt} SimpleFluxFSAWithFiniteDifferenceApprox.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWcvsRoberts_FSA_dns)
add_executable(${tgt} cvsRoberts_FSA_dns_cw.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWReplaceEqnInCvsRoberts_FSA_dns_cw)
add_executable(${tgt} ReplaceEqnInCvsRoberts_FSA_dns_cw.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWVenkatramanFSAWithFiniteDifferenceApprox)
add_executable(${tgt} VenkatramanFSAWithFiniteDifferenceApprox.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWVenkatramanAdjointWithFiniteDifferenceApprox)
add_executable(${tgt} VenkatramanAdjointWithFiniteDifferenceApprox.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})

set(tgt CWSimpleFluxFSAWithFDATake2)
add_executable(${tgt} SimpleFluxFSAWithFDATake2.c)
target_link_libraries(${tgt} PRIVATE ${SUNDIALS_LIBRARIES})
target_include_directories(${tgt} PRIVATE ${SUNDIALS_INCLUDE_DIRS})














