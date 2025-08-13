enable_language(Fortran)
file(GLOB tRecX_zslapack_SRC "./zlapack/*.f")
add_library(tRecX_zslapack SHARED ${tRecX_zslapack_SRC})
#SET_TARGET_PROPERTIES(tRecX_zslapack PROPERTIES LINKER_LANGUAGE Fortran)
