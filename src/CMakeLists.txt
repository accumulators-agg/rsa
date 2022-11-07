add_executable(
  main

  main.cpp
)

target_link_libraries(
  main

  gmp         # gmp
  gmpxx       # gmp
  crypto      # openssl
  fmt # https://github.com/fmtlib/fmt
)
