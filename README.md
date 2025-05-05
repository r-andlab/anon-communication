# Anonymous Communication

## Build

Build Dependencies:

- C++ Compiler (C++17 Standard)
- libBSD (`libbsd-dev` on Ubuntu)
- Boost (`libboost-all-dev` on Ubuntu)

Build using CMake:

1. `cmake -S . -B build`
2. `cmake --build build`

Inside `build` will be `anoncomm_client` and `anoncomm_server`.

## Manual testing instructions:

- `anoncomm_server -o 0 annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
- `anoncomm_server -o 1 localhost annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
- `anoncomm_server -o 2 localhost localhost annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
