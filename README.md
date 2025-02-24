# Anonymous Communication

## Non-docker instructions:

Build with `make`.  The build has been tested on Ubuntu 20.04 and Ubuntu
22.04.  You'll need libbsd-dev and libboost-all-dev.

- `./prac -o 0 annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
- `./prac -o 1 localhost annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
- `./prac -o 2 localhost localhost annoncomm -m 10 -d 10 -i 5 -e 1 -opt 0 -s 0`
