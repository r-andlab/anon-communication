#include <bsd/stdlib.h> // arc4random_buf

#include "online.hpp"

#include "anoncomm.hpp"
#include "shapes.hpp"

void online_main(MPCIO &mpcio, const PRACOptions &opts, char **args) {
  if (!*args) {
    std::cerr
        << "Mode is required as the first argument when not preprocessing.\n";
    return;
  } else if (!strcmp(*args, "annoncomm")) {
    ++args;
    AnnonComm(mpcio, opts, args);
  } else if (!strcmp(*args, "test")) {
    ++args;
    TestMode(mpcio, opts, args);
  } else {
    std::cerr << "Unknown mode " << *args << "\n";
  }
}
