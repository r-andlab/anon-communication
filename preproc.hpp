#ifndef __PREPROC_HPP__
#define __PREPROC_HPP__

#include "mpcio.hpp"
#include "options.hpp"

void preprocessing_comp(MPCIO &mpcio, const PRACOptions &opts, char **args);
void preprocessing_server(MPCServerIO &mpcio, const PRACOptions &opts, char **args);

#endif
