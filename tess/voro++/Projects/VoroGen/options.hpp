#ifndef options_h
#define options_h

#include <cxxopts.hpp>
#include "rve.hpp"

namespace options {
    
  RVE::files get_opts(int argc, char* argv[]);

  void populate_opts(cxxopts::Options& opts);

  void check_result(const cxxopts::ParseResult& result);

  RVE::files parameterize(const cxxopts::ParseResult& result); 

}

#endif
