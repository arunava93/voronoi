#include "rve.hpp"
#include "options.hpp"

int main(int argc, char* argv[]) 
{
  RVE::files io_files=options::get_opts(argc,argv);
  
  RVE rve{io_files};

  rve.gen_tesselation();
  rve.update_Lloyd();

  rve.add_cuts();

  rve.write_output();
  rve.write_hull();

  return 0;
}