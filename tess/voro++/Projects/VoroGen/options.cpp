#include "options.hpp"
#include "helper.hpp"

#include <iostream>

using namespace options;

RVE::files options::get_opts(int argc, char* argv[])
{
  cxxopts::Options opts("A Voronoi generator","Based on Voro++ library");
  populate_opts(opts);

  auto result = opts.parse(argc, argv);
  check_result(result);

  return parameterize(result);
}

void options::populate_opts(cxxopts::Options& opts) 
{
  opts.add_options()
        ("input", "Input file name", cxxopts::value<std::string>())
        ("output", "Output file name", cxxopts::value<std::string>())
        ("h,help", "Print usage");

  opts.parse_positional({"input", "output"});
}

void options::check_result(const cxxopts::ParseResult& result) 
{
  if (result.count("input")<1 || result.count("output")<1) 
  {
    std::cout<< "voro++: Please provide input and output file names." << std::endl;
    exit(1);
  } else 
  {
    if (!file_exists(result["input"].as<std::string>())) 
    {
      std::cout<< "voro++: " << result["input"].as<std::string>() << ": No such file or directory" << std::endl;
      exit(1);
    }
  }
}

RVE::files options::parameterize(const cxxopts::ParseResult& result) 
{
  RVE::files files;

  files.input = result["input"].as<std::string>();
  files.output= result["output"].as<std::string>();

  return files;
}