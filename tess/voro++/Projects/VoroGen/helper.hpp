#ifndef helper_h
#define helper_h

#include <fstream>
#include <sstream>
#include <cmath>

inline double rndoff(double const num, int const digits) 
{
  return std::round(num*pow(10.0,digits))/pow(10.0,digits);
}

inline bool file_exists(std::string fname) 
{
  std::ifstream f(fname.c_str());
  return f ? true : false;
}

inline void skip_lines(std::ifstream& fp, int lines_to_skip)
{
  for (int i=0;i<lines_to_skip;i++) fp.ignore(1000,'\n');
}

inline std::istringstream get_line_stream(std::ifstream& fp)
{
  std::string line;
  std::getline(fp,line);
  std::istringstream istream(line);

  return istream;
}

#endif