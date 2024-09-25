#ifndef rve_h
#define rve_h

#include <voro++.hh>
#include <vector>
#include <memory>

#include "helper.hpp"

class RVE {

  public:

    struct files 
    {
      std::string input, output;
    };

    struct parameters 
    {
      int    lloyd;
      std::vector< std::vector<double> > cuts;
      std::vector<double> sizes;
      std::vector<bool> period;
    };

  public:

    RVE(const files io_files);
    ~RVE() {}

    void gen_tesselation();
    void update_Lloyd();
    void add_cuts();

    void write_output();
    void write_hull();

  private:
    
    static const int lines_to_skip=2;
    
    static const int inline c_N[3]={6,6,6};
    static const int c_MAX_PARTICLES=10000;
    static const int c_INIT_MEM=8;

    files m_files;
    parameters m_parameters;

    int m_particles=0;
    double m_sumz=0.0;
    double m_coords[c_MAX_PARTICLES][3];
    double m_radius[c_MAX_PARTICLES];

    std::unique_ptr<voro::container_poly> p_container;
    std::vector< std::shared_ptr<voro::wall_plane> > p_planes;

  private:

    void read_particles();
    void read_parameters();

    void add_cut(std::vector<double> cut);
    void update_locations(double x_cent[][3]);
    
    void find_periodic_hull(std::vector<std::vector<double>> &all_faces,std::vector<std::vector<int>> &hull,std::vector<std::vector<double>> &hull_dist);
    void write_hull_to_file(std::vector<std::vector<double>> &all_faces,std::vector<std::vector<int>> hull, std::vector<std::vector<double>> hull_dist);
    
    std::vector<std::vector<double>> get_faces();
    void get_face_cent(int n_faces,std::vector<int> &face_ind,std::vector<double> &face_vert,std::vector<std::vector<double>> &all_faces);
    void get_cent(double x_cent[][3]);

    std::vector<std::vector<double>> unique_sort(std::vector<std::vector<double>> &vect);
    
};

#endif