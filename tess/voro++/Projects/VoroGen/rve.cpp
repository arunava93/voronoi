#include "rve.hpp"

#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

RVE::RVE(const files io_files) : m_files(io_files)
{ 
  read_parameters();

  p_container=std::unique_ptr<voro::container_poly>(new voro::container_poly(
                  -m_parameters.sizes[0]/2,m_parameters.sizes[0]/2,-m_parameters.sizes[1]/2,m_parameters.sizes[1]/2,-m_parameters.sizes[2]/2,m_parameters.sizes[2]/2,
                   c_N[0],c_N[1],c_N[2],m_parameters.period[0],m_parameters.period[1],m_parameters.period[2],c_INIT_MEM));
}

void RVE::gen_tesselation() 
{
  if (m_particles<1) read_particles();

  p_container->clear();

  for(int i=0;i<m_particles;i++) 
  {
    p_container->put(i,m_coords[i][0],m_coords[i][1],m_coords[i][2],m_radius[i]);
  }
}

void RVE::write_output() 
{
  p_container->print_custom("%i %v %C %s %t %l %w %P",m_files.output.c_str());
}

void RVE::update_Lloyd() 
{
  double x_cent[m_particles][3];

  if (m_parameters.lloyd<1001) for (int j=0;j<m_parameters.lloyd;j++) 
  {
    get_cent(x_cent);
    update_locations(x_cent);
    gen_tesselation();
  }
}

void RVE::add_cut(std::vector<double> cut) 
{
  std::shared_ptr<voro::wall_plane> plane(new voro::wall_plane(cut[0],cut[1],cut[2],cut[3]));
  p_container->add_wall(plane.get());
  p_planes.push_back(plane);
}

void RVE::add_cuts()
{
  for (std::vector<double> cut:m_parameters.cuts)
  {
    add_cut(cut);
  }
}

void RVE::write_hull() 
{
  std::vector<std::vector<int>> hull(7);
  std::vector<std::vector<double>> hull_dist(7);

  std::vector<std::vector<double>> all_faces=get_faces();

  find_periodic_hull(all_faces,hull,hull_dist);
  write_hull_to_file(all_faces,hull,hull_dist);
}

void RVE::read_particles() 
{
  std::ifstream fp(m_files.input.c_str());

  double x,y,z,r;
  int particles=0;
  int j;
  int id;

  if (fp)
  {
    skip_lines(fp,lines_to_skip);
    while (fp)
    {
      if (fp >> id >> x >> y >> z >> r) 
      {
        m_coords[particles][0]=x;
        m_coords[particles][1]=y;
        m_coords[particles][2]=z;
        m_radius[particles]=r;

        m_sumz=m_sumz+fabs(z);
        particles++;
      } 
      else if (!fp.eof()) 
      {
        voro::voro_fatal_error("File import error",VOROPP_FILE_ERROR);
      }
    }
  } else 
  {
    std::cerr << "voro++: Unable to open file " << m_files.input.c_str() << std::endl;
    exit(VOROPP_FILE_ERROR);
  }
  m_particles=particles;
}

void RVE::read_parameters() 
{
  std::ifstream fp(m_files.input.c_str());

  int    lloyd;
  double cut_z;
  double size_x,size_y,size_z;
  bool period_x,period_y,period_z;
  std::vector< std::vector<double> > cuts;

  if (fp)
  {
    std::istringstream line_stream;

    line_stream=get_line_stream(fp);
    line_stream >> lloyd >> size_x >> size_y >> size_z >> period_x >> period_y >> period_z;

    line_stream=get_line_stream(fp);
    double norm_x,norm_y,norm_z,depth;
    while (line_stream >> norm_x >> norm_y >> norm_z >> depth)
    {
      std::vector<double> cut{norm_x,norm_y,norm_z,depth};
      cuts.push_back(cut);
    }
  } else 
  {
    std::cerr << "voro++: Unable to open file " << m_files.input.c_str() << std::endl;
    exit(VOROPP_FILE_ERROR);
  }
  
  m_parameters.lloyd=lloyd;
  m_parameters.sizes=std::vector<double>{size_x,size_y,size_z};
  m_parameters.period=std::vector<bool>{period_x,period_y,period_z};
  m_parameters.cuts=cuts;
}

void RVE::get_cent(double x_cent[][3]) 
{
  voro::voronoicell c;
  voro::c_loop_all cl(*p_container.get());

  if(cl.start()) do if(p_container->compute_cell(c,cl)) 
  {
    int id=cl.pid();
    c.centroid(x_cent[id][0],x_cent[id][1],x_cent[id][2]);
  } while (cl.inc());
}

void RVE::update_locations(double x_cent[][3]) 
{
  for(int i=0;i<m_particles;i++) 
  {
    m_coords[i][0]=m_coords[i][0]+x_cent[i][0];
    m_coords[i][1]=m_coords[i][1]+x_cent[i][1];
    m_coords[i][2]=(m_sumz>1e-10) ? m_coords[i][2]+x_cent[i][2] : 0.0;
  }
}

std::vector<std::vector<double>> RVE::get_faces() 
{
  voro::voronoicell c;
  voro::c_loop_all cl(*p_container);

  std::vector<std::vector<double>> all_faces;
  
  if(cl.start()) do if(p_container->compute_cell(c,cl)) 
  {
    double x,y,z;
    
    cl.pos(x,y,z);
    std::vector<int> f;
    std::vector<double> v;
    
    c.face_vertices(f);
    c.vertices(x,y,z,v);
    get_face_cent(c.number_of_faces(),f,v,all_faces);
      
  } while (cl.inc());

  return unique_sort(all_faces);
}

void RVE::get_face_cent(int n_faces,std::vector<int> &face_ind,std::vector<double> &face_vert,std::vector<std::vector<double>> &all_faces) 
{
  int cnt=0;
  for (int i=0;i<n_faces;i++) 
  {
    int n=face_ind[cnt];
    std::vector<double> face_cent={0.0,0.0,0.0};
    for (int j=0;j<n;j++) 
    {
      cnt++;
      face_cent[0]+=face_vert[3*face_ind[cnt]];
      face_cent[1]+=face_vert[3*face_ind[cnt]+1];
      face_cent[2]+=face_vert[3*face_ind[cnt]+2];
    }
    cnt++;
    all_faces.push_back({rndoff(face_cent[0]/n,8),rndoff(face_cent[1]/n,8),rndoff(face_cent[2]/n,8)});
  }
}

void RVE::find_periodic_hull(std::vector<std::vector<double>> &all_faces,std::vector<std::vector<int>> &hull,std::vector<std::vector<double>> &hull_dist) 
{
  for (int i=0;i<all_faces.size();i++) 
  {
    std::vector<double> face_a=all_faces[i];
    for (int j=i+1;j<all_faces.size();j++) 
    {
      std::vector<double> face_b=all_faces[j];
      std::vector<int> inds(3),chks(3);

      for (int k=0;k<3;k++) 
      {
        inds[k]=std::abs(face_a[k]-face_b[k]) == m_parameters.sizes[k];
        chks[k]=inds[k] || std::abs(face_a[k]-face_b[k]) == 0.0d;
      }

      if (chks[0]+chks[1]+chks[2]==3) 
      {
        int loc=1*inds[0]+2*inds[1]+4*inds[2]-1;
        hull[loc].push_back(i);hull[loc].push_back(j);
        for (int k=0;k<3;k++) hull_dist[loc].push_back(face_a[k]-face_b[k]);
      }
    }
  }
}

void RVE::write_hull_to_file(std::vector<std::vector<double>> &all_faces,std::vector<std::vector<int>> hull, std::vector<std::vector<double>> hull_dist) 
{
  std::ofstream fout;fout.precision(8);
  fout.open(m_files.output+".hull");

  for (std::vector<double> face:all_faces) 
  {
    fout<<"["<<face[0]<<' '<<face[1]<<' '<<face[2]<<"] ";
  }
  fout<<'\n';

  for (std::vector<int> h:hull) 
  {
    if (h.size()==0) 
    {
      fout<<"[-1 -1]";
    } else 
    {
      for (int i=0;i<h.size()/2;i++) fout<<"["<<h[i*2]<<' '<<h[i*2+1]<<"]"<<' ';
    }
    fout<<'\n';
  }
  for (std::vector<double> h:hull_dist) {
    if (h.size()==0) 
    {
      fout<<"[0 0 0]";
    } else 
    {
      for (int i=0;i<h.size()/3;i++) fout<<"["<<h[i*3]<<' '<<h[i*3+1]<<' '<<h[i*3+2]<<"]"<<' ';
    }
    fout<<'\n';
  }
  fout.close();
}

std::vector<std::vector<double>> RVE::unique_sort(std::vector<std::vector<double>> &vect) 
{
  auto sort_helper=[] (const std::vector<double> &a,const std::vector<double> &b) -> bool 
  {
    if (a[0]!=b[0]) return a[0]<b[0];
    if (a[1]!=b[1]) return a[1]<b[1];
    return a[2]<b[2];
  };

  std::sort(vect.begin(),vect.end(),sort_helper);
  auto last=std::unique(vect.begin(),vect.end());
  vect.erase(last,vect.end());

  return vect;
}