#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/squared_distance_3.h>
//#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h> 
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface;
typedef CGAL::Point_3<Kernel> Point;
typedef CGAL::Segment_3<Kernel> Segment;
namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::Segment_3 Segment_3;
typedef GT::FT FT;

typedef FT (*Function) (Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;


Surface surface;
static double d;

double lattice_func(Point_3 p){
  Surface::Edge_iterator e_it = surface.edges_begin();
  //double sign= -1;
  double value = 100;
  for(; e_it!=surface.edges_end(); e_it++){
    Point p1 = e_it->vertex()->point();
    Point p2 = e_it->opposite()->vertex()->point();
    Segment seg(p1, p2);
    double dist = sqrt(CGAL::squared_distance(Point(p.x(), p.y(), p.z()), seg)) - d;
    //double tmp_sign;
    if(value > dist){      
      value = dist;
      //sign = tmp_sign;
    }
  }
  return value;
};



void edge2lattice(Surface &surface, std::vector<double> &volume, double l, std::vector<double> center){

}

int main( int argc, char** argv )
{

  d = 0.02;

  std::ifstream is(argv[1]); 
  is >> surface;
  std::cout << "mesh readed" << std::endl;
   // compute bounding box size  
  surface.
  Surface::Edge_iterator e_it = surface.edges_begin();   
  std::cout << "surface.edges_begin()" << std::endl;
  Point point = e_it->vertex()->point();
  std::cout << "e_it->vertex()->point()" << std::endl;  
  std::vector<double> min(3), max(3);
  min[0] = max[0] = point.x();
  min[1] = max[1] = point.y();
  //min[2] = max[2] = point.z();   

  for(; e_it != surface.edges_end(); e_it++){
    Point p = e_it->vertex()->point();
    if(p.x() < min[0]){
      min[0] = p.x();
    }
    else if(p.x() > max[0]){
      max[0] = p.x();
    }
    if(p.y() < min[1]){
      min[1] = p.y();
    }
    else if(p.y() > max[1]){
      max[1] = p.y();
    }
    /*
    if(p.z() < min[2]){
      min[2] = p.z();
    }
    else if(p.z() > max[2]){
      max[2] = p.z();
    } 
    */       
   }

   std::cout << min[0] << " " << min[1] << std::endl;
   std::cout << max[0] << " " << max[1] << std::endl;   

   std::vector<double> size(3), center(3);
   size[0] = max[0] - min[0];
   size[1] = max[1] - min[1];
   //size[2] = max[2] - min[2];      
   std::cout << size[0] << " " << size[1] << " " << size[2] << std::endl; 
   center[0] = (max[0] - min[0])/2;
   center[1] = (max[1] - min[1])/2;
   //center[2] = (max[2] - min[2])/2;      

   //double max_size = *std::max_element(size.begin(), size.end());
   //std::vector<double> volume;
   //edge2lattice(surface, volume, max_size, center);

   // --- polygonize implicit function ---
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Surface_3 surf(lattice_func,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 5.)); // bounding sphere
  // Note that "2." above is the *squared* radius of the bounding sphere!

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.005,  // radius bound
                                                     0.005); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surf, criteria, CGAL::Manifold_tag());

  std::ofstream ofs("lattice.off");
  //ofs << tr;
  CGAL::output_surface_facets_to_off(ofs, c2t3);
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

   return 0;
}
