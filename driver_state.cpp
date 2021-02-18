#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
	pixel * alloc_color = new pixel[width*height];//example 8
	//pixel temp_pix = make_pixel(0,0,0);
	//alloc_color = &temp_pix;
	for(int x = 0; x < width*height; x++){ //example 0-7
		alloc_color[x] = make_pixel(0,0,0);
	}
    state.image_width=width;
    state.image_height=height;
    state.image_color=alloc_color;
    state.image_depth=0;
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
   // std::cout<<"TODO: implement rendering."<<std::endl;
    data_geometry * alloc_geometry = new data_geometry[state.num_vertices];
    data_vertex input_for_vertex_shader;
     
	for(int x = 0; x < state.num_vertices*state.floats_per_vertex;x = x + state.floats_per_vertex){
		alloc_geometry[x/state.floats_per_vertex].data = new float[MAX_FLOATS_PER_VERTEX];  //allocate data_geometry data memeber
		input_for_vertex_shader.data = &state.vertex_data[x]; //here we make the data_vertex point to the relevant vertex data
		state.vertex_shader(input_for_vertex_shader, alloc_geometry[x/state.floats_per_vertex],state.vertex_data); //this fills in the data	
		std::cout << alloc_geometry[x/state.floats_per_vertex].gl_Position[0] << "--" << alloc_geometry[x/state.floats_per_vertex].gl_Position[1] << std::endl << x << std:: endl;
		
	}
    
    switch(type){
	case  render_type::indexed:
	break;
	case render_type::triangle:
	for (int i = 0; i < state.num_vertices; i = i + 3)
		rasterize_triangle(state, alloc_geometry[i],alloc_geometry[i+1], alloc_geometry[i+2]);
	break;
	case render_type::fan:
	break;
	case render_type::strip:
	break;
	case render_type::invalid:
	break;

	}
    
}
 

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.

double area_calc(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y){
	double d = 0.5*((p1_x*p2_y - p2_x*p1_y) + (p2_x*p0_y - p0_x*p2_y) + (p0_x*p1_y - p1_x*p0_y));
	return d;
}

void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
	
	vec4 v0_Position = v0.gl_Position/v0.gl_Position[3];
	if(v0.gl_Position[3] == 0)
		v0_Position = v0.gl_Position;
	vec4 v1_Position = v1.gl_Position/v1.gl_Position[3];
	if(v1.gl_Position[3] == 0)
		v1_Position = v1.gl_Position;
	vec4 v2_Position = v2.gl_Position/v2.gl_Position[3];
	if(v2.gl_Position[3] == 0)
		v2_Position = v2.gl_Position;
	
	
	double v0_NDC_x = v0_Position[0]*(state.image_width-1)/2 + (state.image_width-1)/2 + 1/2;
	double v0_NDC_y = v0_Position[1]*(state.image_height-1)/2 + (state.image_height-1)/2 + 1/2;
	double v1_NDC_x = v1_Position[0]*(state.image_width-1)/2 + (state.image_width-1)/2 + 1/2;
	double v1_NDC_y = v1_Position[1]*(state.image_height-1)/2 + (state.image_height-1)/2 + 1/2;
	double v2_NDC_x = v2_Position[0]*(state.image_width-1)/2 + (state.image_width-1)/2 + 1/2;
	double v2_NDC_y = v2_Position[1]*(state.image_height-1)/2 + (state.image_height-1)/2 + 1/2;
	
	double total_area = area_calc(v0_NDC_x, v0_NDC_y, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y);
	double alpha;
	double beta;
	double gamma;
	
	std::cout << std::endl << v0_NDC_x
	<< std::endl << v0_NDC_y 
	<< std::endl << v1_NDC_x 
 << std::endl << v1_NDC_y 
	<< std::endl << v2_NDC_x 
	<< std::endl <<v2_NDC_y ;
	
	
	for(int y = 0; y < state.image_height; y++)
		for(int x = 0; x < state.image_width;x++){
			alpha = area_calc(x+0.5, y+0.5, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y)/total_area;
			beta = area_calc(x+0.5, y+0.5, v2_NDC_x, v2_NDC_y, v0_NDC_x, v0_NDC_y)/total_area;
			gamma = area_calc(x+0.5, y+0.5, v0_NDC_x, v0_NDC_y, v1_NDC_x, v1_NDC_y)/total_area;
			if(alpha > 0 && beta > 0 && gamma > 0){
				
				state.image_color[x + y*(state.image_width)] = make_pixel(255,255,255);
			}
			
		}
	
    std::cout<<"TODO: implement rasterization"<<std::endl;
}




