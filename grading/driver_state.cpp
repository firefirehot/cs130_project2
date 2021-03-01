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
	pixel * alloc_color = new pixel[width*height];
	float * alloc_depth = new float[width*height];
	for(int x = 0; x< width*height; x++)
		alloc_depth[x] = 1.1;
	for(int x = 0; x< width*height; x++)
		alloc_color[x] = make_pixel(0,0,0);
    state.image_width=width;
    state.image_height=height;
    state.image_color=alloc_color;
    state.image_depth=alloc_depth;
   // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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
	
    //std::cout<<"TODO: implement rendering."<<std::endl;
    data_geometry * alloc_geometry = new data_geometry[state.num_vertices];
    data_vertex input_for_vertex_shader;
    
    for(int x = 0; x < state.num_vertices * state.floats_per_vertex; x = x + state.floats_per_vertex){
		alloc_geometry[x/state.floats_per_vertex].data = new float[MAX_FLOATS_PER_VERTEX];
		input_for_vertex_shader.data = &state.vertex_data[x];
		
		for(int f = 0; f < state.floats_per_vertex;f++)//I don't know if this has to be done but it seems right
			alloc_geometry[x/state.floats_per_vertex].data[f] = state.vertex_data[x+f];//same here
		
		state.vertex_shader(input_for_vertex_shader, alloc_geometry[x/state.floats_per_vertex],state.uniform_data);
		//std::cout << alloc_geometry[x/state.floats_per_vertex].gl_Position[0] <<  " __" << alloc_geometry[x/state.floats_per_vertex].gl_Position[1] << std::endl;
		//std::cout << input_for_vertex_shader.data[0] <<  " __" << input_for_vertex_shader.data[1] << std::endl;
		//std::cout << alloc_geometry[x/state.floats_per_vertex].data[0] <<  " __" << alloc_geometry[x/state.floats_per_vertex].data[1] << std::endl;
}
	switch(type){
	case render_type::indexed:
	break;
	case render_type::triangle:
	for(int i = 0; i < state.num_vertices; i = i+3)
		rasterize_triangle(state,alloc_geometry[i],alloc_geometry[i+1],alloc_geometry[i+2]);
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
	/*clipping goes before divide by w*/
	
	
vec4 v0_Position = v0.gl_Position/v0.gl_Position[3];
if(v0.gl_Position[3] == 0)
	v0_Position = v0.gl_Position;
vec4 v1_Position = v1.gl_Position/v1.gl_Position[3];
if(v1.gl_Position[3] == 0)
	v1_Position = v1.gl_Position;
vec4 v2_Position = v2.gl_Position/v2.gl_Position[3];
if(v2.gl_Position[3] == 0)
	v2_Position = v2.gl_Position;


//.75*100/2 + 100/2
//NDC*Width/2 + width/2 + 1/2

double v0_NDC_x = v0_Position[0]*(state.image_width)/2 + (state.image_width-1)/2;
double v0_NDC_y = v0_Position[1]*(state.image_height)/2 + (state.image_height-1)/2;
double v1_NDC_x = v1_Position[0]*(state.image_width)/2 + (state.image_width-1)/2;
double v1_NDC_y = v1_Position[1]*(state.image_height)/2 + (state.image_height-1)/2;
double v2_NDC_x = v2_Position[0]*(state.image_width)/2 + (state.image_width-1)/2;
double v2_NDC_y = v2_Position[1]*(state.image_height)/2 + (state.image_height-1)/2;


double total_area = area_calc(v0_NDC_x, v0_NDC_y, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y);
double alpha;
double beta;
double gamma;

data_fragment input;
data_output out;
float array[MAX_FLOATS_PER_VERTEX];
input.data = array;
/*
for(int i = 0; i < state.floats_per_vertex; i++){
	if(state.interp_rules[i] == interp_type::flat){
		input.data[i] = v0.data[i];
	//	std:: cout << "-- " << v0.data[i];
	}
	else if (state.interp_rules[i] == interp_type::smooth){
		//input.data[i] = ;
		
	}
	else if (state.interp_rules[i] == interp_type::noperspective ){
		input.data[i] = ;
		
	}
	else
		std::cout << "error in data fragment calculation,driverstate.cpp/rasterize_triangle" << std::endl;
}
 */

/*std::cout << std::endl << v0_NDC_x
<< std::endl << v0_NDC_y
<< std::endl << v1_NDC_x
 << std::endl << v1_NDC_y
<< std::endl << v2_NDC_x
<< std::endl <<v2_NDC_y ;*/


for(int y = 0; y < state.image_height; y++)
	for(int x = 0; x < state.image_width;x++){
		alpha = area_calc(x, y, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y)/total_area;
		//PBC
		beta = area_calc(x, y, v2_NDC_x, v2_NDC_y, v0_NDC_x, v0_NDC_y)/total_area;
		//PCA
		gamma = area_calc(x, y, v0_NDC_x, v0_NDC_y, v1_NDC_x, v1_NDC_y)/total_area;
		//PAB
		
		for(int i = 0; i < state.floats_per_vertex; i++){
			if(state.interp_rules[i] == interp_type::flat){
				input.data[i] = v0.data[i];
			//	std:: cout << "-- " << v0.data[i];
			}
			else if (state.interp_rules[i] == interp_type::smooth){
				//input.data[i] = ;
			
			}
			else if (state.interp_rules[i] == interp_type::noperspective ){
				input.data[i] = alpha*v0.data[i] + beta*v1.data[i] + gamma*v2.data[i];
			}	
			else
				std::cout << "error in data fragment calculation,driverstate.cpp/rasterize_triangle" << std::endl;
		}
		
		
		if(alpha >= 0 && beta >= 0 && gamma >= 0){
			
		state.fragment_shader(input, out, state.uniform_data);
		if(state.image_depth[x + y*(state.image_width)] > alpha*v0_Position[2] + beta*v1_Position[2] + gamma*v2_Position[2]){
			state.image_depth[x + y*(state.image_width)] = alpha*v0_Position[2] + beta*v1_Position[2] + gamma*v2_Position[2];
			state.image_color[x + y*(state.image_width)] = make_pixel(out.output_color[0]*255,out.output_color[1]*255,out.output_color[2]*255);
		}
		}
	}

   // std::cout<<"TODO: implement rasterization"<<std::endl;
}

