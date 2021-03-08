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
    data_geometry * alloc_geometry = new data_geometry[state.num_vertices];//here I create a array for the vertexes
    data_vertex input_for_vertex_shader;
    
    for(int x = 0; x < state.num_vertices * state.floats_per_vertex; x = x + state.floats_per_vertex){
		alloc_geometry[x/state.floats_per_vertex].data = new float[MAX_FLOATS_PER_VERTEX]; //I give each vertex's data MAX_FLOATS... number of spaces
		input_for_vertex_shader.data = &state.vertex_data[x];
		
		for(int f = 0; f < state.floats_per_vertex;f++)
			alloc_geometry[x/state.floats_per_vertex].data[f] = state.vertex_data[x+f];
		
		state.vertex_shader(input_for_vertex_shader, alloc_geometry[x/state.floats_per_vertex],state.uniform_data);
}
	switch(type){
	case render_type::indexed:
	break;
	case render_type::triangle:
	for(int i = 0; i < state.num_vertices; i = i+3){
		clip_triangle(state,alloc_geometry[i],alloc_geometry[i+1],alloc_geometry[i+2],0);
		//rasterize_triangle(state,alloc_geometry[i],alloc_geometry[i+1],alloc_geometry[i+2]);
	}
	break;
	case render_type::fan:
	if(state.num_vertices > 2)
		clip_triangle(state,alloc_geometry[0],alloc_geometry[1],alloc_geometry[2],0);
	for(int i = 2; i < state.num_vertices-1; i++)
		clip_triangle(state,alloc_geometry[0],alloc_geometry[i],alloc_geometry[i+1],0);
	break;
	case render_type::strip:
	for(int i = 0; i < state.num_vertices-2; i++){
		if(i%2 == 0)
			clip_triangle(state,alloc_geometry[i],alloc_geometry[i+1],alloc_geometry[i+2],0);
		else
			clip_triangle(state,alloc_geometry[i],alloc_geometry[i+2],alloc_geometry[i+1],0);
	}
	break;
	case render_type::invalid:
	break;
	}
	for(int x = 0; x < state.num_vertices; x++){
		delete [] alloc_geometry[x].data;	
	}
	delete [] alloc_geometry;
	//dealocate at end of render
}


data_geometry find_point(const data_geometry& v0, const data_geometry& v1, int face,driver_state& state){
	data_geometry * return_geo = new data_geometry;
	return_geo->data = new float[MAX_FLOATS_PER_VERTEX];
	float alpha;
	switch(face){
		
	case 0://right
	alpha = (v1.gl_Position[3] - v1.gl_Position[0])/(v1.gl_Position[3] - v1.gl_Position[0] + v0.gl_Position[0] - v0.gl_Position[3]);
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;
		
	case 1://left
		alpha = (v1.gl_Position[3] + v1.gl_Position[0])/(v1.gl_Position[3] + v1.gl_Position[0] - v0.gl_Position[3] - v0.gl_Position[0]);
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;	
		
	case 2://up
		alpha = (v1.gl_Position[3] - v1.gl_Position[1])/(v1.gl_Position[3] - v1.gl_Position[1] + v0.gl_Position[1] - v0.gl_Position[3]);
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;	
		
	case 3://down
			alpha = (v1.gl_Position[3] + v1.gl_Position[1])/(v1.gl_Position[3] + v1.gl_Position[1] - v0.gl_Position[1] - v0.gl_Position[3]);
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;	
		
	case 4://far
			alpha = (v1.gl_Position[3] - v1.gl_Position[2])/(v1.gl_Position[3] - v1.gl_Position[2] + v0.gl_Position[2] - v0.gl_Position[3]);
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;	
		
	case 5://near
			alpha = (v1.gl_Position[3] + v1.gl_Position[2])/(v1.gl_Position[3] + v1.gl_Position[2] - v0.gl_Position[2] - v0.gl_Position[3] );
	return_geo->gl_Position = alpha*v0.gl_Position + (1-alpha)*v1.gl_Position;
	break;		
	default:
	std::cout << "-error-";
	break;
		
	}
	for (int i = 0; i < state.floats_per_vertex; i ++){
	return_geo->data[i] = alpha*v0.data[i] + (1-alpha)*v1.data[i];	
	}
	return *return_geo;
	
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.

void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
	float vert_0 = v0.gl_Position[3];
	//float vert_0 = 1;
	float vert_1 = v1.gl_Position[3];
	//float vert_1 = 1;
	float vert_2 = v2.gl_Position[3];
	//float vert_2 = 1;
	int ind = face / 2;
	int neg = 1;
	if(face % 2 == 1){
		neg = -1;
		//ind++;
	}
	
	if(face > 5){
		rasterize_triangle(state, v0,v1,v2);
		return;
	}
	
	if(neg * v0.gl_Position[ind] > vert_0 || neg * v1.gl_Position[ind] > vert_1 || neg * v2.gl_Position[ind] > vert_2){
		//v0 inside
		if(neg * v0.gl_Position[ind] < vert_0 && neg * v1.gl_Position[ind] > vert_1 && neg * v2.gl_Position[ind] > vert_2)
		{
			clip_triangle(state, v0, find_point(v0,v1,face,state), find_point(v0,v2,face,state),face + 1);
		}
		//v1 inside
		else if(neg * v0.gl_Position[ind] > vert_0 && neg * v1.gl_Position[ind] < vert_1 && neg * v2.gl_Position[ind] > vert_2)
		{
			clip_triangle(state, v1, find_point(v1,v2,face,state), find_point(v1,v0,face,state),face + 1);
		}
		//v2 inside
		else if(neg * v0.gl_Position[ind] > vert_0 && neg * v1.gl_Position[ind] > vert_1 && neg * v2.gl_Position[ind] < vert_2)
		{
			clip_triangle(state, v2, find_point(v2,v0,face,state), find_point(v2,v1,face,state),face + 1);
		}
		//v0 and v1 inside
		else if(neg * v0.gl_Position[ind] < vert_0 && neg * v1.gl_Position[ind] < vert_1 && neg * v2.gl_Position[ind] > vert_2)
		{
			clip_triangle(state, v0, v1,find_point(v0,v2,face,state),face + 1);
			clip_triangle(state, v1, find_point(v1,v2,face,state), find_point(v0,v2,face,state),face + 1);
		}
		//v0 and v2 inside
		else if(neg * v0.gl_Position[ind] < vert_0 && neg * v1.gl_Position[ind] > vert_1 && neg * v2.gl_Position[ind] < vert_2)
		{
			clip_triangle(state, v2, v0,find_point(v2,v1,face,state),face + 1);
			clip_triangle(state, v0, find_point(v0,v1,face,state), find_point(v2,v1,face,state),face + 1);
		}
		//v1 and v2 inside
		else if(neg * v0.gl_Position[ind] > vert_0 && neg * v1.gl_Position[ind] < vert_1 && neg * v2.gl_Position[ind] < vert_2)
		{
			clip_triangle(state, v1, v2,find_point(v1,v0,face,state),face + 1);
			clip_triangle(state, v2, find_point(v2,v0,face,state), find_point(v1,v0,face,state),face + 1);
		}
		
	}
	else{
		clip_triangle(state, v0, v1, v2,face + 1);
	}
	
	
	
	return;
}





// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
double area_calc(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y){
double d = 0.5*((p1_y - p2_y)*p0_x + (p2_y - p0_y)*p1_x + (p0_y - p1_y)*p2_x);
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
double alphaP;
double betaP;
double gammaP;

data_fragment input;
data_output out;
float array[MAX_FLOATS_PER_VERTEX];
input.data = array;

/*std::cout << std::endl << v0_NDC_x
<< std::endl << v0_NDC_y
<< std::endl << v1_NDC_x
 << std::endl << v1_NDC_y
<< std::endl << v2_NDC_x
<< std::endl <<v2_NDC_y ;*/

float wp;
for(int y = 0; y < state.image_height; y++)
	for(int x = 0; x < state.image_width;x++){
		alpha = area_calc(x, y, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y)/total_area;
		//alpha = area_calc(x, y, v1_NDC_x, v1_NDC_y, v2_NDC_x, v2_NDC_y)/total_area;
		//PBC
		beta = area_calc(v0_NDC_x, v0_NDC_y, x, y, v2_NDC_x, v2_NDC_y)/total_area;
		//beta = area_calc(x, y,  v2_NDC_x, v2_NDC_y, v0_NDC_x, v0_NDC_y)/total_area;
		//PCA
		gamma = area_calc(v0_NDC_x,v0_NDC_y, v1_NDC_x, v1_NDC_y,  x, y)/total_area;
		//gamma = area_calc(x, y, v1_NDC_x, v1_NDC_y, v0_NDC_x, v0_NDC_y )/total_area;
		//PAB
		
		for(int i = 0; i < state.floats_per_vertex; i++){
			if(state.interp_rules[i] == interp_type::flat){
				input.data[i] = v0.data[i];
			//	std:: cout << "-- " << v0.data[i];
			}
			else if (state.interp_rules[i] == interp_type::smooth){
				wp = 1/(alpha/v0.gl_Position[3] + beta/v1.gl_Position[3] + gamma/v2.gl_Position[3]);
				alphaP = (alpha*wp)/v0.gl_Position[3];
				betaP = (beta*wp)/v1.gl_Position[3];
				gammaP = (gamma*wp)/v2.gl_Position[3];
				input.data[i] = alphaP*v0.data[i] + betaP*v1.data[i] + gammaP*v2.data[i];
			
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

