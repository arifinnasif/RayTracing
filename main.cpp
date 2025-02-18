#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h>
#include "shapes.cpp"
#include "utils.cpp"
#include "bitmap_image.hpp"

#include <iostream>

using namespace std;


extern bool isTexture;





// camera global variables
struct point eye = {.x=0, .y=100.0, .z=50};
struct point center = {.x=0, .y=0, .z=0};
struct point up = {.x=0, .y=0, .z=100};

vec l, perp_to_l_hor, perp_to_l_ver; // calculate immediately after data load





// fustrum global variables
double near_dist;
double far_dist;
double fov_y;
double aspect_ratio;





// ray tracing global variables
int recursion_level;
int image_width;
int image_height;
point** point_buffer;
color** color_buffer;






// objects global variables
int number_of_objects;
int number_of_normal_lights;
int number_of_spot_lights;
double checkerboard_cell_width;
vector<Shape *> objects;
vector<NormalLight *> normal_lights;
vector<SpotLight *> spot_lights;
CheckerBoard* checkerboard;




// check if normal light illuminates hit point
bool is_illuminated(point hit_point, NormalLight* normal_light) {
    vec temp = normal_light->position - hit_point;
    vec direction = scale_to_r(temp, 1);
    Ray ray{
        .origin = hit_point,
        .direction = direction
    };
    for (int i = 0; i < objects.size(); i++)
    {
        double distance_to_light_sqrd = dot_product(temp, temp);
        HitInfo hit_info = objects[i]->getHitInfo(ray);
        if (hit_info.isHit && hit_info.t * hit_info.t < distance_to_light_sqrd) { // check t is less than distance to light source
            // cout<<"end   "<<hit_info.hitObjectType<<endl;
            return false;
        }
    }
    // cout<<"end   light"<<endl;
    return true;
}

// check if spot light illuminates hit point
bool is_illuminated(point hit_point, SpotLight* spot_light) {
    vec temp = spot_light->position - hit_point;
    vec direction = scale_to_r(temp, 1);
    Ray ray{
        .origin = hit_point,
        .direction = direction
    };
    for (int i = 0; i < objects.size(); i++)
    {
        double distance_to_light_sqrd = dot_product(temp, temp);
        HitInfo hit_info = objects[i]->getHitInfo(ray);
        if (hit_info.isHit && hit_info.t * hit_info.t < distance_to_light_sqrd) {
            return false;
        }
    }
    return spot_light->isPointInCone(hit_point);
}





void initialize_buffers() {
    point_buffer = new point*[image_height];
    color_buffer = new color*[image_height];
    for (int i = 0; i < image_height; i++)
    {
        point_buffer[i] = new point[image_width];
        color_buffer[i] = new color[image_width];
    }

    double vert = near_dist * tan(deg2rad(fov_y/2));
    double hor = vert * aspect_ratio;

    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            point temp = eye + scale_to_r(l,near_dist)
            + scale_to_r(perp_to_l_hor, hor * (2.0*(j+0.5)/image_width - 1))
            + scale_to_r(perp_to_l_ver, vert * (1.0 - 2.0*(i+0.5)/image_height)); // 0.5 is added to get middle of pixel

            point_buffer[i][j].x = temp.x;
            point_buffer[i][j].y = temp.y;
            point_buffer[i][j].z = temp.z;

            // color buffer to background color
            color_buffer[i][j].r = 0;
            color_buffer[i][j].g = 0;
            color_buffer[i][j].b = 0;
        }
    }

    cout<<"Point buffer generation done"<<endl;
    
}


void clear_buffers() {
    for (int i = 0; i < image_height; i++)
    {
        delete[] point_buffer[i];
        delete[] color_buffer[i];
    }
    delete[] point_buffer;
    delete[] color_buffer;
}

color ray_tracing(Ray ray, int arg_recursion_level) {
    if (arg_recursion_level == 0) return {0,0,0};
    color ret{0, 0, 0};
    // int min_index = -1;
    HitInfo min_hit_info;
    min_hit_info.t = MAXFLOAT;
    for (int i = 0; i < objects.size(); i++)
    {
        HitInfo hit_info = objects[i]->getHitInfo(ray);
        if (hit_info.isHit && hit_info.t < min_hit_info.t) {
            // min_index = i;
            min_hit_info = hit_info;
        }
    }
    if (!min_hit_info.isHit) return {0,0,0}; // return background color
    // phong model
    // ambient
    ret = {
        .r = min_hit_info.hitColor.r * min_hit_info.ambientCoeff,
        .g = min_hit_info.hitColor.g * min_hit_info.ambientCoeff,
        .b = min_hit_info.hitColor.b * min_hit_info.ambientCoeff
    };

    // diffuse and specular
    /**
     * @todo: add spotlight
     */
    double phong, lambert;
    phong = 0;
    lambert = 0;
    // point hit_point_moved_towards_eye = min_hit_info.hitPoint + scale_to_r(eye - min_hit_info.hitPoint, EPSILON*2.0);
    // point hit_point_moved_towards_normal = min_hit_info.hitPoint + scale_to_r(min_hit_info.normal, EPSILON);
    // cout<<"debug"<<endl;
    for (int i=0; i < normal_lights.size(); i++) {
        // cout<<"checking illumination"<<endl;
        // cout<<"start "<<min_hit_info.hitObjectType<<endl;
        if(!is_illuminated(min_hit_info.hitPoint, normal_lights[i])) continue;
        // cout<<"illuminated"<<endl;
        vec temp = normal_lights[i]->position - min_hit_info.hitPoint;
        double distance_sqrd_to_source = dot_product(temp, temp);
        vec source_direction_from_hit_point = scale_to_r(temp, 1);
        double scaling_factor = exp(-distance_sqrd_to_source*normal_lights[i]->fallOff);

        
        lambert += dot_product(source_direction_from_hit_point, min_hit_info.normal)*scaling_factor;
        phong += pow(dot_product(source_direction_from_hit_point, min_hit_info.reflectionDirection), min_hit_info.shininess)*scaling_factor;
    }

    // if ( min_hit_info.hitObjectType == "sphere")  {
    //     cout<<lambert<<endl;
    //     cout<<min_hit_info.diffuseCoeff<<endl;
    // }

    for (int i=0; i < spot_lights.size(); i++) {
        // cout<<"checking illumination"<<endl;
        // cout<<"start "<<min_hit_info.hitObjectType<<endl;
        if(!is_illuminated(min_hit_info.hitPoint, spot_lights[i])) continue;
        // cout<<"illuminated"<<endl;
        vec temp = spot_lights[i]->position - min_hit_info.hitPoint;
        double distance_sqrd_to_source = dot_product(temp, temp);
        vec source_direction_from_hit_point = scale_to_r(temp, 1);
        double scaling_factor = exp(-distance_sqrd_to_source*spot_lights[i]->fallOff);
        
        lambert += dot_product(source_direction_from_hit_point, min_hit_info.normal)*scaling_factor;
        phong += pow(dot_product(source_direction_from_hit_point, min_hit_info.reflectionDirection), min_hit_info.shininess)*scaling_factor;
    }

    // if (min_hit_info.hitObjectType == "sphere") cout<<lambert<<endl;

    

    // diffuse
    ret = {
        .r = ret.r + min_hit_info.diffuseCoeff*lambert*min_hit_info.hitColor.r,
        .g = ret.g + min_hit_info.diffuseCoeff*lambert*min_hit_info.hitColor.g,
        .b = ret.b + min_hit_info.diffuseCoeff*lambert*min_hit_info.hitColor.b,
    };

    // cout<<min_hit_info.diffuseCoeff*lambert*min_hit_info.hitColor.r<<endl;
    // cout<<lambert<<endl;

    // specular
    ret = {
        .r = ret.r + min_hit_info.specularCoeff*phong*min_hit_info.hitColor.r,
        .g = ret.g + min_hit_info.specularCoeff*phong*min_hit_info.hitColor.g,
        .b = ret.b + min_hit_info.specularCoeff*phong*min_hit_info.hitColor.b,
    };
    
    Ray reflected_ray{
        .origin = min_hit_info.hitPoint,
        .direction = min_hit_info.reflectionDirection
    };
    color reflected_color = ray_tracing(reflected_ray, arg_recursion_level - 1);
    ret = {
        .r = ret.r + reflected_color.r * min_hit_info.reflectionCoeff,
        .g = ret.g + reflected_color.g * min_hit_info.reflectionCoeff,
        .b = ret.b + reflected_color.b * min_hit_info.reflectionCoeff
    };
    return ret;
}

void read_description_file() {
    // read description file in c++
    ifstream description_file("description.txt");
    string line;
    // read until file ends
    cout<<"reading description file"<<endl;
    double checkerboard_cell_width;
    double checkerboard_ambient_coeff;
    double checkerboard_diffuse_coeff;
    double checkerboard_reflection_coeff;
    description_file>>near_dist>>far_dist>>fov_y>>aspect_ratio;
    description_file>>recursion_level;
    description_file>>image_width;
    image_height = image_width / aspect_ratio;
    description_file>>checkerboard_cell_width;
    description_file>>checkerboard_ambient_coeff>>checkerboard_diffuse_coeff>>checkerboard_reflection_coeff;

    checkerboard = new CheckerBoard(
        checkerboard_cell_width,
        checkerboard_ambient_coeff,
        checkerboard_diffuse_coeff,
        checkerboard_reflection_coeff
    );
    objects.push_back(checkerboard);
    // point p1 = {80,0,0};
    // point p2 = {80,80,0};
    // point p3 = {80,0,80};
    // point sqr = {0,0,80};
    // color c = {.9,0,.5};

    // objects.push_back(new Square(p1,p2,sqr,c,0.2,0.3,0.1,0.3,30));
    // objects.push_back(new Triangle(p1,p2,p3,c,0.2,0.3,0.1,0.3,30));
    description_file>>number_of_objects;

    for (int i = 0; i < number_of_objects; i++)
    {
        description_file>>line;
        if (line == "sphere") {
            point center;
            double radius;
            double r, g, b;
            double ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff;
            double shininess;

            description_file>>center.x>>center.y>>center.z;
            description_file>>radius;
            description_file>>r>>g>>b;
            description_file>>ambient_coeff>>diffuse_coeff>>specular_coeff>>reflection_coeff;
            description_file>>shininess;

            cout<<"sphere"<<endl;
            objects.push_back(new Sphere(
                center, radius,
                r, g, b,
                ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff,
                shininess
                ));
            

        } else if (line == "cube") {
            point bottom_left;
            double side;
            color cube_color;
            double ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff;
            double shininess;

            description_file>>bottom_left.x>>bottom_left.y>>bottom_left.z;
            description_file>>side;
            description_file>>cube_color.r>>cube_color.g>>cube_color.b;
            description_file>>ambient_coeff>>diffuse_coeff>>specular_coeff>>reflection_coeff;
            description_file>>shininess;

            objects.push_back(new Cube(
                bottom_left, side,
                cube_color,
                ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff,
                shininess
            ));


            
            cout<<"cube"<<endl;
        } else if (line == "pyramid") {
            point lowest_point;
            double width, height;
            color pyramid_color;
            double ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff;
            double shininess;

            description_file>>lowest_point.x>>lowest_point.y>>lowest_point.z;
            description_file>>width>>height;
            description_file>>pyramid_color.r>>pyramid_color.g>>pyramid_color.b;
            description_file>>ambient_coeff>>diffuse_coeff>>specular_coeff>>reflection_coeff;
            description_file>>shininess;

            objects.push_back(new Pyramid(
                lowest_point,
                width, height,
                pyramid_color,
                ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff,
                shininess
            ));
            
            cout<<"pyramid"<<endl;
        } else {
            cout<<"unknown object"<<endl;
        }
        

        // std::istringstream iss(line);
        // int a, b;
        // if (!(iss >> a >> b)) { break; } // error

        // process pair (a,b)
        // cout<<line<<endl;
        // cout<<"---"<<endl;
    }

    description_file>>number_of_normal_lights;
    
    for(int i = 0; i < number_of_normal_lights; i++) {
        double x, y, z;
        double fall_off_rate;

        description_file>>x>>y>>z;
        description_file>>fall_off_rate;

        normal_lights.push_back(new NormalLight(x,y,z,fall_off_rate));
        
        cout<<"normal light"<<endl;
    }

    description_file>>number_of_spot_lights;

    for(int i = 0; i < number_of_spot_lights; i++) {
        double x, y, z;
        double fall_off_rate;
        double look_at_x, look_at_y, look_at_z;
        double cut_off_angle;

        description_file>>x>>y>>z;
        description_file>>fall_off_rate;
        description_file>>look_at_x>>look_at_y>>look_at_z;
        description_file>>cut_off_angle;

        spot_lights.push_back(new SpotLight(x,y,z,fall_off_rate, look_at_x, look_at_y, look_at_z, cut_off_angle));
        
        cout<<"spot light"<<endl;
    }

    cout<<"Done taking input"<<endl;

}


void initialize_parameters() {
    l=scale_to_r(center-eye, 1);
    perp_to_l_hor = scale_to_r(cross_product(l,up), 1);
    perp_to_l_ver = scale_to_r(cross_product(cross_product(l,up), l), 1);

}

// we will perform ray tracing here
void rt_handler(){
    cout<<"handling ray tracing"<<endl;
    initialize_buffers();
    int total_pixels = image_height * image_width;
    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            Ray ray{
                .origin = point_buffer[i][j],
                .direction = scale_to_r(point_buffer[i][j] - eye, 1)
            };
            color color = ray_tracing(ray, recursion_level);
            color_buffer[i][j].r = color.r;
            color_buffer[i][j].g = color.g;
            color_buffer[i][j].b = color.b;

            int done_pixels = i*image_width + j + 1;
            if (done_pixels % (total_pixels/10) == 0) {
                cout<<"Rendering "<<done_pixels*10/total_pixels*10<<"% complete"<<endl;
            }
        }
        
    }
    cout<<"Rendering "<<100<<"% complete"<<endl;

    // write to image
    bitmap_image image(image_width, image_height);
    for (int i = 0; i < image_height; i++)
    {
        for (int j = 0; j < image_width; j++)
        {
            image.set_pixel(j, i, min(color_buffer[i][j].r*255, 255.0), min(color_buffer[i][j].g*255, 255.0), min(color_buffer[i][j].b*255, 255.0));
        }
    }
    image.save_image("out.bmp");
    cout<<"Image Saved"<<endl;
    clear_buffers();
    cout<<"ray tracing done"<<endl;
}

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {
    double v = 1;
    vec l_prev;
    l_prev.x = center.x - eye.x;
    l_prev.y = center.y - eye.y;
    l_prev.z = center.z - eye.z;

    l = scale_to_r(l_prev,1);
    perp_to_l_hor = scale_to_r(cross_product(l_prev,up), 1);
    perp_to_l_ver = cross_product(cross_product(l,up), l);
    perp_to_l_ver = scale_to_r(perp_to_l_ver, 1);
    if(dot_product(perp_to_l_ver, up) > 0) up=perp_to_l_ver;
    else if(dot_product(perp_to_l_ver, up) < 0) {
        up.x=-perp_to_l_ver.x;
        up.y=-perp_to_l_ver.y;
        up.z=-perp_to_l_ver.z;
    }
    // print_vec(up);

    switch (key) {
        case 'w':
            eye.y+=v;
            break;

        case 's':
            eye.y-=v;
            break;

        case '0':
            // handle ray tracing
            rt_handler();

            break;

        case '1':
            center.x-=perp_to_l_hor.x*v;
            center.y-=perp_to_l_hor.y*v;
            center.z-=perp_to_l_hor.z*v;
            break;

        case '2':
            center.x+=perp_to_l_hor.x*v;
            center.y+=perp_to_l_hor.y*v;
            center.z+=perp_to_l_hor.z*v;
            break;

        case '3':
            center.x+=perp_to_l_ver.x*v;
            center.y+=perp_to_l_ver.y*v;
            center.z+=perp_to_l_ver.z*v;
            break;

        case '4':
            center.x-=perp_to_l_ver.x*v;
            center.y-=perp_to_l_ver.y*v;
            center.z-=perp_to_l_ver.z*v;
            break;

        case '5':
            up.x+=perp_to_l_hor.x*v;
            up.y+=perp_to_l_hor.y*v;
            up.z+=perp_to_l_hor.z*v;
            break;

        case '6':
            up.x-=perp_to_l_hor.x*v;
            up.y-=perp_to_l_hor.y*v;
            up.z-=perp_to_l_hor.z*v;
            break;
        
        case ' ':
            isTexture = !isTexture;
            break;

        case 'q':
            for(int i = 0; i < objects.size(); i++) {
                delete objects[i];
            }

            for(int i = 0; i < normal_lights.size(); i++) {
                delete normal_lights[i];
            }

            for(int i = 0; i < spot_lights.size(); i++) {
                delete spot_lights[i];
            }
            exit(0);   // exit
            break;
        default:
            return;


    }
    glutPostRedisplay();    // Post a paint request to activate display()
}



/* Callback handler for special-key event */
void specialKeyListener(int key, int x,int y) {
    double v = 1;
    vec l_prev;
    l_prev.x = center.x - eye.x;
    l_prev.y = center.y - eye.y;
    l_prev.z = center.z - eye.z;

    l = scale_to_r(l_prev,1);
    perp_to_l_hor = scale_to_r(cross_product(l,up), 1);
    perp_to_l_ver = cross_product(cross_product(l,up), l);
    perp_to_l_ver = scale_to_r(perp_to_l_ver, 1);
    if(dot_product(perp_to_l_ver, up) > 0) up=perp_to_l_ver;
    else if(dot_product(perp_to_l_ver, up) < 0) {
        up.x=-perp_to_l_ver.x;
        up.y=-perp_to_l_ver.y;
        up.z=-perp_to_l_ver.z;
    }
    // print_vec(up);

    switch (key) {

    case GLUT_KEY_UP:
        eye.x+=l.x*v;
        eye.y+=l.y*v;
        eye.z+=l.z*v;

        center.x+=l.x*v;
        center.y+=l.y*v;
        center.z+=l.z*v;
        break;
    case GLUT_KEY_DOWN:
        eye.x-=l.x*v;
        eye.y-=l.y*v;
        eye.z-=l.z*v;

        center.x-=l.x*v;
        center.y-=l.y*v;
        center.z-=l.z*v;
        break;
    case GLUT_KEY_RIGHT:
        eye.x+=perp_to_l_hor.x*v;
        eye.y+=perp_to_l_hor.y*v;
        eye.z+=perp_to_l_hor.z*v;

        center.x+=perp_to_l_hor.x*v;
        center.y+=perp_to_l_hor.y*v;
        center.z+=perp_to_l_hor.z*v;
        // print_vec(perp_to_l_hor);
        break;
    case GLUT_KEY_LEFT:
        eye.x-=perp_to_l_hor.x*v;
        eye.y-=perp_to_l_hor.y*v;
        eye.z-=perp_to_l_hor.z*v;

        center.x-=perp_to_l_hor.x*v;
        center.y-=perp_to_l_hor.y*v;
        center.z-=perp_to_l_hor.z*v;
        // print_vec(perp_to_l_hor);
        break;
    case GLUT_KEY_PAGE_UP:
        eye.x+=up.x*v;
        eye.y+=up.y*v;
        eye.z+=up.z*v;

        center.x+=up.x*v;
        center.y+=up.y*v;
        center.z+=up.z*v;
        // print_vec(up);
        break;
    case GLUT_KEY_PAGE_DOWN:
        eye.x-=up.x*v;
        eye.y-=up.y*v;
        eye.z-=up.z*v;

        center.x-=up.x*v;
        center.y-=up.y*v;
        center.z-=up.z*v;
        break;
    
    default:
        break;
    }
    // cout<<"---"<<endl;
    // cout<<"eye : ";
    // print_vec(eye);
    // cout<<"center : ";
    // print_vec(center);
    // cout<<"up : ";
    // print_vec(up);
    // cout<<"l : ";
    // print_vec(l);
    // cout<<"perp_to_l_hor : ";
    // print_vec(perp_to_l_hor);
    // cout<<"perp_to_l_ver : ";
    // print_vec(perp_to_l_ver);
    // cout<<"---"<<endl;
    glutPostRedisplay();    // Post a paint request to activate display()
}


void drawAxes() {
    glLineWidth(20);
    double length = 50;
    glBegin(GL_LINES);
        glColor3f(1,0,0);   // Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(length,0,0);

        glColor3f(0,1,0);   // Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,length,0);

        glColor3f(0,0,1);   // Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,length);


    glEnd();
}


/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    gluLookAt(SCALE*eye.x, SCALE*eye.y, SCALE*eye.z,
              SCALE*center.x, SCALE*center.y, SCALE*center.z,
              up.x, up.y, up.z);
    // draw
    drawAxes();
    
    // ************************************************************************
    // *                                                                      *
    // *                             DRAW HERE                                *
    // *                                                                      *
    // ************************************************************************
    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }

    for (int i = 0; i < normal_lights.size(); i++)
    {
        normal_lights[i]->draw();
    }

    for (int i = 0; i < spot_lights.size(); i++)
    {
        spot_lights[i]->draw();
    }

    
    // checkerboard->draw();

    glutSwapBuffers();  // Render now
}


/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
}


/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(fov_y, aspect_ratio, near_dist, far_dist);
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    read_description_file();
    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitWindowSize(image_width, image_height);               // Set the window's initial width & height
    glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("Ray Tracing");      // Create a window with the given title
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event
    initGL();                                   // Our own OpenGL initialization
    // initialization codes
    initialize_parameters();
    glutMainLoop();                             // Enter the event-processing loop
    return 0;
}
