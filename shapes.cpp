#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>

#include <iostream>
#include <bits/stdc++.h>
#include "utils.cpp"
#include "bitmap_image.hpp"

using namespace std;

bool isTexture=false;


#define SCALE 1.00


typedef struct {
    bool isHit = false;
    double t;
    color hitColor; // color of the object that was hit
    point hitPoint;
    vec normal;
    vec reflectionDirection;
    string hitObjectType;
    double ambientCoeff, diffuseCoeff, specularCoeff, reflectionCoeff, shininess;
} HitInfo;

typedef struct {
    point origin;
    vec direction;
} Ray;





class Shape
{
public:
    // double r, g, b;
    double ambientCoeff, diffuseCoeff, specularCoeff, reflectionCoeff;
    double shininess;
    Shape(double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess) {
        // this->r = r;
        // this->g = g;
        // this->b = b;
        this->ambientCoeff = ambient_coeff;
        this->diffuseCoeff = diffuse_coeff;
        this->specularCoeff = specular_coeff;
        this->reflectionCoeff = reflection_coeff;
        this->shininess = shininess;
    }
    virtual string getType() = 0;
    virtual HitInfo getHitInfo(Ray) = 0;
    virtual void draw() = 0;
    
};

class Triangle : public Shape
{
public:
    point p1, p2, p3;
    color triangleColor;
    Triangle(point p1, point p2, point p3, color triangle_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess)  : Shape(ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, shininess) {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        this->triangleColor = triangle_color;

        // cout<<"triangle color: "<<triangle_color.r<<" "<<triangle_color.g<<" "<<triangle_color.b<<endl;
        // cout<<"triangle color: "<<this->triangleColor.r<<" "<<this->triangleColor.g<<" "<<this->triangleColor.b<<endl;
        // cout<<"triange ambient coeff: "<<ambientCoeff<<endl;
        // cout<<"triange diffuse coeff: "<<diffuseCoeff<<endl;
        // cout<<"triange specular coeff: "<<specularCoeff<<endl;
        // cout<<"triange reflection coeff: "<<reflectionCoeff<<endl;
        // cout<<"triange shininess: "<<shininess<<endl;


    }

    // check if ray intersects with triangle
    HitInfo getHitInfo(Ray ray) {
        HitInfo ret;
        double beta, gamma, t;
        bool hasSolve = solve_simultaneous_equations(
            p1.x - p2.x, p1.x - p3.x, ray.direction.x,      p1.x - ray.origin.x,
            p1.y - p2.y, p1.y - p3.y, ray.direction.y,      p1.y - ray.origin.y,
            p1.z - p2.z, p1.z - p3.z, ray.direction.z,      p1.z - ray.origin.z,
            &beta,       &gamma,       &t
            );
        
        

        if(hasSolve && beta >= 0 && gamma >= 0 && beta + gamma <= 1 && t > 0) {
            ret.isHit = true;
            ret.hitObjectType="triangle";
            ret.t = t;
            

            vec normal = cross_product(p2 - p1, p3 - p1);
            if (dot_product(normal, ray.direction) > 0) {
                normal = {-normal.x, -normal.y, -normal.z};
            }

            normal = scale_to_r(normal, 1);

            ret.hitPoint = ray.origin + scale_to_r(ray.direction, t) + scale_to_r(normal, EPSILON);
            ret.reflectionDirection = get_reflection_direction(ray.direction, normal);
            ret.normal = scale_to_r(normal,1);
            ret.ambientCoeff = ambientCoeff;
            ret.diffuseCoeff = diffuseCoeff;
            ret.specularCoeff = specularCoeff;
            ret.reflectionCoeff = reflectionCoeff;
            ret.shininess = shininess;
            ret.hitColor = triangleColor;


        }
        
        
        return ret;
    }

    string getType() {
        return "triangle";
    }

    void draw() {
        glColor3d(triangleColor.r, triangleColor.g, triangleColor.b);
        glBegin(GL_TRIANGLES);
            glVertex3d(SCALE*p1.x, SCALE*p1.y, SCALE*p1.z);
            glVertex3d(SCALE*p2.x, SCALE*p2.y, SCALE*p2.z);
            glVertex3d(SCALE*p3.x, SCALE*p3.y, SCALE*p3.z);
        glEnd();
    }
};


class Square : public Shape
{
public:
    point p1, p2, p3, p4;
    color squareColor;
    Square(point p1, point p2, vec sense, color square_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess)  : Shape(ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, shininess) {
        this->p1 = p1;
        this->p2 = p2;
        double arm_len = mod(p2-p1);
        this->p3 = p2+scale_to_r(sense, arm_len);
        this->p4 = p1+scale_to_r(sense, arm_len);

        this->squareColor = square_color;

        

    }

    // check if ray intersects with square
    HitInfo getHitInfo(Ray ray) {
        HitInfo ret;
        double beta, gamma, t;
        bool hasSolve = solve_simultaneous_equations(
            p1.x - p2.x, p1.x - p4.x, ray.direction.x,      p1.x - ray.origin.x,
            p1.y - p2.y, p1.y - p4.y, ray.direction.y,      p1.y - ray.origin.y,
            p1.z - p2.z, p1.z - p4.z, ray.direction.z,      p1.z - ray.origin.z,
            &beta,       &gamma,       &t
            );
        
        

        if(hasSolve && beta >= 0 && gamma >= 0 && beta <= 1 && gamma <= 1 && t > 0) {
            ret.isHit = true;
            ret.hitObjectType="square";
            ret.t = t;
            

            vec normal = cross_product(p2 - p1, p3 - p1);
            if (dot_product(normal, ray.direction) > 0) {
                normal = {-normal.x, -normal.y, -normal.z};
            }

            normal = scale_to_r(normal, 1);

            ret.hitPoint = ray.origin + scale_to_r(ray.direction, t) + scale_to_r(normal, EPSILON);
            ret.reflectionDirection = get_reflection_direction(ray.direction, normal);
            ret.normal = normal;
            ret.ambientCoeff = ambientCoeff;
            ret.diffuseCoeff = diffuseCoeff;
            ret.specularCoeff = specularCoeff;
            ret.reflectionCoeff = reflectionCoeff;
            ret.shininess = shininess;
            ret.hitColor = squareColor;


        }
        
        
        return ret;
    }

    string getType() {
        return "triangle";
    }

    void draw() {
        glColor3d(squareColor.r, squareColor.g, squareColor.b);
        glBegin(GL_QUADS);
            glVertex3d(SCALE*p1.x, SCALE*p1.y, SCALE*p1.z);
            glVertex3d(SCALE*p2.x, SCALE*p2.y, SCALE*p2.z);
            glVertex3d(SCALE*p3.x, SCALE*p3.y, SCALE*p3.z);
            glVertex3d(SCALE*p4.x, SCALE*p4.y, SCALE*p4.z);
        glEnd();
    }
};


void drawSubsphere(double radius, int subdivisions) {
    long long num_points_per_row = 1<<subdivisions +1;
    struct point points[num_points_per_row][num_points_per_row];

    vec n1, n2;
    double theta1, theta2;

    for(int i = 0; i < num_points_per_row; i++) {
        theta1 = deg2rad(-45.0+90.0*i/(num_points_per_row-1));
        n1.x=-sin(theta1); n1.y=cos(theta1); n1.z=0;
        for(int j=0; j < num_points_per_row; j++) {
            theta2 = deg2rad(-45.0+90.0*j/(num_points_per_row-1));
            n2.x=sin(theta2); n2.y=0; n2.z=cos(theta2);

            points[i][j] = scale_to_r(cross_product(n1, n2), radius);
        }
    }

    // double c = 2.0/3;

    glBegin(GL_QUADS);
        for(int i = 0; i < num_points_per_row - 1; i++) {
            for(int j = 0; j < num_points_per_row - 1; j++) {
                // glColor3d(c,c,c);
                glVertex3d(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3d(points[i+1][j].x, points[i+1][j].y, points[i+1][j].z);
                glVertex3d(points[i+1][j+1].x, points[i+1][j+1].y, points[i+1][j+1].z);
                glVertex3d(points[i][j+1].x, points[i][j+1].y, points[i][j+1].z);
                // c+=1.0/3.0/((num_points_per_row - 1) * (num_points_per_row - 1));
            }
        }
    glEnd();
}

class CheckerBoard : public Shape
{
public:
    double cellWidth;
    
    bitmap_image blackTextureImage;
    bitmap_image whiteTextureImage;
    double maxX, maxY, minX, minY;
    CheckerBoard(double cell_width, double ambient_coeff, double diffuse_coeff, double reflection_coeff)  : Shape(ambient_coeff, diffuse_coeff, 0, reflection_coeff, 0) {
        this->cellWidth = cell_width;
        
        
        maxX = 1000;
        maxY = 1000;
        minX = -1000;
        minY = -1000;

        // this->blackTextureImage = new bitmap_image("texture_b.bmp");
        // this->whiteTextureImage = new bitmap_image("texture_w.bmp");
        this->blackTextureImage = bitmap_image("texture_b.bmp");
        this->whiteTextureImage = bitmap_image("texture_w.bmp");
        // this->isTexture = true;
    }

    color getWhiteTexture(double u, double v) {
        int x_pixel = u * (whiteTextureImage.width());
        int y_pixel = v * (whiteTextureImage.height());

        // cout<<"x_pixel: "<<x_pixel<<endl;
        // cout<<"y_pixel: "<<y_pixel<<endl;
        // cout<<"width: "<<whiteTextureImage.width()<<endl;
        // cout<<"height: "<<whiteTextureImage.height()<<endl;

        unsigned char r, g, b;
        whiteTextureImage.get_pixel(x_pixel, y_pixel, r, g, b);

        // cout<<"r: "<<(int)r<<endl;
        // cout<<"g: "<<(int)g<<endl;
        // cout<<"b: "<<(int)b<<endl;
        
        return {(double)r/255.0, (double)g/255.0, (double)b/255.0};
    }

    color getBlackTexture(double u, double v) {
        int x_pixel = u * (blackTextureImage.width());
        int y_pixel = v * (blackTextureImage.height());

        // cout<<"x_pixel: "<<x_pixel<<endl;
        // cout<<"y_pixel: "<<y_pixel<<endl;
        // cout<<"width: "<<blackTextureImage.width()<<endl;
        // cout<<"height: "<<blackTextureImage.height()<<endl;

        unsigned char r, g, b;
        blackTextureImage.get_pixel(x_pixel, y_pixel, r, g, b);

        // cout<<"r: "<<(int)r<<endl;
        // cout<<"g: "<<(int)g<<endl;
        // cout<<"b: "<<(int)b<<endl;
        
        return {(double)r/255.0, (double)g/255.0, (double)b/255.0};
    }

    string getType() {
        return "checkerboard";
    }

    // check if ray intersects with checkerboard
    HitInfo getHitInfo(Ray ray) {
        HitInfo ret;
        // plane eqn: z = 0
        // ray.origin.z + t* ray.direction.z = 0

        if (ray.direction.z <= EPSILON && ray.direction.z >= -EPSILON) return ret;

        double t = - ray.origin.z / ray.direction.z;
        if(t < 0) return ret;

        
        point hit_point = ray.origin + scale_to_r(ray.direction, t);
        vec normal{.x=0, .y=0, .z=1};
        if (dot_product(normal, ray.direction) > 0) {
            normal = {-normal.x, -normal.y, -normal.z};
        }

        // print_vec(normal);

        ret.isHit = true;
        ret.hitObjectType="checkerboard";
        ret.t = t;
        ret.hitPoint = hit_point+scale_to_r(normal, EPSILON);
        ret.reflectionDirection = get_reflection_direction(ray.direction, normal);
        ret.normal = normal;
        ret.ambientCoeff = ambientCoeff;
        ret.diffuseCoeff = diffuseCoeff;
        ret.specularCoeff = specularCoeff;
        ret.reflectionCoeff = reflectionCoeff;
        ret.shininess = shininess;
        ret.hitColor = getTextureAt(hit_point.x, hit_point.y);
        
        
        return ret;
    }

    color getTextureAt(double x, double y) {
        int i = (x - minX) / cellWidth;
        int j = (y - minY) / cellWidth;
        double u = x/cellWidth - floor(x / cellWidth);
        double v = y/cellWidth - floor(y / cellWidth);
        if ((i + j) % 2 == 0) {
            if(isTexture) { 
                return getWhiteTexture(u, v);
            }
            return {1, 1, 1};
        } else {
            if(isTexture) { 
                return getBlackTexture(u, v);
            }
            return {0, 0, 0};
        }
    }

    color getColorAt(double x, double y) {
        int i = (x - minX) / cellWidth;
        int j = (y - minY) / cellWidth;
        if ((i + j) % 2 == 0) {
            return {1, 1, 1};
        } else {
            return {0, 0, 0};
        }
    }

    void draw() {
        for (double x = minX; x < maxX; x += cellWidth) {
            for (double y = minY; y < maxY; y += cellWidth) {
                color c = getColorAt(x, y);
                glColor3d(c.r, c.g, c.b);
                glBegin(GL_QUADS);
                    glVertex3d(SCALE*x, SCALE*y, 0);
                    glVertex3d(SCALE*(x + cellWidth), SCALE*y, 0);
                    glVertex3d(SCALE*(x + cellWidth), SCALE*(y + cellWidth), 0);
                    glVertex3d(SCALE*x, SCALE*(y + cellWidth), 0);
                glEnd();
            }
        }
        
    }
};

class Sphere : public Shape
{
public:
    double radius;
    point center;
    color sphereColor;
    Sphere(point center, double radius, double r, double g, double b, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess) : Shape(ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, shininess) {
        this->center=center;
        this->radius = radius;

        this->sphereColor.r = r;
        this->sphereColor.g = g;
        this->sphereColor.b = b;
    }

    string getType() {
        return "sphere";
    }

    HitInfo getHitInfo(Ray ray) {
        HitInfo ret;
        vec ro = ray.origin - center;
        double tp = - dot_product(ray.direction, ro);
        double d2 = dot_product(ro, ro) - tp*tp;
        if(d2 > radius*radius) {
            return ret;
        }
        double t_prime = sqrt(radius*radius - d2);
        double t;
        vec normal;
        point hit_point;
        if(dot_product(ro, ro) < radius*radius) {
            // inside sphere
            t = tp + t_prime;
            hit_point = ray.origin + scale_to_r(ray.direction, t);
            normal = scale_to_r(center - hit_point, 1);
        } else {
            // outside sphere
            t = tp - t_prime;
            hit_point = ray.origin + scale_to_r(ray.direction, t);
            normal = scale_to_r(hit_point - center, 1);
        }

        if(t < 0) {
            return ret;
        }

        ret.isHit = true;
        ret.hitObjectType="sphere";
        ret.t = t;
        ret.hitPoint = hit_point+scale_to_r(normal, EPSILON);
        ret.reflectionDirection = get_reflection_direction(ray.direction, normal);
        ret.normal = normal;
        ret.ambientCoeff = ambientCoeff;
        ret.diffuseCoeff = diffuseCoeff;
        ret.specularCoeff = specularCoeff;
        ret.reflectionCoeff = reflectionCoeff;
        ret.shininess = shininess;
        ret.hitColor = sphereColor;
        return ret;
    }

    void draw() {
        // cout<<"drawing sphere"<<endl;
        glColor3d(sphereColor.r, sphereColor.g, sphereColor.b);
        
        glPushMatrix();
            glTranslated(center.x*SCALE, center.y*SCALE, center.z*SCALE);
            for(int i = 0; i < 4; i++) {
                glPushMatrix();
                    // glTranslated(x, y, z);
                    glRotated(90*i, 0, 0, 1);
                    drawSubsphere(radius*SCALE, 5);
                glPopMatrix();
            }
            glRotated(90, 0, 1, 0);
            drawSubsphere(radius*SCALE, 5);
            glRotated(-180, 0, 1, 0);
            drawSubsphere(radius*SCALE, 5);
        glPopMatrix();

    }
};


class Cube : public Shape
{
public:
    double side;
    // point topUpperLeft, topUpperRight; // upper -> away from eye or towards +ve y
    // point topLowerLeft, topLowerRight; // lower -> near from eye or towards -ve y

    
    // point bottomUpperLeft, bottomUpperRight;
    // point bottomLowerLeft, bottomLowerRight;

    point points[2][2][2];

    color cubeColor;
    vector<Shape *> subShapes;
    Cube(point bottom_lower_left, double side, color cube_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess) : Shape(ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, shininess) {
        this->side = side;

        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                for (int k=0; k<2; k++) {
                    points[i][j][k] = bottom_lower_left + (vec){side*i,side*j,side*k};
                }
            }
        }

        // this->p   = bottom_lower_left + (vec){0.0,0.0,0.0}; // to match my definition of bottom left with that of in the spec
        // this->bottomLowerRight  = bottom_lower_left + (vec){side,0.0,0.0};
        // this->bottomUpperLeft   = bottom_lower_left + (vec){0.0,side,0.0};
        // this->bottomUpperRight  = bottom_lower_left + (vec){side,side,0.0};

        // this->topLowerLeft      = bottom_lower_left + (vec){0.0,0.0,side};
        // this->topLowerRight     = bottom_lower_left + (vec){side,0.0,side};
        // this->topUpperLeft      = bottom_lower_left + (vec){0.0,side,side};
        // this->topUpperRight     = bottom_lower_left + (vec){side,side,side};
        
        

        this->cubeColor = cube_color;

        subShapes.push_back(new Square(points[0][0][0], points[1][0][0],
                                    (vec){0.0, 1.0, 0.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Square(points[0][0][0], points[1][0][0],
                                    (vec){0.0, 0.0, 1.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Square(points[0][1][1], points[1][1][1],
                                    (vec){0.0, -1.0, 0.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Square(points[0][1][1], points[1][1][1],
                                    (vec){0.0, 0.0, -1.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Square(points[1][0][0], points[1][1][0],
                                    (vec){0.0, 0.0, 1.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Square(points[0][0][0], points[0][1][0],
                                    (vec){0.0, 0.0, 1.0},
                                    cube_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));
    }

    ~Cube() {
        for (int i=0; i < subShapes.size(); i++) {
            delete subShapes[i];
        }
    }

    string getType() {
        return "cube";
    }

    HitInfo getHitInfo(Ray ray) {
        HitInfo min_hit_info;
        min_hit_info.t=MAXFLOAT;

        for (int i=0; i < subShapes.size(); i++) {
            HitInfo temp = subShapes[i]->getHitInfo(ray);
            if(temp.isHit && temp.t < min_hit_info.t) min_hit_info = temp;
        }
        
        return min_hit_info;
    }

    void draw() {
        // cout<<"drawing cube"<<endl;
        glColor3d(cubeColor.r, cubeColor.g, cubeColor.b);

        // front face
        glBegin(GL_QUADS);
            glVertex3d(points[0][0][0].x, points[0][0][0].y, points[0][0][0].z);
            glVertex3d(points[1][0][0].x, points[1][0][0].y, points[1][0][0].z);
            glVertex3d(points[1][0][1].x, points[1][0][1].y, points[1][0][1].z);
            glVertex3d(points[0][0][1].x, points[0][0][1].y, points[0][0][1].z);
        glEnd();

        // glColor3d(0.1,0.1,0.1);

        // back face
        glBegin(GL_QUADS);
            glVertex3d(points[0][1][0].x, points[0][1][0].y, points[0][1][0].z);
            glVertex3d(points[1][1][0].x, points[1][1][0].y, points[1][1][0].z);
            glVertex3d(points[1][1][1].x, points[1][1][1].y, points[1][1][1].z);
            glVertex3d(points[0][1][1].x, points[0][1][1].y, points[0][1][1].z);
        glEnd();

        // glColor3d(0.2,0.2,0.2);

        // right face
        glBegin(GL_QUADS);
            glVertex3d(points[1][0][0].x, points[1][0][0].y, points[1][0][0].z);
            glVertex3d(points[1][1][0].x, points[1][1][0].y, points[1][1][0].z);
            glVertex3d(points[1][1][1].x, points[1][1][1].y, points[1][1][1].z);
            glVertex3d(points[1][0][1].x, points[1][0][1].y, points[1][0][1].z);
        glEnd();

        // glColor3d(0.3,0.3,0.3);

        // left face
        glBegin(GL_QUADS);
            glVertex3d(points[0][0][0].x, points[0][0][0].y, points[0][0][0].z);
            glVertex3d(points[0][1][0].x, points[0][1][0].y, points[0][1][0].z);
            glVertex3d(points[0][1][1].x, points[0][1][1].y, points[0][1][1].z);
            glVertex3d(points[0][0][1].x, points[0][0][1].y, points[0][0][1].z);
        glEnd();

        // glColor3d(0.4,0.4,0.4);

        // top face
        glBegin(GL_QUADS);
            glVertex3d(points[0][0][1].x, points[0][0][1].y, points[0][0][1].z);
            glVertex3d(points[1][0][1].x, points[1][0][1].y, points[1][0][1].z);
            glVertex3d(points[1][1][1].x, points[1][1][1].y, points[1][1][1].z);
            glVertex3d(points[0][1][1].x, points[0][1][1].y, points[0][1][1].z);
        glEnd();

        // glColor3d(0.5,0.5,0.5);

        // bottom face
        glBegin(GL_QUADS);
            glVertex3d(points[0][0][0].x, points[0][0][0].y, points[0][0][0].z);
            glVertex3d(points[1][0][0].x, points[1][0][0].y, points[1][0][0].z);
            glVertex3d(points[1][1][0].x, points[1][1][0].y, points[1][1][0].z);
            glVertex3d(points[0][1][0].x, points[0][1][0].y, points[0][1][0].z);
        glEnd();

    }
};


class Pyramid : public Shape
{
public:
    double width, height;


    point basePoints[2][2];
    point topPoint;

    color pyramidColor;
    vector<Shape *> subShapes;
    Pyramid(point lowest_point, double width, double height, color pyramid_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double shininess) : Shape(ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, shininess) {
        point mid_point = lowest_point + (vec){width/2.0, width/2.0, 0};
        this->width = width;
        this->height = height;

        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                basePoints[i][j] = mid_point + (vec){width/2.0*(-1.0 + 2.0*i), width/2.0*(-1.0 + 2.0*j), 0};
                // print_vec(basePoints[i][j]);
            }
        }

        this->topPoint = mid_point+(vec){0.0, 0.0, height};
        this->pyramidColor = pyramid_color;

        subShapes.push_back(new Square(basePoints[0][0], basePoints[1][0],
                                    (vec){0.0, 1.0, 0.0},
                                    pyramid_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Triangle(basePoints[0][0], basePoints[1][0], topPoint,
                                    pyramid_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Triangle(basePoints[1][0], basePoints[1][1], topPoint,
                                    pyramid_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Triangle(basePoints[1][1], basePoints[0][1], topPoint,
                                    pyramid_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));

        subShapes.push_back(new Triangle(basePoints[0][1], basePoints[0][0], topPoint,
                                    pyramid_color,
                                    ambient_coeff,
                                    diffuse_coeff,
                                    specular_coeff,
                                    reflection_coeff,
                                    shininess));
    }

    ~Pyramid() {
        for (int i=0; i < subShapes.size(); i++) {
            delete subShapes[i];
        }
    }

    string getType() {
        return "pyramid";
    }

    HitInfo getHitInfo(Ray ray) {
        HitInfo min_hit_info;
        min_hit_info.t=MAXFLOAT;

        for (int i=0; i < subShapes.size(); i++) {
            HitInfo temp = subShapes[i]->getHitInfo(ray);
            if(temp.isHit && temp.t < min_hit_info.t) min_hit_info = temp;
        }
        
        return min_hit_info;
    }

    void draw() {
        glColor3d(pyramidColor.r, pyramidColor.g, pyramidColor.b);
        glBegin(GL_QUADS);
            glVertex3d(basePoints[0][0].x, basePoints[0][0].y, basePoints[0][0].z);
            glVertex3d(basePoints[1][0].x, basePoints[1][0].y, basePoints[1][0].z);
            glVertex3d(basePoints[1][1].x, basePoints[1][1].y, basePoints[1][1].z);
            glVertex3d(basePoints[0][1].x, basePoints[0][1].y, basePoints[0][1].z);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3d(basePoints[0][0].x, basePoints[0][0].y, basePoints[0][0].z);
            glVertex3d(basePoints[1][0].x, basePoints[1][0].y, basePoints[1][0].z);
            glVertex3d(topPoint.x, topPoint.y, topPoint.z);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3d(basePoints[1][0].x, basePoints[1][0].y, basePoints[1][0].z);
            glVertex3d(basePoints[1][1].x, basePoints[1][1].y, basePoints[1][1].z);
            glVertex3d(topPoint.x, topPoint.y, topPoint.z);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3d(basePoints[1][1].x, basePoints[1][1].y, basePoints[1][1].z);
            glVertex3d(basePoints[0][1].x, basePoints[0][1].y, basePoints[0][1].z);
            glVertex3d(topPoint.x, topPoint.y, topPoint.z);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3d(basePoints[0][1].x, basePoints[0][1].y, basePoints[0][1].z);
            glVertex3d(basePoints[0][0].x, basePoints[0][0].y, basePoints[0][0].z);
            glVertex3d(topPoint.x, topPoint.y, topPoint.z);
        glEnd();
    }
};



class NormalLight
{
public:
    point position;
    double fallOff;
    NormalLight(double x, double y, double z, double fallOff) {
        this->position.x = x;
        this->position.y = y;
        this->position.z = z;
        this->fallOff = fallOff;
    }

    void draw() {
        // cout<<"drawing sphere"<<endl;
        double radius = 5;
        // glColor3d(1, 1, 1);
        
        glPushMatrix();
            glTranslated(position.x, position.y, position.z);
            for(int i = 0; i < 4; i++) {
                glPushMatrix();
                    // glTranslated(x, y, z);
                    glRotated(90*i, 0, 0, 1);
                    glColor3d(0.4+0.1*i, 0.4+0.1*i, 0.4+0.1*i);
                    drawSubsphere(radius, 5);
                glPopMatrix();
            }
            glRotated(90, 0, 1, 0);
            glColor3d(0.8, 0.8, 0.8);
            drawSubsphere(radius, 5);
            glColor3d(0.9, 0.9, 0.9);
            glRotated(-180, 0, 1, 0);
            drawSubsphere(radius, 5);
        glPopMatrix();

    }

};

class SpotLight
{
public:
    point position;
    double fallOff;
    vec lookAt;
    double cutOffAngle;
    SpotLight(double x, double y, double z, double fallOff, double lookAtX, double lookAtY, double lookAtZ, double cutOffAngle) {
        this->position.x = x;
        this->position.y = y;
        this->position.z = z;
        this->fallOff = fallOff;
        // this->lookAt.x = lookAtX;
        // this->lookAt.y = lookAtY;
        // this->lookAt.z = lookAtZ;
        this->lookAt = scale_to_r((point){lookAtX, lookAtY, lookAtZ} - this->position, 1);
        this->cutOffAngle = cutOffAngle;
    }

    void draw() {
        // drawing a cone

        double height = 30;
        double radius = 5;

        glPushMatrix();
            glTranslated(position.x, position.y, position.z);
            vec rotation_axis =  cross_product(lookAt, (vec){0, 0, 1});
            glRotated(-rad2deg(acos(dot_product(lookAt, (vec){0, 0, 1}))), rotation_axis.x, rotation_axis.y, rotation_axis.z);

            glBegin(GL_POLYGON);
                for (int i = 0; i < 360; i=i+10)
                {
                    double theta = deg2rad((double)i);; //get the current angle

                    double x = radius * cos(theta); //calculate the x component
                    double y = radius * sin(theta); //calculate the y component

                    glVertex3d(x, y, 0.0); //output vertex

                }
            glEnd();

            for(int i = 0; i < 360; i=i+10) {
                double theta = deg2rad((double)i); //get the current angle

                double x1 = radius * cos(theta); //calculate the x component
                double y1 = radius * sin(theta); //calculate the y component

                double x2 = radius * cos(deg2rad((double)i+10.0)); //calculate the x component
                double y2 = radius * sin(deg2rad((double)i+10.0)); //calculate the y component

                double c = 0.5 + 0.5 * cos(theta);

                glColor3d(c, c, c);

                glBegin(GL_TRIANGLES);
                    glVertex3d(x1, y1, 0);
                    glVertex3d(x2, y2, 0);
                    glVertex3d(0, 0, height);
                glEnd();
            }

        glPopMatrix();

        
    }

    bool isPointInCone(point p) {
        vec v = scale_to_r(p - position, 1);

        double angle = rad2deg(acos(dot_product(v, lookAt)));
        return angle <= cutOffAngle;
    }

};
