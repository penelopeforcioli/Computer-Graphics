#define _CRT_SECURE_NO_WARNINGS 1
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <cfloat>
#include <tuple>
#include <omp.h>
#include <list>
#include <chrono>
#include <iostream>

#include "stb_image_write.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0.0, 1.0);


/* --------------------------------- AUX FUNCTIONS AND CLASSES --------------------------------------------------------------------------*/

double sqr(double x) {
 return x*x;
}

void BoxMuller(double stdev, double& x, double& y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = stdev * sqrt(-2*log(r1)) * cos(2*M_PI*r2);
    y = stdev * sqrt(-2*log(r1)) * sin(2*M_PI*r2);
}

// from the pastebin ressource on moodle
class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};



/* --------------------------------- VECTOR CLASS AS WELL AS HELPER FUNCTIONS FOR DEALING WITH VECTORS ------------------------------------*/
class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator-(const Vector& a) {
 return Vector(-a[0], -a[1], -a[2]);
}

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}




/* --------------------------------- INTERSECTION STRUCTURE ------------------------------------------------------------------------------*/

struct Intersection {
    Vector P;
    Vector N;
    Vector color;
    double t = 0.;
    bool reflective;
    bool intersects = false;

    Intersection( Vector color = Vector(0., 0., 0.),  bool reflects = false){
        this->color = color;
       
        this->reflective = reflects;
    }
};



/* --------------------------------- RAY CLASS --------------------------------------------------------------------------------------------*/
class Ray {
public:
	Vector O,u;
	Ray( const Vector& o,  const Vector& U):O(o), u(U) {}
};





/*  --------------------------------- GEOMETRY --------------------------------------------------------------------------------------------*/

class Geometry {
public:
    virtual Intersection intersect(const Ray& ray) = 0;
    Vector albedo;
    bool reflective;
    bool transparent;  
};

/*  --------------------------------- BOUNDINGBOX ------------------------------------------------------------------------------------------------*/

class BoundingBox {
public:
    Vector b_min;
    Vector b_max;

    BoundingBox(Vector min = Vector(), Vector max = Vector()) {
        b_min = min;
        b_max = max;
    }

    
    // helped by a friend here @lasocki but changed it for my own code structure
    bool intersect(const Ray &ray,double* t){
        double t0[3],t1[3];
        for (int i = 0; i<3; i++){
            t0[i] = std::min( (b_min[i] - ray.O[i])/ ray.u[i], (b_max[i] - ray.O[i]) / ray.u[i] );;
            t1[i] = std::max( (b_min[i] - ray.O[i])/ ray.u[i], (b_max[i] - ray.O[i]) / ray.u[i] ); }
        double max, min;
        max = std::max(t0[0],std::max(t0[1],t0[2]));
        min = std::min(t1[0],std::min(t1[1],t1[2]));
        if (min > max and max > 0){
            *t= max;
            return true;}
        *t= max;
        return false;
    }

};

bool Möller_Trumbore(const Ray& ray, const Vector& A, const Vector& B, const Vector& C, double &alpha, double &beta, double &gamma, double &t, Vector &N){
	Vector e1,e2;
	e1 = B-A;
	e2 = C-A;
	N = cross(e1,e2);

	double dv = dot(ray.u,N);

	beta = dot(e2,cross((A - ray.O),ray.u)) / dv;
	gamma = - dot(e1,cross((A - ray.O),ray.u)) /dv;
	alpha = 1-beta-gamma;

	t = dot(A-ray.O,N)/dot(ray.u,N);

	if (alpha<0 or alpha>1 or beta<0 or beta>1 or gamma<0 or gamma>1 or t<0){
		return false;
	}
	return true;
}

struct Node {
    BoundingBox bbox;
    int starting_triangle;
    int ending_triangle;
    Node* child_left;
    Node* child_right;
};

/*  --------------------------------- MESH ------------------------------------------------------------------------------------------------*/

class TriangleMesh : public Geometry {


public:

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    Node* root;
    double scaling_factor;
    Vector translation;

    ~TriangleMesh() {}
    TriangleMesh(double scaling_factor, Vector translation, Vector color = Vector(0., 0., 0.), bool reflective = false, bool transparent = false) {
        this->scaling_factor = scaling_factor;
        this->translation = translation;
        this->albedo = color;
        this->reflective = reflective;
        this->transparent = transparent;
        this->root = new Node;
    };



    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);
        this->buildBVH(root, 0, indices.size());
    }
  

BoundingBox bounding_box(int starting_triangle, int ending_triangle) {
    // got advised to use these macros by @Lasocki
    double minX = DBL_MAX, minY = DBL_MAX, minZ = DBL_MAX;
    double maxX = DBL_MIN, maxY = DBL_MIN, maxZ = DBL_MIN;
    
    for (int i = starting_triangle; i < ending_triangle; i++) {
        std::vector<Vector> triangle_vertices;
        triangle_vertices.push_back(this->vertices[this->indices[i].vtxi]);
        triangle_vertices.push_back(this->vertices[this->indices[i].vtxj]);
        triangle_vertices.push_back(this->vertices[this->indices[i].vtxk]);
        
        for (size_t j = 0; j < triangle_vertices.size(); j++) {
            Vector Vertex = scaling_factor * triangle_vertices[j] + translation;
            if (Vertex[0] < minX) minX = Vertex[0];
            if (Vertex[0] > maxX) maxX = Vertex[0];
            if (Vertex[1] < minY) minY = Vertex[1];
            if (Vertex[1] > maxY) maxY = Vertex[1];
            if (Vertex[2] < minZ) minZ = Vertex[2];
            if (Vertex[2] > maxZ) maxZ = Vertex[2];
        }
    }
    return BoundingBox(Vector(minX, minY, minZ), Vector(maxX, maxY, maxZ));
}


    Vector BaryCenter(int i) {
        Vector res = (scaling_factor * this->vertices[this->indices[i].vtxi]) + translation;
        res = res + (scaling_factor * this->vertices[this->indices[i].vtxj]) + translation;
         res = res +(scaling_factor * this->vertices[this->indices[i].vtxk]) + translation;
        return res /3.;
    }

    
   void buildBVH(Node *node, int starting_triangle, int ending_triangle) {
    node->bbox = this->bounding_box(starting_triangle, ending_triangle);
    node->starting_triangle = starting_triangle;
    node->ending_triangle = ending_triangle;

    Vector diag = node->bbox.b_max - node->bbox.b_min;
    Vector middle_diag = node->bbox.b_min + diag*0.5;
    int longest_axis;
    if (abs(diag[0]) > abs(diag[1]) && abs(diag[0]) > abs(diag[2])) {
        longest_axis = 0; } else {
        if (abs(diag[1]) > abs(diag[2])) {
            longest_axis = 1;} else {
            longest_axis = 2;
        }
    }

    int pivot_index = starting_triangle;

    for (int i = starting_triangle; i < ending_triangle; i++) {
        Vector barycenter = this->BaryCenter(i);
        if (barycenter[longest_axis] < middle_diag[longest_axis]) {
            std::swap(indices[i], indices[pivot_index]);
            pivot_index++;
        }
    }

    if (pivot_index <= starting_triangle || pivot_index >= ending_triangle - 1 || ending_triangle-starting_triangle < 5) {
        return;
    }

    node->child_left = new Node();
    node->child_right = new Node();
    this->buildBVH(node->child_left, starting_triangle, pivot_index);
    this->buildBVH(node->child_right, pivot_index, ending_triangle);
}

    
    Intersection intersect(const Ray &ray) override {

        Intersection I(this->albedo, this->reflective);
        double t;
        double t_min{DBL_MAX};
        if (!this->root->bbox.intersect(ray, &t))
            return Intersection();

        std::list<Node*> nodes_to_visit;
        nodes_to_visit.push_front(root);
        while (!nodes_to_visit.empty()) {

            Node *curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();

            if (curNode->child_left) {
                if (curNode->child_left->bbox.intersect(ray,  &t))
                    if (t < t_min)
                        nodes_to_visit.push_back(curNode->child_left);
                if (curNode->child_right->bbox.intersect(ray, &t))
                    if (t < t_min)
                        nodes_to_visit.push_back(curNode->child_right);
            } else {
                Vector A, B, C, e1, e2, N;
                for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                    TriangleIndices triangle_indices = this->indices[i];
                    A = scaling_factor * vertices[triangle_indices.vtxi] + translation;
                    B = scaling_factor * vertices[triangle_indices.vtxj] + translation;
                    C = scaling_factor * vertices[triangle_indices.vtxk] + translation;
                    e1 = B - A;
                    e2 = C - A;
                    N = cross(e1, e2);

                    double beta = dot(cross(A - ray.O, ray.u), e2) / dot(ray.u, N);
                    double gamma = - dot(cross(A - ray.O, ray.u), e1) / dot(ray.u, N);
                    double alpha = 1. - beta - gamma;
                    double t = dot(A - ray.O, N) / dot(ray.u, N);

                    if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < t_min) {
                        t_min = t;
                        I.intersects = true;
                        I.t = t;
                        I.P = A + beta * e1 + gamma * e2;
                        I.N = N;
                    }
                }
            }
        }
        return I;
    }


    Intersection intersect_withMuller(const Ray &ray)  {
        int ID=-1;
        Intersection intersection;
		intersection.t = 1E9;
        intersection.intersects = false;
        for (int i =0; i < indices.size(); i++){
            Intersection intersection2;
            double alpha,beta,gamma,t;
            Vector N;
            intersection2.intersects =  Möller_Trumbore(ray,vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], alpha,beta,gamma,t,N);
            intersection2.P= ray.O + t * ray.u;
            intersection2.N= N;
            intersection2.t=t;
            if (intersection2.intersects){
				if (intersection2.t < intersection.t){
					ID = i;
					intersection.t = intersection2.t;
					intersection.P = intersection2.P;
					intersection.N = intersection2.N;
					intersection.intersects  = true;
				}
			}

        }
        return intersection;

    }



};

/* --------------------------------- SPHERE CLASS ----------------------------------------------------------------------------------------*/

class Sphere : public Geometry {
public:
        Vector C;
        double radius;  
       
       
         Sphere(Vector C, double radius, Vector color, bool reflective = false,  bool transparent = false) {
            this->C = C;
            this->radius = radius;
            this->albedo = color;
            this->reflective = reflective;
            this->transparent = transparent;
        }

        Intersection intersect(const Ray& ray) override {
            Intersection I(this->albedo,  this->reflective);
            Vector origin_center = ray.O - this->C;
            double delta = pow(dot(ray.u, origin_center), 2) - (dot(origin_center, origin_center) - pow(radius, 2));
            if (delta >= 0.) {
            double t1 = dot(ray.u, -1.*origin_center) - sqrt(delta);
            double t2 = dot(ray.u, -1.*origin_center) + sqrt(delta);
            if (t1 > 0) {
                 I.t = t1;
            } else {
            if (t2 > 0) {
                    I.t = t2;
            } else {
                 I.t = 0.0;
                }
            }
            if (t2 < 0.) {
        I.intersects = false;
        } else {
          I.intersects = true;
         }
}
            I.P = ray.O + (I.t * ray.u);
            I.N = (I.P - this->C);
            I.N.normalize();
            return I;
        }
};

/* ----------------------------- SCENE CLASS ----------------------------------------------------------------------------------------- */

class Scene {
       
    public:
        std::vector<Geometry*> geometries;
        Vector L; // position of the light
        double I; // intensity of the light
        explicit Scene(Vector light, double intensity) { L = light; I =intensity; }

        void addGeometry(Geometry* geometry) { geometries.push_back(geometry); }

    Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double cos_theta = sqrt(1 - r2);
    double sin_theta = sqrt(r2);
    double x = cos( 2 * M_PI * r1) * sin_theta;
    double y = sin( 2 * M_PI * r1) * sin_theta;
    double z = cos_theta;
    Vector T1, T2;
    if (fabs(N[0]) <= fabs(N[1]) && fabs(N[0]) <= fabs(N[2])) {
        T1 = Vector(0, N[2], -N[1]);
    } else if (fabs(N[1]) <= fabs(N[0]) && fabs(N[1]) <= fabs(N[2])) {
        T1 = Vector(N[2], 0, -N[0]);
    } else {
        T1 = Vector(N[1], -N[0], 0);
    }
    T2 = cross(N, T1);
    T1.normalize();
    T2.normalize();
    return x * T1 + y * T2 + z * N;
}

Intersection intersect(const Ray& ray, int &ID){
		Intersection intersection;
		intersection.t = 1E9;
		intersection.intersects = false;
		for (int i =0; i < geometries.size(); i++){
			Intersection intersection2;
			intersection2 =geometries[i]->intersect(ray) ;
			if (intersection2.intersects){
				if (intersection2.t < intersection.t){
					ID = i;
					intersection.t = intersection2.t;
					intersection.P = intersection2.P;
					intersection.N = intersection2.N;
					intersection.intersects  = true;
				}
			}
		}

		return intersection;};


Vector getColor(const Ray& ray, int bounce){
	
		Vector color(0,0,0);
		int ID;
		double n1 =1;
		double n2 = 1.5;


		if (bounce < 0){
            return color ;} 
		Intersection inter1= intersect(ray, ID);

	
		if (inter1.intersects){
			double eps = 0.001 ; // 1e-10;
            Vector localP = inter1.P + (eps * inter1.N);
            Vector localN = inter1.N;


			if (geometries[ID]->reflective){
				Ray reflected_ray = Ray(localP, ray.u - (2*dot(ray.u, localN) * localN));
				return getColor(reflected_ray,bounce-1);
			}

			if (geometries[ID]->transparent) {
			
				if (dot(ray.u, localN) > 0){
					std::swap(n1,n2);
					localN = - localN;
					}
			
				if ( 1-sqr(n1/n2) * (1-sqr(dot(ray.u, localN))) <0){
					Ray reflected_ray = Ray(localP, ray.u - (2*dot(inter1.N, ray.u) * inter1.N));
					return getColor(reflected_ray ,bounce-1);}
				
				Vector wt_T = (n1/n2) * (ray.u - dot(ray.u, localN) * localN);
				Vector wt_N = -1. * localN * sqrt(1. - pow(n1/n2, 2) * (1 - pow(dot(ray.u, localN), 2)));
				Vector w_t = wt_T + wt_N;
				Ray refracted(localP + 0.001 * localN, w_t);
    			return getColor(refracted, bounce-1);}	
				
                double d = (L - localP).norm();
                Vector w_i = L - localP;
				w_i.normalize();

				int IDlight;
                Intersection Interlight = intersect(Ray(L, w_i*(-1.)), IDlight); 

                double visibility;
                if (!Interlight.intersects || Interlight.t > d) {
                            visibility = 1.0;
                        } else {
                visibility = 0.0;
                    }

                color = I / (4*M_PI*d*d) * geometries[ID]->albedo / M_PI * visibility * std::max(0., dot(w_i, localN));

                Ray random_ray = Ray(localP, random_cos(localN));
                color = color + geometries[ID]->albedo * getColor(random_ray, bounce - 1);
			}
  			return color;
 }

  			

};

Vector computeraydirection(Vector Camera, int x, int y, double fov, int W, int H) {
    // Convert pixel indices to coordinates 
	Vector u(Camera[0] + y + 0.5 - W/2, - (Camera[1] + x + 0.5 - H/2), Camera[2] -W / (2 * tan( fov / 2)));
	u.normalize();
    return u;
}



int main() {

    
    auto start = std::chrono::high_resolution_clock::now(); 

    int W = 512;
	int H = 512;
	double fov = 60 * M_PI / 180.;
	Vector Q(0,0,55); 
    int max_depth = 5;
    int K = 64;

    Scene scene = Scene(Vector(-10, 20, 40), 2E10); // light position , light intensity 
   
    Sphere S1 = Sphere(Vector(-20, 0, 0), 8, Vector(0.3, 0.2, 0.6), false,  false);
    Sphere S2 = Sphere(Vector(0, 0, 0), 8, Vector(0.3, 0.2, 0.6), true,  false);
    Sphere S3 = Sphere(Vector(20, 0, 0), 8, Vector(0.3, 0.2, 0.6), false,  true);
   
    Sphere ceiling = Sphere(Vector(0,1000,0),  940, Vector(0.7,0.1,0.25));
    Sphere floor = Sphere(Vector(0,-1000,0), 990, Vector(0.7,0.2,0.25));
    Sphere front = Sphere(Vector(0,0,-1000), 940, Vector(0.7,0.3,0.25));
    Sphere back = Sphere(Vector(0,0,1000),  940, Vector(0.7,0.3,0.25));
    Sphere left = Sphere(Vector(-1000,0,0), 940, Vector(0.7,0.3,0.25));
    Sphere right = Sphere(Vector(1000,0,0),  940, Vector(0.7,0.1,0.25));

   
    scene.addGeometry(&ceiling);
    scene.addGeometry(&floor);
    scene.addGeometry(&front);
    scene.addGeometry(&back);
    scene.addGeometry(&left);
    scene.addGeometry(&right);


    TriangleMesh cat = TriangleMesh(0.6, Vector(0, -10, 0), Vector(0.3, 0.2, 0.6));
    cat.readOBJ("cat/cat.obj");
    scene.addGeometry(&cat);
   
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++){
        for (int j = 0; j < W; j++) {
            Vector color = Vector(0., 0., 0.);
            double x, y;
            for (int k = 0; k < K; k++) {
        
                BoxMuller(0.5, x, y);
				Vector u = computeraydirection(Q,i+y,j+x,  fov,  W,  H);
                Ray ray(Q,u);
                color = color + scene.getColor(ray,5); 
            }

            image[(i * W + j) * 3 + 0] = std::min(255., pow(color[0]/K, 0.45));
            image[(i * W + j) * 3 + 1] = std::min(255., pow(color[1]/K, 0.45));
            image[(i * W + j) * 3 + 2] = std::min(255., pow(color[2]/K, 0.45));
        }
        }

    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    std::cout << "Time for Rendering " << duration.count() << " milliseconds." << std::endl;
    
    return 0;
}

