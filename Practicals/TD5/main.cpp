#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <algorithm>
#include <random>
#include <vector>
#include <iostream>

std::default_random_engine generator;
std::normal_distribution<double> N01(0., 1.);



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
    Vector normalize() {
        double n = norm();
        Vector temp;
        temp[0] = data[0]/n;
        temp[1] = data[1]/n;
        temp[2] = data[2]/n;
        return temp;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

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


Vector random_direction() {
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    double r1 = dis(gen);
    double r2 = dis(gen);
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    
    return Vector(x, y, z);
}


void color_matching(int W, int H,   unsigned char* I, unsigned char* M) {
	

    std::vector<std::pair<int, int>> projI(W*H);
    std::vector<std::pair<int, int>> projM(W*H);

    
    Vector image_pixel, color_pixel, v;


	for (size_t k = 0; k < 50; k++) {
		v = random_direction();

		for (size_t i = 0; i < W*H; i++) {
			unsigned char *Out_I = I + i * 3;
			unsigned char *Out_M = M + i * 3;
			image_pixel = Vector(Out_I[0], Out_I[1], Out_I[2]);
			color_pixel = Vector(Out_M[0], Out_M[1], Out_M[2]);
			projI[i] = std::make_pair(dot(image_pixel, v), i);
			projM[i] = std::make_pair(dot(color_pixel, v), i);}
			

        std::sort(projI.begin(), projI.end());
        std::sort(projM.begin(), projM.end());


		for (size_t i = 0; i < W*H; i++) {
			unsigned char *Out = I + projI[i].second * 3;
			image_pixel = Vector(Out[0], Out[1], Out[2]) + (projM[i].first - projI[i].first) * v;
			for (size_t j = 0; j < 3; j++){
				Out[j] = image_pixel[j];}
		
	}
}}


int main() {

	int I_W, I_H, I_C;
	int M_W, M_H, M_C;
	
	unsigned char *I = stbi_load("8733654151_b9422bb2ec_k.jpg", &I_W, &I_H, &I_C,STBI_rgb);
	unsigned char *M = stbi_load("redim.jpg", &M_W, &M_H, &M_C ,STBI_rgb);


   

	color_matching(I_W, I_H,  I, M);

	stbi_write_png("imgA.jpg", I_W, I_H, 3, &I[0], 0);

	std::cout << " Done"   << std::endl;

	return 0;
}

