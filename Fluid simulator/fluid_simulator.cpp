#include <string>
#include "vector.cpp"
#include <chrono>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <iostream>

#include "lbfgs/lbfgs.h"
#include "lbfgs/lbfgs.c"


// ---------------------------------- POLYGON CLASS ----------------------------------------------------------------------------------------------------------------
// This part was implemented following the live coding in TD 
// Some changes were made afterwards with the help of Jasmine Watissee and Victor Lasocki 

class Polygon {  
public:
	double area() {
        int N = vertices.size() ; 
		if (N<= 2) {
            return 0;}
		double area = 0;
		for (int i = 0; i < N; i++){
			area += vertices[i][0] * vertices[(i+1)%N][1] - vertices[(i+1)%N][0] * vertices[i][1];}
		return std::abs(area / 2) ;
	}


    // Formula at the end of page 103 of the lecture notes
	double discreteintegral(const Vector& P) {
        int N = vertices.size() ; 

        // these base cases were seen in class
		if (N <= 2) return 0;
		if (area() == 0) return 0;

		double result = 0;
		for (int t = 1; t < N-1; t++){
            // power cell is triangualted
			Vector c[3] = {vertices[0], vertices[t], vertices[t+1]};

			double sumdotproduct = 0;
			for (int k = 0; k < 3; k++){
				for (int l= k; l< 3; l++){
					sumdotproduct += dot(c[k]-P, c[l]-P);
				}
			}
			result += sumdotproduct *std::abs((c[2][1]-c[0][1]) * (c[1][0]-c[0][0]) - (c[2][0]-c[0][0]) * (c[1][1]-c[0][1])) * .5 / 6 ;
		}

		return result;
	}

   

   Vector centroid(){
		int N = vertices.size();
		// Got told by Jasmine Watissee to add a base case
		if (N <= 2) return Vector(0, 0, 0);
		

		// this follows the known formula for a centroid
		// https://mathworld.wolfram.com/PolygonCentroid.html
		Vector centroid;

		for (int i = 0; i < N; i++){
			centroid = centroid + Vector((vertices[i][0] + vertices[(i + 1) % N][0]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]), (vertices[i][1] + vertices[(i + 1) % N][1]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]));
		}
		return centroid / (6* area());
	}




	

 

	std::vector<Vector> vertices;
};	

// clipping one 
// Follows p 88 of the lecture notes
Polygon clipping(const Polygon &cell,   const Vector& P0, const double& weight0, const Vector& Pi,const double& weighti ) {
	Polygon clippedpolygon;
	int vertexCount = cell.vertices.size();
	for (int i=0; i< vertexCount; i++) {
		const Vector& previousVertex = cell.vertices[((i - 1) >= 0) ? (i - 1) : (cell.vertices.size() - 1)];
		const Vector& currentVertex = cell.vertices[i];
		if ((currentVertex - P0).norm2() - weight0 <= (currentVertex - Pi).norm2()-weighti) { 
			if ((previousVertex - P0).norm2() -weight0 > (previousVertex - Pi).norm2() - weighti) { 
				Vector midPoint = (P0 + Pi) * 0.5;
				midPoint= midPoint + (weight0-weighti) * (Pi-P0)/ (2* (P0-Pi).norm2());
				double param = dot(midPoint - previousVertex, Pi - P0) / dot(currentVertex - previousVertex, Pi - P0);
				Vector intersection = previousVertex + param * (currentVertex - previousVertex);
				clippedpolygon.vertices.push_back(intersection);
			}
			clippedpolygon.vertices.push_back(currentVertex);
		} else {
			if ((previousVertex - P0).norm2() - weight0 <= (previousVertex - Pi).norm2()- weighti) { 
				Vector midPoint = (P0 + Pi) * 0.5;
				midPoint= midPoint + (weight0-weighti) * (Pi-P0)/ (2* (P0-Pi).norm2());
				double param = dot(midPoint - previousVertex, Pi - P0) / dot(currentVertex - previousVertex, Pi - P0);
				Vector intersection = previousVertex + param * (currentVertex - previousVertex);
				clippedpolygon.vertices.push_back(intersection);
			}
		}
	}
	return clippedpolygon;
};



// clipping for disks 
Polygon clipping_disk(const Polygon &cell, const Vector &u, const Vector &v) {
	Vector normal(v[1]-u[1], u[0]-v[0], 0);
	Polygon result;
	int N = cell.vertices.size();

	for (int i=0; i< N; i++){
		const Vector& A = cell.vertices[i==0?(N-1): i-1];
		const Vector& B = cell.vertices[i];
		Vector P = A + (dot(u-A, normal)/dot(B-A, normal)) * (B-A);

		if (dot(B-u, normal) <= 0) { 
			if (dot(A-u, normal) > 0){ 
				result.vertices.push_back(P);
			}
			result.vertices.push_back(B);
		} else if (dot(A-u, normal) <= 0) { 
			result.vertices.push_back(P);
		}
	}

	return result;
}



// ---------------------------------- FUNCTIONS TO SAVE RESULTS -------------------------------------------------------------------------------------------------

void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
		FILE* f = fopen(filename.c_str(), "w+"); 
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
		    fprintf(f, "<g>\n");
		    fprintf(f, "<polygon points = \""); 
		    for (int j = 0; j < polygons[i].vertices.size(); j++) {
			    fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		    }
		    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		    fprintf(f, "</g>\n");
        }
		fprintf(f, "</svg>\n");
		fclose(f);
	}





void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
		FILE* f;
		if (frameid == 0) {
			f = fopen(filename.c_str(), "w+");
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
			fprintf(f, "<g>\n");
		} else {
			f = fopen(filename.c_str(), "a+");
		}
		fprintf(f, "<g>\n");
		for (int i = 0; i < polygons.size(); i++) {
			fprintf(f, "<polygon points = \""); 
			for (int j = 0; j < polygons[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
		}
		fprintf(f, "<animate\n");
		fprintf(f, "	id = \"frame%u\"\n", frameid);
		fprintf(f, "	attributeName = \"display\"\n");
		fprintf(f, "	values = \"");
		for (int j = 0; j < nbframes; j++) {
			if (frameid == j) {
				fprintf(f, "inline");
			} else {
				fprintf(f, "none");
			}
			fprintf(f, ";");
		}
		fprintf(f, "none\"\n	keyTimes = \"");
		for (int j = 0; j < nbframes; j++) {
			fprintf(f, "%2.3f", j / (double)(nbframes));
			fprintf(f, ";");
		}
		fprintf(f, "1\"\n	dur = \"5s\"\n");
		fprintf(f, "	begin = \"0s\"\n");
		fprintf(f, "	repeatCount = \"indefinite\"/>\n");
		fprintf(f, "</g>\n");
		if (frameid == nbframes - 1) {
			fprintf(f, "</g>\n");
			fprintf(f, "</svg>\n");
		}
		fclose(f);
	}


void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
        int W = 1000, H = 1000;
        std::vector<unsigned char> image(W*H * 3, 255);
		#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < cells.size(); i++) {
 
            double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
            for (int j = 0; j < cells[i].vertices.size(); j++) {
                bminx = std::min(bminx, cells[i].vertices[j][0]);
                bminy = std::min(bminy, cells[i].vertices[j][1]);
                bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
                bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
            }
            bminx = std::min(W-1., std::max(0., W * bminx));
            bminy = std::min(H-1., std::max(0., H * bminy));
            bmaxx = std::max(W-1., std::max(0., W * bmaxx));
            bmaxy = std::max(H-1., std::max(0., H * bmaxy));
 
            for (int y = bminy; y < bmaxy; y++) {
                for (int x = bminx; x < bmaxx; x++) {
                    int prevSign = 0;
                    bool isInside = true;
                    double mindistEdge = 1E9;
                    for (int j = 0; j < cells[i].vertices.size(); j++) {
                        double x0 = cells[i].vertices[j][0] * W;
                        double y0 = cells[i].vertices[j][1] * H;
                        double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                        double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                        double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
						int sign;
						if (det>0) sign= 1;
						if 	(det<0) sign=  -1;
						if 	(det==0) sign=  0;

                        if (prevSign == 0) prevSign = sign; else
                            if (sign == 0) sign = prevSign; else
                            if (sign != prevSign) {
                                isInside = false;
                                break;
                            }
                        prevSign = sign;

                        double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                        double distEdge = std::abs(det)/ edgeLen;
                        double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                        if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                        mindistEdge = std::min(mindistEdge, distEdge);
                    }
                    if (isInside) {
              
                          image[((H - y - 1)*W + x) * 3] = 0;
                          image[((H - y - 1)*W + x) * 3 + 1] = 0;
                          image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    
                        if (mindistEdge <= 2) {
                            image[((H - y - 1)*W + x) * 3] = 0;
                            image[((H - y - 1)*W + x) * 3 + 1] = 0;
                            image[((H - y - 1)*W + x) * 3 + 2] = 0;
                        }
 
                    }
                    
                }
            }
        }
        std::ostringstream os;
        os << filename << frameid << ".png";
        stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
    }


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------





// --------------------------------------VORONOI DIAGRAM------------------------------------------------------------------------------------------------------------------------------------

class VoronoiDiagram {
public:
	VoronoiDiagram(){
		for (int i=0; i < sides_disk; i++) {
			disk[i][0] = cos(i * 2 * M_PI / sides_disk);
			disk[i][1] = sin(i * 2 * M_PI / sides_disk);
			disk[i][2] = 0;
			}
	}

	VoronoiDiagram(std::vector<Vector> Points, std::vector<Polygon> Cells, std::vector<double> Weights){
		this->points = Points; 
		this->cells = Cells; 
		this->weights = Weights; 
		int N = cells.size();

		
		for (int i=0; i < sides_disk; i++) {
			disk[i][0] = cos(i * 2 * M_PI / sides_disk);
			disk[i][1] = sin(i * 2 * M_PI / sides_disk);
			disk[i][2] = 0;}
	}



	void compute() {
		int N = points.size();
		Polygon intial_square;
		intial_square.vertices.push_back(Vector(0,0,0)); 
		intial_square.vertices.push_back(Vector(1,0,0));
		intial_square.vertices.push_back(Vector(1,1,0));
		intial_square.vertices.push_back(Vector(0,1,0));

		cells.resize(N);
		#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			Polygon cell = intial_square;
			for (int j = 0; j < N; j++) {
				if (i != j){
				cell = clipping(cell, points[i],  weights[i],points[j], weights[j]);}}
			for (int j=0; j< sides_disk; j++) {
				double radius = sqrt(weights[i] - weights[N]);
				Vector u = disk[j] * radius + points[i];
				Vector v = disk[(j+1)% sides_disk]*radius + points[i];
				cell = clipping_disk(cell, u, v);}
			cells[i] = cell;
		}
	}

    // Below is the compute function before we started visualising fluids
    // I coded this at the same time as the professor in class
    /* void compute() {
		Polygon initialSquare;
		initialSquare.vertices.push_back(Vector(0, 0, 0));
		initialSquare.vertices.push_back(Vector(1, 0, 0));
		initialSquare.vertices.push_back(Vector(1, 1, 0));
		initialSquare.vertices.push_back(Vector(0, 1, 0));

		cells.resize(points.size());
		#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			Polygon cell = initialSquare;
			for (int j = 0; j < points.size(); j++) {
				if (i != j) {
				cell = clipping(cell, points[i], weights[i],  points[j], weights[j]);
				}
			}
			cells[i] = cell;
		}
	}; */

	std::vector<Vector> points;
	std::vector<Polygon> cells;
	std::vector<double> weights;
	static const int sides_disk = 200;
	Vector disk[sides_disk];

};
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// ---------------------------------SEMI DISCRETE OPTIMAL TRANSPORT---------------------------------------------------------------------------------------------------------------------------------------

class SemiDiscreteOT {
public:

	// changes the proportion of air / fluid
	double volume_air = 0.4;
	double volume_fluid = 0.6;
	
	void optimize() {
		int N = diagram.points.size();
		diagram.weights.resize(N+1);

		for (int i=0; i < N; i++) {
			diagram.weights[i] = volume_fluid / N; 
		}
		diagram.weights[N] = volume_air;
		
		double objectivefct = -1;

		// I added this after looking up on the internet why my lbfgs function was not terminating
		// I don't remember the particular website where I found that you can use parameters

		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

		int ret = lbfgs(N+1, &diagram.weights[0], &objectivefct, _evaluate, NULL, this, &param);
		std::cout << "L-BFGS optimization terminated with status code = %d\n" << ret << std::endl;
		
		diagram.compute();
	}
    // Below is the function optimize before we started visualing the fluid

    /*void optimize(){
		int N = diagram.points.size();
		diagram.weights.resize(N);
		
		for (int i=0; i < N; i++){
			diagram.weights[i] = 1/N;
		}		

		double objectivefct = -1;
		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
		
		int ret = lbfgs(N, &updated_weights[0], &objectivefct, _evaluate, NULL, this, &param);
		std::cout << "L-BFGS optimization terminated with status code = %d\n" << ret<< std::endl;

		diagram.compute();
	}*/

	VoronoiDiagram diagram;

protected:

	static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
		diagram.compute();

        lbfgsfloatval_t fx = 0.0;
		
		double lambda = volume_air / (n-1);
		double lambda_air = volume_fluid;
		double sum_area = 0;

        for (int i = 0; i < n-1;i++) {

			sum_area += diagram.cells[i].area();
			fx += diagram.cells[i].discreteintegral(diagram.points[i]) + lambda * diagram.weights[i];
			fx -=  diagram.weights[i] * diagram.cells[i].area();
			g[i] = diagram.cells[i].area() - lambda;
        }

		double estimated_air = 1 - sum_area;
		g[n-1] = estimated_air - lambda_air;
		fx += diagram.weights[n-1] * - g[n-1];
		return -fx;
    }

    // Below is the function evaluate before we started visualing the field
    /*
    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {	
		diagram.compute();
        lbfgsfloatval_t fx = 0.0;
        double lambda = 1./n;
        double sum_area =0;

        for (int i = 0; i < n;i++) {
			sum_area += diagram.cells[i].area();
			fx += diagram.cells[i].discreteintegral(diagram.points[i]) + lambda * diagram.weights[i];
			fx -=  diagram.weights[i] * diagram.cells[i].area();
			g[i] = diagram.cells[i].area() - lambda;
        }
        return -fx;
    }
    */

};



int main(){
	auto start = std::chrono::high_resolution_clock::now();


	int N = 100; // Number of particles 
	double dt =0.01; // Time step
	int steps= 50; // Number of steps
	double mi = 50; // mass of a particle
	double epsilon = 0.004 * 0.004; 
	Vector g(0,-9.81,0);


	SemiDiscreteOT ot;
	ot.diagram.points.resize(N);
	std::vector<Vector> velocity(N, Vector(0, 0, 0)); // velocity of the particles
	std::vector<Vector> positions(N, Vector(0, 0, 0)); // positions of the particles


	// Particles are scattered using randomness
	for (int i=0; i < N; i++) {
		positions[i] = Vector(rand()/double(RAND_MAX), rand()/double(RAND_MAX), 0);
		}
	
	for (int i = 0; i< steps; i++){

		
		ot.diagram.points = positions;
		ot.optimize();

		// Gallouët Mérigot scheme p 122
		for (int i=0; i< N; i++) {
            
			Vector Fi_spring = 1./epsilon*(ot.diagram.cells[i].centroid() - positions[i]);
			Vector Fi= Fi_spring + mi * g;
			velocity[i] =velocity[i] + dt/mi * Fi;
			positions[i] = positions[i] + dt * velocity[i];

		// Jasmine Watissee told me to add this
			positions[i][0] = std::min(1-1E-9, std::max(1E-9, positions[i][0]));
			positions[i][1] = std::min(1-1E-9, std::max(1E-9, positions[i][1]));
		}

		save_frame(ot.diagram.cells, "./fluid_frames/fluid", i);
		
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

	return 0;

}



