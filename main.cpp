#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

long long amtRaysDrawn = 0;
int windowWidth = 800;
int windowHeight = 600;
int maxDepth = 4;
float *framebuffer;
GLFWwindow *pWindow;
glm::vec3 eyept(0.0f, 0.0f, 1000.0f);
int scene = -1;
int aa = -1;

// FUNCTIONS
glm::vec3 normalized(glm::vec3 v) {
	float denom = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
	return v / denom;
}
glm::vec4 normalized(glm::vec4 v) {
	float denom = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z) + (v.w * v.w));
	return v / denom;
}
glm::vec3 reflected(glm::vec3 v, glm::vec3 normal) {
	return v - (2.0f * glm::dot(normal, v) * normal);
}
bool quadraticFormula(const float &a, const float &b, const float &c, float &t0, float &t1) {
	float discrem = (b * b) - (4 * a * c);
	if (discrem < 0) {
		//std::cout << "h" << std::endl;
		return false;
	}
	else if (discrem == 0) {
		t0 = t1 = -0.5f * (b / a);
	}
	else {
		float q = (b > 0) ? 0.5f * (-b + sqrt(discrem)) : 0.5f * (-b - sqrt(discrem));
		t0 = q / a;
		t1 = c / q;
	}
	return true;
}
template <typename T> T max(T a, T b) { return a >= b ? a : b; }
bool hasValid(std::vector<float> nums) { for (auto& n : nums) { if (n > 0) { return true; } } return false; }
bool shadowRay(glm::vec3& origin, glm::vec3 dir);
float clamp(float min, float max, float n) { if (n > max) { return max; } else if (n < min) { return min; } else { return n; } }

// CLASSES

struct Ray {
	glm::vec3 p;					// starting point
	glm::vec3 v;					// direction
	bool isShadowRay;
	Ray() {}
	Ray(glm::vec3 pp, glm::vec3 vv) : p(pp), v(vv) {}
	glm::vec3 evaluate(float t) {
		glm::vec3 result;
		result.x = p.x + v.x * t;
		result.y = p.y + v.y * t;
		result.z = p.z + v.z * t;
		return result;
	}
};

class ImagePlane {
	private:
		double width;
		double height;
		glm::vec3 center;
		int curX;
		int curY;
	public:
		ImagePlane() : width(windowWidth), height(windowHeight), center(glm::vec3(0.0f, 0.0f, 0.0f)), curX(0), curY(0) {}
		glm::vec3 nextpixel() {
			double xdist = width * ((curX+1) / width - 0.5);
			double ydist = height * ((curY+1) / height - 0.5);
			++curX;
			if (curX == width) {
				curX = 0;
				++curY;
			}
			return glm::vec3(xdist, ydist, 0.0);
		}
};

struct Light {
	float Ia;		// ambient intensity
	float Id;		// diffuse intensity
	float Is;		// specular intensity
	glm::vec3 pos;
	Light(float a, float d, float s, glm::vec3 p)
		: Ia(a), Id(d), Is(s), pos(p) {}
};

std::vector< Light > lights;

class Object {
	private:
		glm::vec3 ka;				// ambient component
		glm::vec3 kd;				// diffuse component
		glm::vec3 ks;				// specular component
		glm::vec3 k_reflection;		// reflection component
		glm::vec3 k_refraction;		// refraction component
		glm::vec3 basecolor;		// color of object
		float ns;					// specularity exponent
		float ior;					// index of refraction
		bool hasTexture;
	public:
		// SETTERS
		void set_ka(glm::vec3 kaa) { ka = kaa; }
		void set_kd(glm::vec3 kdd) { kd = kdd; }
		void set_ks(glm::vec3 kss) { ks = kss; }
		void set_ns(float nss) { ns = nss; }
		void set_material(glm::vec3 kaa, glm::vec3 kdd, glm::vec3 kss, float nss) {
			ka = kaa;
			kd = kdd;
			ks = kss;
			ns = nss;
		}
		void set_reflection(glm::vec3 reflect) { k_reflection = reflect; }
		void set_refraction(glm::vec3 refract, float iorr) { k_refraction = refract; ior = iorr; }
		void set_basecolor(glm::vec3 c) { basecolor = c; }
		// GETTERS
		virtual glm::vec3 getNormal(glm::vec3 p) = 0;
		Ray getReflectionRay(glm::vec3 direction, glm::vec3 position) {
			float c1 = -1.0f * glm::dot(position, direction);
			glm::vec3 r1 = direction + (2.0f * position * c1);

			return Ray(position, r1);
		}
		Ray getRefractionRay(glm::vec3 direction, glm::vec3 position, bool outside) {
			float cosI = clamp(-1.0f, 1.0f, glm::dot(direction, position));
			float currentIOR = 1;
			float nextIOR = ior;
			if (cosI < 0) { cosI = -1.0f * cosI; } else { std::swap(currentIOR, nextIOR); position = -1.0f * position; }
			float rateIOR = currentIOR / nextIOR;
			float k = 1 - rateIOR * rateIOR * (1 - cosI * cosI);
			glm::vec3 res = (rateIOR * (direction + (rateIOR * cosI - sqrtf(k)) * position));
			return Ray(outside ? position + (.001f * direction) : position + (-.001f * direction), k < 0 ? glm::vec3(0.0f) : res);
		}
		glm::vec3 getBaseColor() { return basecolor; }
		glm::vec3 getReflectivity() { return k_reflection; }
		glm::vec3 getRefractivity() { return k_refraction; }

		virtual std::vector<float> intersect(const Ray& ray) = 0;
		glm::vec3 getColor(glm::vec3 position) {

			glm::vec3 color(0.0f, 0.0f, 0.0f);
			glm::vec3 pnormal = getNormal(position);

			// calculate ambient component
			glm::vec3 ambient = lights[0].Ia * ka;
			for (auto& light : lights) {
				glm::vec3 direction = normalized(light.pos - position);
				if (!shadowRay(position, light.pos)) {
					// calculate diffuse component
					glm::vec3 diffuse = light.Id * kd * max(glm::dot(pnormal, direction), 0.0f);

					// calculate specular component
					glm::vec3 viewDir = normalized(eyept - position);
					glm::vec3 reflectDir = reflected(-1.f * direction, pnormal);
					glm::vec3 specular(0.0f, 0.0f, 0.0f);
					if (glm::dot(viewDir, reflectDir) > 0) {
						float spec = pow(max(glm::dot(viewDir, reflectDir), 0.0f), ns);
						specular = light.Is * ks * spec;
					}
					glm::vec3 result(0.0f, 0.0f, 0.0f);
					result = (diffuse)* basecolor;
					color += result;
					color += specular;
				}
			}
			color += ambient * basecolor;

			return color;
		}
		bool isReflective() {
			return (k_reflection.x > 0.0) && (k_reflection.y > 0.0) && (k_reflection.x > 0.0);
		}
		bool isRefractive() {
			return (k_refraction.x > 0.0) && (k_refraction.y > 0.0) && (k_refraction.x > 0.0);
		}
};

std::vector< Object * > objs;

class Plane : public Object
{
	private:
		glm::vec4 equation;			// Ax + By + Cz + D = 0
		glm::vec3 pt;
	public:
		Plane(glm::vec4 e) {
			equation = e;
		}
		Plane(glm::vec3 normal, glm::vec3 pt) {
				this->pt = pt;
			float D = -1.0f * ((normal.x * pt.x) + (normal.y * pt.y) + (normal.z * pt.z));
			equation.x = normal.x;
			equation.y = normal.y;
			equation.z = normal.z;
			equation.w = D;
		}

		// SETTERS
		void setEquation(glm::vec4 eq) { equation = eq; }

		// VIRTUAL FUNCTIONS

		// for a plane with equation of form Ax + By + Cz + D = 0, the
		// normal is defined as (A, B, C)
		virtual glm::vec3 getNormal(glm::vec3 p = { 0, 0, 0 }) {
			return normalized(glm::vec3(equation.x, equation.y, equation.z));
		}
		virtual std::vector<float> intersect(const Ray& ray) {
			std::vector<float> result;
			float denom = (equation.x * ray.v.x) + (equation.y * ray.v.y) + (equation.z * ray.v.z);
			float numer = glm::dot(glm::vec3(equation.x, equation.y, equation.z), ray.p) + equation.w;
			if (denom != 0.0f) {
				result.push_back(-1.0f * (numer / denom));
			}
			return result;
		}


};

class Sphere : public Object
{
	private:
		glm::vec3 center;
		float radius;
	public:
		Sphere(glm::vec3 c, float r) : center(c), radius(r) {}
		// SETTERS
		void setCenter(glm::vec3 c) { center = c; }
		void setRadius(float rad) { radius = rad; }

		// VIRTUAL FUNCTIONS

		// get vector from center of sphere to the point you're calculating normal of,
		// divide by the radius of the sphere to normalize the vector
		virtual glm::vec3 getNormal(glm::vec3 p) {
			return glm::vec3(p.x - center.x, p.y - center.y, p.z - center.z) / radius;
		}
		// will return two t values, eventually find the min positive
		virtual std::vector<float> intersect(const Ray& ray) {
			float t0, t1; std::vector<float> results;
			glm::vec3 l = ray.p - center;
			float a = glm::dot(ray.v, ray.v);
			float b = 2.0f * glm::dot(l, ray.v);
			float c = glm::dot(l, l) - (radius * radius);
			if (quadraticFormula(a, b, c, t0, t1)) {
				results.push_back(t0);
				results.push_back(t1);
			}
			return results;
		}
};

// Draws the scene
void drawit(void)
{
	glDrawPixels(windowWidth, windowHeight, GL_RGB, GL_FLOAT, framebuffer);
}

// Clears framebuffer to black
void clearFramebuffer()
{
	memset(framebuffer, 0, windowWidth * windowHeight * sizeof(float));
}

// Origin is lower left corner, x positive is right side, y positive is upwards.
void setFramebuffer(int x, int y, float R, float G, float B)
{
	framebuffer[(y * windowWidth + x) * 3] = R;
	framebuffer[(y * windowWidth + x) * 3 + 1] = G;
	framebuffer[(y * windowWidth + x) * 3 + 2] = B;
}

void setFramebuffer(int x, int y, glm::vec3 color)
{
	framebuffer[(y * windowWidth + x) * 3] = color.r;
	framebuffer[(y * windowWidth + x) * 3 + 1] = color.g;
	framebuffer[(y * windowWidth + x) * 3 + 2] = color.b;
}

void sceneSetup(int scene)
{
	if (scene == 0) {
		glm::vec3 gray(0.7f, 0.7f, 0.7f);
		Plane* p1 = new Plane(glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 150.0f, 0.0f));
		p1->set_basecolor(gray);
		p1->set_material(glm::vec3(0.7f, 0.7f, 0.7f), glm::vec3(0.5f, 0.5f, 0.5f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		Plane* p2 = new Plane(glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(150.0f, 0.0f, 0.0f));
		p2->set_basecolor(glm::vec3(0.0f, 1.0f, 0.0f));
		p2->set_material(glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		Plane* p5 = new Plane(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(-150.0f, 0.0f, 0.0f));
		p5->set_basecolor(glm::vec3(1.0f, 0.0f, 0.0f));
		p5->set_material(glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		Plane* p3 = new Plane(glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 50.0f));
		p3->set_basecolor(gray);
		p3->set_material(glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		Plane* p4 = new Plane(glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, -150.0f, 0.0f));
		p4->set_basecolor(gray);
		p4->set_material(glm::vec3(0.6f, 0.6f, 0.6f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);

		Sphere* s1 = new Sphere(glm::vec3(0.0f, -1.0f, 200.0f), 50.0f);
		s1->set_basecolor(glm::vec3(0.0f, 0.0f, 1.0f));
		s1->set_material(glm::vec3(0.4f, 0.4f, 0.4f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		s1->set_reflection(glm::vec3(0.5f, 0.5f, 0.5f));
		s1->set_refraction(glm::vec3(0.5f, 0.5f, 0.5f), 0.5f);
		Sphere* s2 = new Sphere(glm::vec3(50.0f, -1.0f, 400.0f), 50.0f);
		s2->set_basecolor(gray);
		s2->set_material(glm::vec3(0.4f, 0.4f, 0.4f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
		s2->set_reflection(glm::vec3(0.5f, 0.5f, 0.5f));
		s2->set_refraction(glm::vec3(0.5f, 0.5f, 0.5f), 1.31f);

		objs.push_back(p1);
		objs.push_back(p2);
		objs.push_back(p3);
		objs.push_back(p4);
		objs.push_back(p5);
		objs.push_back(s1);
		objs.push_back(s2);

		// add lights
		Light pointlight1(0.51f, 0.51f, 0.51f, glm::vec3(0.0f, 130.f, 300.0f));
		Light pointlight2(0.1f, 0.15f, 0.11f, glm::vec3(20.0f, 120.0f, 300.0f));

		lights.push_back(pointlight1);
		lights.push_back(pointlight2);
	}
	else if (scene == 1) {
		float y = 0.0f;
		for (int i = 0; i < 5; ++i) {
			for (int j = 0; j < 10; ++j) {
				Sphere * s = new Sphere(glm::vec3(j * 20.0f, y, 200.0f), 5.0f);
				s->set_basecolor(glm::vec3(0.0f, 1.0f, 0.0f));
				s->set_material(glm::vec3(0.4f, 0.4f, 0.4f), glm::vec3(0.8f, 0.8f, 0.8f), glm::vec3(0.2f, 0.2f, 0.2f), 32.f);
				objs.push_back(s);
			}
			y += 25.0f;
		} 
		Light pointlight1(0.4f, 0.4f, 0.4f, glm::vec3(0.0f, 150.0f, 400.0f));
		lights.push_back(pointlight1);
	}
}

float getLowestTval(std::vector<float> tvals) {
	float mint = 10000.0f;
	for (auto& t : tvals) {
		if (t < mint && t > 0) {
			mint = t;
		}
	}
	return mint;
}

bool trace(float& mintval, int& minobj, Ray& ray) {
	for (int i = 0; i < objs.size(); ++i) {
		std::vector<float> tvals = objs[i]->intersect(ray);
		if (tvals.size() > 0) {
			if (getLowestTval(tvals) < mintval) {
				if (getLowestTval(tvals) > 0.0f && ( (ray.isShadowRay && getLowestTval(tvals) < 1.0f) || !ray.isShadowRay )) {
					mintval = getLowestTval(tvals);
					minobj = i;
				}
			}
		}
	}
	return minobj != -1;
}

bool shadowRay(glm::vec3& origin, glm::vec3 dir) {
	float mintval = 100000.0f;
	Ray curRay(origin, dir - origin);
	origin = curRay.evaluate(.001f);
	curRay.p = origin;
	curRay.isShadowRay = true;
	int minobj = -1;
	if (trace(mintval, minobj, curRay)) {
		return true;
	}
	else {
		return false;
	}
}

glm::vec3 castRay(Ray ray, const int& depth) {
	if (depth > maxDepth) {
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}
	++amtRaysDrawn;
	glm::vec3 hitColor(0.0f, 0.0f, 0.0f);
	int N = 128;
	float mintval = 1000000.0f;
	int minobj = -1;
	if (trace(mintval, minobj, ray)) {
		glm::vec3 phit = ray.evaluate(mintval);

		// compute direct lighting
		glm::vec3 color = objs[minobj]->getColor(phit);

		// compute indirect lighting
		glm::vec3 indirectColor(0.0f, 0.0f, 0.0f);
		glm::vec3 reflectColor(0.0f);
		glm::vec3 refractColor(0.0f);
		if (objs[minobj]->isReflective()) {
			Ray reflectedRay = objs[minobj]->getReflectionRay(ray.v, phit + (.001f * (-1.f * ray.v)));
			reflectColor += castRay(reflectedRay, depth + 1);
		}
		 
		if (objs[minobj]->isRefractive()) {
			bool outside = glm::dot(ray.v, phit) < 0;
			Ray refractedRay = objs[minobj]->getRefractionRay(ray.v, phit, outside);
			refractColor += castRay(refractedRay, depth + 1);
		}

		// combine color
		hitColor = color + (reflectColor * objs[minobj]->getReflectivity()) + (refractColor * objs[minobj]->getRefractivity());
	}
	else {
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}
	return hitColor;

}
 
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	// set up scene's contents
	objs.clear();
	sceneSetup(scene);

	// create image plane
	ImagePlane ip;
	auto start = std::chrono::high_resolution_clock::now();

	// RENDER PIXELS
	for (int i = 0; i < windowHeight; ++i) {
		for (int j = 0; j < windowWidth; ++j) {
			glm::vec3 pixelColor(0.0f);
			if (aa == 1) {
				for (int samples = 0; samples < 4; ++samples) {
					float jitterMatrix[4 * 2] = {
						-1.0 / 4.0,  3.0 / 4.0,
						3.0 / 4.0,  1.0 / 3.0,
						-3.0 / 4.0, -1.0 / 4.0,
						1.0 / 4.0, -3.0 / 4.0,
					};
					double xdist = ((double)windowWidth) * ((j + jitterMatrix[2 * samples]) / ((double)windowWidth) - 0.5);
					double ydist = ((double)windowHeight) * ((i + jitterMatrix[2 * samples + 1]) / ((double)windowHeight) - 0.5);
					Ray primaryRay(eyept, glm::vec3(xdist, ydist, 0.0f) - eyept);
					pixelColor += castRay(primaryRay, 0);
				}
				setFramebuffer(j, i, (pixelColor / 4.0f));
			}
			else {
				glm::vec3 currpixel = ip.nextpixel();
				Ray primaryRay(eyept, currpixel - eyept);
				pixelColor += castRay(primaryRay, 0);
				setFramebuffer(j, i, pixelColor);
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "Done. Took " << duration.count() * 1e-6 << " seconds to render. " << std::endl;
	std::cout << "Amount of rays drawn: " << amtRaysDrawn << std::endl;
}

void KeyCallbackFunc(GLFWwindow* pwnd, int key, int scancode, int action, int mode)
{
}

void init()
{
	framebuffer = new float[3 * windowWidth * windowHeight];
	glfwInit();
	glfwSetTime(0.0);
	if (pWindow)
	{
		glfwDestroyWindow(pWindow);
	}
	pWindow = glfwCreateWindow(windowWidth, windowHeight, "Assignment 6 - Aksel Taylan", NULL, NULL);
	glfwMakeContextCurrent(pWindow);
	glfwSetKeyCallback(pWindow, KeyCallbackFunc);
	glewExperimental = true;
	glewInit();
	glViewport(0, 0, windowWidth, windowHeight);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	clearFramebuffer();
}

int main()
{
	std::cout << "First enter which scene, 0 or 1. Next, determine if you want anti-aliasing (1 is yes). " << std::endl;
	std::cin >> scene >> aa;

	init();
	sceneSetup(scene);

	display();
	while (!glfwWindowShouldClose(pWindow)) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		drawit();
		glFlush();
		glfwSwapBuffers(pWindow);
		glfwPollEvents();
	}
	delete[] framebuffer;
	return 0;
}