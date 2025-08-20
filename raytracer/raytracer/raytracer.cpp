#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

#if defined __linux__ || defined __APPLE__
//These values are already present on Linux and Apple
#else
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0)
        {
            T invNor = 1 /sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    //Operator overloads
    Vec3<T> operator * (const T &f) const {return Vec3<T>(x * f, y * f, z * f);}
    Vec3<T> operator * (const Vec3<T> &v) const {return Vec3<T>(x * v.x, y * v.y, z * v.z);}
    T dot(const Vec3<T> &v) const {return x * v.x+y * v.y+z * v.z;}
    Vec3<T> operator - (const Vec3<T> &v) const {return Vec3<T>(x - v.x, y - v.y, z - v.z);}
    Vec3<T> operator + (const Vec3<T> &v) const {return Vec3<T>(x + v.x, y + v.y, z + v.z);}
    Vec3<T>& operator += (const Vec3<T> &v){x += v.x, y += v.y, z += v.z; return *this;}
    Vec3<T>& operator *= (const Vec3<T> &v){x *= v.x, y *= v.y, z *= v.z; return *this;}
    Vec3<T> operator - () const {return Vec3<T>(-x, -y, -z);}
    T length2() const {return x*x + y*y + z*z;}
    T length() const { return sqrt(length2());}
    friend std::ostream& operator << (std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};
typedef Vec3<float> Vec3f;
class Sphere
{
    public:
    Vec3f center;                       //Also the pivot point of the sphere
    float radius, radius2;
    Vec3f surfaceColor, emissionColor;
    float transparency, reflection;

    //[comment]
    // Note that all values are passed by reference for efficiency (avoids copying) but introduces the risk of modifying
    // the value, this is why const is used.
    // May be unnecessary in the case of floats though.
    //[/comment]
    Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vec3f &ec = 0
    ) : center(c), radius (r), radius2(r * r), surfaceColor(sc), emissionColor(ec), transparency(transp), reflection(refl){}

    //[comment]
    // Compute a ray-sphere intersection using the geometric solution
    //[/comment]
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false; //if this is true then it's "behind the ray origin"
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false; //the ray never meets the sphere
        float thc = sqrt(radius2 - d2);
        //resulting intersection points
        t0 = tca - thc; //entry point
        t1 = tca + thc; //exit point
        
        return true;
    }
};

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1-mix);
}

Vec3f trace(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Sphere> &spheres,
    const int &depth)
{
    float tnear = INFINITY;
    const Sphere* sphere = NULL; //the sphere we are considering for rendering (?)
    //Iterate on the spheres. Find the intersections (hits) with each of them.
    for (unsigned i = 0; i < spheres.size(); ++i){
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayorig, raydir, t0, t1)){
            if (t0<0) t0 = t1;
            if (t0 < tnear){
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    
    if (!sphere) return Vec3f(2); //If there is no intersection you return black OR the BG color
    Vec3f surfaceColor = 0;
    Vec3f phit = rayorig + raydir * tnear;
    Vec3f nhit = phit - sphere->center;
    nhit.normalize(); //evaluate the normal for the hit.
    float bias = 1e-4;
    bool inside = false;
    if (raydir.dot(nhit)>0) nhit = -nhit, inside = true; //we are inside the sphere (?)

    //Sphere is reflactive and/or transparent (can be reflective but not transparent but not transparent
    // and not reflective) we check first for reflections.
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -raydir.dot(nhit);
        float fresneleffect = mix(pow(1-facingratio,3), 1, 0.1); //evaluate refraction of light with fresnel effect
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        refldir.normalize();
        Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
        Vec3f refraction = 0;
        //Sphere is also transparent. Compute refraction ray (transmission)
        if (sphere->transparency) {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }
        //Mix reflection and refraction for the result
        surfaceColor = (
            reflection * fresneleffect +
            refraction * (1-fresneleffect)*sphere->transparency) * sphere->surfaceColor;
    }else //Sphere is opaque (or depth exceedes MAX_RAY_DEPTH)
    {
        //Cut down some calculations by not raytracing recursively
        for (unsigned i = 0; i < spheres.size(); ++i){
            if (spheres[i].emissionColor.x > 0){
                //This is a light
                Vec3f transmission = 1;
                Vec3f lightDirection = spheres[i].center - phit;
                lightDirection.normalize();
                for (unsigned j = 0; j < spheres.size(); ++j){
                    if (i!=j){
                        float t0, t1;
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)){
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += sphere->surfaceColor * transmission *
                    std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
            }
        }
    }
    return surfaceColor + sphere->emissionColor;
}


/// Main rendering function. We compute a camera ray for each pixel of the image
/// trace it and return a color. If the ray hits a sphere, we return the color of the
/// sphere at the intersection point, else we return the background color.
/// @param spheres - the spheres to render.
void render(const std::vector<Sphere> &spheres)
{
    float width = 640.0f, height = 480.0f;
    Vec3f *image = new Vec3f[width*height], *pixel = image;
    float invWidth = 1.0f/float(width), invHeight = 1.0f/float(height);
    float fov = 30, aspectratio = width/float(height);
    float angle = tan(M_PI*0.5*fov/180.0);
    //trace the rays for each pixel that has to be on screen
    for (unsigned y = 0; y < height; ++y){
        for (unsigned x = 0; x < width; ++x, ++pixel){
            //Mapping coordinates to camera space
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            //Building the ray direction
            Vec3f raydir(xx,yy,-1);
            raydir.normalize();
            *pixel = trace(Vec3f(0),raydir, spheres, 0);
        }
    }
    //Savign the result to a PPM image
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;
}
    
//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
    std::vector<Sphere> spheres; //Store the spheres in a vector
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0,     -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0)); //The plane holding the spheres
    spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
    // light
    spheres.push_back(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    render(spheres);
    
    return 0;
}