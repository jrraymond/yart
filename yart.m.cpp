#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

using dbl = double;
using i32 = int32_t;

constexpr dbl EPSILON = 0.00001;

struct Color
{
    Color() = default;
    Color(dbl r, dbl g, dbl b);
    Color(const Color&) = default;
    Color(Color&&) = default;
    Color& operator=(const Color&) = default;
    Color& operator=(Color&&) = default;
    ~Color() = default;

    Color& operator*=(const Color& rhs);
    Color& operator*=(dbl s);
    Color& operator+=(const Color& rhs);
    bool operator==(const Color& rhs) const;
    bool operator!=(const Color& rhs) const;

    dbl red = 0;
    dbl green = 0;
    dbl blue = 0;
};

Color operator*(Color lhs, const Color& s);
Color operator*(Color lhs, dbl s);
Color operator+(Color lhs, const Color& rhs);
std::ostream& operator<<(std::ostream& out, const Color& c);

struct Vec
{
    Vec() = default;
    Vec(dbl x, dbl y, dbl z);
    Vec(const Vec&) = default;
    Vec(Vec&&) = default;
    Vec& operator=(const Vec&) = default;
    Vec& operator=(Vec&&) = default;
    ~Vec() = default;

    dbl magnitude2() const;
    dbl magnitude() const;
    Vec& normalize();
    
    Vec operator-() const;
    Vec& operator*=(dbl s);
    Vec& operator+=(const Vec& rhs);
    Vec& operator-=(const Vec& rhs);

    bool operator==(const Vec& rhs) const;
    bool operator!=(const Vec& rhs) const;

    dbl x = 0;
    dbl y = 0;
    dbl z = 0;
};

Vec operator*(Vec v, dbl s);
Vec operator+(Vec lhs, const Vec& rhs);
Vec operator-(Vec lhs, const Vec& rhs);
std::ostream& operator<<(std::ostream& out, const Vec& v);

struct Ray
{
  Ray() = default;
  Ray(Vec base, Vec dir);
  Ray(const Ray&) = default;
  Ray(Ray&&) = default;
  Ray& operator=(const Ray&) = default;
  Ray& operator=(Ray&&) = default;
  ~Ray() = default;

  Vec base;
  Vec dir;
};

struct Material
{
  Material(Color c);
  Material(Color c, dbl k_a, dbl k_d);
  Material(Color c, dbl k_a, dbl k_d, dbl k_s, dbl phong);

  dbl k_ambient;
  dbl k_diffuse;
  dbl k_specular;
  dbl phong_exp;
  Color ambient;
  Color diffuse;
  Color specular;
};

class Object;

struct HitRecord
{
  Vec pt;
  Ray ray;
  const Object* object;
};

struct Sphere 
{
    Sphere() = default;
    Sphere(Vec c, dbl radius);
    Sphere(const Sphere&) = default;
    Sphere(Sphere&&) = default;
    Sphere& operator=(const Sphere&) = default;
    Sphere& operator=(Sphere&&) = default;
    ~Sphere() = default;

    bool hits(const Ray& ray, dbl mn, dbl* mx, Vec* pt) const;
    Vec normal(const Vec& pt) const;

    Vec center;
    dbl radius;
};


class Object
{
  public:
    enum type { SPHERE = 0 };

    Object(Sphere s, Material m);

    bool hits(const Ray& ray, dbl mn, dbl* mx, HitRecord* hr) const;
    Vec normal(const Vec& pt) const;

    Material material;

  private:
    type d_tag;
    union Shape {
      Sphere sphere;
      Shape(Sphere&& s) : sphere(s) {}
    };
    Shape d_shape;
};


struct Eye
{
  Eye() = default;
  Eye(Vec pos, const Vec& lookAt, const Vec& up);
  Vec pos;
  Vec u;
  Vec v;
  Vec w;
};

struct Frame
{
  dbl width;
  dbl height;
  dbl distance;
};


struct Light
{
  Vec pos;
  Color color;
};

struct Scene
{
  Eye eye;
  Frame frame;
  i32 img_wd;
  i32 img_ht;
  Color ambient;
  std::vector<Light> lights;

  Ray get_ray(i32 i, i32 j) const;

  Color color(const HitRecord& hr) const;

  template <typename F>
  Color trace(const Ray& ray, F begin, F end) const;
};



// =============================================================================
// MAIN
// =============================================================================
int run_all_tests();

int main(int argc, char** argv)
{
  if (argc > 1)
    return run_all_tests();
  Material blue(   Color(0, 0, 1), 0.4, 0.2, 0.1, 0.0001);
  Material blue2(  Color(0, 0, 1), 0.4, 0.2, 0.1, 0.001);
  Material red(    Color(1, 0, 0), 0.4, 0.2, 0.1, 0.01);
  Material red2(   Color(1, 0, 0), 0.4, 0.2, 0.1, 0.11);
  Material green(  Color(0, 1, 0), 0.4, 0.2, 0.1, 1);
  Material green2( Color(0, 1, 0), 0.4, 0.2, 0.1, 10);
  Material white(  Color(1, 1, 1), 0.4, 0.2, 0.1, 100);
  Material yellow( Color(1, 1, 0), 0.4, 0.2, 0.1, 1000);
  Material yellow2(Color(1, 1, 0), 0.4, 0.2, 0.1, 10000);

  i32 wd = 800;
  i32 ht = 600;

  constexpr i32 NUM_OBJS = 9;
  const std::array<Object, NUM_OBJS> objects = {{
    Object{Sphere{Vec{}, 1}, white},
    Object{Sphere{Vec{2, 0, 0}, 1}, blue},
    Object{Sphere{Vec{4, 0, 0}, 1}, blue2},
    Object{Sphere{Vec{0, 2, 0}, 1}, red},
    Object{Sphere{Vec{0, 4, 0}, 1}, red2},
    Object{Sphere{Vec{-2, 0, 0}, 1}, green},
    Object{Sphere{Vec{-4, 0, 0}, 1}, green2},
    Object{Sphere{Vec{0, -2, 0}, 1}, yellow},
    Object{Sphere{Vec{0, -4, 0}, 1}, yellow2}
  }};
  std::vector<Light> lights = {
    Light{Vec{10, 10, 10}, Color{0.3, 0.3, 0.3}}
  };
  const Eye eye{Vec{0, 0, 10}, Vec{}, Vec{0, 1, 0}};
  const Frame frame{static_cast<dbl>(wd)/100, static_cast<dbl>(ht)/100, 6};
  const Scene scene{eye, frame, wd, ht, Color{1, 1, 1}, std::move(lights)};

  std::ofstream out;
  out.open("out.ppm");
  out << "P3\n" << wd << ' ' << ht << " 255\n";
  for (i32 j = ht - 1; j >= 0; --j) {
    for (i32 i = 0; i < wd; ++i) {
      const Ray ray = scene.get_ray(i, j);
      const Color c = scene.trace(ray, objects.cbegin(), objects.cend());
      out << c << ' ';
    }
    out << '\n';
  }
  out.close();
  return 0;
}



// -----------------------------------------------------------------------------
// Color function definitions
// -----------------------------------------------------------------------------
Color::Color(dbl r, dbl g, dbl b)
  : red(r), green(g), blue(b) {}

Color& Color::operator*=(const Color& rhs)
{
  red *= rhs.red;
  green *= rhs.green;
  blue *= rhs.blue;
  return *this;
}

Color& Color::operator+=(const Color& rhs)
{
  red += rhs.red;
  green += rhs.green;
  blue += rhs.blue;
  return *this;
}

Color& Color::operator*=(dbl s)
{
  red *= s;
  green *= s;
  blue *= s;
  return *this;
}

bool Color::operator==(const Color& rhs) const
{
  dbl epsilon = std::numeric_limits<dbl>::epsilon();
  return std::abs(red - rhs.red) <= epsilon
    && std::abs(green - rhs.green) <= epsilon
    && std::abs(blue - rhs.blue) <= epsilon;
}

bool Color::operator!=(const Color& rhs) const
{
  return !(*this == rhs);
}

Color operator*(Color lhs, const Color& rhs)
{
  return lhs *= rhs;
}

Color operator+(Color lhs, const Color& rhs)
{
  return lhs += rhs;
}

Color operator*(Color lhs, dbl s)
{
  lhs *= s;
  return lhs;
}

std::ostream& operator<<(std::ostream& out, const Color& c)
{
  i32 r = std::clamp(c.red, 0.0, 1.0) * 255.0;
  i32 g = std::clamp(c.green, 0.0, 1.0) * 255.0;
  i32 b = std::clamp(c.blue, 0.0, 1.0) * 255.0;
  return out << r << ' ' << g << ' ' << b;
}

// -----------------------------------------------------------------------------
// Vec function definitions
// -----------------------------------------------------------------------------
dbl dot(const Vec& u, const Vec& v);

Vec::Vec(dbl x, dbl y, dbl z)
  : x(x), y(y), z(z) {}

Vec Vec::operator-() const
{
  return *this * -1;
}

Vec& Vec::operator*=(dbl s)
{
  x *= s;
  y *= s;
  z *= s;
  return *this;
}

Vec& Vec::operator+=(const Vec& rhs)
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

Vec& Vec::operator-=(const Vec& rhs)
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}

Vec operator+(Vec lhs, const Vec& rhs)
{
  lhs += rhs;
  return lhs;
}

Vec operator-(Vec lhs, const Vec& rhs)
{
  lhs -= rhs;
  return lhs;
}

Vec operator*(Vec v, dbl s)
{
  v *= s;
  return v;
}

bool Vec::operator==(const Vec& rhs) const
{
  dbl epsilon = std::numeric_limits<dbl>::epsilon();
  return std::abs(x - rhs.x) <= epsilon 
      && std::abs(y - rhs.y) <= epsilon 
      && std::abs(z - rhs.z) <= epsilon;
}

bool Vec::operator!=(const Vec& rhs) const
{
  return !(*this == rhs);
}

dbl Vec::magnitude2() const
{
  return dot(*this, *this);
}

dbl Vec::magnitude() const
{
  return sqrt(magnitude2());
}

Vec& Vec::normalize()
{
  dbl l = magnitude();
  x /= l;
  y /= l;
  z /= l;
  return *this;
}

Vec cross(const Vec& u, const Vec& v)
{
  return Vec{u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x};
}

dbl dot(const Vec& u, const Vec& v)
{
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

std::ostream& operator<<(std::ostream& out, const Vec& v)
{
  return out << v.x << ' ' << v.y << ' ' << v.z;
}

// -----------------------------------------------------------------------------
// Ray function definitions
// -----------------------------------------------------------------------------
Ray::Ray(Vec base, Vec dir)
  : base(std::move(base)), dir(std::move(dir))
{ }

// -----------------------------------------------------------------------------
// Material function defintions
// -----------------------------------------------------------------------------
Material::Material(Color c)
  : k_ambient(1), k_diffuse(0), k_specular(0), phong_exp(0),
  ambient(c), diffuse(), specular()
{}

Material::Material(Color c, dbl k_a, dbl k_d)
  : k_ambient(k_a), k_diffuse(k_d), k_specular(0), phong_exp(0),
  ambient(c), diffuse(c), specular()
{}

Material::Material(Color c, dbl k_a, dbl k_d, dbl k_s, dbl phong)
  : k_ambient(k_a), k_diffuse(k_d), k_specular(k_s), phong_exp(phong),
  ambient(c), diffuse(c), specular(Color{1, 1, 1})
{}

// -----------------------------------------------------------------------------
// Sphere function definitions
// -----------------------------------------------------------------------------
Sphere::Sphere(Vec center, dbl radius)
  : center(std::move(center)), radius(std::move(radius)) 
{}

bool Sphere::hits(const Ray& ray, dbl mn, dbl* mx, Vec* pt) const
{
  Vec bc = ray.base - center;
  dbl b = dot(bc, ray.dir);
  dbl c = dot(bc, bc) - radius * radius;
  if (b > 0 && c > 0)
    return false;
  dbl disc = b * b - c;
  if (disc < 0)
    return false;
  dbl t = -b - sqrt(disc);
  if (t <= mn || t >= *mx)
    return false;
  *mx = t;
  *pt = ray.base + ray.dir * t;
  return true;
}

Vec Sphere::normal(const Vec& pt) const
{
  return (pt - center).normalize();
}

// -----------------------------------------------------------------------------
// Object function definitions
// -----------------------------------------------------------------------------
Object::Object(Sphere s, Material m)
  : material(std::move(m)), d_tag(SPHERE), d_shape(std::move(s))
{}

bool Object::hits(const Ray& ray, dbl mn, dbl* mx, HitRecord* hr) const
{
  switch (d_tag) {
    case SPHERE:
      if (d_shape.sphere.hits(ray, mn, mx, &hr->pt)) {
        hr->ray = ray;
        hr->object = this;
        return true;
      }
  }
  return false;
}

Vec Object::normal(const Vec& pt) const
{
  switch (d_tag) {
    case SPHERE:
      return d_shape.sphere.normal(pt);
  }
  assert(false);
  return Vec{};
}


// -----------------------------------------------------------------------------
// Eye function definitions
// -----------------------------------------------------------------------------
Eye::Eye(Vec pos, const Vec& look_at, const Vec& up)
  : pos(std::move(pos))
{
    w = (pos - look_at).normalize(),
    u = cross(up, w).normalize();
    v = cross(w, u);
}

// -----------------------------------------------------------------------------
// Scene function defintions
// -----------------------------------------------------------------------------
Color Scene::color(const HitRecord& hr) const
{
  Color c(hr.object->material.ambient * hr.object->material.k_ambient);
  Vec normal = hr.object->normal(hr.pt);
  for (const auto& light : lights) {
    dbl l_dot_n = dot(light.pos, normal);
    if (l_dot_n >= 0)
      c += hr.object->material.diffuse * hr.object->material.k_diffuse * l_dot_n;
    Vec l = (light.pos - hr.pt).normalize();
    Vec h = (-hr.ray.dir + l).normalize();
    dbl h_dot_n = std::max(0.0, dot(h, normal));
    c += hr.object->material.specular * hr.object->material.k_specular
      * std::pow(h_dot_n, hr.object->material.phong_exp);
    c *= light.color;
  }
  return c;
}

template <typename F>
Color Scene::trace(const Ray& ray, F begin, F end) const
{
  bool hit = false;
  HitRecord hr;
  dbl mn = EPSILON;
  dbl mx = std::numeric_limits<dbl>::max();
  for (; begin != end; ++begin)
    hit = begin->hits(ray, mn, &mx, &hr) || hit;
  if (!hit) return Color();
  return color(hr);
}

Ray Scene::get_ray(i32 i, i32 j) const
{
  const dbl s_u = (i + 0.5) / img_wd * frame.width - frame.width / 2;
  const dbl s_v = (j + 0.5) / img_ht * frame.height - frame.height / 2;
  const dbl s_w = frame.distance;
  const Vec u = eye.u * s_u;
  const Vec v = eye.v * s_v;
  const Vec w = eye.w * s_w;
  const Vec pt = u + v + w;
  Vec dir = pt - eye.pos;
  return Ray{eye.pos, dir.normalize()};
}

// -----------------------------------------------------------------------------
// TESTS
// -----------------------------------------------------------------------------
bool test_sphere_hit();
bool test_get_ray();
bool test_eye();
bool test_trace();
bool test_get_ray_and_trace();


#define TEST_ASSERT(cond) \
  if (!cond) { \
    std::cout << "Expected " << #cond << "\n"; \
    return false; \
  } \


#define TEST_EQ(a, b) \
  if (a != b) { \
    std::cout << #a << ": " << a << " != " << #b ": " << b << "\n"; \
    return false; \
  } \


#define TEST_DOUBLE_EQ(a, b) \
  if (a - b >= std::numeric_limits<dbl>::epsilon()) { \
    std::cout << #a << ": " << a << " != " << #b ": " << b << "\n"; \
    return false; \
  } \


#define RUN_TEST(test) \
  std::cout << "Running " << #test << std::endl; \
  if (test()) \
    std::cout << ". . . PASS\n"; \
  else \
    std::cout << ". . . FAIL\n"; \


bool test_sphere_hit()
{
  // IF
  Sphere s{Vec{}, 1};
  Ray r{Vec{0, 0, 10}, Vec{0, 0, -1}};
  dbl mx = std::numeric_limits<dbl>::max();
  Vec pt;
  // WHEN
  bool hit = s.hits(r, EPSILON, &mx, &pt);
  // THEN
  TEST_ASSERT(hit);
  TEST_EQ(pt, (Vec{0, 0, 1}));
  return true;
}

bool test_get_ray()
{
  // IF
  i32 wd = 800;
  i32 ht = 600;
  Eye eye{Vec{0, 0, 10}, Vec{}, Vec{0, 1, 0}};
  Frame frame{static_cast<dbl>(wd)/100, static_cast<dbl>(ht)/100};
  Scene scene{eye, frame, wd, ht};
  // WHEN
  Ray ray = scene.get_ray(wd/2, ht/2);
  // THEN
  TEST_EQ(ray.base, eye.pos);
  TEST_EQ(ray.dir, (Vec{0.0005, 0.0005, -1}));
  return true;
}

bool test_eye()
{
  // IF
  Vec pt{0, 0, 10};
  Vec look_at;
  Vec up{0, 1, 0};
  // WHEN
  Eye eye{pt, look_at, up};
  // THEN
  TEST_EQ(eye.pos, pt);
  TEST_EQ(eye.u, (Vec{1, 0, 0}));
  TEST_EQ(eye.v, (Vec{0, 1, 0}));
  TEST_EQ(eye.w, (Vec{0, 0, 1}));
  return true;
}

bool test_object_hits()
{
  // IF
  Color green{0, 1, 0};
  Object obj{Sphere{Vec{}, 1}, Material{green}};
  Ray ray{Vec{0, 0, 10}, {0, 0, -1}};
  dbl mx = std::numeric_limits<dbl>::max();
  HitRecord hr;
  // WHEN
  bool hit = obj.hits(ray, EPSILON, &mx, &hr);
  // THEN
  TEST_ASSERT(hit);
  return true;
}

bool test_trace()
{
  // IF
  Color green{0, 1, 0};
  const std::array<Object, 1> objs = {{
    Object{Sphere{Vec{}, 1}, Material{green}}
  }};
  Ray ray{Vec{0, 0, 10}, {0, 0, -1}};
  Scene scene;
  // WHEN
  Color c = scene.trace(ray, objs.cbegin(), objs.cend());
  // THEN
  TEST_EQ(c, green);
  return true;
}

bool test_get_ray_and_trace()
{
  // IF
  i32 wd = 800, ht = 600;
  Color green{0, 1, 0};
  const std::array<Object, 1> objs = {{
    Object{Sphere{Vec{}, 1}, Material{green}}
  }};
  Eye eye{Vec{0, 0, 10}, Vec{}, Vec{0, 1, 0}};
  Frame frame{8, 6, 6};
  Scene scene{eye, frame, wd, ht};
  // WHEN
  Ray ray = scene.get_ray(wd/2, ht/2);
  Color c = scene.trace(ray, objs.cbegin(), objs.cend());
  // THEN
  TEST_EQ(ray.base, eye.pos);
  TEST_EQ(ray.dir, (Vec{0, 0, -1}));
  TEST_EQ(c, green);
  return true;
}

int run_all_tests()
{
  RUN_TEST(test_sphere_hit);
  RUN_TEST(test_eye);
  RUN_TEST(test_get_ray);
  RUN_TEST(test_trace);
  RUN_TEST(test_object_hits);
  RUN_TEST(test_get_ray_and_trace);
  return 0;
}
