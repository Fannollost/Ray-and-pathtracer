#include "precomp.h"
/*
class material {
public:
	virtual bool scatter(
		Ray& ray, float3& att, Ray& scattered, float3 normal
	) {};
};

class diffuse : public material {
public:
	diffuse(float3 a) : albedo(a) {}

	virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 normal) const {
		float3 dir = normal + 
		UnitVector();
		if (isZero(dir)) dir = normal;
		scattered = Ray(ray.IntersectionPoint(), dir, ray.color);
		att = albedo;
		return true;
	}

public:
	float3 albedo;

};	 */
