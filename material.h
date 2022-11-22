/*#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H

class material
{
public: 
	virtual bool scatter(
		Ray& ray, float3& att, Ray& scattered, float3 normal
	);
};

class diffuse : public material {
public:
	diffuse(float3 a) : albedo(a) {}
	virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 normal) const;
	float3 albedo;
};						 
#endif // !MATERIAL_H 	*/

