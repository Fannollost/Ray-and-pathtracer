#pragma once
// -----------------------------------------------------------
// scene.h
// Simple test scene for ray tracing experiments. Goals:
// - Super-fast scene intersection
// - Easy interface: scene.FindNearest / IsOccluded
// - With normals and albedo: GetNormal / GetAlbedo
// - Area light source (animated), for light transport
// - Primitives can be hit from inside - for dielectrics
// - Can be extended with other primitives and/or a BVH
// - Optionally animated - for temporal experiments
// - Not everything is axis aligned - for cache experiments
// - Can be evaluated at arbitrary time - for motion blur
// - Has some high-frequency details - for filtering
// Some speed tricks that severely affect maintainability
// are enclosed in #ifdef SPEEDTRIX / #endif. Mind these
// if you plan to alter the scene in any way.
// -----------------------------------------------------------

#define SPEEDTRIX

#define PLANE_X(o,i) {if((t=-(ray.O.x+o)*ray.rD.x)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Y(o,i) {if((t=-(ray.O.y+o)*ray.rD.y)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Z(o,i) {if((t=-(ray.O.z+o)*ray.rD.z)<ray.t)ray.t=t,ray.objIdx=i;}

namespace Tmpl8 {
	class material;
	class diffuse;
	class metal;
__declspec(align(64)) class Ray
{
public:
	Ray() = default;
	Ray( float3 origin, float3 direction, float3 color ,float distance = 1e34f)
	{
		O = origin, D = direction, t = distance;
		// calculate reciprocal ray direction for triangles and AABBs
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
	#ifdef SPEEDTRIX
		d0 = d1 = d2 = 0;
	#endif
	}
	float3 IntersectionPoint() { return O + t * D; }
	material* GetMaterial() { return m; }
	void SetInside(float3 normal) {
		inside = dot(normalize(D), normal) < 0;
		hitNormal = inside ? normal : -normal;
	}
	// ray data
#ifndef SPEEDTRIX
	float3 O, D, rD;
#else
	union { struct { float3 O; float d0; }; __m128 O4; };
	union { struct { float3 D; float d1; }; __m128 D4; };
	union { struct { float3 rD; float d2; }; __m128 rD4; };
#endif
	float t = 1e34f;
	int objIdx = -1;
	bool inside = false; // true when in medium
	float3 color = 0;
	float3 hitNormal;
	material* m;
};


// -----------------------------------------------------------
// Triangle Primitive
// 
// 
// -----------------------------------------------------------
class Triangle {
public:
	Triangle() = default;
	Triangle(int idx, float3 ver0, float3 ver1, float3 ver2, material* m) : objIdx(idx), v0(ver0), v1(ver1), v2(ver2), mat(m) {
		e1 = v1 - v0;
		e2 = v2 - v0;
		N = cross(e1, e2);
	}
	void Intersect(Ray& ray, float t_min) const {		 //scratchapixel implementation
		float area = length(N);
		float NdotRayDir = dot(N, ray.D);
		if (fabs(NdotRayDir) < t_min) return;
		float d = -dot(N, v0);
		float t = -(dot(N, ray.O) + d) / NdotRayDir;
		if (t < 0) return;
		float3 p = ray.O + t * ray.D;
		float3 c;

		float3 vp0 = p - v0;
		c = cross(e1, vp0);
		if (dot(N, c) < 0) return;
		float3 vp1 = p - v1;
		float3 e3 = v2 - v1;
		c = cross(e3, vp1);
		if (dot(N, c) < 0) return;
		float3 e4 = v0 - v2;
		float3 vp2 = p - v2;
		c = cross(e4, vp2);
		if (dot(N, c) < 0) return;

		ray.t = t, ray.objIdx = objIdx, ray.m = mat, ray.color = col;
			//ray.SetInside(GetNormal(ray.IntersectionPoint()));

		
	}
	float3 GetNormal(const float3 I) const { return N; }
	float3 v0, v1, v2, e1, e2, N;
	int objIdx = -1;
	float3 col;
	material* mat;
};

class Light {
public: 
	Light() = default;
	Light(int idx, float3 p, float str, float3 c, float3 n) : objIdx(idx), pos(p), strength(str), col(c), normal(n) {}
	float3 GetNormal() { return normal; }
	virtual float3 GetLightPosition() { return pos; }
	float3 GetLightColor() { return col; }
	virtual float3 GetLightIntensityAt(float3 p, float3 n, const material& m) { return 1; }
	float3 pos;
	float3 col;
	float strength;
	int objIdx;
	float3 normal;
};

class AreaLight : public Light {
public:
	AreaLight() = default;
	AreaLight(int idx, float3 p, float str, float3 c, float r, float3 n, int s) : Light(idx, p, str, c, n) {
		radius = r;
		samples = s;
	}
	float3 GetLightIntensityAt(float3 p, float3 n, const material& m) override {
	/*	float3 lightRayDirection = GetLightPosition() - p;*/
		float dis = abs(length(pos - p));

		//if (m.type == "diff") {
		//	Ray r = Ray(p, normalize(lightRayDirection), m.col, dis);
		//}				   
		float relStr =  1 / (dis) * strength;
		float totCol = 0; 
		float3 dir = pos - p;
		float str = dot(n, normalize(dir));
		if (str < 0.0f) str = 0.0f;

		for (int i = 0; i < samples; i++) {
		 //lerp between values then divide by samples  
		}
		return relStr * str * GetLightColor();
	}
	float3 GetLightPosition() override {
		float newRad = radius * sqrt(RandomFloat());
		float theta = RandomFloat() * 2 * PI;
		//return float3(pos.x + newRad * cos(theta), pos.y + newRad * sin(theta), pos.z);
	   	return pos;
	}
	int samples;
	float radius;
};

class DirectionalLight : public Light {
public:
	DirectionalLight() = default;
	DirectionalLight(int idx, float3 p, float str, float3 c, float3 n, int s, float bs) : Light(idx, p, str, c, n) {
		samples = s;
		beamSize = bs;
	}
	float3 GetLightPosition() override {
		//float3 r = dot(pos), -normal) / length(normal);
		//pos - r;
	}
	float3 GetLightIntensityAt(float3 p, float3 n, const material& m) override {
		float3 dir = p - pos;
		float dis = length(dir);
		float str = abs(dot(n, dir) / (length(n) * length(dir)));
		return (1 / dis * str * strength) * GetLightColor();

		/*float3 str = dot(normal, n);
		float distDiff = length(str - p);
		if(distDiff > 0)
		return 0.0f;*/
		//float3 diff = p - pos;
		//float angle = asin(cross(normal, diff) / (length(normal) * length(diff))) * 180.0f / PI;
		//float dis = abs(length(pos - p));
		//float relStr = 1 / (dis)*strength;
		//return relStr * sqrt(fmaxf(normal / 2.0f) - angle) / (normal / 2), 0)));
	}

	int samples;
	float beamSize;
};
// -----------------------------------------------------------
// Sphere primitive
// Basic sphere, with explicit support for rays that start
// inside it. Good candidate for a dielectric material.
// -----------------------------------------------------------
class Sphere
{
public:
	Sphere() = default;
	Sphere( int idx, float3 p, float r, material* m) : 
		pos( p ), r2( r* r ), invr( 1 / r ), objIdx( idx ),  mat(m) {}
	void Intersect( Ray& ray , float t_min) const
	{
		float3 oc = ray.O - this->pos;
		float b = dot( oc, ray.D );
		float c = dot( oc, oc ) - this->r2;
		float t, d = b * b - c;
		if (d <= 0) return;
		d = sqrtf( d ), t = -b - d;
		if (t < ray.t && t > t_min)
		{
			ray.t = t, ray.objIdx = objIdx, ray.m = mat;
			ray.SetInside(GetNormal(ray.IntersectionPoint()));
			return;
		}
		t = d - b;
		if (t < ray.t && t > t_min)
		{
			ray.t = t, ray.objIdx = objIdx, ray.m = mat;
			ray.SetInside(GetNormal(ray.IntersectionPoint()));
			return;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		return (I - this->pos) * invr;
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 0.93f );
	}
	float3 pos = 0;
	float r2 = 0, invr = 0;
	int objIdx = -1;
	material* mat;
};

// -----------------------------------------------------------
// Plane primitive
// Basic infinite plane, defined by a normal and a distance
// from the origin (in the direction of the normal).
// -----------------------------------------------------------
class Plane
{
public:
	Plane() = default;
	Plane( int idx, float3 normal, float dist, material* m) : N( normal ), d( dist ), objIdx( idx ), mat(m) {}
	void Intersect( Ray& ray, float t_min) const
	{
		float t = -(dot( ray.O, this->N ) + this->d) / (dot( ray.D, this->N ));
		if (t < ray.t && t > t_min) ray.t = t, ray.objIdx = objIdx, ray.m = mat,
			ray.SetInside(N);
	}
	float3 GetNormal( const float3 I ) const
	{
		return N;
	}
	float3 GetAlbedo( const float3 I ) const
	{
		if (N.y == 1)
		{
			// floor albedo: checkerboard
			int ix = (int)(I.x * 2 + 96.01f);
			int iz = (int)(I.z * 2 + 96.01f);
			// add deliberate aliasing to two tile
			if (ix == 98 && iz == 98) ix = (int)(I.x * 32.01f), iz = (int)(I.z * 32.01f);
			if (ix == 94 && iz == 98) ix = (int)(I.x * 64.01f), iz = (int)(I.z * 64.01f);
			return float3( ((ix + iz) & 1) ? 1 : 0.3f );
		}
		else if (N.z == -1)
		{
			// back wall: logo
			static Surface logo( "assets/logo.png" );
			int ix = (int)((I.x + 4) * (128.0f / 8));
			int iy = (int)((2 - I.y) * (64.0f / 3));
			uint p = logo.pixels[(ix & 127) + (iy & 63) * 128];
			uint3 i3( (p >> 16) & 255, (p >> 8) & 255, p & 255 );
			return float3( i3 ) * (1.0f / 255.0f);
		}
		return float3( 0.93f );
	}
	float3 N;
	float d;
	int objIdx = -1;
	material* mat;
};

// -----------------------------------------------------------
// Cube primitive
// Oriented cube. Unsure if this will also work for rays that
// start inside it; maybe not the best candidate for testing
// dielectrics.
// -----------------------------------------------------------
class Cube
{
public:
	Cube() = default;
	Cube( int idx, float3 pos, float3 size, material* m, mat4 transform = mat4::Identity() )
	{
		objIdx = idx;
		b[0] = pos - 0.5f * size, b[1] = pos + 0.5f * size;
		mat = m;
		M = transform, invM = transform.FastInvertedTransformNoScale();
	}
	void Intersect( Ray& ray, float t_min ) const
	{
		// 'rotate' the cube by transforming the ray into object space
		// using the inverse of the cube transform.
		float3 O = TransformPosition( ray.O, invM );
		float3 D = TransformVector( ray.D, invM );
		float rDx = 1 / D.x, rDy = 1 / D.y, rDz = 1 / D.z;
		int signx = D.x < 0, signy = D.y < 0, signz = D.z < 0;
		float tmin = (b[signx].x - O.x) * rDx;
		float tmax = (b[1 - signx].x - O.x) * rDx;
		float tymin = (b[signy].y - O.y) * rDy;
		float tymax = (b[1 - signy].y - O.y) * rDy;
		if (tmin > tymax || tymin > tmax) return;
		tmin = max( tmin, tymin ), tmax = min( tmax, tymax );
		float tzmin = (b[signz].z - O.z) * rDz;
		float tzmax = (b[1 - signz].z - O.z) * rDz;
		if (tmin > tzmax || tzmin > tmax) return;
		tmin = max( tmin, tzmin ), tmax = min( tmax, tzmax );
		if (tmin > t_min)
		{
			if (tmin < ray.t) ray.t = tmin, ray.objIdx = objIdx,  ray.m = mat,
				ray.SetInside(GetNormal(ray.IntersectionPoint()));
		}
		else if (tmax > t_min)
		{
			if (tmax < ray.t) ray.t = tmax, ray.objIdx = objIdx,  ray.m = mat,
				ray.SetInside(GetNormal(ray.IntersectionPoint()));
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		// transform intersection point to object space
		float3 objI = TransformPosition( I, invM );
		// determine normal in object space
		float3 N = float3( -1, 0, 0 );
		float d0 = fabs( objI.x - b[0].x ), d1 = fabs( objI.x - b[1].x );
		float d2 = fabs( objI.y - b[0].y ), d3 = fabs( objI.y - b[1].y );
		float d4 = fabs( objI.z - b[0].z ), d5 = fabs( objI.z - b[1].z );
		float minDist = d0;
		if (d1 < minDist) minDist = d1, N.x = 1;
		if (d2 < minDist) minDist = d2, N = float3( 0, -1, 0 );
		if (d3 < minDist) minDist = d3, N = float3( 0, 1, 0 );
		if (d4 < minDist) minDist = d4, N = float3( 0, 0, -1 );
		if (d5 < minDist) minDist = d5, N = float3( 0, 0, 1 );
		// return normal in world space
		return TransformVector( N, M );
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 1, 1, 1 );
	}
	float3 b[2];
	mat4 M, invM;
	int objIdx = -1; 
	material* mat;
};

// -----------------------------------------------------------
// Quad primitive
// Oriented quad, intended to be used as a light source.
// -----------------------------------------------------------
class Quad
{
public:
	Quad() = default;
	Quad(int idx, float s, material* m, mat4 transform = mat4::Identity())
	{
		objIdx = idx;
		size = s * 0.5f;
		mat = m;
		T = transform, invT = transform.FastInvertedTransformNoScale();
	}
	void Intersect( Ray& ray, float t_min ) const
	{
		const float3 O = TransformPosition( ray.O, invT );
		const float3 D = TransformVector( ray.D, invT );
		const float t = O.y / -D.y;
		if (t < ray.t && t > t_min)
		{
			float3 I = O + t * D;
			if (I.x > -size && I.x < size && I.z > -size && I.z < size)
				ray.t = t, ray.objIdx = objIdx, ray.m = mat, 
			ray.SetInside(GetNormal(ray.IntersectionPoint()));
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		// TransformVector( float3( 0, -1, 0 ), T ) 
		return float3( -T.cell[1], -T.cell[5], -T.cell[9] );
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 10 );
	}
	float size;
	mat4 T, invT;
	int objIdx = -1;
	float3 col; 
	material* mat;
};

class material {
public:
	material(float3 c, float spec) : col(c), specularity(spec) {}
	virtual bool scatter(
		Ray& ray, float3& att, Ray& scattered, float3 normal, float3& energy
	) {
		return false;
	}
	float3 col;
	string type;
	float specularity;
};

class diffuse : public material {
public:
	diffuse(float3 a, float3 c, float spec) : albedo(a), material(c, spec) { type = "diff"; }

	virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 normal, float3& energy) override {
		
		float3 dir = ray.IntersectionPoint() + normal; //+ RandomInHemisphere(normal);
		if (isZero(dir)) dir = normal;
		scattered = Ray(ray.IntersectionPoint(), normalize(dir - ray.IntersectionPoint()), ray.color);
		att = albedo;
		float3 retention = float3(1) - albedo;
		float3 newEnergy(energy - retention);
		energy = newEnergy.x > 0 ? newEnergy : 0;
		return true;
	}

public:
	float3 albedo;
};

class metal : public material {
public:
	metal(float3 a, float f, float3 c, float spec) : albedo(a), fuzzy(f < 1 ? f : 1), material(c, spec) { type = "metal"; }
	virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 normal, float3& energy) override {
		float3 dir = reflect(ray.D, normal); //add fuzzy
		scattered = Ray(ray.IntersectionPoint(), dir, ray.color);
		att = albedo;
		energy = energy;
		return dot(scattered.D, normal) > 0;
	}

public:
	float3 albedo;
	float fuzzy;
};

class glass : public material {
public:
	glass(float refIndex, float3 c, float spec, float3 a) : ir(refIndex), absorption(a), material(c, spec) { type = "glass"; }
	 virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 normal, float3& energy) override {
		 att = float3(1.0f);
		 float refrRatio = ray.inside ? (1.0 / ir) : ir;
		 float3 uDir = UnitVector(ray.D);

		 float ctheta = fmin(dot(-uDir, ray.hitNormal), 1.0);			//maybe use normal
		 float stheta = sqrt(1.0 - ctheta * ctheta);
		 float3 refDir;
		 if ((refrRatio * stheta) > 1.0) {
			 refDir = reflect(uDir, ray.hitNormal);
		 }
		 else {
			 refDir = refractRay(uDir, ray.hitNormal, refrRatio);
		 };

		 if (ray.inside)
		 {
			 energy.x *= exp(-absorption.x * ray.t);
			 energy.y *= exp(-absorption.y * ray.t);
			 energy.z *= exp(-absorption.z * ray.t);
		 }
		 else {
			 energy = energy;
		 }
		 scattered = Ray(ray.IntersectionPoint(),refDir, ray.color);
		 return true;
	 }
	 float3 refractRay(float3 oRayDir, float3 normal, float refRatio) {
		 float theta = fmin(dot(-oRayDir, normal), 1.0);
		 float3 perpendicular = refRatio * (oRayDir + theta * normal);
		 float3	parallel = -sqrt(fabs(1.0 - pow(length(perpendicular), 2))) * normal;
		 return perpendicular + parallel;
	 }
	 float ir;

	 float3 fresnel(float3 I, float3 normal, float ior, float &kr) {
		 float cosi = clamp(-1.0f, 1.0f, dot(I, normal));
		 float etai = 1, etat = ior;
		 if (cosi > 0) { std::swap(etai, etat); }
		 float sint = etai / etat * sqrtf(std::max(0.0f, 1 - cosi * cosi));
		 if (sint >= 1) {
			 kr = 1;
		 }
		 else {
			 float cost = sqrtf(std::max(0.0f, 1 - sint * sint));
			 cosi = fabsf(cosi);
			 float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));	    
			 float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
			 kr = (Rs * Rs + Rp * Rp) / 2;
		 }

	 }//scratchapixel

	 float3 absorption;
};

// -----------------------------------------------------------
// Scene class
// We intersect this. The query is internally forwarded to the
// list of primitives, so that the nearest hit can be returned.
// For this hit (distance, obj id), we can query the normal and
// albedo.
// -----------------------------------------------------------
class Scene
{
public:
	Scene()
	{
		float3 white = float3(1.0, 1.0, 1.0);
		float3 red = float3(1.0, 0, 0);
		float3 blue = float3(0, 1.0, 0);
		float3 green = float3(0, 0, 1.0);
		// we store all primitives in one continuous buffer
		quad = Quad(0, 1, new diffuse(float3(0.8f), white, 0));									// 0: light source
		light[0] = new AreaLight(11, float3(0.1f, 1, 0), 1.0f,  white, 0.1f, float3(0, -1, 0), 4);			//DIT FF CHECKEN!
		light[1] = new AreaLight(12, float3(0.1f, -1, 0), 1.0f,  white, 0.1f, float3(0, -1, 0), 4);			//DIT FF CHECKEN!
		sphere = Sphere( 1, float3( 0 ), 0.5f, new glass(1.5f, red, 1.0f, 0.0f));				// 1: bouncing ball
		//sphere = Sphere( 1, float3( 0 ), 0.5f, new metal(0.8f, 1.0f, blue, 0.5f));				// 1: bouncing ball
		sphere2 = Sphere( 2, float3( 1,-0.1f,0 ), 0.2f, new metal(1.0f, 1.0f, white, 0.8f));				// 1: bouncing ball
		//sphere3 = Sphere( 3, float3( -0.3f,-0.1f,-0.2f ), 0.2f, white, new glass(0.1f));				// 1: bouncing ball
		//sphere2 = Sphere( 2, float3( 0, 2.5f, -3.07f ), 8, blue, new diffuse(float3(0.2f)));	// 2: rounded corners
		//cube = Cube( 3, float3( 0 ), float3( 1.15f ) , green, new diffuse(float3(0.5f)));		// 3: cube
		plane[0] = Plane( 4, float3( 1, 0, 0 ), 3, new diffuse(0.8f, blue, 0.0f));			// 4: left wall
		plane[1] = Plane( 5, float3( -1, 0, 0 ), 2.99f, new diffuse(0.8f, red, 0.0f));		// 5: right wall
		plane[2] = Plane( 6, float3( 0, 1, 0 ), 1, new diffuse(0.8f, white, 0.0f));			// 6: floor
		plane[3] = Plane( 7, float3( 0, -1, 0 ), 2, new diffuse(0.8f, white, 0.0f));			// 7: ceiling
		plane[4] = Plane( 8, float3( 0, 0, 1 ), 3, new diffuse(0.8f, red, 0.0f));				// 8: front wall
		plane[5] = Plane( 9, float3( 0, 0, -1 ), 3.99f, new diffuse(0.8f, green, 0.0f));		// 9: back wall
		triangle = Triangle(10, float3(0.0f, 0.0f, 1.0f), float3(0.2f, 0, 1.0f), float3(0.1f, 0.2f, 1.0f), new diffuse(0.8f, blue, 0.0f));
		SetTime( 0 );
		// Note: once we have triangle support we should get rid of the class
		// hierarchy: virtuals reduce performance somewhat.
	}
	void SetTime( float t )
	{
		// default time for the scene is simply 0. Updating/ the time per frame 
		// enables animation. Updating it per ray can be used for motion blur.
		animTime = t;
		// light source animation: swing
		mat4 M1base = mat4::Translate( float3( 0, 2.6f, 2 ) );
		mat4 M1 = M1base * mat4::RotateZ( sinf( animTime * 0.6f ) * 0.1f ) * mat4::Translate( float3( 0, -0.9, 0 ) );
		quad.T = M1, quad.invT = M1.FastInvertedTransformNoScale();
		// cube animation: spin
		mat4 M2base = mat4::RotateX( PI / 4 ) * mat4::RotateZ( PI / 4 );
		mat4 M2 = mat4::Translate( float3( 1.4f, 0, 2 ) ) * mat4::RotateY( animTime * 0.5f ) * M2base;
		cube.M = M2, cube.invM = M2.FastInvertedTransformNoScale();
		// sphere animation: bounce
		float tm = 1 - sqrf( fmodf( animTime, 2.0f ) - 1 );
		sphere.pos = float3( -1.4f, -0.5f + tm, 2 );
	}
	float3 GetLightPos() const
	{
		// light point position is the middle of the swinging quad
		float3 corner1 = TransformPosition( float3( -0.5f, 0, -0.5f ), quad.T );
		float3 corner2 = TransformPosition( float3( 0.5f, 0, 0.5f ), quad.T );
		return (corner1 + corner2) * 0.5f - float3( 0, 0.01f, 0 );
	}
	float3 GetLightColor() const
	{
		return quad.col;
	}
	void FindNearest( Ray& ray, float t_min) const
	{
		// room walls - ugly shortcut for more speed
		float t;
		for (int i = 0; i < size(plane); i++) plane[i].Intersect(ray, t_min);
		//if (ray.D.x < 0) PLANE_X( 3, 4 ) else PLANE_X( -2.99f, 5 );
		//if (ray.D.y < 0) PLANE_Y( 1, 6 ) else PLANE_Y( -2, 7 );
		//if (ray.D.z < 0) PLANE_Z( 3, 8 ) else PLANE_Z( -3.99f, 9 );
		quad.Intersect( ray , t_min);
		sphere.Intersect( ray, t_min );
		sphere2.Intersect( ray, t_min );
		triangle.Intersect(ray, t_min);
		//cube.Intersect( ray, t_min);
		//triangle.Intersect(ray, t_min);
		//triangle2.Intersect(ray, t_min);
	}
	bool IsOccluded( Ray& ray, float t_min) const
	{
		float rayLength = ray.t;
		// skip planes: it is not possible for the walls to occlude anything
		quad.Intersect( ray, t_min );
		sphere.Intersect( ray, t_min);
		sphere2.Intersect( ray, t_min);
		triangle.Intersect(ray, t_min);

		//cube.Intersect( ray , t_min);
		//triangle.Intersect(ray, t_min);
		//triangle2.Intersect(ray, t_min);

		return ray.t < rayLength;
		// technically this is wasteful: 
		// - we potentially search beyond rayLength
		// - we store objIdx and t when we just need a yes/no
		// - we don't 'early out' after the first occlusion
	}
	float3 GetNormal( int objIdx, float3 I, float3 wo ) const
	{
		// we get the normal after finding the nearest intersection:
		// this way we prevent calculating it multiple times.
		if (objIdx == -1) return float3( 0 ); // or perhaps we should just crash
		float3 N;
		if (objIdx == 0) N = quad.GetNormal(I);
		else if (objIdx == 1) N = sphere.GetNormal(I);
		else if (objIdx == 2) N = sphere2.GetNormal(I);
		//else if (objIdx == 3) N = cube.GetNormal(I);
		else if (objIdx == 10) N = triangle.GetNormal(I);
		else 
		{
			// faster to handle the 6 planes without a call to GetNormal
			N = float3( 0 );
			N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
		}
		if (dot( N, wo ) > 0) N = -N; // hit backside / inside
		return N;
	}
	float3 GetAlbedo( int objIdx, float3 I ) const
	{
		if (objIdx == -1) return float3( 0 ); // or perhaps we should just crash
		if (objIdx == 0) return quad.GetAlbedo( I );
		if (objIdx == 1) return sphere.GetAlbedo( I );
		if (objIdx == 2) return sphere2.GetAlbedo( I );
		//if (objIdx == 3) return cube.GetAlbedo( I );
		return plane[objIdx - 4].GetAlbedo( I );
		// once we have triangle support, we should pass objIdx and the bary-
		// centric coordinates of the hit, instead of the intersection location.
	}
	float GetReflectivity( int objIdx, float3 I ) const
	{
		if (objIdx == 1 /* ball */) return 1;
		if (objIdx == 6 /* floor */) return 0.3f;
		return 0;
	}
	float GetRefractivity( int objIdx, float3 I ) const
	{
		return objIdx == 3 ? 1.0f : 0.0f;
	}
	__declspec(align(64)) // start a new cacheline here
	float animTime = 0;

	Quad quad;
	Light* light[2];
	Sphere sphere;
	Sphere sphere2;
	Sphere sphere3;
	Cube cube;
	Plane plane[6];
	Triangle triangle;
	Triangle triangle2;
	int aaSamples = 1;
	float mediumIr = 1.0f;
};
}