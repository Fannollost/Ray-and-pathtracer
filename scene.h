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

namespace Tmpl8 {

__declspec(align(64)) class Ray
{
public:
	Ray() = default;
	Ray( float3 origin, float3 direction, float distance = 1e34f )
	{
		O = origin, D = direction, t = distance;
		// calculate reciprocal ray direction for triangles and AABBs
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
	#ifdef SPEEDTRIX
		d0 = d1 = d2 = 0;
	#endif
	}
	float3 IntersectionPoint() { return O + t * D; }
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
};

class Primitive
{
public:
	virtual void Intersect( const int idx, Ray& ray ) const = 0;
	virtual float3 GetNormal( const float3 I ) const = 0;
	virtual float3 GetAlbedo( const float3 I ) const = 0;
};

// -----------------------------------------------------------
// Sphere primitive
// Basic sphere, with explicit support for rays that start
// inside it. Good candidate for a dielectric material.
// -----------------------------------------------------------
class Sphere : public Primitive
{
public:
	Sphere() = default;
	Sphere( float3 p, float r ) : pos( p ), r( r ), invr( 1 / r ) {}
	void Intersect( const int idx, Ray& ray ) const
	{
		float3 oc = ray.O - this->pos;
		float b = dot( oc, ray.D );
		float c = dot( oc, oc ) - this->r * this->r;
		float t, d = b * b - c;
		if (d <= 0) return;
		d = sqrtf( d ), t = -b - d;
		if (t < ray.t && t > 0)
		{
			ray.t = t, ray.objIdx = idx;
			return;
		}
		t = d - b;
		if (t < ray.t && t > 0)
		{
			ray.t = t, ray.objIdx = idx;
			return;
		}
	}
	float3 GetNormal( const float3 I ) const
	{
		return (I - this->pos) * invr;
	}
	float3 GetAlbedo( const float3 I ) const
	{
		return float3( 1, 0.2f, 0.2f );
	}
	float3 pos;
	float r, invr;
};

// -----------------------------------------------------------
// Plane primitive
// Basic infinite plane, defined by a normal and a distance
// from the origin (in the direction of the normal).
// -----------------------------------------------------------
class Plane : public Primitive
{
public:
	Plane() = default;
	Plane( float3 normal, float dist ) : N( normal ), d( dist ) {}
	void Intersect( const int idx, Ray& ray ) const
	{
		float t = -(dot( ray.O, this->N ) + this->d) / (dot( ray.D, this->N ));
		if (t < ray.t && t > 0) ray.t = t, ray.objIdx = idx;
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
};

// -----------------------------------------------------------
// Cube primitive
// Oriented cube. Unsure if this will also work for rays that
// start inside it; maybe not the best candidate for testing
// dielectrics.
// -----------------------------------------------------------
class Cube : public Primitive
{
public:
	Cube() = default;
	Cube( float3 pos, float3 size, mat4 transform = mat4::Identity() )
	{
		b[0] = pos - 0.5f * size, b[1] = pos + 0.5f * size;
		M = transform, invM = transform.FastInvertedTransformNoScale();
	}
	void Intersect( const int idx, Ray& ray ) const
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
		if (tmin > 0)
		{
			if (tmin < ray.t) ray.t = tmin, ray.objIdx = idx;
		}
		else if (tmax > 0)
		{
			if (tmax < ray.t) ray.t = tmax, ray.objIdx = idx;
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
		return float3( 0.2f, 1, 0.2f );
	}
	float3 b[2];
	mat4 M, invM;
};

// -----------------------------------------------------------
// Quad primitive
// Oriented quad, intended to be used as a light source.
// -----------------------------------------------------------
class Quad : public Primitive
{
public:
	Quad() = default;
	Quad( float s, mat4 transform = mat4::Identity() )
	{
		size = s * 0.5f;
		T = transform, invT = transform.FastInvertedTransformNoScale();
	}
	void Intersect( const int idx, Ray& ray ) const
	{
		const float3 O = TransformPosition( ray.O, invT );
		const float3 D = TransformVector( ray.D, invT );
		const float t = O.y / -D.y;
		if (t < ray.t && t > 0)
		{
			float3 I = O + t * D;
			if (I.x > -size && I.x < size && I.z > -size && I.z < size)
				ray.t = t, ray.objIdx = idx;
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
		// we store all primitives in one continuous buffer
		quad = Quad( 1 );
		sphere = Sphere( float3( 0 ), 0.5f );
		sphere2 = Sphere( float3( 0, 2.5f, -3.07f ), 8 );
		cube = Cube( float3( 0 ), float3( 1.15f ) );
		plane[0] = Plane( float3( 1, 0, 0 ), 3 ); // left wall
		plane[1] = Plane( float3( -1, 0, 0 ), 2.99f ); // right wall
		plane[2] = Plane( float3( 0, 1, 0 ), 1 ); // floor
		plane[3] = Plane( float3( 0, -1, 0 ), 2 ); // ceiling
		plane[4] = Plane( float3( 0, 0, 1 ), 3 ); // front wall
		plane[5] = Plane( float3( 0, 0, -1 ), 3.99f ); // back wall
		primitive.push_back( &quad ); // light
		primitive.push_back( &sphere ); // sphere
		primitive.push_back( &sphere2 ); // room corner detail
		primitive.push_back( &cube ); // cube
		for (int i = 0; i < 6; i++) primitive.push_back( &plane[i] );
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
		static mat4 M1base = mat4::Translate( float3( 0, 2.6f, 2 ) );
		mat4 M1 = M1base * mat4::RotateZ( sinf( animTime * 0.6f ) * 0.1f ) * mat4::Translate( float3( 0, -0.9, 0 ) );
		((Quad*)primitive[0])->T = M1;
		((Quad*)primitive[0])->invT = M1.FastInvertedTransformNoScale();
		// cube animation: spin
		mat4 M2base = mat4::RotateX( PI / 4 ) * mat4::RotateZ( PI / 4 );
		mat4 M2 = mat4::Translate( float3( 1.4f, 0, 2 ) ) * mat4::RotateY( animTime * 0.5f ) * M2base;
		((Cube*)primitive[3])->M = M2;
		((Cube*)primitive[3])->invM = M2.FastInvertedTransformNoScale();
		// sphere animation: bounce
		float tm = 1 - sqrf( fmodf( animTime, 2.0f ) - 1 );
		((Sphere*)primitive[1])->pos = float3( -1.4f, -0.5f + tm, 2 );
	}
	void FindNearest( Ray& ray )
	{
	#ifdef SPEEDTRIX
		size_t s = 4;
		// room walls - ugly shortcut for more speed
	#if 0
		static __m128 offs0 = _mm_setr_ps( 3, 1, 3, 0 );
		static __m128 offs1 = _mm_setr_ps( -2.99f, -2, -3.99f, 0 );
		static __m128 far4 = _mm_set1_ps( 1e34f );
		static __m128 zero4 = _mm_setzero_ps();
		__m128 offs = _mm_blendv_ps( offs1, offs0, _mm_cmplt_ps( ray.D4, zero4 ) );
		__m128 t4 = _mm_sub_ps( zero4, _mm_mul_ps( _mm_add_ps( ray.O4, offs ), ray.rD4 ) );
		if (t4.m128_f32[0] < ray.t) ray.t = t4.m128_f32[0], ray.objIdx = ray.D.x < 0 ? 4 : 5;
		if (t4.m128_f32[1] < ray.t) ray.t = t4.m128_f32[1], ray.objIdx = ray.D.y < 0 ? 6 : 7;
		if (t4.m128_f32[2] < ray.t) ray.t = t4.m128_f32[2], ray.objIdx = ray.D.z < 0 ? 8 : 9;
	#else
		// somehow this is faster than the sse code?
		if (ray.D.x < 0)
		{
			float t = -(ray.O.x + 3) * ray.rD.x;
			if (t < ray.t) ray.t = t, ray.objIdx = 4;
		}
		else
		{
			float t = -(ray.O.x - 2.99f) * ray.rD.x;
			if (t < ray.t) ray.t = t, ray.objIdx = 5;
		}
		if (ray.D.y < 0)
		{
			float t = -(ray.O.y + 1) * ray.rD.y;
			if (t < ray.t) ray.t = t, ray.objIdx = 6;
		}
		else
		{
			float t = -(ray.O.y - 2) * ray.rD.y;
			if (t < ray.t) ray.t = t, ray.objIdx = 7;
	}
		if (ray.D.z < 0)
		{
			float t = -(ray.O.z + 3) * ray.rD.z;
			if (t < ray.t) ray.t = t, ray.objIdx = 8;
		}
		else
		{
			float t = -(ray.O.z - 3.99f) * ray.rD.z;
			if (t < ray.t) ray.t = t, ray.objIdx = 9;
		}
	#endif
	#else
		// process all entries in the primitive vector
		size_t s = primitive.size();
	#endif
		for (size_t idx = 0; idx < s; idx++)
		{
			primitive[idx]->Intersect( (int)idx, ray );
		}
}
	bool IsOccluded( Ray& ray )
	{
		float rayLength = ray.t;
	#ifndef SPEEDTRIX
		FindNearest( ray );
	#else
		// it is not possible for the walls to occlude anything
		for (size_t idx = 0; idx < 3; idx++)
			primitive[idx]->Intersect( (int)idx, ray );
	#endif
		return ray.t < rayLength;
		// technically this is wasteful: 
		// - we potentially search beyond rayLength
		// - we store objIdx and t when we just need a yes/no
		// - we don't 'early out' after the first occlusion
	}
	float3 GetNormal( int objIdx, float3 I, float3 wo )
	{
		// we get the normal after finding the nearest intersection:
		// this way we prevent calculating it multiple times.
		if (objIdx == -1) return float3( 0 ); // or perhaps we should just crash
		float3 N;
	#ifndef SPEEDTRIX
		N = primitive[objIdx]->GetNormal( I );
	#else
		if (objIdx < 4) N = primitive[objIdx]->GetNormal( I ); else
		{
			// faster to handle the 6 planes without a call to GetNormal
			N = float3( 0 );
			N[(objIdx - 4) / 2] = 1 - 2 * (float)(objIdx & 1);
		}
	#endif
		if (dot( N, wo ) > 0) N = -N; // hit backside / inside
		return N;
	}
	float3 GetAlbedo( int objIdx, float3 I )
	{
		if (objIdx == -1) return float3( 0 ); // or perhaps we should just crash
		return primitive[objIdx]->GetAlbedo( I );
		// once we have triangle support, we should pass objIdx and the bary-
		// centric coordinates of the hit, instead of the intersection location.
	}
	__declspec(align(64)) // start a new cacheline here
		float animTime = 0;
	vector<Primitive*> primitive;
	Quad quad;
	Sphere sphere;
	Sphere sphere2;
	Cube cube;
	Plane plane[6];
};

}