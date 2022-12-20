#pragma once
#define TRIANGLES 0
#define SPHERES 1
namespace Tmpl8{
	class Scene;
	class Ray;
	class DataCollector;
	class Mesh;
	class Triangle;
struct BVHNode
{
	union
	{
		struct { float3 aabbMin; uint leftFirst; };
		__m128 aabbMin4;
	};
	union
	{
		struct { float3 aabbMax; uint primCount; };
		__m128 aabbMax4;
	};
	bool isLeaf() { return primCount > 0; }
};

struct aabb
{
	float3 bmin = 1e30f, bmax = -1e30f;
	void grow(float3 p) { bmin = fminf(bmin, p); bmax = fmaxf(bmax, p); }
	void grow(aabb& b) { if (b.bmin.x != 1e30f) { grow(b.bmin); grow(b.bmax); } }
	float area()
	{
		float3 e = bmax - bmin; // box extent
		return e.x * e.y + e.y * e.z + e.z * e.x;
	}
};

enum SplitMethod {
	BINNEDSAH = 0,
	MIDDLE = 1,
	LONGESTAXIS = 2,
	SAH = 3
};

class bvh
{
	public:
		bvh(Scene* s);
		bvh(Mesh* m);

		void Build();
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint rootNodeIdx);
		void Intersect(Ray& ray);
		static float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		float IntersectAABB_SSE(const Ray& ray, const __m128 bmin4, const __m128 bmax4);
		float EvaluateSAH(BVHNode &node, int axis, float pos);
		float CalculateNodeCost(BVHNode& node);
		float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
		
		void Split(uint nodeIdx);
		void SubdividePrim(uint rootNodeIdx);
		bool IsOccluded(Ray& ray);
		void Refit();
		Triangle getTriangle(uint idx);
	public:
		uint rootNodeIdx = 0, nodesUsed = 2, NTri = 0, NSph = 0, NPla = 0, N = 0;
		uint* primitiveIdx;
		class Scene* scene;
		BVHNode* bvhNode; //- 1];
		Mesh* mesh;
		mat4 invTransform;
		aabb bounds;
		class DataCollector* dataCollector;
		int splitMethod;
		
};

struct Bin { aabb bounds; int primCount = 0; };

}

