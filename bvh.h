#pragma once

namespace Tmpl8{
	class Scene;
	class Ray;
class bvh
{
	public:
		bvh(Scene* s);
	
		struct BVHNode
		{
			union
			{
				struct { float3 aabbMin; uint leftFirst; };
				__m128 aabbMin4;
			};
			union
			{
				struct { float3 aabbMax; uint triCount; };
				__m128 aabbMax4;
			};
			bool isLeaf() { return triCount > 0; }
		};
	
		void Build();
		void UpdateNodeBounds(uint nodeIdx);
		float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
		float CalculateNodeCost(BVHNode& node);
		void Subdivide(uint rootNodeIdx);
		void Intersect(Ray& ray);
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		float IntersectAABB_SSE(const Ray& ray, const __m128 bmin4, const __m128 bmax4);
		float EvaluateSAH(BVHNode &node, int axis, float pos);
	public:
		uint rootNodeIdx = 0, nodesUsed = 2, N = 12582;
		uint* triIdx;
		class Scene* scene;
		BVHNode bvhNode[2 * 12582 - 1];
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

	struct Bin { aabb bounds; int triCount = 0; };
}

