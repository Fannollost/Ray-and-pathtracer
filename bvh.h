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
	
		void BuildBVH();
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint rootNodeIdx);
		void IntersectBVH(Ray& ray);
		float IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax);
		float IntersectAABB_SSE(const Ray& ray, const __m128 bmin4, const __m128 bmax4);
		float EvaluateSAH(BVHNode &node, int axis, float pos);
	public:
		uint rootNodeIdx = 0, nodesUsed = 2, N = 12583;
		class Scene* scene;
		BVHNode bvhNode[2 * 12583 - 1];
	};
}

