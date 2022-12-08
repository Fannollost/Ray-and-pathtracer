#pragma once

namespace Tmpl8{
	class Scene;
class bvh
{
	public:
		bvh(Scene* s);
	
		struct BVHNode {
			float3 aabbMin, aabbMax;
			uint leftNode, firstTriIdx, triCount;
			bool isLeaf() { return triCount > 0; };
			uint firstPrim, primCount;
		};
	
		void BuildBVH();
		void UpdateNodeBounds(uint nodeIdx);
		void Subdivide(uint rootNodeIdx);
	public:
		uint rootNodeIdx = 0, nodesUsed = 1, N = 0;
		class Scene* scene;
		BVHNode bvhNode[2 * 20 - 1];
	};
}

