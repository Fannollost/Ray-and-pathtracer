#pragma once

namespace Tmpl8{
	class bvh;
	struct BVHNode;
	struct aabb;
class DataCollector
{
	public:
		DataCollector();
		void ResetDataCollector();
		void UpdateNodeCount(int nc);
		void UpdateSummedArea(float3 aabbMin, float3 aabbMax);
		void UpdateAverageTraversalSteps(float ats);
		void UpdateIntersectedPrimitives(float ip);
		void UpdateTreeDepth(bool isLeaf);
		int CalculateDepth(BVHNode& node) {}
		void ExportResults();
		int GetTreeDepth() { return maxTreeDepth; }
		int GetSummedNodeArea() { return summedNodeArea; }
	private:
		int nodeCount;
		float summedNodeArea;
		float averageTraversalStepsPerRay;
		int intersectedPrimitiveCountPerIteration;
		int maxTreeDepth = 0, currDepth = 0;

	//nodeCount, summed node area, traversal steps, intersected primitive count, tree depth;
};
}

