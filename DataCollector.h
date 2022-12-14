#pragma once

namespace Tmpl8{
	class bvh;
	struct BVHNode;
class DataCollector
{
	public:
		DataCollector();
		void ResetDataCollector();
		void UpdateNodeCount(int nc);
		void UpdateSummedArea(float sna);
		void UpdateAverageTraversalSteps(float ats);
		void UpdateIntersectedPrimitives(float ip);
		void UpdateTreeDepth(const bvh& b);
		int CalculateDepth(BVHNode& node) {}
		void ExportResults();
	private:
		int nodeCount;
		float summedNodeArea;
		float averageTraversalStepsPerRay;
		int intersectedPrimitiveCountPerIteration;
		int treeDepth;

	//nodeCount, summed node area, traversal steps, intersected primitive count, tree depth;
};
}

