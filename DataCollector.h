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
		void UpdateAverageTraversalSteps(int ats);
		void UpdateIntersectedPrimitives();
		void UpdateTreeDepth(bool isLeaf);
		void UpdateBuildTime(float bt);
		void UpdateFPS(float fps);
		int CalculateDepth(BVHNode& node) {}

		float GetAverageTraversalSteps(int frameNr);
		float GetIntersectedPrimitives(int frameNumber);
		float GetAverageFPS(int frameNumber);
		float GetBuildTime() { return bvhBuildTime; }
		int GetNodeCount() { return nodeCount; }
		int GetTreeDepth() { return maxTreeDepth; }
		int GetSummedNodeArea() { return summedNodeArea; }
	

	private:
		int nodeCount;
		float summedNodeArea;
		long traversalStepsPerIteration;
		long intersectedPrimitiveCountPerIteration;
		float bvhBuildTime;
		float averageFPS;
		float averagePrimitivePerScreen = 0;
		float averageTraversalStepsPerScreen = 0;
		int maxTreeDepth = 0, currDepth = 0;

	//nodeCount, summed node area, traversal steps, intersected primitive count, tree depth;
};
}

