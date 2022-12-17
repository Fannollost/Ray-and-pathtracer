#include "precomp.h"
#include "DataCollector.h"


DataCollector::DataCollector() {
	ResetDataCollector();
}

void DataCollector::ResetDataCollector() {
	DataCollector::nodeCount = 0;
	DataCollector::summedNodeArea = 0;
	DataCollector::averageTraversalStepsPerRay = 0;
	DataCollector::intersectedPrimitiveCountPerIteration = 0;
	DataCollector::maxTreeDepth = 0;
	DataCollector::currDepth = 0;
}

void DataCollector::UpdateNodeCount(int nc) {
	nodeCount += nc;
}

void DataCollector::UpdateSummedArea(float3 aabbMin, float3 aabbMax) {
	float3 diagnal = aabbMax - aabbMin;
	summedNodeArea += diagnal.x * diagnal.y + diagnal.y * diagnal.z + diagnal.z * diagnal.x;
}

void DataCollector::UpdateAverageTraversalSteps(float ats) {

}

void DataCollector::UpdateIntersectedPrimitives(float ip) {

}


void DataCollector::UpdateTreeDepth(bool isLeaf) {
	currDepth++;
	if (isLeaf) {
		maxTreeDepth = currDepth > maxTreeDepth ? currDepth : maxTreeDepth;
		currDepth = 0;
	}
	/*BVHNode* node = &b.bvhNode[b.rootNodeIdx], * stack[64];
	uint stackPtr = 0;
	int depth = 0, maxDepth = 0;
	while (1) {
		if (node->isLeaf()) {
			if(stackPtr == 0) treeDepth = maxDepth;
			else {
				node = stack[--stackPtr];
				maxDepth = depth > maxDepth ? depth : maxDepth;
			}
		}	

		BVHNode* c1 = &b.bvhNode[node->leftFirst];
		BVHNode* c2 = &b.bvhNode[node->leftFirst + 1];
	} */

}

void DataCollector::ExportResults() {

}
