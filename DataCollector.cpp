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
	DataCollector::treeDepth = 0;
}

void DataCollector::UpdateNodeCount(int nc) {
	nodeCount += nc;
}

void DataCollector::UpdateSummedArea(float sna) {

}

void DataCollector::UpdateAverageTraversalSteps(float ats) {

}

void DataCollector::UpdateIntersectedPrimitives(float ip) {

}


void DataCollector::UpdateTreeDepth(const bvh& b) {
	BVHNode* node = &b.bvhNode[b.rootNodeIdx], * stack[64];
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
	}

}

void DataCollector::ExportResults() {

}
