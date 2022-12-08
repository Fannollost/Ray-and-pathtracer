#include "precomp.h"
#include "bvh.h"

bvh::bvh(Scene* s) { scene = s; };

void bvh::BuildBVH() {
	int count = 0;
	for (int i = 0; i < size(scene->tri); ++i) {
		scene->tri[i].centroid = (scene->tri[i].v0 +
				scene->tri[i].v1 + scene->tri[i].v2) * 0.3333f;
			count++;
			scene->triIdx[i] = i;
	}
	N = count;
	BVHNode& root = bvhNode[rootNodeIdx];
	root.triCount = count;
	root.leftNode = root.firstTriIdx = 0;
	root.firstPrim = 0, root.primCount = N;
	UpdateNodeBounds(rootNodeIdx);
	Subdivide(rootNodeIdx);
}

void bvh::UpdateNodeBounds(uint nodeIdx) {
	BVHNode& node = bvhNode[nodeIdx];
	node.aabbMin = float3(1e30f);
	node.aabbMax = float3(-1e30f);
	for (uint first = node.firstTriIdx, i = 0; i < node.triCount; i++) {
		Triangle& leafTri = scene->tri[first + i];
		node.aabbMin = fminf(node.aabbMin, leafTri.v0);
		node.aabbMin = fminf(node.aabbMin, leafTri.v1);
		node.aabbMin = fminf(node.aabbMin, leafTri.v2);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v0);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v1);
		node.aabbMax = fmaxf(node.aabbMax, leafTri.v2);
	}
}

void bvh::Subdivide(uint nodeIdx) {
	BVHNode& node = bvhNode[nodeIdx];
	float3 extent = node.aabbMax - node.aabbMin;
	int axis = 0;
	if (extent.y > extent.x) axis = 1;
	if (extent.z > extent[axis]) axis = 2;
	float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
	int i = node.firstPrim;
	int j = i + node.primCount - 1;
	while (i <= j) {
		if (scene->tri[i].centroid[axis] < splitPos) i++;
		else swap(scene->tri[i], scene->tri[j--]);
	}
	int leftCount = i - node.firstPrim;
	if (leftCount == 0 || leftCount == node.primCount) return;
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	node.leftNode = leftChildIdx;
	bvhNode[leftChildIdx].firstPrim = node.firstPrim;
	bvhNode[leftChildIdx].primCount = leftCount;
	bvhNode[rightChildIdx].firstPrim = i;
	bvhNode[rightChildIdx].primCount = node.primCount - leftCount;
	node.primCount = 0;

}