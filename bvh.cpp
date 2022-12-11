#include "precomp.h"
#include "bvh.h"
//#define USE_SSE			//FASTER without? very weird

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
	root.leftFirst = 0;
	//root.firstPrim = 0, root.primCount = N;
	UpdateNodeBounds(rootNodeIdx);
	Subdivide(rootNodeIdx);
}

void bvh::UpdateNodeBounds(uint nodeIdx) {
	BVHNode& node = bvhNode[nodeIdx];
	node.aabbMin = float3(1e30f);
	node.aabbMax = float3(-1e30f);
	for (uint first = node.leftFirst, i = 0; i < node.triCount; i++) {
		uint leafTriIdx = scene->triIdx[first + i];
		Triangle& leafTri = scene->tri[leafTriIdx];
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
	// determine split axis using SAH
	int bestAxis = -1;
	float bestPos = 0, bestCost = 1e30f;
	for (int axis = 0; axis < 3; axis++) for (uint i = 0; i < node.triCount; i++)
	{
		Triangle& triangle = scene->tri[scene->triIdx[node.leftFirst + i]];
		float candidatePos = triangle.centroid[axis];
		float cost = EvaluateSAH(node, axis, candidatePos);
		if (cost < bestCost)
			bestPos = candidatePos, bestAxis = axis, bestCost = cost;
	}
	int axis = bestAxis;
	float splitPos = bestPos;
	float3 e = node.aabbMax - node.aabbMin; // extent of parent
	float parentArea = e.x * e.y + e.y * e.z + e.z * e.x;
	float parentCost = node.triCount * parentArea;
	if (bestCost >= parentCost) return;
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		if (scene->tri[scene->triIdx[i]].centroid[axis] < splitPos)
			i++;
		else
			swap(scene->triIdx[i], scene->triIdx[j--]);
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == node.triCount) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	bvhNode[leftChildIdx].leftFirst = node.leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	UpdateNodeBounds(leftChildIdx);
	UpdateNodeBounds(rightChildIdx);
	// recurse
	Subdivide(leftChildIdx);
	Subdivide(rightChildIdx);
}

float bvh::EvaluateSAH(BVHNode& node, int axis, float pos)
{
	// determine triangle counts and bounds for this split candidate
	aabb leftBox, rightBox;
	int leftCount = 0, rightCount = 0;
	for (uint i = 0; i < node.triCount; i++)
	{
		Triangle& triangle = scene->tri[scene->triIdx[node.leftFirst + i]];
		if (triangle.centroid[axis] < pos)
		{
			leftCount++;
			leftBox.Grow(triangle.v0);
			leftBox.Grow(triangle.v1);
			leftBox.Grow(triangle.v2);
		}
		else
		{
			rightCount++;
			rightBox.Grow(triangle.v0);
			rightBox.Grow(triangle.v1);
			rightBox.Grow(triangle.v2);
		}
	}
	float cost = leftCount * leftBox.Area() + rightCount * rightBox.Area();
	return cost > 0 ? cost : 1e30f;
}
void bvh::IntersectBVH(Ray& ray) {
	float t_min = 0.0001f;
	BVHNode* node = &bvhNode[rootNodeIdx], *stack[64];
	uint stackPtr = 0;
	while(1){
		//if (!IntersectAABB(ray, node->aabbMin, node->aabbMax)) return;
		if (node->triCount > 0) {
			for (uint i = 0; i < node->triCount; i++) {
				scene->tri[scene->triIdx[node->leftFirst + i]].Intersect(ray, t_min);
			}
			if (stackPtr == 0) break; else node = stack[--stackPtr];
			continue;
		}

		BVHNode* c1 = &bvhNode[node->leftFirst];
		BVHNode* c2 = &bvhNode[node->leftFirst + 1];
#ifdef USE_SSE
		float dist1 = IntersectAABB_SSE(ray, c1->aabbMin4, c1->aabbMax4);
		float dist2 = IntersectAABB_SSE(ray, c2->aabbMin4, c2->aabbMax4);
#else
		float dist1 = IntersectAABB(ray, c1->aabbMin, c1->aabbMax);
		float dist2 = IntersectAABB(ray, c2->aabbMin, c2->aabbMax);
#endif
		if (dist1 > dist2) { swap(dist1, dist2); swap(c1, c2); }
		if (dist1 == 1e30f) {
			if (stackPtr == 0) break; else node = stack[--stackPtr];
		}
		else {
			node = c1;
			if (dist2 != 1e30f) stack[stackPtr++] = c2;
		}
	}
}


float bvh::IntersectAABB_SSE(const Ray& ray, const __m128 bmin4, const __m128 bmax4)
{
	static __m128 mask4 = _mm_cmpeq_ps(_mm_setzero_ps(), _mm_set_ps(1, 0, 0, 0));
	__m128 t1 = _mm_mul_ps(_mm_sub_ps(_mm_and_ps(bmin4, mask4), ray.O4), ray.rD4);
	__m128 t2 = _mm_mul_ps(_mm_sub_ps(_mm_and_ps(bmax4, mask4), ray.O4), ray.rD4);
	__m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
	float tmax = min(vmax4.m128_f32[0], min(vmax4.m128_f32[1], vmax4.m128_f32[2]));
	float tmin = max(vmin4.m128_f32[0], max(vmin4.m128_f32[1], vmin4.m128_f32[2]));
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}

float bvh::IntersectAABB(const Ray& ray, const float3 bmin, const float3 bmax)
{
	float tx1 = (bmin.x - ray.O.x) * ray.rD.x, tx2 = (bmax.x - ray.O.x) * ray.rD.x;
	float tmin = min(tx1, tx2), tmax = max(tx1, tx2);
	float ty1 = (bmin.y - ray.O.y) * ray.rD.y, ty2 = (bmax.y - ray.O.y) * ray.rD.y;
	tmin = max(tmin, min(ty1, ty2)), tmax = min(tmax, max(ty1, ty2));
	float tz1 = (bmin.z - ray.O.z) * ray.rD.z, tz2 = (bmax.z - ray.O.z) * ray.rD.z;
	tmin = max(tmin, min(tz1, tz2)), tmax = min(tmax, max(tz1, tz2));
	if (tmax >= tmin && tmin < ray.t && tmax > 0) return tmin; else return 1e30f;
}