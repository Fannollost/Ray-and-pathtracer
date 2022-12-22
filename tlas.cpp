#include "precomp.h"

tlas::tlas(bvhInstance* bvhList, int N)
{
    // copy a pointer to the array of bottom level accstruc instances
    blas = bvhList;
    blasCount = N;
    // allocate TLAS nodes
    tlasNode = (TLASNode*)_aligned_malloc(sizeof(TLASNode) * 2 * N, 64);
    nodesUsed = 2;
}

void tlas::build()
{
    // assign a TLASleaf node to each BLAS
    int nodeIdx[256], nodeIndices = blasCount;
    nodesUsed = 1;
    for (uint i = 0; i < blasCount; i++)
    {
        nodeIdx[i] = nodesUsed;
        tlasNode[nodesUsed].aabbMin = blas[i].bounds.bmin;
        tlasNode[nodesUsed].aabbMax = blas[i].bounds.bmax;
        tlasNode[nodesUsed].BLAS = i;
        tlasNode[nodesUsed++].leftRight = 0; // makes it a leaf
    }

    // use agglomerative clustering to build the TLAS
    int A = 0, B = FindBestMatch(nodeIdx, nodeIndices, A);
    while (nodeIndices > 1)
    {
        int C = FindBestMatch(nodeIdx, nodeIndices, B);
        if (A == C)
        {
            int nodeIdxA = nodeIdx[A], nodeIdxB = nodeIdx[B];
            TLASNode& nodeA = tlasNode[nodeIdxA];
            TLASNode& nodeB = tlasNode[nodeIdxB];
            TLASNode& newNode = tlasNode[nodesUsed];
            newNode.leftRight = nodeIdxA + (nodeIdxB << 16);
            newNode.aabbMin = fminf(nodeA.aabbMin, nodeB.aabbMin);
            newNode.aabbMax = fmaxf(nodeA.aabbMax, nodeB.aabbMax);
            nodeIdx[A] = nodesUsed++;
            nodeIdx[B] = nodeIdx[nodeIndices - 1];
            B = FindBestMatch(nodeIdx, --nodeIndices, A);
        }
        else A = B, B = C;
    }
    tlasNode[0] = tlasNode[nodeIdx[A]];
}

int tlas::FindBestMatch(int* list, int N, int A)
{
    float smallest = 1e30f;
    int bestB = -1;
    for (int B = 0; B < N; B++) if (B != A)
    {
        float3 bmax = fmaxf(tlasNode[list[A]].aabbMax, tlasNode[list[B]].aabbMax);
        float3 bmin = fminf(tlasNode[list[A]].aabbMin, tlasNode[list[B]].aabbMin);
        float3 e = bmax - bmin;
        float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
        if (surfaceArea < smallest) smallest = surfaceArea, bestB = B;
    }
    return bestB;
}

void tlas::Intersect(Ray& ray)
{
    TLASNode* node = &tlasNode[0], * stack[64];
    uint stackPtr = 0;
    while (1)
    {
        if (node->isLeaf())
        {
            blas[node->BLAS].BIntersect(ray);
            if (stackPtr == 0) break; else node = stack[--stackPtr];
            continue;
        }
        TLASNode* child1 = &tlasNode[node->leftRight & 0x0000FFFF];
        TLASNode* child2 = &tlasNode[node->leftRight >> 16];
        float dist1 = bvh::IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
        float dist2 = bvh::IntersectAABB(ray, child2->aabbMin, child2->aabbMax);
        if (dist1 > dist2) { swap(dist1, dist2); swap(child1, child2); }
        if (dist1 == 1e30f)
        {
            if (stackPtr == 0) break; else node = stack[--stackPtr];
        }
        else
        {
            node = child1;
            if (dist2 != 1e30f) stack[stackPtr++] = child2;
        }
    }
}

bool tlas::IsOccluded(Ray& ray)
{
    TLASNode* node = &tlasNode[0], * stack[64];
    uint stackPtr = 0;
    while (1)
    {
        if (node->isLeaf())
        {
            if(blas[node->BLAS].IsOccluded(ray)) return true;
            if (stackPtr == 0) break; else node = stack[--stackPtr];
            continue;
        }
        TLASNode* child1 = &tlasNode[node->leftRight & 0x0000FFFF];
        TLASNode* child2 = &tlasNode[node->leftRight >> 16];
        float dist1 = bvh::IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
        float dist2 = bvh::IntersectAABB(ray, child2->aabbMin, child2->aabbMax);
        if (dist1 > dist2) { swap(dist1, dist2); swap(child1, child2); }
        if (dist1 == 1e30f)
        {
            if (stackPtr == 0) break; else node = stack[--stackPtr];
        }
        else
        {
            node = child1;
            if (dist2 != 1e30f) stack[stackPtr++] = child2;
        }
    }
    return false;
}