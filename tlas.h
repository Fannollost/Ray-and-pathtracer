
#pragma once

struct TLASNode
{
    float3 aabbMin;
    uint leftRight; // 2x16 bits
    float3 aabbMax;
    uint BLAS;
    bool isLeaf() { return leftRight == 0; }
};

class tlas
{
public:
    tlas(bvhInstance* bvhList, int N);

    void tlas::build();
    int tlas::FindBestMatch(int* list, int N, int A);
    void Intersect(Ray& ray,  bool debug = false);
    bool IsOccluded(Ray& ray);
public:
    TLASNode* tlasNode;
    uint nodesUsed = 0;
    bvhInstance* blas;
    uint blasCount;

    
};
