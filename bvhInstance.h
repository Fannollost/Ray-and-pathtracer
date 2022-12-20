
#pragma once

// instance of a BVH, with transform and world bounds
class bvhInstance
{
public:
    bvhInstance() = default;
    bvhInstance(bvh* blas) : bvh(blas) { SetTransform(mat4()); bounds = bvh->bounds; }
    void SetTransform(mat4& transform);
    void Intersect(Ray& ray);
private:
    mat4 invTransform; // inverse transform
    mat4 matTransform;
public:
    bvh* bvh = 0;
    Tmpl8::aabb bounds; // in world space
};

