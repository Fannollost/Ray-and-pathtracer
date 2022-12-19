#include "precomp.h"

void bvhInstance::Intersect(Ray& ray)
{
    // backup ray and transform original
    Ray backupRay = ray;
    ray.O = TransformPosition(ray.O, invTransform);
    ray.D = TransformVector(ray.D, invTransform);
    ray.rD = float3(1 / ray.D.x, 1 / ray.D.y, 1 / ray.D.z);
    // trace ray through BVH
    bvh->Intersect(ray);
    backupRay.m = ray.m;
    backupRay.t = ray.t;
    backupRay.objIdx = ray.objIdx;
    backupRay.hitNormal = normalize(TransformVector(ray.hitNormal, matTransform));
    //restore ray origin and direction
    ray = backupRay;
}

void bvhInstance::SetTransform(mat4& transform) {
    invTransform = transform.Inverted();
    matTransform = transform;
    // calculate world-space bounds using the new matrix
    float3 bmin = bvh->bounds.bmin, bmax = bvh->bounds.bmax;
    for (int i = 0; i < 8; i++)
        bounds.grow(TransformPosition(float3(i & 1 ? bmax.x : bmin.x, i & 2 ? bmax.y : bmin.y, i & 4 ? bmax.z : bmin.z), transform));
}