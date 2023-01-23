#include "precomp.h"
#include "QTable.h"

void QTable::GeneratePoints(Scene& s) {
	tempBounces = maxBounces;
	for(int i = 0; i < emittedRays; i++){
		float3 emittedDir = RandomVectorInUnitSphere();
		Ray emitted = Ray(emitterPos, emittedDir, float3(0));
		Bounce(s, emitted);
		tempBounces = maxBounces;
	}
	s.rebuildBVH();
	cout << "NODES IN TREE: " << kdTree->count << endl;
}

void QTable::Bounce(Scene& s, Ray& emitted) {

	if (tempBounces <= 0) return;
	s.FindNearest(emitted, 0.001f);

	if (emitted.objIdx == -1) return;

	if (emitted.objIdx >= 11 && emitted.objIdx < 11 + size(s.lights)) return;
	if (emitted.GetMaterial()->type == DIFFUSE )
		//&& distance(nearestPoint, emitted.IntersectionPoint() > rejectRadius)
	{
		if (kdTree->getNearestDist(kdTree->rootNode, emitted.IntersectionPoint(), 3) >= rejectRadius) {
			kdTree->insert(NULL, emitted.IntersectionPoint());
			//cout << emitted.IntersectionPoint().x << ", " << emitted.IntersectionPoint().y << ", " << emitted.IntersectionPoint().z << endl;
			s.instantiateDebugPoint(emitted.IntersectionPoint(), emitted.hitNormal);
		}
			
	}
	tempBounces--;
	float3 emittedDir = RandomInHemisphere(emitted.hitNormal);
	Ray bounce = Ray(emitted.IntersectionPoint(), emittedDir, float3(0));
	Bounce(s, bounce);
}