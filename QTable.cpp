#include "precomp.h"
#include "QTable.h"

void QTable::GeneratePoints(const Scene& s) {
	tempBounces = maxBounces;
	for(int i = 0; i < emittedRays; i++){
		float3 emittedDir = RandomVectorInUnitSphere();
		Ray emitted = Ray(emitterPos, emittedDir, float3(0));
		Bounce(s, emitted);
		tempBounces = maxBounces;
	}
	cout << "NODES IN TREE: " << kdTree->count << endl;
}

void QTable::Bounce(const Scene& s, Ray& emitted) {

	if (tempBounces <= 0) return;
	s.FindNearest(emitted, 0.001f);

	if (emitted.objIdx == -1) return;

	if (emitted.objIdx >= 11 && emitted.objIdx < 11 + size(s.lights)) return;
	if (emitted.GetMaterial()->type == DIFFUSE )
		//&& distance(nearestPoint, emitted.IntersectionPoint() > rejectRadius)
	{
		cout << emitted.IntersectionPoint().x << ", " << emitted.IntersectionPoint().y << ", " << emitted.IntersectionPoint().z << endl;
		if(kdTree->getNearestDist(kdTree->rootNode,emitted.IntersectionPoint(), 3) >= rejectRadius)
			kdTree->insert(NULL, emitted.IntersectionPoint());
	}
	tempBounces--;
	float3 emittedDir = RandomVectorInUnitSphere();
	Ray bounce = Ray(emitted.IntersectionPoint(), emittedDir, float3(0));
	Bounce(s, bounce);
}