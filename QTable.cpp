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
		if(kdTree->getNearestDist(kdTree->rootNode,emitted.IntersectionPoint(), 3) >= rejectRadius)
			kdTree->insert(NULL, emitted.IntersectionPoint());
	}
	tempBounces--;
	float3 emittedDir = RandomInHemisphere(emitted.hitNormal);
	Ray bounce = Ray(emitted.IntersectionPoint(), emittedDir, float3(0));
	Bounce(s, bounce);
}

void QTable::Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, const material& m) {
	int idx = kdTree->findNearest(kdTree->rootNode, origin, 3);
	HemisphereMapping& value = table.at(idx);
	float val = value.getValue(wIndex);
	float3 dir = value.getDir(wIndex);

	float qUpdate = (1.0f - lr) * val + lr * (length(irradiance) + ApproxIntegral(hitPoint, dir, m));
	value.updateByIndex(wIndex, qUpdate);
}

float QTable::ApproxIntegral(const float3 hitPoint, const float3& w, const material& m) {
	auto r = table.insert({1, HemisphereMapping(resx,resy) });  ////????????
	auto& mapping = r.first->second;

	float sum = 0.0f;
	float3 normal(0, 0, 1.0f); //need to fix;

	for (int i = 0; i < (int)mapping.size(); i++) {
		auto Qy = mapping.getValue(i);
		float3 wi = mapping.getDir(i);
		float cosAngle = max(dot(wi, normal), 0.f);
		float3 fs = 0.7f; //--> BRDF calculate (pass as arg?)
		sum += length(fs) * Qy * cosAngle;

	}

	return sum / (float)mapping.size();
}

