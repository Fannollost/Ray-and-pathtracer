#include "precomp.h"
#include "QTable.h"

void QTable::GeneratePoints(const Scene& s) {
	tempBounces = maxBounces;
	for(int i = 0; i < emittedRays; i++){
		float3 emittedDir = float3(RandomFloat() * 2 - 1, RandomFloat() * 2 - 1, RandomFloat() * 2 - 1);
		Ray emitted = Ray(emitterPos, emittedDir, float3(0));
		Bounce(s, emitted);
		tempBounces = maxBounces;
	}
	cout << "NODES IN TREE: " << kdTree->count << endl;
}

void QTable::SampleDirection(const int i, HemisphereMapping::Sample& s) {
	auto r = table.insert({ i,HemisphereMapping(explorationRate,resx,resy) });
	auto& it = r.first;
	it->second.SampleDirection(s, float3(0), trainingPhase);
}

void QTable::Bounce(const Scene& s, Ray& emitted) {

	if (tempBounces <= 0) return;
	s.FindNearest(emitted, 0.001f);

	if (emitted.objIdx == -1) return;

	if (emitted.objIdx >= 11 && emitted.objIdx < 11 + size(s.lights)) return;
	if (emitted.GetMaterial()->type == DIFFUSE )
		//&& distance(nearestPoint, emitted.IntersectionPoint() > rejectRadius)
	{
		float d = kdTree->getNearestDist(kdTree->rootNode, emitted.IntersectionPoint(), 0);
		//cout << d << endl;
		KDTree::Node *nearestNode = kdTree->nearestNode;
		float weight;
		if (nearestNode == nullptr) weight = 1; else weight = dot(emitted.hitNormal, nearestNode->normal);
		if(d >= rejectRadius * weight)													 {
			kdTree->insert(kdTree->rootNode, emitted.IntersectionPoint(), emitted.hitNormal);
			if (emitted.objIdx == 7) cout << "sampling point on wahed cirkel " << emitted.hitNormal.x << ", " << emitted.hitNormal.y << ", " << emitted.hitNormal.z << endl;
		}
	}
	tempBounces--;
	float3 emittedDir = RandomInHemisphere(emitted.hitNormal);
	Ray bounce = Ray(emitted.IntersectionPoint(), emittedDir, float3(0));
	Bounce(s, bounce);
}

void QTable::Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, const Ray& r, float3 BRDF) {

	float dist = kdTree->getNearestDist(kdTree->rootNode, origin, 0);
	int idx = kdTree->nearestNode->idx;
	
	float distHit = kdTree->getNearestDist(kdTree->rootNode, hitPoint, 0);
	int hitIdx = kdTree->nearestNode->idx;

	auto mapping = table.at(idx); 
	float val = mapping.getValue(wIndex);
	float3 dir = mapping.getDir(wIndex);
	
	float qUpdate = (1.0f - lr) * val + lr * (length(irradiance) + ApproxIntegral(hitIdx, dir, r, BRDF));
	mapping.updateByIndex(wIndex, qUpdate);
}

float QTable::ApproxIntegral(const int idx, const float3& w, const Ray& ray, float3 BRDF) {
	auto r = table.insert({idx, HemisphereMapping(explorationRate, resx,resy) });
	auto& mapping = r.first->second;

	float sum = 0.0f;

	for (int i = 0; i < (int)mapping.size(); i++) {
		auto Qy = mapping.getValue(i);
		float3 wi = mapping.getDir(i);
		float cosAngle = max(dot(wi, ray.hitNormal), 0.f); ///////////// EVALUATE BRDF IN Wi DIRECTION!
		//float3 fs = 0.7f; //--> BRDF calculate (pass as arg?)
		sum += length(BRDF) * Qy * cosAngle;
	}

	return sum / (float)mapping.size();
}

