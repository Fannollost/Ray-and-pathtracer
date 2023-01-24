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
	auto& r = table.insert({ i,HemisphereMapping(explorationRate,resx,resy) });
	auto& it = r.first;
	it->second.SampleDirection(s, float3(0), trainingPhase);
}

void QTable::Bounce(const Scene& s, Ray& emitted) {

	if (tempBounces <= 0) return;
	s.FindNearest(emitted, 0.001f);

	if (emitted.objIdx == -1) return;

	
	if (!(emitted.objIdx >= 11 && emitted.objIdx < 11 + size(s.lights))) {
		if (emitted.GetMaterial()->type == DIFFUSE)
		{
			float d = kdTree->getNearestDist(kdTree->rootNode, emitted.IntersectionPoint(), 0);
			KDTree::Node* nearestNode = kdTree->nearestNode;
			float weight;
			if (nearestNode == nullptr) weight = 1; else weight = dot(emitted.hitNormal, nearestNode->normal);
			if (d >= rejectRadius * weight) {
				kdTree->insert(kdTree->rootNode, emitted.IntersectionPoint(), emitted.hitNormal);
			}
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

	auto& mapping = table.at(idx); 

	float val = mapping.getValue(wIndex);
	float3 dir = mapping.getDir(wIndex);
	
	float qUpdate = (1.0f - lr) * val + lr * (length(irradiance) + ApproxIntegral(hitIdx, dir, r, BRDF));
	mapping.updateByIndex(wIndex, qUpdate);
}

float QTable::ApproxIntegral(const int idx, const float3& w, const Ray& ray, float3 BRDF) {
	auto& r = table.insert({idx, HemisphereMapping(explorationRate, resx,resy) });
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


QTable* QTable::parseQTable(string path) {
	QTable* res = new QTable(8, 5, 0.25f, float3(0, 0, 0), 5, 0.2f);
	string line;
	ifstream file(path, ios::in);
	while (getline(file, line)) {
		int idx;
		float px, py, pz, nx, ny, nz;
		const char* constL = line.c_str();
		sscanf(constL, "%i/%f/%f/%f/%f/%f/%f", &idx, &px, &py, &pz, &nx, &ny, &nz);
		Tmpl8::KDTree::Node* node = kdTree->insert(kdTree->rootNode, float3(px, py, pz), float3(nx, ny, nz));
		for (int i = 0; i < 40; i++) {
			auto& r = res->table.insert({ idx, HemisphereMapping(explorationRate, resx,resy) });
			float val;
			getline(file, line);
			const char* constL = line.c_str();
			sscanf(constL, "%f", &val);
			res->table.at(node->idx).updateByIndex(i, val);
		}
		
		
	}
	return res;
}

string ToString(Tmpl8::KDTree::Node* node) {
	return node->idx + "/" + to_string(node->point.x) + "/" + to_string(node->point.y) + "/" + to_string(node->point.z) + "/" + to_string(node->normal.x) + "/" + to_string(node->normal.y) + "/" + to_string(node->normal.z);
}

void QTable::writeQTable(string exportFile, Tmpl8::KDTree::Node* node) {
	std::ofstream myFile(exportFile);

	if (node == NULL) return;

	myFile << ToString(node) << "\n";
	for (int i = 0; i < 40; i++) {
		myFile << table.at(node->idx).getValue(i) << "\n";
	}
	writeQTable(exportFile, node->left);
	writeQTable(exportFile, node->right);
}

void QTable::ToString(string exportFile) {
	writeQTable(exportFile, kdTree->rootNode);
}