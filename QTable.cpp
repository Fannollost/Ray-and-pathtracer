#include "precomp.h"
#include "QTable.h"

void QTable::GeneratePoints(Scene& s) {
	tempBounces = maxBounces;
	for(int i = 0; i < emittedRays; i++){
		float3 emittedDir = float3(RandomFloat() * 2 - 1, RandomFloat() * 2 - 1, RandomFloat() * 2 - 1);
		Ray emitted = Ray(emitterPos, emittedDir, float3(0));
		Bounce(s, emitted);
		tempBounces = maxBounces;
	}
	//s.rebuildBVH();
	cout << "NODES IN TREE: " << kdTree->count << endl;
}

void QTable::SampleDirection(const int i, HemisphereMapping::Sample& s) {
	auto& r = table.insert({ i,HemisphereMapping(explorationRate) });
	auto& it = r.first;
	it->second.SampleDirection(s, float3(0), trainingPhase);
}

void QTable::Bounce(Scene& s, Ray& emitted) {

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
				KDTree::Node* node = kdTree->insert(kdTree->rootNode, emitted.IntersectionPoint() + emitted.hitNormal * 0.001f, emitted.hitNormal);
				table.insert({ kdTree->lastInsertedIdx,HemisphereMapping(explorationRate) });
				//s.instantiateDebugPoint(emitted.IntersectionPoint() + emitted.hitNormal * 0.001f, emitted.hitNormal, kdTree->count-1);
			}
		}
	}

	tempBounces--;
	float3 emittedDir = RandomInHemisphere(emitted.hitNormal);
	Ray bounce = Ray(emitted.IntersectionPoint() + emittedDir * 0.0001f, emittedDir, float3(0));
	Bounce(s, bounce);
}

void QTable::Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, Ray& r, float3 BRDF, Scene& s) {

	float dist = kdTree->getNearestDist(kdTree->rootNode, origin, 0);
	int idx = kdTree->nearestNode->idx;

	float3 tempIr = irradiance;

	float distHit = kdTree->getNearestDist(kdTree->rootNode, hitPoint, 0);
	int hitIdx = kdTree->nearestNode->idx;

	auto& mapping = table.at(idx); 

	float val = mapping.getValue(wIndex);
	float3 dir = mapping.getDir(wIndex);
	float sbeh = ApproxIntegral(hitIdx, dir, r, BRDF);

	float qUpdate = (1.0f - lr) * val + lr * (length(tempIr) + ApproxIntegral(hitIdx, dir, r, BRDF));
	mapping.updateByIndex(wIndex, qUpdate);
	//s.updateQDebug(idx, wIndex, qUpdate);
}

float QTable::ApproxIntegral(const int idx, const float3& w, const Ray& ray, float3 BRDF) {
	auto& r = table.insert({idx, HemisphereMapping(explorationRate) });
	auto& mapping = r.first->second;

	float sum = 0.0f;

	for (int i = 0; i < (int)mapping.size(); i++) {
		auto Qy = mapping.getValue(i);
		float3 wi = mapping.getDir(i);
		float cosAngle = max(dot(wi, ray.hitNormal), 0.f); 
		sum += length(BRDF) * Qy * cosAngle;
	}

	return sum / (float)mapping.size();
}


void QTable::parseQTable(string path, Scene& s) {
	trainingPhase = false;
	tempBounces = maxBounces;
	string line;
	ifstream file(path, ios::in);
	while (getline(file, line)) {
		int idx;
		float px, py, pz, nx, ny, nz;
		const char* constL = line.c_str();
		sscanf(constL, "%i/%f/%f/%f/%f/%f/%f", &idx, &px, &py, &pz, &nx, &ny, &nz);
		kdTree->insert(kdTree->rootNode, float3(px, py, pz), float3(nx, ny, nz), idx);
		//s.instantiateDebugPoint(float3(px, py, pz), float3(nx, ny, nz), kdTree->lastInsertedIdx);
		for (int i = 0; i < 40; i++) {
			auto& r = table.insert({ idx, HemisphereMapping(explorationRate) });
			float val;
			getline(file, line);
			const char* constL = line.c_str();
			sscanf(constL, "%f", &val);
			table.at(kdTree->lastInsertedIdx).updateByIndex(i, val);
			//s.updateQDebug(kdTree->lastInsertedIdx, i, val);
		}		
	}
	//s.rebuildBVH();
	return;
}



void QTable::writeQTable(string exportFile, KDTree::Node* node) {
	std::ofstream myFile;
	myFile.open(exportFile, std::ios_base::app);

	if (node == NULL) return;

	myFile << KDTree::ToString(node) << "\n";
	for (int i = 0; i < 40; i++) {
		myFile << table.at(node->idx).getValue(i) << "\n";
	}
	myFile.close();
	writeQTable(exportFile, node->left);
	writeQTable(exportFile, node->right);
}

void QTable::exportQTable(string exportFile) {
	std::ofstream myFile;
	string sceneTrainingPath = exportFile;
	myFile.open(sceneTrainingPath, std::ios_base::trunc);
	myFile.close();
	
	writeQTable(sceneTrainingPath, kdTree->rootNode);
}