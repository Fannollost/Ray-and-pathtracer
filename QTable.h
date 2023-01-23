#pragma once
#include <map>
#include "bvh.h"
namespace Tmpl8 {

	struct BVHNode;
	class material;
	class Scene;
	class KDTree;
	class QTable
	{
	public:
		QTable(int x, int y, float lr, float3 emitterPos, int maxBounces, float rejectRadius, float explorationRate = 0.2f) 
			: resx(x), resy(y), lr(lr), emitterPos(emitterPos), maxBounces(maxBounces), rejectRadius(rejectRadius) {
			tableSize = resx * resy;
		}
		void GeneratePoints(const Scene& s);
		void Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, const Ray& r, float3 BRDF);
		void SampleDirection(const int i, HemisphereMapping::Sample& s);
		void ToString(string exportFile);
		string ToString(Tmpl8::KDTree::Node* node);
		QTable* parseQTable(string path);
		void writeQTable(string exportFile, Tmpl8::KDTree::Node* node);
	
		KDTree* kdTree = new KDTree();
		bool trainingPhase = true;
		int tableSize;
		float3 nextDir;
	private:
		void Bounce(const Scene& s, Ray& r);
		float ApproxIntegral(const int idx, const float3& w, const Ray& r, float3 BRDF);
		std::map<int, HemisphereMapping> table;   //Int is placeholder!
		float3 emitterPos;
		int resx, resy, tempBounces, maxBounces, emittedRays = 1000;
		float lr, rejectRadius, explorationRate = 0.2f;
	};
}

