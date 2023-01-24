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
			:lr(lr), emitterPos(emitterPos), maxBounces(maxBounces), rejectRadius(rejectRadius) {
			tableSize = 40;
		}
		void GeneratePoints(Scene& s);
		void Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, const Ray& r, float3 BRDF, Scene& s);
		void SampleDirection(const int i, HemisphereMapping::Sample& s);
		void exportQTable(string exportFile);
		void parseQTable(string path, Scene& s);
		void writeQTable(string exportFile, Tmpl8::KDTree::Node* node);
	
		KDTree* kdTree = new KDTree();
		bool trainingPhase = true;
		int tableSize;
		float3 nextDir;
	private:
		void Bounce(Scene& s, Ray& r);
		float ApproxIntegral(const int idx, const float3& w, const Ray& r, float3 BRDF);
		std::map<int, HemisphereMapping> table;   //Int is placeholder!
		float3 emitterPos;
		int tempBounces, maxBounces, emittedRays = 20;
		float lr, rejectRadius, explorationRate = 0.2f;
	};
}

