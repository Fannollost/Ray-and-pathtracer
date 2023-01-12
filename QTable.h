#pragma once
#include <map>
#include "bvh.h"
namespace Tmpl8 {

	struct BVHNode;
	class Scene;
	class KDTree;
	class QTable
	{
	public:
		QTable(int x, int y, float lr, float3 emitterPos, int maxBounces, float rejectRadius) 
			: resx(x), resy(y), lr(lr), emitterPos(emitterPos), maxBounces(maxBounces), rejectRadius(rejectRadius) {
		}

		void GeneratePoints(const Scene& s);
		void Update(const BVHNode* nOrig, const BVHNode* nHit, int wIndex, const float3& irradiance, const Ray& r);
		
	
	private:
		void Bounce(const Scene& s, Ray& r);

		KDTree* kdTree = new KDTree();
		std::map<const BVHNode*, int> table;   //Int is placeholder!
		float3 emitterPos;
		int resx, resy, tempBounces, maxBounces, emittedRays = 10;
		float lr, rejectRadius;
	};
}

