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
		QTable(int x, int y, float lr, float3 emitterPos, int maxBounces, float rejectRadius) 
			: resx(x), resy(y), lr(lr), emitterPos(emitterPos), maxBounces(maxBounces), rejectRadius(rejectRadius) {
		}

		void GeneratePoints(const Scene& s);
		void Update(const float3 origin, const float3 hitPoint, int wIndex, const float3& irradiance, const material& m);
		
	
	private:
		void Bounce(const Scene& s, Ray& r);
		float ApproxIntegral(const float3 hitPoint, const float3& w, const material& m);
		KDTree* kdTree = new KDTree();
		std::map<int, HemisphereMapping> table;   //Int is placeholder!
		float3 emitterPos;
		int resx, resy, tempBounces, maxBounces, emittedRays = 10;
		float lr, rejectRadius;
	};
}

