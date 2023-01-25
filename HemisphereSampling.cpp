#include "precomp.h"
#include "HemisphereSampling.h"
#include <numeric>

void HemisphereSampling::SampleDirection(Sample& s, float3 normal, bool trainingPhase) const {
	s.dir = RandomInHemisphere(normal);
	s.prob = INV2PI;
}

float3 CosineWeightedSampling::cosineWeightedSample(float u, float v) const {
	float r = sqrt(u);
	float t = TWOPI * v;

	float x = r * cos(t);
	float y = r * sin(t);
	return float3(x, y, sqrt(max(0.0f, 1 - u)));
}

void CosineWeightedSampling::SampleDirection(Sample& s, float3 normal, bool trainingPhase) const {
	s.dir = cosineWeightedSample(RandomFloat(), RandomFloat());
	float cosAngle = max(dot(float3(0, 0, 1.0f), s.dir), 0.0f);
	s.prob = cosAngle * INVPI;
}


void HemisphereMapping::SampleDirection(Sample& s, float3 normal, bool training) const {
	float r = RandomFloat();
	std::vector<float> normalizedGrid(grid.size());
	//if(!training) {
		float sum = 0.0f;
		for (int i = 0; i < grid.size(); i++) {
			sum += grid[i]; 
		}


		for (int i = 0; i < grid.size(); i++){	
			normalizedGrid[i] = grid[i] / sum;
		}

		std::vector<float> cum;
		float totSum = 0.0f;
		for (size_t i = 0; i < grid.size(); i++)
		{
			totSum += normalizedGrid[i];
			cum.push_back(totSum);
		}

		//cout << cum[39] << endl;
		//std::partial_sum(normalizedGrid.begin(), normalizedGrid.end(), normalizedGrid.begin());   //WHY IS THIS 1???????????????????????????
		s.idx = grid.size() - 1;
		s.prob = 1 / grid.size();
		float f = 0;
		for (int i = 0; i < (int)cum.size(); i++) {	 
			r = RandomFloat();
			if (r < explorationRate && training)
			{
				s.idx = ((int)grid.size() * RandomFloat());
				s.prob = 1 / (grid.size() * 2 * PI);
				break;
			}
			if (r < normalizedGrid[i] && sum != 0) {
				s.idx = i;
				s.prob = normalizedGrid[s.idx];
				break;
			}
		}
		s.dir = mapIndexToDirection(s.idx);
}

float3 HemisphereMapping::mapIndexToDirection(int dirIdx) const {
	float3 vertices[26] = { float3(0.276388,-0.850649,0.44722),float3(0.688189,-0.499997,0.525737),float3(0.162456,-0.499997,0.850653),float3(-0.723608,-0.525724,0.44722),float3(-0.262868,-0.809013,0.525736),float3(-0.425324,-0.309012,0.850653),float3(-0.723608,0.525724,0.44722),float3(-0.850648,0,0.525736),float3(-0.425323,0.309011,0.850654),float3(0.276388,0.850649,0.447221),float3(-0.262868,0.809011,0.525739),float3(0.162456,0.499995,0.850654),float3(0.894426,0,0.447217),float3(0.688189,0.499997,0.525737),float3(0.525731,0,0.850651),float3(0,0,1),float3(0.951058,0.309013,0),float3(0.587788,0.809015,0),float3(0,1,0),float3(-0.587788,0.809015,0),float3(-0.951058,0.309013,0),float3(-0.951058,-0.309013,0),float3(-0.587788,-0.809015,0),float3(0,-1,0),float3(0.587788,-0.809015,0),float3(0.951058,-0.309013,0) };
	int3 faces[40] = { int3(15, 12, 16), int3(12, 9, 16), int3(9, 6, 16), int3(6, 3, 16), int3(3, 15, 16), int3(3, 2, 15), int3(2, 13, 15), int3(13, 14, 15), int3(15, 14, 12), int3(14, 10, 12), int3(10, 11, 12), int3(12, 11, 9), int3(11, 7, 9), int3(7, 8, 9), int3(9, 8, 6), int3(8, 4, 6), int3(4, 5, 6), int3(6, 5, 3), int3(5, 1, 3), int3(1, 2, 3), int3(25, 2, 1), int3(25, 26, 2), int3(26, 13, 2), int3(26, 17, 13), int3(17, 14, 13), int3(17, 18, 14), int3(18, 10, 14), int3(18, 19, 10), int3(19, 11, 10), int3(19, 20, 11), int3(20, 7, 11), int3(20, 21, 7), int3(21, 8, 7), int3(21, 22, 8), int3(22, 4, 8), int3(22, 23, 4), int3(23, 5, 4), int3(23, 24, 5), int3(24, 1, 5), int3(24, 25, 1) };
	
	float randomShiftX = RandomFloat();
	float randomShiftY = RandomFloat();

	float3 v0 = vertices[faces[dirIdx].x];
	float3 v1 = vertices[faces[dirIdx].y];
	float3 v2 = vertices[faces[dirIdx].z];

	float3 res = normalize(v0 + (v0 - v1) * randomShiftX + (v0 - v2) * randomShiftY);
	cout << dirIdx << ": " << res.x << ", " << res.y << ", " << res.z << " length: " << endl;
	return res;
}

