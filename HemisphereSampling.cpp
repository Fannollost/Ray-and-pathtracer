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
	if(!training || training){
		std::vector<float> normalizedGrid(grid.size());
		float sum = 0.0f;
		for (std::size_t i = 0; i < grid.size(); i++)
			sum += grid[i];
		for (std::size_t i = 0; i < grid.size(); i++)
			normalizedGrid[i] = grid[i] / sum;

		std::vector<float> cum(grid.size());
		std::partial_sum(normalizedGrid.begin(), normalizedGrid.end(), cum.begin());
		s.idx = grid.size() - 1;
		for (int i = 0; i < (int)cum.size(); i++) {
			if (r < cum[i]) {
				s.idx = i;
				break;
			}
		}
		s.dir = normalize(mapIndexToDirection(s.idx));
		s.prob = ((float)grid.size() * normalizedGrid[s.idx]) * INV2PI;	
	}
	/*else {
		s.idx = ((int)(resX * resY) * RandomFloat());
		s.prob = 1 / (resX * resY) * RandomFloat();
	} */

}

float3 HemisphereMapping::mapIndexToDirection(int dirIdx) const {
	int idxX = dirIdx % resX;
	int idxY = dirIdx / resX;

	float sizeX = 1.0f / (float)resX;
	float sizeY = 1.0f / (float)resY;

	//float randomShiftX = RandomFloat();
	//float randomShiftY = RandomFloat();

	float x = ((float)idxX * sizeX); //+ randomShiftX) * sizeX;
		float y = ((float)idxY * sizeY);//+ randomShiftY) * sizeY;

	return simpleMap(x, y);
}

float3 HemisphereMapping::simpleMap(float x, float y) const {
	return HemisphereSampling::uniformSample(x, y);
}

