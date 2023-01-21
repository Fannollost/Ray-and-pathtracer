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
		for (int i = 0; i < grid.size(); i++)
			sum += grid[i]; 

		for (int i = 0; i < grid.size(); i++){																	   //FIX DIT! 
			//if(grid[i] != 0) cout << "WTF:" << grid[i] << " SUM " << sum << " DIE DING " << grid[i] / sum << endl;
			normalizedGrid[i] = grid[i] / sum;
		}

		std::vector<float> cum;
		float totSum = 0.0f;
		for (size_t i = 0; i < grid.size(); i++)
		{
			totSum += normalizedGrid[i];
			cum.push_back(totSum);
		}

		//std::partial_sum(normalizedGrid.begin(), normalizedGrid.end(), normalizedGrid.begin());   //WHY IS THIS 1???????????????????????????
		s.idx = grid.size() - 1;
		for (int i = 0; i < (int)cum.size(); i++) {	   
			//r = RandomFloat();  // welicht eruit gooien?
				if (r > cum[i]) {
					s.idx = i;
					s.prob = normalizedGrid[s.idx]  * grid.size() * INV2PI;

					break;
				}
		}

		s.dir = normalize(mapIndexToDirection(s.idx));


}

float3 HemisphereMapping::mapIndexToDirection(int dirIdx) const {
	int idxX = dirIdx % resX;
	int idxY = dirIdx / resX;

	float sizeX = 1.0f / (float)resX;
	float sizeY = 1.0f / (float)resY;

	float randomShiftX = RandomFloat();
	float randomShiftY = RandomFloat();

	float x = ((float)idxX + randomShiftX) * sizeX;
	float y = ((float)idxY + randomShiftY) * sizeY;

	return simpleMap(x, y);
}

float3 HemisphereMapping::simpleMap(float x, float y) const {
	return HemisphereSampling::uniformSample(x, y);
}

