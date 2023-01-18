#pragma once

class HemisphereSampling
{
public:
	struct Sample {
		float3 dir;
		int idx;
		float prob;
	};

	//float3 uniformSample(float u, float v) const;
	virtual void SampleDirection(Sample& s, float3 normal, bool trainingPhase) const;
	float3 uniformSample(float u1, float u2) const {
		float r = std::sqrt(1.0f - u1 * u1);
		float phi = 2 * PI * u2;

		return float3(std::cos(phi) * r, std::sin(phi) * r, u1);
	}
};

class CosineWeightedSampling : public HemisphereSampling {
public:
	float3 cosineWeightedSample(float u, float v) const;
	void SampleDirection(Sample& s, float3 normal, bool trainingPhase) const override;
};

class HemisphereMapping : public HemisphereSampling {
public:
	HemisphereMapping(int resX = 8, int resY = 5) : resX(resX), resY(resY), grid(resX*resY){}

	std::size_t size() {
		return grid.size();
	}

	const float getValue(int idx) {
		return grid[idx];
	}

	float3 getDir(int idx) {
		return normalize(mapIndexToDirection(idx));
	}

	void updateByIndex(int idx, float update) {
		grid[idx] = update;
	}

	void SampleDirection(Sample& s, float3 normal, bool training) const override;
private:
	float3 mapIndexToDirection(int dirIdx) const;
	float3 simpleMap(float x, float y) const;
	int resX, resY;
	std::vector<float> grid;
};
