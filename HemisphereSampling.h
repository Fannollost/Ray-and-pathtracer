#pragma once

class HemisphereSampling
{
public:
	struct Sample {
		float3 dir;
		int idx;
		float prob;
	};
	virtual void SampleDirection(Sample& s, float3 normal, bool trainingPhase) const;
};

class CosineWeightedSampling : public HemisphereSampling {
public:
	float3 cosineWeightedSample(float u, float v) const;
	void SampleDirection(Sample& s, float3 normal, bool trainingPhase) const override;
};

class HemisphereMapping : public HemisphereSampling {
public:
	HemisphereMapping(float er = 0.2f) : grid(40, ((float)1/40)), explorationRate(er){}

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
	std::vector<float> grid;
	float explorationRate = 0.05f;
};
