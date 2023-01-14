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
	virtual void SampleDirection(Sample& s, float3 normal) const;
};

class CosineWeightedSampling : public HemisphereSampling {
public:
	float3 cosineWeightedSample(float u, float v) const;
	void SampleDirection(Sample& s, float3 normal) const override;
};
