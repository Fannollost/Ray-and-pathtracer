#include "precomp.h"
#include "HemisphereSampling.h"


void HemisphereSampling::SampleDirection(Sample& s, float3 normal) const {
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

void CosineWeightedSampling::SampleDirection(Sample& s, float3 normal) const {
	s.dir = cosineWeightedSample(RandomFloat(), RandomFloat());
	float cosAngle = max(dot(float3(0, 0, 1.0f), s.dir), 0.0f);
	s.prob = cosAngle * INVPI;
}