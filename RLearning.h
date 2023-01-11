#pragma once
#include <map>
#include <array>
#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>

const int dim_action_space = 72;
using Key = std::array<float, 6>;
using QValue = std::array<float, dim_action_space + 1>;
using StateAction = std::array<float, 7>;
using StateActionCount = float;

struct StatesStruct {
	std::array<float, 6> oldState;
	int oldAction;
	float probability;
	int oldId;
};

class RLearning
{
public:
	int training, test, ActionSpaceMode;
	virtual void UpdateQTable(Key& state, Key& nextState, Ray& hit, std::map<Key, QValue>* dict, std::map<int, float3>* dictAction,
		int& old_action, float& BRDF, float3& nl, float3& x, float prob, std::map<StateAction, StateActionCount>* dictStateActionCount) = 0;
	virtual float3 SamplingScattering(std::map<Key, QValue>* dict, std::map<int, float3>* dictAction, int& id, float3& x, float3& nl,
		StatesStruct& states_rec, std::map<StateAction, StateActionCount>* dictStateActionCount, float& epsilon, const Scene& scene);
	virtual std::map<Key, QValue>* LoadWeights(string fileName);
	virtual void ComputeQProb(const int& action, const std::array<float, dim_action_space + 1>& qvalue, StatesStruct& states_rec, const float& total);
	virtual void SampleSphereCoordinates(float3& sphereCoord, const float3& pointOldCoord);
	virtual void UpdateStateActionCount(const Key& state, const int& action, std::map<StateAction, StateActionCount>* dictStateActionCount,
		bool sequence_complete);

	float epsilon, pathLength = 0;
	float* ptrPathLength = &pathLength;
	int counter_states = 0;
	int counter_red = 0;
	int counter = 0;
	int ssp_train = 100;
	string weights_path = "gewichten";
};

class QLearning : public RLearning {
public:
	QLearning(int t, int tes, int am) { training = t, test = tes, ActionSpaceMode = am; };
	void UpdateQTable(Key& state, Key& nextState, Ray& hit, std::map<Key, QValue>* dict, std::map<int, float3>* dictAction,
		int& old_action, float& BRDF, float3& nl, float3& x, float prob, std::map<StateAction, StateActionCount>* dictStateActionCount);
};

