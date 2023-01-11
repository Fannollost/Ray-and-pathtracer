#include "precomp.h"
#include "RLearning.h"

std::map<Key, QValue>* RLearning::LoadWeights(string fileName) {
	std::ifstream in("weights/" + fileName + ".txt");
	string line;
	vector<string> results; 
	std::array<float, 6> temp_state;
	std::array<float, dim_action_space + 1> temp_qvalue;
	std::map<Key, QValue>* dict = new std::map<Key, QValue>;
	while (getline(in, line)) {
		istringstream iss(line);
		while (getline(iss, line, ',')) {
			results.push_back(line);
		}

		for (int i = 0; i < 6; ++i)
			temp_state[i] = (float)std::stod(results[i]);

		for (int i = 0; i < dim_action_space + 1; ++i)
			temp_qvalue[i] = (float)std::stod(results[i + 6]);

		std::map<Key, QValue>& addrDict = *dict;
		addrDict[temp_state] = temp_qvalue;
		results.clear();
	}
	cout << "Weight Loading done!" << endl;
	return dict;
}

float3 RLearning::SamplingScattering(std::map<Key, QValue>* dict, std::map<int, float3>* dictAction, int& id, float3& x, float3& nl,
	StatesStruct& states_rec, std::map<StateAction, StateActionCount>* dictStateActionCount, float& epsilon, const Scene& scene) {
	return float3(0);
}
void RLearning::ComputeQProb(const int& action, const std::array<float, dim_action_space + 1>& qvalue, StatesStruct& states_rec, const float& total) {}
void RLearning::SampleSphereCoordinates(float3& sphereCoord, const float3& pointOldCoord){}
void RLearning::UpdateStateActionCount(const Key& state, const int& action, std::map<StateAction, StateActionCount>* dictStateActionCount,
	bool sequence_complete){}

void QLearning::UpdateQTable(Key& state, Key& nextState, Ray& hit, std::map<Key, QValue>* dict, std::map<int, float3>* dictAction,
	int& old_action, float& BRDF, float3& nl, float3& x, float prob, std::map<StateAction, StateActionCount>* dictStateActionCount){}
