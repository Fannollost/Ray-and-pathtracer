#pragma once
#include <iostream>
namespace Tmpl8
{

	class Renderer : public TheApp
	{
	public:
		// game flow methods
		void Init();
		float3 Trace(Ray& ray, int depth, float3 energy);
		float3 Sample(Ray& ray, int depth, float3 energy, const int sampleIdx);
		float3 Debug(Ray& ray, float3 totCol);
		HemisphereSampling::Sample SampleDirection(const Ray& hitPoint);
		void Tick(float deltaTime);
		void Shutdown() { /* implement if you want to do something on exit */ }
		// input handling
		void MouseUp(int button) {
			mousePressed = false;
		}
		void MouseDown(int button) { mousePressed = true; }
		void MouseMove(int x, int y) {
			if (mousePressed) {
				camera.RotateScreenY((float(x) - float(mousePos.x)) / SCRHEIGHT);
				camera.RotateScreenX((float(y) - float(mousePos.y)) / SCRWIDTH);
			}
			mousePos.x = x, mousePos.y = y;

		}
		void MouseWheel(float y) {
			y > 0 ? camera.speed += 0.1f : camera.speed -= 0.1f;
		}
		void KeyUp(int key) {
			switch (key) {
			case KEYBOARD_W:
				camera.MoveCameraZY(-1,majPressed);
				break;
			case KEYBOARD_S:
				camera.MoveCameraZY(1,majPressed);
				break;
			case KEYBOARD_D:
				camera.MoveCameraX(-1);
				break;
			case KEYBOARD_A:
				camera.MoveCameraX(1);
				break;
			case KEYBOARD_PLUS:
				camera.FOV(-1.f);
				break;
			case KEYBOARD_MINUS:
				camera.FOV(1.f);
				break;
			case KEYBOARD_MAJ:
				majPressed = false;
				camera.transferYZ();
				break;
			}
		}
		void KeyDown(int key) {
			switch (key) {
			case KEYBOARD_W:
				camera.MoveCameraZY(1, majPressed);
				break;
			case KEYBOARD_S:
				camera.MoveCameraZY(-1, majPressed);
				break;
			case KEYBOARD_D:
				camera.MoveCameraX(1);
				break;
			case KEYBOARD_A:
				camera.MoveCameraX(-1);
				break;
			case KEYBOARD_SPACE:
				camera.TogglePause();
				break;
			case KEYBOARD_PLUS:
				camera.FOV(1.f);
				break;
			case KEYBOARD_MINUS:
				camera.FOV(-1.f);
				break;
			case KEYBOARD_F:
				camera.ToogleFisheye();
				break;
			case KEYBOARD_P:
				camera.SetChange(true);
				scene.toogleRaytracer();
				break;
			case KEYBOARD_MAJ:
				majPressed = true;
				camera.transferYZ();
				break;
			case KEYBOARD_B:
				scene.ExportData();
			}
			/* implement if you want to handle keys */
		}
		// data members
		int2 mousePos;
		bool mousePressed = false;
		bool learningEnabled = false;
		int maxRayDepth = 5;
		float4* accumulator;
		Scene scene;
		Camera camera;
		float2 xBox = float2(-1, 1), yBox = float2(-1, 1), zBox = float2(-1, 1);	//makeboudningbox
		QTable* qTable = new QTable(8,5,0.25f, float3(0,1,1), 7, 0.15f);
		bool majPressed = false;
		enum UserInput {
			KEYBOARD_B = 66,
			KEYBOARD_W = 87,
			KEYBOARD_D = 68,
			KEYBOARD_S = 83,
			KEYBOARD_A = 65,
			KEYBOARD_F = 70,
			KEYBOARD_P = 80,
			KEYBOARD_SPACE = 32,
			KEYBOARD_PLUS = 334,
			KEYBOARD_MINUS = 333,
			KEYBOARD_MAJ = 340
		};
	};
};