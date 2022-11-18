#pragma once
#include <iostream>

using namespace std;
namespace Tmpl8
{

class Renderer : public TheApp
{
public:
	// game flow methods
	void Init();
	float3 Trace( Ray& ray, int depth );
	void Tick( float deltaTime );
	void Shutdown() { /* implement if you want to do something on exit */ }
	// input handling
	void MouseUp( int button ) { 
		/* implement if you want to detect mouse button presses */ }
	void MouseDown( int button ) { /* implement if you want to detect mouse button presses */ }
	void MouseMove( int x, int y ) { mousePos.x = x, mousePos.y = y; }
	void MouseWheel(float y) { 
		y > 0 ? camera.speed += 0.1f : camera.speed -= 0.1f;
	}
	void KeyUp( int key ) {
		
		/* implement if you want to handle keys */ }
	void KeyDown( int key ) { 
		switch (key) {
			case KEYBOARD_W:
				camera.MoveCameraY(1);
				break;
			case KEYBOARD_S:
				camera.MoveCameraY(-1);
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
		}
		/* implement if you want to handle keys */ }
	// data members
	int2 mousePos;
	float4* accumulator;
	Scene scene;
	Camera camera;

	enum UserInput {
		KEYBOARD_W = 87,
		KEYBOARD_D = 68,
		KEYBOARD_S = 83,
		KEYBOARD_A = 65,
		KEYBOARD_SPACE = 32
	};
};

} // namespace Tmpl8