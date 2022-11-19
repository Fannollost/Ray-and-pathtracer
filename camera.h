#pragma once

// default screen resolution
#define SCRWIDTH	1280
#define SCRHEIGHT	720
// #define FULLSCREEN
// #define DOUBLESIZE

namespace Tmpl8 {

class Camera
{
public:
	Camera()
	{
		// setup a basic view frustum
		camPos = float3( 0, 0, -2 );
		topLeft = float3( -aspect, 1, 0 );
		topRight = float3( aspect, 1, 0 );
		bottomLeft = float3( -aspect, -1, 0 );
		speed = 0.1f;
	}
	Ray GetPrimaryRay( const int x, const int y )
	{
		// calculate pixel position on virtual screen plane
		const float u = (float)x * (1.0f / SCRWIDTH);
		const float v = (float)y * (1.0f / SCRHEIGHT);
		const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
		return Ray( camPos, normalize( P - camPos ), float3(0) );
	}
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float3 camPos;
	float3 topLeft, topRight, bottomLeft;
	float speed; 
	bool paused = false;

	void MoveCameraY(int dir) {
		camPos += float3(0, speed * dir, 0);
		topLeft += float3(0, speed * dir, 0);
		topRight += float3(0, speed * dir, 0);
		bottomLeft += float3(0, speed * dir, 0);
	}

	void MoveCameraX(int dir) {
		camPos += float3(speed * dir, 0, 0);
		topLeft += float3(speed * dir, 0, 0);
		topRight += float3(speed * dir, 0, 0);
		bottomLeft += float3(speed * dir, 0, 0);
	}

	void TogglePause() {
		paused = !paused;
	}
};

}