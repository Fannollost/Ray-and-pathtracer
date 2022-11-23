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

	float3 center() {
		return topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
	}

	void MoveCameraY(int dir) {
		camPos += float3(0, speed * dir, 0);
		topLeft += float3(0, speed * dir, 0);
		topRight += float3(0, speed * dir, 0);
		bottomLeft += float3(0, speed * dir, 0);
	}

	void MoveCameraX(int dir) {
		camPos += speed * dir * (topRight - topLeft);
		topLeft += speed * dir * (topRight - topLeft);
		topRight += speed * dir * (topRight - topLeft);
		bottomLeft += speed * dir * (topRight - topLeft);
	}

	void RotateScreenX(float theta) {
		topLeft = RotateX(topLeft, camPos, theta);
		topRight = RotateX(topRight, camPos, theta);
		bottomLeft = RotateX(bottomLeft, camPos, theta);
	}

	void RotateScreenY(float theta) {
		topLeft = RotateY(topLeft, camPos,theta);
		topRight = RotateY(topRight, camPos, theta);
		bottomLeft = RotateY(bottomLeft, camPos, theta);
	}

	float3 RotateY(float3 p, float3 center, float theta) {
		double c = cos(theta);
		double s = sin(theta);
		float3 res = float3(0.f);

		float3 vect = p - center;
		float3 xTransform = float3(c, 0, -s);
		float3 zTransform = float3(s, 0, c);

		res[0] = dot(vect, xTransform);
		res[1] = vect[1];
		res[2] = dot(vect, zTransform);

		return res + center;
	}

	float3 RotateX(float3 p, float3 center, float theta) {
		double c = cos(theta);
		double s = sin(theta);
		float3 res = float3(0.f);

		float3 vect = p - center;
		float3 zTransform = float3(0, -s, c);
		float3 yTransform = float3(0, c, s);

		res[0] = vect[0];
		res[1] = dot(vect, yTransform);
		res[2] = dot(vect, zTransform);

		return res + center;
	}

	void FOV(float x) {
		const float3 c = center();
		if (length(c - camPos) > speed || x > 0) {
			topLeft += (c - camPos) * speed * x;
			topRight += (c - camPos) * speed * x;
			bottomLeft += (c - camPos) * speed * x;
		}
	}



	void TogglePause() {
		paused = !paused;
	}
};

}