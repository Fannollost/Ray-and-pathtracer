#pragma once

// default screen resolution
#define SCRWIDTH	600//1280
#define SCRHEIGHT	400//720
// #define FULLSCREEN
// #define DOUBLESIZE

namespace Tmpl8 {

class Camera
{
public:
	Camera()
	{
		// setup a basic view frustum
		camPos = float3( 0, 1, -2 );
		topLeft = float3( -aspect, 2, 0 );
		topRight = float3( aspect, 2, 0 );
		bottomLeft = float3( -aspect, 0, 0 );
		screenCenter = topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
		speed = 0.1f;
	}
	Ray GetPrimaryRay(const int x, const int y)
	{
		if (fishEye) {
			float theta = viewAngle * PI;
			const float u = (float)(x - SCRWIDTH / 2) * (aspect*viewAngle / SCRWIDTH);
			const float v = (float)(y - SCRWIDTH / 2) * (viewAngle / SCRHEIGHT);
			float3 newRay = RotateX(RotateY(normalize(screenCenter - camPos), camPos, -u), camPos, -v);
			return Ray(camPos, normalize(newRay), float3(0));
		}
		else {
			// calculate pixel position on virtual screen plane
			const float u = (float)x * (1.0f / SCRWIDTH);
			const float v = (float)y * (1.0f / SCRHEIGHT);
			const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
			return Ray(camPos, normalize(P - camPos), float3(0));
		}

	}
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float viewAngle = 0.25;
	float3 camPos;
	float3 topLeft, topRight, bottomLeft, screenCenter;
	float speed;
	float yAngle = 0;
	float3 mov = float3(0.f);
	float fovChange = 0.f;
	bool paused = false;
	bool fishEye = false;
	bool changed = false;

	void ToogleFisheye() {
		changed = true;
		fishEye = !fishEye;
	}
	
	void SetChange(bool s) { changed = s; }
	bool GetChange() { return changed; }

	void MoveTick() {
		screenCenter = topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
		float3 velocity = speed * ((mov[2] * normalize(screenCenter - camPos)) + (mov[0] * normalize(topRight - topLeft)) + float3(0, 1,0) * mov[1]);
		if (length(velocity) > 0) {
			camPos += velocity;
			topLeft += velocity;
			topRight += velocity;
			bottomLeft += velocity;
			changed = true;
		}
		screenCenter = topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
	}

	void FOVTick() {
		if (fishEye) {
			if(viewAngle>0.11f && viewAngle<0.99)viewAngle -= 0.01f * fovChange; changed = true;
		}
		else {
			screenCenter = topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
			if (fovChange != 0 && (length(screenCenter - camPos) > 0.1f || fovChange > 0)) {
				topLeft += normalize(screenCenter - camPos) * 0.1f * fovChange;
				topRight += normalize(screenCenter - camPos) * 0.1f * fovChange;
				bottomLeft += normalize(screenCenter - camPos) * 0.1f * fovChange;
				changed = true;
			}
			screenCenter = topLeft + .5f * (topRight - topLeft) + .5f * (bottomLeft - topLeft);
		}
	}


	void MoveCameraZY(int dir, bool maj) {
		if (maj) {
			mov[2] = 0;
			mov[1] += dir;
		} else {
			mov[1] = 0;
			mov[2] += dir;
		}
	}

	void transferYZ() {
		float tmp = mov[1];
		mov [1] = mov[2];
		mov[2] = tmp;
	}

	void MoveCameraX(int dir) {
		mov[0] += dir;
	}

	void RotateScreenX(float theta) {
		topLeft = RotateX(topLeft, camPos, theta);
		topRight = RotateX(topRight, camPos, theta);
		bottomLeft = RotateX(bottomLeft, camPos, theta);
		changed = true;
	}

	void RotateScreenY(float theta) {
		yAngle += theta;
		topLeft = RotateY(topLeft, camPos,theta);
		topRight = RotateY(topRight, camPos, theta);
		bottomLeft = RotateY(bottomLeft, camPos, theta);
		changed = true;
	}

	float3 RotateZ(float3 p, float3 center, float theta) {
		double c = cos(theta);
		double s = sin(theta);
		float3 res = float3(0.f);

		float3 vect = p - center;
		float3 xTransform = float3(c, -s, 0);
		float3 yTransform = float3(s, c, 0);

		res[0] = dot(vect, xTransform);
		res[1] = dot(vect, yTransform);
		res[2] = vect[2];

		return res + center;
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
		vect = RotateY(vect, float3(0.f), -yAngle);
		float3 zTransform = float3(0, -s, c);
		float3 yTransform = float3(0, c, s);

		res[0] = vect[0];
		res[1] = dot(vect, yTransform);
		res[2] = dot(vect, zTransform);

		res = RotateY(res, float3(0.f), yAngle);

		return res + center;
	}

	void FOV(float x) {
		fovChange += x;
	}



	void TogglePause() {
		paused = !paused;
	}
};

}