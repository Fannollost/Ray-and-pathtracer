#include "precomp.h"
// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( accumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );

}

enum MAT_TYPE {
	DIFFUSE = 1,
	METAL = 2,
	GLASS = 3,
};
// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace(Ray& ray, int depth, float3 energy)
{
	if (depth <= 0) return float3(0, 0, 0);
	float t_min = 0.001f;
	scene.FindNearest(ray, t_min);
	if (ray.objIdx == -1) return 0; // or a fancy sky color	
	float3 totCol = float3(0);
	material* m = ray.GetMaterial();
	float3 f = m->col;
	float3 I = ray.O + ray.t * ray.D;
	float3 N = ray.hitNormal;// scene.GetNormal(ray.objIdx, I, ray.D);

	if (!scene.raytracer) {
		double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
		if (depth < 5 || !p) {
			if (RandomFloat() < p) {
				f = f * (1 / p);
			}
			else {
				return totCol;
			}
		}
	}
	switch (m->type) {
	case GLASS:	{
		glass* g = (glass*)m;
		float3 refractionColor = 0;
		// compute fresnel
		float kr;
		g->fresnel(normalize(ray.D), normalize(ray.hitNormal), g->ir, kr);
		bool outside = dot(ray.D, ray.hitNormal) < 0;
		float3 bias = 0.0001f * ray.hitNormal;
		float3 norm = outside ? ray.hitNormal : -ray.hitNormal;
		float r = !outside ? g->ir : (1/g->ir);
		// compute refraction if it is not a case of total internal reflection
		if (outside)
		{
			energy.x *= exp(g->absorption.x * -ray.t);
			energy.y *= exp(g->absorption.y * -ray.t);
			energy.z *= exp(g->absorption.z * -ray.t);
		}
		else {
			energy = energy;
		}
		if (kr < 1) {
			float3 refractionDirection = normalize(g->RefractRay(ray.D, norm,r));
			float3 refractionRayOrig = outside ? ray.IntersectionPoint() - bias : ray.IntersectionPoint() + bias;
			Ray refrRay = Ray(refractionRayOrig, refractionDirection, ray.color);
			refractionColor = g->col * Trace(refrRay, depth - 1, energy);
		}

		float3 reflectionDirection = normalize(reflect(ray.D, norm));
		float3 reflectionRayOrig = outside ? ray.IntersectionPoint() + bias : ray.IntersectionPoint() - bias;
		Ray reflRay = Ray(reflectionRayOrig, reflectionDirection, ray.color);
		float3 reflectionColor = g->col * Trace(reflRay, depth - 1, energy);

		// mix the two
		totCol += reflectionColor * kr + refractionColor * (1 - kr);
		break;
	}
	case METAL:	{
		Ray reflected;
		((metal*)m)->scatter(ray, reflected, N, energy);
		totCol += m->col * Trace(reflected, depth - 1, energy) * energy;
		break;
	}
	case DIFFUSE:
		for (int i = 0; i < sizeof(scene.light) / sizeof(scene.light[0]); i++)
		{
			Ray scattered;
			float3 attenuation;
			float3 lightRayDirection = scene.light[i]->GetLightPosition() - ray.IntersectionPoint();
			float len2 = dot(lightRayDirection, lightRayDirection);
			lightRayDirection = normalize(lightRayDirection);
			Ray r = Ray(ray.IntersectionPoint() + lightRayDirection * 1e-4f ,lightRayDirection, ray.color, sqrt(len2));
			((diffuse*)m)->scatter(ray, attenuation, scattered, normalize(lightRayDirection),
				scene.light[i]->GetLightIntensityAt(ray.IntersectionPoint(), N, *m), N, energy);
			if (scene.IsOccluded(r, t_min)) {
				 if (r.m->type != GLASS)
					continue;
			}

			if (((diffuse*)m)->shinieness != 0)
				totCol += ((diffuse*)m)->shinieness * m->col * Trace(Ray(ray.IntersectionPoint(), reflect(ray.D, ray.hitNormal), ray.color), depth - 1, energy) * energy;

			totCol += (1- ((diffuse*)m)->shinieness) * m->col * attenuation * energy;
		}
	}
	
	return totCol;
}


/*
float3 Renderer::Sample(Ray& ray, int depth, float3 energy) {
	float t;
	float t_min = 0.001f;
	scene.FindNearest(ray, t_min);
	if (ray.objIdx == -1 || depth > 4) return float3(0);
	float3 totCol;
	material* m = ray.GetMaterial();
	float3 objCol = ray.m->col;
	float3 f = m->col;
	float p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
	if (depth < 5) {
		if (RandomFloat() < 1) {
			f = f * (1 / p);
		}
		else {
			return objCol;
		}
	}
	//cout << p;
	//float3 nl = dot(ray.D, ray.hitNormal) < 0 ? ray.hitNormal : ray.hitNormal * -1;
	if (m->type == DIFFUSE) {
		Ray scattered;
		float3 attenuation;
		int randomLight = (int)(Rand(size(scene.light)) - 0.0001f);
		Light* light = scene.light[randomLight];
		float3 lightRayDirection = light->GetLightPosition() - ray.IntersectionPoint();
		float len2 = dot(lightRayDirection, lightRayDirection);
		lightRayDirection = normalize(lightRayDirection);
		Ray r = Ray(ray.IntersectionPoint() + lightRayDirection * 1e-4f, lightRayDirection, ray.color, sqrt(len2));
		((diffuse*)m)->scatter(ray, attenuation, scattered, normalize(lightRayDirection),
			light->GetLightIntensityAt(ray.IntersectionPoint(), ray.hitNormal, *m), ray.hitNormal, energy);
		if (scene.IsOccluded(r, t_min)) {
				return float3(0);
		} 

		//if (!scene.raytracer)
		//	totCol += Trace(scattered, depth - 1, energy) * energy;
		  
		return m->emission + m->col * attenuation * Sample(scattered,depth+1,energy) * energy;
		//Ray r = Ray(ray.IntersectionPoint(), RandomInHemisphere(ray.hitNormal), float3(1));
		//return m->emission + f * Sample(r, depth++, energy);
	}
	else if (m->type == METAL) {
		Ray r = Ray(ray.IntersectionPoint(), reflect(ray.D, ray.hitNormal), float3(1));
		return m->emission + f * Sample(r, depth++, energy);
	}
	else if (m->type == GLASS) {
		/*glass* g = (glass*)m;
		float kr;
		float3 refractionColor;
		g->fresnel(ray.IntersectionPoint(), ray.hitNormal, g->ir, kr);
		bool outside = dot(ray.D, ray.hitNormal) < 0;
		float3 bias = 0.001f * ray.hitNormal;
		float3 absorbance = (1 - g->col) * g->absorption * -ray.t;
		if (outside)
		{
			energy.x *= exp(g->absorption.x * -ray.t);
			energy.y *= exp(g->absorption.y * -ray.t);
			energy.z *= exp(g->absorption.z * -ray.t);
		}
		else {
			energy = energy;
		}
		float odds = 0.25 + kr * 0.5;
		if (RandomFloat() < odds) {
			float3 refractionDirection = normalize(g->RefractRay(ray.D, ray.hitNormal, g->ir));
			float3 refractionRayOrig = outside ? ray.hitNormal - bias : ray.hitNormal + bias;
			Ray refrRay = Ray(refractionRayOrig, refractionDirection, ray.color);
			refractionColor = f * Sample(refrRay, depth++, energy);
			return m->emission + refractionColor;
		} else  {
			float3 reflectionDirection = normalize(reflect(ray.D, ray.hitNormal));
			float3 reflectionRayOrig = outside ? ray.hitNormal + bias : ray.hitNormal - bias;
			Ray reflRay = Ray(reflectionRayOrig, reflectionDirection, ray.color);
			float3 reflectionColor = f * Sample(reflRay, depth++, energy);
			return m->emission + reflectionColor;
		}										 *//*
		glass* g = (glass*)m;
		float3 refractionColor = 0;
		// compute fresnel
		float kr;
		g->fresnel(normalize(ray.D), normalize(ray.hitNormal), g->ir, kr);
		bool outside = dot(ray.D, ray.hitNormal) < 0;
		float3 bias = 0.0001f * ray.hitNormal;
		float3 norm = outside ? ray.hitNormal : -ray.hitNormal;
		float r = !outside ? g->ir : (1/g->ir);
		// compute refraction if it is not a case of total internal reflection
		if (outside)
		{
			energy.x *= exp(g->absorption.x * -ray.t);
			energy.y *= exp(g->absorption.y * -ray.t);
			energy.z *= exp(g->absorption.z * -ray.t);
		}
		else {
			energy = energy;
		}
		float odds = 0.25 + kr * 0.5;
		if (RandomFloat() < odds) {
			float3 refractionDirection = normalize(g->RefractRay(ray.D, norm,r));
			float3 refractionRayOrig = outside ? ray.IntersectionPoint() - bias : ray.IntersectionPoint() + bias;
			Ray refrRay = Ray(refractionRayOrig, refractionDirection, ray.color);
			refractionColor = g->col * Trace(refrRay, depth - 1, energy);
			return refractionColor;
		}
		else{
			float3 reflectionDirection = normalize(reflect(ray.D, norm));
			float3 reflectionRayOrig = outside ? ray.IntersectionPoint() + bias : ray.IntersectionPoint() - bias;
			Ray reflRay = Ray(reflectionRayOrig, reflectionDirection, ray.color);
			float3 reflectionColor = g->col * Trace(reflRay, depth - 1, energy);
			return reflectionColor;
		}
	}
}
*/

float3 Renderer::Sample(Ray& ray, int depth, float3 energy) {
	float t_min = 0.01f;
	float eps = 0.0001f;
	float3 totCol = 0;
	//float3 hitCol = Trace(ray, 1, energy);
	if (depth < 0) return float3(1.0, 1.0, 1.0);
	scene.FindNearest(ray, t_min);
	float3 hemiDir = RandomInHemisphere(ray.hitNormal);
	Ray rayToHemisphere = Ray(ray.IntersectionPoint() + hemiDir * eps, hemiDir, ray.color);
	for (int i = 0; i < size(scene.light); i++)
	{
		scene.light[i]->Intersect(rayToHemisphere, t_min);
		if (rayToHemisphere.objIdx == 11 || rayToHemisphere.objIdx == 12) {
			float3 BRDF = ray.m->albedo * INVPI;
			float cos_i = dot(hemiDir, ray.D);
			return 2.0f * PI * BRDF * scene.light[i]->GetLightColor() * cos_i;
		}
	}

	return float3(0);
}
// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime, int frameNr)
{
	
	// animation
	if (!camera.paused && scene.raytracer) {
		static float animTime = 0;
		scene.SetTime(animTime += deltaTime * 0.002f);
	}

	if(scene.raytracer){
		camera.MoveTick();
		camera.FOVTick();
		camera.aspectTick();
	}
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++) {
			float3 totCol = float3(0);				//antialiassing
			for (int s = 0; s < scene.aaSamples; ++s) {
				if (scene.raytracer) {
					float newX = x + RandomFloat(); //+ random(-1.0f, 1.0f);
					float newY = y + RandomFloat(); //+ random(-1.0f, 1.0f);
					totCol += Trace(camera.GetPrimaryRay(newX, newY), 4, float3(1));
					accumulator[x + y * SCRWIDTH] = (totCol / scene.aaSamples);
				}
				else {
					float newX = x + random(-1.0f, 1.0f);
					float newY = y + random(-1.0f, 1.0f);
					totCol += Sample(camera.GetPrimaryRay(newX, newY),4, float3(1));
					accumulator[x + y * SCRWIDTH] += totCol;
					//cout << totCol.x;
				}
			}
			//accumulator[x + y * SCRWIDTH] = totCol;
		}
		if(scene.raytracer)
			frameNr = 1;
		// translate accumulator contents to rgb32 pixels
		for (int dest = y * SCRWIDTH, x = 0; x < SCRWIDTH; x++)	{
			float4 acc = accumulator[x + y * SCRWIDTH] / frameNr; /// iteration;
			screen->pixels[dest + x] = (RGBF32_to_RGB8(&acc));///iteration ;
		}
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	printf( "%5.2fms (%.1ffps) - %.1fMrays/s %.1fCameraSpeed\n", avg, fps, rps / 1000000, camera.speed );
}

