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
	float3 attenuation;
	Ray scattered;

	material* m = ray.GetMaterial();
	float3 I = ray.O + ray.t * ray.D;
	float3 N = ray.hitNormal;// scene.GetNormal(ray.objIdx, I, ray.D);
	bool s = m->scatter(ray, attenuation, scattered, N, energy);
	for(int i = 0; i < sizeof(scene.light) / sizeof(scene.light[0]); i++){
		if(m->specularity != 1){
			float3 lightRayDirection = scene.light[i]->GetLightPosition() - ray.IntersectionPoint();
			float len = length(lightRayDirection);
			Ray r = Ray(ray.IntersectionPoint(), normalize(lightRayDirection), ray.color, len);
			if (scene.IsOccluded(r, t_min)) continue;
			totCol += (1 - m->specularity) * m->col * scene.light[i]->GetLightIntensityAt(ray.IntersectionPoint(), N, *m) * energy;
		}
	}

	if(m->specularity != 0)
		totCol += m->specularity * attenuation * (Trace(scattered, depth - 1, energy)) * energy;

	return totCol;
}
	//float3 albedo = scene.GetAlbedo(ray.objIdx, I);

	//RandomInHemisphere(N));
	//Ray r = Ray(ray.IntersectionPoint(), normalize(target - ray.IntersectionPoint()), ray.color);

	//return 0.5 * Trace(r, depth - 1);

	//GENERATES SKYBOX
	/*float3 unit_dir = UnitVector(ray.D);
	auto t = 0.5 * (unit_dir.y + 1.0);
	return (1.0 - t) * float3(1.0, 1.0, 1.0) + t * float3(0.5, 0.7, 1.0);  */

	/*scene.FindNearest(r, t_min);
	if (r.t < length(normalize(r.D))) {
		return float3(0, 0, 0);
	}	*/					 


	//return ray.color;
	/* visualize normal */// return (N + 1) * 0.5f;

	/* visualize distance */ // return 0.1f * float3( ray.t, ray.t, ray.t );
	/* visualize albedo */ // return albedo;


// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick( float deltaTime )
{
	
	// animation
	if (!camera.paused) {
		static float animTime = 0;
		scene.SetTime(animTime += deltaTime * 0.002f);
	}
	camera.MoveTick();
	camera.FOVTick();
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
				float newX = x + RandomFloat();
				float newY = y + RandomFloat();
				totCol += Trace(camera.GetPrimaryRay(newX, newY), 6, float3(1));
			}
			accumulator[x + y * SCRWIDTH] = totCol / scene.aaSamples;
		}
		// translate accumulator contents to rgb32 pixels
		for (int dest = y * SCRWIDTH, x = 0; x < SCRWIDTH; x++)
			screen->pixels[dest + x] = 
				RGBF32_to_RGB8( &accumulator[x + y * SCRWIDTH] );
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000 / avg, rps = (SCRWIDTH * SCRHEIGHT) * fps;
	printf( "%5.2fms (%.1fps) - %.1fMrays/s %.1fCameraSpeed\n", avg, fps, rps / 1000000, camera.speed );
}

