#pragma once
// -----------------------------------------------------------
// scene.h
// Simple test scene for ray tracing experiments. Goals:
// - Super-fast scene intersection
// - Easy interface: scene.FindNearest / IsOccluded
// - With normals and albedo: GetNormal / GetAlbedo
// - Area light source (animated), for light transport
// - Primitives can be hit from inside - for dielectrics
// - Can be extended with other primitives and/or a BVH
// - Optionally animated - for temporal experiments
// - Not everything is axis aligned - for cache experiments
// - Can be evaluated at arbitrary time - for motion blur
// - Has some high-frequency details - for filtering
// Some speed tricks that severely affect maintainability
// are enclosed in #ifdef SPEEDTRIX / #endif. Mind these
// if you plan to alter the scene in any way.
// -----------------------------------------------------------

#define SPEEDTRIX

#define PLANE_X(o,i) {if((t=-(ray.O.x+o)*ray.rD.x)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Y(o,i) {if((t=-(ray.O.y+o)*ray.rD.y)<ray.t)ray.t=t,ray.objIdx=i;}
#define PLANE_Z(o,i) {if((t=-(ray.O.z+o)*ray.rD.z)<ray.t)ray.t=t,ray.objIdx=i;}

namespace Tmpl8 {
	class material;
	class diffuse;
	class metal;
	enum MAT_TYPE {
		DIFFUSE = 1,
		METAL = 2,
		GLASS = 3,
	};
	__declspec(align(64)) class Ray
	{
	public:
		Ray() = default;
		Ray(float3 origin, float3 direction, float3 color, float distance = 1e34f)
		{
			exists = true;
			O = origin, D = direction, t = distance;
			// calculate reciprocal ray direction for triangles and AABBs
			rD = float3(1 / D.x, 1 / D.y, 1 / D.z);
#ifdef SPEEDTRIX
			d0 = d1 = d2 = 0;
#endif
		}
		float3 IntersectionPoint() { return O + t * D; }
		void SetMaterial(material* mat) { m = mat; }
		material* GetMaterial() { return m; }
		void SetInside(float3 normal) {
			inside = dot(normalize(D), normal) < 0;
			hitNormal = inside ? normal : -normal;
		}
		// ray data
#ifndef SPEEDTRIX
		float3 O, D, rD;
#else
		union { struct { float3 O; float d0; }; __m128 O4; };
		union { struct { float3 D; float d1; }; __m128 D4; };
		union { struct { float3 rD; float d2; }; __m128 rD4; };
#endif
		float t = 1e34f;
		int objIdx = -1;
		bool inside = false; // true when in medium
		bool exists = false;
		float3 color = 0;
		float3 hitNormal;
		material* m;
	};

	class Light {
	public:
		Light() = default;
		Light(int idx, float3 p, float str, float3 c, float3 n, bool rt)
			: objIdx(idx), pos(p), strength(str), col(c), normal(n), raytracer(rt) {}
		float3 GetNormal() { return normal; }
		virtual float3 GetLightPosition() { return pos; }
		float3 GetLightColor() { return col; }
		virtual float3 GetLightIntensityAt(float3 p, float3 n, const material& m) { return 1; }
		virtual void Intersect(Ray& ray, float t_min) { return; }
		float3 pos;
		bool raytracer;
		float3 col;
		float strength;
		int objIdx;
		float3 normal;
	};

	class AreaLight : public Light {
	public:
		AreaLight() = default;
		AreaLight(int idx, float3 p, float str, float3 c, float r, float3 n, int s, bool rt)
			: Light(idx, p, str, c, n, rt) {
			radius = r;
			radius2 = r * r;
			samples = s;
			area = 2 * radius2 * PI;
		}

		void Intersect(Ray& ray, float t_min) override {
			float d = dot(normal, ray.D);
			//cout << ray.O.x << ", " << ray.O.y << ", " << ray.O.z << endl;
			////cout << "normal" <<  ray.D.x << ", " << ray.D.y << ", " << ray.D.z << endl;
			float3 dir = pos - ray.O;
			float t = dot(dir, normal) / d;
				if (t >= t_min) {
					float3 intersection = ray.O +ray.D * t;
					float3 v = intersection - pos;
					float dis2 = dot(v, v);
					if (sqrtf(dis2) <= radius) {
						ray.t = t, ray.SetInside(normal), ray.color = col;
							ray.objIdx = objIdx;
					}
			}


		}
		float3 GetLightIntensityAt(float3 p, float3 n, const material& m) override {
			float dis = abs(length(pos - p));
			float relStr = 1 / (dis) * strength;
			float3 dir = pos - p;
			float str = dot(n, normalize(dir));
			if (str < 0.0f) str = 0.0f;

			for (int i = 0; i < samples; i++) {
				//lerp between values then divide by samples  
			}
			return relStr * str * GetLightColor();
		}
		float3 GetLightPosition() override {
			if (raytracer) return pos;
			float newRad = radius * sqrt(RandomFloat());
			float theta = random(-1.0f, 1.0f) * 2 * PI;
			return float3(pos.x + newRad * cos(theta), pos.y + newRad * sin(theta), pos.z);
		}
		int samples;
		float radius, radius2;
		float area;
	};

	class DirectionalLight : public Light {
	public:
		DirectionalLight() = default;
		DirectionalLight(int idx, float3 p, float str, float3 c, float3 n, float r, bool rt) : Light(idx, p, str, c, n, rt) {
			sinAngle = sin(r * PI / 2);
		}
		void Intersect(Ray& ray, float t_min) override {
			return;
		}
		float3 GetLightPosition() override {
			return pos;
		}
		float3 GetLightIntensityAt(float3 p, float3 n, const material& m) override {
			float3 dir = p - pos;
			float sTheta = length(cross(dir, normal)) / length(dir) * length(normal);
			if (dot(dir, normal) < 0) {
				return 0;
			}
			float dis = length(dir);
			float str = sinAngle - sTheta > 0 ? asin(sinAngle) - asin(sTheta) : 0;
			
			//return str > 0 ? 1 : 0;

			return 1 / dis * str * strength;
		}

		float sinAngle;
	};

	class Object {
	public:
		Object() = default;
		Object(int idx, material* m) : objIdx(idx), mat(m) {}
		float3 GetIndex() { return objIdx; }
		material* GetMaterial() { return mat; }
		virtual void Intersect(Ray& ray, float t_min) { return; }
		virtual float3 GetNormal(const float3 I) { return float3(); }
		virtual float3 GetAlbedo(const float3 I) { return float3(); }
		int objIdx = -1;
		material* mat;
	};

	// -----------------------------------------------------------
	// Triangle Primitive
	// 
	// 
	// -----------------------------------------------------------
	class Triangle : public Object {
	public:
		Triangle() = default;
		Triangle(int idx, material* m, float3 ver0, float3 ver1, float3 ver2) : objIdx(idx), v0(ver0), v1(ver1), v2(ver2), mat(m) {
			e1 = v1 - v0;
			e2 = v2 - v0;
			N = normalize(cross(e1, e2));
		}
		void Intersect(Ray& ray, float t_min) override {		 //scratchapixel implementation
			float NdotRayDir = dot(N, ray.D);
			if (fabs(NdotRayDir) < t_min) return;
			float d = -dot(N, v0);
			float t = -(dot(N, ray.O) + d) / NdotRayDir;
			if (t < 0) return;
			float3 p = ray.O + t * ray.D;
			float3 c ;

			float3 vp0 = p - v0;
			c = cross(e1, vp0);
			if (dot(N, c) < 0) return ;
			float3 vp1 = p - v1;
			float3 e3 = v2 - v1;
			c = cross(e3, vp1);
			if (dot(N, c) < 0) return;
			float3 e4 = v0 - v2;
			float3 vp2 = p - v2;
			c = cross(e4, vp2);
			if (dot(N, c) < 0) return;

			if (t < ray.t && t > t_min) {
				ray.t = t, ray.objIdx = objIdx, ray.m = mat,
					ray.SetInside(GetNormal(ray.IntersectionPoint()));
			}


		}
		float3 GetNormal(const float3 I) const { return N; }
		float3 v0, v1, v2, e1, e2, N;
		int objIdx = -1;
		float3 col;
		material* mat;
	};

	// -----------------------------------------------------------
	// Sphere primitive
	// Basic sphere, with explicit support for rays that start
	// inside it. Good candidate for a dielectric material.
	// -----------------------------------------------------------
	class Mesh : public Object {
	public:
		Mesh() = default;
		Mesh(int idx,  material* m, const char* path, float3 pos, float scale) : objIdx(idx), mat(m), meshPath(path) {
			ifstream file(meshPath, ios::in);
			if (!file)
			{
				std::cerr << "Cannot open " << meshPath << std::endl;
				exit(1);
			}
			string line;
			float x, y, z;
			while (getline(file, line))
			{
				if (line.substr(0, 2) == "v ") {
					istringstream v(line.substr(2));
					v >> x; v >> y; v >> z;
					vertices.push_back(float3(x*scale+pos.x, y * scale +pos.y, z * scale +pos.z));
				}
				else if (line.substr(0, 2) == "f ") {
					int v0, v1, v2;
					int temp;
					const char* constL = line.c_str();
					sscanf(constL, "f %i//%i %i//%i %i//%i", &v0, &temp, &v1, &temp, &v2, &temp);
					faces.push_back(int3(v0, v1, v2));
				}
			}
			for (int i = 0; i < faces.size(); i++) {
				triangles.push_back(Triangle(objIdx*1000+i, mat, vertices[(faces[i] - 1).x] , vertices[(faces[i] - 1).y] , vertices[(faces[i] - 1).z]));
			}
		}

		void Intersect(Ray& ray, float t_min) override {
			for (int i = 0; i < triangles.size(); i++) {
				triangles[i].Intersect(ray, t_min);
			}
			//triangles[0].Intersect(ray, t_min);
		}
		float3 GetNormal(const float3 I) const
		{
			float3 res;
			int n = 0;
			for (int i = 0; i < triangles.size(); i++) {
				float3 p0 = I - triangles[i].v0;
				float e1Coord = dot(triangles[i].e1, p0);
				float e2Coord = dot(triangles[i].e2, p0);
				if (fabs(dot(triangles[i].N, p0)) < 1 / LARGE_FLOAT && e1Coord > 0 && e1Coord < 1 && e2Coord > 0 && e2Coord < 1) {
					res = triangles[i].N;
					n++;
				}	
			}
			return res;
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(1.f);
		}

		const char* meshPath;
		float3 col = 0;
		int objIdx = -1;
		material* mat;
		int size = 0;
		int vertexNb = 0;
		int facesNb = 0;
		vector<float3> vertices;
		vector<int3> faces;
		vector<Triangle> triangles;
	};

	// -----------------------------------------------------------
	// Sphere primitive
	// Basic sphere, with explicit support for rays that start
	// inside it. Good candidate for a dielectric material.
	// -----------------------------------------------------------
	class Sphere : public Object {
	public:
		Sphere() = default;
		Sphere(int idx, material* m, float3 p, float r) : pos(p), r2(r* r), invr(1 / r), objIdx(idx), mat(m) {}
		void Intersect(Ray& ray, float t_min) override {
			float3 oc = ray.O - this->pos;
			float b = dot(oc, ray.D);
			float c = dot(oc, oc) - this->r2;
			float t, d = b * b - c;
			if (d <= 0) return;
			d = sqrtf(d), t = -b - d;
			if (t < ray.t && t > t_min)
			{
				ray.t = t, ray.objIdx = objIdx, ray.m = mat;
				ray.SetInside(GetNormal(ray.IntersectionPoint()));
				return;
			}
			t = d - b;
			if (t < ray.t && t > t_min)
			{
				ray.t = t, ray.objIdx = objIdx, ray.m = mat;
				ray.SetInside(GetNormal(ray.IntersectionPoint()));
				return;
			}
		}
		float3 GetNormal(const float3 I) const
		{
			return (I - this->pos) * invr;
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(0.93f);
		}
		float3 pos = 0;
		float r2 = 0, invr = 0;
		int objIdx = -1;
		material* mat;
	};

	// -----------------------------------------------------------
	// Plane primitive
	// Basic infinite plane, defined by a normal and a distance
	// from the origin (in the direction of the normal).
	// -----------------------------------------------------------
	class Plane : public Object {
	public:
		Plane() = default;
		Plane(int idx, material* m, float3 normal, float dist) : N(normal), d(dist), objIdx(idx), mat(m) {}
		void Intersect(Ray& ray, float t_min) const
		{
			float t = -(dot(ray.O, this->N) + this->d) / (dot(ray.D, this->N));
			if (t < ray.t && t > t_min) ray.t = t, ray.objIdx = objIdx, ray.m = mat,
				ray.SetInside(N);
		}
		float3 GetNormal(const float3 I) const
		{
			return N;
		}
		float3 GetAlbedo(const float3 I) const
		{
			if (N.y == 1)
			{
				// floor albedo: checkerboard
				int ix = (int)(I.x * 2 + 96.01f);
				int iz = (int)(I.z * 2 + 96.01f);
				// add deliberate aliasing to two tile
				if (ix == 98 && iz == 98) ix = (int)(I.x * 32.01f), iz = (int)(I.z * 32.01f);
				if (ix == 94 && iz == 98) ix = (int)(I.x * 64.01f), iz = (int)(I.z * 64.01f);
				return float3(((ix + iz) & 1) ? 1 : 0.3f);
			}
			else if (N.z == -1)
			{
				// back wall: logo
				static Surface logo("assets/logo.png");
				int ix = (int)((I.x + 4) * (128.0f / 8));
				int iy = (int)((2 - I.y) * (64.0f / 3));
				uint p = logo.pixels[(ix & 127) + (iy & 63) * 128];
				uint3 i3((p >> 16) & 255, (p >> 8) & 255, p & 255);
				return float3(i3) * (1.0f / 255.0f);
			}
			return float3(0.93f);
		}
		float3 N;
		float d;
		int objIdx = -1;
		material* mat;
	};

	// -----------------------------------------------------------
	// Cube primitive
	// Oriented cube. Unsure if this will also work for rays that
	// start inside it; maybe not the best candidate for testing
	// dielectrics.
	// -----------------------------------------------------------
	class Cube : public Object {
	public:
		Cube() = default;
		Cube(int idx, material* m, float3 pos, float3 size, mat4 transform = mat4::Identity())
		{
			objIdx = idx;
			b[0] = pos - 0.5f * size, b[1] = pos + 0.5f * size;
			mat = m;
			M = transform, invM = transform.FastInvertedTransformNoScale();
		}
		void Intersect(Ray& ray, float t_min) override {
			// 'rotate' the cube by transforming the ray into object space
			// using the inverse of the cube transform.
			float3 O = TransformPosition(ray.O, invM);
			float3 D = TransformVector(ray.D, invM);
			float rDx = 1 / D.x, rDy = 1 / D.y, rDz = 1 / D.z;
			int signx = D.x < 0, signy = D.y < 0, signz = D.z < 0;
			float tmin = (b[signx].x - O.x) * rDx;
			float tmax = (b[1 - signx].x - O.x) * rDx;
			float tymin = (b[signy].y - O.y) * rDy;
			float tymax = (b[1 - signy].y - O.y) * rDy;
			if (tmin > tymax || tymin > tmax) return;
			tmin = max(tmin, tymin), tmax = min(tmax, tymax);
			float tzmin = (b[signz].z - O.z) * rDz;
			float tzmax = (b[1 - signz].z - O.z) * rDz;
			if (tmin > tzmax || tzmin > tmax) return;
			tmin = max(tmin, tzmin), tmax = min(tmax, tzmax);
			if (tmin > t_min)
			{
				if (tmin < ray.t) ray.t = tmin, ray.objIdx = objIdx, ray.m = mat,
					ray.SetInside(GetNormal(ray.IntersectionPoint()));
			}
			else if (tmax > t_min)
			{
				if (tmax < ray.t) ray.t = tmax, ray.objIdx = objIdx, ray.m = mat,
					ray.SetInside(GetNormal(ray.IntersectionPoint()));
			}
		}
		float3 GetNormal(const float3 I) const
		{
			// transform intersection point to object space
			float3 objI = TransformPosition(I, invM);
			// determine normal in object space
			float3 N = float3(-1, 0, 0);
			float d0 = fabs(objI.x - b[0].x), d1 = fabs(objI.x - b[1].x);
			float d2 = fabs(objI.y - b[0].y), d3 = fabs(objI.y - b[1].y);
			float d4 = fabs(objI.z - b[0].z), d5 = fabs(objI.z - b[1].z);
			float minDist = d0;
			if (d1 < minDist) minDist = d1, N.x = 1;
			if (d2 < minDist) minDist = d2, N = float3(0, -1, 0);
			if (d3 < minDist) minDist = d3, N = float3(0, 1, 0);
			if (d4 < minDist) minDist = d4, N = float3(0, 0, -1);
			if (d5 < minDist) minDist = d5, N = float3(0, 0, 1);
			// return normal in world space
			return TransformVector(N, M);
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(1, 1, 1);
		}
		float3 b[2];
		mat4 M, invM;
		int objIdx = -1;
		material* mat;
	};

	// -----------------------------------------------------------
	// Quad primitive
	// Oriented quad, intended to be used as a light source.
	// -----------------------------------------------------------
	class Quad {
	public:
		Quad() = default;
		Quad(int idx, material* m, float s, mat4 transform = mat4::Identity())
		{
			objIdx = idx;
			size = s * 0.5f;
			mat = m;
			T = transform, invT = transform.FastInvertedTransformNoScale();
		}
		void Intersect(Ray& ray, float t_min) const {
			const float3 O = TransformPosition(ray.O, invT);
			const float3 D = TransformVector(ray.D, invT);
			const float t = O.y / -D.y;
			if (t < ray.t && t > t_min)
			{
				float3 I = O + t * D;
				if (I.x > -size && I.x < size && I.z > -size && I.z < size)
					ray.t = t, ray.objIdx = objIdx, ray.m = mat,
					ray.SetInside(GetNormal(ray.IntersectionPoint()));
			}
		}
		float3 GetNormal(const float3 I) const
		{
			// TransformVector( float3( 0, -1, 0 ), T ) 
			return float3(-T.cell[1], -T.cell[5], -T.cell[9]);
		}
		float3 GetAlbedo(const float3 I) const
		{
			return float3(10);
		}
		float size;
		mat4 T, invT;
		int objIdx = -1;
		material* mat;
	};

	class material {
	public:
		material(float3 c, bool rt) : col(c), raytracer(rt) {}

		void SetColor(float3 c) { col = c; }
		float3 col;
		int type;
		bool raytracer;
	};

	class diffuse : public material {
	public:
		diffuse(float3 a = 0, float3 c = 0, float ks = 0.2, float kd = 0.8, int n = 2, bool rt = true)
			: specu(ks), diffu(kd), N(n), albedo(a), material(c, rt) {
			type = DIFFUSE;
		}
		void SetSpecularity(float ks) { specu = ks; }
		void SetDiffuse(float kd) { diffu = kd; }
		void SetN(int n) { N = n; }
		virtual bool scatter(Ray& ray, float3& att, Ray& scattered, float3 lightDir, float3 lightIntensity, float3 normal, float3& energy) {
			float3 reflectionDirection = reflect(-lightDir, normal);
			float3 specularColor, lightAttenuation;
			specularColor = powf(fmax(0.0f, -dot(reflectionDirection, ray.D)), N) * lightIntensity;
			lightAttenuation = lightIntensity;

			att = albedo * lightAttenuation * diffu + specularColor * specu;
			float3 dir;
			if (raytracer) {
				dir = ray.IntersectionPoint() + normal;
			}
			else {
				dir = ray.IntersectionPoint() + normal + RandomInHemisphere(normal);
			}
			if (isZero(dir)) dir = normal;
			scattered = Ray(ray.IntersectionPoint(), normalize(dir - ray.IntersectionPoint()), ray.color);
			//att = albedo;  */
			float3 retention = float3(1) - albedo;
			float3 newEnergy(energy - retention);
			energy = newEnergy.x > 0 ? newEnergy : 0;
			return true;
		}

	public:
		float3 albedo;
		float specu, diffu;
		int N;
	};

	class metal : public material {
	public:
		metal(float f, float3 c, bool rt) : fuzzy(f < 1 ? f : 1), material(c, rt) { type = METAL; }
		virtual bool scatter(Ray& ray, Ray& reflected, float3 normal, float3& energy) {
			float3 dir = reflect(ray.D, normal);
			reflected = Ray(ray.IntersectionPoint(), dir, ray.color * col);
			energy = energy;
			return dot(reflected.D, normal) > 0;
		}

	public:
		float fuzzy;
	};

	class glass : public material {
	public:
		glass(float refIndex, float3 c, float3 a, bool rt) : ir(refIndex), absorption(a), material(c, rt) { type = GLASS; }
		virtual bool scatter(Ray& ray, Ray& scattered, Ray& refracted, float3 normal, float3& energy) {
			float kR;
			float refrRatio = ray.inside ? (1.0 / ir) : ir;
			fresnel(ray.IntersectionPoint(), normal, refrRatio, kR);
			float3 uDir = UnitVector(ray.D);
			float ctheta = fmin(dot(-uDir, ray.hitNormal), 1.0);			//maybe use normal
			float stheta = sqrt(1.0 - ctheta * ctheta);
			float3 refrDir, reflDir;
			bool reflected = (refrRatio * stheta) > 1.0;

			float3 absorbance = (1 - col) * absorption * -ray.t;
			if (ray.inside)
			{
				energy.x *= exp(absorbance.x);
				energy.y *= exp(absorbance.y);
				energy.z *= exp(absorbance.z);
			}
			else {
				energy = energy;
			}
			if (kR < 1) {
				refrDir = normalize(refractRay(uDir, ray.hitNormal, refrRatio));
				float3 refrOrig = !ray.inside ? ray.IntersectionPoint() - 0.001f : ray.IntersectionPoint() + 0.001f;
				refracted = Ray(refrOrig, refrDir, col); //check for color of sphere itself
			}
			kr = kR;
			reflDir = normalize(reflect(uDir, ray.hitNormal));
			float3 reflOrig = ray.inside ? ray.IntersectionPoint() - 0.001f : ray.IntersectionPoint() + 0.001f;
			scattered = Ray(reflOrig, reflDir, col);
			return true;
		}
		float3 refractRay(float3 oRayDir, float3 normal, float refRatio) {
			float theta = fmin(dot(-oRayDir, normal), 1.0);
			float3 perpendicular = refRatio * (oRayDir + theta * normal);
			float3	parallel = -sqrt(fabs(1.0 - pow(length(perpendicular), 2))) * normal;
			return perpendicular + parallel;
		} 
			/*float3 refractRay(float3 point, float3 normal, float ir) {
			float cosi = clamp(dot(point, normal ),-1.0f, 1.0f);
			float etai = 1, etat = ir;
			float3 n = normal;
			if (cosi < 0) { cosi = -cosi; }
			else { std::swap(etai, etat); n = -normal; }
			float eta = etai / etat;
			float k = 1 - eta * eta * (1 - cosi * cosi);
			return k < 0 ? 0 : eta * point + (eta * cosi - sqrtf(k)) * n;
		}  */
		float ir;
		float kr;

		//gotta check this/
		void fresnel(float3 I, float3 normal, float ior, float& kr) {
			float cosi = clamp(-1.0f, 1.0f, dot(I, normal));
			float etai = 1, etat = ior;
			if (cosi > 0) { std::swap(etai, etat); }
			float sint = etai / etat * sqrtf(std::max(0.0f, 1 - cosi * cosi));
			if (sint >= 1) {
				kr = 1;
			}
			else {
				float cost = sqrtf(std::max(0.0f, 1 - sint * sint));
				cosi = fabsf(cosi);
				float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
				float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
				kr = (Rs * Rs + Rp * Rp) / 2;
			}

		}//scratchapixel

		float3 absorption;

	};

	// -----------------------------------------------------------
	// Scene class
	// We intersect this. The query is internally forwarded to the
	// list of primitives, so that the nearest hit can be returned.
	// For this hit (distance, obj id), we can query the normal and
	// albedo.
	// -----------------------------------------------------------
	class Scene
	{
	public:
		Scene()
		{
			float3 white = float3(1.0, 1.0, 1.0);
			float3 red = float3(199, 70, 123) / 255;
			float3 blue = float3(112, 66, 219) / 255;
			float3 babyblue = float3(0, 0, 1.0f);
			float3 green = float3(105, 250, 144) / 255;
			diffuse* blueDiff = new diffuse(float3(0.8f), blue, 0.3f, 0.7f, 1200);
			diffuse* standardDiff = new diffuse(float3(0.8f), white, 0.2, 0.8f, 2, raytracer);
			glass* standardGlass = new glass(1.5f, white, float3(0.0f), raytracer);
			glass* blueGlass = new glass(1.5f, blue, float3(0.00f), raytracer);
			glass* diamond = new glass(2.4f, white, float3(0.00f), raytracer);
			diffuse* specularDiff = new diffuse(float3(0.8f), red, 0.3f, 0.7f, 7, raytracer);
			diffuse* greenDiff = new diffuse(float3(0.8f), green, 0.6f, 0.4f, 1250, raytracer);
			metal* standardMetal = new metal(0.7f, white, raytracer);
			// we store all primitives in one continuous buffer

			//light[0] = new DirectionalLight(11, float3(0, 2, 0), 8.0f, white, float3(0, -1, 1), 0.9, raytracer);			//DIT FF CHECKEN!
			light[1] = new AreaLight(12, float3(0,1,0), 3.0f, white, 0.1f, float3(0, -1, 0), 4, raytracer);
			light[0] = new AreaLight(11, float3(0.1f,-0.9f, 0), 3.0f, white, 0.1f, float3(0, 1, 0), 4, raytracer);			//DIT FF CHECKEN!
			//light[2] = new AreaLight(13, float3(0.1f, -1, 0), 2.0f, white, 0.1f, float3(0, -1, 0), 4, raytracer);			//DIT FF CHECKEN!

			plane[0] = Plane(0, specularDiff, float3(1, 0, 0), 3);			// 0: left wall
			plane[1] = Plane(1, new diffuse(0.8f, red, 0), float3(-1, 0, 0), 2.99f);		// 1: right wall
			plane[2] = Plane(2, new diffuse(0.8f, white, 0), float3(0, 1, 0), 1);			// 2: floor
			plane[3] = Plane(3, new diffuse(0.8f, white, 0), float3(0, -1, 0), 2);			// 3: ceiling
			plane[4] = Plane(4, new diffuse(0.8f, red, 0), float3(0, 0, 1), 3);			// 4: front wall
			plane[5] = Plane(5, new diffuse(0.8f, green, 0), float3(0, 0, -1), 3.99f);		// 5: back wall
			//quad = Quad(6, new diffuse(0.8f, white, 0), 1);							// 6: light source

			obj[0] = new Sphere(7, standardMetal, float3(0), 0.5f);			// 1: bouncing ball
			//obj[0] = new Sphere(7, red, new metal(1.0f, 1.0f), float3(-1.5f, 0, 2), 0.5f);		// 1: static ball => set animOn to false
			obj[1] = new Sphere(8, specularDiff, float3(0, 2.5f, -3.07f), 8);		// 2: rounded corners
			//obj[2] = new Sphere(9, white, new glass(0.1f), float3(1.5f, 0, 2), 0.5f);			// 3: static glass sphere => set animOn to false
			obj[2] = new Cube(9, blueDiff, float3(0), float3(1.15f));		// 3: spinning cube
			obj[3] = new Mesh(10, standardMetal, "shape.obj", float3(0,0,2), 0.5f);
			
			//obj[3] = new Triangle(10, new diffuse(0.8f, blue, 0), float3(0.0f, 0.0f, 1.0f), float3(0.2f, 0, 1.0f), float3(0.2f, 0.2f, 1.0f));	// 4: Triangle

			SetTime(0);
			// Note: once we have triangle support we should get rid of the class
			// hierarchy: virtuals reduce performance somewhat.
		}
		void SetTime(float t)
		{
			// default time for the scene is simply 0. Updating/ the time per frame 
			// enables animation. Updating it per ray can be used for motion blur.
			// light source animation: swing
			animTime = t;

			mat4 M1base = mat4::Translate(float3(0, 2.6f, 2));
			mat4 M1 = M1base * mat4::RotateZ(sinf(animTime * 0.6f) * 0.1f) * mat4::Translate(float3(0, -0.9, 0));
			quad.T = M1, quad.invT = M1.FastInvertedTransformNoScale();

			if (animOn) {
				// cube animation: spin
				mat4 M2base = mat4::RotateX(PI / 4) * mat4::RotateZ(PI / 4);
				mat4 M2 = mat4::Translate(float3(1.4f, 0, 2)) * mat4::RotateY(animTime * 0.5f) * M2base;
				if (size(obj) >= 3)((Cube*)obj[2])->M = M2, ((Cube*)obj[2])->invM = M2.FastInvertedTransformNoScale();

				// sphere animation: bounce
				float tm = 1 - sqrf(fmodf(animTime, 2.0f) - 1);
				if (size(obj) >= 1)((Sphere*)obj[0])->pos = float3(-1.4f, -0.5f + tm, 2);
			}

		}
		float3 GetLightPos() const
		{
			// light point position is the middle of the swinging quad
			float3 corner1 = TransformPosition(float3(-0.5f, 0, -0.5f), quad.T);
			float3 corner2 = TransformPosition(float3(0.5f, 0, 0.5f), quad.T);
			return (corner1 + corner2) * 0.5f - float3(0, 0.01f, 0);
		}
		float3 GetLightColor() const
		{
			return quad.mat->col;
		}
		void FindNearest(Ray& ray, float t_min) const
		{

			for (int i = 0; i < size(plane); i++) plane[i].Intersect(ray, t_min);

			for (int i = 0; i < size(obj); i++) obj[i]->Intersect(ray, t_min);
		}
		bool IsOccluded(Ray& ray, float t_min) const
		{
			float rayLength = ray.t;
			// skip planes: it is not possible for the walls to occlude anything
			quad.Intersect(ray, t_min);
			for (int i = 0; i < size(obj); i++) obj[i]->Intersect(ray, t_min);

			return ray.t < rayLength;
			// technically this is wasteful: 
			// - we potentially search beyond rayLength
			// - we store objIdx and t when we just need a yes/no
			// - we don't 'early out' after the first occlusion
		}
		float3 GetNormal(int objIdx, float3 I, float3 wo) const
		{
			// we get the normal after finding the nearest intersection:
			// this way we prevent calculating it multiple times.
			if (objIdx == -1) return float3(0); // or perhaps we should just crash
			float3 N;
			if (objIdx == 6) N = quad.GetNormal(I);
			else if (objIdx >= 7 && objIdx < 7 + size(obj)) N = obj[objIdx - 7]->GetNormal(I);
			else
			{
				// faster to handle the 6 planes without a call to GetNormal
				N = float3(0);
				N[objIdx / 2] = 1 - 2 * (float)(objIdx + 4);
			}
			if (dot(N, wo) > 0) N = -N; // hit backside / inside
			return N;
		}
		float3 GetAlbedo(int objIdx, float3 I) const
		{
			if (objIdx == -1) return float3(0); // or perhaps we should just crash
			if (objIdx == 0) return quad.GetAlbedo(I);
			if (objIdx >= 7 && objIdx < 7 + size(obj)) return obj[objIdx - 7]->GetAlbedo(I);
			return plane[objIdx].GetAlbedo(I);
			// once we have triangle support, we should pass objIdx and the bary-
			// centric coordinates of the hit, instead of the intersection location.
		}
		float GetReflectivity(int objIdx, float3 I) const
		{
			if (objIdx == 1 /* ball */) return 1;
			if (objIdx == 6 /* floor */) return 0.3f;
			return 0;
		}
		float GetRefractivity(int objIdx, float3 I) const
		{
			return objIdx == 3 ? 1.0f : 0.0f;
		}
		__declspec(align(64)) // start a new cacheline here
			float animTime = 0;

		Light* light[2];
		Object* obj[4];
		Quad quad;
		Plane plane[6];
		int aaSamples = 1;
		bool raytracer = false;
		float mediumIr = 1.0f;
		bool animOn = true; // set to false while debugging to prevent some cast error from primitive object type
	};
}