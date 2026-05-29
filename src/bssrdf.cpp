#include "bssrdf.h"



enum BSSRDFSampleAxis {
	UAxis,
	VAxis,
	NAxis
};

/*
	3/9/2026:

	1.
	assume we want to sample points (x, y) with prob. pdf(x, y) ~ exp(-falloff * sqrt(x^2 + y^2))   (2d gaussian centered at origin)

	need a normalization c s.t it is a pdf: pdf(x, y) = c * exp(-falloff * sqrt(x^2 + y^2))

	2.
	switch from cartesian to polar corrdinates:
	Jacobian determinant: |J| = r

	pdf is now: pdf(r, theta) = c * r * exp(-falloff * r)

	3. integrate over r and theta to get c: c = falloff / pi

	pdf is now pdf(x, y) = a/pi * exp(-falloff * (x^2 + y^2))

	or

	pdf(r, theta) = a / pi * r * exp(-falloff * r * r)

	4. marginal pdf of r: pdf(r) = 2 * a * r * exp(-falloff * r * r)

	   marginal pdf of theta: pdf(theta) = 1 / (2 * pi)

	5. cdf and inverse transform sampling....

	6. r = sqrt(-ln(u1) / falloff).
	   theta = 2 * pi * u2
*/

// generate a random sample from 2d Gaussian distribution
// This produces samples on an infinite plane.
Vector2 gaussSample2d(const Vector2& rnd_param_uv, Real falloff) {
	float r = sqrt(-log(rnd_param_uv[0]) / falloff);
	float theta = 2 * c_PI * rnd_param_uv[1];
	return Vector2{ r * cos(theta), r * sin(theta) };
}

// with integral upper bound = Rmax
Vector2 gaussSample2d(const Vector2& rnd_param_uv, Real falloff, Real Rmax) {
	float r = log(1 - rnd_param_uv[0] * (1 - exp(-falloff * Rmax * Rmax))) / (-falloff);
	r = sqrt(r);

	float theta = 2 * c_PI * rnd_param_uv[1];
	return Vector2{ r * cos(theta), r * sin(theta) };
}

// pdf(x, y) = a/pi * exp(-falloff * (x^2 + y^2))
// project the 3d point to tangent plane, then eval
Real gaussianSample2dPdf(const Vector3& pCenter, const Vector3& pSample, const Vector3& N, Real falloff) {
	Vector3 d = pSample - pCenter;
	Vector3 proj = d - dot(d, N) * N;  // project to tangent plane
	return falloff / c_PI * exp(-falloff * dot(proj, proj));
}

Real gaussianSample2dPdf(const Vector3& pCenter, const Vector3& pSample, const Vector3& N, Real falloff, Real Rmax) {
	return gaussianSample2dPdf(pCenter, pSample, N, falloff) /
		(1.0f - exp(-falloff * Rmax * Rmax));
}

Real gaussianSample2dPdf(Vector2 sample, Real falloff, Real Rmax) {
	Real deno = 1 / c_PI * falloff * exp(-falloff * (sample.x * sample.x + sample.y * sample.y));
	return deno / (1.0f - exp(-falloff * Rmax * Rmax));
}

Real expSample(Real u, Real falloff) {
	return log(u) / (-falloff);
}

Real expSamplePdf(Real x, Real falloff) {
	return falloff * exp(-falloff * x);
}


Spectrum henyeyGreenstein_phase(Real g, Vector3 &dir_in, Vector3 &dir_out){
	return make_const_spectrum(c_INVFOURPI *
		(1 - g * g) /
		(pow((1 + g * g + 2 * g * dot(dir_in, dir_out)), Real(3) / Real(2))));
}

// constructor
BSSRDF::BSSRDF(Texture<Spectrum> base_color, Texture<Spectrum> normMap, Spectrum sigma_a, Spectrum sigma_s, Real g, Real eta, bool hN)
	: base_color(base_color), normalMap(normMap), sigma_a_(sigma_a), sigma_s_(sigma_s), g_(g), eta_(eta) {
	hasN = hN;

	simgma_s_prime_ = sigma_s_ * (make_const_spectrum(1) - make_const_spectrum(g_));
	sigma_t_prime_ = sigma_a_ + simgma_s_prime_;
	alpha_prime_ = simgma_s_prime_ / sigma_t_prime_;
	sigma_tr_ = sqrt(make_const_spectrum(3) * sigma_a_ * sigma_t_prime_);
	D_ = make_const_spectrum(1) / (make_const_spectrum(3) * sigma_t_prime_);


	// Rd ~ exp(-sigma_tr * r*r) 
	// at center, Rd ~ exp(0) = 1
	// relative contribution compared to center is Rd(r) / Rd(0) = exp(-sigma_tr * r*r)
	// exp(-sigma_tr * r*r) = skip ratio, solve for r: 
	Real skip_ratio = 0.01f;  // skip points that contribute less than 1% of the center point

	Real lum_sigmaTr =
		Real(0.2126) * sigma_tr_[0] +
		Real(0.7152) * sigma_tr_[1] +
		Real(0.0722) * sigma_tr_[2];

	Rmax_ = sqrt(log(skip_ratio) / -lum_sigmaTr);
}

// the average diffuse Fresnel reflectance. It can be solved like this.
Real BSSRDF::F_dr(Real eta) const {
	// Donner. C 2006 Chapter 5
	if (eta < 1.0f) {
		return -0.4399f + 0.7099f / eta - 0.3319f / (eta * eta) +
			0.0636f / (eta * eta * eta);
	}
	else {
		return -1.4399f / (eta * eta) + 0.7099f / eta + 0.6681f +
			0.0636f * eta;
	}
}


// d2: the radical distance squared between 2 points on the surface.
// dipole approximation evaluation
// instead of evaluating S(),  Assume the BSSRDF is radially symmetric and mostly a function of distance: Rd()
Spectrum BSSRDF::Rd(Real d2) const {
	// Donner. C 2006 Chapter 5

	Real fdr = F_dr(eta_);
	Real A = (1 + fdr) / (1 - fdr);  // accounts for Fresnel internal reflection

	Spectrum zr = make_const_spectrum(1) / sigma_t_prime_;  // the distance beneath the surface to the positive dipole light
	Spectrum zv = zr + (make_const_spectrum(4) * A * D_);   // the distance above the surface to the negative dipole light

	Spectrum dr = sqrt(d2 + zr * zr);
	Spectrum dv = sqrt(d2 + zv * zv);

	Spectrum sigma_tr_dr = sigma_tr_ * dr;
	Spectrum sigma_tr_dv = sigma_tr_ * dv;

	// calculate Rd
	Spectrum coef = alpha_prime_ / (4 * c_PI);
	Spectrum exp_trdr = exp(-sigma_tr_dr);
	Spectrum exp_trdv = exp(-sigma_tr_dv);
	Spectrum rd = zr * (sigma_tr_dr + make_const_spectrum(1)) * exp_trdr / (dr * dr * dr) +
		zv * (sigma_tr_dv + make_const_spectrum(1)) * exp_trdv / (dv * dv * dv);
	rd *= coef;

	return rd;
}

// sample a probe ray
// return which axis to sample, the ray sampled, and pdf of this sampled ray
std::tuple<BSSRDFSampleAxis, Ray, Real> sampleProbeRay(const PathVertex& vertex, pcg32_state& rng, Real sigmaTr, Real Rmax, const Scene& scene) {
	Vector2 rnd_uv_gaussian = { next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
	Real rnd_u_axis = next_pcg32_real<Real>(rng);
	Frame frame(vertex.shading_frame.n);			// Culprit for discontinuity triangles, see notes in the end of the file

	// disk distribution sampling
	Vector2 p_sample = gaussSample2d(rnd_uv_gaussian, sigmaTr, Rmax);
	Real pdf_sample = gaussianSample2dPdf(p_sample, sigmaTr, Rmax);

	// probe ray start at a sphere of radius Rmax centered at vertex.position
	// the len of the 3rd dim is calculated here
	Real r2 = p_sample.x * p_sample.x + p_sample.y * p_sample.y;
	Real len_3rd = sqrt(Rmax * Rmax - r2);		// this is half len: on semisphere
	len_3rd = max(len_3rd, Real(0));

	// pick which axis to sample: N , U, or V,  probability ratio: 0.5, 0.25, 0.25
	BSSRDFSampleAxis axis;

	// use -len_3rd such that the ray dir is always pointing along the frame dir
	Vector3 point;
	Vector3 dir;

	Real pdf = pdf_sample;
	// 50% sample N axis
	if (rnd_u_axis <= 0.5) {
		point = vertex.position + to_world(frame, Vector3{ p_sample.x, p_sample.y, -len_3rd });
		axis = NAxis;
		dir = frame.n;  // dir along normal
		pdf *= 0.5;
	}
	// 25% sample U axis
	else if (rnd_u_axis <= 0.75) {
		point = vertex.position + to_world(frame, Vector3{ -len_3rd, p_sample.x, p_sample.y });
		dir = frame.x;  // dir along U axis
		axis = UAxis;
		pdf *= 0.25;
	}
	// 25% sample V axis
	else {
		point = vertex.position + to_world(frame, Vector3{ p_sample.x, -len_3rd, p_sample.y });
		dir = frame.y;  // dir along V axis
		axis = VAxis;
		pdf *= 0.25;
	}

	Ray ray{ point, dir, get_shadow_epsilon(scene), 2 * len_3rd };		// need to tune tnear

	return { axis, ray, pdf };
}

Real MISWeight(const PathVertex& vertex, const Vector3& pIn, const Vector3& nIn, BSSRDFSampleAxis mainAxis, Real sigmaTr, Real pdf, Real Rmax) {
	Vector3 pWo = vertex.position;
	Real weight = 0;

	Frame frame(vertex.shading_frame.n);
	Vector3 N = frame.n;
	Vector3 U = frame.x;
	Vector3 V = frame.y;

	switch (mainAxis)
	{
	case NAxis:
	{
		// dot is the jacobian: from projected disk space to surface area space
		Real uPdf = 0.25 * gaussianSample2dPdf(pWo, pIn, U, sigmaTr, Rmax) * abs(dot(U, nIn));
		Real vPdf = 0.25 * gaussianSample2dPdf(pWo, pIn, V, sigmaTr, Rmax) * abs(dot(V, nIn));
		Real nu = pdf * pdf;
		weight = nu / (nu + uPdf * uPdf + vPdf * vPdf);

		break;
	}
	case VAxis:
	{
		Real nPdf = 0.5 * gaussianSample2dPdf(pWo, pIn, N, sigmaTr, Rmax) * abs(dot(N, nIn));
		Real uPdf2 = 0.25 * gaussianSample2dPdf(pWo, pIn, U, sigmaTr, Rmax) * abs(dot(U, nIn));
		Real nu = pdf * pdf;
		weight = nu / (nPdf * nPdf + uPdf2 * uPdf2 + nu);

		break;
	}
	case UAxis:
	{
		Real nPdf = 0.5 * gaussianSample2dPdf(pWo, pIn, N, sigmaTr, Rmax) * abs(dot(N, nIn));
		Real vPdf = 0.25 * gaussianSample2dPdf(pWo, pIn, V, sigmaTr, Rmax) * abs(dot(V, nIn));
		Real nu = pdf * pdf;
		weight = nu / (nPdf * nPdf + nu + vPdf * vPdf);

		break;
	}
	default:
		break;
	}
	return weight;
}

// single scattering contrib. The light pdf in this estimator is in solid angle space. So no G term inside
// G_term is not that G
/*
	1. when intersect bssrdf mtl, go inside along the refraction dir
	2. sample distance t like in vol_path
	3. sample light source at 2. and test if occluded from other obj
	4. eval
*/
Spectrum L_single(const Scene& scene, const PathVertex& vertex, const BSSRDF& bssrdf, const Vector3& wo, pcg32_state& rng) {
	Vector3 pVertex = vertex.position;
	Vector3 nVertex = vertex.shading_frame.n;
	Real cos_out = abs(dot(wo, nVertex));
	
	Real eta = bssrdf.eta_;	// inter_eta / world_eta
	Real f_trans = 1 - fresnel_dielectric(cos_out, eta);
	Real lum_sigmaTr = 0.2126 * bssrdf.sigma_tr_[0] + 0.7152 * bssrdf.sigma_tr_[1] + 0.0722 * bssrdf.sigma_tr_[2];

	Spectrum sigma_s = bssrdf.sigma_s_;
	Spectrum sigma_t = sigma_s + bssrdf.sigma_a_;

	// refraction
	Vector3 dir_refract = -wo / eta + (fabs(cos_out) / eta - cos_out) * nVertex;
	dir_refract = normalize(dir_refract);

	Spectrum contrib = make_zero_spectrum();

	for (int i = 0; i < BSSRDF_SAMPLE_NUM; i++) {
		Real t = expSample(next_pcg32_real<Real>(rng), lum_sigmaTr);
		Vector3 pSample = pVertex + t * dir_refract;
		Real samplePdf = expSamplePdf(t, lum_sigmaTr);

		// sample light
		Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
		Real light_w = next_pcg32_real<Real>(rng);
		Real shape_w = next_pcg32_real<Real>(rng);
		int light_id = sample_light(scene, light_w);
		const Light& light = scene.lights[light_id];
		PointAndNormal point_on_light =
			sample_point_on_light(light, pSample, light_uv, shape_w, scene);
		Real pdf_light = light_pmf(scene, light_id) *
			pdf_point_on_light(light, point_on_light, pSample, scene);

		Vector3 dir_light;

		Ray shadow_ray;
		if (is_envmap(light)) {
			dir_light = -point_on_light.normal;
			shadow_ray = { pSample, dir_light,
								   get_shadow_epsilon(scene),
								   infinity<Real>() };

			std::optional<PathVertex> isect = intersect(scene, shadow_ray, RayDifferential());
			if (!isect) continue;

			PathVertex inter_vertex = isect.value();		// this is p_wi
			Material mtl = scene.materials[inter_vertex.material_id];

			//if (!std::holds_alternative<BSSRDF>(mtl) || inter_vertex.shape_id != vertex.shape_id) continue;
			if (!std::holds_alternative<BSSRDF>(mtl)) continue;

			// update shadow ray from p_wi
			Vector3 p_wi = inter_vertex.position;
			Vector3 ni = inter_vertex.shading_frame.n;

			Real distanceLeft = infinity<Real>();
			shadow_ray.org = p_wi;
			shadow_ray.tfar = distanceLeft;			// might be wrong here

			// shadow ray
			isect = intersect(scene, shadow_ray, RayDifferential());
			if (isect) continue;

			Spectrum p = henyeyGreenstein_phase(bssrdf.g_, -dir_refract, dir_light);
			Real G_term = abs(dot(ni, dir_refract)) / abs(dot(ni, dir_light));

			Real Fti = 1 - fresnel_dielectric(abs(dot(ni, dir_light)), eta);
			Spectrum sigmaTC = sigma_t + G_term * sigma_t;
			Real di = length(p_wi - pSample);
			Real et = 1.0 / eta;
			Real cosi = abs(dot(dir_light, ni));
			Real diPrime = di * cosi / sqrt(1 - et * et * (1 - cosi * cosi));

			Vector3 expPt1 = exp(-diPrime * sigma_t);
			Vector3 expPt2 = exp(-t * sigma_t);
			Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
			contrib += (f_trans * Fti * p * sigma_s / sigmaTC) * expPt1 * expPt2 * L / (pdf_light * samplePdf);
		}
		// non env map light source
		else {
			dir_light = normalize(point_on_light.position - pSample);
			shadow_ray = { pSample, dir_light,
								   get_shadow_epsilon(scene),
								   (1 - get_shadow_epsilon(scene)) *
									   distance(point_on_light.position, pSample) };

			// convert light pdf to solid angle space:
			pdf_light *= distance_squared(point_on_light.position, pSample) / abs(dot(-dir_light, point_on_light.normal));

			std::optional<PathVertex> isect = intersect(scene, shadow_ray, RayDifferential());
			if (!isect) continue;

			PathVertex inter_vertex = isect.value();		// this is p_wi
			Material mtl = scene.materials[inter_vertex.material_id];

			// if (!std::holds_alternative<BSSRDF>(mtl) || inter_vertex.shape_id != vertex.shape_id) continue;
			if (!std::holds_alternative<BSSRDF>(mtl)) continue;

			// update shadow ray from p_wi
			Vector3 p_wi = inter_vertex.position;
			Vector3 ni = inter_vertex.shading_frame.n;
			
			Real distanceLeft = distance(point_on_light.position, p_wi);
			shadow_ray.org = p_wi;
			shadow_ray.tfar = distanceLeft;			// might be wrong here

			// shadow ray
			isect = intersect(scene, shadow_ray, RayDifferential());
			if (isect) continue;

			Spectrum p = henyeyGreenstein_phase(bssrdf.g_, -dir_refract, dir_light);
			Real G_term = abs(dot(ni, dir_refract)) / abs(dot(ni, dir_light));

			Real Fti = 1 - fresnel_dielectric(abs(dot(ni, dir_light)), eta);
			Spectrum sigmaTC = sigma_t + G_term * sigma_t;
			Real di = length(p_wi - pSample);
			Real et = 1.0 / eta;
			Real cosi = abs(dot(dir_light, ni));
			Real diPrime = di * cosi / sqrt(1 - et * et * (1 - cosi * cosi));

			Vector3 expPt1 = exp(-diPrime * sigma_t);
			Vector3 expPt2 = exp(-t * sigma_t);
			Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
			contrib += (f_trans * Fti * p * sigma_s / sigmaTC) * expPt1 * expPt2 * L / (pdf_light * samplePdf);
		}
	}

	contrib = contrib / (Real)BSSRDF_SAMPLE_NUM;
	return contrib;
}

Spectrum L_diffusion(const Scene& scene, const PathVertex& vertex, const BSSRDF& bssrdf, const Vector3& wo, pcg32_state& rng) {
	Real Rmax = bssrdf.Rmax_;
	Real cos_out = abs(dot(wo, vertex.shading_frame.n));
	Real F_trans = 1 - fresnel_dielectric(cos_out, bssrdf.eta_);

	// luminance == scalar proxy
	Real lum_sigmaTr = 0.2126 * bssrdf.sigma_tr_[0] + 0.7152 * bssrdf.sigma_tr_[1] + 0.0722 * bssrdf.sigma_tr_[2];

	Spectrum contrib = make_const_spectrum(0);
	for (int i = 0; i < BSSRDF_SAMPLE_NUM; i++) {
		auto [spAxis, probRay, pdf_disk_sample] = sampleProbeRay(vertex, rng, lum_sigmaTr, Rmax, scene);

		std::optional<PathVertex> isect = intersect(scene, probRay, RayDifferential());
		if (!isect) continue;

		PathVertex probVertex = isect.value();


		Material mtl = scene.materials[probVertex.material_id];
		// only calc the contrib from the same object
		if (!std::holds_alternative<BSSRDF>(mtl) || probVertex.shape_id != vertex.shape_id) continue;
		//if (!std::holds_alternative<BSSRDF>(mtl)) continue;				// an biased hack

		Vector3 pProbVertex = probVertex.position;
		Spectrum val_Rd = bssrdf.Rd(length_squared(pProbVertex - vertex.position));

		// sample light: need to treat env_map differently with light source

		Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
		Real light_w = next_pcg32_real<Real>(rng);
		Real shape_w = next_pcg32_real<Real>(rng);
		int light_id = sample_light(scene, light_w);
		const Light& light = scene.lights[light_id];
		PointAndNormal point_on_light =
			sample_point_on_light(light, probVertex.position, light_uv, shape_w, scene);
		Real pdf_light = light_pmf(scene, light_id) *
			pdf_point_on_light(light, point_on_light, probVertex.position, scene);

		// shadow ray
		Vector3 dir_light;
		Ray shadow_ray;

		if (is_envmap(light)) {
			dir_light = -point_on_light.normal;

			shadow_ray = { pProbVertex, dir_light,
			   get_shadow_epsilon(scene),
			   infinity<Real>() };
		}
		else {
			dir_light = point_on_light.position - pProbVertex;
			dir_light = normalize(dir_light);
			shadow_ray = { pProbVertex, dir_light,
			   get_shadow_epsilon(scene),
			   (1 - get_shadow_epsilon(scene)) *
				   distance(point_on_light.position, pProbVertex) };
		}
		
		if (occluded(scene, shadow_ray)) {
			continue;
		}

		Vector3 ni = probVertex.shading_frame.n;
		Real cos_in = abs(dot(ni, dir_light));
		Real cos_theta_prime = abs(dot(-dir_light, point_on_light.normal));
		Real Geo;
		if (is_envmap(light)) {
			Geo = 1;	// pdf_light is in solid angle space
		}
		else{
			Geo = cos_in * cos_theta_prime / length_squared(point_on_light.position - pProbVertex);
		}
		
		Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);

		Real F_ti = 1 - fresnel_dielectric(cos_in, bssrdf.eta_);

		// MIS for axis samplings
		Real pdf = pdf_disk_sample * abs(dot(probRay.dir, probVertex.shading_frame.n));		// surface area space
		Real w = MISWeight(vertex, pProbVertex, probVertex.shading_frame.n, spAxis, lum_sigmaTr, pdf, Rmax);

		//Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

		contrib += (w / c_PI * F_trans * F_ti * val_Rd * L * Geo) / (pdf * pdf_light);
	}

	contrib = contrib / (Real)BSSRDF_SAMPLE_NUM;
	return contrib;
}



Spectrum L_bssrdf(const Scene& scene, const PathVertex& vertex, const BSSRDF& bssrdf, const Vector3& wo, pcg32_state& rng) {
	Spectrum final = make_zero_spectrum();

#if BSSRDF_SINGLE == 1
	Spectrum single = L_single(scene, vertex, bssrdf, wo, rng);
	final += single;
#endif 

#if BSSRDF_MULTI == 1
	Spectrum multi = L_diffusion(scene, vertex, bssrdf, wo, rng);
	final += multi;
#endif 

	return final;
}





/*

For BSSRDF diffusion sampling, the probe distribution is radially symmetric, so the tangent orientation should be arbitrary and 
independent of the mesh UV parameterization. 

Using vertex.shading_frame.x/y ties the probe ray placement to per-triangle dpdu/dpdv, 
which can rotate or flip across triangle boundaries even when the shading normal is smooth. 

Rebuilding a frame from only vertex.shading_frame.n keeps the probe plane aligned with the smooth normal 
while removing triangle-local tangent discontinuities from the sampling pattern.

*/