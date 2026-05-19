#pragma once

#include "scene.h"
#include "pcg.h"

#include "microfacet.h"

#define BSSRDF_SAMPLE_NUM 64

#define BSSRDF_SINGLE 1
#define BSSRDF_MULTI 1
#define BSSRDF_SPEC 1

Spectrum L_bssrdf(const Scene& scene, const PathVertex& vertex, const BSSRDF& bssrdf, const Vector3& wo, pcg32_state& rng);
